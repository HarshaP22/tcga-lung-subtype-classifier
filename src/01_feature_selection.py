"""Step 1: five feature-selection methods + overlap.

Reproduces the paper's selection stage: DGE, PCA, mRMR, XGBoost, Lasso, each
producing an independent ranked/selected gene set; then the Venn intersection
-> genes in >=2 methods (paper: 131) and >=3 methods (paper: 17).

A single fixed patient-level 7:3 split is created here and reused by every
downstream step (paper: "we split the training sets and testing sets the same
way for all validations").

Outputs (data/processed, results/tables)
-------
processed/split.json          {train:[samples], test:[samples]}
processed/dge_stats.csv        per-gene logFC (LUSC vs LUAD), p, FDR
results/tables/feature_sets.json   {method: [genes]}, plus overlap131/overlap17
results/tables/overlap_membership.csv  gene x method 0/1 + n_methods
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent))
import config as C
import utils as U


def patient_split(labels: pd.DataFrame, samples, seed=C.SEED, test_size=C.TEST_SIZE):
    """Stratified-by-histology, grouped-by-patient train/test split."""
    lab = labels.set_index("sample").loc[list(samples)]
    rng = np.random.default_rng(seed)
    train, test = [], []
    for hist, grp in lab.groupby("histology"):
        patients = grp["patient"].unique()
        rng.shuffle(patients)
        n_test = int(round(len(patients) * test_size))
        test_pat = set(patients[:n_test])
        for s, row in grp.iterrows():
            (test if row["patient"] in test_pat else train).append(s)
    return sorted(train), sorted(test)


def run_dge(counts, logcpm, y):
    """edgeR-style DGE approximation: Welch t-test on log-CPM + BH FDR + logFC.

    logFC defined as mean(LUSC) - mean(LUAD) in log2 space (paper convention).
    Returns DataFrame indexed by gene with columns logFC, pval, fdr.
    """
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    luad = logcpm.loc[:, y[y == "LUAD"].index]
    lusc = logcpm.loc[:, y[y == "LUSC"].index]
    logfc = lusc.mean(axis=1) - luad.mean(axis=1)
    t, p = stats.ttest_ind(lusc, luad, axis=1, equal_var=False)
    fdr = multipletests(p, method="fdr_bh")[1]
    return pd.DataFrame({"logFC": logfc.values, "pval": p, "fdr": fdr},
                        index=logcpm.index)


def select_dge(dge, k=C.N_DGE):
    sig = dge[(dge.fdr < C.DGE_FDR) & (dge.logFC.abs() > C.DGE_LOGFC)]
    return list(sig.sort_values("pval").head(k).index)


def select_pca(X_train, k=C.N_PCA):
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    Z = StandardScaler().fit_transform(X_train.values)
    pca = PCA(n_components=10, random_state=C.SEED).fit(Z)
    load = pd.Series(np.abs(pca.components_[0]), index=X_train.columns)
    return list(load.sort_values(ascending=False).head(k).index)


def select_mrmr(X_train, y_train, k=C.N_MRMR):
    """Fast vectorized mRMR (MID scheme) over the FULL gene set.

    relevance = ANOVA F-stat; redundancy = mean |Pearson r| to already-selected
    genes. Greedy: pick argmax(relevance - mean_redundancy). Implemented with a
    standardized matrix so each step is a single matrix-vector product, giving genuine
    minimum-redundancy selection at ~O(K * n_genes * n_samples)."""
    import numpy as np
    from sklearn.feature_selection import f_classif
    yc = (y_train == C.POS_LABEL).astype(int).values
    genes = list(X_train.columns)
    F, _ = f_classif(X_train.values, yc)
    rel = np.nan_to_num(F).astype(np.float64)
    rel = (rel - rel.min()) / (rel.max() - rel.min() + 1e-12)  # scale to [0,1]
    # so the redundancy penalty (mean |corr|, also in [0,1]) actually competes
    # standardized gene matrix (genes x samples) for fast correlation
    Z = X_train.values.T.astype(np.float64)          # (G, n)
    Z = Z - Z.mean(axis=1, keepdims=True)
    sd = Z.std(axis=1, keepdims=True); sd[sd == 0] = 1.0
    Z = Z / sd
    n = Z.shape[1]
    G = Z.shape[0]
    selected = [int(np.argmax(rel))]
    red = np.zeros(G)
    chosen = np.zeros(G, dtype=bool); chosen[selected[0]] = True
    for _ in range(min(k, G) - 1):
        last = selected[-1]
        red += np.abs(Z @ Z[last] / n)               # |corr| to newest pick
        score = rel - red / len(selected)
        score[chosen] = -np.inf
        nxt = int(np.argmax(score))
        selected.append(nxt); chosen[nxt] = True
    return [genes[i] for i in selected]


def select_xgb(X_train, y_train):
    from xgboost import XGBClassifier
    yb = (y_train == C.POS_LABEL).astype(int)
    clf = XGBClassifier(n_estimators=400, max_depth=4, learning_rate=0.1,
                        subsample=0.8, colsample_bytree=0.5, n_jobs=4,
                        eval_metric="logloss", random_state=C.SEED)
    clf.fit(X_train, yb)
    imp = pd.Series(clf.feature_importances_, index=X_train.columns)
    imp = imp[imp > 0].sort_values(ascending=False)
    # cap to the paper's XGBoost selection size (148); 400 shallow trees spread
    # negligible gain across a long tail, so the top slice holds the real signal.
    return list(imp.head(C.XGB_TOP_N).index), clf


def select_lasso(X_train, y_train):
    """L1-logistic feature selection using the 1-SE rule (glmnet-style): pick the
    most-regularized lambda whose CV accuracy is within 1 SE of the best. This
    yields a sparse, stable gene set rather than the barely-penalized lambda.min."""
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_score
    from sklearn.preprocessing import StandardScaler
    yb = (y_train == C.POS_LABEL).astype(int).values
    Xs = StandardScaler().fit_transform(X_train.values)
    Cs = np.logspace(-3, 1, 20)
    skf = StratifiedKFold(C.CV_FOLDS, shuffle=True, random_state=C.SEED)
    means, ses = [], []
    for c in Cs:
        s = cross_val_score(LogisticRegression(penalty="l1", solver="liblinear",
                            C=c, max_iter=3000, random_state=C.SEED), Xs, yb, cv=skf, scoring="accuracy")
        means.append(s.mean()); ses.append(s.std() / np.sqrt(len(s)))
    means, ses = np.array(means), np.array(ses)
    best = means.argmax()
    c_1se = Cs[means >= means[best] - ses[best]].min()   # strongest penalty within 1SE
    clf = LogisticRegression(penalty="l1", solver="liblinear", C=c_1se,
                             max_iter=3000, random_state=C.SEED).fit(Xs, yb)
    coef = pd.Series(clf.coef_.ravel(), index=X_train.columns)
    return list(coef[coef != 0].index)


CACHE = C.PROC / "fs_cache"


def _load_common():
    logcpm = pd.read_parquet(C.PROC / "expr_logcpm.parquet")
    counts = pd.read_parquet(C.PROC / "counts_filtered.parquet")
    labels = pd.read_csv(C.PROC / "labels.csv")
    lab_t = labels[labels.sample_type == "tumor"].set_index("sample")
    samples = [s for s in logcpm.columns if s in lab_t.index]
    y = lab_t.loc[samples, "histology"]
    if (C.PROC / "split.json").exists():
        sp = json.load(open(C.PROC / "split.json"))
        train, test = sp["train"], sp["test"]
    else:
        train, test = patient_split(labels, samples)
        json.dump({"train": train, "test": test}, open(C.PROC / "split.json", "w"))
    return logcpm, counts, y, train, test


def run_method(name):
    """Run one selector and cache its gene list (resumable across calls)."""
    import time, os
    CACHE.mkdir(exist_ok=True)
    out = CACHE / f"{name}.json"
    if out.exists() and os.environ.get("FS_FORCE") != "1":
        print(f"  {name}: cached ({len(json.load(open(out)))} genes)"); return
    logcpm, counts, y, train, test = _load_common()
    X = logcpm.T
    X_train, y_train = X.loc[train], y.loc[train]
    t0 = time.time()
    if name == "DGE":
        dge = run_dge(counts, logcpm, y_train); dge.to_csv(C.PROC / "dge_stats.csv")
        genes = select_dge(dge)
    elif name == "PCA":
        genes = select_pca(X_train)
    elif name == "mRMR":
        genes = select_mrmr(X_train, y_train)
    elif name == "XGBoost":
        genes, _ = select_xgb(X_train, y_train)
    elif name == "Lasso":
        genes = select_lasso(X_train, y_train)
    else:
        raise ValueError(name)
    json.dump(list(genes), open(out, "w"))
    print(f"  {name}: {len(genes)} genes  ({time.time()-t0:.0f}s)")


def combine():
    """Build the overlap sets once all five method caches exist."""
    sets = {m: json.load(open(CACHE / f"{m}.json"))
            for m in ["DGE", "PCA", "mRMR", "XGBoost", "Lasso"]}

    # membership matrix
    all_genes = sorted(set().union(*sets.values()))
    M = pd.DataFrame(0, index=all_genes, columns=list(sets))
    for m, gs in sets.items():
        M.loc[[g for g in gs if g in M.index], m] = 1
    M["n_methods"] = M.sum(axis=1)
    M.sort_values("n_methods", ascending=False).to_csv(
        C.TAB / "overlap_membership.csv")

    overlap131 = sorted(M.index[M.n_methods >= C.OVERLAP_CLASSIFY])
    overlap17 = sorted(M.index[M.n_methods >= C.OVERLAP_BIOMARKER])
    print(f"\nOverlap >=2 methods: {len(overlap131)} genes (paper 131)")
    print(f"Overlap >=3 methods: {len(overlap17)} genes (paper 17)")

    # how many of the paper's 17 did we recover?
    rec = sorted(set(overlap17) & set(C.PAPER_17_BIOMARKERS))
    print(f"Recovered {len(rec)}/{len(C.PAPER_17_BIOMARKERS)} paper biomarkers: {rec}")

    out = {**sets, "overlap131": overlap131, "overlap17": overlap17}
    json.dump(out, open(C.TAB / "feature_sets.json", "w"), indent=1)
    print("Wrote feature_sets.json, overlap_membership.csv, dge_stats.csv")


if __name__ == "__main__":
    import sys
    args = sys.argv[1:] or ["DGE", "PCA", "XGBoost", "Lasso", "mRMR", "combine"]
    for a in args:
        combine() if a == "combine" else run_method(a)
