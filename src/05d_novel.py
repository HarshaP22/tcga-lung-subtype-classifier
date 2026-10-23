"""Extension D: novel exploration ("find something new").

1. Stability selection: bootstrap the training set B times, re-run Lasso +
   XGBoost, and record how often each gene is picked. Stable-but-unreported
   genes are candidate biomarkers the single-run paper could have missed.
2. Under-reported candidates: genes in the 131-overlap with high single-gene
   AUC that are NOT among the paper's highlighted genes.
3. Confound audit: association of each biomarker with sex / tumour stage
   (needs phenotype); flags markers whose signal may be confounded.
4. Unsupervised substructure: cluster tumours within LUAD and within LUSC on
   the 131 panel; report silhouette and whether clean subclusters emerge.

Outputs
-------
results/tables/stability_selection.csv, novel_candidates.csv,
                confound_audit.csv, subclusters.csv
results/figures/stability_hist.png, subcluster_<hist>.png
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parent))
import config as C


def _find(patterns, root):
    for p in patterns:
        h = sorted(root.glob(p))
        if h:
            return h[0]
    return None


def stability(X, y, B=30, pool=2000):
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.feature_selection import f_classif
    from xgboost import XGBClassifier
    # restrict to top-`pool` genes by relevance for tractable resampling
    F = np.nan_to_num(f_classif(X.values, y.values)[0])
    top = pd.Series(F, index=X.columns).sort_values(ascending=False).head(pool).index
    X = X[top]
    rng = np.random.default_rng(C.SEED)
    counts = pd.Series(0, index=X.columns, dtype=float)
    n = len(X)
    for b in range(B):
        idx = rng.integers(0, n, n)
        Xb, yb = X.iloc[idx], y.iloc[idx]
        sc = StandardScaler().fit(Xb.values)
        l1 = LogisticRegression(penalty="l1", solver="liblinear", C=0.1, random_state=b,
                                max_iter=1000).fit(sc.transform(Xb.values), yb)
        sel = set(X.columns[(l1.coef_.ravel() != 0)])
        xgb = XGBClassifier(n_estimators=150, max_depth=3, learning_rate=0.1,
                            subsample=0.8, colsample_bytree=0.3, n_jobs=4,
                            eval_metric="logloss", random_state=b).fit(Xb, yb)
        sel |= set(X.columns[xgb.feature_importances_ > 0])
        counts[list(sel)] += 1
    return (counts / B).sort_values(ascending=False)


def main():
    logcpm = pd.read_parquet(C.PROC / "expr_logcpm.parquet")
    labels = pd.read_csv(C.PROC / "labels.csv").set_index("sample")
    sets = json.load(open(C.TAB / "feature_sets.json"))
    split = json.load(open(C.PROC / "split.json"))
    from sklearn.metrics import roc_auc_score
    X = logcpm.T
    y = labels.loc[X.index, "histology"]
    yb = (y == C.POS_LABEL).astype(int)
    genes131 = [g for g in sets["overlap131"] if g in X.columns]

    # 1. stability selection over full gene space
    print("Stability selection …")
    freq = stability(X.loc[split["train"]], yb.loc[split["train"]], B=30)
    freq.to_csv(C.TAB / "stability_selection.csv", header=["selection_freq"])
    plt.figure(figsize=(6, 4)); freq.head(40).iloc[::-1].plot.barh(fontsize=6)
    plt.xlabel("selection frequency"); plt.title("Top stable genes (Lasso+XGB, 50 boots)")
    plt.tight_layout(); plt.savefig(C.FIG / "stability_hist.png", dpi=150, bbox_inches="tight"); plt.close()

    # 2. under-reported candidates: high AUC, in 131, not paper-highlighted
    reported = set(C.PAPER_17_BIOMARKERS) | {"TP63", "NAPSA", "MLPH", "DSC3",
                 "NKX2-1", "SOX2", "MUC1", "KRT5", "KRT6A", "KRT14", "KRT17", "PVRL1"}
    pos = (y == C.POS_LABEL)
    auc = {}
    for g in genes131:
        a = roc_auc_score(pos.astype(int), logcpm.loc[g]); auc[g] = max(a, 1 - a)
    aucs = pd.Series(auc).sort_values(ascending=False)
    novel = aucs[[g for g in aucs.index if g not in reported]].head(25)
    novel_df = pd.DataFrame({"auc": novel, "stability": freq.reindex(novel.index).fillna(0)})
    novel_df.to_csv(C.TAB / "novel_candidates.csv")
    print("Top under-reported high-AUC candidates:\n",
          novel_df.head(12).round(3).to_string())

    # 3. confound audit (needs phenotype)
    ph = _find(["*clinical*", "*phenotype*"], C.RAW)
    if ph is not None:
        rows = []
        pfiles = list(C.RAW.glob("*clinical*")) + list(C.RAW.glob("*phenotype*"))
        phen = pd.concat([pd.read_csv(f, sep="\t", low_memory=False) for f in pfiles],
                         ignore_index=True)
        idc = phen.columns[0]
        phen = phen.set_index(idc)
        sexcol = next((c for c in phen.columns if "gender" in c.lower()), None)
        from scipy import stats
        hist = labels["histology"]
        for g in sets["overlap17"]:
            if g not in logcpm.index:
                continue
            expr = logcpm.loc[g]
            rec = {"gene": g}
            if sexcol:
                sx = phen[sexcol].reindex(expr.index)
                # test sex association WITHIN each histology (hold histology fixed),
                # otherwise sex just proxies the LUAD/LUSC male/female imbalance
                for h in ["LUAD", "LUSC"]:
                    m = (hist.reindex(expr.index) == h)
                    e, s = expr[m], sx[m]
                    grp = [e[s == v] for v in s.dropna().unique()]
                    if len(grp) == 2 and min(len(grp[0]), len(grp[1])) >= 5:
                        rec[f"sex_p_{h}"] = stats.mannwhitneyu(grp[0], grp[1]).pvalue
            rows.append(rec)
        ca = pd.DataFrame(rows)
        ca.to_csv(C.TAB / "confound_audit.csv", index=False)
        from statsmodels.stats.multitest import multipletests
        for h in ["LUAD", "LUSC"]:
            col = f"sex_p_{h}"
            if col in ca:
                sig = (multipletests(ca[col].dropna(), method="fdr_bh")[1] < 0.05).sum()
                print(f"  sex-confounded biomarkers within {h} (FDR<0.05): "
                      f"{sig}/{ca[col].notna().sum()}")
        print("Confound audit written (sex association of 17 biomarkers).")
    else:
        print("No phenotype file: confound audit skipped.")

    # 4. unsupervised substructure within each histology
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.metrics import silhouette_score
    from sklearn.preprocessing import StandardScaler
    sub_rows = []
    for h in ["LUAD", "LUSC"]:
        samp = y[y == h].index
        Z = StandardScaler().fit_transform(X.loc[samp, genes131].values)
        for k in (2, 3, 4):
            cl = AgglomerativeClustering(n_clusters=k).fit_predict(Z)
            sub_rows.append({"histology": h, "k": k, "silhouette": silhouette_score(Z, cl)})
    pd.DataFrame(sub_rows).to_csv(C.TAB / "subclusters.csv", index=False)
    print("Subcluster silhouettes:\n", pd.DataFrame(sub_rows).to_string(index=False))


if __name__ == "__main__":
    main()
