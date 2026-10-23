"""Step 6: overlap sensitivity analysis.

The paper reports 131 genes in >=2 methods and 17 in >=3, but several selector
thresholds are unspecified. This step asks: how sensitive is the overlap to those
choices, and do the conclusions survive? It compares configurations:

  principled   : our default (PCA on z-scored genes, Lasso 1-SE rule, XGB top-148)
  paper_sizes  : every selector matched to the paper's set sizes
                 (DGE/PCA/mRMR 500, XGBoost 148, Lasso ~68 via lambda search)
  paper_sizes_pca_centered : as above, PCA on centred (not scaled) log-CPM

For each: overlap sizes, recovery of the paper's biomarkers, and a fixed random
forest's test accuracy/AUROC on the >=2 and >=3 panels.

Outputs
-------
results/tables/sensitivity_overlap.csv
results/tables/feature_sets_paper_sizes.json
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent))
import config as C
import importlib
fs = importlib.import_module("01_feature_selection")


def lasso_target(X_train, y_train, target=68):
    """L1-logistic at the lambda whose non-zero count is closest to `target`."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    yb = (y_train == C.POS_LABEL).astype(int).values
    Xs = StandardScaler().fit_transform(X_train.values)
    best = None
    for c in np.logspace(-2.5, 0, 26):
        clf = LogisticRegression(penalty="l1", solver="liblinear", C=c,
                                 max_iter=3000, random_state=C.SEED).fit(Xs, yb)
        nz = int((clf.coef_ != 0).sum())
        if best is None or abs(nz - target) < abs(best[1] - target):
            best = (c, nz, clf.coef_.ravel().copy())
        if nz > target * 1.5:
            break
    coef = pd.Series(best[2], index=X_train.columns)
    return list(coef[coef != 0].index)


def pca_centered(X_train, k=C.N_PCA):
    from sklearn.decomposition import PCA
    Z = X_train.values - X_train.values.mean(axis=0)
    pca = PCA(n_components=5, random_state=C.SEED).fit(Z)
    load = pd.Series(np.abs(pca.components_[0]), index=X_train.columns)
    return list(load.sort_values(ascending=False).head(k).index)


def overlap(sets):
    all_g = sorted(set().union(*map(set, sets.values())))
    n = pd.Series({g: sum(g in set(s) for s in sets.values()) for g in all_g})
    return sorted(n.index[n >= 2]), sorted(n.index[n >= 3]), sorted(n.index[n == 5])


def rf_score(X, y, train, test, genes):
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import roc_auc_score, accuracy_score
    g = [x for x in genes if x in X.columns]
    rf = RandomForestClassifier(n_estimators=500, max_features="sqrt",
                                n_jobs=4, random_state=C.SEED)
    rf.fit(X.loc[train, g], y.loc[train])
    p = rf.predict_proba(X.loc[test, g])[:, list(rf.classes_).index(C.POS_LABEL)]
    return (accuracy_score(y.loc[test], rf.predict(X.loc[test, g])),
            roc_auc_score((y.loc[test] == C.POS_LABEL).astype(int), p))


def main():
    logcpm, counts, y, train, test = fs._load_common()
    X = logcpm.T
    Xtr, ytr = X.loc[train], y.loc[train]
    cache = {m: json.load(open(fs.CACHE / f"{m}.json"))
             for m in ["DGE", "PCA", "mRMR", "XGBoost", "Lasso"]}

    print("Lasso matched to ~68 genes …", flush=True)
    lf = fs.CACHE / "Lasso68_seeded.json"
    if lf.exists():
        lasso68 = json.load(open(lf))
    else:
        lasso68 = lasso_target(Xtr, ytr, 68); json.dump(lasso68, open(lf, "w"))
    print(f"  -> {len(lasso68)} genes", flush=True)
    pca_c = pca_centered(Xtr)

    configs = {
        "principled": cache,
        "paper_sizes": {**cache, "Lasso": lasso68},
        "paper_sizes_pca_centered": {**cache, "Lasso": lasso68, "PCA": pca_c},
    }
    rows = []
    paper17 = set(C.PAPER_17_BIOMARKERS) | {"PVRL1"}
    for name, sets in configs.items():
        o2, o3, o5 = overlap(sets)
        a2, u2 = rf_score(X, y, train, test, o2)
        a3, u3 = rf_score(X, y, train, test, o3)
        rec = sorted((set(o3) & paper17))
        rows.append({"config": name,
                     **{f"n_{m}": len(v) for m, v in sets.items()},
                     "overlap_ge2": len(o2), "overlap_ge3": len(o3), "in_all5": len(o5),
                     "acc_ge2": a2, "auroc_ge2": u2, "acc_ge3": a3, "auroc_ge3": u3,
                     "paper_biomarkers_in_ge3": len(rec),
                     "recovered": ";".join(rec), "all5_genes": ";".join(o5)})
        print(f"{name:26s} >=2:{len(o2):4d} >=3:{len(o3):4d} all5:{len(o5):3d} | "
              f"acc>=2 {a2:.3f}  acc>=3 {a3:.3f} | paper markers in >=3: {len(rec)}",
              flush=True)
        if name == "paper_sizes_pca_centered":
            json.dump({**sets, "overlap_ge2": o2, "overlap_ge3": o3},
                      open(C.TAB / "feature_sets_paper_sizes.json", "w"), indent=1)
    pd.DataFrame(rows).to_csv(C.TAB / "sensitivity_overlap.csv", index=False)
    print("Wrote sensitivity_overlap.csv")

    # ---- B. overlap size vs top-N of the three ranked selectors -------------
    # DGE (by p), PCA (|PC1|) and mRMR (greedy order) are all rankings, so a
    # top-N prefix is well defined. XGBoost/Lasso held at paper-matched sizes.
    dge = pd.read_csv(C.PROC / "dge_stats.csv", index_col=0)
    sig = dge[(dge.fdr < C.DGE_FDR) & (dge.logFC.abs() > C.DGE_LOGFC)]
    dge_rank = list(sig.sort_values("pval").index)
    curve = []
    for N in [50, 75, 100, 125, 150, 200, 250, 300, 400, 500]:
        sets = {"DGE": dge_rank[:N], "PCA": cache["PCA"][:N], "mRMR": cache["mRMR"][:N],
                "XGBoost": cache["XGBoost"], "Lasso": lasso68}
        o2, o3, _ = overlap(sets)
        a2, u2 = rf_score(X, y, train, test, o2)
        curve.append({"top_N": N, "overlap_ge2": len(o2), "overlap_ge3": len(o3),
                      "acc_ge2": a2, "auroc_ge2": u2,
                      "paper_markers_in_ge3": len(set(o3) & paper17)})
        print(f"  top-N={N:3d}  >=2:{len(o2):4d}  >=3:{len(o3):4d}  acc>=2 {a2:.3f}",
              flush=True)
    cv = pd.DataFrame(curve); cv.to_csv(C.TAB / "sensitivity_topN.csv", index=False)

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(cv.top_N, cv.overlap_ge2, "o-", label=">=2 methods")
    ax[0].plot(cv.top_N, cv.overlap_ge3, "s-", label=">=3 methods")
    ax[0].axhline(131, ls="--", c="C0", lw=.8); ax[0].axhline(17, ls="--", c="C1", lw=.8)
    ax[0].text(55, 138, "paper 131", color="C0", fontsize=8)
    ax[0].text(55, 24, "paper 17", color="C1", fontsize=8)
    ax[0].set_xlabel("top-N kept from DGE / PCA / mRMR"); ax[0].set_ylabel("# overlap genes")
    ax[0].set_title("Overlap size depends on selector depth"); ax[0].legend()
    ax[1].plot(cv.top_N, cv.acc_ge2, "o-", label="accuracy")
    ax[1].plot(cv.top_N, cv.auroc_ge2, "s-", label="AUROC")
    ax[1].set_ylim(0.9, 1.0); ax[1].set_xlabel("top-N"); ax[1].set_title(">=2 panel performance is flat")
    ax[1].legend()
    plt.tight_layout(); plt.savefig(C.FIG / "sensitivity_overlap.png", dpi=150); plt.close()
    print("Wrote sensitivity_topN.csv + sensitivity_overlap.png")


if __name__ == "__main__":
    main()
