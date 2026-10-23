"""Extension B: SHAP interpretability + feature-count ablation.

(1) SHAP on an XGBoost LUAD-vs-LUSC model built from the overlap-131 panel;
    reports the genes that most drive predictions and whether they match known
    LUAD (NKX2-1/NAPSA/MUC1) vs LUSC (KRT5, TP63, SOX2, KRT6A) markers.
(2) Ablation: AUROC as a function of the number of genes (ranked by XGBoost
    gain): how few genes still reach >=0.95 AUROC?

Outputs
-------
results/figures/shap_summary.png, shap_bar.png
results/figures/auc_vs_ngenes.png
results/tables/ablation_auc.csv, shap_top_genes.csv
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


def main():
    logcpm = pd.read_parquet(C.PROC / "expr_logcpm.parquet")
    labels = pd.read_csv(C.PROC / "labels.csv").set_index("sample")
    sets = json.load(open(C.TAB / "feature_sets.json"))
    split = json.load(open(C.PROC / "split.json"))
    train, test = split["train"], split["test"]
    X = logcpm.T
    y = (labels.loc[X.index, "histology"] == C.POS_LABEL).astype(int)

    import shap
    from xgboost import XGBClassifier
    from sklearn.metrics import roc_auc_score

    genes131 = [g for g in sets["overlap131"] if g in X.columns]
    clf = XGBClassifier(n_estimators=400, max_depth=4, learning_rate=0.08,
                        subsample=0.8, colsample_bytree=0.7, n_jobs=4,
                        eval_metric="logloss", random_state=C.SEED)
    clf.fit(X.loc[train, genes131], y.loc[train])

    # exact tree SHAP from XGBoost natively (robust across xgboost/shap versions)
    import xgboost as xgb
    dtest = xgb.DMatrix(X.loc[test, genes131])
    contribs = clf.get_booster().predict(dtest, pred_contribs=True)
    sv = contribs[:, :-1]     # drop bias column
    shap.summary_plot(sv, X.loc[test, genes131], show=False, max_display=20)
    plt.title("SHAP: LUAD(+) vs LUSC(-), overlap-131"); plt.tight_layout()
    plt.savefig(C.FIG / "shap_summary.png", dpi=150, bbox_inches="tight"); plt.close()

    mean_abs = pd.Series(np.abs(sv).mean(axis=0), index=genes131).sort_values(ascending=False)
    mean_abs.head(25).to_csv(C.TAB / "shap_top_genes.csv")
    plt.figure(figsize=(6, 7))
    mean_abs.head(20).iloc[::-1].plot.barh()
    plt.xlabel("mean |SHAP|"); plt.title("Top SHAP genes"); plt.tight_layout()
    plt.savefig(C.FIG / "shap_bar.png", dpi=150, bbox_inches="tight"); plt.close()

    known_luad = {"NKX2-1", "NAPSA", "MUC1", "MLPH", "TTF1"}
    known_lusc = {"KRT5", "KRT6A", "KRT14", "KRT17", "TP63", "SOX2", "TP53", "DSC3"}
    top = set(mean_abs.head(20).index)
    print("Top-20 SHAP genes overlapping known LUSC markers:", sorted(top & known_lusc))
    print("Top-20 SHAP genes overlapping known LUAD markers:", sorted(top & known_luad))

    # ---- ablation: AUROC vs #genes (rank by XGB gain on full training) -----
    full = XGBClassifier(n_estimators=400, max_depth=4, learning_rate=0.08,
                         subsample=0.8, colsample_bytree=0.3, n_jobs=4,
                         eval_metric="logloss", random_state=C.SEED)
    full.fit(X.loc[train], y.loc[train])
    ranked = pd.Series(full.feature_importances_, index=X.columns).sort_values(ascending=False)
    ks = [2, 5, 10, 17, 25, 50, 100, 131, 250, 500]
    rows = []
    from sklearn.ensemble import RandomForestClassifier
    for k in ks:
        gg = list(ranked.head(k).index)
        rf = RandomForestClassifier(n_estimators=400, n_jobs=4, random_state=C.SEED)
        rf.fit(X.loc[train, gg], y.loc[train])
        p = rf.predict_proba(X.loc[test, gg])[:, 1]
        rows.append({"n_genes": k, "auroc": roc_auc_score(y.loc[test], p)})
    ab = pd.DataFrame(rows); ab.to_csv(C.TAB / "ablation_auc.csv", index=False)
    plt.figure(figsize=(6.5, 4.5))
    plt.plot(ab.n_genes, ab.auroc, "o-")
    plt.axhline(0.95, color="r", ls="--", lw=.8, label="0.95")
    plt.xscale("log"); plt.xlabel("# genes (log)"); plt.ylabel("test AUROC")
    plt.title("Feature-count ablation"); plt.legend(); plt.tight_layout()
    plt.savefig(C.FIG / "auc_vs_ngenes.png", dpi=150, bbox_inches="tight"); plt.close()
    print(ab.to_string(index=False))


if __name__ == "__main__":
    main()
