"""Step 2: random-forest validation of every gene set (paper Table of stats).

For each selected gene set (5 methods + overlap131 + overlap17):
 - Random Forest (ranger analogue), 5-fold CV on train for n_estimators/depth,
 - evaluate on the locked test set, LUAD=positive / LUSC=negative,
 - accuracy (10k-bootstrap 95% CI), sensitivity, specificity, precision,
   F-measure, AUROC.
Also: confusion matrices + clustered heatmaps for the 131- and 17-gene panels.

Outputs
-------
results/tables/classification_metrics.csv
results/figures/heatmap_131.png, heatmap_17.png
results/figures/confusion_131.png, confusion_17.png
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(str(Path(__file__).resolve().parent))
import config as C
import utils as U


def fit_eval(X_train, y_train, X_test, y_test):
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import GridSearchCV
    from sklearn.metrics import roc_auc_score
    grid = {"n_estimators": [300, 500], "max_features": ["sqrt", 0.1],
            "min_samples_leaf": [1, 3]}
    gs = GridSearchCV(RandomForestClassifier(random_state=C.SEED, n_jobs=4),
                      grid, cv=C.CV_FOLDS, scoring="accuracy", n_jobs=4)
    gs.fit(X_train, y_train)
    clf = gs.best_estimator_
    y_pred = clf.predict(X_test)
    proba = clf.predict_proba(X_test)[:, list(clf.classes_).index(C.POS_LABEL)]
    stats = U.binary_stats(y_test, y_pred, C.POS_LABEL)
    stats["auroc"] = roc_auc_score((y_test == C.POS_LABEL).astype(int), proba)
    acc, lo, hi = U.bootstrap_ci(y_test.values, y_pred, C.N_BOOTSTRAP, C.SEED, C.POS_LABEL)
    stats.update(acc_boot=acc, acc_lo=lo, acc_hi=hi, cv_acc=gs.best_score_)
    return stats, clf, y_pred


def heatmap(logcpm, genes, y, path, title):
    sub = logcpm.loc[[g for g in genes if g in logcpm.index]]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1) + 1e-9, axis=0)
    order = y.sort_values().index
    z = z[order]
    col_colors = y.loc[order].map({"LUAD": "#2c7fb8", "LUSC": "#de2d26"})
    g = sns.clustermap(z, col_cluster=False, cmap="RdBu_r", vmin=-2.5, vmax=2.5,
                       col_colors=col_colors.values, xticklabels=False,
                       yticklabels=(len(genes) <= 20), figsize=(9, max(4, len(genes) * 0.12)))
    g.fig.suptitle(title, y=1.02)
    g.savefig(path, dpi=150, bbox_inches="tight")
    plt.close("all")


CACHE = C.PROC / "clf_cache"
ORDER = ["DGE", "PCA", "mRMR", "XGBoost", "Lasso", "overlap131", "overlap17"]


def run_one(name):
    """Fit+evaluate one gene set, cache metrics + test predictions (resumable)."""
    import os
    CACHE.mkdir(exist_ok=True)
    out = CACHE / f"{name}.json"
    if out.exists() and os.environ.get("CLF_FORCE") != "1":
        print(f"  {name}: cached"); return
    logcpm = pd.read_parquet(C.PROC / "expr_logcpm.parquet")
    labels = pd.read_csv(C.PROC / "labels.csv").set_index("sample")
    sets = json.load(open(C.TAB / "feature_sets.json"))
    split = json.load(open(C.PROC / "split.json"))
    X = logcpm.T
    y = labels.loc[X.index, "histology"]
    genes = [g for g in sets[name] if g in X.columns]
    st, _, y_pred = fit_eval(X.loc[split["train"], genes], y.loc[split["train"]],
                             X.loc[split["test"], genes], y.loc[split["test"]])
    st = {"method": name, "n_genes": len(genes), **st}
    st = {k: (int(v) if isinstance(v, np.integer) else
              float(v) if isinstance(v, np.floating) else v)
          for k, v in st.items()}
    st["pred"] = list(map(str, y_pred))
    json.dump(st, open(out, "w"))
    print(f"{name:11s} n={len(genes):4d}  acc={st['accuracy']:.3f} "
          f"[{st['acc_lo']:.3f},{st['acc_hi']:.3f}]  sens={st['sensitivity']:.3f} "
          f"spec={st['specificity']:.3f} F={st['f_measure']:.3f} AUROC={st['auroc']:.3f}")


def finalize():
    from sklearn.metrics import ConfusionMatrixDisplay, confusion_matrix
    logcpm = pd.read_parquet(C.PROC / "expr_logcpm.parquet")
    labels = pd.read_csv(C.PROC / "labels.csv").set_index("sample")
    sets = json.load(open(C.TAB / "feature_sets.json"))
    split = json.load(open(C.PROC / "split.json"))
    test = split["test"]; y_test = labels.loc[test, "histology"]
    rows, preds = [], {}
    for name in ORDER:
        st = json.load(open(CACHE / f"{name}.json"))
        preds[name] = pd.Series(st.pop("pred"), index=test)
        rows.append(st)
    pd.DataFrame(rows).to_csv(C.TAB / "classification_metrics.csv", index=False)
    labels_order = [C.POS_LABEL, C.NEG_LABEL]
    for tag, genes in [("131", sets["overlap131"]), ("17", sets["overlap17"])]:
        gg = genes if len(genes) <= 60 else genes[:60]
        heatmap(logcpm[test], gg, y_test, C.FIG / f"heatmap_{tag}.png",
                f"Overlap-{tag} panel ({len(genes)} genes): test samples")
        cm = confusion_matrix(y_test, preds[f"overlap{tag}"], labels=labels_order)
        ConfusionMatrixDisplay(cm, display_labels=labels_order).plot(cmap="Blues")
        plt.title(f"Overlap-{tag}: test confusion")
        plt.savefig(C.FIG / f"confusion_{tag}.png", dpi=150, bbox_inches="tight")
        plt.close()
    print("Wrote classification_metrics.csv + heatmaps + confusion matrices")


def main():
    for name in ORDER:
        run_one(name)
    finalize()


if __name__ == "__main__":
    import sys
    args = sys.argv[1:]
    if not args:
        main()
    else:
        for a in args:
            finalize() if a == "finalize" else run_one(a)
