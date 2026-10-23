"""Extension A: three-class LUAD / LUSC / Normal (the source guide's Project 1).

Adds solid-tissue normals as a third class and evaluates multi-class
classification with the overlap-131 panel vs the full gene set.

Outputs
-------
results/tables/threeclass_metrics.csv
results/figures/confusion_threeclass.png
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


def three_class_labels(labels):
    def lab(row):
        if row.sample_type == "normal":
            return "Normal"
        if row.sample_type == "tumor":
            return row.histology
        return None
    return labels.assign(y3=labels.apply(lab, axis=1)).dropna(subset=["y3"])


def grouped_split(df, seed=C.SEED, test_size=C.TEST_SIZE):
    rng = np.random.default_rng(seed)
    train, test = [], []
    for _, grp in df.groupby("y3"):
        pats = grp["patient"].unique(); rng.shuffle(pats)
        tst = set(pats[:int(round(len(pats) * test_size))])
        for s, r in grp.set_index("sample").iterrows():
            (test if r["patient"] in tst else train).append(s)
    return train, test


def main():
    allx = pd.read_parquet(C.PROC / "expr_logcpm_all.parquet")
    labels = pd.read_csv(C.PROC / "labels.csv")
    sets = json.load(open(C.TAB / "feature_sets.json"))
    lab3 = three_class_labels(labels)
    lab3 = lab3[lab3["sample"].isin(allx.columns)]
    print("class counts:\n", lab3["y3"].value_counts().to_string())

    X = allx.T.loc[lab3["sample"]]
    y = lab3.set_index("sample")["y3"]
    train, test = grouped_split(lab3)
    train = [s for s in train if s in X.index]; test = [s for s in test if s in X.index]

    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import (f1_score, roc_auc_score, confusion_matrix,
                                 ConfusionMatrixDisplay, accuracy_score)
    rows = []
    for name, genes in [("overlap131", sets["overlap131"]), ("all_genes", list(X.columns))]:
        gg = [g for g in genes if g in X.columns]
        clf = RandomForestClassifier(n_estimators=500, n_jobs=4, random_state=C.SEED)
        clf.fit(X.loc[train, gg], y.loc[train])
        pred = clf.predict(X.loc[test, gg])
        proba = clf.predict_proba(X.loc[test, gg])
        macro_f1 = f1_score(y.loc[test], pred, average="macro")
        try:
            auroc = roc_auc_score(pd.get_dummies(y.loc[test]), proba, multi_class="ovr", average="macro")
        except Exception:
            auroc = np.nan
        acc = accuracy_score(y.loc[test], pred)
        rows.append({"features": name, "n_genes": len(gg), "accuracy": acc,
                     "macro_f1": macro_f1, "macro_auroc": auroc})
        print(f"{name:11s} n={len(gg):5d}  acc={acc:.3f} macroF1={macro_f1:.3f} AUROC={auroc:.3f}")
        if name == "overlap131":
            cm = confusion_matrix(y.loc[test], pred, labels=clf.classes_)
            ConfusionMatrixDisplay(cm, display_labels=clf.classes_).plot(cmap="Blues")
            plt.title("Three-class (overlap-131)")
            plt.savefig(C.FIG / "confusion_threeclass.png", dpi=150, bbox_inches="tight")
            plt.close()
    pd.DataFrame(rows).to_csv(C.TAB / "threeclass_metrics.csv", index=False)
    print("Wrote threeclass_metrics.csv + confusion_threeclass.png")


if __name__ == "__main__":
    main()
