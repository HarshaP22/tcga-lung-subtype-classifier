"""Step 7: one-page summary figure for the README.

Six panels: (A) classifier performance per gene set, (B) biomarker AUC vs the
paper, (C) TCGA vs external (GSE28582) AUC, (D) feature-count ablation,
(E) overlap sensitivity to selector depth, (F) top pathways up in LUSC.

Output: results/figures/summary_figure.png
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parent))
import config as C

PAPER_AUC = {"KRT5": 0.973, "ARHGEF38": 0.957, "KRT6A": 0.953, "PERP": 0.953,
             "NECTIN1": 0.951}
ALIAS = {"NECTIN1": "PVRL1"}  # TCGA/gencode v22 still uses the old symbol


def main():
    T = C.TAB
    m = pd.read_csv(T / "classification_metrics.csv")
    t = pd.read_csv(T / "biomarker_auc_tcga.csv", index_col=0)
    e = pd.read_csv(T / "biomarker_auc_external.csv", index_col=0)
    ab = pd.read_csv(T / "ablation_auc.csv")
    sv = pd.read_csv(T / "sensitivity_topN.csv")
    up = pd.read_csv(T / "enrichment_up.csv")

    fig, ax = plt.subplots(2, 3, figsize=(17, 9.5))
    fig.suptitle("Overlapping feature selection for LUAD vs LUSC: reproduction of "
                 "Chen & Dhahbi (Sci Rep 2021) on TCGA (n = 1,029 tumours)",
                 fontsize=13, fontweight="bold")

    # A: classifier performance
    a = ax[0, 0]
    lbl = m.method.replace({"overlap131": "Overlap ≥2", "overlap17": "Overlap ≥3"})
    x = np.arange(len(m))
    col = ["#bbbbbb"] * 5 + ["#2c7fb8", "#7fcdbb"]
    a.bar(x, m.accuracy, color=col,
          yerr=[m.accuracy - m.acc_lo, m.acc_hi - m.accuracy], capsize=3)
    a.plot(x, m.auroc, "kD", ms=5, label="AUROC")
    for i, n in enumerate(m.n_genes):
        a.text(i, 0.905, f"n={n}", ha="center", fontsize=7, rotation=90, color="white")
    a.set_xticks(x); a.set_xticklabels(lbl, rotation=35, ha="right", fontsize=8)
    a.set_ylim(0.9, 1.0); a.set_ylabel("test accuracy (95% bootstrap CI)")
    a.set_title("A. Random-forest validation per gene set"); a.legend(loc="upper left", fontsize=8)

    # B: biomarker AUC vs paper
    b = ax[0, 1]
    genes = list(PAPER_AUC)
    ours = [t.loc[ALIAS.get(g, g), "auc"] for g in genes]
    xx = np.arange(len(genes))
    b.bar(xx - 0.2, list(PAPER_AUC.values()), 0.4, label="paper", color="#fdae61")
    b.bar(xx + 0.2, ours, 0.4, label="this work", color="#2c7fb8")
    for i, v in enumerate(ours):
        b.text(i + 0.2, v + 0.002, f"{v:.3f}", ha="center", fontsize=7)
    b.set_xticks(xx); b.set_xticklabels(genes, fontsize=9)
    b.set_ylim(0.9, 0.99); b.set_ylabel("single-gene AUC (LUAD vs LUSC)")
    b.set_title("B. Biomarker AUCs reproduce the paper"); b.legend(fontsize=8)

    # C: external validation
    c = ax[0, 2]
    comp = t[["auc"]].join(e, how="inner")
    c.scatter(comp.auc, comp.auc_external, s=12, alpha=.6, color="#2c7fb8")
    offs = {"KRT5": (-28, -14), "KRT6A": (-40, 6), "DSC3": (6, -12),
            "DSG3": (6, 4), "TRIM29": (-48, -4), "PVRL1": (6, -14)}
    for g, o in offs.items():
        if g in comp.index:
            xy = (comp.loc[g, "auc"], comp.loc[g, "auc_external"])
            c.scatter(*xy, s=22, color="#de2d26", zorder=3)
            c.annotate(g, xy, xytext=o, textcoords="offset points", fontsize=7,
                       arrowprops=dict(arrowstyle="-", lw=.5, color="grey"))
    c.plot([0.8, 1], [0.8, 1], "k--", lw=.7)
    c.set_xlabel("AUC (TCGA RNA-seq)"); c.set_ylabel("AUC (GSE28582 microarray)")
    c.set_title(f"C. External cohort GSE28582\n(mean AUC {comp.auc_external.mean():.3f} "
                f"over {len(comp)} genes)")

    # D: ablation
    d = ax[1, 0]
    d.plot(ab.n_genes, ab.auroc, "o-", color="#2c7fb8")
    d.axhline(0.95, ls="--", c="r", lw=.8)
    d.set_xscale("log"); d.set_ylim(0.94, 1.0)
    d.set_xlabel("# genes (ranked by XGBoost gain, log scale)"); d.set_ylabel("test AUROC")
    d.set_title("D. A handful of genes is enough")

    # E: sensitivity
    s = ax[1, 1]
    s.plot(sv.top_N, sv.overlap_ge2, "o-", label="≥2 methods")
    s.plot(sv.top_N, sv.overlap_ge3, "s-", label="≥3 methods")
    s.axhline(131, ls="--", c="C0", lw=.8); s.axhline(17, ls="--", c="C1", lw=.8)
    s.text(330, 138, "paper: 131", color="C0", fontsize=8)
    s.text(330, 24, "paper: 17", color="C1", fontsize=8)
    s2 = s.twinx(); s2.plot(sv.top_N, sv.acc_ge2, "k:", lw=1); s2.set_ylim(0.9, 1.0)
    s2.set_ylabel("≥2-panel accuracy (dotted)", fontsize=8)
    s.set_xlabel("top-N kept from DGE / PCA / mRMR"); s.set_ylabel("# overlap genes")
    s.set_title("E. Overlap size is threshold-driven; accuracy isn't"); s.legend(fontsize=8, loc="upper left")

    # F: pathways
    f = ax[1, 2]
    pc = "Adjusted P-value"
    top = up.sort_values(pc).head(8).iloc[::-1]
    names = top.Term.str.replace(r"\s*\(GO:\d+\)|\s*R-HSA-\d+", "", regex=True).str.slice(0, 40)
    f.barh(names, -np.log10(top[pc].astype(float)), color="#de2d26")
    f.set_xlabel("-log10 adjusted p"); f.tick_params(axis="y", labelsize=8)
    f.set_title("F. Pathways up in LUSC (keratinization)")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = C.FIG / "summary_figure.png"
    plt.savefig(out, dpi=150); plt.close()
    print("Wrote", out)


if __name__ == "__main__":
    main()
