"""Step 3: biomarker ROC/AUC (TCGA) + external validation on GSE28582.

For the >=3-method overlap genes (the paper's 17 biomarkers), compute the
single-gene AUC for discriminating LUAD vs LUSC in TCGA, then re-check the same
genes in the independent GSE28582 microarray cohort (50 LUAD + 28 LUSC).

Outputs
-------
results/tables/biomarker_auc_tcga.csv
results/tables/biomarker_auc_external.csv   (if GSE28582 present)
results/figures/roc_biomarkers_up.png, roc_biomarkers_down.png
results/figures/auc_tcga_vs_external.png
"""
from __future__ import annotations
import sys, json, gzip, re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parent))
import config as C
import utils as U
from sklearn.metrics import roc_auc_score, roc_curve


def single_gene_auc(expr_row, y_pos_mask):
    """AUC using one gene's expression as the score (orient to max AUC)."""
    a = roc_auc_score(y_pos_mask.astype(int), expr_row.values)
    return max(a, 1 - a)


def tcga_auc(logcpm, y, genes):
    pos = (y == C.POS_LABEL)          # LUAD positive
    rows = []
    for g in genes:
        if g in logcpm.index:
            rows.append({"gene": g, "auc": single_gene_auc(logcpm.loc[g], pos),
                         "logFC_LUSCvsLUAD": logcpm.loc[g, y[y=="LUSC"].index].mean()
                                             - logcpm.loc[g, y[y=="LUAD"].index].mean()})
    return pd.DataFrame(rows).set_index("gene")


# --------------------------- GSE28582 parsing ------------------------------
def _find(patterns, root):
    for pat in patterns:
        hits = sorted(root.glob(pat))
        if hits:
            return hits[0]
    return None


def load_gse28582():
    f = _find(["*GSE28582*GPL570*matrix*.txt.gz", "*GSE28582*matrix*.txt*",
               "*GSE28582*.txt.gz"], C.EXTERNAL)
    if f is None:
        return None
    op = gzip.open if str(f).endswith(".gz") else open
    meta, data_lines, in_tbl = {}, [], False
    with op(f, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_tbl = True; continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_tbl:
                data_lines.append(line)
            elif line.startswith("!Sample_"):
                key = line.split("\t", 1)[0]
                vals = [v.strip().strip('"') for v in line.rstrip("\n").split("\t")[1:]]
                meta.setdefault(key, []).append(vals)
    from io import StringIO
    expr = pd.read_csv(StringIO("".join(data_lines)), sep="\t", index_col=0)
    expr.index = expr.index.astype(str).str.strip('"')
    # derive histology label from any characteristics/source/title line
    labels = _label_from_meta(meta, expr.columns)
    return expr, labels


def _label_from_meta(meta, cols):
    text = {c: "" for c in cols}
    for key, blocks in meta.items():
        if any(k in key.lower() for k in ["characteristics", "source", "title", "description"]):
            for vals in blocks:
                for c, v in zip(cols, vals):
                    text[c] += " " + v.lower()
    lab = {}
    for c, t in text.items():
        if "squamous" in t: lab[c] = "LUSC"
        elif "adeno" in t:  lab[c] = "LUAD"
        else:               lab[c] = "other"
    return pd.Series(lab)


def map_probes_to_symbols(expr):
    """Collapse GPL570 probes to gene symbols using GPL570 annotation if present."""
    ann = _find(["GPL570*annot*.gz", "GPL570*.txt*", "*GPL570*annot*"], C.RAW)
    if ann is None:
        print("  WARNING: GPL570 annotation not found: external step will map by "
              "probe only (skipped).")
        return None
    op = gzip.open if str(ann).endswith(".gz") else open
    rows, header, started = [], None, False
    with op(ann, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("^") or line.startswith("!") or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts; started = True; continue
            if started:
                rows.append(parts)
    ann_df = pd.DataFrame(rows, columns=header)
    idcol = ann_df.columns[0]
    symcol = next((c for c in ann_df.columns if "gene symbol" in c.lower()
                   or c.lower() == "gene_symbol"), None)
    if symcol is None:
        return None
    m = ann_df.set_index(idcol)[symcol]
    expr = expr.copy()
    expr["_sym"] = m.reindex(expr.index.astype(str)).values
    expr = expr.dropna(subset=["_sym"])
    expr["_sym"] = expr["_sym"].str.split("///").str[0].str.strip()
    return expr.groupby("_sym").max()


def main():
    logcpm = pd.read_parquet(C.PROC / "expr_logcpm.parquet")
    labels = pd.read_csv(C.PROC / "labels.csv").set_index("sample")
    sets = json.load(open(C.TAB / "feature_sets.json"))
    y = labels.loc[logcpm.columns, "histology"]
    genes = sets["overlap17"]

    tcga = tcga_auc(logcpm, y, genes).sort_values("auc", ascending=False)
    tcga.to_csv(C.TAB / "biomarker_auc_tcga.csv")
    print("TCGA single-gene AUCs (17 biomarkers):")
    print(tcga.round(3).to_string())

    # ROC curves split by up/down in LUSC
    up = tcga[tcga.logFC_LUSCvsLUAD > 0].index
    dn = tcga[tcga.logFC_LUSCvsLUAD <= 0].index
    pos = (y == C.POS_LABEL)
    for grp, tag, ttl in [(up, "up", "Up in LUSC"), (dn, "down", "Up in LUAD")]:
        grp = grp[:12]          # top-12 by AUC per direction (readable ROC panel)
        plt.figure(figsize=(6, 6))
        for g in grp:
            fpr, tpr, _ = roc_curve(pos.astype(int), logcpm.loc[g])
            if roc_auc_score(pos.astype(int), logcpm.loc[g]) < 0.5:
                fpr, tpr, _ = roc_curve((~pos).astype(int), logcpm.loc[g])
            plt.plot(fpr, tpr, lw=1, label=f"{g} ({tcga.loc[g,'auc']:.3f})")
        plt.plot([0, 1], [0, 1], "k--", lw=.7)
        plt.xlabel("1 - specificity"); plt.ylabel("sensitivity")
        plt.title(f"Biomarker ROC: {ttl} (TCGA)"); plt.legend(fontsize=7)
        plt.savefig(C.FIG / f"roc_biomarkers_{tag}.png", dpi=150, bbox_inches="tight")
        plt.close()

    # external validation
    ext = load_gse28582()
    if ext is None:
        print("GSE28582 not found in data/external: external validation skipped.")
        return
    expr, elab = ext
    expr_sym = map_probes_to_symbols(expr)
    keep = elab.isin(["LUAD", "LUSC"])
    print(f"GSE28582: {(elab=='LUAD').sum()} LUAD, {(elab=='LUSC').sum()} LUSC "
          f"(paper 50/28)")
    if expr_sym is None:
        print("No probe->symbol mapping; external AUC skipped.")
        return
    ex = expr_sym.loc[:, keep[keep].index]
    ey = elab.loc[ex.columns]
    epos = (ey == C.POS_LABEL)
    rows = []
    for g in genes:
        if g in ex.index:
            rows.append({"gene": g, "auc_external": single_gene_auc(ex.loc[g].astype(float), epos)})
    edf = pd.DataFrame(rows).set_index("gene")
    edf.to_csv(C.TAB / "biomarker_auc_external.csv")
    comp = tcga[["auc"]].join(edf, how="inner")
    plt.figure(figsize=(6, 6))
    plt.scatter(comp["auc"], comp["auc_external"])
    for g, r in comp.iterrows():
        plt.annotate(g, (r["auc"], r["auc_external"]), fontsize=7)
    plt.plot([0.5, 1], [0.5, 1], "k--", lw=.7)
    plt.xlabel("AUC (TCGA)"); plt.ylabel("AUC (GSE28582)")
    plt.title("Biomarker AUC: TCGA vs external")
    plt.savefig(C.FIG / "auc_tcga_vs_external.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Wrote external AUCs + comparison figure.")


if __name__ == "__main__":
    main()
