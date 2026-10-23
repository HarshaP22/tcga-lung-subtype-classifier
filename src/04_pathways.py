"""Step 4: pathway enrichment (GO / KEGG / Reactome) on the >=2-method overlap.

Splits the >=2-method overlap into up/down in LUSC vs LUAD (paper: 74 up, 57
down) and runs offline over-representation analysis against local Enrichr .gmt
libraries, with the ~14k analysed genes as background (Fisher's exact, like the
paper's DAVID step).

Requires .gmt libraries in data/raw/genesets/ (see docs/DATA.md §D). If
absent, the step is skipped with a message.

Outputs
-------
results/tables/enrichment_up.csv, enrichment_down.csv
results/figures/pathways_up.png, pathways_down.png
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

GENESET_DIR = C.RAW / "genesets"


def read_gmt(path):
    """Parse a .gmt into {term: [genes]} (col0=term, col1=desc, col2+=genes)."""
    d = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                d[parts[0]] = [g for g in parts[2:] if g]
    return d


def run_enrich(genes, background, gmts, tag):
    import gseapy as gp
    frames = []
    for gmt in gmts:
        try:
            enr = gp.enrich(gene_list=list(genes), gene_sets=read_gmt(gmt),
                            background=list(background), outdir=None)
            res = enr.results.copy()
            res["library"] = gmt.stem
            frames.append(res)
        except Exception as e:
            print(f"  [{tag}] {gmt.name}: {e}")
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    pcol = "Adjusted P-value" if "Adjusted P-value" in out else "P-value"
    return out.sort_values(pcol)


def barplot(df, path, title, top=10):
    if df.empty:
        return
    pcol = "Adjusted P-value" if "Adjusted P-value" in df else "P-value"
    d = df.head(top).iloc[::-1]
    plt.figure(figsize=(8, 4.5))
    plt.barh(d["Term"].str.slice(0, 45), -np.log10(d[pcol].astype(float) + 1e-300))
    plt.xlabel("-log10 adj.P"); plt.title(title); plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight"); plt.close()


def main():
    logcpm = pd.read_parquet(C.PROC / "expr_logcpm.parquet")
    dge = pd.read_csv(C.PROC / "dge_stats.csv", index_col=0)
    sets = json.load(open(C.TAB / "feature_sets.json"))
    genes131 = [g for g in sets["overlap131"] if g in dge.index]
    background = list(logcpm.index)

    up = [g for g in genes131 if dge.loc[g, "logFC"] > 0]     # up in LUSC
    dn = [g for g in genes131 if dge.loc[g, "logFC"] <= 0]    # up in LUAD
    print(f"131 overlap: {len(up)} up in LUSC / {len(dn)} up in LUAD "
          f"(paper 74 / 57)")

    if not GENESET_DIR.exists() or not list(GENESET_DIR.glob("*.gmt")):
        print(f"No .gmt libraries in {GENESET_DIR}: pathway step skipped "
              "(see docs/DATA.md §D).")
        # still record the gene split for the report
        pd.Series(up, name="gene").to_csv(C.TAB / "genes_up_in_LUSC.csv", index=False)
        pd.Series(dn, name="gene").to_csv(C.TAB / "genes_up_in_LUAD.csv", index=False)
        return
    gmts = sorted(GENESET_DIR.glob("*.gmt"))
    up_res = run_enrich(up, background, gmts, "up")
    dn_res = run_enrich(dn, background, gmts, "down")
    up_res.to_csv(C.TAB / "enrichment_up.csv", index=False)
    dn_res.to_csv(C.TAB / "enrichment_down.csv", index=False)
    barplot(up_res, C.FIG / "pathways_up.png", "Pathways enriched in genes up in LUSC")
    barplot(dn_res, C.FIG / "pathways_down.png", "Pathways enriched in genes up in LUAD")
    if not up_res.empty:
        print("Top up-in-LUSC terms:\n", up_res.head(8)[["Term", "library"]].to_string())
    if not dn_res.empty:
        print("Top up-in-LUAD terms:\n", dn_res.head(8)[["Term", "library"]].to_string())


if __name__ == "__main__":
    main()
