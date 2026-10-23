"""Step 0: preprocessing (memory-lean).

Reads TCGA-LUAD / TCGA-LUSC Xena gene-expression matrices, inverts the
log2(count+1) transform back to raw counts, merges the two cohorts, applies the
paper's low-count filter (drop genes with CPM<1 in >600 tumor samples) at the
Ensembl level, maps the survivors to gene symbols, then TMM + log2-CPM
normalizes.

Order note: the CPM filter is applied before Ensembl->symbol collapse (keeps
peak memory well under the sandbox's limit). Effect on the gene list is limited
to the handful of Ensembl IDs sharing a symbol; documented deviation.

Outputs (data/processed/)
-------
expr_logcpm.parquet        genes(symbols) x TUMOR samples, TMM+log2-CPM  [main task]
expr_logcpm_all.parquet    genes(symbols) x ALL samples (incl. normal)   [3-class]
counts_filtered.parquet    raw filtered counts, tumor samples (genes x samples)
labels.csv                 sample, patient, histology, sample_type
"""
from __future__ import annotations
import sys, gc
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent))
import config as C
import utils as U


def _find(patterns):
    for pat in patterns:
        hits = sorted(C.RAW.glob(pat))
        if hits:
            return hits[0]
    return None


def _load_cohort(hist: str) -> pd.DataFrame:
    f = _find([f"TCGA-{hist}.star_counts.tsv.gz", f"TCGA-{hist}.star_counts.tsv",
               f"TCGA-{hist}.htseq_counts.tsv.gz", f"TCGA-{hist}.htseq_counts.tsv",
               f"*{hist}*counts*.tsv*"])
    if f is None:
        raise FileNotFoundError(f"No expression matrix for {hist} in {C.RAW}")
    print(f"  [{hist}] {f.name}", flush=True)
    df = pd.read_csv(f, sep="\t", index_col=0)
    df.index.name = "gene_id"
    # invert log2(count+1) -> counts, store float32 to save memory
    df = np.rint(np.power(2.0, df) - 1.0).clip(lower=0).astype("float32")
    return df


def _symbol_map():
    pm_file = _find(["gencode*.probeMap", "*.probeMap", "*probemap*"])
    if pm_file is None:
        return None
    sym = U.load_probemap(pm_file)
    sym.index = sym.index.astype(str).str.split(".").str[0]
    return sym[~sym.index.duplicated()]


def _collapse_to_symbols(counts, sym):
    stripped = counts.index.astype(str).str.split(".").str[0]
    mapped = pd.Series(sym.reindex(stripped).values, index=counts.index)
    unmapped = mapped.isna().mean() * 100
    out = counts.assign(_sym=mapped.values).dropna(subset=["_sym"])
    out = out.groupby("_sym").sum()
    out.index.name = "gene"
    print(f"  mapped to {out.shape[0]} unique symbols ({unmapped:.1f}% ids unmapped)",
          flush=True)
    return out.astype("float32")


def main():
    print("Loading cohorts …", flush=True)
    luad, lusc = _load_cohort("LUAD"), _load_cohort("LUSC")
    genes = luad.index.union(lusc.index)
    counts = pd.concat([luad.reindex(genes).fillna(0),
                        lusc.reindex(genes).fillna(0)], axis=1).astype("float32")
    counts = counts.loc[:, ~counts.columns.duplicated()]
    luad_cols = set(luad.columns)
    del luad, lusc; gc.collect()
    print(f"  merged: {counts.shape[0]} genes x {counts.shape[1]} samples", flush=True)

    samples = counts.columns
    labels = pd.DataFrame({
        "sample": samples,
        "patient": [U.tcga_patient(s) for s in samples],
        "histology": ["LUAD" if s in luad_cols else "LUSC" for s in samples],
        "sample_type": [U.tcga_sample_type(s) for s in samples],
    })
    n_lt = ((labels.sample_type == "tumor") & (labels.histology == "LUAD")).sum()
    n_st = ((labels.sample_type == "tumor") & (labels.histology == "LUSC")).sum()
    print(f"  tumor LUAD={n_lt}, tumor LUSC={n_st}, normal="
          f"{(labels.sample_type=='normal').sum()} (paper: 529 / 498 tumor)", flush=True)

    # ---- low-count filter (paper: drop genes CPM<1 in >600 tumor samples) ---
    # done without materializing a full float CPM matrix
    tumor_cols = labels.loc[labels.sample_type == "tumor", "sample"]
    tumor_cols = [c for c in tumor_cols if c in counts.columns]
    ct = counts[tumor_cols]
    lib = ct.sum(axis=0)                      # per-sample library size
    thresh = lib / 1e6                        # counts equivalent to CPM = 1
    below = ct.lt(thresh, axis=1)             # bool: CPM < 1
    n_low = below.sum(axis=1)
    keep = n_low <= C.MIN_SAMPLES_ABOVE_CPM
    del ct, below; gc.collect()
    counts = counts.loc[keep]
    print(f"  gene filter: {counts.shape[0]} genes retained (paper 14,010)", flush=True)

    # ---- map Ensembl -> symbol (now on the small filtered matrix) -----------
    sym = _symbol_map()
    if sym is not None:
        counts = _collapse_to_symbols(counts, sym)
    else:
        print("  WARNING: no probeMap; keeping Ensembl IDs", flush=True)

    # ---- TMM + log2-CPM ----------------------------------------------------
    print("TMM + log2-CPM …", flush=True)
    logcpm_tumor = U.tmm_logcpm(counts[tumor_cols])
    logcpm_all = U.tmm_logcpm(counts)

    logcpm_tumor.to_parquet(C.PROC / "expr_logcpm.parquet")
    logcpm_all.to_parquet(C.PROC / "expr_logcpm_all.parquet")
    counts[tumor_cols].to_parquet(C.PROC / "counts_filtered.parquet")
    labels.to_csv(C.PROC / "labels.csv", index=False)
    print(f"Done → data/processed/  (tumor {logcpm_tumor.shape}, all {logcpm_all.shape})",
          flush=True)


if __name__ == "__main__":
    main()
