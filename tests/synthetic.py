"""Build a small synthetic cohort in the same file formats as the real downloads.

Layout written under <root>/data: Xena-style log2(count+1) STAR matrices with
versioned Ensembl IDs and TCGA barcodes, a gencode probeMap, clinical and survival
tables, a GEO series matrix with a GPL570 annotation table, and .gmt gene sets.

A known LUSC-up signal is planted on genes named after real squamous markers, so
the tests can check that the pipeline recovers it.
"""
from __future__ import annotations
import gzip
from pathlib import Path
import numpy as np
import pandas as pd

SIGNAL = ["KRT5", "KRT6A", "DSG3", "DSC3", "TP63", "PERP", "TRIM29", "KRT17"]
N_GENES = 1500


def build(root: Path, seed: int = 0) -> None:
    rng = np.random.default_rng(seed)
    raw, ext = root / "data" / "raw", root / "data" / "external"
    (raw / "genesets").mkdir(parents=True, exist_ok=True)
    ext.mkdir(parents=True, exist_ok=True)

    # gene ids: counts use one version suffix, the probeMap another (like v36 vs v22)
    acc = [f"ENSG{100000 + i:011d}" for i in range(N_GENES)]
    ids_counts = [f"{a}.{(i % 7) + 2}" for i, a in enumerate(acc)]
    ids_map = [f"{a}.1" for a in acc]
    symbols = SIGNAL + [f"GENE{i}" for i in range(N_GENES - len(SIGNAL))]
    signal_idx = list(range(len(SIGNAL))) + list(range(20, 50))   # named + extra signal

    def cohort(proj, n_tumour, n_normal, lusc):
        cols = ([f"TCGA-{proj}-{1000 + i:04d}-01A" for i in range(n_tumour)] +
                [f"TCGA-{proj}-{1000 + i:04d}-11A" for i in range(n_normal)])
        base = rng.negative_binomial(5, 0.02, size=(N_GENES, len(cols))).astype(float)
        base[60:80, :n_tumour] *= 3                  # tumour-vs-normal genes
        if lusc:
            base[signal_idx, :n_tumour] *= 6         # squamous signal
        return pd.DataFrame(np.log2(base + 1), index=ids_counts, columns=cols)

    luad, lusc = cohort("05", 60, 10, False), cohort("22", 50, 10, True)
    for name, df in [("LUAD", luad), ("LUSC", lusc)]:
        df.to_csv(raw / f"TCGA-{name}.star_counts.tsv.gz", sep="\t",
                  index_label="Ensembl_ID", compression="gzip")
        n = df.shape[1]
        pd.DataFrame({"sample": df.columns, "OS.time": rng.integers(30, 3000, n),
                      "OS": rng.integers(0, 2, n), "_PATIENT": [c[:12] for c in df.columns]}
                     ).to_csv(raw / f"TCGA-{name}.survival.tsv.gz", sep="\t", index=False)
        pd.DataFrame({"sample": df.columns,
                      "gender.demographic": rng.choice(["male", "female"], n),
                      "ajcc_pathologic_stage.diagnoses": rng.choice(["Stage I", "Stage II"], n)}
                     ).to_csv(raw / f"TCGA-{name}.clinical.tsv.gz", sep="\t", index=False)

    pd.DataFrame({"id": ids_map, "gene": symbols, "chrom": "chr1", "chromStart": 0,
                  "chromEnd": 1, "strand": "+"}).to_csv(
        raw / "gencode.v22.annotation.gene.probeMap", sep="\t", index=False)

    with open(raw / "genesets" / "SYNTH.gmt", "w") as f:
        f.write("squamous_signal\tNA\t" + "\t".join(SIGNAL + [symbols[i] for i in range(20, 50)]) + "\n")
        f.write("background_set\tNA\t" + "\t".join(symbols[200:260]) + "\n")

    # GEO: GPL570-style annotation (plain-text platform table) + series matrix
    probes = [f"{i}_at" for i in range(300)]
    psym = symbols[:300]
    with open(raw / "GPL570.annot.txt", "w") as f:
        f.write("^PLATFORM = GPL570\n!platform_table_begin\nID\tGene Symbol\n")
        for p, g in zip(probes, psym):
            f.write(f"{p}\t{g}\n")
        f.write("!platform_table_end\n")
    gsm = [f"GSM{i}" for i in range(78)]
    hist = ["adenocarcinoma"] * 50 + ["squamous cell carcinoma"] * 28
    expr = rng.normal(6, 1, size=(300, 78))
    for i in signal_idx:
        if i < 300:
            expr[i, 50:] += 3
    with gzip.open(ext / "GSE28582-GPL570_series_matrix.txt.gz", "wt") as f:
        f.write("!Sample_title\t" + "\t".join(f'"{g}"' for g in gsm) + "\n")
        f.write("!Sample_source_name_ch1\t" + "\t".join(f'"NSCLC {h}"' for h in hist) + "\n")
        f.write("!series_matrix_table_begin\n")
        f.write('"ID_REF"\t' + "\t".join(f'"{g}"' for g in gsm) + "\n")
        for p, row in zip(probes, expr):
            f.write(f'"{p}"\t' + "\t".join(f"{v:.3f}" for v in row) + "\n")
        f.write("!series_matrix_table_end\n")
