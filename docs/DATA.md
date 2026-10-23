# Data

All inputs are public and need no login or controlled access. Nothing is committed to the
repository. One command fetches everything (about 250 MB):

```bash
bash download_data.sh
```

The script skips files you already have, tries fallback URLs when a source has renamed a
file, checks every gzip, and prints a summary. The table below says what each file is and
why the pipeline needs it, so you can also fetch any file by hand.

## Files

| # | File | Saved to | Source | Used for |
|---|---|---|---|---|
| 1 | `TCGA-LUAD.star_counts.tsv.gz` | `data/raw/` | UCSC Xena, GDC hub | LUAD expression, the core input |
| 2 | `TCGA-LUSC.star_counts.tsv.gz` | `data/raw/` | UCSC Xena, GDC hub | LUSC expression, the core input |
| 3 | `gencode.v22.annotation.gene.probeMap` | `data/raw/` | [ucscXena/xena-GDC-ETL](https://github.com/ucscXena/xena-GDC-ETL) | Ensembl ID to gene symbol |
| 4 | `TCGA-LUAD.clinical.tsv.gz` | `data/raw/` | UCSC Xena, GDC hub | Sex and stage, for the confound audit |
| 5 | `TCGA-LUSC.clinical.tsv.gz` | `data/raw/` | UCSC Xena, GDC hub | Sex and stage, for the confound audit |
| 6 | `TCGA-LUAD.survival.tsv.gz` | `data/raw/` | UCSC Xena, GDC hub | Overall survival |
| 7 | `TCGA-LUSC.survival.tsv.gz` | `data/raw/` | UCSC Xena, GDC hub | Overall survival |
| 8 | `GSE28582-GPL570_series_matrix.txt.gz` | `data/external/` | NCBI GEO | External validation cohort |
| 9 | `GPL570.annot.txt` | `data/raw/` | NCBI GEO | Affymetrix probe to gene symbol |
| 10-12 | `GO_Biological_Process_2021.gmt`, `KEGG_2021_Human.gmt`, `Reactome_2022.gmt` | `data/raw/genesets/` | Enrichr | Pathway enrichment |

Files 1-3 are enough to run the core reproduction (steps 0-3). The others enable the
external validation, pathways, survival and confound analyses. Any step whose inputs are
missing is skipped with a message.

Cohort pages, if you'd rather browse:
- [GDC TCGA LUAD on UCSC Xena](https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Lung%20Adenocarcinoma%20(LUAD)&host=https%3A%2F%2Fgdc.xenahubs.net)
- [GDC TCGA LUSC on UCSC Xena](https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Lung%20Squamous%20Cell%20Carcinoma%20(LUSC)&host=https%3A%2F%2Fgdc.xenahubs.net)
- [GSE28582 on GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28582)

## Things worth knowing about these sources

- **Xena stores counts as `log2(count + 1)`.** `00_preprocess.py` inverts this back to raw
  counts before TMM normalization.
- **STAR, not HTSeq.** The paper used GDC HTSeq counts, which the GDC retired in 2024. STAR
  counts (gencode v36) are the current equivalent. The loader also accepts HTSeq files if
  you have them.
- **Gencode versions differ.** The counts use gencode v36 IDs while the probeMap uses v22.
  Only about 53% of versioned IDs match directly, but 93% match once the version suffix is
  stripped. The loader always matches on unversioned IDs.
- **GSE28582 is a SuperSeries.** It holds a SNP-array subseries and an expression
  subseries. Only the GPL570 expression matrix is used, and large-cell carcinomas are
  excluded, which leaves 50 LUAD and 28 LUSC samples.
- **Tumour vs normal comes from the TCGA barcode.** Sample code `-01` is primary tumour,
  `-11` is solid-tissue normal, and anything else (recurrent, metastatic) is excluded.
