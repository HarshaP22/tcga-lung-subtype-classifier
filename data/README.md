# data/

Nothing in this folder is committed except this file. Run `bash download_data.sh` from the
repository root to populate it (about 250 MB, all public). See [docs/DATA.md](../docs/DATA.md)
for what each file is.

```
data/
├── raw/          TCGA expression, clinical and survival files, probe maps, genesets/
├── external/     GSE28582 expression matrix (external validation)
└── processed/    written by the pipeline: normalized matrices, split, stage caches
```
