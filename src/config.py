"""Central configuration: paths, constants, and paper-matched parameters.

Reproduction of Chen & Dhahbi (2021), Sci Rep 11:13323,
"LUAD and LUSC cancer classification, biomarker identification, and gene
expression analysis using overlapping feature selection methods."
DOI: 10.1038/s41598-021-92725-8
"""
import os
from pathlib import Path

# ---- paths -----------------------------------------------------------------
# repo root; LUNG_ROOT lets the test suite point the pipeline at synthetic data
ROOT = Path(os.environ.get("LUNG_ROOT", Path(__file__).resolve().parents[1]))
RAW = ROOT / "data" / "raw"
EXTERNAL = ROOT / "data" / "external"
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
FIG = RESULTS / "figures"
TAB = RESULTS / "tables"
for _p in (PROC, FIG, TAB):
    _p.mkdir(parents=True, exist_ok=True)

# ---- reproducibility -------------------------------------------------------
SEED = 42

# ---- preprocessing (paper Methods: "RNA-Seq data processing") --------------
# "Genes < 1 CPM in over 600 samples were considered noise and discarded to
#  obtain 14,010 genes."  ~529 LUAD + 498 LUSC = 1027 samples in the paper.
CPM_MIN = 1.0
MIN_SAMPLES_ABOVE_CPM = 600          # keep gene if CPM>=1 in >= (N-600) samples
# implemented as: discard if (CPM<1) in > 600 samples  -> keep otherwise

# ---- feature-selection sizes (paper Results: "Selection of genes") ---------
N_DGE = 500          # top 500 DEGs by lowest p-value
N_PCA = 500          # top 500 genes by |PC1 loading|
N_MRMR = 500         # top 500 genes by mRMR ranking
# XGBoost / Lasso: number is data-driven (paper got 148 / 68 at threshold 0.5)
XGB_IMPORTANCE_THRESH = 0.5   # relative-to-max importance cut (see script notes)
XGB_TOP_N = 148               # cap XGBoost genes to the paper's selection size
LASSO_SELECTION = "nonzero"   # genes with non-zero L1 coefficient

# DGE filter (paper: edgeR, FDR<5e-2 and |logFC|>0.5 -> 4702 DEGs)
DGE_FDR = 5e-2
DGE_LOGFC = 0.5

# ---- classification / validation (paper Methods) ---------------------------
TEST_SIZE = 0.30              # 7:3 train/test
CV_FOLDS = 5                 # 5-fold CV for hyperparameters
N_BOOTSTRAP = 10000          # bootstrap replicates for accuracy CI
# In the paper LUAD = positive, LUSC = negative for binary metrics.
POS_LABEL = "LUAD"
NEG_LABEL = "LUSC"

# overlap thresholds
OVERLAP_CLASSIFY = 2         # genes found by >=2 methods  -> ~131 genes
OVERLAP_BIOMARKER = 3        # genes found by >=3 methods  -> ~17 genes
# NOTE: feature_sets.json keys "overlap131" / "overlap17" are named after the
# paper's set sizes. In this reproduction they hold 339 and 140 genes.

# ---- the paper's headline gene sets (for validation/comparison) ------------
# 17 candidate biomarkers (>=3 methods) reported in the paper. NOTE: reconstructed
# from the article main text (the exact 17 live in a supplementary table we don't
# have); used only to score how many our pipeline independently recovers.
PAPER_17_BIOMARKERS = [
    "NECTIN1", "PERP", "KRT5", "KRT6A", "ARHGEF38", "KRT17", "S100A2",
    "TRIM29", "REPS1", "GPC1", "ELFN2", "ARHGAP12", "QSOX1", "CELSR2",
    "KRT14", "TUBA1C", "MUC1",
]
# Five claimed-novel biomarkers.
PAPER_NOVEL = ["ARHGAP12", "ARHGEF38", "ELFN2", "NECTIN1", "REPS1"]
