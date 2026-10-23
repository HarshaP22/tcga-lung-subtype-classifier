"""Shared helpers: IO, normalization, TCGA barcode parsing, metrics."""
from __future__ import annotations
import gzip
import numpy as np
import pandas as pd
from pathlib import Path


# ---------------------------------------------------------------- barcode ----
def tcga_sample_type(barcode: str) -> str:
    """Map a TCGA aliquot barcode to a coarse sample type from the sample code.

    Barcode: TCGA-XX-XXXX-01A-...  positions 14-15 = sample type code.
      01 Primary Tumor, 02 Recurrent, 06 Metastatic, 11 Solid Tissue Normal.
    """
    parts = barcode.split("-")
    if len(parts) < 4:
        return "unknown"
    code = parts[3][:2]
    return {
        "01": "tumor", "02": "recurrent", "03": "tumor", "06": "metastatic",
        "11": "normal",
    }.get(code, "other")


def tcga_patient(barcode: str) -> str:
    """Patient-level id (first 3 fields): TCGA-XX-XXXX."""
    return "-".join(barcode.split("-")[:3])


# ------------------------------------------------------------------- IO ------
def read_xena_matrix(path: Path) -> pd.DataFrame:
    """Read a Xena gene-expression matrix (genes x samples), first col = gene id.

    Xena GDC files are log2(count+1). We return the raw matrix as-is; callers
    decide whether to invert. Index = gene id (Ensembl, possibly versioned).
    """
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index.name = "gene_id"
    return df


def load_probemap(path: Path) -> pd.Series:
    """Xena probeMap -> Series mapping gene id (Ensembl, versioned) to symbol."""
    pm = pd.read_csv(path, sep="\t")
    idcol = pm.columns[0]
    symcol = "gene" if "gene" in pm.columns else pm.columns[1]
    return pm.set_index(idcol)[symcol]


def invert_log2p1(df: pd.DataFrame) -> pd.DataFrame:
    """Invert Xena's log2(x+1) transform back to (approximately integer) counts."""
    return np.rint(np.power(2.0, df) - 1.0).clip(lower=0)


# ----------------------------------------------------------- normalization ---
def tmm_logcpm(counts: pd.DataFrame, prior_count: float = 1.0) -> pd.DataFrame:
    """edgeR-style TMM normalization + log2-CPM. counts: genes x samples.

    Uses `conorm.tmm` for the TMM factors, then CPM on effective library sizes.
    Returns log2-CPM (genes x samples).
    """
    import conorm
    norm = conorm.tmm(counts)              # TMM-normalized counts (genes x samples)
    lib = norm.sum(axis=0)
    cpm = norm.divide(lib, axis=1) * 1e6
    return np.log2(cpm + prior_count)


def cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0)
    return counts.divide(lib, axis=1) * 1e6


# --------------------------------------------------------------- metrics -----
def binary_stats(y_true, y_pred, pos_label="LUAD"):
    """Accuracy/sens/spec/precision/F1 with LUAD=positive (paper convention)."""
    from sklearn.metrics import (accuracy_score, precision_score, recall_score,
                                 f1_score, confusion_matrix)
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    pos = pos_label
    tp = np.sum((y_true == pos) & (y_pred == pos))
    tn = np.sum((y_true != pos) & (y_pred != pos))
    fp = np.sum((y_true != pos) & (y_pred == pos))
    fn = np.sum((y_true == pos) & (y_pred != pos))
    sens = tp / (tp + fn) if (tp + fn) else np.nan     # recall / sensitivity
    spec = tn / (tn + fp) if (tn + fp) else np.nan
    prec = tp / (tp + fp) if (tp + fp) else np.nan
    f1 = 2 * prec * sens / (prec + sens) if (prec + sens) else np.nan
    acc = (tp + tn) / len(y_true)
    return dict(accuracy=acc, sensitivity=sens, specificity=spec,
                precision=prec, f_measure=f1, tp=tp, tn=tn, fp=fp, fn=fn)


def bootstrap_ci(y_true, y_pred, n=10000, seed=42, pos_label="LUAD"):
    """Bootstrap 95% CI for accuracy over paired (y_true,y_pred)."""
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    m = len(y_true)
    accs = np.empty(n)
    for i in range(n):
        idx = rng.integers(0, m, m)
        accs[i] = np.mean(y_true[idx] == y_pred[idx])
    return float(np.mean(accs)), float(np.percentile(accs, 2.5)), float(np.percentile(accs, 97.5))
