"""End-to-end smoke test on synthetic data (no downloads).

Runs the real pipeline scripts against a small synthetic cohort (tests/synthetic.py)
by pointing config.ROOT at a temporary directory via LUNG_ROOT, then checks that
each stage produced its outputs and that the planted squamous signal is recovered.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent))
import synthetic  # noqa: E402

SRC = Path(__file__).resolve().parents[1] / "src"
STEPS = ["00_preprocess.py", "01_feature_selection.py", "02_classification.py",
         "03_biomarkers_roc.py", "04_pathways.py", "05a_three_class.py", "05c_survival.py"]


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    root = tmp_path_factory.mktemp("synthetic_repo")
    synthetic.build(root)
    env = {**os.environ, "LUNG_ROOT": str(root), "MPLBACKEND": "Agg"}
    for step in STEPS:
        res = subprocess.run([sys.executable, str(SRC / step)], env=env,
                             capture_output=True, text=True)
        assert res.returncode == 0, f"{step} failed:\n{res.stdout[-2000:]}\n{res.stderr[-2000:]}"
    return root


def test_preprocess_outputs(run):
    proc = run / "data" / "processed"
    expr = pd.read_parquet(proc / "expr_logcpm.parquet")
    labels = pd.read_csv(proc / "labels.csv")
    assert expr.shape[1] == 110                          # tumours only (60 + 50)
    assert (labels.sample_type == "normal").sum() == 20
    assert "KRT5" in expr.index                          # versioned IDs mapped to symbols


def test_split_is_patient_level(run):
    split = json.load(open(run / "data" / "processed" / "split.json"))
    patient = lambda s: "-".join(s.split("-")[:3])
    assert not {patient(s) for s in split["train"]} & {patient(s) for s in split["test"]}


def test_overlap_recovers_planted_signal(run):
    sets = json.load(open(run / "results" / "tables" / "feature_sets.json"))
    for m in ["DGE", "PCA", "mRMR", "XGBoost", "Lasso", "overlap131", "overlap17"]:
        assert len(sets[m]) > 0, m
    assert {"KRT5", "DSG3", "DSC3"} <= set(sets["overlap17"])


def test_classifier_beats_chance(run):
    m = pd.read_csv(run / "results" / "tables" / "classification_metrics.csv")
    assert (m.accuracy > 0.9).all()


def test_external_validation_ran(run):
    ext = pd.read_csv(run / "results" / "tables" / "biomarker_auc_external.csv", index_col=0)
    assert ext.loc["KRT5", "auc_external"] > 0.9


def test_downstream_steps_wrote_outputs(run):
    tabs = run / "results" / "tables"
    for f in ["enrichment_up.csv", "threeclass_metrics.csv", "survival_biomarkers.csv"]:
        assert (tabs / f).exists(), f
