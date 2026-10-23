"""Extension C: survival / prognosis (replaces the paper's KM-Plotter step).

For each of the 17 biomarkers, split patients at the median expression within
each histology and test overall-survival difference (log-rank) + univariate
Cox. This reproduces, on TCGA clinical data, the paper's Kaplan-Meier finding
that several genes are prognostic, particularly in LUAD.

Requires TCGA survival files in data/raw/ (see docs/DATA.md §B). Skipped if
absent.

Outputs
-------
results/tables/survival_biomarkers.csv
results/figures/km_<gene>_<hist>.png   for significant genes
"""
from __future__ import annotations
import sys, json, glob
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parent))
import config as C
import utils as U


def load_survival():
    frames = []
    for f in sorted(C.RAW.glob("*survival*.tsv*")):
        df = pd.read_csv(f, sep="\t")
        frames.append(df)
    if not frames:
        return None
    surv = pd.concat(frames, ignore_index=True)
    # Xena survival cols: sample, OS, OS.time, _PATIENT
    cols = {c.lower(): c for c in surv.columns}
    s = pd.DataFrame({
        "sample": surv[cols.get("sample", "sample")],
        "OS": surv[cols.get("os", "OS")],
        "OS_time": surv[cols.get("os.time", cols.get("os_time", "OS.time"))],
    })
    return s.dropna()


def main():
    logcpm = pd.read_parquet(C.PROC / "expr_logcpm.parquet")
    labels = pd.read_csv(C.PROC / "labels.csv").set_index("sample")
    sets = json.load(open(C.TAB / "feature_sets.json"))
    genes = sets["overlap17"]
    surv = load_survival()
    if surv is None:
        print("No survival files in data/raw: survival extension skipped.")
        return

    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test

    surv = surv.set_index("sample")
    common = [s for s in logcpm.columns if s in surv.index]
    hist = labels.loc[common, "histology"]
    rows = []
    for h in ["LUAD", "LUSC"]:
        samp = [s for s in common if hist.loc[s] == h]
        sd = surv.loc[samp]
        for g in genes:
            if g not in logcpm.index:
                continue
            expr = logcpm.loc[g, samp]
            hi = expr >= expr.median()
            lr = logrank_test(sd.OS_time[hi], sd.OS_time[~hi],
                              sd.OS[hi], sd.OS[~hi])
            # Cox (continuous)
            cph_df = pd.DataFrame({"T": sd.OS_time.values, "E": sd.OS.values,
                                   "x": expr.values})
            try:
                cph = CoxPHFitter().fit(cph_df, "T", "E")
                hr = float(np.exp(cph.params_["x"])); cp = float(cph.summary.p["x"])
            except Exception:
                hr, cp = np.nan, np.nan
            rows.append({"gene": g, "histology": h, "logrank_p": lr.p_value,
                         "cox_HR": hr, "cox_p": cp,
                         "direction": "high=worse" if hr > 1 else "high=better"})
            if lr.p_value < 0.05:
                kmf = KaplanMeierFitter()
                plt.figure(figsize=(5, 4))
                for mask, lbl in [(hi, "high"), (~hi, "low")]:
                    kmf.fit(sd.OS_time[mask], sd.OS[mask], label=f"{g} {lbl}")
                    kmf.plot(ci_show=False)
                plt.title(f"{g} in {h}  (log-rank p={lr.p_value:.3g})")
                plt.xlabel("days"); plt.ylabel("OS")
                plt.savefig(C.FIG / f"km_{g}_{h}.png", dpi=150, bbox_inches="tight")
                plt.close()
    res = pd.DataFrame(rows).sort_values("logrank_p")
    res.to_csv(C.TAB / "survival_biomarkers.csv", index=False)
    sig = res[res.logrank_p < 0.05]
    print(f"Prognostic (log-rank p<0.05): {len(sig)} gene-histology pairs")
    print(sig.to_string(index=False))


if __name__ == "__main__":
    main()
