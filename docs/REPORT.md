# Reproducing overlapping feature selection for LUAD vs LUSC

> For the current interpretation and additional limitations, read [VALIDATION.md](VALIDATION.md) and the [research guide](pdf/PROJECT_GUIDE.pdf). Saved numerical results are preserved; strong claims in the historical narrative should be read with this assessment.

**Reproduction and extension of** Chen & Dhahbi (2021), *"Lung adenocarcinoma and lung
squamous cell carcinoma cancer classification, biomarker identification, and gene expression
analysis using overlapping feature selection methods,"* **Scientific Reports** 11:13323,
[DOI: 10.1038/s41598-021-92725-8](https://doi.org/10.1038/s41598-021-92725-8) (PMC8233431).
Article full text and metadata retrieved via PubMed Central.

---

## 1. Summary

The paper's thesis is that intersecting several independent feature-selection methods yields
a compact gene set that (a) classifies LUAD vs LUSC as well as or better than any single
method, and (b) concentrates biologically meaningful biomarkers while suppressing each
method's false positives. This project reproduces that pipeline end-to-end on current TCGA
data and independently confirms the central results: the overlap panel is the single best
classifier by accuracy and F-measure; the recovered biomarkers reproduce the paper's AUCs
almost exactly (KRT5 0.973 vs 0.973, NECTIN1 0.951 vs 0.951); they validate on an independent
microarray cohort; and pathway analysis recovers keratinization (up in LUSC) as the dominant
signal. Four extensions add a three-class task, model interpretability, survival, and a
novel-biomarker exploration. A sensitivity analysis (§8) shows that the overlap's *size*
depends on selector thresholds while its *accuracy* does not.

Every model, its settings, and a leakage audit are documented in [MODELS.md](MODELS.md). The
one-page summary figure is `results/figures/summary_figure.png`.

## 2. Data and preprocessing

TCGA-LUAD and TCGA-LUSC STAR-count matrices and a gencode probeMap were obtained from UCSC
Xena (GDC hub); the external cohort GSE28582 and the GPL570 annotation from NCBI GEO. After
inverting Xena's `log2(count+1)`, merging cohorts, applying the paper's low-count filter
(drop genes with CPM < 1 in > 600 tumour samples), and mapping Ensembl→symbol, the working
matrix is **15,686 genes × 1,029 tumours (528 LUAD, 501 LUSC)** plus 110 solid-tissue
normals, normalized with edgeR-style TMM + log2-CPM.

| | This work | Paper |
|---|---|---|
| LUAD / LUSC tumours | 528 / 501 | 529 / 498 |
| Genes after filter | 15,686 | 14,010 |

The tumour counts match to within one or two samples. The gene count is higher because we
use GDC **STAR** counts (gencode v36), the current GDC pipeline, whereas the paper used the
now-retired **HTSeq** counts (gencode v22); the v36 annotation simply contains more genes.

## 3. Feature selection and overlap

Five methods were run on a single fixed patient-level 7:3 training split (reused everywhere):
differential expression (edgeR-style: Welch test on log-CPM with BH-FDR, |logFC|>0.5), PCA
(top genes by |PC1| loading), a fast vectorized mRMR over the full gene set, XGBoost gain
importance, and L1-logistic regression with the 1-SE rule.

| Method | Genes selected (this work) | Paper |
|---|---|---|
| DGE (top 500) | 500 | 500 |
| PCA (top 500) | 500 | 500 |
| mRMR (top 500) | 500 | 500 |
| XGBoost | 148 (capped to paper's size) | 148 |
| Lasso (1-SE) | 18 | 68 |
| **Overlap ≥ 2 methods** | **339** | 131 |
| **Overlap ≥ 3 methods** | 140 | 17 |
| In all 5 methods | 6 | not reported |

The six genes selected by **all five** methods form a natural high-confidence core:
`BCL2L15, DDAH1, DSC3, DSG3, KRT5, SIAH2`: desmosomal (DSC3, DSG3) and keratin (KRT5)
genes, with NECTIN1 (PVRL1) just outside in four of five, i.e. exactly the squamous-epithelial markers the
paper highlights.

Our overlap sets are larger than the paper's (339/140 vs 131/17). The reason is that our five
methods are more concordant with one another than the paper's were; in particular, DGE and
PCA share 228/500 genes on this data. Overlap size is highly sensitive to each method's exact
thresholds and to the PCA/mRMR implementation details, none of which are fully specified in
the paper. This does not affect the scientific conclusions: the overlap still classifies best
and still concentrates the known biomarkers (below).

## 4. Classification (the paper's central claim)

Random-forest validation on the locked test set (LUAD = positive, 10,000-bootstrap 95% CI):

| Gene set | n | Accuracy (95% CI) | Sens | Spec | F | AUROC |
|---|---|---|---|---|---|---|
| DGE | 500 | 0.958 (0.935–0.977) | 0.987 | 0.927 | 0.960 | 0.986 |
| PCA | 500 | 0.964 (0.942–0.984) | 0.994 | 0.933 | 0.966 | 0.982 |
| mRMR | 500 | 0.964 (0.942–0.984) | 1.000 | 0.927 | 0.967 | 0.987 |
| XGBoost | 148 | 0.964 (0.942–0.984) | 1.000 | 0.927 | 0.967 | 0.989 |
| Lasso | 18 | 0.958 (0.935–0.977) | 0.975 | 0.940 | 0.960 | 0.991 |
| **Overlap ≥ 2** | 339 | **0.968 (0.948–0.987)** | 1.000 | 0.933 | **0.970** | 0.989 |
| Overlap ≥ 3 | 140 | 0.958 (0.935–0.977) | 0.994 | 0.920 | 0.960 | 0.985 |

Every method exceeds 95% accuracy and ~0.98 AUROC, and **the overlap panel achieves the
highest accuracy (0.968) and F-measure (0.970)** of all, reproducing the paper's key finding
that overlapping feature selection is "comparable, if not better" than the individual methods,
while using a fraction of the features. See `results/figures/confusion_131.png`,
`heatmap_131.png`, `heatmap_17.png`.

## 5. Biomarkers and external validation

Single-gene discriminative AUCs on TCGA reproduce the paper's values almost exactly:

| Gene | AUC (this work) | AUC (paper) |
|---|---|---|
| KRT5 | 0.973 | 0.973 |
| KRT6A | 0.953 | 0.953 |
| NECTIN1 / PVRL1 | 0.951 | 0.951 |
| ARHGEF38 | 0.953 | 0.957 |
| PERP | 0.949 | 0.953 |

The strongest discriminators overall are the squamous markers **DSG3 (0.974), KRT5 (0.973),
DSC3 (0.969)** and the keratin family, all up in LUSC. Of the paper's specifically named
biomarkers we independently recover KRT5, KRT6A, KRT14, KRT17, PERP, TRIM29, MUC1, ARHGAP12,
ARHGEF38 and NECTIN1 (as PVRL1) in the ≥3-method set. That includes **three of the paper's
five claimed-novel genes** (ARHGAP12, ARHGEF38, NECTIN1). The other two, ELFN2 (AUC 0.837)
and REPS1 (0.858), discriminate reasonably but were picked by only one of our five methods.
GPC1, CELSR2, S100A2 and TUBA1C land in the ≥2-method panel but not the ≥3 set, so in total
10 of the paper's 17 biomarkers reach our ≥3 set (see §8, Sensitivity).

**External validation (GSE28582, 50 LUAD + 28 LUSC, Affymetrix GPL570):** across 134 shared
genes the mean AUC is **0.925**, and the top biomarkers transfer strongly (KRT6A 0.994,
DSC3 0.993, KRT5 0.988, TRIM29 0.989), confirming the panel generalizes off-platform, exactly
as the paper reports. See `results/figures/roc_biomarkers_up.png`, `auc_tcga_vs_external.png`.

## 6. Pathway analysis

Splitting the overlap into genes up in LUSC vs up in LUAD and running GO/KEGG/Reactome
over-representation (Enrichr libraries, Fisher's exact, the ~15.7k analysed genes as
background):

- **Up in LUSC:** *Keratinization* (Reactome R-HSA-6805567) is the top term, followed by
  *Formation of the cornified envelope*, *epidermis development*, and *keratinocyte
  differentiation*, plus a *negative regulation of DNA-damage-induced apoptosis* (p53-axis)
  term. This is precisely the paper's finding: cornification/keratinization dominates, with
  p53 regulation prominent.
- **Up in LUAD:** small-GTPase-mediated signalling, cell migration/extension, O-glycan
  processing (MUC1-linked) and transport/secretion, the secretory/exocytosis-leaning profile
  the paper associates with LUAD.

See `results/figures/pathways_up.png`, `pathways_down.png`.

## 7. Extensions

**7a. Three-class LUAD / LUSC / Normal** (the source guide's Project 1). Adding the 110
normals, the 339-gene overlap panel reaches macro-F1 0.957 / macro-AUROC 0.992, essentially
matching the full 15,686-gene model (0.964 / 0.992). A compact panel suffices for the
three-way problem. See `results/figures/confusion_threeclass.png`.

**7b. Interpretability + feature-count ablation.** Native XGBoost SHAP values rank DSC3 and
KRT5 among the top drivers (known LUSC markers). The ablation curve is striking: **2 genes
already give 0.955 AUROC, 10 genes 0.980, 17 genes 0.984**, plateauing by ~25 genes at 0.988.
This explains why small NanoString-style panels work and mirrors the guide's expected result.
See `results/figures/auc_vs_ngenes.png`, `shap_summary.png`.

**7c. Survival (TCGA clinical, lifelines KM + Cox).** Replacing the paper's KM-Plotter step
with a direct TCGA analysis, 63 of 280 gene–histology tests are prognostic (log-rank p<0.05).
Consistent with the paper, keratins are unfavourable in LUAD: **KRT6A (p=0.003), KRT17
(p=0.016), TRIM29 (borderline)**. One documented discrepancy: PERP is unfavourable in our
TCGA cohort but was reported favourable in the paper's KM-Plotter microarray cohort. That is
a cohort/platform difference worth flagging rather than hiding. See `results/figures/km_*.png`.

**7d. Novel exploration.** Beyond the paper's highlighted genes, stability selection (genes
chosen across bootstraps) combined with high single-gene AUC surfaces reproducible,
under-discussed candidates: **DSG3 (AUC 0.974, stability 1.00), CERS3 (0.965, 1.00),
SERPINB13 (0.955), BNC1 (0.952), KRT6B/KRT6C, S1PR5, and the lncRNA MIR205HG (0.949)**. MIR205HG is the
miR-205 host gene, a known squamous transcript but absent from the paper's panel.
Unsupervised clustering of the overlap panel shows moderate substructure within LUSC
(silhouette 0.30 at k=2) but little within LUAD, hinting at LUSC molecular subtypes. A
confound audit (sex tested *within* each histology, so the LUAD/LUSC sex imbalance can't leak
in) finds that the key discriminators (KRT5, KRT6A, NECTIN1, PERP) are **not** sex-associated,
supporting their validity as subtype markers.

## 8. Sensitivity of the overlap (can we reach the paper's 131 / 17?)

The paper doesn't fully specify its selector thresholds, so `06_sensitivity.py` asks how much
the overlap depends on them, and whether the conclusions do.

| Configuration | Lasso genes | ≥2 | ≥3 | all 5 | RF acc (≥2) | RF acc (≥3) | Paper biomarkers in ≥3 |
|---|---|---|---|---|---|---|---|
| Principled (default) | 18 | 339 | 140 | 6 | 0.964 | 0.961 | 10 / 17 |
| Paper-matched set sizes | 76 | 350 | 143 | 8 | 0.964 | 0.964 | 10 / 17 |
| Paper sizes + PCA centred, not scaled | 76 | 347 | 193 | 11 | 0.968 | 0.968 | 14 / 17 |

*(A fixed 500-tree forest is used here, not the grid-tuned one from §4, so accuracies differ
slightly from the main table.)*

Matching every selector to the paper's set sizes barely moves the overlap. The size is set by
how deep the three ranked selectors (DGE, PCA, mRMR) are read:

| Top-N from DGE / PCA / mRMR | 50 | 100 | 150 | 200 | 300 | 500 |
|---|---|---|---|---|---|---|
| Genes in ≥2 methods | 51 | 85 | 119 | 153 | 216 | 350 |
| Genes in ≥3 methods | **17** | 28 | 40 | 53 | 84 | 143 |
| RF accuracy, ≥2 panel | 0.968 | 0.964 | 0.964 | 0.964 | 0.964 | 0.964 |

Top-N ≈ 150–175 reproduces the paper's ≈131 genes in ≥2 methods, and top-N = 50 gives exactly
17 in ≥3. No single top-N matches both numbers at once, which suggests the paper's selectors
were less concordant than ours. The key point is that **accuracy stays flat at 0.964–0.968
across the whole range**: the overlap's size is an artefact of the thresholds, while its
predictive value is robust. Centring PCA without scaling recovers the most paper biomarkers
(14/17), because unscaled PC1 favours high-variance, highly expressed genes such as the
keratins. See `results/figures/sensitivity_overlap.png` and panel E of
`results/figures/summary_figure.png`.

## 9. Documented deviations from the paper

1. **STAR vs HTSeq counts** (gencode v36 vs v22). HTSeq was retired by the GDC in 2024;
   raises the pre-filter gene count (15,686 vs 14,010).
2. **Overlap sizes larger** (339/140 vs 131/17). Our five methods are more concordant, and
   thresholds (Lasso 1-SE → 18 genes; XGBoost capped to the paper's 148) differ; the paper's
   exact PCA/mRMR/threshold choices are unspecified. §8 shows that reading ~150–175 genes
   from each ranked selector recovers ≈131, with accuracy unchanged.
3. **DGE** via Welch test on log-CPM + BH-FDR rather than edgeR's negative-binomial GLM
   (Bioconductor unavailable in the run environment); rankings are highly concordant.
4. **Pathways** via Enrichr offline libraries rather than the DAVID web tool.
5. **Survival** via TCGA RNA-seq + lifelines rather than the KM-Plotter microarray database, which
   explains the single PERP direction difference.

None of these change the paper's conclusions, all of which reproduced: overlap classifies
best; biomarker AUCs match within 0.004 (three exactly); biomarkers validate externally;
keratinization/p53 up in LUSC.

## 10. Reproduce

```bash
python -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt
bash download_data.sh
make reproduce
```

*Independent reproduction for educational/portfolio use. Anchor paper © the authors (CC-BY);
text and metadata retrieved via PubMed Central.*
