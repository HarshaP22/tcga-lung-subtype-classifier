# Validation assessment

This portfolio is an adapted reproduction and extension of Chen and Dhahbi
(2021), DOI [10.1038/s41598-021-92725-8](https://doi.org/10.1038/s41598-021-92725-8).
It demonstrates public-data engineering, feature selection, classification and
biological interpretation. The saved results support strong retrospective
subtype discrimination; they do not establish a clinical assay or demonstrate
that the overlap panel outperforms alternatives beyond sampling variation.

Read the [project research guide](pdf/PROJECT_GUIDE.pdf) for a detailed explanation
of the data, models, findings and limitations. Older REPORT, MODELS and BLOG
documents and their PDFs retain the original analysis narrative. This assessment
qualifies their stronger claims. Numerical result tables have not been changed.

## What the results measure

- The 339-gene overlap panel classified 299 of 309 binary test samples correctly:
  accuracy 0.9676, F1 0.9695 and AUROC 0.9887. PCA, mRMR and XGBoost panels each
  classified 298/309 correctly. L1 selection has the highest saved AUROC, 0.9913.
  A one-sample accuracy advantage does not establish superiority.
- Historical identifiers `overlap131` and `overlap17` contain 339 and 140 genes
  respectively. Their names refer to the source paper, not the current sizes.
- The external mean AUC of 0.9253 is the arithmetic mean of 134 **single-gene**
  AUCs in GSE28582. The complete multigene classifier was not transferred or
  evaluated externally. The code independently chooses favorable AUC direction
  in each cohort; a stronger test would freeze direction using training data.
- TCGA single-gene AUCs are evaluated on the full discovery cohort, including
  samples used for selection, and should not be described as held-out estimates.

## Evaluation boundaries

1. The binary outer split groups patients within histology, but filtering and
   normalization precede splitting. Test distributions influence preprocessing.
2. Feature selection is performed once before random-forest cross validation.
   Inner validation folds therefore influence the selected inputs. L1 scaling is
   also fitted before its cross validation. Fit these steps within each fold.
3. Inner cross validation is sample based. If patients have multiple samples,
   patient-grouped inner folds are needed.
4. Multiple panels, selector thresholds and ablation sizes are examined on the
   same test set. Those comparisons are development evidence; use nested
   evaluation and a new untouched validation cohort for final claims.
5. The three-class extension creates a different split after binary feature
   selection. Its split groups within each class independently, so tumour and
   normal samples from the same patient can cross the train/test boundary.
   This also affects the all-gene baseline. Assign each patient globally before
   splitting classes or selecting features.
6. Accuracy bootstrap intervals condition on fixed test predictions. They do
   not represent uncertainty from the full selection and tuning process.

## Survival and exploratory biology

The saved survival table has 280 gene/histology tests and 63 unadjusted log-rank
p-values below 0.05. Applying Benjamini-Hochberg once across those 280 saved
log-rank p-values gives 11 comparisons with q < 0.05. This is a calculation from
the saved table, not a rerun of survival models. The testing family should be
specified in advance for a confirmatory analysis.

| LUAD gene | Log-rank p | BH q across 280 | Univariate Cox HR |
| --- | ---: | ---: | ---: |
| PERP | 0.000208 | 0.0109 | 1.356 |
| KRT6A | 0.002771 | 0.0560 | 1.135 |
| KRT17 | 0.015520 | 0.1403 | 1.098 |

Cox p-values are separate tests and require their own prespecified multiplicity
handling. Current models do not adjust for stage, age or treatment and do not
report proportional-hazards checks. Survival associations are exploratory.

Agreement among feature selectors does not provide independent biological
confirmation or prove that false positives cancel out. Pathway enrichment is
association, SHAP explains model predictions rather than causation, and candidate
stability does not establish novelty. Clustering requires batch and purity checks
before assigning new molecular subtypes. Absence of detected sex association does
not establish absence of confounding.

## Differences from the source study

The implementation uses STAR inputs, a Welch-test differential-expression
approximation instead of edgeR, a custom ANOVA/correlation mRMR variant, and
different selection thresholds. XGBoost is capped at 148 genes; the L1 panel has
18 genes. Enrichr libraries and TCGA survival analysis replace DAVID and
KM-Plotter steps. These are material adaptations, not exact reproductions.

## Reproduction and next steps

Use an isolated environment, record exact input checksums and patient IDs, and
distinguish regenerated stages from reused caches. Synthetic tests check software
behavior; they do not independently reproduce the real-data results. Raw inputs
and processed caches are intentionally absent from Git.

Before stronger claims, rerun a globally patient-grouped, nested pipeline; fit
preprocessing and selection within training folds; freeze panel, transformation,
direction and threshold; and evaluate the complete classifier on an independent
cohort. Validate prognosis separately with clinical covariates and multiplicity
control.
