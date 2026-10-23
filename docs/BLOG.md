# What five feature selectors agree on: reproducing a lung-cancer biomarker study

> For the current interpretation and additional limitations, read [VALIDATION.md](VALIDATION.md) and the [research guide](pdf/PROJECT_GUIDE.pdf). Saved numerical results are preserved; strong claims in the historical narrative should be read with this assessment.

Lung cancer is not one disease. Its two dominant non-small-cell subtypes, adenocarcinoma
(LUAD) and squamous cell carcinoma (LUSC), look different under the microscope, arise from
different cells, carry different mutations, and increasingly get different drugs. Yet they are
still often lumped together, and clean molecular markers that separate them are scarcer than
you'd expect. A 2021 *Scientific Reports* paper by Chen and Dhahbi proposed a simple, almost
old-fashioned idea to find such markers: instead of trusting any one machine-learning feature
selector, run several and keep only the genes they *agree* on. I reproduced their pipeline
from scratch on current public data. It held up remarkably well, and it taught me a few things
about doing genomics ML carefully.

## The idea: intersection as a noise filter

Every feature-selection method has a characteristic way of being wrong. Differential
expression testing is sensitive to outliers and throws off false positives. Gradient-boosted
trees (XGBoost) overfit. LASSO is unstable when genes are correlated: it picks one of a
correlated cluster almost at random. PCA finds directions of maximum variance that may have
nothing to do with your biological question and are hard to interpret. mRMR greedily balances
relevance against redundancy but is, well, greedy.

The paper's move is to treat these weaknesses as independent noise. A false positive that
fools XGBoost is unlikely to also fool LASSO *and* differential expression *and* mRMR. So if
you take five methods and keep the genes selected by several of them, the survivors have
passed multiple, differently-shaped filters. The intersection should be enriched for real
signal and depleted of each method's private mistakes. Genes found by ≥2 methods go into the
classifier and pathway analysis; genes found by ≥3 become candidate biomarkers.

It's a satisfying idea because it needs no new algorithm, just discipline.

## Getting the data right is most of the work

The dataset is ~1,000 tumours from The Cancer Genome Atlas: about 500 LUAD and 500 LUSC,
each measured across ~20,000 genes by RNA-seq. Two preprocessing decisions matter more than
any modeling choice that follows.

**Normalization.** RNA-seq gives you counts, and counts are not comparable across samples
until you correct for sequencing depth and composition. It's tempting to reach for TPM
(transcripts per million), because that's what everyone quotes. But TPM is designed for
comparing genes *within* a sample, not the same gene *across* samples for machine learning;
it leaves a mean–variance dependency that models latch onto as fake signal. The right tool is
library-size normalization with a between-sample factor (edgeR's TMM) followed by a log
transform. I used TMM + log2-CPM, matching the paper. (VST from DESeq2 is an equally valid
choice; TPM is the trap.)

**Splitting.** TCGA has multiple aliquots per patient, and a naive random split can put two
samples from the *same person* in both train and test, a subtle leak that inflates every
metric. So the split has to be by *patient*, not by sample, and stratified by subtype. Do
this once, save it, and reuse the identical split for every method so the comparison is fair.
These two habits, proper normalization refit on training data only and per-patient splits,
are the difference between a portfolio project that looks honest to a reviewer and one that
doesn't.

One wrinkle: the paper used the GDC's HTSeq counts, which the GDC retired in 2024. I used the
current STAR counts instead and documented the swap. This is the kind of thing that quietly
breaks "just re-run the old code" reproductions, and it's worth being explicit about.

## Did it reproduce? Yes, strikingly so

Five selectors, one shared training split, and then the intersection. A random forest scored
each gene set on the untouched test set, with 10,000-bootstrap confidence intervals.

Every method cleared 95% accuracy and ~0.98 AUROC. LUAD and LUSC really are that separable
from bulk expression. The headline result reproduced cleanly: **the overlap panel had the
highest accuracy (0.968) and F-measure (0.970) of all**, edging out every individual method
while using far fewer genes. Agreement really does beat any single selector.

The biomarkers were the most convincing part. The paper reports single-gene AUCs for its top
markers; mine landed almost on top of theirs: KRT5 at 0.973 (paper: 0.973), KRT6A 0.953
(0.953), NECTIN1 0.951 (0.951), ARHGEF38 0.953 (0.957). When independent code on refreshed
data lands on the same third decimal, you believe the finding. And the genes make biological
sense: the strongest discriminators are keratins (KRT5/6A/14/17), desmosomal cadherins (DSC3,
DSG3) and nectins, the machinery of squamous epithelium. LUSC lights them up; LUAD doesn't.
The model isn't finding a statistical fluke; it's rediscovering what a pathologist sees.

To make sure this wasn't a TCGA-specific artifact, I validated the markers on a completely
separate microarray cohort (GSE28582, 78 patients, a different measurement technology
entirely). Mean AUC across shared genes was 0.925, with the top keratins still above 0.99.
The signal survives a platform change, the strongest possible sign it's real.

Finally, pathway enrichment on the overlap genes put *keratinization* and *cornified envelope
formation* at the very top of the LUSC-up list, with a p53/DNA-damage term close behind,
matching the paper, and matching the biology (TP53 mutations are more common in LUSC). The
LUAD side leaned toward GTPase signalling, cell migration and glycosylation.

## Where I pushed further

Reproducing is table stakes; a portfolio wants a little more. Four extensions:

**Three classes.** Adding normal lung tissue and doing LUAD-vs-LUSC-vs-normal, the compact
339-gene panel reached macro-F1 0.957, basically tying the full 15,000-gene model. You don't
need the whole transcriptome.

**How few genes?** I ranked genes by importance and asked how AUROC grows with panel size. Two
genes already give 0.955; ten give 0.980; seventeen give 0.984, and it plateaus by
twenty-five. This is exactly why clinical assays like NanoString panels can get away with
50–200 genes: the discriminative information is concentrated in a handful.

**Survival.** Using TCGA clinical data with Kaplan–Meier and Cox models, several markers were
prognostic: high keratin expression (KRT6A, KRT17) tracked with worse outcomes in LUAD, as in
the paper. One gene, PERP, pointed the opposite way from the paper, a reminder that survival
associations are cohort-dependent (they used a different microarray database), and worth
reporting rather than quietly dropping.

**Something new?** Stability selection (keeping genes chosen consistently across
resampling) combined with high discrimination surfaced reproducible candidates the paper
didn't spotlight: DSG3, CERS3, SERPINB13, BNC1, and notably MIR205HG, the long non-coding
host gene of miR-205, a known squamous transcript. I also checked whether the markers were
secretly just tracking patient sex (LUSC skews male); testing *within* each subtype, the key
markers were clean. And LUSC showed modest molecular substructure that LUAD didn't, a thread
worth pulling on another day.

## What I'd tell someone starting this

Three things carried the project. First, spend your care budget on data hygiene:
normalization and per-patient splitting decide whether your impressive numbers are real.
Second, a strong, boring baseline (a random forest on well-chosen genes) is astonishingly
hard to beat here, and saying so plainly is a feature, not a failure. Third, always ask what
the model learned: when the top features turn out to be the exact keratins a pathologist would
name, you've done more than fit a classifier; you've recovered biology. The overlapping-
selection idea works not because it's clever, but because it's honest about the fact that no
single method deserves your full trust.

*Full methods, tables, and figures are in the project's REPORT.md and results/ folder.
Reproduction of Chen & Dhahbi, Sci Rep 2021 (DOI 10.1038/s41598-021-92725-8), retrieved via
PubMed Central.*
