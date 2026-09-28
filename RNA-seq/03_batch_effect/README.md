# Variance decomposition with Bioconductor variancePartition

Added at revision, answering the editor —

> "Please quantify the experiment/batch effect in the RNA-seq dataset, for
> example using PCA or MDS, and show the relative contribution of experiment,
> treatment and time point."

— and Reviewer 4, comment 16:

> "a PCA or MDS with variance explained per axis is a minimum requirement."

## Method

`variance_partition.R` uses **variancePartition** as the sole
variance-decomposition engine (Hoffman & Schadt, *BMC Bioinformatics* 2016;
17:483; <https://bioconductor.org/packages/variancePartition>). Four steps, all
package-native:

1. **`canCorPairs()`** — canonical correlations between trial, time point and
   condition. All zero here (largest 6e-17), so the trial is *crossed* with the
   design factors, not confounded with them. This is what makes the trial
   effect separable regardless of its size, and the package asks you to check
   it before decomposing anything.
2. **`fitExtractVarPartModel()` on genes** — the package's primary product: the
   share of each gene's expression variance attributable to each factor.
   Summarised as median, mean, variance-weighted mean and quartiles;
   visualised with `plotVarPart()`.
3. **`fitExtractVarPartModel()` on principal component scores** — the same
   estimator applied axis by axis. Combined with the variance each axis
   carries, this is the "PCA with variance explained per axis" the reviewer
   asked for. Treating component scores as the features is legitimate: the
   model and the estimator are unchanged, only the response differs.
4. **`plotPercentBars()`** — worked examples: the ten genes most affected by
   trial, and the ten most affected by time point.

variancePartition reports no p-values, so one is added from the *same* model:
limma's moderated *F*-test across the coefficients coding each term, giving the
number of genes for which that factor explains more than noise. Effect size
from variancePartition, evidence from the F-test, one model fit.

One gene set throughout: `filterByExpr` against the full design, 18,425 genes,
no most-variable-gene selection and no `ntop` cut-off.

## Sum-to-zero contrasts are mandatory here

**This is the one thing to know before editing the script.**

For a fixed-effect term, variancePartition reports `var(X_j * beta_j)` as a
fraction of the total. Under R's **default treatment contrasts**, the columns
coding a main effect and those coding its interaction are correlated: the
shared variance gets counted twice, the denominator inflates, and the split
between main effect and interaction is wrong.

Measured on this dataset, the default coding attributes **0.01%** of PC1 to
time point. The correct value is **36.8%**.

`options(contrasts = c("contr.sum", "contr.poly"))` makes those columns
orthogonal for a balanced design, and the decomposition becomes the unique
orthogonal partition. The script sets this at the top, asserts the design is
balanced (orthogonality depends on it), and verifies the per-axis result
against a classical balanced ANOVA before writing anything —
`check_orthogonal_partition.tsv` records the agreement, currently 1.7e-14
percentage points. If that check ever fails, the contrast option has been
changed or the design is no longer balanced.

Random effects are **not** used: Experiment and Time have two levels and
Condition three, and `lme4` shrinks two-level variance components to
essentially zero (Time comes out at ~1e-31). variancePartition's own guidance
is that categorical variables with few levels enter as fixed effects.

## Outputs (`results/batch_effect/`)

| File | Contents |
|---|---|
| `table1_canCorPairs.tsv` | correlations between the experimental factors — evidence of no confounding |
| `table2_variancePartition_summary.tsv` | per-gene medians, means, variance-weighted means, quartiles |
| `table3_variancePartition_per_gene.tsv` | per-gene values, all 18,425 genes |
| `table4_variance_per_axis.tsv` | PC1–PC10: % variance and attribution — the reviewer's minimum requirement |
| `table5_significance_per_term.tsv` | genes with a significant effect per factor (moderated F-test, FDR < 0.05) |
| `check_orthogonal_partition.tsv` | agreement with the balanced-ANOVA partition |
| `design_balance.txt` | proof the design is balanced (5 per cell) |
| `Fig_PCA.pdf` | Figure 2d |
| `Fig_variance_per_axis.pdf` | Figure 2e |
| `Fig_variancePartition_violin.pdf` | Figure 2c |
| `Fig_canCorPairs.pdf` | factor-correlation matrix, supplementary |
| `Fig_percentBars_top_genes.pdf` | worked per-gene examples, supplementary |
| `figure2_panels_cd.rds` | the panel objects, read by both `04_figures/` builders |
| `sessionInfo.txt` | package versions for this run |

## Headline numbers

Per gene:

| Source | Median gene | Variance-weighted | IQR |
|---|---:|---:|---:|
| Time point | 11.7% | 17.6% | 3.0 – 27.3% |
| Trial (Experiment) | 9.9% | 13.1% | 2.5 – 23.7% |
| Trial × Time point | 9.2% | 17.4% | 2.3 – 22.2% |
| Condition (treatment) | 1.6% | 2.3% | 0.6 – 3.3% |
| Residual (between-plant) | 50.1% | 49.5% | 33.3 – 67.8% |

Per axis: PC1 (35.9% of total variance) is 36.8% time point, 19.1% trial,
35.4% trial × time point, 0.1% treatment. Trial, time and their interaction
account for 91% of PC1 and 97% of PC2.

Trial plus its interaction with time point: **30.6%** variance-weighted, and
27.2% for the median gene. (The variance-weighted figures are weighted means
and so add across terms; the medians do not — 27.2% is the median of the
per-gene sum, computed from `table3`, not 9.9 + 9.2.) Treatment never exceeds a
few percent under any summary.

Significance, by moderated F-test at FDR < 0.05: time point 12,333 of 18,425
genes (66.9%), trial 11,859 (64.4%), trial × time point 11,471 (62.3%),
treatment **0 (0.0%)**.

## Why no p-values on the per-axis table

Deliberate. Principal components are ordered by the variance they carry, so
testing them is selection-affected. On this dataset doing so manufactures a
nominally significant treatment effect on PC2 — 0.49% of that axis, 0.07% of
total variance, BH-adjusted FDR = 0.047 across the 40 term-by-axis tests —
which the gene-level test contradicts with zero significant genes. Significance
belongs at the gene level, where there is no selection step; "variance
explained per axis" is an effect-size question and the per-axis table answers
it as one.
