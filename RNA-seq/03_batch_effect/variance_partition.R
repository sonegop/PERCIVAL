# ===========================================================================
# PERCIVAL — contribution of experimental trial (batch), time point and
# treatment to RNA-seq transcriptome variance, using Bioconductor
# variancePartition as the sole variance-decomposition engine.
#
# https://bioconductor.org/packages/variancePartition
# Hoffman & Schadt (2016) BMC Bioinformatics 17:483
#
# Answers:
#   Editor: "Please quantify the experiment/batch effect in the RNA-seq
#   dataset, for example using PCA or MDS, and show the relative contribution
#   of experiment, treatment and time point."
#   Reviewer 4, comment 16: "a PCA or MDS with variance explained per axis is
#   a minimum requirement."
#
# Four outputs, all from variancePartition:
#   1. canCorPairs()            — are the factors confounded with one another?
#   2. fitExtractVarPartModel() on genes — the per-gene decomposition, the
#      package's primary product; violin plot via plotVarPart().
#   3. fitExtractVarPartModel() on principal component scores — the same
#      decomposition applied axis by axis, which together with the variance
#      each axis carries answers the reviewer's per-axis requirement.
#   4. plotPercentBars()        — worked examples for individual genes.
#
# variancePartition reports no p-values. One is added, from the same model:
# limma's moderated F-test per term, giving the number of genes for which each
# factor explains more than noise. See section 4b for why this is done at the
# gene level and deliberately not per principal component.
#
# ---------------------------------------------------------------------------
# CONTRAST CODING IS NOT OPTIONAL HERE. READ BEFORE EDITING.
#
# For a fixed-effect term, variancePartition reports var(X_j * beta_j) as a
# fraction of the total. Under R's DEFAULT treatment contrasts the columns
# coding a main effect and those coding its interaction are correlated, so the
# shared variance is counted twice, the denominator inflates, and the split
# between main effect and interaction is wrong. On this dataset the default
# coding attributes 0.01% of PC1 to time point; the correct value is 36.8%.
#
# Sum-to-zero coding (contr.sum) makes those columns orthogonal for a balanced
# design, and the decomposition then becomes the unique orthogonal partition —
# verified in this script against a classical balanced ANOVA, which it
# reproduces to numerical precision (see check_orthogonal_partition.tsv).
#
# The design is fully crossed and balanced: 3 Condition x 2 Time x
# 2 Experiment x 5 biological replicates = 60 libraries, five per cell.
# ---------------------------------------------------------------------------
#
# Inputs  : metadata/targets.txt, data/RNAseq_rawcounts.tsv
# Outputs : results/batch_effect/
# Usage   : Rscript 03_batch_effect/variance_partition.R
# ===========================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(DESeq2)
  library(variancePartition)
  library(ggplot2)
})

## Locate config/paths.R by walking up from this script (or the working dir).
.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.d <- if (length(.f)) dirname(normalizePath(sub("^--file=", "", .f[1]))) else normalizePath(getwd())
while (!file.exists(file.path(.d, "config", "paths.R")) && dirname(.d) != .d) .d <- dirname(.d)
if (!file.exists(file.path(.d, "config", "paths.R")))
  stop("Run from inside the repository, or set PERCIVAL_ROOT.", call. = FALSE)
source(file.path(.d, "config", "paths.R")); rm(.f, .d)

# See the header block: this line is load-bearing, not cosmetic.
options(contrasts = c("contr.sum", "contr.poly"))

OUT <- file.path(RES_DIR, "batch_effect")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(..., "\n", sep = "")

# ---------------------------------------------------------------------------
# 1. Metadata and counts
# ---------------------------------------------------------------------------
require_inputs(TARGETS_FILE, COUNTS_FILE)

targets <- read.delim(TARGETS_FILE, header = TRUE, stringsAsFactors = FALSE)
counts  <- as.matrix(read.delim(COUNTS_FILE, header = TRUE, row.names = 1,
                                check.names = FALSE))
stopifnot(all(c("Sample", "Condition", "Time", "Experiment") %in% colnames(targets)))
stopifnot(setequal(colnames(counts), targets$Sample))
counts <- counts[, targets$Sample]

targets$Condition  <- factor(targets$Condition, levels = c("Cntr", "E1_D3", "E1_D4"))
targets$Time       <- factor(targets$Time, levels = c("T1", "T2"))
targets$Experiment <- factor(paste0("Exp", targets$Experiment))
rownames(targets)  <- targets$Sample

balance <- with(targets, table(Condition, Time, Experiment))
capture.output({
  cat("Samples per Condition x Time x Experiment\n\n")
  print(ftable(balance))
  cat("\nAll cells equal to 5 replicates: ", all(balance == 5), "\n", sep = "")
  cat("Total libraries: ", sum(balance), "\n", sep = "")
}, file = file.path(OUT, "design_balance.txt"))
msg("Design balanced (all cells = 5): ", all(balance == 5))
stopifnot(all(balance == 5))     # sum-to-zero orthogonality relies on balance

# ---------------------------------------------------------------------------
# 2. Filtering, normalisation, voom
#    One gene set throughout; no most-variable-gene selection anywhere.
# ---------------------------------------------------------------------------
design_full <- model.matrix(~ Experiment + Time * Condition, data = targets)
dge  <- DGEList(counts = counts)
keep <- filterByExpr(dge, design = design_full)
dge  <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])     # TMM
msg("Genes retained by filterByExpr: ", sum(keep), " of ", length(keep))

vobj <- voom(dge, design_full, plot = FALSE)

# ---------------------------------------------------------------------------
# 3. canCorPairs — are the factors confounded with one another?
#
#    variancePartition's own diagnostic. Off-diagonal entries are canonical
#    correlations between the metadata variables. Zeros mean the factors are
#    mutually orthogonal, i.e. trial is CROSSED with condition and time point,
#    not confounded with them. This is the single most important thing for a
#    reuser to know about the batch, and it is the first thing the package
#    asks you to check.
# ---------------------------------------------------------------------------
cc <- canCorPairs(~ Experiment + Time + Condition, targets)
write.table(data.frame(Variable = rownames(cc), round(cc, 6), check.names = FALSE),
            file.path(OUT, "table1_canCorPairs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
msg("\n--- canCorPairs: correlation between the experimental factors ---")
print(round(cc, 4))
msg("Maximum off-diagonal correlation: ",
    sprintf("%.2e", max(abs(cc[upper.tri(cc)]))),
    "  (0 = fully crossed, no confounding)")

# plotCorrMatrix draws via ComplexHeatmap/image and needs a roomy device.
pdf(file.path(OUT, "Fig_canCorPairs.pdf"), width = 7, height = 6)
plotCorrMatrix(cc, sort = FALSE)
dev.off()

# ---------------------------------------------------------------------------
# 4. Per-gene variance decomposition — the package's primary product
#
#    Experiment and Time have two levels, Condition three. variancePartition
#    requires categorical variables with few levels to be modelled as FIXED
#    effects: as random effects, lme4 shrinks two-level components to zero
#    (Time comes out at ~1e-31), which is plainly wrong.
#
#    Experiment:Time is retained because it is large; see the summary below.
# ---------------------------------------------------------------------------
form <- ~ Experiment * Time + Condition

vp   <- fitExtractVarPartModel(vobj, form, targets)
vp_s <- sortCols(vp)
vpm  <- as.matrix(vp_s)

# Two summaries of the same per-gene distribution:
#   MedianPct  — the typical gene; unweighted.
#   VarWeightedPct — genes weighted by how much expression variance they carry,
#                    so highly variable genes count more. This is the summary
#                    comparable to a whole-transcriptome figure.
gene_var <- apply(vobj$E, 1, var)
wts      <- gene_var / sum(gene_var)

vp_summary <- data.frame(
  Term            = colnames(vp_s),
  MedianPct       = 100 * apply(vpm, 2, median),
  MeanPct         = 100 * colMeans(vpm),
  VarWeightedPct  = 100 * colSums(vpm * wts),
  Q1Pct           = 100 * apply(vpm, 2, quantile, 0.25),
  Q3Pct           = 100 * apply(vpm, 2, quantile, 0.75),
  row.names = NULL, check.names = FALSE
)
write.table(vp_summary, file.path(OUT, "table2_variancePartition_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(Gene = rownames(vp_s), round(100 * as.data.frame(vp_s), 4),
                       check.names = FALSE),
            file.path(OUT, "table3_variancePartition_per_gene.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
msg("\n--- Per-gene variance explained (%) ---")
print(vp_summary, digits = 4)

p_violin <- plotVarPart(vp_s) +
  ggtitle("Variance explained per gene",
          subtitle = sprintf("%s genes; sum-to-zero contrasts",
                             format(nrow(vp_s), big.mark = ",")))
pdf(file.path(OUT, "Fig_variancePartition_violin.pdf"), width = 6.5, height = 4.6)
print(p_violin)
dev.off()

# ---------------------------------------------------------------------------
# 4b. How much vs how often: a significance test for the same model
#
#     variancePartition answers "how much of the variance does each factor
#     explain", but reports no p-value. The natural companion — not a bolt-on
#     from a different framework — is limma's moderated F-test on the SAME
#     model that variancePartition decomposes: for each term, an omnibus test
#     across the coefficients coding that term, asking whether the factor
#     explains more of a gene's expression than noise. Effect size from
#     variancePartition, evidence from the F-test, one model.
#
#     voom weights are taken from the fuller factorial design (conservative for
#     estimating the mean-variance trend); the linear model itself uses exactly
#     the variancePartition formula, so `assign` maps coefficients to terms.
#
#     Deliberately NOT done: p-values per principal component. The components
#     are ordered by variance, so testing them is selection-affected, and on
#     this dataset it manufactures a nominally significant treatment effect on
#     PC2 (0.49% of that axis, BH-adjusted FDR = 0.047) that the gene-level
#     test flatly contradicts (0 genes). Significance belongs at the gene
#     level; the per-axis table reports effect size only. See
#     docs/ for the full reasoning.
# ---------------------------------------------------------------------------
design_vp <- model.matrix(form, data = targets)
term_of   <- attr(design_vp, "assign")
term_lab  <- attr(terms(form), "term.labels")

fit <- eBayes(lmFit(vobj, design_vp), robust = TRUE)

signif_tbl <- do.call(rbind, lapply(seq_along(term_lab), function(i) {
  tt <- topTable(fit, coef = which(term_of == i), number = Inf, sort.by = "none")
  n  <- sum(tt$adj.P.Val < 0.05)
  data.frame(Term = term_lab[i],
             GenesTested = nrow(tt),
             Genes_FDR05 = n,
             PctGenes_FDR05 = round(100 * n / nrow(tt), 2),
             row.names = NULL, check.names = FALSE)
}))
signif_tbl$MedianPct <- vp_summary$MedianPct[match(signif_tbl$Term, vp_summary$Term)]
write.table(signif_tbl, file.path(OUT, "table5_significance_per_term.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
msg("\n--- Moderated F-test: genes for which each factor beats noise ---")
print(signif_tbl, digits = 4)

# ---------------------------------------------------------------------------
# 5. Per-axis decomposition — the reviewer's minimum requirement
#
#    A PCA supplies the axes and the variance each carries; variancePartition
#    then attributes each axis to the experimental factors, exactly as it does
#    for a gene. Treating the component scores as the features is legitimate
#    because the model and the estimator are unchanged; only the response
#    differs.
# ---------------------------------------------------------------------------
dds     <- DESeqDataSetFromMatrix(counts[keep, , drop = FALSE], targets, ~ 1)
vst_mat <- assay(vst(dds, blind = TRUE))          # blind: unsupervised
stopifnot(identical(colnames(vst_mat), targets$Sample))

pca      <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
var_prop <- pca$sdev^2 / sum(pca$sdev^2)
n_pc     <- 10

pc_mat <- t(pca$x[, seq_len(n_pc)])
rownames(pc_mat) <- colnames(pca$x)[seq_len(n_pc)]

vp_pc <- as.data.frame(fitExtractVarPartModel(pc_mat, form, targets))
vp_pc <- vp_pc[rownames(pc_mat), , drop = FALSE]

pc_tbl <- data.frame(PC = rownames(vp_pc),
                     PctVariance = round(100 * var_prop[seq_len(n_pc)], 2),
                     round(100 * vp_pc, 2),
                     row.names = NULL, check.names = FALSE)
colnames(pc_tbl) <- c("PC", "PctVariance", paste0("R2pct_", colnames(vp_pc)))
write.table(pc_tbl, file.path(OUT, "table4_variance_per_axis.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
msg("\n--- Variance per axis, and attribution of each axis ---")
print(pc_tbl, digits = 4)

# ---------------------------------------------------------------------------
# 5b. Verification that the coding above gives the orthogonal partition
#
#     For a balanced factorial design the classical ANOVA sums of squares are
#     the unique orthogonal decomposition. Under contr.sum, variancePartition
#     must reproduce them. If this check ever fails, the contrast option at the
#     top of the script has been changed or the design is no longer balanced.
# ---------------------------------------------------------------------------
anova_ref <- t(vapply(seq_len(n_pc), function(j) {
  ss <- anova(lm(y ~ Experiment * Time + Condition,
                 data = cbind(targets, y = pca$x[, j])))
  setNames(ss[["Sum Sq"]] / sum(ss[["Sum Sq"]]), rownames(ss))
}, numeric(ncol(vp_pc))))
rownames(anova_ref) <- rownames(vp_pc)
anova_ref <- anova_ref[, colnames(vp_pc), drop = FALSE]

max_dev <- max(abs(100 * (as.matrix(vp_pc) - anova_ref)))
check <- data.frame(PC = rownames(vp_pc),
                    MaxAbsDeviation_pct = round(apply(
                      abs(100 * (as.matrix(vp_pc) - anova_ref)), 1, max), 8),
                    row.names = NULL, check.names = FALSE)
write.table(check, file.path(OUT, "check_orthogonal_partition.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
msg("\nAgreement with the balanced-ANOVA orthogonal partition: max deviation ",
    sprintf("%.2e", max_dev), " percentage points")
stopifnot(max_dev < 1e-6)

# ---------------------------------------------------------------------------
# 6. Figures
# ---------------------------------------------------------------------------
pc_lab <- function(k) sprintf("PC%d (%.1f%%)", k, 100 * var_prop[k])

p_pca <- ggplot(data.frame(targets, PC1 = pca$x[, 1], PC2 = pca$x[, 2]),
                aes(PC1, PC2, colour = Time, shape = Experiment)) +
  geom_point(size = 2.8, alpha = 0.9) +
  facet_wrap(~ Condition) +
  labs(x = pc_lab(1), y = pc_lab(2),
       title = "VST PCA of 60 tomato leaf RNA-seq libraries",
       subtitle = sprintf("%s genes after filterByExpr; blind VST; no most-variable-gene selection",
                          format(nrow(vst_mat), big.mark = ","))) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "Fig_PCA.pdf"), p_pca, width = 9, height = 4.2)

keep_as_is <- c("Experiment", "Time", "Condition", "Experiment:Time", "Residuals")
lev <- c("Residuals", "Condition", "Experiment:Time", "Time", "Experiment")
r2_long <- do.call(rbind, lapply(seq_len(n_pc), function(k) {
  data.frame(PC    = factor(rownames(vp_pc)[k], levels = rownames(vp_pc)),
             Term  = colnames(vp_pc),
             R2pct = 100 * as.numeric(vp_pc[k, ]),
             PCvar = 100 * var_prop[k])
}))
r2_long$Term <- factor(r2_long$Term, levels = lev)

p_axis <- ggplot(r2_long, aes(PC, R2pct, fill = Term)) +
  geom_col(width = 0.72) +
  geom_text(data = unique(r2_long[, c("PC", "PCvar")]),
            aes(PC, 104, label = sprintf("%.1f%%", PCvar)),
            inherit.aes = FALSE, size = 2.9, vjust = 0) +
  scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 25)) +
  labs(y = "Variance of the component explained (%)", x = NULL,
       title = "Attribution of each principal component to the experimental factors",
       subtitle = "variancePartition applied axis by axis. Above each bar: share of total variance carried by that axis") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "Fig_variance_per_axis.pdf"), p_axis, width = 9, height = 5)

# The three ggplot objects behind Figure 2c, 2d and 2e are kept so that
# the 04_figures/ builders can compose the published figure without
# recomputing the decomposition.
saveRDS(list(p_pca = p_pca, p_axis = p_axis, p_violin = p_violin),
        file.path(OUT, "figure2_panels_cd.rds"))

# plotPercentBars: the genes most affected by trial, and by time point
top_batch <- head(order(vpm[, "Experiment"], decreasing = TRUE), 10)
top_time  <- head(order(vpm[, "Time"], decreasing = TRUE), 10)
pdf(file.path(OUT, "Fig_percentBars_top_genes.pdf"), width = 7.5, height = 4.4)
print(plotPercentBars(vp_s[top_batch, ]) +
        ggtitle("Ten genes most affected by experimental trial"))
print(plotPercentBars(vp_s[top_time, ]) +
        ggtitle("Ten genes most affected by time point"))
dev.off()

# ---------------------------------------------------------------------------
# 7. Headline
# ---------------------------------------------------------------------------
batch_terms <- grep("Experiment", vp_summary$Term, value = TRUE)
msg("\n=== HEADLINE ===")
msg("Trial and its interaction with time point, variance-weighted: ",
    sprintf("%.1f%%", sum(vp_summary$VarWeightedPct[vp_summary$Term %in% batch_terms])))
msg("Trial and its interaction with time point, median gene:      ",
    sprintf("%.1f%%", sum(vp_summary$MedianPct[vp_summary$Term %in% batch_terms])))
msg("Treatment (Condition), median gene:                          ",
    sprintf("%.1f%%", vp_summary$MedianPct[vp_summary$Term == "Condition"]))
msg("Genes with a significant treatment effect (FDR<0.05):        ",
    signif_tbl$Genes_FDR05[signif_tbl$Term == "Condition"], " of ",
    signif_tbl$GenesTested[signif_tbl$Term == "Condition"])

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
msg("\nOutputs written to: ", OUT)
