#!/usr/bin/env Rscript
# ===========================================================================
# Supplementary Table 6 — variance decomposition of the RNA-seq dataset.
#
# Assembles the five tables written by 03_batch_effect/variance_partition.R
# into one workbook, one sheet per table plus a legend. This script composes
# only: it reads the deposited TSVs and recomputes nothing, so the workbook
# cannot disagree with the analysis behind it.
#
# The model term names are relabelled for readers ("Experiment" -> "Trial",
# "Condition" -> "Treatment"); the legend sheet records the mapping so the
# workbook can still be matched against the code and the TSVs.
#
# Output: results/batch_effect/Supplementary_Table_6_variance_decomposition.xlsx
#
# Usage: Rscript 03_batch_effect/build_supplementary_table.R
#        (run after 03_batch_effect/variance_partition.R)
# ===========================================================================

suppressPackageStartupMessages(library(openxlsx))

## Locate config/paths.R by walking up from this script (or the working dir).
.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.d <- if (length(.f)) dirname(normalizePath(sub("^--file=", "", .f[1]))) else normalizePath(getwd())
while (!file.exists(file.path(.d, "config", "paths.R")) && dirname(.d) != .d) .d <- dirname(.d)
if (!file.exists(file.path(.d, "config", "paths.R")))
  stop("Run from inside the repository, or set PERCIVAL_ROOT.", call. = FALSE)
source(file.path(.d, "config", "paths.R")); rm(.f, .d)

IN  <- file.path(RES_DIR, "batch_effect")
OUT <- file.path(IN, "Supplementary_Table_6_variance_decomposition.xlsx")

SRC <- c(canCor      = "table1_canCorPairs.tsv",
         summary     = "table2_variancePartition_summary.tsv",
         per_gene    = "table3_variancePartition_per_gene.tsv",
         per_axis    = "table4_variance_per_axis.tsv",
         significance = "table5_significance_per_term.tsv")
require_inputs(file.path(IN, SRC))

read_tbl <- function(f) read.delim(file.path(IN, f), check.names = FALSE,
                                   stringsAsFactors = FALSE)
tbl <- lapply(SRC, read_tbl)

# ---------------------------------------------------------------------------
# Reader-facing labels. Everything downstream goes through these two maps, so
# the code names and the published names can never drift apart silently.
# ---------------------------------------------------------------------------
TERM <- c(Experiment        = "Trial",
          Time              = "Time point",
          Condition         = "Treatment",
          "Experiment:Time" = "Trial x time point",
          Residuals         = "Residual (between plants)")

relabel <- function(x) ifelse(x %in% names(TERM), TERM[x], x)

rename_cols <- function(df, prefix = "") {
  nm <- sub(paste0("^", prefix), "", names(df))
  names(df) <- ifelse(nm %in% names(TERM), TERM[nm], names(df))
  df
}

## 1 — canonical correlations between the three design factors
s_cancor <- tbl$canCor
s_cancor$Variable <- relabel(s_cancor$Variable)
s_cancor <- rename_cols(s_cancor)

## 2 — per-principal-component attribution
s_axis <- rename_cols(tbl$per_axis, prefix = "R2pct_")
names(s_axis)[names(s_axis) == "PctVariance"] <- "Total variance carried by the axis (%)"

## The summary and significance tables are written at full precision by
## variance_partition.R; two decimals is what a reader needs here.
round_numeric <- function(df, digits = 2) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], round, digits)
  df
}

## 3 — per-gene summary across all genes
s_summary <- round_numeric(tbl$summary)
s_summary$Term <- relabel(s_summary$Term)
names(s_summary) <- c("Term", "Median (%)", "Mean (%)",
                      "Variance-weighted mean (%)", "1st quartile (%)",
                      "3rd quartile (%)")

## 4 — how often, not just how much
s_signif <- round_numeric(tbl$significance)
s_signif$Term <- relabel(s_signif$Term)
names(s_signif) <- c("Term", "Genes tested", "Genes at FDR < 0.05",
                     "Percent of genes tested", "Median variance explained (%)")

## 5 — the full per-gene decomposition
s_genes <- rename_cols(tbl$per_gene)

legend <- data.frame(
  Sheet = c("1_canonical_correlations", "2_per_axis", "3_per_gene_summary",
            "4_significance_per_term", "5_per_gene"),
  Contents = c(
    "Canonical correlations between trial, time point and treatment. Values of zero (to numerical precision) show that the three factors are mutually orthogonal, so the trial is crossed with the design rather than confounded with it.",
    "Share of each of the first ten principal components attributed to each model term, and the share of total transcriptome variance carried by that component. Computed on the variance-stabilised matrix.",
    "Distribution of the per-gene variance explained, summarised across all genes tested. The variance-weighted mean weights each gene by the expression variance it carries.",
    "Number of genes with a significant effect for each term, from a moderated F-test (limma) on the same model, Benjamini-Hochberg FDR < 0.05.",
    "Per-gene variance explained (%) for every gene tested, one row per gene. Rows sum to 100%."
  ),
  stringsAsFactors = FALSE
)

notes <- data.frame(
  Item = c("Dataset", "Genes tested", "Model", "Contrast coding", "Software",
           "Term names", "Source", "Code"),
  Value = c(
    "60 tomato leaf RNA-seq libraries: 3 treatments x 2 time points x 5 biological replicates x 2 independent trials.",
    sprintf("%s of %s annotated genes, retained by edgeR::filterByExpr. No selection of most-variable genes.",
            format(nrow(s_genes), big.mark = ","), "43,752"),
    "~ trial * time point + treatment",
    "Sum-to-zero (contr.sum). The design is balanced, so this yields the unique orthogonal partition; treatment contrasts would split the main effects and their interaction incorrectly.",
    paste0("variancePartition ", tryCatch(as.character(utils::packageVersion("variancePartition")), error = function(e) "1.42.0"),
           ", limma, DESeq2. Full session details in results/batch_effect/sessionInfo.txt."),
    "Trial = Experiment, Treatment = Condition in the deposited TSVs and in the analysis code.",
    paste(SRC, collapse = ", "),
    "https://github.com/sonegop/PERCIVAL (10.5281/zenodo.22232657)"
  ),
  stringsAsFactors = FALSE
)

wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", valign = "top")
add <- function(name, df, widths = "auto", wrap = FALSE) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hdr)
  if (wrap)
    addStyle(wb, name, createStyle(wrapText = TRUE, valign = "top"),
             rows = 2:(nrow(df) + 1), cols = 2, gridExpand = TRUE)
  setColWidths(wb, name, cols = seq_along(df), widths = widths)
  freezePane(wb, name, firstRow = TRUE)
}

addWorksheet(wb, "0_legend")
writeData(wb, "0_legend", "Supplementary Table 6. Variance decomposition of the tomato leaf RNA-seq dataset.",
          startRow = 1)
addStyle(wb, "0_legend", createStyle(textDecoration = "bold"), rows = 1, cols = 1)
writeData(wb, "0_legend", legend, startRow = 3, headerStyle = hdr)
writeData(wb, "0_legend", notes, startRow = nrow(legend) + 5, headerStyle = hdr)
addStyle(wb, "0_legend", createStyle(wrapText = TRUE, valign = "top"),
         rows = 4:(nrow(legend) + nrow(notes) + 6), cols = 2, gridExpand = TRUE)
setColWidths(wb, "0_legend", cols = 1:2, widths = c(26, 110))

add("1_canonical_correlations", s_cancor)
add("2_per_axis",               s_axis)
add("3_per_gene_summary",       s_summary)
add("4_significance_per_term",  s_signif)
add("5_per_gene",               s_genes)

saveWorkbook(wb, OUT, overwrite = TRUE)
message("Wrote ", OUT)
message(sprintf("  %d sheets; %s gene rows", length(names(wb)),
                format(nrow(s_genes), big.mark = ",")))
