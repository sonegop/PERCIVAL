#!/usr/bin/env Rscript
# ===========================================================================
# Figure 2 — tall alternative layout.
#
#   a  read retention through the pipeline        02_rnaseq_qc/plot_read_retention.R
#   b  sample-to-sample distance heatmap          02_rnaseq_qc/plot_sample_distance_heatmap.R
#   c  PCA of the VST matrix, with the attribution of each principal
#      component to the experimental factors beneath it
#   d  per-gene variance explained                03_batch_effect/variance_partition.R
#
# Kept as an alternative to the figure submitted with the revision, which is
# built by 04_figures/build_figure2_pagefit.R and letters the panels
# differently (violins c, PCA d, per-axis bars e). The letters below are drawn
# into this script's images and do not match the manuscript. Each producing
# script saves its panel object, so this script only
# composes — it recomputes nothing and cannot drift from the panels the
# analysis scripts wrote.
#
# Per-panel titles and subtitles are dropped here: in the assembled figure
# that information belongs in the legend, not repeated inside every panel.
#
# Outputs (results/figures/):
#   Figure2_complete.{pdf,png}      all four panels
#   Figure2_panels_c_d.{pdf,png}    the new material on its own
#
# Usage: Rscript 04_figures/build_figure2.R
#        (run after 02_rnaseq_qc/* and 03_batch_effect/variance_partition.R)
# ===========================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(grid)
})

## Locate config/paths.R by walking up from this script (or the working dir).
.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.d <- if (length(.f)) dirname(normalizePath(sub("^--file=", "", .f[1]))) else normalizePath(getwd())
while (!file.exists(file.path(.d, "config", "paths.R")) && dirname(.d) != .d) .d <- dirname(.d)
if (!file.exists(file.path(.d, "config", "paths.R")))
  stop("Run from inside the repository, or set PERCIVAL_ROOT.", call. = FALSE)
source(file.path(.d, "config", "paths.R"))
source(file.path(.d, "02_rnaseq_qc", "distance_heatmap_functions.R")); rm(.f, .d)

QC_DIR  <- file.path(RES_DIR, "rnaseq_qc")
BE_DIR  <- file.path(RES_DIR, "batch_effect")
OUT_DIR <- file.path(RES_DIR, "figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PANEL_A  <- file.path(QC_DIR, "figure2_panel_a.rds")
PANEL_B  <- file.path(QC_DIR, "figure2_panel_b.rds")
PANEL_CD <- file.path(BE_DIR, "figure2_panels_cd.rds")
require_inputs(PANEL_A, PANEL_B, PANEL_CD)

# ---------------------------------------------------------------------------
# Panel geometry, in inches. The width is that of the submitted Figure 2, so
# panels a and b render at exactly the size already accepted; the heights of
# a and b are likewise unchanged, and the new row is added beneath them.
# ---------------------------------------------------------------------------
FIG_WIDTH   <- 21
H_A         <- 9
H_B         <- 19
H_CD        <- 11
W_LEFT      <- 0.62      # PCA + per-axis bars share the left column of row c/d
H_PCA       <- 0.46      # ... split between the scatter and the bars
LABEL_SIZE  <- 34

# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------
bare <- function(p) p + labs(title = NULL, subtitle = NULL)

p_a  <- bare(readRDS(PANEL_A))
cd   <- readRDS(PANEL_CD)
p_pca    <- bare(cd$p_pca)
p_axis   <- bare(cd$p_axis)
p_violin <- bare(cd$p_violin)

## The Heatmap object carries closures that do not serialise, so it is rebuilt
## from the stored 60 x 60 distance matrix and captured as a grob.
d    <- readRDS(PANEL_B)
ht   <- build_distance_heatmap(d$dist_mat, d$targets)
## grabbing the grob needs an open device; pdf(NULL) keeps it off disk, which
## is what otherwise leaves a stray Rplots.pdf behind.
pdf(NULL)
g_b  <- grid.grabExpr(ComplexHeatmap::draw(ht, heatmap_legend_side = "right",
                                           annotation_legend_side = "right"))
invisible(dev.off())
p_b  <- ggdraw() + draw_grob(g_b)

# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------
col_c <- plot_grid(p_pca, p_axis, ncol = 1, rel_heights = c(H_PCA, 1 - H_PCA))

row_cd <- plot_grid(col_c, p_violin, nrow = 1,
                    rel_widths = c(W_LEFT, 1 - W_LEFT),
                    labels = c("c", "d"), label_size = LABEL_SIZE, vjust = 1.1)

fig2 <- plot_grid(p_a, p_b, row_cd, ncol = 1,
                  rel_heights = c(H_A, H_B, H_CD),
                  labels = c("a", "b", ""), label_size = LABEL_SIZE, vjust = 1.1)

save_both <- function(plot, stem, width, height) {
  ggsave(file.path(OUT_DIR, paste0(stem, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(OUT_DIR, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 150, limitsize = FALSE)
  message(sprintf("  %-24s %4.1f x %4.1f in", stem, width, height))
}

message("Writing to ", OUT_DIR)
save_both(fig2,   "Figure2_complete",   FIG_WIDTH, H_A + H_B + H_CD)

## The new material on its own, for the response letter and for anyone who
## wants the added panels without the two that were already published.
row_cd_alone <- plot_grid(col_c, p_violin, nrow = 1,
                          rel_widths = c(W_LEFT, 1 - W_LEFT),
                          labels = c("c", "d"), label_size = 22, vjust = 1.1)
save_both(row_cd_alone, "Figure2_panels_c_d", 14, 7.5)
