#!/usr/bin/env Rscript
# ===========================================================================
# Figure 2 — page-fit alternative layout.
#
# Same five panels as 04_figures/build_figure2.R, re-laid out to fit a single
# journal page at the Springer Nature double-column width (183 mm max). The
# device is 180 x 202 mm; every panel is sized to be legible at that width
# rather than at the 21 x 39 in canvas the tall version uses.
#
#   a          full width     read retention through the pipeline
#   b | c      65% | 35%      sample-distance heatmap | per-gene variance
#   d | e      55% | 45%      PCA of the VST matrix  | variance per component
#
# Panel letters differ from the tall layout, where the PCA and the per-axis
# bars share panel c and the violins are panel d:
#
#   tall layout        page-fit layout
#   c (upper)  PCA           -> d
#   c (lower)  per-axis bars -> e
#   d          violins       -> c
#
# Prose that cites panel letters (Technical Validation, response letter) is
# written against the tall layout and would need remapping if this version is
# adopted for the manuscript.
#
# This script composes only: each panel object is read from the .rds its
# producing script wrote, so nothing is recomputed and the two layouts cannot
# disagree about the numbers.
#
# Outputs (results/figures/):
#   Figure2_complete_pagefit.pdf   vector, 180 x 202 mm
#   Figure2_complete_pagefit.png   600 dpi, 4252 x 4773 px
#
# Usage: Rscript 04_figures/build_figure2_pagefit.R
#        (run after 02_rnaseq_qc/* and 03_batch_effect/variance_partition.R)
# ===========================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
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
# Geometry, in millimetres. BASE is the body font size in points at final
# size; Springer Nature accepts 5-7 pt, so 7 pt leaves one point of headroom
# for a production editor who scales the figure down.
# ---------------------------------------------------------------------------
FIG_W  <- 180
FIG_H  <- 202
BASE   <- 7
TAG_PT <- 11

ROW_H  <- c(0.30, 0.52, 0.18)   # a / b+c / d+e
W_B    <- 13 / 20               # heatmap share of the middle row
W_D    <- 11 / 20               # PCA share of the bottom row

## The heatmap grob keeps absolute units once captured, so its legend and
## annotation sizes are tuned against the cell it lands in, not against a
## default 7-inch device.
CELL_W <- FIG_W * W_B
CELL_H <- FIG_H * ROW_H[2]

# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------
## In an assembled figure the per-panel titles belong in the legend.
bare <- function(p) p + labs(title = NULL, subtitle = NULL)

## Size deltas only, never a complete theme: plotVarPart() sets the rotated
## axis labels of panel c through the theme, and a complete theme would
## silently drop them.
shrink <- function(p, base = BASE) {
  p + theme(
    text            = element_text(size = base),
    axis.text       = element_text(size = base - 1),
    axis.title      = element_text(size = base),
    strip.text      = element_text(size = base - 1,
                                   margin = margin(1.2, 1.2, 1.2, 1.2, "pt")),
    legend.text     = element_text(size = base - 1),
    legend.title    = element_text(size = base, face = "plain"),
    legend.key.size = unit(3, "mm"),
    legend.margin   = margin(0, 0, 0, 0),
    legend.box.spacing = unit(2, "pt"),
    plot.margin     = margin(2, 2, 2, 2, "pt")
  )
}

cd <- readRDS(PANEL_CD)

## a — the 60 rotated library names would fall to ~3 pt at 180 mm, so they are
## dropped here. The bar-to-library mapping is in the supplementary QC table.
p_a <- shrink(bare(readRDS(PANEL_A))) +
  theme(axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        axis.title.x    = element_blank(),
        legend.position = "bottom",
        legend.title    = element_blank())

## c — five term names on a 63 mm axis, so they are angled rather than shrunk
## further.
p_c <- shrink(bare(cd$p_violin)) +
  theme(axis.text.x = element_text(size = BASE - 2, angle = 35, hjust = 1))

p_d <- shrink(bare(cd$p_pca)) +
  theme(legend.position = "bottom", legend.box = "horizontal")

## e — the axis title and the five term names have to fit 81 mm. The per-axis
## percentages above the bars are a geom_text layer, whose size is in mm and
## therefore does not follow the theme: it is set on the layer itself.
p_e <- shrink(bare(cd$p_axis)) +
  labs(y = "Variance explained (%)") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = BASE - 2),
        legend.key.size = unit(2.4, "mm"),
        axis.text.x     = element_text(size = BASE - 2, angle = 90,
                                       hjust = 1, vjust = 0.5))
p_e$layers[[2]]$aes_params$size <- 1.5   # mm; was 2.9 for the standalone figure

## b — rebuilt from the stored 60 x 60 distance matrix: the Heatmap object
## carries closures that do not serialise. Both margin label sets are dropped
## for the same reason as in panel a; group identity is carried by the
## annotation tracks.
d  <- readRDS(PANEL_B)
ht <- build_distance_heatmap(d$dist_mat, d$targets,
                             base_size         = 6,
                             show_sample_names = FALSE,
                             anno_size         = unit(2, "mm"),
                             anno_gap          = unit(0.6, "mm"),
                             legend_height     = unit(25, "mm"))

## grid.grabExpr needs an open device, and ComplexHeatmap makes layout
## decisions against its size, so the device matches the destination cell.
## pdf(NULL) keeps the capture off disk.
pdf(NULL, width = CELL_W / 25.4, height = CELL_H / 25.4)
g_b <- grid.grabExpr(
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right",
                       annotation_legend_side = "right",
                       padding = unit(c(1, 1, 1, 1), "mm")))
invisible(dev.off())
p_b <- wrap_elements(full = g_b)

# ---------------------------------------------------------------------------
# Assembly. patchwork binds the design letters to the plots in the order they
# are added, so the order below must stay aligned with the design string.
# ---------------------------------------------------------------------------
design <- "
AAAAAAAAAAAAAAAAAAAA
BBBBBBBBBBBBBCCCCCCC
DDDDDDDDDDDEEEEEEEEE
"

fig2 <- p_a + p_b + p_c + p_d + p_e +
  plot_layout(design = design, heights = ROW_H) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.margin = margin(2, 2, 2, 2, "pt"))
  ) &
  theme(plot.tag = element_text(size = TAG_PT, face = "bold", hjust = 0, vjust = 1))

message("Writing to ", OUT_DIR)
ggsave(file.path(OUT_DIR, "Figure2_complete_pagefit.pdf"), fig2,
       width = FIG_W, height = FIG_H, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Figure2_complete_pagefit.png"), fig2,
       width = FIG_W, height = FIG_H, units = "mm", dpi = 600)
message(sprintf("  Figure2_complete_pagefit  %d x %d mm", FIG_W, FIG_H))
