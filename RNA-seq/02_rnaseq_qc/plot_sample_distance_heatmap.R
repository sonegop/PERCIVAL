#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Figure 2b — sample-to-sample distance heatmap.
#
# Blind variance-stabilising transformation (DESeq2) of the raw gene-level
# count matrix, Euclidean distances between the 60 libraries, drawn as a
# symmetric heatmap with rows and columns ordered by complete-linkage
# hierarchical clustering of the distance matrix.
#
# Input  : data/RNAseq_rawcounts.tsv, metadata/targets.txt
# Output : results/rnaseq_qc/sample_distance_heatmap.{pdf,png}
#          results/rnaseq_qc/figure2_panel_b.rds  (distances, for Figure 2)
# Usage  : Rscript 02_rnaseq_qc/plot_sample_distance_heatmap.R
# ---------------------------------------------------------------------------

## Locate config/paths.R by walking up from this script (or the working dir).
.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.d <- if (length(.f)) dirname(normalizePath(sub("^--file=", "", .f[1]))) else normalizePath(getwd())
while (!file.exists(file.path(.d, "config", "paths.R")) && dirname(.d) != .d) .d <- dirname(.d)
if (!file.exists(file.path(.d, "config", "paths.R")))
  stop("Run from inside the repository, or set PERCIVAL_ROOT.", call. = FALSE)
source(file.path(.d, "config", "paths.R"))
source(file.path(.d, "02_rnaseq_qc", "distance_heatmap_functions.R")); rm(.f, .d)

require_inputs(COUNTS_FILE, TARGETS_FILE)
OUT_DIR <- file.path(RES_DIR, "rnaseq_qc")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

d  <- compute_sample_distances(COUNTS_FILE, TARGETS_FILE)
ht <- build_distance_heatmap(d$dist_mat, d$targets)

message(sprintf("%d libraries, %s genes, maximum Euclidean distance %.1f",
                nrow(d$dist_mat), format(d$n_genes, big.mark = ","), max(d$dist_mat)))

draw_ht  <- function() draw(ht, heatmap_legend_side = "right",
                            annotation_legend_side = "right")
size     <- distance_heatmap_size(nrow(d$dist_mat))
pdf_dev  <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

pdf_dev(file.path(OUT_DIR, "sample_distance_heatmap.pdf"),
        width = size[["width"]], height = size[["height"]])
draw_ht(); dev.off()

png(file.path(OUT_DIR, "sample_distance_heatmap.png"),
    width = size[["width"]] * 100, height = size[["height"]] * 100, res = 150)
draw_ht(); dev.off()

## Only the 60 x 60 distance matrix and the sample sheet are kept: the Heatmap
## object itself carries closures that do not serialise.
saveRDS(d[c("dist_mat", "targets")], file.path(OUT_DIR, "figure2_panel_b.rds"))

message("Wrote sample_distance_heatmap.{pdf,png} and figure2_panel_b.rds to ", OUT_DIR)
