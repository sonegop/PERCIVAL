# ---------------------------------------------------------------------------
# Figure 2b — pure builders for the sample-to-sample distance heatmap.
#
# Sourced by 02_rnaseq_qc/plot_sample_distance_heatmap.R (which draws and
# saves it) and by 04_figures/build_figure2.R (which places it in Figure 2).
# Nothing here reads or writes a file except through its arguments.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

## Blind VST of the raw counts, then Euclidean distances between libraries.
## No gene filtering: the panel is a global library-similarity check, so the
## transformation and the distances use all annotated genes.
compute_sample_distances <- function(counts_file, targets_file) {
  targets <- read.delim(targets_file, stringsAsFactors = FALSE)
  targets$Condition  <- factor(targets$Condition, levels = c("Cntr", "E1_D3", "E1_D4"))
  targets$Time       <- factor(targets$Time,      levels = c("T1", "T2"))
  targets$Experiment <- factor(targets$Experiment)

  counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
  counts <- counts[, targets$Sample]              # enforce the sample-sheet order
  stopifnot(identical(colnames(counts), targets$Sample))

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData   = targets,
                                design    = ~ Condition)
  vsd <- vst(dds, blind = TRUE)
  colnames(vsd) <- targets$Sample_name            # readable labels on both margins

  list(dist_mat = as.matrix(dist(t(assay(vsd)))),
       targets  = targets,
       n_genes  = nrow(counts))
}

## Discrete annotation palette for one factor.
annotation_palette <- function(x, palette_name) {
  lvls <- levels(factor(x))
  stats::setNames(
    colorRampPalette(brewer.pal(max(3L, min(length(lvls), 8L)), palette_name))(length(lvls)),
    lvls
  )
}

## Rows and columns are ordered by complete-linkage clustering of the distance
## matrix itself, not of the rows of that matrix.
##
## The arguments after `base_size` exist only so that the same heatmap can be
## drawn into a small page-fit panel (04_figures/build_figure2_pagefit.R), where
## 60 margin labels would fall below 3 pt. Their defaults reproduce the panel as
## published, so callers that do not set them are unaffected.
build_distance_heatmap <- function(dist_mat, targets, base_size = 8,
                                   show_sample_names = TRUE,
                                   anno_size     = unit(5, "mm"),
                                   anno_gap      = unit(1.5, "mm"),
                                   legend_height = unit(3, "cm")) {
  col_fun <- colorRamp2(
    breaks = seq(0, max(dist_mat), length.out = 255),
    colors = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  )

  ann_cols <- list(
    Condition  = annotation_palette(targets$Condition,  "Set2"),
    Experiment = annotation_palette(targets$Experiment, "Set1"),
    Time       = annotation_palette(targets$Time,       "Oranges")
  )
  leg_par <- list(title_gp = gpar(fontsize = base_size + 1),
                  labels_gp = gpar(fontsize = base_size))

  ha_col <- HeatmapAnnotation(
    Condition = targets$Condition, Experiment = targets$Experiment,
    Time = targets$Time, col = ann_cols,
    annotation_legend_param = list(Condition = leg_par, Experiment = leg_par, Time = leg_par),
    annotation_name_gp = gpar(fontsize = base_size), gap = anno_gap,
    simple_anno_size = anno_size
  )
  ha_row <- rowAnnotation(
    Condition = targets$Condition, Experiment = targets$Experiment,
    Time = targets$Time, col = ann_cols,
    show_legend = FALSE, show_annotation_name = FALSE, gap = anno_gap,
    simple_anno_size = anno_size
  )

  Heatmap(
    dist_mat,
    name                        = "Distance",
    col                         = col_fun,
    top_annotation              = ha_col,
    left_annotation             = ha_row,
    clustering_distance_rows    = function(m) as.dist(m),
    clustering_distance_columns = function(m) as.dist(m),
    clustering_method_rows      = "complete",
    clustering_method_columns   = "complete",
    show_row_names              = show_sample_names,
    show_column_names           = show_sample_names,
    show_row_dend               = FALSE,
    show_column_dend            = FALSE,
    row_names_gp                = gpar(fontsize = base_size),
    column_names_gp             = gpar(fontsize = base_size),
    column_names_rot            = 90,
    heatmap_legend_param        = list(
      title         = "Euclidean\ndistance",
      title_gp      = gpar(fontsize = base_size + 1),
      labels_gp     = gpar(fontsize = base_size),
      legend_height = legend_height
    )
  )
}

## Canvas size scales with the number of libraries so margin labels stay legible.
distance_heatmap_size <- function(n) c(width = max(7, n * 0.3) + 3,
                                       height = max(7, n * 0.3) + 1)
