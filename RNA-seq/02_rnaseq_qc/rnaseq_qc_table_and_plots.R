# ============================================================================
# RNA-seq supplementary table + Technical Validation QC plots
# Scientific Data data descriptor - PERCIVAL multi-omics tomato biowaste
#
# Inputs:
#   - data/multiqc_data/multiqc_data.json
#   - data/RNAseq_rawcounts.tsv
#   - metadata/targets.txt
#
# Outputs (results/rnaseq_qc/):
#   Supplementary_Table_RNAseq_QC.tsv         - sample-level QC stats
#   Supplementary_Table_RNAseq_QC.xlsx        - same as XLSX
#   figures_QC/Fig_2A_read_retention.pdf      - stacked bar read retention
#   figures_QC/Fig_2B_mapping_rate_boxplot.pdf
#   figures_QC/Fig_2C_replicate_correlation_heatmap.pdf
#   figures_QC/Fig_2D_PCA.pdf
#   figures_QC/Fig_2E_marker_gene_heatmap.pdf
#   figures_QC/Fig_2F_gene_detection_saturation.pdf
#
# Format follows Supplementary Table 1 of Zhao et al., Sci Data 2024
# (41597_2024_3183 - Spermophilus alashanicus): columns Sample / Raw reads /
# Clean reads / Error rate(%) / Q20(%) / Q30(%) / GC content(%) / Dup_pair /
# rRNA ratio(%). rRNA ratio is reported as N.A. - no SortMeRNA / RSeQC step
# was run in the pipeline; rRNA depletion was performed at library prep
# (poly-A selection with KAPA Stranded mRNA-Seq Kit).
# ============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(ComplexHeatmap)
  library(matrixStats)
  library(openxlsx)
  library(ggrepel)
  library(circlize)
})

## Locate config/paths.R by walking up from this script (or the working dir).
.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.d <- if (length(.f)) dirname(normalizePath(sub("^--file=", "", .f[1]))) else normalizePath(getwd())
while (!file.exists(file.path(.d, "config", "paths.R")) && dirname(.d) != .d) .d <- dirname(.d)
if (!file.exists(file.path(.d, "config", "paths.R")))
  stop("Run from inside the repository, or set PERCIVAL_ROOT.", call. = FALSE)
source(file.path(.d, "config", "paths.R")); rm(.f, .d)

# ---- Paths (all resolved through config/paths.R) -------------------------- #
MQC_JSON <- file.path(DATA_DIR, "multiqc_data", "multiqc_data.json")
TARGETS  <- TARGETS_FILE                     # metadata/targets.txt
                                             # COUNTS_FILE comes from paths.R
OUT_DIR  <- file.path(RES_DIR, "rnaseq_qc")
FIG_DIR  <- file.path(OUT_DIR, "figures_QC")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
require_inputs(MQC_JSON, TARGETS, COUNTS_FILE)

# ---- 1. Load metadata ---------------------------------------------------- #
meta <- read_tsv(TARGETS, show_col_types = FALSE) %>%
  mutate(Group = paste(Condition, Time, Experiment, sep = "_"))

stopifnot(nrow(meta) == 60)

# Helper to strip the long sequencing-lane suffix from MultiQC sample IDs
# Handles fastp / fastqc / featurecounts (`_MKDL...`) and STAR (`___STARpass1`)
# and samtools (`__stats`) suffixes.
short_id <- function(x) {
  x <- sub("___STARpass[0-9]+$", "", x)
  x <- sub("__stats$", "", x)
  x <- sub("_MKDL[0-9].*$", "", x)
  x
}

# ---- 2. Parse MultiQC JSON ----------------------------------------------- #
mqc <- fromJSON(MQC_JSON, simplifyDataFrame = FALSE)

# fastp - raw / clean reads, Q20, Q30, GC, duplication
fastp_records <- mqc$report_saved_raw_data$multiqc_fastp
fastp_df <- bind_rows(lapply(names(fastp_records), function(s) {
  r <- fastp_records[[s]]
  tibble(
    fastp_sample = s,
    raw_reads    = r$summary$before_filtering$total_reads,
    clean_reads  = r$summary$after_filtering$total_reads,
    q20_pct      = round(100 * r$summary$after_filtering$q20_rate, 2),
    q30_pct      = round(100 * r$summary$after_filtering$q30_rate, 2),
    gc_pct       = round(100 * r$summary$after_filtering$gc_content, 2),
    dup_pair     = round(r$duplication$rate, 4)
  )
})) %>%
  mutate(Sample = short_id(fastp_sample))

`%||%` <- function(a, b) if (is.null(a)) b else a

# STAR - total/uniquely mapped/multimapped + per-base mismatch rate
# (used as 'Error rate(%)' in the supplementary table because the
# samtools_stats error_rate field was 0.0 across all samples in the
# MultiQC export; STAR mismatch_rate is the equivalent quantity reported
# in % directly.)
star_records <- mqc$report_saved_raw_data$multiqc_star
star_df <- bind_rows(lapply(names(star_records), function(s) {
  r <- star_records[[s]]
  tibble(
    star_sample = s,
    star_total_reads        = r$total_reads,
    star_uniquely_mapped    = r$uniquely_mapped,
    star_uniquely_mapped_pct= r$uniquely_mapped_percent,
    star_multimapped        = r$multimapped,
    star_multimapped_pct    = r$multimapped_percent,
    error_rate_pct          = round(r$mismatch_rate %||% NA_real_, 4)
  )
})) %>%
  mutate(Sample = short_id(star_sample))

# featureCounts - assigned + percent_assigned
fc_records <- mqc$report_saved_raw_data$multiqc_featurecounts
fc_df <- bind_rows(lapply(names(fc_records), function(s) {
  r <- fc_records[[s]]
  tibble(
    fc_sample          = s,
    fc_assigned        = r$Assigned,
    fc_percent_assigned= r$percent_assigned
  )
})) %>%
  mutate(Sample = short_id(fc_sample))

# ---- 3. Build supplementary QC table ------------------------------------- #
qc_tbl <- meta %>%
  left_join(fastp_df, by = "Sample") %>%
  left_join(star_df,  by = "Sample") %>%
  left_join(fc_df,    by = "Sample") %>%
  transmute(
    Sample            = Sample,
    Experiment        = Experiment,
    Condition         = Condition,
    Time              = Time,
    `Raw reads`       = raw_reads,
    `Clean reads`     = clean_reads,
    `Error rate(%)`   = error_rate_pct,
    `Q20(%)`          = q20_pct,
    `Q30(%)`          = q30_pct,
    `GC content(%)`   = gc_pct,
    `Dup_pair`        = dup_pair,
    `rRNA ratio(%)`   = "N.A."
  )

# Order by sample number (PERCIVAL_n_m)
qc_tbl <- qc_tbl %>%
  mutate(.n = as.integer(sub("^PERCIVAL_([0-9]+)_.*", "\\1", Sample)),
         .e = as.integer(sub("^PERCIVAL_[0-9]+_([12])$", "\\1", Sample))) %>%
  arrange(.e, .n) %>%
  select(-.n, -.e)

write_tsv(qc_tbl, file.path(OUT_DIR, "Supplementary_Table_RNAseq_QC.tsv"))

wb <- createWorkbook()
addWorksheet(wb, "RNAseq_QC")
writeData(wb, "RNAseq_QC", qc_tbl, headerStyle = createStyle(textDecoration = "bold"))
setColWidths(wb, "RNAseq_QC", cols = 1:ncol(qc_tbl), widths = "auto")
saveWorkbook(wb, file.path(OUT_DIR, "Supplementary_Table_RNAseq_QC.xlsx"), overwrite = TRUE)

cat("Supplementary Table written:\n  ",
    file.path(OUT_DIR, "Supplementary_Table_RNAseq_QC.tsv"), "\n  ",
    file.path(OUT_DIR, "Supplementary_Table_RNAseq_QC.xlsx"), "\n")

# ---- 4. Read retention table for Fig 2A ---------------------------------- #
retention <- meta %>%
  left_join(fastp_df %>% select(Sample, raw_reads, clean_reads), by = "Sample") %>%
  left_join(star_df  %>% select(Sample, star_uniquely_mapped),   by = "Sample") %>%
  left_join(fc_df    %>% select(Sample, fc_assigned),            by = "Sample") %>%
  mutate(
    raw_reads_pe          = raw_reads / 2,
    clean_reads_pe        = clean_reads / 2,
    `Trimmed-out`         = raw_reads_pe - clean_reads_pe,
    `Multi/unmapped`      = clean_reads_pe - star_uniquely_mapped,
    `Unique not assigned` = star_uniquely_mapped - fc_assigned,
    `Assigned to feature` = fc_assigned
  ) %>%
  select(Sample, Experiment, Condition, Time,
         `Trimmed-out`, `Multi/unmapped`,
         `Unique not assigned`, `Assigned to feature`) %>%
  pivot_longer(cols = c(`Trimmed-out`, `Multi/unmapped`,
                        `Unique not assigned`, `Assigned to feature`),
               names_to = "Stage", values_to = "Reads") %>%
  mutate(Stage = factor(Stage,
           levels = c("Trimmed-out","Multi/unmapped",
                      "Unique not assigned","Assigned to feature")),
         Group_label = paste(Condition, Time, sep="_"))

p2A <- ggplot(retention,
              aes(x = reorder(Sample, as.integer(sub("PERCIVAL_([0-9]+)_.*","\\1", Sample))),
                  y = Reads / 1e6, fill = Stage)) +
  geom_col(width = 0.85) +
  facet_grid(. ~ Experiment, scales = "free_x", space = "free_x",
             labeller = labeller(Experiment = function(x) paste0("Experiment ", x))) +
  scale_fill_manual(values = c("Trimmed-out"        = "#9ecae1",
                               "Multi/unmapped"     = "#fdae6b",
                               "Unique not assigned"= "#bdbdbd",
                               "Assigned to feature"= "#31a354")) +
  labs(x = NULL, y = "Read pairs (millions)",
       title = "Figure 2A - Read retention through the RNA-seq pipeline",
       subtitle = "raw -> fastp -> STAR uniquely mapped -> featureCounts assigned",
       fill = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey92"))

ggsave(file.path(FIG_DIR, "Fig_2A_read_retention.pdf"), p2A,
       width = 11, height = 5)

# ---- 5. Mapping-rate boxplot (Fig 2B) ------------------------------------ #
# geom_jitter draws from the RNG, so the seed is fixed here: without it the
# point positions, and therefore the output file, change on every run.
set.seed(42)
map_df <- meta %>%
  left_join(star_df %>% select(Sample, star_uniquely_mapped_pct), by = "Sample") %>%
  mutate(Treatment = factor(Condition, levels = c("Cntr","E1_D3","E1_D4")),
         Time      = factor(Time,      levels = c("T1","T2")),
         Experiment= factor(Experiment))

p2B <- ggplot(map_df, aes(x = Treatment, y = star_uniquely_mapped_pct,
                          fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.85) +
  facet_grid(Experiment ~ Time, labeller = label_both) +
  geom_hline(yintercept = 85, linetype = "dashed", colour = "red") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "STAR uniquely-mapped reads (%)",
       title = "Figure 2B - Per-sample uniquely-mapped read fraction",
       subtitle = "dashed line = 85% acceptance threshold") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "Fig_2B_mapping_rate_boxplot.pdf"), p2B,
       width = 7, height = 5)

# ---- 6. Counts-based analyses (correlation, PCA, markers, saturation) ----- #
counts <- read.table(COUNTS_FILE, header = TRUE, sep = "\t",
                     row.names = 1, check.names = FALSE)
stopifnot(all(meta$Sample %in% colnames(counts)))
counts <- as.matrix(counts[, meta$Sample])

# Filter very low-expressed genes for stability
keep <- rowSums(edgeR::cpm(counts) >= 1) >= 5
counts_f <- counts[keep, ]

# DESeq2 VST for PCA + correlation
dds <- DESeqDataSetFromMatrix(countData = counts_f,
                              colData   = meta,
                              design    = ~ 1)
vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)

# ---- 6a. Replicate Pearson correlation heatmap (Fig 2C) ------------------ #
cor_mat <- cor(vst_mat, method = "pearson")

ann_df <- meta %>%
  transmute(Experiment = factor(Experiment),
            Time       = factor(Time),
            Treatment  = factor(Condition, levels = c("Cntr","E1_D3","E1_D4"))) %>%
  as.data.frame()
rownames(ann_df) <- meta$Sample

ann <- HeatmapAnnotation(df = ann_df,
                         col = list(
                           Experiment = c(`1` = "#1b9e77", `2` = "#d95f02"),
                           Time       = c(T1 = "#7570b3", T2 = "#e7298a"),
                           Treatment  = c(Cntr = "#66c2a5", E1_D3 = "#fc8d62", E1_D4 = "#8da0cb")
                         ))

pdf(file.path(FIG_DIR, "Fig_2C_replicate_correlation_heatmap.pdf"),
    width = 11, height = 9)
draw(Heatmap(cor_mat,
             name = "Pearson r",
             col = colorRamp2(c(0.85, 0.95, 1.0),
                              c("#2c7bb6", "#ffffbf", "#d7191c")),
             top_annotation = ann,
             cluster_rows = TRUE, cluster_columns = TRUE,
             show_row_names = TRUE, show_column_names = TRUE,
             row_names_gp = gpar(fontsize = 6),
             column_names_gp = gpar(fontsize = 6),
             column_title = "Figure 2C - Pearson correlation (VST counts) of all 60 RNA-seq libraries"))
dev.off()

# ---- 6b. PCA (Fig 2D) ---------------------------------------------------- #
top_var <- order(rowVars(vst_mat), decreasing = TRUE)[1:2000]
pca <- prcomp(t(vst_mat[top_var, ]), center = TRUE, scale. = FALSE)
ve  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_df <- as_tibble(pca$x[, 1:3]) %>%
  mutate(Sample = colnames(vst_mat)) %>%
  left_join(meta, by = "Sample") %>%
  mutate(Treatment = factor(Condition, levels = c("Cntr","E1_D3","E1_D4")),
         Experiment= factor(Experiment),
         Time      = factor(Time))

p2D <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Treatment, shape = Time)) +
  geom_point(size = 3, alpha = 0.85) +
  facet_wrap(~ Experiment, labeller = labeller(Experiment = function(x) paste0("Experiment ", x))) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Figure 2D - PCA on top-2000 variable genes (VST)",
       x = sprintf("PC1 (%.1f%%)", ve[1]),
       y = sprintf("PC2 (%.1f%%)", ve[2])) +
  theme_bw(base_size = 10)

ggsave(file.path(FIG_DIR, "Fig_2D_PCA.pdf"), p2D, width = 9, height = 5)

# ---- 6c. Marker-gene heatmap (Fig 2E) ------------------------------------ #
# A small panel of housekeeping + photosynthesis + nitrogen-assimilation genes
# (ITAG5.0 IDs). If a gene is not present in the count matrix it is skipped.
marker_table <- tibble::tribble(
  ~id,                         ~symbol,
  "Solyc06T000004.1.ITAG5.0",  "EF1a",
  "Solyc01T000516.1.ITAG5.0",  "ACTIN",
  "Solyc10T000146.1.ITAG5.0",  "UBQ",
  "Solyc04T000300.1.ITAG5.0",  "GAPDH",
  "Solyc02T000799.1.ITAG5.0",  "RBCS",
  "Solyc03T000076.1.ITAG5.0",  "LHCB",
  "Solyc01T002336.1.ITAG5.0",  "NIA",
  "Solyc01T002190.1.ITAG5.0",  "GS1",
  "Solyc03T000684.1.ITAG5.0",  "PAL",
  "Solyc05T002477.1.ITAG5.0",  "CHS"
)
marker_table <- marker_table[marker_table$id %in% rownames(vst_mat), ]
markers <- marker_table$id

if (length(markers) >= 3) {
  m <- vst_mat[markers, , drop = FALSE]
  rownames(m) <- paste0(marker_table$symbol, " (", marker_table$id, ")")
  m <- t(scale(t(m)))
  cairo_pdf(file.path(FIG_DIR, "Fig_2E_marker_gene_heatmap.pdf"),
            width = 12, height = max(3.5, 0.45 * length(markers) + 2))
  draw(Heatmap(m, name = "z-score (VST)",
               top_annotation = ann,
               cluster_rows = TRUE, cluster_columns = TRUE,
               row_names_gp = gpar(fontsize = 9),
               column_names_gp = gpar(fontsize = 6),
               column_title = paste0("Figure 2E - Marker-gene expression across the 60 RNA-seq libraries (",
                                     length(markers), " markers)")))
  dev.off()
} else {
  message("Marker-gene heatmap skipped - fewer than 3 markers matched ITAG5.0 IDs.")
}

# ---- 6d. Gene-detection saturation (Fig 2F) ------------------------------ #
# Subsample each library at increasing depths and count detected genes (>=1 read)
set.seed(42)
depths_pct <- c(5, 10, 25, 50, 75, 100)
sat <- bind_rows(lapply(meta$Sample, function(s) {
  v <- counts[, s]
  total <- sum(v)
  if (total == 0) return(NULL)
  do.call(rbind, lapply(depths_pct, function(p) {
    n_keep <- round(total * p / 100)
    if (n_keep < 100) return(NULL)
    sub <- rmultinom(1, n_keep, prob = v / total)[, 1]
    data.frame(Sample = s, Depth_pct = p,
               Detected = sum(sub >= 1), stringsAsFactors = FALSE)
  }))
})) %>%
  left_join(meta, by = "Sample")

p2F <- ggplot(sat, aes(x = Depth_pct, y = Detected, group = Sample,
                       colour = Condition)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 0.8, alpha = 0.7) +
  facet_wrap(~ Experiment,
             labeller = labeller(Experiment = function(x) paste0("Experiment ", x))) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Figure 2F - Gene-detection saturation",
       subtitle = "Detected genes (>=1 count) at increasing fractions of library depth",
       x = "Library depth (%)", y = "Genes detected") +
  theme_bw(base_size = 10)

ggsave(file.path(FIG_DIR, "Fig_2F_gene_detection_saturation.pdf"), p2F,
       width = 9, height = 5)

cat("All Technical-Validation figures written to:\n  ", FIG_DIR, "\n")
