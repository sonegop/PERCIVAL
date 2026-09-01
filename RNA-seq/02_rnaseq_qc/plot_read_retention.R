#!/usr/bin/env Rscript
## Read retention through the RNA-seq pipeline
## raw -> fastp -> STAR uniquely mapped -> featureCounts assigned

library(tidyverse)

## Locate config/paths.R by walking up from this script (or the working dir).
.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.d <- if (length(.f)) dirname(normalizePath(sub("^--file=", "", .f[1]))) else normalizePath(getwd())
while (!file.exists(file.path(.d, "config", "paths.R")) && dirname(.d) != .d) .d <- dirname(.d)
if (!file.exists(file.path(.d, "config", "paths.R")))
  stop("Run from inside the repository, or set PERCIVAL_ROOT.", call. = FALSE)
source(file.path(.d, "config", "paths.R")); rm(.f, .d)

IN_FILE <- file.path(DATA_DIR, "read_retention_merged.tsv")
OUT_DIR <- file.path(RES_DIR, "rnaseq_qc")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
require_inputs(IN_FILE)
df <- read_tsv(IN_FILE, show_col_types = FALSE)

## Order samples within each experiment by their numeric suffix (PERCIVAL_<n>_<exp>)
df <- df %>%
  mutate(
    sample_num = as.integer(str_extract(Sample, "(?<=_)[0-9]+(?=_)")),
    Experiment_lab = paste("Experiment", Experiment)
  ) %>%
  arrange(Experiment, sample_num) %>%
  mutate(Sample_name = factor(Sample_name, levels = unique(Sample_name)))

## Build the four stacked categories from the cumulative retention columns
plot_df <- df %>%
  transmute(
    Sample_name,
    Experiment_lab,
    `Assigned to feature`   = Successfully.assigned.alignments...featureCounts.v.2.1.1.,
    `Unique not assigned`   = Uniquely.mapped.reads.number..STAR.v2.7.11b. - Successfully.assigned.alignments...featureCounts.v.2.1.1.,
    `Multi/unmapped`        = Total.sequences.after.preprocessing..fastp.v0.23.4. - Uniquely.mapped.reads.number..STAR.v2.7.11b.,
    `Trimmed-out`           = Total.Sequences - Total.sequences.after.preprocessing..fastp.v0.23.4.
  ) %>%
  pivot_longer(
    cols = c(`Assigned to feature`, `Unique not assigned`, `Multi/unmapped`, `Trimmed-out`),
    names_to = "Category",
    values_to = "Reads"
  ) %>%
  mutate(
    Category = factor(Category,
                       levels = c("Trimmed-out", "Multi/unmapped", "Unique not assigned", "Assigned to feature")),
    Reads_millions = Reads / 1e6
  )

pal <- c(
  "Trimmed-out"          = "#A6D3E8",
  "Multi/unmapped"       = "#F5A85A",
  "Unique not assigned"  = "#BDBDBD",
  "Assigned to feature"  = "#2E9E4F"
)

p <- ggplot(plot_df, aes(x = Sample_name, y = Reads_millions, fill = Category)) +
  geom_col(width = 0.75) +
  facet_grid(~ Experiment_lab, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal, name = NULL) +
  labs(
    title = "Read retention through the RNA-seq pipeline",
    subtitle = "raw -> fastp -> STAR uniquely mapped -> featureCounts assigned",
    x = NULL,
    y = "Read pairs (millions)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "plain", size = 11),
    legend.position = "bottom",
    plot.title = element_text(face = "plain", size = 16),
    plot.subtitle = element_text(size = 11, color = "grey30")
  )

ggsave(file.path(OUT_DIR, "read_retention.png"), p, width = 14, height = 6, dpi = 150)
ggsave(file.path(OUT_DIR, "read_retention.pdf"), p, width = 14, height = 6, device = cairo_pdf)

## Kept so that 04_figures/build_figure2.R can compose Figure 2 panel a
## without rebuilding the plot.
saveRDS(p, file.path(OUT_DIR, "figure2_panel_a.rds"))
