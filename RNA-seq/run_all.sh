#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Runs the whole downstream (post-alignment) analysis in dependency order.
# The HPC steps in 01_preprocessing_alignment/ are NOT run here: they need a
# SLURM cluster and the raw FASTQ files. See that directory's README.
# ---------------------------------------------------------------------------
set -euo pipefail

PERCIVAL_ROOT="${PERCIVAL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
export PERCIVAL_ROOT
echo "Repository root: $PERCIVAL_ROOT"

run() {
  echo; echo "=============================================================="
  echo ">>> $1"
  echo "=============================================================="
  Rscript "$PERCIVAL_ROOT/$1"
}

run 02_rnaseq_qc/plot_read_retention.R
run 02_rnaseq_qc/plot_sample_distance_heatmap.R
run 02_rnaseq_qc/rnaseq_qc_table_and_plots.R
run 03_batch_effect/variance_partition.R
run 03_batch_effect/build_supplementary_table.R
run 04_figures/build_figure2.R
run 04_figures/build_figure2_pagefit.R

echo; echo "Done. Outputs in ${PERCIVAL_RESULTS:-$PERCIVAL_ROOT/results}"
