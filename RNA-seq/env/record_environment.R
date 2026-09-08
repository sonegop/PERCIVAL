#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Records the R environment this repository was run in.
#
# Writes env/r_packages.tsv (the packages the analysis scripts load, with the
# installed version) and env/sessionInfo.txt (the full session, after loading
# them, so the transitive dependency versions are captured too).
#
# Usage: Rscript env/record_environment.R
# ---------------------------------------------------------------------------

PKGS <- c(
  # variance decomposition and differential expression
  "variancePartition", "limma", "edgeR", "DESeq2", "matrixStats",
  # figures
  "ggplot2", "cowplot", "patchwork", "ComplexHeatmap", "circlize", "RColorBrewer", "ggrepel",
  # data handling and I/O
  "tidyverse", "dplyr", "tidyr", "readr", "openxlsx", "jsonlite"
)

.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.d <- if (length(.f)) dirname(normalizePath(sub("^--file=", "", .f[1]))) else normalizePath(getwd())
OUT <- .d

versions <- vapply(PKGS, function(p) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  if (is.na(v)) "not installed" else v
}, character(1))

write.table(
  data.frame(package = PKGS, installed = versions != "not installed", version = versions),
  file.path(OUT, "r_packages.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

invisible(lapply(PKGS[versions != "not installed"],
                 function(p) suppressPackageStartupMessages(
                   library(p, character.only = TRUE))))

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
message("Wrote r_packages.tsv and sessionInfo.txt to ", OUT)
