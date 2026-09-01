# ---------------------------------------------------------------------------
# PERCIVAL — central path configuration
#
# All scripts in this repository source this file. No absolute path appears
# anywhere else. Resolution order for the repository root:
#   1. environment variable PERCIVAL_ROOT
#   2. the parent of the directory holding this file (config/paths.R), which
#      is the repository root by construction
#   3. the current working directory
#
# Note that `here::here()` is deliberately NOT used: it anchors on the nearest
# .Rproj / .git above the working directory, which is the wrong root whenever
# this repository sits inside a larger analysis project.
#
# Large input matrices (RNAseq_rawcounts.tsv) are NOT versioned in this
# repository. Download them from Zenodo (see data/README.md) into data/ , or
# point PERCIVAL_DATA at a directory that already holds them.
# ---------------------------------------------------------------------------

percival_root <- function() {
  env <- Sys.getenv("PERCIVAL_ROOT", unset = NA_character_)
  if (!is.na(env) && nzchar(env)) return(normalizePath(env, mustWork = TRUE))
  # `ofile` is set by source() in the calling frame and points at this file.
  for (i in rev(seq_len(sys.nframe()))) {
    of <- try(get("ofile", envir = sys.frame(i), inherits = FALSE), silent = TRUE)
    if (!inherits(of, "try-error") && is.character(of) && length(of) == 1L)
      return(dirname(dirname(normalizePath(of, mustWork = TRUE))))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

ROOT     <- percival_root()
DATA_DIR <- Sys.getenv("PERCIVAL_DATA", unset = file.path(ROOT, "data"))
META_DIR <- file.path(ROOT, "metadata")
RES_DIR  <- Sys.getenv("PERCIVAL_RESULTS", unset = file.path(ROOT, "results"))

COUNTS_FILE  <- file.path(DATA_DIR, "RNAseq_rawcounts.tsv")
TARGETS_FILE <- file.path(META_DIR, "targets.txt")

dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

# Fail early and legibly rather than deep inside an analysis.
require_inputs <- function(...) {
  files <- c(...)
  missing <- files[!file.exists(files)]
  if (length(missing)) {
    stop("Missing required input file(s):\n  ",
         paste(missing, collapse = "\n  "),
         "\nSee data/README.md for how to obtain them.", call. = FALSE)
  }
  invisible(TRUE)
}
