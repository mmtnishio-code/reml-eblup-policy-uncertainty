# ============================================================================
# Validate completed production outputs against the archived manuscript values.
# This script NEVER edits analysis output.  A mismatch is reported as an error.
# Base R only.
#
# Usage:
#   Rscript validation/validate_all_production_outputs.R \
#     --main path/to/Simulation_II_method_performance.csv \
#     --truevc path/to/TrueVC_method_performance.csv \
#     --rlup path/to/RLUP_prior_sensitivity_MASTER.rds \
#     --matched path/to/RAMN_matched_output_dir \
#     --pedigree path/to/pedigree_structure_sensitivity_FINAL_comparison.csv \
#     --reselection path/to/reselection_bootstrap_master_only_ALL.rds
# ============================================================================

.this <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.vdir <- if (!is.null(.this)) dirname(normalizePath(.this, winslash = "/", mustWork = FALSE)) else file.path(getwd(), "validation")
.root <- normalizePath(file.path(.vdir, ".."), winslash = "/", mustWork = FALSE)

source(file.path(.vdir, "validate_01_main_production.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")
source(file.path(.vdir, "validation_trueVC_reference.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")
source(file.path(.root, "01_GSE_main_and_diagnostics.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")
source(file.path(.root, "02_GSE_RLUP_prior_sensitivity_RNG_FIXED.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")
source(file.path(.vdir, "validate_03_RAMN_matched.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")
source(file.path(.vdir, "validate_04_pedigree_structure.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")
options(gse.validate06.no_dispatch = TRUE)
source(file.path(.vdir, "validate_06_reselection_bootstrap.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")
options(gse.validate06.no_dispatch = NULL)

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--") || i == length(x)) stop("Arguments must be --name value pairs.")
    out[[substring(key, 3L)]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (length(args) == 0L) stop("No production outputs supplied. See header usage.")

if (!is.null(args$main)) validate_GSE_main_simII_production(args$main, stop_on_failure = TRUE)
if (!is.null(args$truevc)) validate_GSE_trueVC_production(args$truevc, stop_on_failure = TRUE)
if (!is.null(args$rlup)) {
  x <- if (grepl("\\.rds$", args$rlup, ignore.case = TRUE)) readRDS(args$rlup) else args$rlup
  validate_RLUP_prior_sensitivity_production(x, stop_on_failure = TRUE)
}
if (!is.null(args$matched)) {
  d <- args$matched
  mt <- if (dir.exists(d)) file.path(d, "RAMN_matched_method_performance.csv") else d
  at <- if (dir.exists(d)) file.path(d, "RAMN_matched_applicability.csv") else NULL
  validate_GSE_RAMN_matched_production(mt, at, stop_on_failure = TRUE)
}
if (!is.null(args$pedigree)) {
  GSE_PEDIGREE_RUN_MODE <- "design_only"
  source(file.path(.root, "04_GSE_pedigree_structure_sensitivity.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")
  validate_pedigree_structure_production(args$pedigree, stop_on_failure = TRUE)
}
if (!is.null(args$reselection)) {
  validate_GSE_reselection_bootstrap_production(
    args$reselection,
    stop_on_failure = TRUE
  )
}

if (!is.null(args$pedigree) && !is.null(args$main) && !is.null(args$truevc)) {
  validate_GSE_pedigree_baseline_against_main(
    pedigree_csv = args$pedigree,
    main_method_csv = args$main,
    trueVC_method_csv = args$truevc,
    tolerance = 1e-8,
    stop_on_failure = TRUE
  )
}

cat("\nAll supplied production-output checks passed.\n")
