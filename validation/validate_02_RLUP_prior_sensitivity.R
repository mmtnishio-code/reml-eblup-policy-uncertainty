# ============================================================================
# Validation driver for 02_GSE_RLUP_prior_sensitivity_RNG_FIXED.R
# ============================================================================
# Run from the GSE_public_code directory:
#
#   Rscript validation/validate_02_RLUP_prior_sensitivity.R selftest
#
# To validate an already completed S=1000 output without rerunning it:
#
#   Rscript validation/validate_02_RLUP_prior_sensitivity.R validate \
#     RLUP_prior_sensitivity_S1000_RNG_FIXED/RLUP_prior_sensitivity_MASTER.rds
#
# Full production run (computationally intensive):
#
#   Rscript validation/validate_02_RLUP_prior_sensitivity.R production \
#     RLUP_prior_sensitivity_S1000_RNG_FIXED
#
# This file does not modify the scientific calculation.  It only calls the
# functions defined in 02 and checks the archived manuscript summaries.
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1L) args[[1]] else "selftest"

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE)),
  error = function(e) file.path(getwd(), "validation")
)
root_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)

source(file.path(root_dir, "01_GSE_main_and_diagnostics.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")
source(file.path(root_dir, "02_GSE_RLUP_prior_sensitivity_RNG_FIXED.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")

check_RLUP_manuscript_configuration()

if (identical(mode, "selftest")) {
  run_RLUP_release_self_test()
  quit(save = "no", status = 0L)
}

if (identical(mode, "validate")) {
  if (length(args) < 2L) {
    stop("validate mode requires the path to RLUP_prior_sensitivity_MASTER.rds")
  }
  x <- readRDS(args[[2]])
  v <- validate_RLUP_prior_sensitivity_production(x, stop_on_failure = TRUE)
  print(v, row.names = FALSE)
  quit(save = "no", status = 0L)
}

if (identical(mode, "production")) {
  outdir <- if (length(args) >= 2L) {
    args[[2]]
  } else {
    file.path(root_dir, "RLUP_prior_sensitivity_S1000_RNG_FIXED")
  }

  run_RLUP_prior_sensitivity_production(
    output_dir = outdir,
    keep_raw = TRUE,
    resume = TRUE
  )
  quit(save = "no", status = 0L)
}

stop("Unknown mode: ", mode, ". Use selftest, validate, or production.")
