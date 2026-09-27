# ============================================================================
# Validation for 06_GSE_procedure_level_reselection_bootstrap.R
# Base R only.
#
# Usage:
#   Rscript validation/validate_06_reselection_bootstrap.R smoke
#   Rscript validation/validate_06_reselection_bootstrap.R validate \
#     paper_output/reselection_bootstrap/reselection_bootstrap_master_only_ALL.rds
#   Rscript validation/validate_06_reselection_bootstrap.R outer \
#     paper_output/bootstrap_diagnostic/Simulation_II_bootstrap_generator_ALL.rds
#
# "smoke" runs only a small new-function test.
# "validate" checks a completed S=1000/B=500 final add-on object against the
# final manuscript values and its embedded exact-reproduction validation table.
# "outer" reconstructs all original outer samples from boot1000 and verifies
# stored Conditional/True-VC summaries; this is much more expensive than smoke.
# ============================================================================

.this <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.vdir <- if (!is.null(.this)) {
  dirname(normalizePath(.this, winslash = "/", mustWork = FALSE))
} else {
  file.path(getwd(), "validation")
}
.root <- normalizePath(file.path(.vdir, ".."), winslash = "/", mustWork = FALSE)

source(file.path(.root, "01_GSE_main_and_diagnostics.R"),
       echo = FALSE, chdir = FALSE, encoding = "UTF-8")
source(file.path(.root, "06_GSE_procedure_level_reselection_bootstrap.R"),
       echo = FALSE, chdir = FALSE, encoding = "UTF-8")

validate_GSE_reselection_bootstrap_production <- function(
    x,
    rounded_tolerance = 0.002,
    stop_on_failure = TRUE) {

  if (is.character(x) && length(x) == 1L) {
    if (dir.exists(x)) {
      x <- file.path(x, "reselection_bootstrap_master_only_ALL.rds")
    }
    if (!file.exists(x)) stop("Reselection-bootstrap result not found: ", x)
    x <- readRDS(x)
  }

  req <- c(
    "settings", "comparison_table", "validation_table",
    "per_outer_diagnostics_table"
  )
  miss <- req[!req %in% names(x)]
  if (length(miss) > 0L) {
    stop("Final reselection object is missing: ", paste(miss, collapse = ", "))
  }

  ok <- TRUE
  messages <- character()

  if (!isTRUE(x$settings$full_production) ||
      as.integer(x$settings$S_reference) != 1000L ||
      as.integer(x$settings$B_inner) != 500L) {
    ok <- FALSE
    messages <- c(messages, "Result is not marked as the full S=1000/B=500 production run.")
  }

  vv <- x$validation_table
  if (nrow(vv) == 0L || any(vv$within_tolerance %in% FALSE) || any(is.na(vv$within_tolerance))) {
    ok <- FALSE
    messages <- c(messages, "Embedded exact-reproduction validation contains FALSE/NA or is empty.")
  }

  cmp <- x$comparison_table
  cmp <- cmp[order(cmp$h2, cmp$n_pheno), , drop = FALSE]

  expected_cmp <- data.frame(
    h2 = rep(c(0.05, 0.20, 0.40), each = 3),
    n_pheno = rep(c(60, 90, 120), times = 3),
    stored_fixed_full_coverage = c(
      0.846, 0.799, 0.776,
      0.653, 0.722, 0.758,
      0.674, 0.764, 0.828
    ),
    paired_fixed_new_coverage = c(
      0.847, 0.806, 0.773,
      0.661, 0.718, 0.763,
      0.670, 0.769, 0.837
    ),
    reselect_coverage = c(
      0.613, 0.645, 0.650,
      0.656, 0.758, 0.810,
      0.753, 0.851, 0.897
    ),
    stored_fixed_full_MSE_ratio = c(
      0.292, 0.373, 0.416,
      0.294, 0.370, 0.442,
      0.285, 0.416, 0.513
    ),
    reselect_MSE_ratio = c(
      0.818, 0.809, 0.773,
      0.693, 0.707, 0.735,
      0.649, 0.791, 0.854
    ),
    row.names = NULL,
    check.names = FALSE
  )

  if (nrow(cmp) != nrow(expected_cmp) ||
      any(abs(cmp$h2 - expected_cmp$h2) > 1e-12) ||
      any(cmp$n_pheno != expected_cmp$n_pheno)) {
    ok <- FALSE
    messages <- c(messages, "Comparison-table scenario rows do not match the 3 x 3 manuscript design.")
  } else {
    cols <- setdiff(names(expected_cmp), c("h2", "n_pheno"))
    for (nm in cols) {
      dd <- max(abs(as.numeric(cmp[[nm]]) - expected_cmp[[nm]]), na.rm = TRUE)
      if (!is.finite(dd) || dd > rounded_tolerance) {
        ok <- FALSE
        messages <- c(
          messages,
          paste0(nm, " differs from the archived manuscript values; max abs diff = ", format(dd, digits = 6))
        )
      }
    }
  }

  s10 <- rbootm_make_supplementary_S10(x)
  s10 <- s10[order(s10$h2, s10$n_pheno), , drop = FALSE]
  expected_s10 <- data.frame(
    h2 = rep(c(0.05, 0.20, 0.40), each = 3),
    n_pheno = rep(c(60, 90, 120), times = 3),
    truth_below_lower = c(
      0.028, 0.015, 0.022,
      0.038, 0.025, 0.026,
      0.072, 0.045, 0.037
    ),
    truth_above_upper = c(
      0.359, 0.340, 0.328,
      0.306, 0.217, 0.164,
      0.175, 0.104, 0.066
    ),
    corr_outer_error_inner_mean_error = c(
      -0.8215388, -0.7693940, -0.7240360,
      -0.7293377, -0.5968656, -0.5147117,
      -0.7065925, -0.4130196, -0.1690633
    ),
    row.names = NULL,
    check.names = FALSE
  )

  if (nrow(s10) != 9L) {
    ok <- FALSE
    messages <- c(messages, "Supplementary S10 reconstruction does not have 9 scenario rows.")
  } else {
    for (nm in c("truth_below_lower", "truth_above_upper", "corr_outer_error_inner_mean_error")) {
      dd <- max(abs(as.numeric(s10[[nm]]) - expected_s10[[nm]]), na.rm = TRUE)
      if (!is.finite(dd) || dd > rounded_tolerance) {
        ok <- FALSE
        messages <- c(messages, paste0("S10 ", nm, " mismatch; max abs diff = ", format(dd, digits = 6)))
      }
    }
  }

  s11 <- rbootm_make_supplementary_S11(x)
  expected_s11 <- data.frame(
    h2 = rep(c(0.05, 0.20, 0.40), each = 5),
    inner_boundary_rate_stratum = rep(
      c("<=0.05", "0.05-0.10", "0.10-0.20", "0.20-0.40", ">0.40"),
      times = 3
    ),
    n = c(
      252, 264, 516, 1719, 249,
      847, 524, 630, 862, 137,
      1906, 401, 329, 304, 60
    ),
    coverage = c(
      0.753968, 0.988636, 1.000000, 0.509599, 0.261044,
      0.896104, 0.994275, 0.974603, 0.371230, 0.072993,
      0.913431, 0.957606, 0.860182, 0.299342, 0.033333
    ),
    mean_inner_error_minus_outer_error = c(
      -0.451767, -0.276802, -0.137636, 0.113803, 0.152428,
      -0.257779, -0.077475, 0.044013, 0.267393, 0.311202,
      -0.121053, 0.081974, 0.195828, 0.408705, 0.498579
    ),
    row.names = NULL,
    check.names = FALSE
  )

  s11 <- s11[order(s11$h2, match(
    s11$inner_boundary_rate_stratum,
    c("<=0.05", "0.05-0.10", "0.10-0.20", "0.20-0.40", ">0.40")
  )), , drop = FALSE]

  if (nrow(s11) != 15L || any(s11$n != expected_s11$n)) {
    ok <- FALSE
    messages <- c(messages, "Supplementary S11 stratum counts do not match the archived manuscript values.")
  } else {
    for (nm in c("coverage", "mean_inner_error_minus_outer_error")) {
      dd <- max(abs(as.numeric(s11[[nm]]) - expected_s11[[nm]]), na.rm = TRUE)
      if (!is.finite(dd) || dd > rounded_tolerance) {
        ok <- FALSE
        messages <- c(messages, paste0("S11 ", nm, " mismatch; max abs diff = ", format(dd, digits = 6)))
      }
    }
  }

  if (ok) {
    cat("Procedure-level reselection-bootstrap production validation: PASS\n")
  } else {
    cat("Procedure-level reselection-bootstrap production validation: FAIL\n")
    cat(paste0("  - ", messages, collapse = "\n"), "\n")
    if (isTRUE(stop_on_failure)) stop("Reselection-bootstrap production validation failed.")
  }

  invisible(list(ok = ok, messages = messages, S10 = s10, S11 = s11))
}

.validate06_no_dispatch <- isTRUE(getOption("gse.validate06.no_dispatch", FALSE))

if (!.validate06_no_dispatch) {
  args <- commandArgs(trailingOnly = TRUE)
  mode <- if (length(args) >= 1L) tolower(args[[1]]) else "smoke"

  if (mode == "smoke") {
    self_test_reselection_bootstrap_master_only(B = 20)
    cat("validate_06 smoke test: PASS\n")
  } else if (mode == "validate") {
    if (length(args) < 2L) {
      stop("Usage: ... validate path/to/reselection_bootstrap_master_only_ALL.rds")
    }
    validate_GSE_reselection_bootstrap_production(
      args[[2]],
      stop_on_failure = TRUE
    )
  } else if (mode == "outer") {
    if (length(args) < 2L) {
      stop("Usage: ... outer path/to/Simulation_II_bootstrap_generator_ALL.rds")
    }
    boot1000 <- readRDS(args[[2]])
    z <- validate_boot1000_outer_reconstruction(
      boot1000 = boot1000,
      tolerance = 1e-8,
      stop_on_failure = TRUE
    )
    stopifnot(isTRUE(z$all_ok))
    cat("Exact outer reconstruction validation: PASS\n")
  } else {
    stop("Unknown mode: ", mode, ". Use smoke, validate, or outer.")
  }
}
