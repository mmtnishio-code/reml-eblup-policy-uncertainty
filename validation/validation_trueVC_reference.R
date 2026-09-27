# ============================================================================
# Validate the True-VC reference used in the manuscript.
# The realised c(Y) is retained; only prediction/PEC variance components are
# replaced by their simulation truth.  Base R only.
# ============================================================================

GSE_expected_trueVC <- function() {
  data.frame(
    h2 = rep(c(0.05, 0.20, 0.40), each = 3),
    n_pheno = rep(c(60, 90, 120), times = 3),
    coverage_expected = c(
      0.949, 0.938, 0.945,
      0.948, 0.942, 0.946,
      0.950, 0.947, 0.944
    ),
    MSE_ratio_expected = c(
      0.971, 0.952, 1.009,
      0.973, 0.969, 0.960,
      0.954, 1.005, 0.964
    )
  )
}

validate_GSE_trueVC_production <- function(
    x,
    rounded_tolerance = 0.0015,
    stop_on_failure = TRUE) {

  tab <- if (is.character(x) && length(x) == 1L) {
    read.csv(x, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.data.frame(x)) {
    x
  } else if (is.list(x) && !is.null(x$method_table)) {
    x$method_table
  } else {
    stop("x must be a True-VC method-table CSV, data.frame, or object with $method_table.")
  }

  req <- c("h2", "n_pheno", "method", "inclusion", "MSE_ratio")
  miss <- setdiff(req, names(tab))
  if (length(miss) > 0L) stop("Missing True-VC columns: ", paste(miss, collapse = ", "))

  accepted_names <- c("true_VC_oracle", "conditional_trueVC_oracle")
  z <- tab[tab$method %in% accepted_names, req, drop = FALSE]
  if (nrow(z) != 9L) {
    stop("Expected exactly 9 True-VC scenario rows; found ", nrow(z), ".")
  }

  exp <- GSE_expected_trueVC()
  obs <- merge(exp, z, by = c("h2", "n_pheno"), all.x = TRUE, sort = FALSE)
  obs$coverage_abs_diff <- abs(obs$inclusion - obs$coverage_expected)
  obs$MSE_ratio_abs_diff <- abs(obs$MSE_ratio - obs$MSE_ratio_expected)
  obs$coverage_ok <- is.finite(obs$coverage_abs_diff) & obs$coverage_abs_diff <= rounded_tolerance
  obs$MSE_ratio_ok <- is.finite(obs$MSE_ratio_abs_diff) & obs$MSE_ratio_abs_diff <= rounded_tolerance
  obs$pass <- obs$coverage_ok & obs$MSE_ratio_ok
  all_ok <- all(obs$pass)

  if (isTRUE(stop_on_failure) && !all_ok) {
    print(obs[!obs$pass, , drop = FALSE], row.names = FALSE)
    stop("True-VC production validation failed. Investigate provenance/settings rather than forcing values.")
  }

  cat("True-VC archived-value validation: ", if (all_ok) "PASS" else "CHECK REQUIRED", "\n", sep = "")
  invisible(obs)
}

run_GSE_trueVC_smoke_test <- function() {
  this <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  vdir <- if (!is.null(this)) dirname(normalizePath(this, winslash = "/", mustWork = FALSE)) else file.path(getwd(), "validation")
  root <- normalizePath(file.path(vdir, ".."), winslash = "/", mustWork = FALSE)
  source(file.path(root, "01_GSE_main_and_diagnostics.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")

  td <- file.path(tempdir(), "gse_trueVC_smoke")
  unlink(td, recursive = TRUE, force = TRUE)
  z <- run_true_VC_oracle_final(
    S = 2,
    h2_values = 0.20,
    info_fractions = c(moderate = 0.70),
    checkpoint_every = 1,
    output_dir = td,
    resume = FALSE,
    keep_raw_in_master = TRUE
  )
  stopifnot(
    nrow(z$method_table) == 2L,
    all(c("conditional", "true_VC_oracle") %in% z$method_table$method),
    all(is.finite(z$method_table$RMSE)),
    all(is.finite(z$method_table$MSE_ratio)),
    all(z$method_table$inclusion >= 0 & z$method_table$inclusion <= 1)
  )
  cat("True-VC smoke test: PASS\n")
  invisible(z)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L || identical(args[[1]], "selftest")) {
    run_GSE_trueVC_smoke_test()
  } else {
    z <- validate_GSE_trueVC_production(args[[1]], stop_on_failure = TRUE)
    print(z, row.names = FALSE)
  }
}
