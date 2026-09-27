# ============================================================================
# Validate archived Simulation-II main production summaries against manuscript
# Tables 2-3 and RAM-N applicability counts.  Base R only.
# ============================================================================

GSE_expected_main_simII <- function() {
  scenarios <- data.frame(
    h2 = rep(c(0.05, 0.20, 0.40), each = 3),
    n_pheno = rep(c(60, 90, 120), times = 3)
  )

  methods <- c("conditional", "RAM_N_policy", "Bayesian_policy", "bootstrap_policy")

  coverage <- rbind(
    c(0.363, 0.836, 0.905, 0.846),
    c(0.454, 0.868, 0.910, 0.799),
    c(0.520, 0.903, 0.931, 0.776),
    c(0.501, 0.953, 0.960, 0.653),
    c(0.646, 0.942, 0.948, 0.722),
    c(0.715, 0.949, 0.952, 0.758),
    c(0.611, 0.981, 0.965, 0.674),
    c(0.736, 0.959, 0.946, 0.764),
    c(0.813, 0.939, 0.942, 0.828)
  )

  ratio <- rbind(
    c(0.180, 0.657, 0.692, 0.292),
    c(0.270, 0.854, 0.810, 0.373),
    c(0.327, 0.993, 0.898, 0.416),
    c(0.209, 1.081, 1.054, 0.294),
    c(0.305, 1.162, 1.038, 0.370),
    c(0.382, 1.168, 1.009, 0.442),
    c(0.215, 1.294, 1.142, 0.285),
    c(0.355, 1.066, 1.028, 0.416),
    c(0.464, 0.962, 0.954, 0.513)
  )

  out <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
    data.frame(
      h2 = scenarios$h2[i],
      n_pheno = scenarios$n_pheno[i],
      method = methods,
      coverage_expected = coverage[i, ],
      MSE_ratio_expected = ratio[i, ],
      stringsAsFactors = FALSE
    )
  }))

  rownames(out) <- NULL
  out
}

GSE_expected_RAMN_applicability <- function() {
  data.frame(
    h2 = rep(c(0.05, 0.20, 0.40), each = 3),
    n_pheno = rep(c(60, 90, 120), times = 3),
    valid_expected = c(402, 456, 506, 635, 764, 851, 835, 931, 972),
    valid_rate_expected = c(0.402, 0.456, 0.506, 0.635, 0.764, 0.851, 0.835, 0.931, 0.972)
  )
}

validate_GSE_main_simII_production <- function(
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
    stop("x must be Simulation_II_method_performance.csv, a data.frame, or an object with $method_table.")
  }

  req <- c("h2", "n_pheno", "method", "inclusion", "MSE_ratio")
  miss <- setdiff(req, names(tab))
  if (length(miss) > 0L) stop("Missing main-output columns: ", paste(miss, collapse = ", "))

  exp <- GSE_expected_main_simII()
  obs <- merge(
    exp,
    tab[, intersect(c(req, "valid", "valid_rate"), names(tab)), drop = FALSE],
    by = c("h2", "n_pheno", "method"),
    all.x = TRUE,
    sort = FALSE
  )

  obs$coverage_abs_diff <- abs(obs$inclusion - obs$coverage_expected)
  obs$MSE_ratio_abs_diff <- abs(obs$MSE_ratio - obs$MSE_ratio_expected)
  obs$coverage_ok <- is.finite(obs$coverage_abs_diff) & obs$coverage_abs_diff <= rounded_tolerance
  obs$MSE_ratio_ok <- is.finite(obs$MSE_ratio_abs_diff) & obs$MSE_ratio_abs_diff <= rounded_tolerance
  obs$pass <- obs$coverage_ok & obs$MSE_ratio_ok

  ramn_check <- NULL
  if (all(c("valid", "valid_rate") %in% names(tab))) {
    r <- tab[tab$method == "RAM_N_policy", c("h2", "n_pheno", "valid", "valid_rate"), drop = FALSE]
    er <- GSE_expected_RAMN_applicability()
    ramn_check <- merge(er, r, by = c("h2", "n_pheno"), all.x = TRUE, sort = FALSE)
    ramn_check$valid_ok <- is.finite(ramn_check$valid) & ramn_check$valid == ramn_check$valid_expected
    ramn_check$valid_rate_ok <- is.finite(ramn_check$valid_rate) &
      abs(ramn_check$valid_rate - ramn_check$valid_rate_expected) <= rounded_tolerance
  }

  all_ok <- all(obs$pass)
  if (!is.null(ramn_check)) all_ok <- all_ok && all(ramn_check$valid_ok & ramn_check$valid_rate_ok)

  if (isTRUE(stop_on_failure) && !all_ok) {
    print(obs[!obs$pass, , drop = FALSE], row.names = FALSE)
    if (!is.null(ramn_check)) print(ramn_check[!(ramn_check$valid_ok & ramn_check$valid_rate_ok), , drop = FALSE], row.names = FALSE)
    stop("Main Simulation-II production validation failed. Investigate settings/RNG/output provenance; do not force values.")
  }

  cat("Main Simulation-II archived-value validation: ", if (all_ok) "PASS" else "CHECK REQUIRED", "\n", sep = "")
  invisible(list(method_check = obs, RAMN_applicability_check = ramn_check))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1L) {
    stop("Usage: Rscript validation/validate_01_main_production.R path/to/Simulation_II_method_performance.csv")
  }
  ans <- validate_GSE_main_simII_production(args[[1]], stop_on_failure = TRUE)
  print(ans$method_check, row.names = FALSE)
  if (!is.null(ans$RAMN_applicability_check)) print(ans$RAMN_applicability_check, row.names = FALSE)
}
