# ============================================================================
# Validation for 03_GSE_RAMN_matched_replicates.R
# Base R only.
# ============================================================================

.this <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.vdir <- if (!is.null(.this)) dirname(normalizePath(.this, winslash = "/", mustWork = FALSE)) else file.path(getwd(), "validation")
.root <- normalizePath(file.path(.vdir, ".."), winslash = "/", mustWork = FALSE)
source(file.path(.root, "03_GSE_RAMN_matched_replicates.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")

GSE_expected_RAMN_matched <- function() {
  scenarios <- data.frame(
    h2 = rep(c(0.05, 0.20, 0.40), each = 3),
    n_pheno = rep(c(60, 90, 120), times = 3),
    applicability = c(0.402, 0.456, 0.506, 0.635, 0.764, 0.851, 0.835, 0.931, 0.972)
  )
  methods <- c("Conditional", "RL-UP", "RAM-N")
  coverage <- rbind(
    c(0.614, 0.826, 0.836),
    c(0.739, 0.853, 0.868),
    c(0.834, 0.901, 0.903),
    c(0.724, 0.951, 0.953),
    c(0.808, 0.944, 0.942),
    c(0.827, 0.957, 0.949),
    c(0.722, 0.969, 0.981),
    c(0.785, 0.947, 0.959),
    c(0.835, 0.941, 0.939)
  )
  ratio <- rbind(
    c(0.189, 0.479, 0.657),
    c(0.304, 0.593, 0.854),
    c(0.393, 0.692, 0.993),
    c(0.292, 0.929, 1.081),
    c(0.413, 1.012, 1.162),
    c(0.481, 1.033, 1.168),
    c(0.314, 1.253, 1.294),
    c(0.413, 1.055, 1.066),
    c(0.498, 0.972, 0.962)
  )
  out <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
    data.frame(
      h2 = scenarios$h2[i],
      n_pheno = scenarios$n_pheno[i],
      RAMN_applicability_expected = scenarios$applicability[i],
      method = methods,
      coverage_expected = coverage[i, ],
      MSE_ratio_expected = ratio[i, ],
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

validate_GSE_RAMN_matched_production <- function(
    method_table,
    applicability_table = NULL,
    rounded_tolerance = 0.0015,
    stop_on_failure = TRUE) {

  mt <- if (is.character(method_table) && length(method_table) == 1L) {
    read.csv(method_table, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.data.frame(method_table)) method_table else stop("method_table must be a CSV path or data.frame.")

  req <- c("h2", "n_pheno", "method", "coverage", "MSE_ratio")
  miss <- setdiff(req, names(mt))
  if (length(miss) > 0L) stop("Missing RAM-N matched columns: ", paste(miss, collapse = ", "))

  exp <- GSE_expected_RAMN_matched()
  obs <- merge(exp, mt[, req, drop = FALSE], by = c("h2", "n_pheno", "method"), all.x = TRUE, sort = FALSE)
  obs$coverage_abs_diff <- abs(obs$coverage - obs$coverage_expected)
  obs$MSE_ratio_abs_diff <- abs(obs$MSE_ratio - obs$MSE_ratio_expected)
  obs$coverage_ok <- is.finite(obs$coverage_abs_diff) & obs$coverage_abs_diff <= rounded_tolerance
  obs$MSE_ratio_ok <- is.finite(obs$MSE_ratio_abs_diff) & obs$MSE_ratio_abs_diff <= rounded_tolerance
  obs$pass <- obs$coverage_ok & obs$MSE_ratio_ok

  app_check <- NULL
  if (!is.null(applicability_table)) {
    at <- if (is.character(applicability_table) && length(applicability_table) == 1L) {
      read.csv(applicability_table, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (is.data.frame(applicability_table)) applicability_table else stop("applicability_table must be CSV path/data.frame.")
    if (!all(c("h2", "n_pheno", "RAMN_applicability") %in% names(at))) stop("RAM-N applicability table lacks required columns.")
    eapp <- unique(exp[, c("h2", "n_pheno", "RAMN_applicability_expected")])
    app_check <- merge(eapp, at[, c("h2", "n_pheno", "RAMN_applicability")], by = c("h2", "n_pheno"), all.x = TRUE, sort = FALSE)
    app_check$abs_diff <- abs(app_check$RAMN_applicability - app_check$RAMN_applicability_expected)
    app_check$pass <- is.finite(app_check$abs_diff) & app_check$abs_diff <= rounded_tolerance
  }

  all_ok <- all(obs$pass) && (is.null(app_check) || all(app_check$pass))
  if (isTRUE(stop_on_failure) && !all_ok) {
    print(obs[!obs$pass, , drop = FALSE], row.names = FALSE)
    if (!is.null(app_check)) print(app_check[!app_check$pass, , drop = FALSE], row.names = FALSE)
    stop("RAM-N matched production validation failed. Investigate source RDS/applicability selection; do not force values.")
  }

  cat("RAM-N matched archived-value validation: ", if (all_ok) "PASS" else "CHECK REQUIRED", "\n", sep = "")
  invisible(list(method_check = obs, applicability_check = app_check))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L || identical(args[[1]], "selftest")) {
    self_test_RAMN_matched_replicate_comparison()
  } else if (identical(args[[1]], "validate")) {
    if (length(args) < 2L) stop("validate mode requires output directory or method-table CSV.")
    p <- args[[2]]
    if (dir.exists(p)) {
      mt <- file.path(p, "RAMN_matched_method_performance.csv")
      at <- file.path(p, "RAMN_matched_applicability.csv")
    } else {
      mt <- p
      at <- if (length(args) >= 3L) args[[3]] else NULL
    }
    ans <- validate_GSE_RAMN_matched_production(mt, at, stop_on_failure = TRUE)
    print(ans$method_check, row.names = FALSE)
    if (!is.null(ans$applicability_check)) print(ans$applicability_check, row.names = FALSE)
  } else {
    stop("Use selftest or validate.")
  }
}
