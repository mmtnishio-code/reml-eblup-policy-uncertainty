# ============================================================================
# GSE: RAM-N MATCHED-REPLICATE COMPARISON
#
# Post-processes existing Simulation-II production RDS files. No new outer
# simulation is generated. Conditional, RL-UP, and RAM-N are compared on the
# same outer replicates for which RAM-N passed the pre-specified diagnostics.
# Base R only.
# ============================================================================

# ------------------------------------------------------------
# 1. Small helpers
# ------------------------------------------------------------

ramn_matched_safe_quantile <- function(x, probs = c(0.025, 0.975)) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(rep(NA_real_, length(probs)))
  }
  unname(
    quantile(
      x,
      probs = probs,
      names = FALSE,
      type = 8
    )
  )
}


ramn_matched_wilson_ci <- function(x, n, level = 0.95) {

  if (
    !is.finite(x) ||
    !is.finite(n) ||
    n <= 0
  ) {
    return(c(lo = NA_real_, hi = NA_real_))
  }

  z <- qnorm(1 - (1 - level) / 2)

  phat <- x / n

  den <- 1 + z^2 / n

  center <- (
    phat +
      z^2 / (2 * n)
  ) / den

  half <- (
    z *
      sqrt(
        phat * (1 - phat) / n +
          z^2 / (4 * n^2)
      )
  ) / den

  c(
    lo = max(0, center - half),
    hi = min(1, center + half)
  )
}


ramn_matched_method_label <- function(x) {

  y <- as.character(x)

  y[y == "conditional"] <- "Conditional"
  y[y == "Bayesian_policy"] <- "RL-UP"
  y[y == "RAM_N_policy"] <- "RAM-N"

  y
}


ramn_matched_extract_method_matrix <- function(
    res,
    method_name,
    idx) {

  x <- lapply(
    res[idx],
    `[[`,
    method_name
  )

  bad <- vapply(
    x,
    is.null,
    logical(1)
  )

  if (any(bad)) {
    stop(
      "Method '",
      method_name,
      "' is NULL in ",
      sum(bad),
      " RAM-N-valid replicate(s). ",
      "The comparison would no longer be a pure RAM-N-applicable matched set."
    )
  }

  mat <- do.call(
    rbind,
    x
  )

  required <- c(
    "mean",
    "sd",
    "lo",
    "hi"
  )

  if (!all(required %in% colnames(mat))) {
    stop(
      "Method '",
      method_name,
      "' does not contain required columns: ",
      paste(required, collapse = ", ")
    )
  }

  mat
}


ramn_matched_method_statistics <- function(
    truth,
    mat,
    target = 0) {

  stopifnot(
    length(truth) == nrow(mat)
  )

  err <- mat[, "mean"] - truth

  mc_mse <- mean(
    err^2
  )

  mean_est_mse <- mean(
    mat[, "sd"]^2
  )

  covered <- (
    truth >= mat[, "lo"] &
      truth <= mat[, "hi"]
  )

  n <- length(truth)

  coverage <- mean(
    covered
  )

  coverage_mcse <- sqrt(
    coverage *
      (1 - coverage) /
      n
  )

  coverage_wilson <- ramn_matched_wilson_ci(
    x = sum(covered),
    n = n
  )

  out <- c(
    n = n,
    bias = mean(err),
    RMSE = sqrt(mc_mse),
    mean_SE = mean(mat[, "sd"]),
    MSE_ratio = mean_est_mse / mc_mse,
    coverage = coverage,
    coverage_MCSE = coverage_mcse,
    coverage_Wilson_lo = unname(coverage_wilson["lo"]),
    coverage_Wilson_hi = unname(coverage_wilson["hi"]),
    width = mean(
      mat[, "hi"] -
        mat[, "lo"]
    )
  )

  if ("p" %in% colnames(mat)) {

    z <- as.numeric(
      truth > target
    )

    out <- c(
      out,
      Brier = mean(
        (mat[, "p"] - z)^2
      ),
      mean_p = mean(
        mat[, "p"]
      ),
      event_rate = mean(z)
    )
  }

  out
}


# ------------------------------------------------------------
# 2. Bootstrap Monte Carlo intervals on the matched outer reps
# ------------------------------------------------------------

ramn_matched_bootstrap_one_method <- function(
    truth,
    mat,
    R_boot = 5000,
    seed = 20260870) {

  n <- length(truth)

  if (n < 2L) {
    return(
      c(
        MSE_ratio_lo = NA_real_,
        MSE_ratio_hi = NA_real_,
        RMSE_lo = NA_real_,
        RMSE_hi = NA_real_,
        coverage_lo = NA_real_,
        coverage_hi = NA_real_
      )
    )
  }

  set.seed(seed)

  b_mse_ratio <- numeric(R_boot)
  b_rmse <- numeric(R_boot)
  b_coverage <- numeric(R_boot)

  for (b in seq_len(R_boot)) {

    ii <- sample.int(
      n,
      size = n,
      replace = TRUE
    )

    tr <- truth[ii]
    mm <- mat[ii, , drop = FALSE]

    err <- mm[, "mean"] - tr

    mc_mse <- mean(
      err^2
    )

    est_mse <- mean(
      mm[, "sd"]^2
    )

    b_mse_ratio[b] <- (
      est_mse /
        mc_mse
    )

    b_rmse[b] <- sqrt(
      mc_mse
    )

    b_coverage[b] <- mean(
      tr >= mm[, "lo"] &
        tr <= mm[, "hi"]
    )
  }

  q_ratio <- ramn_matched_safe_quantile(
    b_mse_ratio
  )

  q_rmse <- ramn_matched_safe_quantile(
    b_rmse
  )

  q_cov <- ramn_matched_safe_quantile(
    b_coverage
  )

  c(
    MSE_ratio_lo = q_ratio[1],
    MSE_ratio_hi = q_ratio[2],
    RMSE_lo = q_rmse[1],
    RMSE_hi = q_rmse[2],
    coverage_lo = q_cov[1],
    coverage_hi = q_cov[2]
  )
}


ramn_matched_bootstrap_pair_difference <- function(
    truth,
    mat_A,
    mat_B,
    R_boot = 5000,
    seed = 20260871) {

  n <- length(truth)

  stopifnot(
    nrow(mat_A) == n,
    nrow(mat_B) == n
  )

  if (n < 2L) {
    return(
      c(
        coverage_difference = NA_real_,
        coverage_difference_lo = NA_real_,
        coverage_difference_hi = NA_real_,
        MSE_ratio_difference = NA_real_,
        MSE_ratio_difference_lo = NA_real_,
        MSE_ratio_difference_hi = NA_real_,
        RMSE_difference = NA_real_,
        RMSE_difference_lo = NA_real_,
        RMSE_difference_hi = NA_real_
      )
    )
  }

  coverage_A <- (
    truth >= mat_A[, "lo"] &
      truth <= mat_A[, "hi"]
  )

  coverage_B <- (
    truth >= mat_B[, "lo"] &
      truth <= mat_B[, "hi"]
  )

  err_A <- (
    mat_A[, "mean"] -
      truth
  )

  err_B <- (
    mat_B[, "mean"] -
      truth
  )

  mse_ratio_A <- (
    mean(mat_A[, "sd"]^2) /
      mean(err_A^2)
  )

  mse_ratio_B <- (
    mean(mat_B[, "sd"]^2) /
      mean(err_B^2)
  )

  rmse_A <- sqrt(
    mean(err_A^2)
  )

  rmse_B <- sqrt(
    mean(err_B^2)
  )

  observed <- c(
    coverage_difference =
      mean(coverage_A - coverage_B),
    MSE_ratio_difference =
      mse_ratio_A - mse_ratio_B,
    RMSE_difference =
      rmse_A - rmse_B
  )

  set.seed(seed)

  b_cov <- numeric(R_boot)
  b_ratio <- numeric(R_boot)
  b_rmse <- numeric(R_boot)

  for (b in seq_len(R_boot)) {

    ii <- sample.int(
      n,
      size = n,
      replace = TRUE
    )

    tr <- truth[ii]

    A <- mat_A[ii, , drop = FALSE]
    B <- mat_B[ii, , drop = FALSE]

    errA <- A[, "mean"] - tr
    errB <- B[, "mean"] - tr

    covA <- (
      tr >= A[, "lo"] &
        tr <= A[, "hi"]
    )

    covB <- (
      tr >= B[, "lo"] &
        tr <= B[, "hi"]
    )

    ratioA <- (
      mean(A[, "sd"]^2) /
        mean(errA^2)
    )

    ratioB <- (
      mean(B[, "sd"]^2) /
        mean(errB^2)
    )

    rmseA <- sqrt(
      mean(errA^2)
    )

    rmseB <- sqrt(
      mean(errB^2)
    )

    b_cov[b] <- mean(
      covA - covB
    )

    b_ratio[b] <- (
      ratioA - ratioB
    )

    b_rmse[b] <- (
      rmseA - rmseB
    )
  }

  q_cov <- ramn_matched_safe_quantile(
    b_cov
  )

  q_ratio <- ramn_matched_safe_quantile(
    b_ratio
  )

  q_rmse <- ramn_matched_safe_quantile(
    b_rmse
  )

  c(
    coverage_difference =
      unname(observed["coverage_difference"]),
    coverage_difference_lo =
      q_cov[1],
    coverage_difference_hi =
      q_cov[2],

    MSE_ratio_difference =
      unname(observed["MSE_ratio_difference"]),
    MSE_ratio_difference_lo =
      q_ratio[1],
    MSE_ratio_difference_hi =
      q_ratio[2],

    RMSE_difference =
      unname(observed["RMSE_difference"]),
    RMSE_difference_lo =
      q_rmse[1],
    RMSE_difference_hi =
      q_rmse[2]
  )
}


# ------------------------------------------------------------
# 3. One scenario
# ------------------------------------------------------------

ramn_matched_analyse_scenario <- function(
    obj,
    file_name = NA_character_,
    target = 0,
    R_boot = 5000,
    seed = 20260872) {

  if (
    is.null(obj$raw) ||
    is.null(obj$summary)
  ) {
    stop(
      "RDS does not contain both $raw and $summary: ",
      file_name
    )
  }

  res <- obj$raw

  total_n <- length(res)

  if (total_n == 0L) {
    stop(
      "No outer replicates in: ",
      file_name
    )
  }

  ram_valid <- !vapply(
    lapply(
      res,
      `[[`,
      "RAM_N_policy"
    ),
    is.null,
    logical(1)
  )

  idx <- which(
    ram_valid
  )

  matched_n <- length(idx)

  if (matched_n == 0L) {
    stop(
      "No RAM-N-valid replicates in: ",
      file_name
    )
  }

  truth <- vapply(
    res[idx],
    `[[`,
    numeric(1),
    "truth"
  )

  methods <- c(
    "conditional",
    "Bayesian_policy",
    "RAM_N_policy"
  )

  mats <- lapply(
    methods,
    function(m) {
      ramn_matched_extract_method_matrix(
        res = res,
        method_name = m,
        idx = idx
      )
    }
  )

  names(mats) <- methods

  h2 <- obj$summary$h2
  information <- obj$summary$information
  n_pheno <- obj$summary$n_pheno

  method_rows <- vector(
    "list",
    length(methods)
  )

  for (j in seq_along(methods)) {

    m <- methods[j]

    st <- ramn_matched_method_statistics(
      truth = truth,
      mat = mats[[m]],
      target = target
    )

    ci <- ramn_matched_bootstrap_one_method(
      truth = truth,
      mat = mats[[m]],
      R_boot = R_boot,
      seed = seed + j
    )

    method_rows[[j]] <- data.frame(
      h2 = h2,
      information = information,
      n_pheno = n_pheno,
      total_outer_replicates = total_n,
      RAMN_matched_n = matched_n,
      RAMN_applicability =
        matched_n / total_n,
      method = ramn_matched_method_label(m),

      bias = unname(st["bias"]),
      RMSE = unname(st["RMSE"]),
      RMSE_MC_lo = unname(ci["RMSE_lo"]),
      RMSE_MC_hi = unname(ci["RMSE_hi"]),

      mean_SE = unname(st["mean_SE"]),

      MSE_ratio = unname(st["MSE_ratio"]),
      MSE_ratio_MC_lo =
        unname(ci["MSE_ratio_lo"]),
      MSE_ratio_MC_hi =
        unname(ci["MSE_ratio_hi"]),

      coverage = unname(st["coverage"]),
      coverage_MCSE =
        unname(st["coverage_MCSE"]),
      coverage_Wilson_lo =
        unname(st["coverage_Wilson_lo"]),
      coverage_Wilson_hi =
        unname(st["coverage_Wilson_hi"]),
      coverage_outer_boot_lo =
        unname(ci["coverage_lo"]),
      coverage_outer_boot_hi =
        unname(ci["coverage_hi"]),

      width = unname(st["width"]),

      stringsAsFactors = FALSE
    )
  }

  method_table <- do.call(
    rbind,
    method_rows
  )

  pair_defs <- list(
    c("conditional", "RAM_N_policy"),
    c("Bayesian_policy", "RAM_N_policy"),
    c("conditional", "Bayesian_policy")
  )

  pair_rows <- vector(
    "list",
    length(pair_defs)
  )

  for (j in seq_along(pair_defs)) {

    pp <- pair_defs[[j]]

    A <- pp[1]
    B <- pp[2]

    dd <- ramn_matched_bootstrap_pair_difference(
      truth = truth,
      mat_A = mats[[A]],
      mat_B = mats[[B]],
      R_boot = R_boot,
      seed = seed + 100 + j
    )

    pair_rows[[j]] <- data.frame(
      h2 = h2,
      information = information,
      n_pheno = n_pheno,
      total_outer_replicates = total_n,
      RAMN_matched_n = matched_n,
      RAMN_applicability =
        matched_n / total_n,

      method_A =
        ramn_matched_method_label(A),
      method_B =
        ramn_matched_method_label(B),

      coverage_difference_A_minus_B =
        unname(dd["coverage_difference"]),
      coverage_difference_MC_lo =
        unname(dd["coverage_difference_lo"]),
      coverage_difference_MC_hi =
        unname(dd["coverage_difference_hi"]),

      MSE_ratio_difference_A_minus_B =
        unname(dd["MSE_ratio_difference"]),
      MSE_ratio_difference_MC_lo =
        unname(dd["MSE_ratio_difference_lo"]),
      MSE_ratio_difference_MC_hi =
        unname(dd["MSE_ratio_difference_hi"]),

      RMSE_difference_A_minus_B =
        unname(dd["RMSE_difference"]),
      RMSE_difference_MC_lo =
        unname(dd["RMSE_difference_lo"]),
      RMSE_difference_MC_hi =
        unname(dd["RMSE_difference_hi"]),

      stringsAsFactors = FALSE
    )
  }

  pair_table <- do.call(
    rbind,
    pair_rows
  )

  applicability <- data.frame(
    h2 = h2,
    information = information,
    n_pheno = n_pheno,
    total_outer_replicates = total_n,
    RAMN_matched_n = matched_n,
    RAMN_applicability =
      matched_n / total_n,
    RAMN_applicability_MCSE =
      sqrt(
        (matched_n / total_n) *
          (1 - matched_n / total_n) /
          total_n
      ),
    stringsAsFactors = FALSE
  )

  list(
    method_table = method_table,
    pair_table = pair_table,
    applicability = applicability,
    matched_indices = idx
  )
}


# ------------------------------------------------------------
# 4. Find production Simulation-II RDS files
# ------------------------------------------------------------

ramn_matched_find_final_rds <- function(input_dir) {

  if (!dir.exists(input_dir)) {
    stop(
      "Input directory does not exist:\n",
      input_dir
    )
  }

  files <- list.files(
    input_dir,
    pattern = "_final\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(files) == 0L) {
    stop(
      "No *_final.rds files found under:\n",
      input_dir
    )
  }

  keep <- logical(
    length(files)
  )

  for (i in seq_along(files)) {

    x <- try(
      readRDS(files[i]),
      silent = TRUE
    )

    keep[i] <- (
      !inherits(x, "try-error") &&
        !is.null(x$raw) &&
        !is.null(x$summary) &&
        !is.null(x$summary$h2) &&
        !is.null(x$summary$information) &&
        !is.null(x$summary$n_pheno)
    )
  }

  files <- files[keep]

  if (length(files) == 0L) {
    stop(
      "No Simulation-II production final RDS with $raw and $summary found."
    )
  }

  files
}


# ------------------------------------------------------------
# 5. Main runner
# ------------------------------------------------------------

run_RAMN_matched_replicate_comparison <- function(
    input_dir = file.path(
      getwd(),
      "simulation_I_II_final_output",
      "Simulation_II"
    ),
    output_dir = file.path(
      getwd(),
      "RAMN_matched_replicate_comparison"
    ),
    target = 0,
    R_boot = 5000,
    seed = 20260873) {

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  files <- ramn_matched_find_final_rds(
    input_dir
  )

  cat(
    "\n============================================================\n",
    "RAM-N MATCHED-REPLICATE COMPARISON\n",
    "============================================================\n",
    "Input directory: ", normalizePath(input_dir), "\n",
    "Number of scenario files: ", length(files), "\n",
    "Outer-bootstrap replicates: ", R_boot, "\n",
    "============================================================\n\n",
    sep = ""
  )

  all_methods <- list()
  all_pairs <- list()
  all_app <- list()
  matched_index_list <- list()

  for (i in seq_along(files)) {

    obj <- readRDS(
      files[i]
    )

    key <- paste0(
      "h2_",
      obj$summary$h2,
      "__",
      obj$summary$information
    )

    cat(
      "===== ",
      key,
      " | n_pheno = ",
      obj$summary$n_pheno,
      " =====\n",
      sep = ""
    )

    ans <- ramn_matched_analyse_scenario(
      obj = obj,
      file_name = basename(files[i]),
      target = target,
      R_boot = R_boot,
      seed = seed + i * 1000L
    )

    all_methods[[key]] <- ans$method_table
    all_pairs[[key]] <- ans$pair_table
    all_app[[key]] <- ans$applicability
    matched_index_list[[key]] <- ans$matched_indices

    cat(
      "RAM-N valid: ",
      ans$applicability$RAMN_matched_n,
      " / ",
      ans$applicability$total_outer_replicates,
      " = ",
      sprintf(
        "%.3f",
        ans$applicability$RAMN_applicability
      ),
      "\n\n",
      sep = ""
    )
  }

  method_table <- do.call(
    rbind,
    all_methods
  )

  pair_table <- do.call(
    rbind,
    all_pairs
  )

  applicability_table <- do.call(
    rbind,
    all_app
  )

  rownames(method_table) <- NULL
  rownames(pair_table) <- NULL
  rownames(applicability_table) <- NULL

  method_table <- method_table[
    order(
      method_table$h2,
      method_table$n_pheno,
      match(
        method_table$method,
        c(
          "Conditional",
          "RL-UP",
          "RAM-N"
        )
      )
    ),
    ,
    drop = FALSE
  ]

  pair_table <- pair_table[
    order(
      pair_table$h2,
      pair_table$n_pheno,
      pair_table$method_A,
      pair_table$method_B
    ),
    ,
    drop = FALSE
  ]

  applicability_table <- applicability_table[
    order(
      applicability_table$h2,
      applicability_table$n_pheno
    ),
    ,
    drop = FALSE
  ]

  write.csv(
    method_table,
    file.path(
      output_dir,
      "RAMN_matched_method_performance.csv"
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    pair_table,
    file.path(
      output_dir,
      "RAMN_matched_pairwise_differences.csv"
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    applicability_table,
    file.path(
      output_dir,
      "RAMN_matched_applicability.csv"
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  out <- list(
    method_table = method_table,
    pairwise_table = pair_table,
    applicability_table = applicability_table,
    matched_indices = matched_index_list,
    settings = list(
      input_dir = normalizePath(input_dir),
      output_dir = normalizePath(output_dir),
      target = target,
      R_boot = R_boot,
      seed = seed,
      definition =
        "Matched set = outer replicates with non-NULL RAM_N_policy"
    )
  )

  saveRDS(
    out,
    file.path(
      output_dir,
      "RAMN_matched_comparison_full.rds"
    )
  )

  cat(
    "\n============================================================\n",
    "Matched-replicate comparison completed.\n",
    "Output directory:\n",
    normalizePath(output_dir),
    "\n============================================================\n\n",
    sep = ""
  )

  print(
    applicability_table,
    row.names = FALSE
  )

  cat(
    "\n--- Matched method performance ---\n"
  )

  print(
    method_table,
    row.names = FALSE
  )

  cat(
    "\n--- Paired differences ---\n"
  )

  print(
    pair_table,
    row.names = FALSE
  )

  invisible(out)
}


# ------------------------------------------------------------
# 6. Compact supplementary-table helpers
# ------------------------------------------------------------

make_RAMN_matched_supplementary_table <- function(x) {

  tab <- x$method_table

  out <- tab[
    ,
    c(
      "h2",
      "information",
      "n_pheno",
      "RAMN_matched_n",
      "RAMN_applicability",
      "method",
      "bias",
      "RMSE",
      "MSE_ratio",
      "MSE_ratio_MC_lo",
      "MSE_ratio_MC_hi",
      "coverage",
      "coverage_MCSE",
      "width"
    ),
    drop = FALSE
  ]

  numeric_cols <- vapply(
    out,
    is.numeric,
    logical(1)
  )

  out[numeric_cols] <- lapply(
    out[numeric_cols],
    function(z) round(z, 3)
  )

  out
}


make_RAMN_matched_key_pairwise_table <- function(x) {

  tab <- x$pairwise_table

  keep <- (
    tab$method_B == "RAM-N" |
      (
        tab$method_A == "Conditional" &
          tab$method_B == "RL-UP"
      )
  )

  out <- tab[
    keep,
    c(
      "h2",
      "information",
      "n_pheno",
      "RAMN_matched_n",
      "method_A",
      "method_B",
      "coverage_difference_A_minus_B",
      "coverage_difference_MC_lo",
      "coverage_difference_MC_hi",
      "MSE_ratio_difference_A_minus_B",
      "MSE_ratio_difference_MC_lo",
      "MSE_ratio_difference_MC_hi"
    ),
    drop = FALSE
  ]

  numeric_cols <- vapply(
    out,
    is.numeric,
    logical(1)
  )

  out[numeric_cols] <- lapply(
    out[numeric_cols],
    function(z) round(z, 3)
  )

  out
}


# ------------------------------------------------------------
# 7. Self-test using synthetic raw results
# ------------------------------------------------------------

self_test_RAMN_matched_replicate_comparison <- function() {

  set.seed(20260874)

  S <- 20

  res <- vector(
    "list",
    S
  )

  for (s in seq_len(S)) {

    truth <- rnorm(1)

    mk <- function(
        mean_shift = 0,
        sd = 0.5) {

      mu <- truth +
        rnorm(1, mean_shift, 0.4)

      c(
        mean = mu,
        sd = sd,
        lo = mu - 1.96 * sd,
        hi = mu + 1.96 * sd,
        p = pnorm(
          mu / sd
        )
      )
    }

    res[[s]] <- list(
      truth = truth,
      conditional = mk(
        mean_shift = 0.05,
        sd = 0.30
      ),
      Bayesian_policy = mk(
        mean_shift = 0.02,
        sd = 0.45
      ),
      RAM_N_policy =
        if (s %% 4 == 0) {
          NULL
        } else {
          mk(
            mean_shift = 0.01,
            sd = 0.48
          )
        }
    )
  }

  obj <- list(
    summary = list(
      h2 = 0.20,
      information = "test",
      n_pheno = 90
    ),
    raw = res
  )

  ans <- ramn_matched_analyse_scenario(
    obj = obj,
    file_name = "synthetic_test.rds",
    R_boot = 200,
    seed = 20260875
  )

  stopifnot(
    ans$applicability$total_outer_replicates == 20,
    ans$applicability$RAMN_matched_n == 15,
    nrow(ans$method_table) == 3,
    all(
      ans$method_table$RAMN_matched_n == 15
    ),
    all(is.finite(ans$method_table$coverage_Wilson_lo)),
    all(is.finite(ans$method_table$coverage_Wilson_hi)),
    all(is.finite(
      ans$pair_table$coverage_difference_A_minus_B
    )),
    all(is.finite(
      ans$pair_table$MSE_ratio_difference_A_minus_B
    )),
    all(is.finite(
      ans$pair_table$RMSE_difference_A_minus_B
    ))
  )

  cat(
    "\nRAM-N matched-replicate self-test passed.\n"
  )

  print(
    ans$method_table,
    row.names = FALSE
  )

  invisible(ans)
}
