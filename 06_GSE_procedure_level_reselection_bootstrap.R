# ============================================================================
# GSE: PROCEDURE-LEVEL RESELECTION BOOTSTRAP
#
# Uses the completed Simulation-II bootstrap master object. For each inner
# sample it repeats REML -> EBLUP -> top-k selection, samples a full breeding-
# value vector conditional on the same phenotype sample, and evaluates the
# procedure-level basic-bootstrap interval. The stored fixed-policy results are
# treated as immutable references for validation. Base R only.
# ============================================================================

RBOOTM_ANALYSIS_VERSION <- "master_only_error_diagnostics_v2_20260925"


# ============================================================================
# 0. Required validated core and small helpers
# ============================================================================

rbootm_require_core <- function() {
  req <- c(
    "fit_reml_animal",
    "predict_ebv_validation",
    "policy_stats_animal",
    "make_policy",
    "simulate_latent_replicate",
    "make_dataset_from_latent",
    "bootdiag_gls_beta",
    "bootdiag_joint_generator",
    "bootdiag_fixed_eta_interval",
    "Vinv_apply",
    "safe_quantile"
  )

  miss <- req[!vapply(req, exists, logical(1), mode = "function")]

  if (length(miss) > 0L) {
    stop(
      "Required validated core functions are missing: ",
      paste(miss, collapse = ", "),
      "\nSource 01_GSE_main_and_diagnostics.R first."
    )
  }

  invisible(TRUE)
}


rbootm_check_master <- function(boot1000) {
  if (is.null(boot1000) || !is.list(boot1000)) {
    stop("boot1000 must be a list read from Simulation_II_bootstrap_generator_ALL.rds.")
  }

  req <- c(
    "design", "information_design", "result", "method_table",
    "comparison_table", "diagnostic_table", "settings"
  )

  miss <- req[!req %in% names(boot1000)]
  if (length(miss) > 0L) {
    stop("boot1000 is missing: ", paste(miss, collapse = ", "))
  }

  st <- boot1000$settings
  req_st <- c(
    "S", "B", "h2_values", "info_fractions", "n_select",
    "sigma2_P", "beta", "seed_info", "seed_outer"
  )

  miss_st <- req_st[!req_st %in% names(st)]
  if (length(miss_st) > 0L) {
    stop("boot1000$settings is missing: ", paste(miss_st, collapse = ", "))
  }

  if (is.null(boot1000$information_design$info_levels) ||
      is.null(boot1000$information_design$base_design)) {
    stop("boot1000$information_design does not contain the expected saved design.")
  }

  invisible(TRUE)
}


rbootm_safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
}


rbootm_dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(path)) stop("Could not create output directory: ", path)
  invisible(path)
}


rbootm_safe_ratio <- function(num, den) {
  if (!is.finite(num) || !is.finite(den) || den <= 0) return(NA_real_)
  num / den
}


rbootm_top_selected <- function(uhat, candidates, n_select) {
  ord <- order(uhat[candidates], decreasing = TRUE)
  candidates[ord[seq_len(n_select)]]
}


rbootm_make_policy_from_uhat <- function(
    uhat,
    n_animals,
    candidates,
    n_select) {

  selected <- rbootm_top_selected(
    uhat = uhat,
    candidates = candidates,
    n_select = n_select
  )

  cstar <- make_policy(
    n_animals = n_animals,
    candidates = candidates,
    selected = selected
  )

  list(selected = selected, cstar = cstar)
}


rbootm_wilson <- function(x, conf = 0.95) {
  x <- as.logical(x)
  x <- x[!is.na(x)]

  n <- length(x)
  if (n == 0L) {
    return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_, n = 0))
  }

  p <- mean(x)
  z <- qnorm(1 - (1 - conf) / 2)
  den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den

  c(
    estimate = p,
    lower = max(0, ctr - half),
    upper = min(1, ctr + half),
    n = n
  )
}


rbootm_pack_method <- function(x) {
  if (is.null(x) || is.null(x$valid) || !isTRUE(x$valid)) return(NULL)

  c(
    mean = as.numeric(x$mean),
    sd = as.numeric(x$sd),
    lo = as.numeric(x$interval[1]),
    hi = as.numeric(x$interval[2]),
    p = as.numeric(x$p)
  )
}


rbootm_get_stored_method <- function(
    boot1000,
    h2,
    information,
    method) {

  z <- boot1000$method_table
  idx <- which(
    abs(z$h2 - h2) < 1e-12 &
      as.character(z$information) == as.character(information) &
      as.character(z$method) == as.character(method)
  )

  if (length(idx) != 1L) {
    stop(
      "Could not uniquely find stored method row: h2=", h2,
      ", information=", information,
      ", method=", method
    )
  }

  z[idx, , drop = FALSE]
}


rbootm_get_stored_diagnostic <- function(
    boot1000,
    h2,
    information) {

  z <- boot1000$diagnostic_table
  idx <- which(
    abs(z$h2 - h2) < 1e-12 &
      as.character(z$information) == as.character(information)
  )

  if (length(idx) != 1L) {
    stop(
      "Could not uniquely find stored diagnostic row: h2=", h2,
      ", information=", information
    )
  }

  z[idx, , drop = FALSE]
}


# ============================================================================
# 1. EXACT reconstruction of the outer RNG schedule
# ============================================================================
#
# IMPORTANT subtlety:
# The validated production code did
#   set.seed(seed_outer)
#   seed_latent <- sample.int(..., S)
#   latent_list <- lapply(seed_latent, simulate_latent_replicate(seed=...))
#   seed_methods <- sample.int(...)
# in this exact order.
#
# simulate_latent_replicate() resets the RNG internally for each latent sample.
# Therefore seed_methods must be generated AFTER all latent samples have been
# regenerated.  This function deliberately preserves that order.
# ============================================================================

rbootm_rebuild_rng_schedule <- function(boot1000) {
  rbootm_require_core()
  rbootm_check_master(boot1000)

  st <- boot1000$settings
  S <- as.integer(st$S)
  h2_values <- as.numeric(st$h2_values)
  info_bundle <- boot1000$information_design

  set.seed(st$seed_outer)

  seed_latent <- sample.int(
    .Machine$integer.max,
    S
  )

  latent_list <- lapply(
    seed_latent,
    function(ss) {
      simulate_latent_replicate(
        info_bundle,
        seed = ss
      )
    }
  )

  n_h <- length(h2_values)
  n_i <- length(info_bundle$info_levels)

  seed_methods <- array(
    sample.int(
      .Machine$integer.max,
      S * n_h * n_i * 3
    ),
    dim = c(S, n_h, n_i, 3)
  )

  list(
    seed_latent = seed_latent,
    latent_list = latent_list,
    seed_methods = seed_methods
  )
}


# ============================================================================
# 2. Reconstruct one OUTER replicate without running any inner bootstrap
# ============================================================================

rbootm_reconstruct_outer_one <- function(
    latent,
    info_bundle,
    info_level,
    h2,
    sigma2_P,
    beta,
    target,
    n_select) {

  base <- info_bundle$base_design

  dat <- make_dataset_from_latent(
    latent = latent,
    info_level = info_level,
    sigma2_P = sigma2_P,
    h2 = h2,
    beta = beta
  )

  fit_sel <- fit_reml_animal(
    y = dat$y,
    prep = info_level$prep
  )

  pred <- predict_ebv_validation(
    eta = fit_sel$eta,
    y = dat$y,
    prep = info_level$prep
  )

  candidates <- base$candidates

  selected <- rbootm_top_selected(
    uhat = pred$uhat,
    candidates = candidates,
    n_select = n_select
  )

  cstar <- make_policy(
    n_animals = nrow(base$ped),
    candidates = candidates,
    selected = selected
  )

  truth <- sum(cstar * dat$u)

  gs_plugin <- policy_stats_animal(
    eta = fit_sel$eta,
    y = dat$y,
    prep = info_level$prep,
    cstar = cstar,
    calc_diag = FALSE
  )

  conditional_plugin <- list(
    valid = TRUE,
    mean = gs_plugin$mu,
    sd = gs_plugin$sd,
    interval = gs_plugin$mu + c(-1, 1) * 1.96 * gs_plugin$sd,
    p = 1 - pnorm(
      target,
      mean = gs_plugin$mu,
      sd = max(gs_plugin$sd, 1e-12)
    )
  )

  eta_true <- log(c(dat$sigma2_A, dat$sigma2_e))

  conditional_trueVC_oracle <- bootdiag_fixed_eta_interval(
    eta = eta_true,
    y = dat$y,
    prep = info_level$prep,
    cstar = cstar,
    target = target
  )

  list(
    dat = dat,
    fit = fit_sel,
    selected = selected,
    cstar = cstar,
    truth = truth,
    conditional_plugin = rbootm_pack_method(conditional_plugin),
    conditional_trueVC_oracle = rbootm_pack_method(conditional_trueVC_oracle),
    h2_hat = fit_sel$theta[1] / sum(fit_sel$theta),
    lower_boundary_outer = as.numeric(
      fit_sel$eta[1] - fit_sel$lower[1] < 0.05
    )
  )
}


# ============================================================================
# 3. Standard method summary -- same definitions as the validated paper code
# ============================================================================

rbootm_summarize_method <- function(
    res,
    method,
    target = 0) {

  x <- lapply(res, `[[`, method)
  ok <- !vapply(x, is.null, logical(1))
  x_ok <- x[ok]

  truth <- vapply(
    res,
    function(z) as.numeric(z$truth),
    numeric(1)
  )

  if (length(x_ok) == 0L) {
    return(c(
      valid = 0,
      bias = NA_real_,
      RMSE = NA_real_,
      mean_SE = NA_real_,
      MSE_ratio = NA_real_,
      inclusion = NA_real_,
      width = NA_real_,
      Brier = NA_real_,
      mean_p = NA_real_,
      event_rate = mean(truth > target),
      coverage_MCSE = NA_real_,
      coverage_Wilson_lo = NA_real_,
      coverage_Wilson_hi = NA_real_
    ))
  }

  truth_ok <- truth[ok]

  grab <- function(nm) {
    vapply(x_ok, function(z) as.numeric(z[[nm]]), numeric(1))
  }

  mu <- grab("mean")
  se <- grab("sd")
  lo <- grab("lo")
  hi <- grab("hi")
  pp <- grab("p")

  err <- mu - truth_ok
  mse <- mean(err^2)
  inside <- truth_ok >= lo & truth_ok <= hi
  coverage <- mean(inside)
  wi <- rbootm_wilson(inside)

  c(
    valid = length(x_ok),
    bias = mean(err),
    RMSE = sqrt(mse),
    mean_SE = mean(se),
    MSE_ratio = if (mse > 0) mean(se^2) / mse else NA_real_,
    inclusion = coverage,
    width = mean(hi - lo),
    Brier = mean((pp - as.numeric(truth_ok > target))^2),
    mean_p = mean(pp),
    event_rate = mean(truth_ok > target),
    coverage_MCSE = sqrt(coverage * (1 - coverage) / length(x_ok)),
    coverage_Wilson_lo = unname(wi["lower"]),
    coverage_Wilson_hi = unname(wi["upper"])
  )
}


rbootm_summary_row <- function(
    res,
    method,
    h2,
    information,
    n_pheno,
    target = 0) {

  z <- rbootm_summarize_method(
    res = res,
    method = method,
    target = target
  )

  data.frame(
    h2 = h2,
    information = information,
    n_pheno = n_pheno,
    method = method,
    as.list(z),
    row.names = NULL,
    check.names = FALSE
  )
}


# ============================================================================
# 4. Validation: does the master RDS reconstruct the original OUTER samples?
# ============================================================================
#
# This uses NO inner bootstrap.  It recomputes only the outer REML/EBLUP
# analyses and compares the 9-scenario summaries against the immutable rows
# stored in boot1000$method_table and boot1000$diagnostic_table.
#
# If this passes, the same outer samples / selection procedure have been
# reconstructed.  Run this once before the expensive new bootstrap.
# ============================================================================

validate_boot1000_outer_reconstruction <- function(
    boot1000,
    tolerance = 1e-8,
    output_dir = NULL,
    stop_on_failure = TRUE) {

  rbootm_require_core()
  rbootm_check_master(boot1000)

  st <- boot1000$settings
  S <- as.integer(st$S)
  h2_values <- as.numeric(st$h2_values)
  n_select <- as.integer(st$n_select)
  sigma2_P <- as.numeric(st$sigma2_P)
  beta <- as.numeric(st$beta)
  target <- if (!is.null(st$target)) as.numeric(st$target) else 0

  info_bundle <- boot1000$information_design
  rng <- rbootm_rebuild_rng_schedule(boot1000)

  if (is.null(output_dir)) {
    output_dir <- file.path(
      getwd(),
      "GSE_reselection_bootstrap_master_only",
      "outer_reconstruction_validation"
    )
  }
  rbootm_dir_create(output_dir)

  validation_rows <- list()
  reconstructed_rows <- list()
  diagnostic_rows <- list()
  rid <- 1L

  compare_cols <- c(
    "valid", "bias", "RMSE", "mean_SE", "MSE_ratio",
    "inclusion", "width", "Brier", "mean_p", "event_rate"
  )

  cat(
    "\n====================================================\n",
    "VALIDATING OUTER RECONSTRUCTION FROM boot1000\n",
    "====================================================\n",
    "No old scenario RDS files are used.\n",
    "S = ", S, " outer replicates per scenario.\n",
    sep = ""
  )

  for (ih in seq_along(h2_values)) {
    h2 <- h2_values[ih]

    for (ii in seq_along(info_bundle$info_levels)) {
      info <- info_bundle$info_levels[[ii]]

      cat(
        "Validating h2=", h2,
        ", information=", info$name,
        ", n_pheno=", info$n_pheno, " ...\n",
        sep = ""
      )

      res <- vector("list", S)
      h2_hat <- numeric(S)
      boundary <- numeric(S)

      for (s in seq_len(S)) {
        z <- rbootm_reconstruct_outer_one(
          latent = rng$latent_list[[s]],
          info_bundle = info_bundle,
          info_level = info,
          h2 = h2,
          sigma2_P = sigma2_P,
          beta = beta,
          target = target,
          n_select = n_select
        )

        res[[s]] <- list(
          truth = z$truth,
          conditional_plugin = z$conditional_plugin,
          conditional_trueVC_oracle = z$conditional_trueVC_oracle
        )

        h2_hat[s] <- z$h2_hat
        boundary[s] <- z$lower_boundary_outer
      }

      for (m in c("conditional_plugin", "conditional_trueVC_oracle")) {
        rec <- rbootm_summary_row(
          res = res,
          method = m,
          h2 = h2,
          information = info$name,
          n_pheno = info$n_pheno,
          target = target
        )

        old <- rbootm_get_stored_method(
          boot1000 = boot1000,
          h2 = h2,
          information = info$name,
          method = m
        )

        diffs <- vapply(
          compare_cols,
          function(nm) abs(as.numeric(rec[[nm]]) - as.numeric(old[[nm]])),
          numeric(1)
        )

        validation_rows[[rid]] <- data.frame(
          h2 = h2,
          information = info$name,
          method = m,
          max_abs_difference = max(diffs, na.rm = TRUE),
          within_tolerance = max(diffs, na.rm = TRUE) <= tolerance,
          row.names = NULL
        )

        reconstructed_rows[[rid]] <- rec
        rid <- rid + 1L
      }

      old_dg <- rbootm_get_stored_diagnostic(
        boot1000 = boot1000,
        h2 = h2,
        information = info$name
      )

      dg_rec <- data.frame(
        h2 = h2,
        information = info$name,
        n_pheno = info$n_pheno,
        mean_outer_h2_hat = mean(h2_hat),
        median_outer_h2_hat = median(h2_hat),
        outer_h2_bias = mean(h2_hat - h2),
        outer_lower_boundary_rate = mean(boundary),
        row.names = NULL,
        check.names = FALSE
      )

      dg_cols <- c(
        "mean_outer_h2_hat", "median_outer_h2_hat",
        "outer_h2_bias", "outer_lower_boundary_rate"
      )

      dg_diff <- vapply(
        dg_cols,
        function(nm) abs(as.numeric(dg_rec[[nm]]) - as.numeric(old_dg[[nm]])),
        numeric(1)
      )

      diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
        dg_rec,
        max_abs_difference = max(dg_diff, na.rm = TRUE),
        within_tolerance = max(dg_diff, na.rm = TRUE) <= tolerance,
        row.names = NULL,
        check.names = FALSE
      )
    }
  }

  validation_table <- do.call(rbind, validation_rows)
  reconstructed_table <- do.call(rbind, reconstructed_rows)
  diagnostic_validation <- do.call(rbind, diagnostic_rows)
  rownames(validation_table) <- NULL
  rownames(reconstructed_table) <- NULL
  rownames(diagnostic_validation) <- NULL

  all_ok <- all(validation_table$within_tolerance) &&
    all(diagnostic_validation$within_tolerance)

  out <- list(
    all_ok = all_ok,
    tolerance = tolerance,
    method_validation = validation_table,
    reconstructed_method_table = reconstructed_table,
    diagnostic_validation = diagnostic_validation,
    immutable_stored_method_table = boot1000$method_table,
    immutable_stored_diagnostic_table = boot1000$diagnostic_table
  )

  saveRDS(
    out,
    file.path(output_dir, "outer_reconstruction_validation.rds")
  )

  write.csv(
    validation_table,
    file.path(output_dir, "outer_reconstruction_method_validation.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    diagnostic_validation,
    file.path(output_dir, "outer_reconstruction_diagnostic_validation.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  cat("\nMethod-summary validation:\n")
  print(validation_table, row.names = FALSE)

  cat("\nOuter REML diagnostic validation:\n")
  print(diagnostic_validation, row.names = FALSE)

  cat("\nOverall outer reconstruction: ", if (all_ok) "PASS" else "FAIL", "\n", sep = "")

  if (!all_ok && isTRUE(stop_on_failure)) {
    stop(
      "Outer reconstruction did not reproduce boot1000 within tolerance. ",
      "Do not run the full reselection bootstrap until this is resolved."
    )
  }

  invisible(out)
}


# ============================================================================
# 5. Full-u* sampling conditional on the EXACT original plug-in y* draws
# ============================================================================
#
# The original fixed-policy bootstrap generated (y*, G_fixed*) jointly using
# bootdiag_joint_generator().  We reproduce those y* draws exactly.
#
# To evaluate c(Y*)'u* after reselection, the full u* vector is required.
# We sample u* | y* with Matheron's conditional-simulation rule:
#
#   draw u0 ~ N(0,G), e0 ~ N(0,R), y0 = Zu0 + e0,
#   u* = u0 + G Z' V^{-1} (y* - X beta - y0).
#
# The original joint generator adds 1e-12 to the y* diagonal as a numerical
# stabilizer, so the same 1e-12 is added to the residual variance in this
# conditional sampler.  This makes the marginal y* covariance identical to the
# original plug-in generator while avoiding a large conditional-covariance
# eigendecomposition for every outer replicate.
#
# This gives three branches from the same inner y*:
#   A. fixed_existing_recomputed:
#        exact original fixed-policy truth G_fixed* from the joint generator;
#   B. fixed_paired_new:
#        c_obs' u* using the newly sampled full u*;
#   C. reselection:
#        c(Y*)' u* with c(Y*) rebuilt from the same y*.
#
# Branch A is used to verify the stored manuscript bootstrap.  Branches B/C are
# the clean paired comparison for the effect of policy re-formation.
# ============================================================================

parametric_bootstrap_reselection_from_master_one <- function(
    y,
    prep,
    candidates,
    n_select,
    cstar_obs,
    B = 500,
    B_reference = B,
    target = 0,
    current_fit = NULL,
    beta_generator = NULL,
    seed_boot,
    A_chol = NULL,
    tolerance = 1e-8,
    return_inner = FALSE) {

  rbootm_require_core()

  B <- as.integer(B)
  B_reference <- as.integer(B_reference)

  if (B < 20L || B_reference < B) {
    stop("Require 20 <= B <= B_reference.")
  }

  if (is.null(current_fit)) {
    current_fit <- fit_reml_animal(y = y, prep = prep)
  }

  sigma2_A <- as.numeric(current_fit$theta[1])
  sigma2_e <- as.numeric(current_fit$theta[2])

  pred_obs <- predict_ebv_validation(
    eta = current_fit$eta,
    y = y,
    prep = prep
  )

  n_animals <- length(pred_obs$uhat)

  pol_obs <- rbootm_make_policy_from_uhat(
    uhat = pred_obs$uhat,
    n_animals = n_animals,
    candidates = candidates,
    n_select = n_select
  )

  cstar_obs <- as.numeric(cstar_obs)
  if (length(cstar_obs) != n_animals) stop("cstar_obs has wrong length.")

  c_diff <- max(abs(cstar_obs - pol_obs$cstar))
  if (!is.finite(c_diff) || c_diff > tolerance) {
    stop(
      "Reconstructed observed policy does not match cstar_obs. max diff=",
      format(c_diff, digits = 16)
    )
  }

  gs_obs <- policy_stats_animal(
    eta = current_fit$eta,
    y = y,
    prep = prep,
    cstar = cstar_obs,
    calc_diag = FALSE
  )
  current <- as.numeric(gs_obs$mu)

  if (is.null(beta_generator)) {
    beta_generator <- bootdiag_gls_beta(
      y = y,
      prep = prep,
      sigma2_A = sigma2_A,
      sigma2_e = sigma2_e
    )
  }
  beta_generator <- as.numeric(beta_generator)

  # EXACT original plug-in fixed-policy joint generator.
  joint <- bootdiag_joint_generator(
    prep = prep,
    cstar = cstar_obs,
    sigma2_A = sigma2_A,
    sigma2_e = sigma2_e,
    beta_generator = beta_generator
  )

  nobs <- prep$nobs
  d_joint <- nobs + 1L

  set.seed(seed_boot)

  # Generate the full B_reference matrix first, exactly as required for the
  # original B=500 production mapping; a pilot uses its first B rows.
  z_joint_reference <- matrix(
    rnorm(B_reference * d_joint),
    nrow = B_reference,
    ncol = d_joint
  )

  z_joint <- z_joint_reference[seq_len(B), , drop = FALSE]

  draws <- sweep(
    z_joint %*% joint$chol,
    2,
    joint$mean,
    "+"
  )

  Ystar <- draws[, seq_len(nobs), drop = FALSE]
  Gfixed_original_star <- as.numeric(draws[, nobs + 1L])

  # ------------------------------------------------------------------------
  # Conditional full-u* generator by Matheron's rule.
  # The +1e-12 exactly matches the y-block stabilizer in
  # bootdiag_joint_generator().
  # ------------------------------------------------------------------------

  A <- prep$A
  ids <- prep$pheno_ids
  sigma2_e_y <- sigma2_e + 1e-12

  if (is.null(A_chol)) {
    A_chol <- chol((A + t(A)) / 2)
  }

  if (!is.matrix(A_chol) ||
      nrow(A_chol) != n_animals ||
      ncol(A_chol) != n_animals) {
    stop("A_chol has incompatible dimensions.")
  }

  # Continue the RNG stream only AFTER the exact original z_joint_reference
  # has been generated.  Therefore these new draws cannot alter the old y*.
  z_u_reference <- matrix(
    rnorm(B_reference * n_animals),
    nrow = B_reference,
    ncol = n_animals
  )

  z_e_reference <- matrix(
    rnorm(B_reference * nobs),
    nrow = B_reference,
    ncol = nobs
  )

  z_u <- z_u_reference[seq_len(B), , drop = FALSE]
  z_e <- z_e_reference[seq_len(B), , drop = FALSE]

  U0 <- sqrt(sigma2_A) * (z_u %*% A_chol)
  E0 <- sqrt(sigma2_e_y) * z_e
  Y0 <- U0[, ids, drop = FALSE] + E0

  mu_y <- as.numeric(joint$mean[seq_len(nobs)])
  innovation <- sweep(Ystar, 2, mu_y, "-") - Y0

  # Vinv_apply accepts a matrix, so all B conditional corrections are
  # calculated in one spectral solve.
  Vi_innovation <- Vinv_apply(
    prep = prep,
    sigma2_A = sigma2_A,
    sigma2_e = sigma2_e_y,
    B = t(innovation)
  )

  correction <- t(
    sigma2_A *
      (A[, ids, drop = FALSE] %*% Vi_innovation)
  )

  Ustar <- U0 + correction

  # Inner outputs.
  err_fixed_existing <- rep(NA_real_, B)
  err_fixed_paired <- rep(NA_real_, B)
  err_reselect <- rep(NA_real_, B)

  ghat_fixed <- rep(NA_real_, B)
  gtrue_fixed_existing <- Gfixed_original_star
  gtrue_fixed_paired <- rep(NA_real_, B)
  ghat_reselect <- rep(NA_real_, B)
  gtrue_reselect <- rep(NA_real_, B)

  overlap_outer <- rep(NA_real_, B)
  exact_outer <- rep(NA_real_, B)
  h2_hat <- rep(NA_real_, B)
  converged <- rep(NA_real_, B)
  boundary <- rep(NA_real_, B)

  selected_obs <- pol_obs$selected

  selected_store <- if (return_inner) {
    matrix(NA_integer_, nrow = B, ncol = n_select)
  } else {
    NULL
  }

  for (b in seq_len(B)) {
    yb <- as.numeric(Ystar[b, ])
    ub <- as.numeric(Ustar[b, ])

    fit_b <- fit_reml_animal(
      y = yb,
      prep = prep
    )

    pred_b <- predict_ebv_validation(
      eta = fit_b$eta,
      y = yb,
      prep = prep
    )

    # Fixed-policy point estimate -- same calculation as the original bootstrap.
    gs_fixed <- policy_stats_animal(
      eta = fit_b$eta,
      y = yb,
      prep = prep,
      cstar = cstar_obs,
      calc_diag = FALSE
    )

    ghat_fixed[b] <- as.numeric(gs_fixed$mu)

    # A. Exact old fixed-policy bootstrap branch.
    err_fixed_existing[b] <-
      ghat_fixed[b] - gtrue_fixed_existing[b]

    # B. New paired fixed-policy branch using the same full u* as reselection.
    gtrue_fixed_paired[b] <- sum(cstar_obs * ub)
    err_fixed_paired[b] <-
      ghat_fixed[b] - gtrue_fixed_paired[b]

    # C. Procedure-level reselection branch.
    pol_b <- rbootm_make_policy_from_uhat(
      uhat = pred_b$uhat,
      n_animals = n_animals,
      candidates = candidates,
      n_select = n_select
    )

    # The reselection point predictor is exactly c(Y*)' uhat*.  Reuse the
    # already-computed EBLUP vector rather than recomputing policy PEC.
    ghat_reselect[b] <- sum(pol_b$cstar * pred_b$uhat)
    gtrue_reselect[b] <- sum(pol_b$cstar * ub)
    err_reselect[b] <- ghat_reselect[b] - gtrue_reselect[b]

    overlap_outer[b] <-
      length(intersect(pol_b$selected, selected_obs)) / n_select

    exact_outer[b] <- as.numeric(
      setequal(pol_b$selected, selected_obs)
    )

    h2_hat[b] <- fit_b$theta[1] / sum(fit_b$theta)
    converged[b] <- as.numeric(fit_b$convergence == 0)
    boundary[b] <- as.numeric(
      fit_b$eta[1] - fit_b$lower[1] < 0.05
    )

    if (return_inner) {
      selected_store[b, ] <- as.integer(pol_b$selected)
    }
  }

  ok <- is.finite(err_fixed_existing) &
    is.finite(err_fixed_paired) &
    is.finite(err_reselect)

  if (sum(ok) < max(20L, ceiling(0.80 * B))) {
    return(list(
      valid = FALSE,
      reason = "too_few_valid_inner_bootstrap_replicates",
      n_valid = sum(ok),
      B = B
    ))
  }

  make_branch <- function(err) {
    e <- err[ok]
    q <- safe_quantile(e, c(0.025, 0.975))
    interval <- c(current - q[2], current - q[1])
    pseudo <- current - e

    list(
      valid = TRUE,
      mean = current - mean(e),
      sd = sd(e),
      interval = interval,
      p = mean(pseudo > target),
      current = current,
      error_mean = mean(e),
      error_sd = sd(e),
      error_E2 = mean(e^2),
      error_q025 = q[1],
      error_q975 = q[2],
      n_valid = length(e),
      valid_fraction = length(e) / B
    )
  }

  fixed_existing <- make_branch(err_fixed_existing)
  fixed_paired <- make_branch(err_fixed_paired)
  reselect <- make_branch(err_reselect)

  out <- list(
    valid = TRUE,
    fixed_existing_recomputed = fixed_existing,
    fixed_paired_new = fixed_paired,
    reselect = reselect,
    generator_sigma2_A = sigma2_A,
    generator_sigma2_e = sigma2_e,
    generator_h2 = sigma2_A / (sigma2_A + sigma2_e),
    beta_generator = beta_generator,
    conditional_sampler = "Matheron_u_given_y",
    conditional_y_residual_variance = sigma2_e_y,
    mean_overlap_with_outer_policy = mean(overlap_outer[ok], na.rm = TRUE),
    exact_outer_policy_rate = mean(exact_outer[ok], na.rm = TRUE),
    policy_switch_rate = 1 - mean(exact_outer[ok], na.rm = TRUE),
    inner_REML_convergence_rate = mean(converged[ok], na.rm = TRUE),
    inner_boundary_rate = mean(boundary[ok], na.rm = TRUE),
    inner_mean_h2_hat = mean(h2_hat[ok], na.rm = TRUE),
    paired_E2_ratio_reselect_over_fixed = rbootm_safe_ratio(
      reselect$error_E2,
      fixed_paired$error_E2
    ),
    paired_error_correlation = suppressWarnings(
      cor(err_fixed_paired[ok], err_reselect[ok])
    )
  )

  if (return_inner) {
    out$error_fixed_existing <- err_fixed_existing
    out$error_fixed_paired <- err_fixed_paired
    out$error_reselect <- err_reselect
    out$Ystar <- Ystar
    out$Ustar <- Ustar
    out$gtrue_fixed_existing <- gtrue_fixed_existing
    out$gtrue_fixed_paired <- gtrue_fixed_paired
    out$gtrue_reselect <- gtrue_reselect
    out$ghat_fixed <- ghat_fixed
    out$ghat_reselect <- ghat_reselect
    out$selected_reselect <- selected_store
    out$overlap_outer <- overlap_outer
    out$exact_outer <- exact_outer
  }

  out
}


# ============================================================================
# 6. Full 3 x 3 runner using ONLY boot1000 + validated core code
# ============================================================================

run_reselection_bootstrap_from_master <- function(
    boot1000,
    outer_indices = NULL,
    B_inner = NULL,
    output_dir = NULL,
    checkpoint_every = 10,
    resume = TRUE,
    tolerance = 1e-8,
    stop_on_full_validation_failure = TRUE,
    keep_inner_vectors = FALSE,
    keep_per_outer_in_master = FALSE) {

  rbootm_require_core()
  rbootm_check_master(boot1000)

  st <- boot1000$settings
  S <- as.integer(st$S)
  B_reference <- as.integer(st$B)
  h2_values <- as.numeric(st$h2_values)
  n_select <- as.integer(st$n_select)
  sigma2_P <- as.numeric(st$sigma2_P)
  beta <- as.numeric(st$beta)
  target <- if (!is.null(st$target)) as.numeric(st$target) else 0

  if (is.null(B_inner)) {
    B <- B_reference
  } else {
    B <- as.integer(B_inner)
  }

  if (B < 20L || B > B_reference) {
    stop("B_inner must be between 20 and the stored B_reference=", B_reference, ".")
  }

  if (is.null(outer_indices)) {
    outer_indices <- seq_len(S)
  }

  outer_indices <- sort(unique(as.integer(outer_indices)))
  if (length(outer_indices) == 0L ||
      any(outer_indices < 1L) ||
      any(outer_indices > S)) {
    stop("outer_indices must be within 1:S.")
  }

  full_outer <- identical(outer_indices, seq_len(S))
  full_B <- identical(B, B_reference)
  full_production <- full_outer && full_B

  if (is.null(output_dir)) {
    output_dir <- file.path(
      getwd(),
      "GSE_reselection_bootstrap_master_only",
      paste0(
        "S", length(outer_indices),
        "_B", B
      )
    )
  }
  rbootm_dir_create(output_dir)

  info_bundle <- boot1000$information_design
  base <- boot1000$design

  # The master saved both objects.  They should refer to the same base design.
  design_check <- isTRUE(all.equal(
    base$A,
    info_bundle$base_design$A,
    tolerance = 0
  )) && identical(
    as.integer(base$candidates),
    as.integer(info_bundle$base_design$candidates)
  )

  if (!design_check) {
    stop("boot1000$design and boot1000$information_design$base_design disagree.")
  }

  candidates <- base$candidates
  A_chol <- chol((base$A + t(base$A)) / 2)
  rng <- rbootm_rebuild_rng_schedule(boot1000)

  cat(
    "\n====================================================\n",
    "PROCEDURE-LEVEL RESELECTION BOOTSTRAP -- MASTER ONLY\n",
    "====================================================\n",
    "Old scenario RDS files: NOT USED\n",
    "Stored manuscript tables: IMMUTABLE REFERENCE\n",
    "Outer replicates: ", length(outer_indices), " / ", S, "\n",
    "Inner B: ", B, " / ", B_reference, "\n",
    "Full production validation enabled: ", full_production, "\n",
    sep = ""
  )

  method_rows <- list()
  comparison_rows <- list()
  validation_rows <- list()
  per_outer_rows <- list()
  scenario_objects <- list()
  row_id <- 1L

  for (ih in seq_along(h2_values)) {
    h2 <- h2_values[ih]

    for (ii in seq_along(info_bundle$info_levels)) {
      info <- info_bundle$info_levels[[ii]]
      key <- paste0("h2_", h2, "__", info$name)

      cat(
        "\n----------------------------------------------------\n",
        key, " | n_pheno=", info$n_pheno, "\n",
        "----------------------------------------------------\n",
        sep = ""
      )

      scenario_tag <- paste0(
        rbootm_safe_name(key),
        "__nOuter", length(outer_indices),
        "__B", B
      )

      checkpoint_file <- file.path(
        output_dir,
        paste0(scenario_tag, "_checkpoint.rds")
      )

      final_file <- file.path(
        output_dir,
        paste0(scenario_tag, "_final.rds")
      )

      if (resume && file.exists(final_file)) {
        done <- readRDS(final_file)

        compatible <- !is.null(done$analysis_version) &&
          identical(as.character(done$analysis_version), RBOOTM_ANALYSIS_VERSION) &&
          !is.null(done$key) &&
          !is.null(done$B) &&
          !is.null(done$outer_indices) &&
          identical(as.character(done$key), key) &&
          identical(as.integer(done$B), B) &&
          identical(as.integer(done$outer_indices), outer_indices)

        if (isTRUE(compatible)) {
          cat("Loading completed compatible scenario.\n")
          method_rows[[row_id]] <- done$method_table
          comparison_rows[[row_id]] <- done$comparison
          validation_rows[[row_id]] <- done$validation
          if (is.null(done$per_outer)) {
            stop("Compatible final scenario is missing per_outer diagnostics: ", final_file)
          }
          po <- done$per_outer
          po$h2 <- h2
          po$information <- info$name
          po$n_pheno <- info$n_pheno
          per_outer_rows[[row_id]] <- po
          scenario_objects[[key]] <- done$summary_object
          row_id <- row_id + 1L
          next
        }
      }

      res <- vector("list", length(outer_indices))
      start_pos <- 1L

      if (resume && file.exists(checkpoint_file)) {
        cp <- readRDS(checkpoint_file)

        compatible_cp <- !is.null(cp$analysis_version) &&
          identical(as.character(cp$analysis_version), RBOOTM_ANALYSIS_VERSION) &&
          !is.null(cp$key) &&
          !is.null(cp$B) &&
          !is.null(cp$outer_indices) &&
          identical(as.character(cp$key), key) &&
          identical(as.integer(cp$B), B) &&
          identical(as.integer(cp$outer_indices), outer_indices)

        if (isTRUE(compatible_cp)) {
          res <- cp$res
          miss <- which(vapply(res, is.null, logical(1)))
          start_pos <- if (length(miss) == 0L) length(res) + 1L else min(miss)
          cat("Resuming at subset position ", start_pos, ".\n", sep = "")
        }
      }

      if (start_pos <= length(outer_indices)) {
        for (pos in seq.int(start_pos, length(outer_indices))) {
          s <- outer_indices[pos]

          outer <- rbootm_reconstruct_outer_one(
            latent = rng$latent_list[[s]],
            info_bundle = info_bundle,
            info_level = info,
            h2 = h2,
            sigma2_P = sigma2_P,
            beta = beta,
            target = target,
            n_select = n_select
          )

          beta_hat_plugin <- bootdiag_gls_beta(
            y = outer$dat$y,
            prep = info$prep,
            sigma2_A = outer$fit$theta[1],
            sigma2_e = outer$fit$theta[2]
          )

          pb <- parametric_bootstrap_reselection_from_master_one(
            y = outer$dat$y,
            prep = info$prep,
            candidates = candidates,
            n_select = n_select,
            cstar_obs = outer$cstar,
            B = B,
            B_reference = B_reference,
            target = target,
            current_fit = outer$fit,
            beta_generator = beta_hat_plugin,
            seed_boot = rng$seed_methods[s, ih, ii, 2],
            A_chol = A_chol,
            tolerance = tolerance,
            return_inner = keep_inner_vectors
          )

          if (is.null(pb$valid) || !isTRUE(pb$valid)) {
            stop(
              "Inner bootstrap failed at ", key,
              ", outer replicate ", s,
              ". Reason: ",
              if (!is.null(pb$reason)) pb$reason else "unknown"
            )
          }

          one <- list(
            outer_replicate = s,
            truth = outer$truth,
            conditional_plugin = outer$conditional_plugin,
            conditional_trueVC_oracle = outer$conditional_trueVC_oracle,
            bootstrap_fixed_recomputed = rbootm_pack_method(
              pb$fixed_existing_recomputed
            ),
            bootstrap_fixed_paired_new = rbootm_pack_method(
              pb$fixed_paired_new
            ),
            bootstrap_reselection = rbootm_pack_method(
              pb$reselect
            ),
            diagnostics = list(
              h2_hat_outer = outer$h2_hat,
              lower_boundary_outer = outer$lower_boundary_outer,

              # Actual outer prediction error for the realised policy.
              outer_point_predictor = pb$reselect$current,
              outer_prediction_error = pb$reselect$current - outer$truth,

              # Original fixed-policy bootstrap error distribution.
              fixed_existing_error_mean =
                pb$fixed_existing_recomputed$error_mean,
              fixed_existing_error_sd =
                pb$fixed_existing_recomputed$error_sd,
              fixed_existing_error_q025 =
                pb$fixed_existing_recomputed$error_q025,
              fixed_existing_error_q975 =
                pb$fixed_existing_recomputed$error_q975,
              fixed_existing_interval_lo =
                pb$fixed_existing_recomputed$interval[1],
              fixed_existing_interval_hi =
                pb$fixed_existing_recomputed$interval[2],

              # Paired fixed-policy branch using the same full u* draws as
              # the reselection branch.
              fixed_paired_error_mean = pb$fixed_paired_new$error_mean,
              fixed_paired_error_sd = pb$fixed_paired_new$error_sd,
              fixed_paired_error_q025 = pb$fixed_paired_new$error_q025,
              fixed_paired_error_q975 = pb$fixed_paired_new$error_q975,
              fixed_paired_interval_lo = pb$fixed_paired_new$interval[1],
              fixed_paired_interval_hi = pb$fixed_paired_new$interval[2],

              # Procedure-level reselection bootstrap error distribution.
              reselect_error_mean = pb$reselect$error_mean,
              reselect_error_sd = pb$reselect$error_sd,
              reselect_error_q025 = pb$reselect$error_q025,
              reselect_error_q975 = pb$reselect$error_q975,
              reselect_interval_lo = pb$reselect$interval[1],
              reselect_interval_hi = pb$reselect$interval[2],
              reselect_error_mean_minus_outer_error =
                pb$reselect$error_mean -
                (pb$reselect$current - outer$truth),

              policy_switch_rate = pb$policy_switch_rate,
              mean_overlap_with_outer_policy = pb$mean_overlap_with_outer_policy,
              exact_outer_policy_rate = pb$exact_outer_policy_rate,
              paired_E2_ratio_reselect_over_fixed =
                pb$paired_E2_ratio_reselect_over_fixed,
              paired_error_correlation = pb$paired_error_correlation,
              inner_REML_convergence_rate = pb$inner_REML_convergence_rate,
              inner_boundary_rate = pb$inner_boundary_rate,
              inner_mean_h2_hat = pb$inner_mean_h2_hat,
              conditional_y_residual_variance =
                pb$conditional_y_residual_variance
            )
          )

          if (keep_inner_vectors) {
            one$inner <- pb
          }

          res[[pos]] <- one

          if (
            pos %% checkpoint_every == 0L ||
            pos == length(outer_indices)
          ) {
            saveRDS(
              list(
                analysis_version = RBOOTM_ANALYSIS_VERSION,
                key = key,
                B = B,
                outer_indices = outer_indices,
                res = res
              ),
              checkpoint_file
            )

            cat(
              "Completed ", pos, " / ", length(outer_indices),
              " (outer replicate ", s, ")\n",
              sep = ""
            )
          }
        }
      }

      if (any(vapply(res, is.null, logical(1)))) {
        stop("Incomplete scenario: ", key)
      }

      # New/current-run summaries.
      method_names <- c(
        "conditional_plugin",
        "conditional_trueVC_oracle",
        "bootstrap_fixed_recomputed",
        "bootstrap_fixed_paired_new",
        "bootstrap_reselection"
      )

      mt <- do.call(
        rbind,
        lapply(
          method_names,
          function(m) {
            rbootm_summary_row(
              res = res,
              method = m,
              h2 = h2,
              information = info$name,
              n_pheno = info$n_pheno,
              target = target
            )
          }
        )
      )
      rownames(mt) <- NULL

      old_fixed <- rbootm_get_stored_method(
        boot1000 = boot1000,
        h2 = h2,
        information = info$name,
        method = "bootstrap_plugin_generator"
      )

      get_new <- function(m) {
        mt[mt$method == m, , drop = FALSE]
      }

      rec_fixed <- get_new("bootstrap_fixed_recomputed")
      paired_fixed <- get_new("bootstrap_fixed_paired_new")
      reselect <- get_new("bootstrap_reselection")

      # Paired outer coverage gain: new paired-fixed vs reselection.
      cov_fixed_pair <- vapply(
        res,
        function(z) {
          q <- z$bootstrap_fixed_paired_new
          as.numeric(z$truth >= q[["lo"]] && z$truth <= q[["hi"]])
        },
        numeric(1)
      )

      cov_re <- vapply(
        res,
        function(z) {
          q <- z$bootstrap_reselection
          as.numeric(z$truth >= q[["lo"]] && z$truth <= q[["hi"]])
        },
        numeric(1)
      )

      delta_cov_pair <- cov_re - cov_fixed_pair
      gain <- mean(delta_cov_pair)
      gain_se <- if (length(delta_cov_pair) >= 2L) {
        sd(delta_cov_pair) / sqrt(length(delta_cov_pair))
      } else {
        NA_real_
      }

      diag_mean <- function(nm) {
        mean(
          vapply(
            res,
            function(z) as.numeric(z$diagnostics[[nm]]),
            numeric(1)
          ),
          na.rm = TRUE
        )
      }

      diag_mean_abs <- function(nm) {
        mean(
          abs(vapply(
            res,
            function(z) as.numeric(z$diagnostics[[nm]]),
            numeric(1)
          )),
          na.rm = TRUE
        )
      }

      comparison <- data.frame(
        h2 = h2,
        information = info$name,
        n_pheno = info$n_pheno,
        n_outer = length(outer_indices),
        B_inner = B,
        full_production = full_production,

        stored_fixed_full_coverage = as.numeric(old_fixed$inclusion),
        stored_fixed_full_MSE_ratio = as.numeric(old_fixed$MSE_ratio),
        stored_fixed_full_bias = as.numeric(old_fixed$bias),
        stored_fixed_full_RMSE = as.numeric(old_fixed$RMSE),
        stored_fixed_full_mean_SE = as.numeric(old_fixed$mean_SE),
        stored_fixed_full_width = as.numeric(old_fixed$width),

        recomputed_fixed_coverage = rec_fixed$inclusion,
        recomputed_fixed_MSE_ratio = rec_fixed$MSE_ratio,
        recomputed_fixed_bias = rec_fixed$bias,
        recomputed_fixed_RMSE = rec_fixed$RMSE,
        recomputed_fixed_mean_SE = rec_fixed$mean_SE,
        recomputed_fixed_width = rec_fixed$width,

        paired_fixed_new_coverage = paired_fixed$inclusion,
        paired_fixed_new_MSE_ratio = paired_fixed$MSE_ratio,

        reselect_coverage = reselect$inclusion,
        reselect_coverage_MCSE = reselect$coverage_MCSE,
        reselect_coverage_Wilson_lo = reselect$coverage_Wilson_lo,
        reselect_coverage_Wilson_hi = reselect$coverage_Wilson_hi,
        reselect_MSE_ratio = reselect$MSE_ratio,
        reselect_bias = reselect$bias,
        reselect_RMSE = reselect$RMSE,
        reselect_mean_SE = reselect$mean_SE,
        reselect_width = reselect$width,

        paired_coverage_gain_reselect_vs_new_fixed = gain,
        paired_coverage_gain_SE = gain_se,
        paired_coverage_gain_lo = if (is.finite(gain_se)) gain - 1.96 * gain_se else NA_real_,
        paired_coverage_gain_hi = if (is.finite(gain_se)) gain + 1.96 * gain_se else NA_real_,

        # Error-distribution diagnostics.  These are descriptive summaries
        # across outer replicates; per-outer values are retained below.
        mean_outer_prediction_error = diag_mean("outer_prediction_error"),
        mean_fixed_existing_boot_error_mean =
          diag_mean("fixed_existing_error_mean"),
        mean_fixed_paired_boot_error_mean =
          diag_mean("fixed_paired_error_mean"),
        mean_reselect_boot_error_mean = diag_mean("reselect_error_mean"),
        mean_abs_reselect_boot_error_mean =
          diag_mean_abs("reselect_error_mean"),
        mean_reselect_boot_error_sd = diag_mean("reselect_error_sd"),
        mean_reselect_boot_error_q025 = diag_mean("reselect_error_q025"),
        mean_reselect_boot_error_q975 = diag_mean("reselect_error_q975"),
        mean_reselect_error_mean_minus_outer_error =
          diag_mean("reselect_error_mean_minus_outer_error"),
        mean_abs_reselect_error_mean_minus_outer_error =
          diag_mean_abs("reselect_error_mean_minus_outer_error"),

        mean_policy_switch_rate = diag_mean("policy_switch_rate"),
        mean_overlap_with_outer_policy = diag_mean("mean_overlap_with_outer_policy"),
        mean_paired_E2_ratio_reselect_over_fixed =
          diag_mean("paired_E2_ratio_reselect_over_fixed"),
        mean_paired_error_correlation = diag_mean("paired_error_correlation"),
        mean_inner_REML_convergence_rate = diag_mean("inner_REML_convergence_rate"),
        mean_inner_boundary_rate = diag_mean("inner_boundary_rate"),
        mean_inner_h2_hat = diag_mean("inner_mean_h2_hat"),
        row.names = NULL,
        check.names = FALSE
      )

      # ----------------------------------------------------------
      # Validation against IMMUTABLE stored manuscript summaries.
      # Exact comparison is required only for full S and full B.
      # ----------------------------------------------------------

      validation <- data.frame(
        h2 = numeric(),
        information = character(),
        check = character(),
        max_abs_difference = numeric(),
        within_tolerance = logical(),
        row.names = NULL,
        check.names = FALSE
      )

      if (full_outer) {
        compare_cols_outer <- c(
          "valid", "bias", "RMSE", "mean_SE", "MSE_ratio",
          "inclusion", "width", "Brier", "mean_p", "event_rate"
        )

        for (m in c("conditional_plugin", "conditional_trueVC_oracle")) {
          rec <- get_new(m)
          old <- rbootm_get_stored_method(
            boot1000, h2, info$name, m
          )

          d <- vapply(
            compare_cols_outer,
            function(nm) abs(as.numeric(rec[[nm]]) - as.numeric(old[[nm]])),
            numeric(1)
          )

          validation <- rbind(
            validation,
            data.frame(
              h2 = h2,
              information = info$name,
              check = paste0("outer_", m),
              max_abs_difference = max(d, na.rm = TRUE),
              within_tolerance = max(d, na.rm = TRUE) <= tolerance,
              row.names = NULL
            )
          )
        }
      }

      if (full_production) {
        compare_cols_boot <- c(
          "valid", "bias", "RMSE", "mean_SE", "MSE_ratio",
          "inclusion", "width", "Brier", "mean_p", "event_rate"
        )

        dboot <- vapply(
          compare_cols_boot,
          function(nm) {
            old_name <- nm
            abs(as.numeric(rec_fixed[[nm]]) - as.numeric(old_fixed[[old_name]]))
          },
          numeric(1)
        )

        validation <- rbind(
          validation,
          data.frame(
            h2 = h2,
            information = info$name,
            check = "fixed_bootstrap_reproduces_stored_manuscript_result",
            max_abs_difference = max(dboot, na.rm = TRUE),
            within_tolerance = max(dboot, na.rm = TRUE) <= tolerance,
            row.names = NULL
          )
        )
      }

      if (nrow(validation) == 0L) {
        validation <- data.frame(
          h2 = h2,
          information = info$name,
          check = "pilot_run_no_exact_stored_summary_check",
          max_abs_difference = NA_real_,
          within_tolerance = NA,
          row.names = NULL
        )
      }

      if (full_production &&
          any(validation$within_tolerance %in% FALSE) &&
          isTRUE(stop_on_full_validation_failure)) {
        saveRDS(
          list(
            method_table = mt,
            comparison = comparison,
            validation = validation,
            res = res
          ),
          file.path(output_dir, paste0(scenario_tag, "_VALIDATION_FAILED.rds"))
        )

        stop(
          "Full-production validation failed for ", key,
          ". Existing manuscript results must NOT be replaced. ",
          "Inspect the saved validation-failed RDS."
        )
      }

      per_outer <- data.frame(
        outer_replicate = vapply(res, `[[`, numeric(1), "outer_replicate"),
        truth = vapply(res, `[[`, numeric(1), "truth"),

        outer_point_predictor = vapply(
          res, function(z) z$diagnostics$outer_point_predictor, numeric(1)
        ),
        outer_prediction_error = vapply(
          res, function(z) z$diagnostics$outer_prediction_error, numeric(1)
        ),

        fixed_existing_error_mean = vapply(
          res, function(z) z$diagnostics$fixed_existing_error_mean, numeric(1)
        ),
        fixed_existing_error_sd = vapply(
          res, function(z) z$diagnostics$fixed_existing_error_sd, numeric(1)
        ),
        fixed_existing_error_q025 = vapply(
          res, function(z) z$diagnostics$fixed_existing_error_q025, numeric(1)
        ),
        fixed_existing_error_q975 = vapply(
          res, function(z) z$diagnostics$fixed_existing_error_q975, numeric(1)
        ),
        fixed_existing_interval_lo = vapply(
          res, function(z) z$diagnostics$fixed_existing_interval_lo, numeric(1)
        ),
        fixed_existing_interval_hi = vapply(
          res, function(z) z$diagnostics$fixed_existing_interval_hi, numeric(1)
        ),

        fixed_paired_error_mean = vapply(
          res, function(z) z$diagnostics$fixed_paired_error_mean, numeric(1)
        ),
        fixed_paired_error_sd = vapply(
          res, function(z) z$diagnostics$fixed_paired_error_sd, numeric(1)
        ),
        fixed_paired_error_q025 = vapply(
          res, function(z) z$diagnostics$fixed_paired_error_q025, numeric(1)
        ),
        fixed_paired_error_q975 = vapply(
          res, function(z) z$diagnostics$fixed_paired_error_q975, numeric(1)
        ),
        fixed_paired_interval_lo = vapply(
          res, function(z) z$diagnostics$fixed_paired_interval_lo, numeric(1)
        ),
        fixed_paired_interval_hi = vapply(
          res, function(z) z$diagnostics$fixed_paired_interval_hi, numeric(1)
        ),

        reselect_error_mean = vapply(
          res, function(z) z$diagnostics$reselect_error_mean, numeric(1)
        ),
        reselect_error_sd = vapply(
          res, function(z) z$diagnostics$reselect_error_sd, numeric(1)
        ),
        reselect_error_q025 = vapply(
          res, function(z) z$diagnostics$reselect_error_q025, numeric(1)
        ),
        reselect_error_q975 = vapply(
          res, function(z) z$diagnostics$reselect_error_q975, numeric(1)
        ),
        reselect_interval_lo = vapply(
          res, function(z) z$diagnostics$reselect_interval_lo, numeric(1)
        ),
        reselect_interval_hi = vapply(
          res, function(z) z$diagnostics$reselect_interval_hi, numeric(1)
        ),
        reselect_error_mean_minus_outer_error = vapply(
          res,
          function(z) z$diagnostics$reselect_error_mean_minus_outer_error,
          numeric(1)
        ),

        fixed_recomputed_cover = vapply(
          res,
          function(z) {
            q <- z$bootstrap_fixed_recomputed
            as.numeric(z$truth >= q[["lo"]] && z$truth <= q[["hi"]])
          },
          numeric(1)
        ),
        fixed_paired_new_cover = cov_fixed_pair,
        reselect_cover = cov_re,
        policy_switch_rate = vapply(
          res,
          function(z) z$diagnostics$policy_switch_rate,
          numeric(1)
        ),
        overlap_with_outer_policy = vapply(
          res,
          function(z) z$diagnostics$mean_overlap_with_outer_policy,
          numeric(1)
        ),
        paired_E2_ratio_reselect_over_fixed = vapply(
          res,
          function(z) z$diagnostics$paired_E2_ratio_reselect_over_fixed,
          numeric(1)
        ),
        paired_error_correlation = vapply(
          res, function(z) z$diagnostics$paired_error_correlation, numeric(1)
        ),
        inner_REML_convergence_rate = vapply(
          res,
          function(z) z$diagnostics$inner_REML_convergence_rate,
          numeric(1)
        ),
        inner_boundary_rate = vapply(
          res, function(z) z$diagnostics$inner_boundary_rate, numeric(1)
        ),
        inner_mean_h2_hat = vapply(
          res, function(z) z$diagnostics$inner_mean_h2_hat, numeric(1)
        ),
        row.names = NULL,
        check.names = FALSE
      )

      summary_object <- list(
        h2 = h2,
        information = info$name,
        n_pheno = info$n_pheno,
        outer_indices = outer_indices,
        B = B,
        method_table = mt,
        comparison = comparison,
        validation = validation
      )

      saveRDS(
        list(
          analysis_version = RBOOTM_ANALYSIS_VERSION,
          key = key,
          B = B,
          outer_indices = outer_indices,
          method_table = mt,
          comparison = comparison,
          validation = validation,
          per_outer = per_outer,
          summary_object = summary_object,
          res = if (keep_per_outer_in_master || keep_inner_vectors) res else NULL
        ),
        final_file
      )

      write.csv(
        mt,
        file.path(output_dir, paste0(scenario_tag, "_methods.csv")),
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )

      write.csv(
        comparison,
        file.path(output_dir, paste0(scenario_tag, "_comparison.csv")),
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )

      write.csv(
        validation,
        file.path(output_dir, paste0(scenario_tag, "_validation.csv")),
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )

      write.csv(
        per_outer,
        file.path(output_dir, paste0(scenario_tag, "_per_outer.csv")),
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )

      method_rows[[row_id]] <- mt
      comparison_rows[[row_id]] <- comparison
      validation_rows[[row_id]] <- validation
      po_master <- per_outer
      po_master$h2 <- h2
      po_master$information <- info$name
      po_master$n_pheno <- info$n_pheno
      per_outer_rows[[row_id]] <- po_master
      scenario_objects[[key]] <- summary_object

      if (keep_per_outer_in_master) {
        scenario_objects[[key]]$per_outer <- per_outer
      }

      row_id <- row_id + 1L
    }
  }

  method_table <- do.call(rbind, method_rows)
  comparison_table <- do.call(rbind, comparison_rows)
  validation_table <- do.call(rbind, validation_rows)
  per_outer_diagnostics_table <- do.call(rbind, per_outer_rows)
  rownames(method_table) <- NULL
  rownames(comparison_table) <- NULL
  rownames(validation_table) <- NULL
  rownames(per_outer_diagnostics_table) <- NULL

  # Put scenario identifiers first in the master per-outer diagnostic table.
  lead_cols <- c("h2", "information", "n_pheno", "outer_replicate")
  other_cols <- setdiff(names(per_outer_diagnostics_table), lead_cols)
  per_outer_diagnostics_table <- per_outer_diagnostics_table[
    , c(lead_cols, other_cols), drop = FALSE
  ]

  write.csv(
    method_table,
    file.path(output_dir, "reselection_bootstrap_new_method_table.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    comparison_table,
    file.path(output_dir, "reselection_bootstrap_main_comparison.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    validation_table,
    file.path(output_dir, "reselection_bootstrap_validation.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    per_outer_diagnostics_table,
    file.path(output_dir, "reselection_bootstrap_error_diagnostics_per_outer.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  out <- list(
    settings = list(
      analysis_version = RBOOTM_ANALYSIS_VERSION,
      S_reference = S,
      B_reference = B_reference,
      outer_indices = outer_indices,
      B_inner = B,
      full_outer = full_outer,
      full_B = full_B,
      full_production = full_production,
      h2_values = h2_values,
      n_select = n_select,
      sigma2_P = sigma2_P,
      beta = beta,
      target = target,
      seed_info = st$seed_info,
      seed_outer = st$seed_outer,
      source_master_estimand = st$estimand,
      output_dir = normalizePath(
        output_dir,
        winslash = "/",
        mustWork = FALSE
      ),
      new_estimand = paste0(
        "procedure-level G_true(Y)=c(Y)'u; inner bootstrap repeats REML, ",
        "EBLUP and top-k policy reconstruction"
      ),
      note = paste0(
        "Existing boot1000 manuscript summaries remain immutable; ",
        "bootstrap_fixed_recomputed is used only as a reproduction check."
      )
    ),
    immutable_stored_method_table = boot1000$method_table,
    immutable_stored_comparison_table = boot1000$comparison_table,
    immutable_stored_diagnostic_table = boot1000$diagnostic_table,
    new_method_table = method_table,
    comparison_table = comparison_table,
    validation_table = validation_table,
    per_outer_diagnostics_table = per_outer_diagnostics_table,
    scenarios = scenario_objects
  )

  saveRDS(
    out,
    file.path(output_dir, "reselection_bootstrap_master_only_ALL.rds")
  )

  cat(
    "\n====================================================\n",
    "RESELECTION BOOTSTRAP COMPLETED\n",
    "====================================================\n",
    sep = ""
  )

  cat("\nMain comparison:\n")
  print(comparison_table, row.names = FALSE)

  cat("\nValidation:\n")
  print(validation_table, row.names = FALSE)

  invisible(out)
}


# ============================================================================
# 7. Convenience printer
# ============================================================================

print_reselection_bootstrap_from_master <- function(x, digits = 3) {
  if (is.null(x$comparison_table)) stop("Object has no $comparison_table.")

  z <- x$comparison_table

  keep <- c(
    "h2", "information", "n_pheno",
    "stored_fixed_full_coverage",
    "recomputed_fixed_coverage",
    "paired_fixed_new_coverage",
    "reselect_coverage",
    "reselect_coverage_Wilson_lo",
    "reselect_coverage_Wilson_hi",
    "stored_fixed_full_MSE_ratio",
    "recomputed_fixed_MSE_ratio",
    "reselect_MSE_ratio",
    "paired_coverage_gain_reselect_vs_new_fixed",
    "mean_policy_switch_rate",
    "mean_paired_E2_ratio_reselect_over_fixed"
  )

  keep <- keep[keep %in% names(z)]
  zz <- z[, keep, drop = FALSE]

  num <- vapply(zz, is.numeric, logical(1))
  zz[num] <- lapply(zz[num], round, digits = digits)

  print(zz, row.names = FALSE)
  invisible(zz)
}


# ============================================================================
# 8. Small smoke test of the NEW inner machinery
# ============================================================================

self_test_reselection_bootstrap_master_only <- function(
    B = 30,
    seed_data = 2026092501,
    seed_boot = 2026092502) {

  rbootm_require_core()

  if (!exists("build_simII_design", mode = "function") ||
      !exists("simulate_animal_data", mode = "function")) {
    stop("The validated self-test helpers are not available.")
  }

  design <- build_simII_design()

  dat <- simulate_animal_data(
    design = design,
    h2 = 0.20,
    sigma2_P = 1,
    beta = 0,
    seed = seed_data
  )

  fit <- fit_reml_animal(
    y = dat$y,
    prep = design$prep
  )

  pred <- predict_ebv_validation(
    eta = fit$eta,
    y = dat$y,
    prep = design$prep
  )

  selected <- rbootm_top_selected(
    pred$uhat,
    design$candidates,
    5
  )

  cstar <- make_policy(
    n_animals = nrow(design$A),
    candidates = design$candidates,
    selected = selected
  )

  z <- parametric_bootstrap_reselection_from_master_one(
    y = dat$y,
    prep = design$prep,
    candidates = design$candidates,
    n_select = 5,
    cstar_obs = cstar,
    B = B,
    B_reference = B,
    target = 0,
    current_fit = fit,
    seed_boot = seed_boot,
    return_inner = TRUE
  )

  stopifnot(
    isTRUE(z$valid),
    all(is.finite(z$fixed_existing_recomputed$interval)),
    all(is.finite(z$fixed_paired_new$interval)),
    all(is.finite(z$reselect$interval)),
    is.finite(z$paired_E2_ratio_reselect_over_fixed)
  )

  cat("\nMaster-only reselection-bootstrap smoke test passed.\n")
  cat("Recomputed old fixed SD: ", z$fixed_existing_recomputed$error_sd, "\n", sep = "")
  cat("Paired fixed SD:         ", z$fixed_paired_new$error_sd, "\n", sep = "")
  cat("Reselection SD:          ", z$reselect$error_sd, "\n", sep = "")
  cat("Policy switch rate:      ", z$policy_switch_rate, "\n", sep = "")
  cat("E2 ratio reselect/fixed: ", z$paired_E2_ratio_reselect_over_fixed, "\n", sep = "")

  invisible(z)
}


# ============================================================================
# 9. Post-processing helpers for manuscript Supplementary tables
# ============================================================================
# These helpers use only the saved per-outer diagnostic table.  They do not
# rerun REML, EBLUP, selection, or bootstrap simulation and do not modify any
# production result.  They reproduce the summaries used for Supplementary
# Tables S10-S11 in the final manuscript.
# ============================================================================

rbootm_get_per_outer_diagnostics <- function(x) {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) stop("File not found: ", x)
    if (grepl("\\.rds$", x, ignore.case = TRUE)) {
      x <- readRDS(x)
    } else {
      x <- read.csv(x, stringsAsFactors = FALSE, check.names = FALSE)
    }
  }

  if (is.data.frame(x)) return(x)
  if (is.list(x) && !is.null(x$per_outer_diagnostics_table)) {
    return(x$per_outer_diagnostics_table)
  }
  stop(
    "Expected a data.frame, a final reselection-bootstrap RDS object, ",
    "or a CSV/RDS path containing per-outer diagnostics."
  )
}

rbootm_make_supplementary_S10 <- function(x) {
  d <- rbootm_get_per_outer_diagnostics(x)
  req <- c(
    "h2", "information", "n_pheno", "truth",
    "reselect_interval_lo", "reselect_interval_hi",
    "reselect_cover", "outer_prediction_error", "reselect_error_mean"
  )
  miss <- req[!req %in% names(d)]
  if (length(miss) > 0L) {
    stop("Missing S10 diagnostic column(s): ", paste(miss, collapse = ", "))
  }

  keys <- unique(d[, c("h2", "information", "n_pheno"), drop = FALSE])
  keys <- keys[order(keys$h2, keys$n_pheno), , drop = FALSE]

  rows <- lapply(seq_len(nrow(keys)), function(i) {
    z <- d[
      abs(d$h2 - keys$h2[i]) < 1e-12 &
        d$information == keys$information[i] &
        d$n_pheno == keys$n_pheno[i],
      , drop = FALSE
    ]

    inside <- as.logical(z$reselect_cover)
    wi <- rbootm_wilson(inside)
    miss_below <- mean(z$truth < z$reselect_interval_lo, na.rm = TRUE)
    miss_above <- mean(z$truth > z$reselect_interval_hi, na.rm = TRUE)
    corr <- suppressWarnings(
      cor(z$outer_prediction_error, z$reselect_error_mean, use = "complete.obs")
    )

    data.frame(
      h2 = keys$h2[i],
      information = keys$information[i],
      n_pheno = keys$n_pheno[i],
      n_outer = sum(!is.na(inside)),
      coverage = mean(inside, na.rm = TRUE),
      coverage_Wilson_lo = unname(wi["lower"]),
      coverage_Wilson_hi = unname(wi["upper"]),
      truth_below_lower = miss_below,
      truth_above_upper = miss_above,
      corr_outer_error_inner_mean_error = corr,
      row.names = NULL,
      check.names = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

rbootm_make_supplementary_S11 <- function(x) {
  d <- rbootm_get_per_outer_diagnostics(x)
  req <- c(
    "h2", "inner_boundary_rate", "reselect_cover",
    "reselect_error_mean_minus_outer_error"
  )
  miss <- req[!req %in% names(d)]
  if (length(miss) > 0L) {
    stop("Missing S11 diagnostic column(s): ", paste(miss, collapse = ", "))
  }

  d$boundary_stratum <- cut(
    d$inner_boundary_rate,
    breaks = c(-Inf, 0.05, 0.10, 0.20, 0.40, Inf),
    labels = c("<=0.05", "0.05-0.10", "0.10-0.20", "0.20-0.40", ">0.40"),
    right = TRUE,
    include.lowest = TRUE
  )

  h2_values <- sort(unique(d$h2))
  strata <- levels(d$boundary_stratum)
  rows <- list()
  rr <- 1L

  for (hh in h2_values) {
    for (ss in strata) {
      z <- d[abs(d$h2 - hh) < 1e-12 & d$boundary_stratum == ss, , drop = FALSE]
      if (nrow(z) == 0L) next

      rows[[rr]] <- data.frame(
        h2 = hh,
        inner_boundary_rate_stratum = ss,
        n = nrow(z),
        coverage = mean(z$reselect_cover, na.rm = TRUE),
        mean_inner_error_minus_outer_error = mean(
          z$reselect_error_mean_minus_outer_error,
          na.rm = TRUE
        ),
        row.names = NULL,
        check.names = FALSE
      )
      rr <- rr + 1L
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

rbootm_write_manuscript_diagnostic_tables <- function(
    x,
    output_dir = "paper_output/reselection_bootstrap") {

  rbootm_dir_create(output_dir)
  s10 <- rbootm_make_supplementary_S10(x)
  s11 <- rbootm_make_supplementary_S11(x)

  write.csv(
    s10,
    file.path(output_dir, "Supplementary_Table_S10_reselection_diagnostics.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  write.csv(
    s11,
    file.path(output_dir, "Supplementary_Table_S11_boundary_strata.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  invisible(list(S10 = s10, S11 = s11))
}
