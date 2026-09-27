# ============================================================================
# GSE: RL-UP PRIOR-SENSITIVITY ANALYSIS
#
# Evaluates sensitivity of RL-UP to the proper prior used with the restricted
# likelihood in Simulation II. Outer data, realised policy, and random numbers
# are shared across priors within each outer replicate.
#
# Main prior:
#   h2 ~ Beta(1,1)
#   psi = log(sigma_A^2 + sigma_e^2) ~ N(0,1^2)
# Sensitivity priors for h2: Beta(2,2) and Beta(0.5,0.5).
#
# Integration is performed on eta_A = log(sigma_A^2), eta_e = log(sigma_e^2).
# The Jacobian for the (h2, psi) -> (eta_A, eta_e) transformation contributes
# h2*(1-h2). The RNG order is intentionally fixed as
# seed_latent -> latent_list -> seed_methods and must not be rearranged.
# ============================================================================

# Source the public core automatically when this file is run from the release
# directory.  If the core is already loaded, no re-source is performed.
.gse_rlup_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.gse_rlup_dir <- if (!is.null(.gse_rlup_script)) {
  dirname(normalizePath(.gse_rlup_script, winslash = "/", mustWork = FALSE))
} else {
  getwd()
}
if (!exists("build_simII_design", mode = "function")) {
  .gse_core <- file.path(.gse_rlup_dir, "01_GSE_main_and_diagnostics.R")
  if (!file.exists(.gse_core)) {
    stop("Cannot find 01_GSE_main_and_diagnostics.R. Source it first or run from the release directory.")
  }
  source(.gse_core, echo = FALSE, chdir = FALSE, encoding = "UTF-8")
}

# ---------------------------------------------------------------------------
# 0. Dependency check
# ---------------------------------------------------------------------------

required_public_functions <- c(
  "build_simII_design",
  "build_information_level_designs",
  "simulate_latent_replicate",
  "make_dataset_from_latent",
  "fit_reml_animal",
  "predict_ebv_validation",
  "make_policy",
  "policy_stats_animal",
  "log_reml_eta",
  "safe_quantile"
)

check_RLUP_public_core <- function() {
  miss <- required_public_functions[
    !vapply(required_public_functions, exists, logical(1), mode = "function")
  ]

  if (length(miss) > 0L) {
    stop(
      paste0(
        "Required functions are missing: ",
        paste(miss, collapse = ", "),
        "\nSource 01_GSE_main_and_diagnostics.R first."
      )
    )
  }

  invisible(TRUE)
}


# ---------------------------------------------------------------------------
# 1. Prior definitions
# ---------------------------------------------------------------------------

default_RLUP_priors <- function(include_scale_wide = FALSE) {

  out <- list(
    beta11 = list(
      label = "h2 ~ Beta(1,1); log(sigmaP2) ~ N(0,1^2)",
      alpha = 1,
      beta = 1,
      logP_mean = 0,
      logP_sd = 1
    ),

    beta22 = list(
      label = "h2 ~ Beta(2,2); log(sigmaP2) ~ N(0,1^2)",
      alpha = 2,
      beta = 2,
      logP_mean = 0,
      logP_sd = 1
    ),

    beta05 = list(
      label = "h2 ~ Beta(0.5,0.5); log(sigmaP2) ~ N(0,1^2)",
      alpha = 0.5,
      beta = 0.5,
      logP_mean = 0,
      logP_sd = 1
    )
  )

  if (isTRUE(include_scale_wide)) {
    out$beta11_scale2 <- list(
      label = "h2 ~ Beta(1,1); log(sigmaP2) ~ N(0,2^2)",
      alpha = 1,
      beta = 1,
      logP_mean = 0,
      logP_sd = 2
    )
  }

  out
}


print_RLUP_priors <- function(priors = default_RLUP_priors()) {

  rows <- lapply(names(priors), function(nm) {
    p <- priors[[nm]]
    q <- qbeta(
      c(0.025, 0.50, 0.975),
      shape1 = p$alpha,
      shape2 = p$beta
    )

    data.frame(
      prior = nm,
      label = p$label,
      h2_prior_mean = p$alpha / (p$alpha + p$beta),
      h2_q025 = q[1],
      h2_median = q[2],
      h2_q975 = q[3],
      log_sigmaP2_mean = p$logP_mean,
      log_sigmaP2_sd = p$logP_sd,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })

  ans <- do.call(rbind, rows)
  print(ans, row.names = FALSE)
  invisible(ans)
}


# ---------------------------------------------------------------------------
# 2. Stable transformations and log prior on eta grid
# ---------------------------------------------------------------------------

logsumexp2_RLUP <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}


h2_from_eta_RLUP <- function(eta_A, eta_e) {
  xi <- eta_A - eta_e
  out <- numeric(length(xi))
  pos <- xi >= 0

  out[pos] <- 1 / (1 + exp(-xi[pos]))
  ex <- exp(xi[!pos])
  out[!pos] <- ex / (1 + ex)

  pmin(pmax(out, 1e-12), 1 - 1e-12)
}


log_prior_RLUP_eta <- function(eta_A, eta_e, prior) {

  h2 <- h2_from_eta_RLUP(eta_A, eta_e)
  psi <- logsumexp2_RLUP(eta_A, eta_e)

  lp_h2 <- dbeta(
    h2,
    shape1 = prior$alpha,
    shape2 = prior$beta,
    log = TRUE
  )

  lp_psi <- dnorm(
    psi,
    mean = prior$logP_mean,
    sd = prior$logP_sd,
    log = TRUE
  )

  # Jacobian for (eta_A, eta_e) -> (h2, psi)
  log_jac <- log(h2) + log1p(-h2)

  lp_h2 + lp_psi + log_jac
}


# ---------------------------------------------------------------------------
# 3. Common adaptive restricted-likelihood x prior grid
# ---------------------------------------------------------------------------

normalize_log_weights_RLUP <- function(z) {
  # Robust log-sum-exp normalization.  Some extreme variance-component
  # grid points can have log posterior = -Inf (or occasionally NA from
  # numerical evaluation).  Such points have zero posterior mass and must
  # not contaminate the normalization of otherwise valid grid points.
  z <- as.numeric(z)

  # If +Inf ever occurs, concentrate mass uniformly on those points.
  pos_inf <- is.infinite(z) & z > 0
  if (any(pos_inf)) {
    w <- numeric(length(z))
    w[pos_inf] <- 1 / sum(pos_inf)
    return(w)
  }

  # Treat NA/NaN/-Inf as invalid zero-mass grid points.
  ok <- is.finite(z)
  if (!any(ok)) {
    stop(
      "RL-UP weight normalization failed: no finite posterior log-weights. " ,
      "This indicates that the adaptive grid contains no numerically valid " ,
      "posterior point for the current outer replicate/prior."
    )
  }

  m <- max(z[ok])
  w <- numeric(length(z))
  w[ok] <- exp(z[ok] - m)

  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    stop(
      "RL-UP weight normalization failed: non-finite or zero weight sum."
    )
  }

  w <- w / sw
  # Remove round-off and force exact normalization.
  w[!is.finite(w) | w < 0] <- 0
  w <- w / sum(w)
  w
}


edge_mass_RLUP <- function(grid, w) {
  a1 <- range(grid$eta_A)
  a2 <- range(grid$eta_e)

  edge <- (
    abs(grid$eta_A - a1[1]) < 1e-12 |
    abs(grid$eta_A - a1[2]) < 1e-12 |
    abs(grid$eta_e - a2[1]) < 1e-12 |
    abs(grid$eta_e - a2[2]) < 1e-12
  )

  sum(w[edge])
}


build_RLUP_prior_common_grid <- function(
    y,
    prep,
    priors = default_RLUP_priors(),
    fit = NULL,
    ncoarse = 31,
    ngrid = 41,
    log_drop = 14,
    pad_coarse_steps = 2) {

  check_RLUP_public_core()

  if (is.null(fit)) {
    fit <- fit_reml_animal(y, prep)
  }

  prior_names <- names(priors)

  # Coarse grid over the same numerical bounds used by REML.
  g1c <- seq(fit$lower[1], fit$upper[1], length.out = ncoarse)
  g2c <- seq(fit$lower[2], fit$upper[2], length.out = ncoarse)

  coarse <- expand.grid(
    eta_A = g1c,
    eta_e = g2c,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  coarse$logL <- vapply(
    seq_len(nrow(coarse)),
    function(i) {
      log_reml_eta(
        c(coarse$eta_A[i], coarse$eta_e[i]),
        y,
        prep
      )
    },
    numeric(1)
  )

  keep_union <- rep(FALSE, nrow(coarse))

  for (nm in prior_names) {
    lp <- log_prior_RLUP_eta(
      coarse$eta_A,
      coarse$eta_e,
      priors[[nm]]
    )

    z <- coarse$logL + lp
    z[!is.finite(z)] <- -Inf
    mx <- max(z)

    if (!is.finite(mx)) {
      stop(
        "No finite coarse-grid posterior point for prior: ", nm
      )
    }

    keep_union <- keep_union | (z >= mx - log_drop)
  }

  if (!any(keep_union)) {
    stop("No posterior region retained on the coarse grid.")
  }

  d1 <- diff(g1c)[1]
  d2 <- diff(g2c)[1]

  lo1 <- max(
    fit$lower[1],
    min(coarse$eta_A[keep_union]) - pad_coarse_steps * d1
  )
  hi1 <- min(
    fit$upper[1],
    max(coarse$eta_A[keep_union]) + pad_coarse_steps * d1
  )
  lo2 <- max(
    fit$lower[2],
    min(coarse$eta_e[keep_union]) - pad_coarse_steps * d2
  )
  hi2 <- min(
    fit$upper[2],
    max(coarse$eta_e[keep_union]) + pad_coarse_steps * d2
  )

  g1 <- seq(lo1, hi1, length.out = ngrid)
  g2 <- seq(lo2, hi2, length.out = ngrid)

  grid <- expand.grid(
    eta_A = g1,
    eta_e = g2,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  grid$logL <- vapply(
    seq_len(nrow(grid)),
    function(i) {
      log_reml_eta(
        c(grid$eta_A[i], grid$eta_e[i]),
        y,
        prep
      )
    },
    numeric(1)
  )

  grid$h2 <- h2_from_eta_RLUP(grid$eta_A, grid$eta_e)
  grid$psi <- logsumexp2_RLUP(grid$eta_A, grid$eta_e)
  grid$sigmaP2 <- exp(grid$psi)

  W <- matrix(
    NA_real_,
    nrow = nrow(grid),
    ncol = length(prior_names),
    dimnames = list(NULL, prior_names)
  )

  edge <- setNames(numeric(length(prior_names)), prior_names)

  for (nm in prior_names) {
    lp <- log_prior_RLUP_eta(
      grid$eta_A,
      grid$eta_e,
      priors[[nm]]
    )

    W[, nm] <- normalize_log_weights_RLUP(grid$logL + lp)
    edge[nm] <- edge_mass_RLUP(grid, W[, nm])
  }

  list(
    fit = fit,
    grid = grid,
    weights = W,
    edge_mass = edge,
    bounds = c(eta_A_lo = lo1, eta_A_hi = hi1, eta_e_lo = lo2, eta_e_hi = hi2),
    ncoarse = ncoarse,
    ngrid = ngrid,
    log_drop = log_drop
  )
}


# ---------------------------------------------------------------------------
# 4. Common-random-number sampling from discrete grid weights
# ---------------------------------------------------------------------------

weighted_index_from_uniform_RLUP <- function(w, u) {
  # Defensive inverse-CDF sampler.  Invalid/negative weights are assigned
  # zero mass; the remaining weights are renormalized before cumsum().
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0

  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    stop(
      "RL-UP discrete sampling failed: all grid weights are invalid or zero."
    )
  }

  w <- w / sw
  cw <- cumsum(w)
  # Protect findInterval() from tiny floating-point non-monotonicity.
  cw <- cummax(cw)
  cw[length(cw)] <- 1

  if (any(!is.finite(cw)) || is.unsorted(cw, strictly = FALSE)) {
    stop(
      "RL-UP discrete sampling failed: cumulative weights are not finite/ordered."
    )
  }

  idx <- findInterval(u, c(0, cw), rightmost.closed = TRUE)
  pmin(pmax(idx, 1L), length(w))
}


# ---------------------------------------------------------------------------
# 5. RL-UP for one realised policy under several priors
# ---------------------------------------------------------------------------

RLUP_prior_one_policy <- function(
    y,
    prep,
    cstar,
    true_h2,
    priors = default_RLUP_priors(),
    M = 2000,
    ncoarse = 31,
    ngrid = 41,
    log_drop = 14,
    target = 0,
    seed_inner = NULL) {

  check_RLUP_public_core()

  if (!is.null(seed_inner)) {
    set.seed(seed_inner)
  }

  fit <- fit_reml_animal(y, prep)

  bg <- build_RLUP_prior_common_grid(
    y = y,
    prep = prep,
    priors = priors,
    fit = fit,
    ncoarse = ncoarse,
    ngrid = ngrid,
    log_drop = log_drop
  )

  prior_names <- names(priors)

  # Common random numbers across priors isolate prior effects.
  u_mix <- runif(M)
  z_pred <- rnorm(M)

  ids_by_prior <- lapply(prior_names, function(nm) {
    weighted_index_from_uniform_RLUP(bg$weights[, nm], u_mix)
  })
  names(ids_by_prior) <- prior_names

  used <- sort(unique(unlist(ids_by_prior, use.names = FALSE)))

  cache_mu <- setNames(numeric(length(used)), as.character(used))
  cache_v <- setNames(numeric(length(used)), as.character(used))

  for (ii in used) {
    eta <- c(bg$grid$eta_A[ii], bg$grid$eta_e[ii])
    ps <- policy_stats_animal(eta, y, prep, cstar)
    cache_mu[as.character(ii)] <- ps$mu
    cache_v[as.character(ii)] <- ps$v
  }

  out <- vector("list", length(prior_names))
  names(out) <- prior_names

  for (nm in prior_names) {
    id <- ids_by_prior[[nm]]

    mu <- unname(cache_mu[as.character(id)])
    v <- unname(cache_v[as.character(id)])

    delta <- mu + sqrt(pmax(v, 0)) * z_pred
    h2_draw <- bg$grid$h2[id]
    sigmaP2_draw <- bg$grid$sigmaP2[id]

    q_delta <- safe_quantile(delta, c(0.025, 0.975))
    q_h2 <- safe_quantile(h2_draw, c(0.025, 0.975))

    out[[nm]] <- list(
      valid = TRUE,
      mean = mean(delta),
      sd = sd(delta),
      lo = q_delta[1],
      hi = q_delta[2],
      p = mean(delta > target),
      posterior_h2_mean = mean(h2_draw),
      posterior_h2_median = median(h2_draw),
      posterior_h2_lo = q_h2[1],
      posterior_h2_hi = q_h2[2],
      posterior_h2_width = diff(q_h2),
      posterior_h2_contains_true = as.numeric(
        true_h2 >= q_h2[1] && true_h2 <= q_h2[2]
      ),
      posterior_sigmaP2_mean = mean(sigmaP2_draw),
      edge_mass = unname(bg$edge_mass[nm])
    )
  }

  list(
    fit = fit,
    posterior = out,
    grid = bg
  )
}


# ---------------------------------------------------------------------------
# 6. One outer scenario: same REML/EBLUP selection for every prior
# ---------------------------------------------------------------------------

one_RLUP_prior_sensitivity_scenario <- function(
    latent,
    info_bundle,
    info_level,
    h2,
    sigma2_P = 1,
    beta = 0,
    target = 0,
    n_select = 5,
    priors = default_RLUP_priors(),
    M = 2000,
    ncoarse = 31,
    ngrid = 41,
    log_drop = 14,
    seed_inner = NULL) {

  check_RLUP_public_core()

  base <- info_bundle$base_design

  dat <- make_dataset_from_latent(
    latent = latent,
    info_level = info_level,
    sigma2_P = sigma2_P,
    h2 = h2,
    beta = beta
  )

  # The realised policy is determined once from the ordinary REML-EBLUP fit.
  fit_sel <- fit_reml_animal(dat$y, info_level$prep)
  pred <- predict_ebv_validation(fit_sel$eta, dat$y, info_level$prep)

  candidates <- base$candidates
  ord <- order(pred$uhat[candidates], decreasing = TRUE)
  selected <- candidates[ord[seq_len(n_select)]]

  cstar <- make_policy(
    n_animals = nrow(base$ped),
    candidates = candidates,
    selected = selected
  )

  truth <- sum(cstar * dat$u)

  rz <- RLUP_prior_one_policy(
    y = dat$y,
    prep = info_level$prep,
    cstar = cstar,
    true_h2 = h2,
    priors = priors,
    M = M,
    ncoarse = ncoarse,
    ngrid = ngrid,
    log_drop = log_drop,
    target = target,
    seed_inner = seed_inner
  )

  out <- list(
    truth = truth,
    h2 = h2,
    information = info_level$name,
    n_pheno = info_level$n_pheno,
    selected = selected,
    reml_eta = fit_sel$eta,
    reml_theta = fit_sel$theta
  )

  for (nm in names(priors)) {
    x <- rz$posterior[[nm]]

    out[[nm]] <- c(
      mean = x$mean,
      sd = x$sd,
      lo = x$lo,
      hi = x$hi,
      p = x$p
    )

    out[[paste0(nm, "_h2")]] <- c(
      mean = x$posterior_h2_mean,
      median = x$posterior_h2_median,
      lo = x$posterior_h2_lo,
      hi = x$posterior_h2_hi,
      width = x$posterior_h2_width,
      contains_true = x$posterior_h2_contains_true,
      sigmaP2_mean = x$posterior_sigmaP2_mean,
      edge_mass = x$edge_mass
    )
  }

  out
}


# ---------------------------------------------------------------------------
# 7. Summary helpers
# ---------------------------------------------------------------------------

summarize_RLUP_prior_method <- function(res, prior_name, target = 0) {

  x <- lapply(res, `[[`, prior_name)
  ok <- !vapply(x, is.null, logical(1))

  if (!any(ok)) {
    return(c(
      valid = 0,
      bias = NA,
      RMSE = NA,
      mean_SE = NA,
      MSE_ratio = NA,
      coverage = NA,
      coverage_MCSE = NA,
      width = NA,
      mean_p = NA,
      event_rate = NA
    ))
  }

  x <- x[ok]
  truth <- vapply(res[ok], `[[`, numeric(1), "truth")
  mat <- do.call(rbind, x)

  err <- mat[, "mean"] - truth
  mc_mse <- mean(err^2)
  mean_est_mse <- mean(mat[, "sd"]^2)
  coverage <- mean(truth >= mat[, "lo"] & truth <= mat[, "hi"])
  n <- length(truth)

  c(
    valid = n,
    bias = mean(err),
    RMSE = sqrt(mc_mse),
    mean_SE = mean(mat[, "sd"]),
    MSE_ratio = mean_est_mse / mc_mse,
    coverage = coverage,
    coverage_MCSE = sqrt(coverage * (1 - coverage) / n),
    width = mean(mat[, "hi"] - mat[, "lo"]),
    mean_p = mean(mat[, "p"]),
    event_rate = mean(truth > target)
  )
}


summarize_RLUP_prior_h2 <- function(res, prior_name, true_h2) {

  H <- do.call(
    rbind,
    lapply(res, `[[`, paste0(prior_name, "_h2"))
  )

  err <- H[, "mean"] - true_h2

  c(
    posterior_h2_mean = mean(H[, "mean"]),
    posterior_h2_median = mean(H[, "median"]),
    h2_bias = mean(err),
    h2_RMSE = sqrt(mean(err^2)),
    h2_interval_coverage = mean(H[, "contains_true"]),
    mean_h2_interval_width = mean(H[, "width"]),
    mean_posterior_sigmaP2 = mean(H[, "sigmaP2_mean"]),
    mean_edge_mass = mean(H[, "edge_mass"]),
    max_edge_mass = max(H[, "edge_mass"])
  )
}


scenario_prior_table_RLUP <- function(
    res,
    h2,
    information,
    n_pheno,
    priors,
    target = 0) {

  rows <- lapply(names(priors), function(nm) {
    m <- summarize_RLUP_prior_method(res, nm, target = target)
    h <- summarize_RLUP_prior_h2(res, nm, true_h2 = h2)

    data.frame(
      h2 = h2,
      information = information,
      n_pheno = n_pheno,
      prior = nm,
      prior_label = priors[[nm]]$label,
      valid = unname(m["valid"]),
      bias = unname(m["bias"]),
      RMSE = unname(m["RMSE"]),
      mean_SE = unname(m["mean_SE"]),
      MSE_ratio = unname(m["MSE_ratio"]),
      coverage = unname(m["coverage"]),
      coverage_MCSE = unname(m["coverage_MCSE"]),
      width = unname(m["width"]),
      posterior_h2_mean = unname(h["posterior_h2_mean"]),
      h2_bias = unname(h["h2_bias"]),
      h2_RMSE = unname(h["h2_RMSE"]),
      h2_interval_coverage = unname(h["h2_interval_coverage"]),
      mean_h2_interval_width = unname(h["mean_h2_interval_width"]),
      mean_posterior_sigmaP2 = unname(h["mean_posterior_sigmaP2"]),
      mean_edge_mass = unname(h["mean_edge_mass"]),
      max_edge_mass = unname(h["max_edge_mass"]),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })

  do.call(rbind, rows)
}


paired_prior_difference_RLUP <- function(
    res,
    prior_name,
    baseline = "beta11") {

  if (prior_name == baseline) {
    return(c(
      mean_point_difference = 0,
      mean_abs_point_difference = 0,
      mean_width_difference = 0,
      coverage_difference = 0,
      MSE_ratio_difference = 0,
      mean_h2_difference = 0,
      mean_abs_h2_difference = 0
    ))
  }

  truth <- vapply(res, `[[`, numeric(1), "truth")
  A <- do.call(rbind, lapply(res, `[[`, baseline))
  B <- do.call(rbind, lapply(res, `[[`, prior_name))
  HA <- do.call(rbind, lapply(res, `[[`, paste0(baseline, "_h2")))
  HB <- do.call(rbind, lapply(res, `[[`, paste0(prior_name, "_h2")))

  mse_ratio_A <- mean(A[, "sd"]^2) / mean((A[, "mean"] - truth)^2)
  mse_ratio_B <- mean(B[, "sd"]^2) / mean((B[, "mean"] - truth)^2)

  cov_A <- mean(truth >= A[, "lo"] & truth <= A[, "hi"])
  cov_B <- mean(truth >= B[, "lo"] & truth <= B[, "hi"])

  d_point <- B[, "mean"] - A[, "mean"]
  d_h2 <- HB[, "mean"] - HA[, "mean"]

  c(
    mean_point_difference = mean(d_point),
    mean_abs_point_difference = mean(abs(d_point)),
    mean_width_difference = mean((B[, "hi"] - B[, "lo"]) - (A[, "hi"] - A[, "lo"])),
    coverage_difference = cov_B - cov_A,
    MSE_ratio_difference = mse_ratio_B - mse_ratio_A,
    mean_h2_difference = mean(d_h2),
    mean_abs_h2_difference = mean(abs(d_h2))
  )
}


# ---------------------------------------------------------------------------
# 8. Outer-replicate bootstrap Monte Carlo intervals
# ---------------------------------------------------------------------------

metric_RLUP_from_indices <- function(res, prior_name, idx) {

  z <- res[idx]
  truth <- vapply(z, `[[`, numeric(1), "truth")
  A <- do.call(rbind, lapply(z, `[[`, prior_name))

  err <- A[, "mean"] - truth

  c(
    coverage = mean(truth >= A[, "lo"] & truth <= A[, "hi"]),
    MSE_ratio = mean(A[, "sd"]^2) / mean(err^2),
    bias = mean(err),
    RMSE = sqrt(mean(err^2)),
    width = mean(A[, "hi"] - A[, "lo"])
  )
}


bootstrap_MC_RLUP_scenario <- function(
    res,
    priors,
    baseline = "beta11",
    R_boot = 2000,
    seed = 20260890) {

  set.seed(seed)

  n <- length(res)
  prior_names <- names(priors)

  point <- lapply(prior_names, function(nm) {
    metric_RLUP_from_indices(res, nm, seq_len(n))
  })
  names(point) <- prior_names

  boot <- lapply(prior_names, function(nm) {
    matrix(NA_real_, nrow = R_boot, ncol = 5,
           dimnames = list(NULL, names(point[[nm]])))
  })
  names(boot) <- prior_names

  diff_boot <- lapply(prior_names, function(nm) {
    matrix(NA_real_, nrow = R_boot, ncol = 5,
           dimnames = list(NULL, c(
             "coverage_difference",
             "MSE_ratio_difference",
             "bias_difference",
             "RMSE_difference",
             "width_difference"
           )))
  })
  names(diff_boot) <- prior_names

  for (b in seq_len(R_boot)) {
    idx <- sample.int(n, n, replace = TRUE)

    met <- lapply(prior_names, function(nm) {
      metric_RLUP_from_indices(res, nm, idx)
    })
    names(met) <- prior_names

    for (nm in prior_names) {
      boot[[nm]][b, ] <- met[[nm]]
      diff_boot[[nm]][b, ] <- c(
        met[[nm]]["coverage"] - met[[baseline]]["coverage"],
        met[[nm]]["MSE_ratio"] - met[[baseline]]["MSE_ratio"],
        met[[nm]]["bias"] - met[[baseline]]["bias"],
        met[[nm]]["RMSE"] - met[[baseline]]["RMSE"],
        met[[nm]]["width"] - met[[baseline]]["width"]
      )
    }
  }

  rows <- list()
  k <- 1L

  for (nm in prior_names) {
    for (q in colnames(boot[[nm]])) {
      ci <- safe_quantile(boot[[nm]][, q], c(0.025, 0.975))
      rows[[k]] <- data.frame(
        prior = nm,
        quantity = q,
        estimate = unname(point[[nm]][q]),
        lo = ci[1],
        hi = ci[2],
        comparison = "absolute",
        stringsAsFactors = FALSE,
        row.names = NULL
      )
      k <- k + 1L
    }

    dpoint <- c(
      coverage_difference = point[[nm]]["coverage"] - point[[baseline]]["coverage"],
      MSE_ratio_difference = point[[nm]]["MSE_ratio"] - point[[baseline]]["MSE_ratio"],
      bias_difference = point[[nm]]["bias"] - point[[baseline]]["bias"],
      RMSE_difference = point[[nm]]["RMSE"] - point[[baseline]]["RMSE"],
      width_difference = point[[nm]]["width"] - point[[baseline]]["width"]
    )

    for (q in colnames(diff_boot[[nm]])) {
      ci <- safe_quantile(diff_boot[[nm]][, q], c(0.025, 0.975))
      rows[[k]] <- data.frame(
        prior = nm,
        quantity = q,
        estimate = unname(dpoint[q]),
        lo = ci[1],
        hi = ci[2],
        comparison = paste0("difference_vs_", baseline),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
      k <- k + 1L
    }
  }

  do.call(rbind, rows)
}


# ---------------------------------------------------------------------------
# 9. File helpers
# ---------------------------------------------------------------------------

RLUP_dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}


RLUP_safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}


write_csv_RLUP <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(path)
}


# ---------------------------------------------------------------------------
# 10. Main 3 x 3 prior-sensitivity runner with checkpoint/resume
# ---------------------------------------------------------------------------

run_RLUP_prior_sensitivity_final <- function(
    S = 1000,
    h2_values = c(0.05, 0.20, 0.40),
    info_fractions = c(low = 0.40, moderate = 0.70, high = 1.00),
    n_select = 5,
    sigma2_P = 1,
    beta = 0,
    target = 0,
    priors = default_RLUP_priors(),
    M = 2000,
    ncoarse = 31,
    ngrid = 41,
    log_drop = 14,
    checkpoint_every = 25,
    R_boot = 2000,
    baseline = "beta11",
    seed_info = 20260850,
    seed_outer = 20260851,
    seed_bootstrap = 20260890,
    output_dir = file.path(getwd(), "RLUP_prior_sensitivity_S1000"),
    keep_raw = TRUE,
    resume = TRUE) {

  check_RLUP_public_core()
  RLUP_dir_create(output_dir)

  if (!baseline %in% names(priors)) {
    stop("baseline prior is not present in priors.")
  }

  base <- build_simII_design()
  info_bundle <- build_information_level_designs(
    base_design = base,
    info_fractions = info_fractions,
    seed_info = seed_info
  )

  cat("\n============================================================\n")
  cat("RL-UP PRIOR SENSITIVITY: realised EBLUP-selection policy\n")
  cat("============================================================\n")
  cat("Outer replicates:", S, "\n")
  cat("Mixture draws per prior:", M, "\n")
  cat("Grid:", ncoarse, "x", ncoarse, "coarse ->", ngrid, "x", ngrid, "fine\n")
  cat("Baseline prior:", baseline, "\n\n")
  print_RLUP_priors(priors)

  # -----------------------------------------------------------------------
  # Reproduce the SAME outer latent seeds as the manuscript's main runner.
  # We also reconstruct the main runner's three method-seed array and use
  # the third component (Bayesian/RL-UP seed) as the common inner seed.
  # -----------------------------------------------------------------------

  set.seed(seed_outer)
  seed_latent <- sample.int(.Machine$integer.max, S)

  # FINAL production RNG order: latent replicates are materialised BEFORE
  # drawing the method-seed array.  Do not move seed_methods above this block.
  latent_list <- lapply(seed_latent, function(s) {
    simulate_latent_replicate(info_bundle, seed = s)
  })

  n_h <- length(h2_values)
  n_i <- length(info_bundle$info_levels)

  seed_methods <- array(
    sample.int(.Machine$integer.max, S * n_h * n_i * 3),
    dim = c(S, n_h, n_i, 3)
  )

  all_results <- list()
  all_tables <- list()
  all_diff_tables <- list()
  all_mc_tables <- list()

  scenario_counter <- 0L

  for (ih in seq_along(h2_values)) {
    h2 <- h2_values[ih]

    for (ii in seq_along(info_bundle$info_levels)) {
      info <- info_bundle$info_levels[[ii]]
      scenario_counter <- scenario_counter + 1L

      key <- paste0("h2_", h2, "__", info$name)

      cat("\n===== ", key, " | n_pheno = ", info$n_pheno, " =====\n", sep = "")

      checkpoint_file <- file.path(
        output_dir,
        paste0(RLUP_safe_name(key), "_checkpoint.rds")
      )

      final_file <- file.path(
        output_dir,
        paste0(RLUP_safe_name(key), "_final.rds")
      )

      res <- vector("list", S)

      scenario_signature <- list(
        S = S,
        key = key,
        prior_names = names(priors),
        prior_parameters = priors,
        M = M,
        ncoarse = ncoarse,
        ngrid = ngrid,
        log_drop = log_drop,
        n_select = n_select,
        sigma2_P = sigma2_P,
        beta = beta,
        target = target,
        seed_info = seed_info,
        seed_outer = seed_outer,
        rng_order = "latent_list_before_seed_methods_v1"
      )

      if (isTRUE(resume) && file.exists(final_file)) {
        done <- readRDS(final_file)

        if (!is.null(done$signature) && identical(done$signature, scenario_signature)) {
          cat("Loading completed scenario from final RDS.\n")
          res <- done$res
        } else {
          cat("Existing final RDS has a different configuration; recomputing scenario.\n")
        }
      }

      if (any(vapply(res, is.null, logical(1)))) {

        start_at <- 1L

        if (isTRUE(resume) && file.exists(checkpoint_file)) {
          cp <- readRDS(checkpoint_file)

          if (!is.null(cp$signature) && identical(cp$signature, scenario_signature)) {
            res <- cp$res
            miss <- which(vapply(res, is.null, logical(1)))
            if (length(miss) == 0L) {
              start_at <- S + 1L
            } else {
              start_at <- min(miss)
            }
            cat("Resuming from replicate", start_at, "\n")
          } else {
            cat("Existing checkpoint has a different configuration; ignoring it.\n")
          }
        }

        if (start_at <= S) {
          for (s in seq.int(start_at, S)) {

            res[[s]] <- tryCatch(
              one_RLUP_prior_sensitivity_scenario(
                latent = latent_list[[s]],
                info_bundle = info_bundle,
                info_level = info,
                h2 = h2,
                sigma2_P = sigma2_P,
                beta = beta,
                target = target,
                n_select = n_select,
                priors = priors,
                M = M,
                ncoarse = ncoarse,
                ngrid = ngrid,
                log_drop = log_drop,
                seed_inner = seed_methods[s, ih, ii, 3]
              ),
              error = function(e) {
                stop(
                  "RL-UP prior sensitivity failed at scenario ", key,
                  ", outer replicate ", s, ": ",
                  conditionMessage(e),
                  call. = FALSE
                )
              }
            )

            if (
              s %% checkpoint_every == 0L ||
              s == S
            ) {
              saveRDS(
                list(
                  signature = scenario_signature,
                  res = res
                ),
                checkpoint_file
              )
              cat("Completed", s, "of", S, "\n")
            }
          }
        }

        saveRDS(
          list(
            signature = scenario_signature,
            res = res
          ),
          final_file
        )
      }

      # Guard against incomplete checkpoint files.
      if (any(vapply(res, is.null, logical(1)))) {
        stop("Scenario contains incomplete outer replicates: ", key)
      }

      tab <- scenario_prior_table_RLUP(
        res = res,
        h2 = h2,
        information = info$name,
        n_pheno = info$n_pheno,
        priors = priors,
        target = target
      )

      diff_rows <- lapply(names(priors), function(nm) {
        d <- paired_prior_difference_RLUP(
          res = res,
          prior_name = nm,
          baseline = baseline
        )

        data.frame(
          h2 = h2,
          information = info$name,
          n_pheno = info$n_pheno,
          prior = nm,
          baseline = baseline,
          t(d),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      })
      diff_tab <- do.call(rbind, diff_rows)

      mc_tab <- bootstrap_MC_RLUP_scenario(
        res = res,
        priors = priors,
        baseline = baseline,
        R_boot = R_boot,
        seed = seed_bootstrap + scenario_counter
      )
      mc_tab$h2 <- h2
      mc_tab$information <- info$name
      mc_tab$n_pheno <- info$n_pheno
      mc_tab <- mc_tab[, c(
        "h2", "information", "n_pheno", "prior", "comparison",
        "quantity", "estimate", "lo", "hi"
      )]

      all_results[[key]] <- if (isTRUE(keep_raw)) res else NULL
      all_tables[[key]] <- tab
      all_diff_tables[[key]] <- diff_tab
      all_mc_tables[[key]] <- mc_tab

      write_csv_RLUP(
        tab,
        file.path(output_dir, paste0(RLUP_safe_name(key), "_prior_summary.csv"))
      )
      write_csv_RLUP(
        diff_tab,
        file.path(output_dir, paste0(RLUP_safe_name(key), "_paired_differences.csv"))
      )
      write_csv_RLUP(
        mc_tab,
        file.path(output_dir, paste0(RLUP_safe_name(key), "_MC_intervals.csv"))
      )
    }
  }

  method_table <- do.call(rbind, all_tables)
  rownames(method_table) <- NULL

  difference_table <- do.call(rbind, all_diff_tables)
  rownames(difference_table) <- NULL

  mc_table <- do.call(rbind, all_mc_tables)
  rownames(mc_table) <- NULL

  # Across-scenario descriptive summary.  This is not a replacement for the
  # condition-specific table; it is only a compact sensitivity summary.
  overall_rows <- lapply(names(priors), function(nm) {
    z <- method_table[method_table$prior == nm, , drop = FALSE]
    d <- difference_table[difference_table$prior == nm, , drop = FALSE]

    data.frame(
      prior = nm,
      prior_label = priors[[nm]]$label,
      mean_abs_bias = mean(abs(z$bias)),
      mean_RMSE = mean(z$RMSE),
      mean_abs_MSE_ratio_minus_1 = mean(abs(z$MSE_ratio - 1)),
      mean_coverage = mean(z$coverage),
      mean_width = mean(z$width),
      mean_posterior_h2 = mean(z$posterior_h2_mean),
      mean_abs_point_difference_vs_baseline = mean(d$mean_abs_point_difference),
      mean_abs_h2_difference_vs_baseline = mean(d$mean_abs_h2_difference),
      max_edge_mass = max(z$max_edge_mass),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })

  overall_table <- do.call(rbind, overall_rows)

  settings <- list(
    S = S,
    h2_values = h2_values,
    info_fractions = info_fractions,
    n_select = n_select,
    sigma2_P = sigma2_P,
    beta = beta,
    target = target,
    priors = priors,
    M = M,
    ncoarse = ncoarse,
    ngrid = ngrid,
    log_drop = log_drop,
    checkpoint_every = checkpoint_every,
    R_boot = R_boot,
    baseline = baseline,
    seed_info = seed_info,
    seed_outer = seed_outer,
    seed_bootstrap = seed_bootstrap,
    rng_order = "latent_list_before_seed_methods_v1"
  )

  out <- list(
    settings = settings,
    method_table = method_table,
    difference_table = difference_table,
    MC_interval_table = mc_table,
    overall_table = overall_table,
    raw = if (isTRUE(keep_raw)) all_results else NULL
  )

  write_csv_RLUP(
    method_table,
    file.path(output_dir, "RLUP_prior_sensitivity_method_table.csv")
  )
  write_csv_RLUP(
    difference_table,
    file.path(output_dir, "RLUP_prior_sensitivity_paired_differences.csv")
  )
  write_csv_RLUP(
    mc_table,
    file.path(output_dir, "RLUP_prior_sensitivity_MC_intervals.csv")
  )
  write_csv_RLUP(
    overall_table,
    file.path(output_dir, "RLUP_prior_sensitivity_overall_summary.csv")
  )

  saveRDS(
    out,
    file.path(output_dir, "RLUP_prior_sensitivity_MASTER.rds")
  )

  cat("\n============================================================\n")
  cat("RL-UP prior-sensitivity analysis completed.\n")
  cat("Output directory:\n", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("============================================================\n\n")

  print(overall_table, row.names = FALSE)

  invisible(out)
}


# ---------------------------------------------------------------------------
# 11. Tables convenient for Supplementary Information
# ---------------------------------------------------------------------------

make_RLUP_supplementary_table <- function(x, digits = 3) {

  z <- x$method_table[, c(
    "h2", "information", "n_pheno", "prior",
    "bias", "RMSE", "MSE_ratio", "coverage", "coverage_MCSE", "width",
    "posterior_h2_mean", "h2_bias", "h2_interval_coverage",
    "mean_edge_mass"
  )]

  num <- vapply(z, is.numeric, logical(1))
  z[num] <- lapply(z[num], round, digits = digits)
  z
}


make_RLUP_supplementary_difference_table <- function(x, digits = 3) {

  z <- x$difference_table[, c(
    "h2", "information", "n_pheno", "prior", "baseline",
    "mean_point_difference", "mean_abs_point_difference",
    "coverage_difference", "MSE_ratio_difference",
    "mean_h2_difference", "mean_abs_h2_difference"
  )]

  num <- vapply(z, is.numeric, logical(1))
  z[num] <- lapply(z[num], round, digits = digits)
  z
}


make_RLUP_weak_information_table <- function(x, digits = 3) {
  z <- make_RLUP_supplementary_table(x, digits = digits)
  z[z$h2 == 0.05, , drop = FALSE]
}


# ---------------------------------------------------------------------------
# 12. Diagnostics / interpretation helpers
# ---------------------------------------------------------------------------

check_RLUP_edge_mass <- function(x, threshold = 0.01) {

  z <- x$method_table[, c(
    "h2", "information", "n_pheno", "prior",
    "mean_edge_mass", "max_edge_mass"
  )]

  bad <- z[z$max_edge_mass > threshold, , drop = FALSE]

  if (nrow(bad) == 0L) {
    cat("No scenario/prior exceeded edge-mass threshold =", threshold, "\n")
  } else {
    cat("WARNING: posterior grid edge mass exceeded threshold.\n")
    print(bad, row.names = FALSE)
    cat("Consider increasing log_drop, padding, or grid resolution for these cases.\n")
  }

  invisible(bad)
}


# ---------------------------------------------------------------------------
# 13. Small self-test
# ---------------------------------------------------------------------------

self_test_RLUP_prior_sensitivity <- function() {

  check_RLUP_public_core()

  tmp <- file.path(tempdir(), "RLUP_prior_sensitivity_selftest")

  z <- run_RLUP_prior_sensitivity_final(
    S = 2,
    h2_values = c(0.05, 0.40),
    info_fractions = c(low = 0.50, high = 1.00),
    priors = default_RLUP_priors(),
    M = 100,
    ncoarse = 15,
    ngrid = 21,
    checkpoint_every = 1,
    R_boot = 100,
    output_dir = tmp,
    keep_raw = TRUE,
    resume = FALSE
  )

  print(make_RLUP_supplementary_table(z))
  print(make_RLUP_supplementary_difference_table(z))
  check_RLUP_edge_mass(z, threshold = 0.10)

  cat("\nRL-UP prior-sensitivity self-test completed.\n")
  invisible(z)
}


# ---------------------------------------------------------------------------
# 14. Production-result validation
# ---------------------------------------------------------------------------

RLUP_expected_production_rows <- function() {
  data.frame(
    h2 = c(
      rep(0.05, 9),
      rep(0.20, 3),
      rep(0.40, 3)
    ),
    n_pheno = c(
      rep(c(60, 90, 120), 3),
      60, 90, 120,
      60, 90, 120
    ),
    prior = c(
      rep("beta11", 3),
      rep("beta22", 3),
      rep("beta05", 3),
      rep("beta11", 3),
      rep("beta11", 3)
    ),
    coverage_expected = c(
      0.905, 0.908, 0.930,
      0.829, 0.855, 0.885,
      0.922, 0.927, 0.936,
      0.957, 0.947, 0.948,
      0.963, 0.946, 0.940
    ),
    MSE_ratio_expected = c(
      0.691, 0.809, 0.898,
      0.483, 0.568, 0.644,
      0.928, 1.089, 1.139,
      1.052, 1.040, 1.008,
      1.141, 1.030, 0.949
    ),
    stringsAsFactors = FALSE
  )
}

validate_RLUP_prior_sensitivity_production <- function(
    x,
    rounded_tolerance = 0.0015,
    edge_mass_limit = 0.10,
    stop_on_failure = TRUE) {

  tab <- if (is.character(x) && length(x) == 1L) {
    read.csv(x, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.list(x) && !is.null(x$method_table)) {
    x$method_table
  } else if (is.data.frame(x)) {
    x
  } else {
    stop("x must be an RL-UP result object, method table, or method-table CSV path.")
  }

  req <- c("h2", "n_pheno", "prior", "coverage", "MSE_ratio", "max_edge_mass")
  miss <- setdiff(req, names(tab))
  if (length(miss) > 0L) stop("Missing RL-UP output columns: ", paste(miss, collapse = ", "))

  exp <- RLUP_expected_production_rows()
  obs <- merge(
    exp,
    tab[, req],
    by = c("h2", "n_pheno", "prior"),
    all.x = TRUE,
    sort = FALSE
  )

  obs$coverage_abs_diff <- abs(obs$coverage - obs$coverage_expected)
  obs$MSE_ratio_abs_diff <- abs(obs$MSE_ratio - obs$MSE_ratio_expected)
  obs$coverage_match <- is.finite(obs$coverage_abs_diff) & obs$coverage_abs_diff <= rounded_tolerance
  obs$MSE_ratio_match <- is.finite(obs$MSE_ratio_abs_diff) & obs$MSE_ratio_abs_diff <= rounded_tolerance
  obs$edge_mass_ok <- is.finite(obs$max_edge_mass) & obs$max_edge_mass <= edge_mass_limit
  global_max_edge_mass <- max(tab$max_edge_mass, na.rm = TRUE)
  global_edge_mass_ok <- is.finite(global_max_edge_mass) && global_max_edge_mass <= edge_mass_limit
  obs$pass <- obs$coverage_match & obs$MSE_ratio_match & obs$edge_mass_ok

  if (isTRUE(stop_on_failure) && (!all(obs$pass) || !global_edge_mass_ok)) {
    print(obs[!obs$pass, , drop = FALSE], row.names = FALSE)
    stop("RL-UP production validation failed. Do not overwrite manuscript numbers; investigate RNG/output provenance.")
  }

  overall_pass <- all(obs$pass) && global_edge_mass_ok
  cat("RL-UP RNG-fixed production validation: ", if (overall_pass) "PASS" else "CHECK REQUIRED", "\n", sep = "")
  cat("Observed maximum edge mass across all scenario/prior rows: ", format(global_max_edge_mass, digits = 6), "\n", sep = "")
  invisible(obs)
}
