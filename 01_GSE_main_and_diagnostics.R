# ============================================================================
# GSE PUBLIC RELEASE -- MAIN ANALYSIS AND DIAGNOSTICS
#
# Manuscript estimand: G = c^T u, where c = c(Y) is the realised REML-EBLUP
# selection-policy contrast formed in each OUTER Monte Carlo replicate.
#
# IMPORTANT historical-name note:
# Some validated internal objects/comments use the legacy label "DeltaG".
# In this release those objects refer to G = c^T u itself; no factor 1/2 is
# applied.  The biological next-generation response G/2 is not the inferential
# target of these simulations.
#
# Scientific-computation policy for this release:
#   * the validated production logic is preserved;
#   * internal names Bayesian_policy / Bayesian_Beta11 are retained because
#     they are the historically validated RL-UP implementation;
#   * only one confirmed overwritten legacy function definition was removed;
#   * no package dependency is introduced (base R only).
# ============================================================================

# ============================================================================
# PUBLIC REPRODUCIBILITY CODE
# Prediction error of genetic gain after EBLUP selection
# Effects of REML variance-component estimation and same-data dependence
#
# Version: 2026-08-26
# Target journal: Genetics Selection Evolution
# Language: R (base R only)
#
# PURPOSE
# -------
# This single standalone file contains the validated simulation and diagnostic
# code used for the manuscript and its Supplementary Information:
#
#   1. Simulation I: fixed-policy calibration
#   2. Simulation II: realised EBLUP-selection policy
#   3. Conditional full-PEC analysis
#   4. RAM-N uncertainty propagation
#   5. RL-UP uncertainty propagation
#   6. Parametric-bootstrap diagnostics
#   7. True-variance-component reference analyses
#   8. Outer-vs-inner second-moment decomposition
#   9. Inner-reselection mechanism diagnostic
#  10. Same-sample vs cross-sample policy diagnostic
#  11. Two-layer common-data / policy-switch decomposition
#  12. Selection-intensity sensitivity for k = 2, 5, 10
#
# IMPORTANT TERMINOLOGY NOTE
# --------------------------
# The validated historical code internally uses the labels
#   "Bayesian_Beta11" and "Bayesian_policy"
# for the grid-based calculation called RL-UP in the manuscript.
# These internal names are intentionally preserved here so that the validated
# computational logic is not altered by cosmetic renaming.
#
# In manuscript terminology, RL-UP uses
#   h2 ~ Beta(1,1)
#   psi = log(sigma_A^2 + sigma_e^2) ~ N(0,1)
# with the Jacobian for transformation to
#   eta = (log sigma_A^2, log sigma_e^2),
# and propagates the resulting restricted-likelihood x proper-prior weights
# through EBLUP and the full PEC.  See the manuscript/Supplementary Information
# for the statistical interpretation and numerical-integration details.
#
# REPRODUCIBILITY SETTINGS USED IN THE PAPER
# ------------------------------------------
# Pedigree seed:          20260816
# Information-level seed: 20260850
# Outer-replicate seed:    20260851
# Simulation-II h2:        0.05, 0.20, 0.40
# Simulation-II n_pheno:   60, 90, 120
# Candidates:              20 latest-generation males
# Main selected set:       top 5 by EBLUP
# Outer replicates:        1000
# RAM-N draws:             2000
# Bootstrap replicates:    500
# RL-UP final grid:        41 x 41 (after adaptive coarse-grid localisation)
#
# The production calculations are computationally intensive.  Run the small
# self-tests first, then a pilot, then the full production analysis.
#
# BASIC USE
# ---------
# 1) Save this text as, for example:
#      01_GSE_main_and_diagnostics.R
#
# 2) Source it in a clean R session:
#      source("01_GSE_main_and_diagnostics.R")
#
# 3) Optional structural self-tests:
#      self_test_bootstrap_generator_diagnostic()
#      self_test_outer_vs_trueVC_inner_decomposition()
#
# 4) Main Simulation I + II production run:
#      final <- run_all_final(
#        S_I = 1000,
#        S_II = 1000,
#        M = 2000,
#        B = 500,
#        ngrid = 41,
#        checkpoint_every = 25,
#        root_output_dir = file.path(getwd(), "paper_output"),
#        resume = FALSE
#      )
#
# 5) Bootstrap generating-model diagnostic:
#      boot1000 <- run_bootstrap_generator_diagnostic_final(
#        S = 1000,
#        B = 500,
#        checkpoint_every = 25,
#        output_dir = file.path(getwd(), "paper_output", "bootstrap_diagnostic"),
#        resume = FALSE
#      )
#
# 6) Outer-vs-inner decomposition:
#      dec_inner1000 <- run_outer_vs_trueVC_inner_decomposition(
#        boot1000 = boot1000,
#        checkpoint_every = 10,
#        resume = FALSE
#      )
#
# 7) Inner-reselection diagnostic:
#      rsel1000 <- run_inner_reselection_VC_diagnostic(
#        boot1000 = boot1000,
#        dec_inner1000 = dec_inner1000,
#        checkpoint_every = 10,
#        resume = FALSE
#      )
#
# 8) Same-sample vs cross-sample diagnostic:
#      cross1000 <- run_same_vs_cross_selection_VC_diagnostic(
#        boot1000 = boot1000,
#        rsel1000 = rsel1000,
#        n_cross_shifts = 5,
#        checkpoint_every = 10,
#        resume = FALSE
#      )
#
# 9) Two-layer decomposition:
#      layer1000 <- run_two_layer_alignment_decomposition(
#        boot1000 = boot1000,
#        cross_reference = cross1000,
#        checkpoint_every = 10,
#        resume = FALSE
#      )
#
# 10) Selection-intensity sensitivity:
#      ksens1000 <- run_selection_intensity_sensitivity(
#        boot1000 = boot1000,
#        layer_reference = layer1000,
#        k_values = c(2, 5, 10),
#        checkpoint_every = 10,
#        resume = FALSE,
#        MC_R = 5000
#      )
#
# OUTPUTS
# -------
# Each production runner writes CSV summaries and RDS objects to a relative
# output directory.  No absolute user-specific paths are required.
#
# PUBLIC-RELEASE CHECKLIST
# ------------------------
# Before depositing in a repository, add:
#   - final manuscript citation / DOI when available;
#   - repository DOI;
#   - sessionInfo() from the environment used for the archived run;
#   - a license file (for example, MIT) if desired.
#
# This code intentionally contains no private data and simulates all phenotypes,
# breeding values, and pedigrees required for the study.
# ============================================================================

# ======================================================================
# NEXT PRIORITY DIAGNOSTIC AFTER TRUE-VC ORACLE
# V3: safe resume; final RDS must match S, B, and scenario key.
# Pilot and production outputs are separated automatically by S and B.
# Simulation II: plug-in vs TRUE-VC bootstrap generating model
# ======================================================================
#
# Main quick check:
#   boot20 <- run_bootstrap_generator_diagnostic_final(
#     S = 20, B = 100, checkpoint_every = 5, resume = FALSE
#   )
#   print(boot20$comparison_table)
#   print(boot20$diagnostic_table)
#
# Main production run:
#   boot1000 <- run_bootstrap_generator_diagnostic_final(
#     S = 1000, B = 500, checkpoint_every = 25, resume = TRUE
#   )
#
# Optional exact reproduction check against previous Simulation II:
#   check <- verify_bootstrap_generator_against_previous(boot1000)
#   print(check)
#
# The file contains the complete earlier production Simulation I/II code
# plus the new diagnostic extension. No other R file needs to be sourced.
# ======================================================================


# ======================================================================
# FINAL PRODUCTION CODE
# Genetic gain uncertainty: Simulation I + Simulation II
# ======================================================================
#
# ONE FILE / ONE SOURCE
#
# After saving this file, only:
#
#
# is needed.
#
# Then:
#
#   final <- run_all_final()
#
# ----------------------------------------------------------------------
# FINAL STUDY STRUCTURE
# ----------------------------------------------------------------------
#
# Simulation I
#   Balanced sire model, FIXED policy.
#   Model:
#     y_ij = mu + 0.5*u_i + e_ij
#     u_i ~ N(0, sigma2_A)
#     e_ij ~ N(0, sigma2_e)
#
#   Scenario A: m=24, k=8, sigma2_A=30, sigma2_e=70
#   Scenario B: m=16, k=5, sigma2_A=10, sigma2_e=90
#
#   Methods:
#     conditional
#     RAM_N
#     bootstrap
#     Bayesian_Beta11
#
# Simulation II
#   Pedigree animal model, realised EBV-selection policy.
#   In each OUTER replicate:
#     REML -> EBLUP -> top-5 selection -> construct c(y)
#     -> hold selected set fixed inside RAM-N/bootstrap/Bayesian
#     -> truth = c(y)'u
#
#   Factorial design:
#     h2 = 0.05, 0.20, 0.40
#       x
#     phenotypic information ~ 60, 90, 120 records
#
#   Methods:
#     conditional
#     conditional_diag
#     RAM_N_policy
#     bootstrap_policy
#     Bayesian_policy
#
#   Main Bayesian prior:
#     h2 ~ Beta(1,1)
#     log(sigma_P^2) ~ N(0, 1^2)
#
# RAM-RL is NOT used by the FINAL runners.
#
# Default production settings:
#   outer replicates = 1000
#   RAM/Bayesian draws = 2000
#   bootstrap replicates = 500
#   Bayesian grid = 41 x 41
#   checkpoint every 25 outer replicates
#
# ======================================================================


# --- Embedded validated animal-model core ---

# ============================================================
# Simulation II: Pedigree-based animal model
# Genetic gain uncertainty with full PEC, RAM-N, RAM-RL,
# parametric bootstrap, and a grid-based Bayesian analysis.
#
# Base R only.
# Main model:
#   y = X beta + Z u + e
#   u ~ N(0, A * sigma2_A)
#   e ~ N(0, I * sigma2_e)
#
# True policy gain:
#   G_true = c' u
#
# Notes:
# - The baseline pedigree uses relatively broad, uniform parent use.
# - The pedigree and A matrix are generated once and then held fixed
#   across repeated Monte Carlo replicates.
# - Early fixed-policy validation helpers use a fixed c to isolate uncertainty.
# - The final Simulation-II production runner reconstructs the realised c(Y)
#   from REML-EBLUP in every outer replicate.
# - RAM-RL is defined here on eta = (log sigma2_A, log sigma2_e):
#     q_RL(eta | y) proportional to exp{ l_R(eta) } d eta
# - The Bayesian analysis uses the same eta grid but multiplies the
#   restricted likelihood by an explicit proper prior on eta.
# ============================================================

# ---------- Basic helpers ----------

rmvn1 <- function(mu, Sigma) {
  R <- chol((Sigma + t(Sigma)) / 2)
  drop(mu + t(R) %*% rnorm(length(mu)))
}

rmvn_rows <- function(n, mu, Sigma) {
  R <- chol((Sigma + t(Sigma)) / 2)
  Z <- matrix(rnorm(n * length(mu)), nrow = n, ncol = length(mu))
  sweep(Z %*% R, 2, mu, "+")
}

safe_quantile <- function(x, probs) {
  unname(quantile(x, probs = probs, names = FALSE, type = 8))
}

# ---------- Pedigree generation ----------

make_moderate_pedigree <- function(
    n_founders = 24,
    n_gen = 3,
    n_per_gen = 40,
    sire_fraction = 0.75,
    dam_fraction = 0.95,
    parent_use = c("uniform", "weighted"),
    seed = 20260816) {

  parent_use <- match.arg(parent_use)

  stopifnot(n_founders >= 8, n_per_gen >= 12)
  set.seed(seed)

  ped <- data.frame(
    id = seq_len(n_founders),
    sire = 0L,
    dam = 0L,
    sex = rep(c("M", "F"), length.out = n_founders),
    generation = 0L,
    stringsAsFactors = FALSE
  )

  next_id <- n_founders + 1L

  for (g in seq_len(n_gen)) {
    prev <- ped[ped$generation == (g - 1L), , drop = FALSE]
    males <- prev$id[prev$sex == "M"]
    females <- prev$id[prev$sex == "F"]

    nsire <- max(3L, min(length(males), round(length(males) * sire_fraction)))
    ndam  <- max(6L, min(length(females), round(length(females) * dam_fraction)))

    used_sires <- sample(males, nsire, replace = FALSE)
    used_dams  <- sample(females, ndam, replace = FALSE)

    # Baseline: uniform parent use.  The weighted option intentionally
    # concentrates contributions and is useful as a high-relatedness
    # sensitivity scenario.
    if (parent_use == "uniform") {
      ps <- rep(1 / nsire, nsire)
      pd <- rep(1 / ndam, ndam)
    } else {
      ps <- seq(nsire, 1) + 1
      ps <- ps / sum(ps)
      pd <- seq(ndam, 1) + 2
      pd <- pd / sum(pd)
    }

    sire_vec <- sample(used_sires, n_per_gen, replace = TRUE, prob = ps)
    dam_vec  <- sample(used_dams,  n_per_gen, replace = TRUE, prob = pd)

    sex_vec <- sample(rep(c("M", "F"), length.out = n_per_gen))

    new <- data.frame(
      id = next_id:(next_id + n_per_gen - 1L),
      sire = as.integer(sire_vec),
      dam = as.integer(dam_vec),
      sex = sex_vec,
      generation = as.integer(g),
      stringsAsFactors = FALSE
    )

    ped <- rbind(ped, new)
    next_id <- next_id + n_per_gen
  }

  rownames(ped) <- NULL
  ped
}

# ---------- Numerator relationship matrix A ----------

make_A <- function(ped) {
  n <- nrow(ped)
  if (!all(ped$id == seq_len(n))) {
    stop("ped$id must be 1:n and parents must precede offspring.")
  }

  A <- matrix(0, n, n)

  for (i in seq_len(n)) {
    s <- ped$sire[i]
    d <- ped$dam[i]

    if (i > 1L) {
      for (j in seq_len(i - 1L)) {
        asj <- if (s == 0L) 0 else A[s, j]
        adj <- if (d == 0L) 0 else A[d, j]
        A[i, j] <- A[j, i] <- 0.5 * (asj + adj)
      }
    }

    if (s == 0L || d == 0L) {
      A[i, i] <- 1
    } else {
      A[i, i] <- 1 + 0.5 * A[s, d]
    }
  }

  A
}

relationship_summary <- function(A, ids = seq_len(nrow(A))) {
  S <- A[ids, ids, drop = FALSE]
  off <- S[upper.tri(S)]
  F <- diag(S) - 1

  c(
    n = length(ids),
    mean_offdiag = mean(off),
    median_offdiag = median(off),
    q90_offdiag = unname(quantile(off, 0.90)),
    max_offdiag = max(off),
    mean_inbreeding = mean(F),
    max_inbreeding = max(F)
  )
}

# ---------- Design / policy ----------

make_policy <- function(n_animals, candidates, selected) {
  if (!all(selected %in% candidates)) stop("selected must be a subset of candidates.")
  cstar <- numeric(n_animals)
  cstar[selected] <- cstar[selected] + 1 / length(selected)
  cstar[candidates] <- cstar[candidates] - 1 / length(candidates)
  stopifnot(abs(sum(cstar)) < 1e-12)
  cstar
}

build_simII_design <- function(
    n_founders = 24,
    n_gen = 3,
    n_per_gen = 40,
    sire_fraction = 0.75,
    dam_fraction = 0.95,
    parent_use = "uniform",
    n_select = 5,
    seed_pedigree = 20260816) {

  ped <- make_moderate_pedigree(
    n_founders = n_founders,
    n_gen = n_gen,
    n_per_gen = n_per_gen,
    sire_fraction = sire_fraction,
    dam_fraction = dam_fraction,
    parent_use = parent_use,
    seed = seed_pedigree
  )

  A <- make_A(ped)

  # Phenotypes: all non-founders.
  pheno_ids <- ped$id[ped$generation >= 1L]

  # Candidate set: males in the latest generation.
  gmax <- max(ped$generation)
  candidates <- ped$id[ped$generation == gmax & ped$sex == "M"]

  if (length(candidates) < n_select + 2L) {
    stop("Too few candidate males; increase n_per_gen or reduce n_select.")
  }

  # Default fixed-policy placeholder used by legacy/fixed-policy validation helpers.
  # The final realised-policy Simulation-II runner ignores this placeholder for
  # selection and reconstructs top-k c(Y) from REML-EBLUP in each outer replicate.
  selected <- candidates[seq_len(n_select)]
  cstar <- make_policy(nrow(ped), candidates, selected)

  X <- matrix(1, nrow = length(pheno_ids), ncol = 1)
  colnames(X) <- "intercept"

  prep <- prepare_animal_model(A, pheno_ids, X)

  list(
    ped = ped,
    A = A,
    pheno_ids = pheno_ids,
    candidates = candidates,
    selected = selected,
    cstar = cstar,
    X = X,
    prep = prep
  )
}

# ---------- Spectral preparation for fast REML / V^{-1} ----------

prepare_animal_model <- function(A, pheno_ids, X) {
  K <- A[pheno_ids, pheno_ids, drop = FALSE]

  eigK <- eigen((K + t(K)) / 2, symmetric = TRUE)
  eigK$values <- pmax(eigK$values, 0)

  qrX <- qr(X)
  p <- qrX$rank
  nobs <- nrow(X)
  if (p >= nobs) stop("X must have rank < number of observations.")

  Qfull <- qr.Q(qrX, complete = TRUE)
  Q <- Qfull[, (p + 1L):nobs, drop = FALSE]

  KR <- crossprod(Q, K %*% Q)
  eigR <- eigen((KR + t(KR)) / 2, symmetric = TRUE)
  eigR$values <- pmax(eigR$values, 0)

  list(
    A = A,
    pheno_ids = pheno_ids,
    X = X,
    K = K,
    eigK = eigK,
    Q = Q,
    eigR = eigR,
    nobs = nobs,
    p = p,
    cAc_cache = new.env(parent = emptyenv())
  )
}

Vinv_apply <- function(prep, sigma2_A, sigma2_e, B) {
  was_vector <- is.null(dim(B))
  if (was_vector) B <- matrix(B, ncol = 1)

  U <- prep$eigK$vectors
  d <- sigma2_A * prep$eigK$values + sigma2_e
  tmp <- crossprod(U, B)
  tmp <- sweep(tmp, 1, d, "/")
  out <- U %*% tmp

  if (was_vector) drop(out) else out
}

# ---------- Restricted likelihood ----------

log_reml_eta <- function(eta, y, prep) {
  sigma2_A <- exp(eta[1])
  sigma2_e <- exp(eta[2])

  yr <- crossprod(prep$Q, y)
  z <- crossprod(prep$eigR$vectors, yr)

  d <- sigma2_A * prep$eigR$values + sigma2_e
  if (any(!is.finite(d)) || any(d <= 0)) return(-Inf)

  # Exact REML log-likelihood up to an additive constant independent of eta.
  -0.5 * sum(log(d) + (z^2) / d)
}

fit_reml_animal <- function(
    y, prep,
    lower = c(log(1e-6), log(1e-6)),
    upper = c(log(10), log(10)),
    start = NULL) {

  if (is.null(start)) {
    vp <- var(y)
    if (!is.finite(vp) || vp <= 0) vp <- 1
    start <- log(c(max(0.20 * vp, 1e-4), max(0.80 * vp, 1e-4)))
  }

  obj <- function(eta) -log_reml_eta(eta, y, prep)

  opt <- optim(
    par = start,
    fn = obj,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper
  )

  H <- optimHess(opt$par, obj)
  Hs <- (H + t(H)) / 2
  ee <- eigen(Hs, symmetric = TRUE)

  Veta <- NULL
  if (all(is.finite(ee$values)) && min(ee$values) > 1e-10) {
    Veta <- ee$vectors %*% diag(1 / ee$values) %*% t(ee$vectors)
  }

  list(
    eta = opt$par,
    theta = exp(opt$par),
    convergence = opt$convergence,
    value = opt$value,
    H = Hs,
    hess_eigen = ee$values,
    Veta = Veta,
    lower = lower,
    upper = upper
  )
}

diagnose_RAM_N <- function(
    fit,
    se_max = 3,
    hess_min = 1e-6,
    boundary_margin = 0.05) {

  reasons <- character()

  if (fit$convergence != 0) reasons <- c(reasons, "REML_nonconvergence")
  if (is.null(fit$Veta)) reasons <- c(reasons, "nonpositive_Hessian")

  if (all(is.finite(fit$hess_eigen)) &&
      min(fit$hess_eigen) <= hess_min) {
    reasons <- c(reasons, "small_Hessian_eigenvalue")
  }

  if (!is.null(fit$Veta)) {
    se <- sqrt(diag(fit$Veta))
    if (any(!is.finite(se)) || max(se) > se_max) {
      reasons <- c(reasons, "large_eta_SE")
    }
  }

  if (any(fit$eta - fit$lower < boundary_margin) ||
      any(fit$upper - fit$eta < boundary_margin)) {
    reasons <- c(reasons, "near_numeric_boundary")
  }

  list(valid = length(reasons) == 0L, reasons = unique(reasons))
}

# ---------- BLUP and policy PEV/PEC ----------

policy_stats_animal <- function(
    eta, y, prep, cstar,
    calc_diag = FALSE) {

  sigma2_A <- exp(eta[1])
  sigma2_e <- exp(eta[2])

  X <- prep$X
  A <- prep$A
  ids <- prep$pheno_ids

  ViX <- Vinv_apply(prep, sigma2_A, sigma2_e, X)
  Viy <- Vinv_apply(prep, sigma2_A, sigma2_e, y)

  XtViX <- crossprod(X, ViX)
  beta_hat <- solve(XtViX, crossprod(X, Viy))

  resid <- drop(y - X %*% beta_hat)
  Vir <- Vinv_apply(prep, sigma2_A, sigma2_e, resid)

  # Cov(y, c'u) / sigma2_A before multiplying by sigma2_A.
  a_obs_c <- drop(A[ids, , drop = FALSE] %*% cstar)

  mu <- sigma2_A * sum(a_obs_c * Vir)

  # Full prediction-error variance:
  # Var(c'u - c'uhat | theta) = c'Gc - h'Ph.
  cAc <- drop(crossprod(cstar, A %*% cstar))
  h <- sigma2_A * a_obs_c

  Vih <- Vinv_apply(prep, sigma2_A, sigma2_e, h)
  Ph <- Vih - ViX %*% solve(XtViX, crossprod(X, Vih))

  v_full <- sigma2_A * cAc - sum(h * Ph)
  v_full <- max(v_full, 0)

  out <- list(
    mu = as.numeric(mu),
    v = as.numeric(v_full),
    sd = sqrt(v_full),
    beta = drop(beta_hat)
  )

  if (calc_diag) {
    idx <- which(abs(cstar) > 0)
    Hobs <- sigma2_A * A[ids, idx, drop = FALSE]
    ViH <- Vinv_apply(prep, sigma2_A, sigma2_e, Hobs)
    PH <- ViH - ViX %*% solve(XtViX, crossprod(X, ViH))

    pev_i <- sigma2_A * diag(A)[idx] - colSums(Hobs * PH)
    pev_i <- pmax(pev_i, 0)

    v_diag <- sum((cstar[idx]^2) * pev_i)

    out$v_diag <- as.numeric(v_diag)
    out$sd_diag <- sqrt(v_diag)
    out$pev_i <- pev_i
  }

  out
}

# ---------- Data simulation ----------

simulate_animal_data <- function(
    design,
    h2 = 0.20,
    sigma2_P = 1,
    beta = 0,
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  if (h2 <= 0 || h2 >= 1) stop("h2 must be between 0 and 1.")

  sigma2_A <- h2 * sigma2_P
  sigma2_e <- (1 - h2) * sigma2_P

  A <- design$A
  n <- nrow(A)

  u <- sqrt(sigma2_A) * rmvn1(rep(0, n), A)
  e <- rnorm(length(design$pheno_ids), 0, sqrt(sigma2_e))

  y <- beta + u[design$pheno_ids] + e

  list(
    y = y,
    u = u,
    truth = sum(design$cstar * u),
    sigma2_A = sigma2_A,
    sigma2_e = sigma2_e
  )
}

# ---------- Conditional EBLUP ----------

conditional_method <- function(y, prep, cstar, target = 0) {
  fit <- fit_reml_animal(y, prep)
  g <- policy_stats_animal(fit$eta, y, prep, cstar, calc_diag = TRUE)

  lohi <- g$mu + c(-1, 1) * 1.96 * g$sd
  p <- 1 - pnorm(target, mean = g$mu, sd = max(g$sd, 1e-12))

  list(
    valid = TRUE,
    mean = g$mu,
    sd = g$sd,
    interval = lohi,
    p = p,
    fit = fit,
    v_diag = g$v_diag,
    sd_diag = g$sd_diag
  )
}

# ---------- RAM-N ----------

RAM_N_animal <- function(
    y, prep, cstar,
    M = 2000,
    target = 0,
    diagnosis = TRUE) {

  fit <- fit_reml_animal(y, prep)

  if (diagnosis) {
    dg <- diagnose_RAM_N(fit)
    if (!dg$valid) {
      return(list(valid = FALSE, reasons = dg$reasons, fit = fit))
    }
  }

  if (is.null(fit$Veta)) {
    return(list(valid = FALSE, reasons = "Veta_unavailable", fit = fit))
  }

  eta_draw <- rmvn_rows(M, fit$eta, fit$Veta)

  mu <- numeric(M)
  v <- numeric(M)

  for (i in seq_len(M)) {
    g <- policy_stats_animal(eta_draw[i, ], y, prep, cstar)
    mu[i] <- g$mu
    v[i] <- g$v
  }

  delta <- mu + sqrt(pmax(v, 0)) * rnorm(M)

  list(
    valid = TRUE,
    mean = mean(mu),
    sd = sqrt(mean(v) + var(mu)),
    interval = safe_quantile(delta, c(0.025, 0.975)),
    p = mean(delta > target),
    delta = delta,
    eta = eta_draw,
    fit = fit
  )
}

# ---------- Adaptive RAM-RL grid ----------

build_RL_grid <- function(
    y, prep, fit = NULL,
    ncoarse = 31,
    ngrid = 51,
    log_drop = 12,
    measure = c("eta", "theta")) {

  measure <- match.arg(measure)
  if (is.null(fit)) fit <- fit_reml_animal(y, prep)

  g1c <- seq(fit$lower[1], fit$upper[1], length.out = ncoarse)
  g2c <- seq(fit$lower[2], fit$upper[2], length.out = ncoarse)
  coarse <- expand.grid(eta_A = g1c, eta_e = g2c)

  coarse$logL <- vapply(
    seq_len(nrow(coarse)),
    function(i) log_reml_eta(as.numeric(coarse[i, 1:2]), y, prep),
    numeric(1)
  )

  if (measure == "theta") {
    coarse$logW <- coarse$logL + coarse$eta_A + coarse$eta_e
  } else {
    coarse$logW <- coarse$logL
  }

  mx <- max(coarse$logW)
  keep <- coarse$logW >= (mx - log_drop)

  d1 <- if (length(g1c) > 1) diff(g1c)[1] else 0
  d2 <- if (length(g2c) > 1) diff(g2c)[1] else 0

  lo1 <- max(fit$lower[1], min(coarse$eta_A[keep]) - 2 * d1)
  hi1 <- min(fit$upper[1], max(coarse$eta_A[keep]) + 2 * d1)
  lo2 <- max(fit$lower[2], min(coarse$eta_e[keep]) - 2 * d2)
  hi2 <- min(fit$upper[2], max(coarse$eta_e[keep]) + 2 * d2)

  g1 <- seq(lo1, hi1, length.out = ngrid)
  g2 <- seq(lo2, hi2, length.out = ngrid)
  grid <- expand.grid(eta_A = g1, eta_e = g2)

  grid$logL <- vapply(
    seq_len(nrow(grid)),
    function(i) log_reml_eta(as.numeric(grid[i, 1:2]), y, prep),
    numeric(1)
  )

  if (measure == "theta") {
    grid$logW <- grid$logL + grid$eta_A + grid$eta_e
  } else {
    grid$logW <- grid$logL
  }

  lw <- grid$logW - max(grid$logW)
  w <- exp(lw)
  w <- w / sum(w)

  list(
    grid = grid,
    weight = w,
    fit = fit,
    measure = measure
  )
}

RAM_RL_animal <- function(
    y, prep, cstar,
    M = 2000,
    ngrid = 51,
    target = 0,
    measure = "eta") {

  fit <- fit_reml_animal(y, prep)
  gr <- build_RL_grid(
    y, prep, fit = fit,
    ngrid = ngrid,
    measure = measure
  )

  id <- sample(
    seq_len(nrow(gr$grid)),
    size = M,
    replace = TRUE,
    prob = gr$weight
  )

  eta_draw <- as.matrix(gr$grid[id, c("eta_A", "eta_e"), drop = FALSE])

  mu <- numeric(M)
  v <- numeric(M)

  # Evaluate repeated grid points only once.
  uid <- unique(id)
  cache_mu <- numeric(length(uid))
  cache_v <- numeric(length(uid))
  names(cache_mu) <- names(cache_v) <- as.character(uid)

  for (j in seq_along(uid)) {
    ii <- uid[j]
    eta <- as.numeric(gr$grid[ii, c("eta_A", "eta_e")])
    gs <- policy_stats_animal(eta, y, prep, cstar)
    cache_mu[j] <- gs$mu
    cache_v[j] <- gs$v
  }

  mu <- cache_mu[as.character(id)]
  v <- cache_v[as.character(id)]

  delta <- mu + sqrt(pmax(v, 0)) * rnorm(M)

  list(
    valid = TRUE,
    mean = mean(mu),
    sd = sqrt(mean(v) + var(mu)),
    interval = safe_quantile(delta, c(0.025, 0.975)),
    p = mean(delta > target),
    delta = delta,
    eta = eta_draw,
    fit = fit,
    grid_obj = gr
  )
}

# ---------- Parametric bootstrap ----------
# Fixed policy; truth is simulated jointly with y under the fitted model.

parametric_bootstrap_animal <- function(
    y, prep, cstar,
    B = 300,
    target = 0) {

  fit <- fit_reml_animal(y, prep)
  sigma2_A <- fit$theta[1]
  sigma2_e <- fit$theta[2]

  gs0 <- policy_stats_animal(fit$eta, y, prep, cstar)
  current <- gs0$mu

  X <- prep$X
  ViX <- Vinv_apply(prep, sigma2_A, sigma2_e, X)
  Viy <- Vinv_apply(prep, sigma2_A, sigma2_e, y)
  beta_hat <- solve(crossprod(X, ViX), crossprod(X, Viy))

  A <- prep$A
  ids <- prep$pheno_ids
  cAc <- drop(crossprod(cstar, A %*% cstar))
  a_obs_c <- drop(A[ids, , drop = FALSE] %*% cstar)

  V <- sigma2_A * prep$K + sigma2_e * diag(prep$nobs)
  h <- sigma2_A * a_obs_c
  var_delta <- sigma2_A * cAc

  Sjoint <- rbind(
    cbind(V, h),
    c(h, var_delta)
  )
  Sjoint <- (Sjoint + t(Sjoint)) / 2
  diag(Sjoint) <- diag(Sjoint) + 1e-12

  mean_joint <- c(drop(X %*% beta_hat), 0)
  R <- chol(Sjoint)

  Zrand <- matrix(
    rnorm(B * length(mean_joint)),
    nrow = B,
    ncol = length(mean_joint)
  )
  draws <- sweep(Zrand %*% R, 2, mean_joint, "+")

  err <- numeric(B)

  for (b in seq_len(B)) {
    yb <- draws[b, seq_len(prep$nobs)]
    truth_b <- draws[b, prep$nobs + 1L]

    fit_b <- fit_reml_animal(yb, prep)
    dhat_b <- policy_stats_animal(fit_b$eta, yb, prep, cstar)$mu
    err[b] <- dhat_b - truth_b
  }

  q <- safe_quantile(err, c(0.025, 0.975))
  interval <- c(current - q[2], current - q[1])

  pseudo_delta <- current - err

  list(
    valid = TRUE,
    mean = current - mean(err),
    sd = sd(err),
    interval = interval,
    p = mean(pseudo_delta > target),
    error = err,
    fit = fit
  )
}

# ---------- Bayesian grid analysis ----------
# Proper prior on eta = log variance components.
# Default is intentionally broad on the sigma_P^2 = 1 simulation scale.

bayesian_grid_animal <- function(
    y, prep, cstar,
    M = 2000,
    ngrid = 51,
    target = 0,
    prior_mean = log(c(0.20, 0.80)),
    prior_sd = c(2, 2),
    grid_obj = NULL) {

  if (is.null(grid_obj)) {
    fit <- fit_reml_animal(y, prep)
    grid_obj <- build_RL_grid(
      y, prep, fit = fit,
      ngrid = ngrid,
      measure = "eta"
    )
  }

  grid <- grid_obj$grid

  logprior <-
    dnorm(grid$eta_A, mean = prior_mean[1], sd = prior_sd[1], log = TRUE) +
    dnorm(grid$eta_e, mean = prior_mean[2], sd = prior_sd[2], log = TRUE)

  logpost <- grid$logL + logprior
  w <- exp(logpost - max(logpost))
  w <- w / sum(w)

  id <- sample(
    seq_len(nrow(grid)),
    size = M,
    replace = TRUE,
    prob = w
  )

  uid <- unique(id)
  cache_mu <- numeric(length(uid))
  cache_v <- numeric(length(uid))
  names(cache_mu) <- names(cache_v) <- as.character(uid)

  for (j in seq_along(uid)) {
    ii <- uid[j]
    eta <- as.numeric(grid[ii, c("eta_A", "eta_e")])
    gs <- policy_stats_animal(eta, y, prep, cstar)
    cache_mu[j] <- gs$mu
    cache_v[j] <- gs$v
  }

  mu <- cache_mu[as.character(id)]
  v <- cache_v[as.character(id)]

  delta <- mu + sqrt(pmax(v, 0)) * rnorm(M)

  list(
    valid = TRUE,
    mean = mean(delta),
    sd = sd(delta),
    interval = safe_quantile(delta, c(0.025, 0.975)),
    p = mean(delta > target),
    delta = delta,
    eta = as.matrix(grid[id, c("eta_A", "eta_e"), drop = FALSE]),
    posterior_weight = w,
    grid_obj = grid_obj
  )
}

# ---------- One full replicate ----------

one_simII <- function(
    design,
    h2 = 0.20,
    sigma2_P = 1,
    beta = 0,
    target = 0,
    M = 1000,
    B = 150,
    ngrid = 41,
    seed = NULL) {

  dat <- simulate_animal_data(
    design,
    h2 = h2,
    sigma2_P = sigma2_P,
    beta = beta,
    seed = seed
  )

  cond <- conditional_method(
    dat$y, design$prep, design$cstar, target = target
  )

  rn <- RAM_N_animal(
    dat$y, design$prep, design$cstar,
    M = M, target = target
  )

  rr <- RAM_RL_animal(
    dat$y, design$prep, design$cstar,
    M = M, ngrid = ngrid, target = target
  )

  pb <- parametric_bootstrap_animal(
    dat$y, design$prep, design$cstar,
    B = B, target = target
  )

  by <- bayesian_grid_animal(
    dat$y, design$prep, design$cstar,
    M = M, ngrid = ngrid, target = target,
    grid_obj = rr$grid_obj
  )

  pack <- function(x) {
    if (is.null(x$valid) || !x$valid) return(NULL)
    c(
      mean = x$mean,
      sd = x$sd,
      lo = unname(x$interval[1]),
      hi = unname(x$interval[2]),
      p = x$p
    )
  }

  out <- list(
    truth = dat$truth,
    true_sigma2_A = dat$sigma2_A,
    true_sigma2_e = dat$sigma2_e,
    conditional = pack(cond),
    RAM_N = pack(rn),
    RAM_RL = pack(rr),
    bootstrap = pack(pb),
    Bayesian = pack(by),
    RAM_N_reasons = if (!is.null(rn$valid) && !rn$valid) rn$reasons else character(),
    conditional_diag = c(
      mean = cond$mean,
      sd = cond$sd_diag,
      lo = cond$mean - 1.96 * cond$sd_diag,
      hi = cond$mean + 1.96 * cond$sd_diag,
      p = 1 - pnorm(target, mean = cond$mean, sd = max(cond$sd_diag, 1e-12))
    )
  )

  out
}

# ---------- Monte Carlo summaries ----------

summarize_simII <- function(res, name, target = 0) {
  x <- lapply(res, `[[`, name)
  ok <- !vapply(x, is.null, logical(1))
  x <- x[ok]

  if (length(x) == 0L) {
    return(c(
      valid = 0, bias = NA, RMSE = NA, mean_SE = NA,
      MSE_ratio = NA, inclusion = NA, width = NA,
      Brier = NA, mean_p = NA, event_rate = NA
    ))
  }

  truth <- vapply(res[ok], `[[`, numeric(1), "truth")
  mat <- do.call(rbind, x)

  err <- mat[, "mean"] - truth
  mc_mse <- mean(err^2)
  mean_est_mse <- mean(mat[, "sd"]^2)
  z <- as.numeric(truth > target)

  c(
    valid = length(x),
    bias = mean(err),
    RMSE = sqrt(mc_mse),
    mean_SE = mean(mat[, "sd"]),
    MSE_ratio = mean_est_mse / mc_mse,
    inclusion = mean(truth >= mat[, "lo"] & truth <= mat[, "hi"]),
    width = mean(mat[, "hi"] - mat[, "lo"]),
    Brier = mean((mat[, "p"] - z)^2),
    mean_p = mean(mat[, "p"]),
    event_rate = mean(z)
  )
}

run_MC_simII <- function(
    S = 50,
    design,
    h2 = 0.20,
    sigma2_P = 1,
    beta = 0,
    target = 0,
    M = 500,
    B = 100,
    ngrid = 41,
    seed = 20260817,
    keep_raw = FALSE) {

  set.seed(seed)
  replicate_seeds <- sample.int(.Machine$integer.max, S)

  res <- vector("list", S)

  for (s in seq_len(S)) {
    res[[s]] <- one_simII(
      design = design,
      h2 = h2,
      sigma2_P = sigma2_P,
      beta = beta,
      target = target,
      M = M,
      B = B,
      ngrid = ngrid,
      seed = replicate_seeds[s]
    )

    if (s %% 10 == 0) {
      cat("Completed", s, "of", S, "replicates\n")
    }
  }

  methods <- c(
    "conditional",
    "conditional_diag",
    "RAM_N",
    "RAM_RL",
    "bootstrap",
    "Bayesian"
  )

  ans <- t(vapply(
    methods,
    function(z) summarize_simII(res, z, target = target),
    numeric(10)
  ))
  rownames(ans) <- methods

  if (keep_raw) {
    list(summary = ans, raw = res)
  } else {
    ans
  }
}


# ---------- Alternative relatedness scenarios ----------

build_simII_design_high_relatedness <- function(
    n_founders = 24,
    n_gen = 3,
    n_per_gen = 40,
    n_select = 5,
    seed_pedigree = 20260816) {

  build_simII_design(
    n_founders = n_founders,
    n_gen = n_gen,
    n_per_gen = n_per_gen,
    sire_fraction = 0.50,
    dam_fraction = 0.80,
    parent_use = "weighted",
    n_select = n_select,
    seed_pedigree = seed_pedigree
  )
}

# ---------- Convenient demo ----------

run_simII_demo <- function() {
  cat("\n--- Building pedigree-based animal model design ---\n")

  design <- build_simII_design()

  cat("\nNumber of animals:", nrow(design$ped), "\n")
  cat("Phenotyped animals:", length(design$pheno_ids), "\n")
  cat("Candidate males:", length(design$candidates), "\n")
  cat("Fixed selected candidates:", paste(design$selected, collapse = ", "), "\n")

  cat("\nRelationship summary among candidate males:\n")
  print(relationship_summary(design$A, design$candidates))

  cat("\n--- One simulated dataset: h2 = 0.20 ---\n")
  set.seed(20260818)
  out <- one_simII(
    design = design,
    h2 = 0.20,
    M = 1000,
    B = 150,
    ngrid = 41,
    seed = 20260818
  )

  cat("True DeltaG =", out$truth, "\n")

  methods <- c("conditional", "conditional_diag", "RAM_N", "RAM_RL", "bootstrap", "Bayesian")
  one_table <- do.call(rbind, lapply(methods, function(m) {
    x <- out[[m]]
    if (is.null(x)) {
      c(mean = NA, sd = NA, lo = NA, hi = NA, p = NA)
    } else {
      x
    }
  }))
  rownames(one_table) <- methods

  print(one_table)

  if (is.null(out$RAM_N)) {
    cat("\nRAM-N diagnostic did not pass. Reasons:\n")
    print(out$RAM_N_reasons)
  }

  invisible(list(design = design, one = out))
}

# ---------- Self-test ----------

self_test_simII <- function() {
  design <- build_simII_design(
    n_founders = 16,
    n_gen = 2,
    n_per_gen = 20,
    n_select = 3,
    seed_pedigree = 123
  )

  stopifnot(
    nrow(design$A) == nrow(design$ped),
    max(abs(design$A - t(design$A))) < 1e-10,
    min(eigen(design$A, symmetric = TRUE, only.values = TRUE)$values) > -1e-8,
    abs(sum(design$cstar)) < 1e-12
  )

  dat <- simulate_animal_data(design, h2 = 0.20, seed = 456)
  fit <- fit_reml_animal(dat$y, design$prep)
  stopifnot(all(is.finite(fit$theta)), all(fit$theta > 0))

  gs <- policy_stats_animal(
    fit$eta, dat$y, design$prep, design$cstar,
    calc_diag = TRUE
  )
  stopifnot(
    is.finite(gs$mu),
    is.finite(gs$v),
    gs$v >= 0,
    is.finite(gs$v_diag),
    gs$v_diag >= 0
  )

  mc <- run_MC_simII(
    S = 2,
    design = design,
    h2 = 0.20,
    M = 100,
    B = 20,
    ngrid = 21,
    seed = 789
  )

  print(mc)
  cat("\nSimulation II self-test completed.\n")
  invisible(TRUE)
}

# ============================================================
# HOW TO RUN
# ============================================================
#
# 1) Quick structural test:
#    self_test_simII()
#
# 2) One full demonstration:
#    demo <- run_simII_demo()
#
# 3) Small Monte Carlo trial:
#    design <- build_simII_design()
#    mc20 <- run_MC_simII(
#      S = 20,
#      design = design,
#      h2 = 0.20,
#      M = 500,
#      B = 100,
#      ngrid = 41,
#      seed = 20260819
#    )
#    print(mc20)
#
# 4) After the code has been validated, larger S can be used.
#    Keep in mind that the bootstrap and RAM-RL grid dominate cost.
#
# ============================================================

# --- Embedded validated no-RAM-RL helpers ---

# ============================================================
# Simulation II validation without RAM-RL
# Pedigree-based animal model
#
# Compare:
#   Baseline 1: Conditional EBLUP + full PEC
#   Baseline 2: Conditional EBLUP + diagonal PEV approximation
#   Method 1:   RAM-N
#   Method 2:   Parametric bootstrap
#   Method 3:   Bayesian
#
# Scenarios:
#   h2 = 0.05, 0.20, 0.40
#   20 outer Monte Carlo replicates per h2 by default
#
# Two experiments:
#   A. Fixed policy
#   B. EBV selection rule: top 5 candidates in each outer replicate
#
# IMPORTANT:
# - RAM-RL is NOT evaluated anywhere in this file.
# - In the selection-rule experiment, the top 5 are reselected in each
#   OUTER replicate, but are held fixed within RAM-N/bootstrap/Bayesian
#   inner calculations. Thus the uncertainty is conditional on the
#   realized EBV-selected set.
# - Bayesian analysis below uses:
#       flat prior for beta,
#       eta = (log sigma2_A, log sigma2_e),
#       independent Normal priors on eta.
#   The adaptive grid is constructed from the posterior itself, not from
#   an RAM-RL normalized likelihood distribution.
#
# REQUIREMENT:
#   First source:
#     Simulation_II_pedigree_animal_model_v2.R
#
# Main commands:
#   fixed_check <- run_fixed_validation_3h2()
#   sel_check   <- run_selection_validation_3h2()
#   all_check   <- run_all_validation_3h2()
# ============================================================


# ------------------------------------------------------------
# Helper: predict EBV for all animals
# ------------------------------------------------------------

predict_ebv_validation <- function(eta, y, prep) {
  sigma2_A <- exp(eta[1])
  sigma2_e <- exp(eta[2])

  X <- prep$X
  A <- prep$A
  ids <- prep$pheno_ids

  ViX <- Vinv_apply(prep, sigma2_A, sigma2_e, X)
  Viy <- Vinv_apply(prep, sigma2_A, sigma2_e, y)

  XtViX <- crossprod(X, ViX)
  beta_hat <- solve(XtViX, crossprod(X, Viy))

  resid <- drop(y - X %*% beta_hat)
  Vir <- Vinv_apply(prep, sigma2_A, sigma2_e, resid)

  uhat <- sigma2_A * drop(
    A[, ids, drop = FALSE] %*% Vir
  )

  list(
    uhat = uhat,
    beta = drop(beta_hat)
  )
}


# ------------------------------------------------------------
# Dedicated Bayesian adaptive grid
# ------------------------------------------------------------

build_bayesian_posterior_grid <- function(
    y, prep,
    fit = NULL,
    ncoarse = 31,
    ngrid = 41,
    log_drop = 14,
    prior_mean = log(c(0.20, 0.80)),
    prior_sd = c(2, 2)) {

  if (is.null(fit)) {
    fit <- fit_reml_animal(y, prep)
  }

  log_prior_eta <- function(eta_A, eta_e) {
    dnorm(
      eta_A,
      mean = prior_mean[1],
      sd = prior_sd[1],
      log = TRUE
    ) +
    dnorm(
      eta_e,
      mean = prior_mean[2],
      sd = prior_sd[2],
      log = TRUE
    )
  }

  # Coarse search over the full REML numerical bounds.
  g1c <- seq(
    fit$lower[1],
    fit$upper[1],
    length.out = ncoarse
  )
  g2c <- seq(
    fit$lower[2],
    fit$upper[2],
    length.out = ncoarse
  )

  coarse <- expand.grid(
    eta_A = g1c,
    eta_e = g2c
  )

  coarse$logL <- vapply(
    seq_len(nrow(coarse)),
    function(i) {
      log_reml_eta(
        as.numeric(
          coarse[i, c("eta_A", "eta_e")]
        ),
        y,
        prep
      )
    },
    numeric(1)
  )

  coarse$logPrior <- log_prior_eta(
    coarse$eta_A,
    coarse$eta_e
  )

  coarse$logPost <- (
    coarse$logL +
    coarse$logPrior
  )

  mx <- max(coarse$logPost)
  keep <- coarse$logPost >= (mx - log_drop)

  if (!any(keep)) {
    stop("Bayesian coarse grid retained no points.")
  }

  d1 <- diff(g1c)[1]
  d2 <- diff(g2c)[1]

  lo1 <- max(
    fit$lower[1],
    min(coarse$eta_A[keep]) - 2 * d1
  )
  hi1 <- min(
    fit$upper[1],
    max(coarse$eta_A[keep]) + 2 * d1
  )

  lo2 <- max(
    fit$lower[2],
    min(coarse$eta_e[keep]) - 2 * d2
  )
  hi2 <- min(
    fit$upper[2],
    max(coarse$eta_e[keep]) + 2 * d2
  )

  g1 <- seq(lo1, hi1, length.out = ngrid)
  g2 <- seq(lo2, hi2, length.out = ngrid)

  grid <- expand.grid(
    eta_A = g1,
    eta_e = g2
  )

  grid$logL <- vapply(
    seq_len(nrow(grid)),
    function(i) {
      log_reml_eta(
        as.numeric(
          grid[i, c("eta_A", "eta_e")]
        ),
        y,
        prep
      )
    },
    numeric(1)
  )

  grid$logPrior <- log_prior_eta(
    grid$eta_A,
    grid$eta_e
  )

  grid$logPost <- (
    grid$logL +
    grid$logPrior
  )

  w <- exp(
    grid$logPost -
    max(grid$logPost)
  )
  w <- w / sum(w)

  list(
    grid = grid,
    weight = w,
    fit = fit,
    prior_mean = prior_mean,
    prior_sd = prior_sd
  )
}


# ----------------------------------------------------------------------
# The earlier validation-era definition of bayesian_animal_validation()
# was removed for the public release.  In the validated cumulative source it
# was overwritten later during source(); only the later production definition
# is retained below.  No production calculation is changed by this removal.
# ----------------------------------------------------------------------

# ------------------------------------------------------------
# Packing and REML diagnostics
# ------------------------------------------------------------

pack_validation_method <- function(x) {
  if (
    is.null(x$valid) ||
    !isTRUE(x$valid)
  ) {
    return(NULL)
  }

  c(
    mean = x$mean,
    sd = x$sd,
    lo = unname(x$interval[1]),
    hi = unname(x$interval[2]),
    p = x$p
  )
}


extract_reml_diagnostics <- function(
    fit,
    true_h2) {

  h2_hat <- (
    fit$theta[1] /
    sum(fit$theta)
  )

  eta_se_A <- NA_real_
  eta_se_e <- NA_real_

  if (!is.null(fit$Veta)) {
    se <- sqrt(diag(fit$Veta))
    eta_se_A <- se[1]
    eta_se_e <- se[2]
  }

  near_lower_A <- (
    fit$eta[1] -
    fit$lower[1]
  ) < 0.05

  c(
    true_h2 = true_h2,
    h2_hat = h2_hat,
    sigma2A_hat = fit$theta[1],
    sigma2e_hat = fit$theta[2],
    hess_min = min(fit$hess_eigen),
    etaA_SE = eta_se_A,
    etae_SE = eta_se_e,
    near_lower_A = as.numeric(
      near_lower_A
    )
  )
}


# ------------------------------------------------------------
# Fixed policy: one replicate
# ------------------------------------------------------------

one_validation_fixed <- function(
    design,
    h2,
    sigma2_P = 1,
    beta = 0,
    target = 0,
    M = 500,
    B = 100,
    ngrid = 41,
    seed_data,
    seed_RAMN,
    seed_boot,
    seed_bayes) {

  dat <- simulate_animal_data(
    design = design,
    h2 = h2,
    sigma2_P = sigma2_P,
    beta = beta,
    seed = seed_data
  )

  fit0 <- fit_reml_animal(
    dat$y,
    design$prep
  )

  cond <- conditional_method(
    dat$y,
    design$prep,
    design$cstar,
    target = target
  )

  set.seed(seed_RAMN)
  rn <- RAM_N_animal(
    dat$y,
    design$prep,
    design$cstar,
    M = M,
    target = target
  )

  set.seed(seed_boot)
  pb <- parametric_bootstrap_animal(
    dat$y,
    design$prep,
    design$cstar,
    B = B,
    target = target
  )

  by <- bayesian_animal_validation(
    dat$y,
    design$prep,
    design$cstar,
    M = M,
    ngrid = ngrid,
    target = target,
    seed_inner = seed_bayes
  )

  list(
    truth = dat$truth,

    conditional =
      pack_validation_method(cond),

    conditional_diag = c(
      mean = cond$mean,
      sd = cond$sd_diag,
      lo = (
        cond$mean -
        1.96 * cond$sd_diag
      ),
      hi = (
        cond$mean +
        1.96 * cond$sd_diag
      ),
      p = 1 - pnorm(
        target,
        mean = cond$mean,
        sd = max(
          cond$sd_diag,
          1e-12
        )
      )
    ),

    RAM_N =
      pack_validation_method(rn),

    bootstrap =
      pack_validation_method(pb),

    Bayesian =
      pack_validation_method(by),

    RAM_N_reasons =
      if (
        !is.null(rn$valid) &&
        !rn$valid
      ) {
        rn$reasons
      } else {
        character()
      },

    reml_diag =
      extract_reml_diagnostics(
        fit0,
        true_h2 = h2
      ),

    bayes_mean_h2 =
      by$mean_h2_post
  )
}


# ------------------------------------------------------------
# Selection rule: one replicate
# ------------------------------------------------------------

one_validation_selection <- function(
    design,
    h2,
    sigma2_P = 1,
    beta = 0,
    target = 0,
    n_select = 5,
    M = 500,
    B = 100,
    ngrid = 41,
    seed_data,
    seed_RAMN,
    seed_boot,
    seed_bayes) {

  dat <- simulate_animal_data(
    design = design,
    h2 = h2,
    sigma2_P = sigma2_P,
    beta = beta,
    seed = seed_data
  )

  fit_sel <- fit_reml_animal(
    dat$y,
    design$prep
  )

  pred <- predict_ebv_validation(
    fit_sel$eta,
    dat$y,
    design$prep
  )

  candidates <- design$candidates

  ord <- order(
    pred$uhat[candidates],
    decreasing = TRUE
  )

  selected <- candidates[
    ord[seq_len(n_select)]
  ]

  cstar <- make_policy(
    n_animals = nrow(design$ped),
    candidates = candidates,
    selected = selected
  )

  truth <- sum(
    cstar * dat$u
  )

  # Oracle top-5 by true breeding value.
  ord_true <- order(
    dat$u[candidates],
    decreasing = TRUE
  )

  oracle_selected <- candidates[
    ord_true[seq_len(n_select)]
  ]

  c_oracle <- make_policy(
    n_animals = nrow(design$ped),
    candidates = candidates,
    selected = oracle_selected
  )

  oracle_truth <- sum(
    c_oracle * dat$u
  )

  overlap_n <- length(
    intersect(
      selected,
      oracle_selected
    )
  )

  overlap_rate <- (
    overlap_n /
    n_select
  )

  jaccard <- (
    length(
      intersect(
        selected,
        oracle_selected
      )
    ) /
    length(
      union(
        selected,
        oracle_selected
      )
    )
  )

  cond <- conditional_method(
    dat$y,
    design$prep,
    cstar,
    target = target
  )

  set.seed(seed_RAMN)
  rn <- RAM_N_animal(
    dat$y,
    design$prep,
    cstar,
    M = M,
    target = target
  )

  set.seed(seed_boot)
  pb <- parametric_bootstrap_animal(
    dat$y,
    design$prep,
    cstar,
    B = B,
    target = target
  )

  by <- bayesian_animal_validation(
    dat$y,
    design$prep,
    cstar,
    M = M,
    ngrid = ngrid,
    target = target,
    seed_inner = seed_bayes
  )

  list(
    truth = truth,

    selected = selected,
    oracle_selected = oracle_selected,

    oracle_truth = oracle_truth,
    oracle_gap = (
      oracle_truth -
      truth
    ),

    overlap_rate = overlap_rate,
    jaccard = jaccard,

    conditional =
      pack_validation_method(cond),

    conditional_diag = c(
      mean = cond$mean,
      sd = cond$sd_diag,
      lo = (
        cond$mean -
        1.96 * cond$sd_diag
      ),
      hi = (
        cond$mean +
        1.96 * cond$sd_diag
      ),
      p = 1 - pnorm(
        target,
        mean = cond$mean,
        sd = max(
          cond$sd_diag,
          1e-12
        )
      )
    ),

    RAM_N =
      pack_validation_method(rn),

    bootstrap =
      pack_validation_method(pb),

    Bayesian =
      pack_validation_method(by),

    RAM_N_reasons =
      if (
        !is.null(rn$valid) &&
        !rn$valid
      ) {
        rn$reasons
      } else {
        character()
      },

    reml_diag =
      extract_reml_diagnostics(
        fit_sel,
        true_h2 = h2
      ),

    bayes_mean_h2 =
      by$mean_h2_post
  )
}


# ------------------------------------------------------------
# Summaries
# ------------------------------------------------------------

summarize_methods_validation <- function(
    res,
    target = 0) {

  methods <- c(
    "conditional",
    "conditional_diag",
    "RAM_N",
    "bootstrap",
    "Bayesian"
  )

  tab <- t(
    vapply(
      methods,
      function(z) {
        summarize_simII(
          res,
          z,
          target = target
        )
      },
      numeric(10)
    )
  )

  rownames(tab) <- methods

  tab <- cbind(
    tab,
    valid_rate = (
      tab[, "valid"] /
      length(res)
    )
  )

  tab
}


summarize_reml_validation <- function(res) {

  D <- do.call(
    rbind,
    lapply(
      res,
      `[[`,
      "reml_diag"
    )
  )

  ramn_valid <- !vapply(
    lapply(
      res,
      `[[`,
      "RAM_N"
    ),
    is.null,
    logical(1)
  )

  c(
    mean_h2_hat =
      mean(D[, "h2_hat"]),

    median_h2_hat =
      median(D[, "h2_hat"]),

    mean_sigma2A_hat =
      mean(D[, "sigma2A_hat"]),

    mean_sigma2e_hat =
      mean(D[, "sigma2e_hat"]),

    mean_hess_min =
      mean(D[, "hess_min"]),

    median_etaA_SE =
      median(
        D[, "etaA_SE"],
        na.rm = TRUE
      ),

    lower_boundary_rate =
      mean(
        D[, "near_lower_A"]
      ),

    RAM_N_valid_rate =
      mean(ramn_valid),

    mean_Bayes_h2 =
      mean(
        vapply(
          res,
          `[[`,
          numeric(1),
          "bayes_mean_h2"
        )
      )
  )
}


summarize_selection_validation <- function(res) {

  c(
    mean_true_gain =
      mean(
        vapply(
          res,
          `[[`,
          numeric(1),
          "truth"
        )
      ),

    mean_oracle_gain =
      mean(
        vapply(
          res,
          `[[`,
          numeric(1),
          "oracle_truth"
        )
      ),

    mean_oracle_gap =
      mean(
        vapply(
          res,
          `[[`,
          numeric(1),
          "oracle_gap"
        )
      ),

    mean_overlap_rate =
      mean(
        vapply(
          res,
          `[[`,
          numeric(1),
          "overlap_rate"
        )
      ),

    mean_jaccard =
      mean(
        vapply(
          res,
          `[[`,
          numeric(1),
          "jaccard"
        )
      )
  )
}


# ------------------------------------------------------------
# Core runner for one h2
# ------------------------------------------------------------

run_one_h2_fixed_validation <- function(
    h2,
    S = 20,
    design = build_simII_design(),
    M = 500,
    B = 100,
    ngrid = 41,
    target = 0,
    seed = 20260830,
    keep_raw = FALSE) {

  set.seed(seed)

  seeds <- matrix(
    sample.int(
      .Machine$integer.max,
      S * 4
    ),
    nrow = S,
    ncol = 4
  )

  res <- vector(
    "list",
    S
  )

  for (s in seq_len(S)) {

    res[[s]] <- one_validation_fixed(
      design = design,
      h2 = h2,
      target = target,
      M = M,
      B = B,
      ngrid = ngrid,
      seed_data = seeds[s, 1],
      seed_RAMN = seeds[s, 2],
      seed_boot = seeds[s, 3],
      seed_bayes = seeds[s, 4]
    )

    if (s %% 5 == 0) {
      cat(
        "Fixed h2 =",
        h2,
        ": completed",
        s,
        "of",
        S,
        "\n"
      )
    }
  }

  out <- list(
    h2 = h2,
    methods =
      summarize_methods_validation(
        res,
        target = target
      ),
    diagnostics =
      summarize_reml_validation(
        res
      )
  )

  if (keep_raw) {
    out$raw <- res
  }

  out
}


run_one_h2_selection_validation <- function(
    h2,
    S = 20,
    design = build_simII_design(),
    n_select = 5,
    M = 500,
    B = 100,
    ngrid = 41,
    target = 0,
    seed = 20260830,
    keep_raw = FALSE) {

  set.seed(seed)

  seeds <- matrix(
    sample.int(
      .Machine$integer.max,
      S * 4
    ),
    nrow = S,
    ncol = 4
  )

  res <- vector(
    "list",
    S
  )

  for (s in seq_len(S)) {

    res[[s]] <- one_validation_selection(
      design = design,
      h2 = h2,
      target = target,
      n_select = n_select,
      M = M,
      B = B,
      ngrid = ngrid,
      seed_data = seeds[s, 1],
      seed_RAMN = seeds[s, 2],
      seed_boot = seeds[s, 3],
      seed_bayes = seeds[s, 4]
    )

    if (s %% 5 == 0) {
      cat(
        "Selection h2 =",
        h2,
        ": completed",
        s,
        "of",
        S,
        "\n"
      )
    }
  }

  out <- list(
    h2 = h2,
    methods =
      summarize_methods_validation(
        res,
        target = target
      ),
    diagnostics =
      summarize_reml_validation(
        res
      ),
    selection =
      summarize_selection_validation(
        res
      )
  )

  if (keep_raw) {
    out$raw <- res
  }

  out
}


# ------------------------------------------------------------
# Three-h2 runners
# ------------------------------------------------------------

run_fixed_validation_3h2 <- function(
    h2_values = c(0.05, 0.20, 0.40),
    S = 20,
    M = 500,
    B = 100,
    ngrid = 41,
    target = 0,
    seed = 20260831,
    keep_raw = FALSE) {

  design <- build_simII_design()

  cat(
    "\n--- FIXED POLICY: candidate relationship summary ---\n"
  )
  print(
    relationship_summary(
      design$A,
      design$candidates
    )
  )

  ans <- vector(
    "list",
    length(h2_values)
  )
  names(ans) <- paste0(
    "h2_",
    h2_values
  )

  for (j in seq_along(h2_values)) {

    h2 <- h2_values[j]

    cat(
      "\n===== FIXED POLICY: h2 =",
      h2,
      "=====\n"
    )

    ans[[j]] <- run_one_h2_fixed_validation(
      h2 = h2,
      S = S,
      design = design,
      M = M,
      B = B,
      ngrid = ngrid,
      target = target,
      seed = seed + j,
      keep_raw = keep_raw
    )

    print(
      ans[[j]]$methods
    )

    cat(
      "\nREML / RAM-N diagnostics:\n"
    )
    print(
      ans[[j]]$diagnostics
    )
  }

  invisible(
    list(
      design = design,
      result = ans
    )
  )
}


run_selection_validation_3h2 <- function(
    h2_values = c(0.05, 0.20, 0.40),
    S = 20,
    n_select = 5,
    M = 500,
    B = 100,
    ngrid = 41,
    target = 0,
    seed = 20260841,
    keep_raw = FALSE) {

  design <- build_simII_design()

  cat(
    "\n--- SELECTION RULE: candidate relationship summary ---\n"
  )
  print(
    relationship_summary(
      design$A,
      design$candidates
    )
  )

  ans <- vector(
    "list",
    length(h2_values)
  )
  names(ans) <- paste0(
    "h2_",
    h2_values
  )

  for (j in seq_along(h2_values)) {

    h2 <- h2_values[j]

    cat(
      "\n===== SELECTION RULE: h2 =",
      h2,
      "=====\n"
    )

    ans[[j]] <- run_one_h2_selection_validation(
      h2 = h2,
      S = S,
      design = design,
      n_select = n_select,
      M = M,
      B = B,
      ngrid = ngrid,
      target = target,
      seed = seed + j,
      keep_raw = keep_raw
    )

    print(
      ans[[j]]$methods
    )

    cat(
      "\nREML / RAM-N diagnostics:\n"
    )
    print(
      ans[[j]]$diagnostics
    )

    cat(
      "\nSelection performance:\n"
    )
    print(
      ans[[j]]$selection
    )
  }

  invisible(
    list(
      design = design,
      result = ans
    )
  )
}


# ------------------------------------------------------------
# Combined convenience runner
# ------------------------------------------------------------

run_all_validation_3h2 <- function(
    h2_values = c(0.05, 0.20, 0.40),
    S = 20,
    n_select = 5,
    M = 500,
    B = 100,
    ngrid = 41,
    target = 0,
    seed_fixed = 20260831,
    seed_selection = 20260841,
    keep_raw = FALSE) {

  fixed <- run_fixed_validation_3h2(
    h2_values = h2_values,
    S = S,
    M = M,
    B = B,
    ngrid = ngrid,
    target = target,
    seed = seed_fixed,
    keep_raw = keep_raw
  )

  selection <- run_selection_validation_3h2(
    h2_values = h2_values,
    S = S,
    n_select = n_select,
    M = M,
    B = B,
    ngrid = ngrid,
    target = target,
    seed = seed_selection,
    keep_raw = keep_raw
  )

  invisible(
    list(
      fixed = fixed,
      selection = selection
    )
  )
}


# ------------------------------------------------------------
# Compact comparison tables after the run
# ------------------------------------------------------------

make_method_table_3h2 <- function(
    validation_object) {

  x <- validation_object$result

  rows <- list()

  k <- 1L

  for (j in seq_along(x)) {

    h2 <- x[[j]]$h2
    tab <- x[[j]]$methods

    for (m in rownames(tab)) {

      rows[[k]] <- data.frame(
        h2 = h2,
        method = m,
        valid = tab[m, "valid"],
        valid_rate = tab[m, "valid_rate"],
        bias = tab[m, "bias"],
        RMSE = tab[m, "RMSE"],
        mean_SE = tab[m, "mean_SE"],
        MSE_ratio = tab[m, "MSE_ratio"],
        inclusion = tab[m, "inclusion"],
        width = tab[m, "width"],
        Brier = tab[m, "Brier"],
        mean_p = tab[m, "mean_p"],
        event_rate = tab[m, "event_rate"],
        row.names = NULL
      )

      k <- k + 1L
    }
  }

  do.call(
    rbind,
    rows
  )
}


make_diagnostic_table_3h2 <- function(
    validation_object) {

  x <- validation_object$result

  do.call(
    rbind,
    lapply(
      x,
      function(z) {
        data.frame(
          h2 = z$h2,
          t(z$diagnostics),
          row.names = NULL
        )
      }
    )
  )
}


make_selection_table_3h2 <- function(
    selection_validation_object) {

  x <- selection_validation_object$result

  do.call(
    rbind,
    lapply(
      x,
      function(z) {
        data.frame(
          h2 = z$h2,
          t(z$selection),
          row.names = NULL
        )
      }
    )
  )
}


# ============================================================
# HOW TO RUN
# ============================================================
#
# 1. Source the base code:
#
#
#    Choose:
#    Simulation_II_pedigree_animal_model_v2.R
#
#
# 2. Source this validation code:
#
#
#    Choose:
#    Simulation_II_validation_3h2_no_RAMRL.R
#
#
# 3A. Fixed policy only:
#
#    fixed_check <- run_fixed_validation_3h2()
#
#    fixed_table <- make_method_table_3h2(fixed_check)
#    fixed_diag  <- make_diagnostic_table_3h2(fixed_check)
#
#    print(fixed_table)
#    print(fixed_diag)
#
#
# 3B. EBV selection rule only:
#
#    sel_check <- run_selection_validation_3h2()
#
#    sel_table <- make_method_table_3h2(sel_check)
#    sel_diag  <- make_diagnostic_table_3h2(sel_check)
#    sel_perf  <- make_selection_table_3h2(sel_check)
#
#    print(sel_table)
#    print(sel_diag)
#    print(sel_perf)
#
#
# 3C. Run both:
#
#    all_check <- run_all_validation_3h2()
#
#    fixed_table <- make_method_table_3h2(
#      all_check$fixed
#    )
#
#    sel_table <- make_method_table_3h2(
#      all_check$selection
#    )
#
#    sel_perf <- make_selection_table_3h2(
#      all_check$selection
#    )
#
#
# Interpretation notes:
#
# - 20 replicates are for code/behavior validation only.
# - Do NOT infer final rankings from 20 replicates.
# - RAM-N performance rows are conditional on diagnostic-valid
#   replicates; always interpret valid_rate together with coverage.
# - conditional vs conditional_diag isolates the role of PEC.
# - Selection-rule intervals remain conditional on the realized
#   top-5 EBV selected set within each outer replicate.
# - For target = 0, Brier score may become uninformative if the
#   event rate is near 0 or 1, especially under selection.
#
# ============================================================

# --- Embedded validated realised-policy 3x3 helpers ---

# ============================================================
# Simulation II-B DEEP DIVE
# Realised EBV-selection policy in a pedigree animal model
#
# PURPOSE
# -------
# In every OUTER Monte Carlo replicate:
#   1. generate genetic values and phenotypes
#   2. estimate variance components by REML
#   3. compute EBLUPs
#   4. select the top 5 candidate males by EBLUP
#   5. define the realised policy c(y)
#   6. hold that realised selected set fixed for uncertainty analysis
#   7. evaluate the true realised gain c(y)'u
#
# This is Simulation II-B:
#   "realised EBV-selection policy"
#
# METHODS
# -------
# Baseline:
#   conditional       = conditional EBLUP + FULL PEC
#   conditional_diag  = conditional EBLUP + diagonal PEV only
#
# Main uncertainty methods:
#   RAM_N_policy
#   bootstrap_policy
#   Bayesian_policy
#
# RAM-RL is NOT used.
# Internal reselection bootstrap is NOT used in the main comparison.
#
#
# FACTORIAL DESIGN
# ----------------
# Heritability:
#   h2 = 0.05, 0.20, 0.40
#
# Information level:
#   low       : all candidate own records + 40% of noncandidate records
#   moderate  : all candidate own records + 70% of noncandidate records
#   high      : all candidate own records + 100% of noncandidate records
#
# With the default pedigree:
#   candidates = 20 males
#   full nonfounder phenotype set = 120 animals
#
# Therefore the default information levels contain approximately:
#   low      = 20 + 40 = 60 records
#   moderate = 20 + 70 = 90 records
#   high     = 20 + 100 = 120 records
#
# IMPORTANT:
# - The same pedigree A is used in all scenarios.
# - The same candidate set is used in all scenarios.
# - Phenotype subsets are NESTED and fixed before simulation.
# - Candidate own records are retained at all information levels.
# - Common latent genetic and residual draws are used across all
#   h2 x information scenarios within an outer replicate.
#
# This paired/common-random-number construction greatly improves
# interpretability of differences among scenarios.
#
#
# REQUIREMENTS
# ------------
# Source these first:
#   1) Simulation_II_pedigree_animal_model_v2.R
#   2) Simulation_II_validation_3h2_no_RAMRL.R
#
#
# MAIN VALIDATION COMMAND
# -----------------------
#   deep20 <- run_realised_policy_deep_validation(S = 20)
#
# ============================================================


# ------------------------------------------------------------
# 1. Build nested information-level designs
# ------------------------------------------------------------

build_information_level_designs <- function(
    base_design = build_simII_design(),
    info_fractions = c(
      low = 0.40,
      moderate = 0.70,
      high = 1.00
    ),
    seed_info = 20260850) {

  if (is.null(names(info_fractions))) {
    stop("info_fractions must be a named numeric vector.")
  }

  if (any(info_fractions <= 0) ||
      any(info_fractions > 1)) {
    stop("All information fractions must be in (0, 1].")
  }

  # Require increasing nested information.
  ord <- order(info_fractions)
  info_fractions <- info_fractions[ord]

  full_ids <- base_design$pheno_ids
  candidates <- base_design$candidates

  # Candidate males retain own phenotype records at every level.
  candidate_records <- intersect(
    candidates,
    full_ids
  )

  noncandidate_pool <- setdiff(
    full_ids,
    candidate_records
  )

  set.seed(seed_info)

  # One fixed random order, independent of all phenotype/BV simulations.
  noncandidate_order <- sample(
    noncandidate_pool,
    length(noncandidate_pool),
    replace = FALSE
  )

  out <- vector(
    "list",
    length(info_fractions)
  )

  names(out) <- names(info_fractions)

  for (j in seq_along(info_fractions)) {

    frac <- info_fractions[j]

    n_non <- round(
      frac *
      length(noncandidate_pool)
    )

    n_non <- max(
      1L,
      min(
        n_non,
        length(noncandidate_pool)
      )
    )

    ids <- c(
      candidate_records,
      noncandidate_order[
        seq_len(n_non)
      ]
    )

    # Put phenotype IDs back into pedigree order.
    ids <- sort(
      unique(ids)
    )

    X <- matrix(
      1,
      nrow = length(ids),
      ncol = 1
    )

    colnames(X) <- "intercept"

    prep <- prepare_animal_model(
      A = base_design$A,
      pheno_ids = ids,
      X = X
    )

    out[[j]] <- list(
      name = names(info_fractions)[j],
      fraction = frac,
      pheno_ids = ids,
      n_pheno = length(ids),
      n_candidate_records =
        length(candidate_records),
      n_noncandidate_records =
        n_non,
      X = X,
      prep = prep
    )
  }

  list(
    base_design = base_design,
    info_levels = out,
    full_pheno_ids = full_ids,
    candidate_records = candidate_records,
    noncandidate_order = noncandidate_order
  )
}


print_information_design <- function(info_bundle) {

  cat("\n--- Information-level design ---\n")

  tab <- do.call(
    rbind,
    lapply(
      info_bundle$info_levels,
      function(z) {
        data.frame(
          information = z$name,
          fraction_noncandidate = z$fraction,
          n_pheno = z$n_pheno,
          candidate_records = z$n_candidate_records,
          noncandidate_records =
            z$n_noncandidate_records,
          row.names = NULL
        )
      }
    )
  )

  print(tab)

  invisible(tab)
}


# ------------------------------------------------------------
# 2. Common latent Monte Carlo replicate
# ------------------------------------------------------------

simulate_latent_replicate <- function(
    info_bundle,
    seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  A <- info_bundle$base_design$A

  # Standardized additive genetic vector:
  # g ~ N(0, A)
  g <- rmvn1(
    rep(0, nrow(A)),
    A
  )

  # Standardized residuals for the largest phenotype set.
  full_ids <- info_bundle$full_pheno_ids

  z_e <- rnorm(
    length(full_ids)
  )

  names(z_e) <- as.character(
    full_ids
  )

  list(
    g = g,
    z_e = z_e
  )
}


make_dataset_from_latent <- function(
    latent,
    info_level,
    sigma2_P = 1,
    h2 = 0.20,
    beta = 0) {

  if (h2 <= 0 || h2 >= 1) {
    stop("h2 must be between 0 and 1.")
  }

  sigma2_A <- (
    h2 *
    sigma2_P
  )

  sigma2_e <- (
    (1 - h2) *
    sigma2_P
  )

  u <- (
    sqrt(sigma2_A) *
    latent$g
  )

  ids <- info_level$pheno_ids

  z <- latent$z_e[
    as.character(ids)
  ]

  if (anyNA(z)) {
    stop(
      "Residual draw could not be matched to phenotype IDs."
    )
  }

  e <- (
    sqrt(sigma2_e) *
    z
  )

  y <- (
    beta +
    u[ids] +
    e
  )

  list(
    y = as.numeric(y),
    u = u,
    sigma2_A = sigma2_A,
    sigma2_e = sigma2_e
  )
}


# ------------------------------------------------------------
# 3. Safe correlations and relationship summaries
# ------------------------------------------------------------

safe_cor <- function(
    x,
    y,
    method = "pearson") {

  if (
    length(x) < 3 ||
    length(y) < 3 ||
    sd(x) == 0 ||
    sd(y) == 0
  ) {
    return(NA_real_)
  }

  suppressWarnings(
    cor(
      x,
      y,
      method = method
    )
  )
}


mean_offdiag_A <- function(
    A,
    ids) {

  if (length(ids) < 2) {
    return(NA_real_)
  }

  S <- A[
    ids,
    ids,
    drop = FALSE
  ]

  mean(
    S[
      upper.tri(S)
    ]
  )
}


# ------------------------------------------------------------
# 4. One realised-policy scenario
# ------------------------------------------------------------

one_realised_policy_scenario <- function(
    latent,
    info_bundle,
    info_level,
    h2,
    sigma2_P = 1,
    beta = 0,
    target = 0,
    n_select = 5,
    M = 500,
    B = 100,
    ngrid = 41,
    seed_RAMN,
    seed_boot,
    seed_bayes) {

  base <- info_bundle$base_design

  dat <- make_dataset_from_latent(
    latent = latent,
    info_level = info_level,
    sigma2_P = sigma2_P,
    h2 = h2,
    beta = beta
  )

  # ----------------------------------------------------------
  # REML + EBLUP selection
  # ----------------------------------------------------------

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

  ord <- order(
    pred$uhat[candidates],
    decreasing = TRUE
  )

  selected <- candidates[
    ord[
      seq_len(n_select)
    ]
  ]

  cstar <- make_policy(
    n_animals = nrow(base$ped),
    candidates = candidates,
    selected = selected
  )

  # ----------------------------------------------------------
  # True realised gain
  # ----------------------------------------------------------

  truth <- sum(
    cstar *
    dat$u
  )

  # ----------------------------------------------------------
  # Oracle benchmark
  # ----------------------------------------------------------

  ord_true <- order(
    dat$u[candidates],
    decreasing = TRUE
  )

  oracle_selected <- candidates[
    ord_true[
      seq_len(n_select)
    ]
  ]

  c_oracle <- make_policy(
    n_animals = nrow(base$ped),
    candidates = candidates,
    selected = oracle_selected
  )

  oracle_truth <- sum(
    c_oracle *
    dat$u
  )

  overlap_n <- length(
    intersect(
      selected,
      oracle_selected
    )
  )

  overlap_rate <- (
    overlap_n /
    n_select
  )

  jaccard <- (
    overlap_n /
    length(
      union(
        selected,
        oracle_selected
      )
    )
  )

  exact_top_set <- as.numeric(
    setequal(
      selected,
      oracle_selected
    )
  )

  # Candidate-level accuracy diagnostics.
  candidate_accuracy <- safe_cor(
    pred$uhat[candidates],
    dat$u[candidates],
    method = "pearson"
  )

  candidate_rank_accuracy <- safe_cor(
    pred$uhat[candidates],
    dat$u[candidates],
    method = "spearman"
  )

  # ----------------------------------------------------------
  # Conditional EBLUP with full PEC
  # ----------------------------------------------------------

  cond <- conditional_method(
    y = dat$y,
    prep = info_level$prep,
    cstar = cstar,
    target = target
  )

  # ----------------------------------------------------------
  # RAM-N: realised-policy conditional
  # ----------------------------------------------------------

  set.seed(seed_RAMN)

  rn <- RAM_N_animal(
    y = dat$y,
    prep = info_level$prep,
    cstar = cstar,
    M = M,
    target = target
  )

  # ----------------------------------------------------------
  # Parametric bootstrap: realised policy fixed internally
  # ----------------------------------------------------------

  set.seed(seed_boot)

  pb <- parametric_bootstrap_animal(
    y = dat$y,
    prep = info_level$prep,
    cstar = cstar,
    B = B,
    target = target
  )

  # ----------------------------------------------------------
  # Bayesian: realised policy fixed internally
  # ----------------------------------------------------------

  by <- bayesian_animal_validation(
    y = dat$y,
    prep = info_level$prep,
    cstar = cstar,
    M = M,
    ngrid = ngrid,
    target = target,
    seed_inner = seed_bayes
  )

  # ----------------------------------------------------------
  # PEC diagnostics
  # ----------------------------------------------------------

  v_full <- (
    cond$sd^2
  )

  v_diag <- (
    cond$sd_diag^2
  )

  pec_ratio_diag_over_full <- if (
    v_full > 0
  ) {
    v_diag / v_full
  } else {
    NA_real_
  }

  # ----------------------------------------------------------
  # Selection structure
  # ----------------------------------------------------------

  selected_mean_relationship <- mean_offdiag_A(
    base$A,
    selected
  )

  selected_mean_inbreeding <- mean(
    diag(base$A)[selected] - 1
  )

  # ----------------------------------------------------------
  # Pack outputs
  # ----------------------------------------------------------

  list(
    truth = truth,

    h2 = h2,
    information = info_level$name,
    n_pheno = info_level$n_pheno,

    selected = selected,
    oracle_selected = oracle_selected,

    oracle_truth = oracle_truth,
    oracle_gap = (
      oracle_truth -
      truth
    ),

    overlap_rate = overlap_rate,
    jaccard = jaccard,
    exact_top_set = exact_top_set,

    candidate_accuracy =
      candidate_accuracy,

    candidate_rank_accuracy =
      candidate_rank_accuracy,

    predicted_gain =
      cond$mean,

    selection_optimism =
      cond$mean -
      truth,

    pec_ratio_diag_over_full =
      pec_ratio_diag_over_full,

    selected_mean_relationship =
      selected_mean_relationship,

    selected_mean_inbreeding =
      selected_mean_inbreeding,

    conditional =
      pack_validation_method(
        cond
      ),

    conditional_diag = c(
      mean = cond$mean,
      sd = cond$sd_diag,
      lo = (
        cond$mean -
        1.96 * cond$sd_diag
      ),
      hi = (
        cond$mean +
        1.96 * cond$sd_diag
      ),
      p = 1 - pnorm(
        target,
        mean = cond$mean,
        sd = max(
          cond$sd_diag,
          1e-12
        )
      )
    ),

    RAM_N_policy =
      pack_validation_method(
        rn
      ),

    bootstrap_policy =
      pack_validation_method(
        pb
      ),

    Bayesian_policy =
      pack_validation_method(
        by
      ),

    RAM_N_reasons =
      if (
        !is.null(rn$valid) &&
        !rn$valid
      ) {
        rn$reasons
      } else {
        character()
      },

    reml_diag =
      extract_reml_diagnostics(
        fit_sel,
        true_h2 = h2
      ),

    bayes_mean_h2 =
      by$mean_h2_post
  )
}


# ------------------------------------------------------------
# 5. Summaries for one h2 x information scenario
# ------------------------------------------------------------

summarize_realised_policy_methods <- function(
    res,
    target = 0) {

  methods <- c(
    "conditional",
    "conditional_diag",
    "RAM_N_policy",
    "bootstrap_policy",
    "Bayesian_policy"
  )

  tab <- t(
    vapply(
      methods,
      function(z) {
        summarize_simII(
          res,
          z,
          target = target
        )
      },
      numeric(10)
    )
  )

  rownames(tab) <- methods

  cbind(
    tab,
    valid_rate =
      tab[, "valid"] /
      length(res)
  )
}


summarize_realised_policy_selection <- function(
    res) {

  get_num <- function(name) {
    vapply(
      res,
      `[[`,
      numeric(1),
      name
    )
  }

  c(
    mean_true_gain =
      mean(
        get_num("truth")
      ),

    mean_predicted_gain =
      mean(
        get_num("predicted_gain")
      ),

    mean_selection_optimism =
      mean(
        get_num("selection_optimism")
      ),

    mean_oracle_gain =
      mean(
        get_num("oracle_truth")
      ),

    mean_oracle_gap =
      mean(
        get_num("oracle_gap")
      ),

    mean_overlap_rate =
      mean(
        get_num("overlap_rate")
      ),

    mean_jaccard =
      mean(
        get_num("jaccard")
      ),

    exact_top5_rate =
      mean(
        get_num("exact_top_set")
      ),

    mean_candidate_accuracy =
      mean(
        get_num("candidate_accuracy"),
        na.rm = TRUE
      ),

    mean_candidate_rank_accuracy =
      mean(
        get_num("candidate_rank_accuracy"),
        na.rm = TRUE
      ),

    mean_PEC_diag_full_ratio =
      mean(
        get_num(
          "pec_ratio_diag_over_full"
        ),
        na.rm = TRUE
      ),

    mean_selected_relationship =
      mean(
        get_num(
          "selected_mean_relationship"
        ),
        na.rm = TRUE
      ),

    mean_selected_inbreeding =
      mean(
        get_num(
          "selected_mean_inbreeding"
        ),
        na.rm = TRUE
      )
  )
}


summarize_realised_policy_reml <- function(
    res) {

  D <- do.call(
    rbind,
    lapply(
      res,
      `[[`,
      "reml_diag"
    )
  )

  ramn_valid <- !vapply(
    lapply(
      res,
      `[[`,
      "RAM_N_policy"
    ),
    is.null,
    logical(1)
  )

  c(
    mean_h2_hat =
      mean(
        D[, "h2_hat"]
      ),

    median_h2_hat =
      median(
        D[, "h2_hat"]
      ),

    mean_sigma2A_hat =
      mean(
        D[, "sigma2A_hat"]
      ),

    mean_sigma2e_hat =
      mean(
        D[, "sigma2e_hat"]
      ),

    mean_hess_min =
      mean(
        D[, "hess_min"]
      ),

    median_etaA_SE =
      median(
        D[, "etaA_SE"],
        na.rm = TRUE
      ),

    lower_boundary_rate =
      mean(
        D[, "near_lower_A"]
      ),

    RAM_N_valid_rate =
      mean(
        ramn_valid
      ),

    mean_Bayes_h2 =
      mean(
        vapply(
          res,
          `[[`,
          numeric(1),
          "bayes_mean_h2"
        )
      )
  )
}


# ------------------------------------------------------------
# 6. Main deep-validation runner
# ------------------------------------------------------------

run_realised_policy_deep_validation <- function(
    S = 20,
    h2_values = c(
      0.05,
      0.20,
      0.40
    ),
    info_fractions = c(
      low = 0.40,
      moderate = 0.70,
      high = 1.00
    ),
    n_select = 5,
    sigma2_P = 1,
    beta = 0,
    target = 0,
    M = 500,
    B = 100,
    ngrid = 41,
    seed_info = 20260850,
    seed_outer = 20260851,
    keep_raw = FALSE) {

  base <- build_simII_design()

  info_bundle <- build_information_level_designs(
    base_design = base,
    info_fractions = info_fractions,
    seed_info = seed_info
  )

  cat(
    "\n====================================================\n"
  )
  cat(
    "Simulation II-B: realised EBV-selection policy\n"
  )
  cat(
    "====================================================\n"
  )

  cat(
    "\nCandidate relationship summary:\n"
  )

  print(
    relationship_summary(
      base$A,
      base$candidates
    )
  )

  print_information_design(
    info_bundle
  )

  # ----------------------------------------------------------
  # Common outer latent replicates
  # ----------------------------------------------------------

  set.seed(seed_outer)

  seed_latent <- sample.int(
    .Machine$integer.max,
    S
  )

  # Separate method seeds for every outer x h2 x info scenario.
  n_h <- length(h2_values)
  n_i <- length(
    info_bundle$info_levels
  )

  seed_methods <- array(
    sample.int(
      .Machine$integer.max,
      S * n_h * n_i * 3
    ),
    dim = c(
      S,
      n_h,
      n_i,
      3
    )
  )

  latent_list <- lapply(
    seed_latent,
    function(s) {
      simulate_latent_replicate(
        info_bundle,
        seed = s
      )
    }
  )

  scenario_results <- list()

  counter <- 1L

  # ----------------------------------------------------------
  # Loop over h2 x information
  # ----------------------------------------------------------

  for (ih in seq_along(h2_values)) {

    h2 <- h2_values[ih]

    for (ii in seq_along(
      info_bundle$info_levels
    )) {

      info <- info_bundle$info_levels[[ii]]

      cat(
        "\n===== h2 =",
        h2,
        "| information =",
        info$name,
        "| n_pheno =",
        info$n_pheno,
        "=====\n"
      )

      res <- vector(
        "list",
        S
      )

      for (s in seq_len(S)) {

        res[[s]] <-
          one_realised_policy_scenario(
            latent =
              latent_list[[s]],

            info_bundle =
              info_bundle,

            info_level =
              info,

            h2 =
              h2,

            sigma2_P =
              sigma2_P,

            beta =
              beta,

            target =
              target,

            n_select =
              n_select,

            M =
              M,

            B =
              B,

            ngrid =
              ngrid,

            seed_RAMN =
              seed_methods[
                s,
                ih,
                ii,
                1
              ],

            seed_boot =
              seed_methods[
                s,
                ih,
                ii,
                2
              ],

            seed_bayes =
              seed_methods[
                s,
                ih,
                ii,
                3
              ]
          )

        if (s %% 5 == 0) {
          cat(
            "Completed",
            s,
            "of",
            S,
            "\n"
          )
        }
      }

      key <- paste0(
        "h2_",
        h2,
        "__",
        info$name
      )

      scenario_results[[key]] <- list(
        h2 = h2,
        information = info$name,
        n_pheno = info$n_pheno,

        methods =
          summarize_realised_policy_methods(
            res,
            target = target
          ),

        selection =
          summarize_realised_policy_selection(
            res
          ),

        diagnostics =
          summarize_realised_policy_reml(
            res
          )
      )

      if (keep_raw) {
        scenario_results[[key]]$raw <- res
      }

      counter <- counter + 1L
    }
  }

  out <- list(
    design = base,
    information_design = info_bundle,
    result = scenario_results,
    settings = list(
      S = S,
      h2_values = h2_values,
      info_fractions =
        info_fractions,
      n_select = n_select,
      sigma2_P = sigma2_P,
      target = target,
      M = M,
      B = B,
      ngrid = ngrid,
      seed_info = seed_info,
      seed_outer = seed_outer
    )
  )

  cat(
    "\n====================================================\n"
  )
  cat(
    "Deep validation completed.\n"
  )
  cat(
    "Use the table functions below for compact output.\n"
  )
  cat(
    "====================================================\n"
  )

  invisible(out)
}


# ------------------------------------------------------------
# 7. Long-format result tables
# ------------------------------------------------------------

make_deep_method_table <- function(
    x) {

  rows <- list()
  k <- 1L

  for (z in x$result) {

    tab <- z$methods

    for (m in rownames(tab)) {

      rows[[k]] <- data.frame(
        h2 =
          z$h2,

        information =
          z$information,

        n_pheno =
          z$n_pheno,

        method =
          m,

        valid =
          tab[
            m,
            "valid"
          ],

        valid_rate =
          tab[
            m,
            "valid_rate"
          ],

        bias =
          tab[
            m,
            "bias"
          ],

        RMSE =
          tab[
            m,
            "RMSE"
          ],

        mean_SE =
          tab[
            m,
            "mean_SE"
          ],

        MSE_ratio =
          tab[
            m,
            "MSE_ratio"
          ],

        inclusion =
          tab[
            m,
            "inclusion"
          ],

        width =
          tab[
            m,
            "width"
          ],

        Brier =
          tab[
            m,
            "Brier"
          ],

        mean_p =
          tab[
            m,
            "mean_p"
          ],

        event_rate =
          tab[
            m,
            "event_rate"
          ],

        row.names = NULL
      )

      k <- k + 1L
    }
  }

  do.call(
    rbind,
    rows
  )
}


make_deep_selection_table <- function(
    x) {

  do.call(
    rbind,
    lapply(
      x$result,
      function(z) {

        data.frame(
          h2 =
            z$h2,

          information =
            z$information,

          n_pheno =
            z$n_pheno,

          t(
            z$selection
          ),

          row.names = NULL
        )
      }
    )
  )
}


make_deep_diagnostic_table <- function(
    x) {

  do.call(
    rbind,
    lapply(
      x$result,
      function(z) {

        data.frame(
          h2 =
            z$h2,

          information =
            z$information,

          n_pheno =
            z$n_pheno,

          t(
            z$diagnostics
          ),

          row.names = NULL
        )
      }
    )
  )
}


# ------------------------------------------------------------
# 8. Focused tables for the most important comparisons
# ------------------------------------------------------------

make_deep_coverage_table <- function(
    x) {

  d <- make_deep_method_table(
    x
  )

  d[
    ,
    c(
      "h2",
      "information",
      "n_pheno",
      "method",
      "valid_rate",
      "bias",
      "RMSE",
      "MSE_ratio",
      "inclusion",
      "width"
    )
  ]
}


make_PEC_comparison_table <- function(
    x) {

  s <- make_deep_selection_table(
    x
  )

  s[
    ,
    c(
      "h2",
      "information",
      "n_pheno",
      "mean_PEC_diag_full_ratio"
    )
  ]
}


# ------------------------------------------------------------
# 9. Optional pairwise information-gain table
# ------------------------------------------------------------

make_information_gain_table <- function(
    x) {

  s <- make_deep_selection_table(
    x
  )

  # Sort information using the actual n_pheno.
  s <- s[
    order(
      s$h2,
      s$n_pheno
    ),
    ,
    drop = FALSE
  ]

  rows <- list()
  k <- 1L

  for (h in unique(s$h2)) {

    d <- s[
      s$h2 == h,
      ,
      drop = FALSE
    ]

    if (nrow(d) < 2) {
      next
    }

    for (j in 2:nrow(d)) {

      a <- d[j - 1, ]
      b <- d[j, ]

      rows[[k]] <- data.frame(
        h2 = h,

        from_information =
          a$information,

        to_information =
          b$information,

        delta_n_pheno =
          b$n_pheno -
          a$n_pheno,

        delta_candidate_accuracy =
          b$mean_candidate_accuracy -
          a$mean_candidate_accuracy,

        delta_rank_accuracy =
          b$mean_candidate_rank_accuracy -
          a$mean_candidate_rank_accuracy,

        delta_true_gain =
          b$mean_true_gain -
          a$mean_true_gain,

        delta_oracle_gap =
          b$mean_oracle_gap -
          a$mean_oracle_gap,

        delta_overlap_rate =
          b$mean_overlap_rate -
          a$mean_overlap_rate,

        row.names = NULL
      )

      k <- k + 1L
    }
  }

  do.call(
    rbind,
    rows
  )
}


# ------------------------------------------------------------
# 10. RAM-N failure reasons
# ------------------------------------------------------------

summarize_RAMN_failure_reasons <- function(
    x) {

  if (is.null(
    x$result[[1]]$raw
  )) {
    stop(
      paste0(
        "Raw results were not retained. ",
        "Re-run with keep_raw = TRUE."
      )
    )
  }

  rows <- list()
  k <- 1L

  for (z in x$result) {

    reason_vec <- unlist(
      lapply(
        z$raw,
        `[[`,
        "RAM_N_reasons"
      )
    )

    if (length(reason_vec) == 0) {

      rows[[k]] <- data.frame(
        h2 = z$h2,
        information =
          z$information,
        reason =
          "none",
        count = 0,
        row.names = NULL
      )

      k <- k + 1L

    } else {

      tt <- table(
        reason_vec
      )

      for (r in names(tt)) {

        rows[[k]] <- data.frame(
          h2 = z$h2,
          information =
            z$information,
          reason = r,
          count =
            as.integer(
              tt[r]
            ),
          row.names = NULL
        )

        k <- k + 1L
      }
    }
  }

  do.call(
    rbind,
    rows
  )
}


# ------------------------------------------------------------
# 11. Small self-test
# ------------------------------------------------------------

self_test_realised_policy_deep <- function() {

  z <- run_realised_policy_deep_validation(
    S = 2,
    h2_values = c(
      0.20
    ),
    info_fractions = c(
      low = 0.50,
      high = 1.00
    ),
    n_select = 5,
    M = 100,
    B = 20,
    ngrid = 21,
    seed_info = 123,
    seed_outer = 456,
    keep_raw = TRUE
  )

  print(
    make_deep_coverage_table(
      z
    )
  )

  print(
    make_deep_selection_table(
      z
    )
  )

  cat(
    "\nRealised-policy deep self-test completed.\n"
  )

  invisible(TRUE)
}


# ============================================================
# HOW TO RUN
# ============================================================
#
# STEP 1: source the animal-model base
#
#
# Choose:
#   Simulation_II_pedigree_animal_model_v2.R
#
#
# STEP 2: source the no-RAM-RL validation helpers
#
#
# Choose:
#   Simulation_II_validation_3h2_no_RAMRL.R
#
#
# STEP 3: source THIS file
#
#
# Choose:
#   Simulation_II_realised_policy_deep_3x3.R
#
#
# STEP 4: optional quick test
#
#   self_test_realised_policy_deep()
#
#
# STEP 5: 3 h2 x 3 information levels x 20 replicates
#
#   deep20 <- run_realised_policy_deep_validation(
#     S = 20,
#     keep_raw = TRUE
#   )
#
#
# STEP 6: compact outputs
#
#   coverage_table <- make_deep_coverage_table(
#     deep20
#   )
#
#   selection_table <- make_deep_selection_table(
#     deep20
#   )
#
#   diagnostic_table <- make_deep_diagnostic_table(
#     deep20
#   )
#
#   pec_table <- make_PEC_comparison_table(
#     deep20
#   )
#
#   info_gain <- make_information_gain_table(
#     deep20
#   )
#
#   ramn_fail <- summarize_RAMN_failure_reasons(
#     deep20
#   )
#
#
#   print(coverage_table)
#   print(selection_table)
#   print(diagnostic_table)
#   print(pec_table)
#   print(info_gain)
#   print(ramn_fail)
#
#
# ============================================================
# INTERPRETATION
# ============================================================
#
# Main scientific questions:
#
# Q1. Does uncertainty calibration deteriorate when h2 is low?
#
# Q2. Does increasing phenotype information improve:
#       candidate EBV accuracy,
#       rank accuracy,
#       true realised gain,
#       top-5 overlap with the oracle,
#       oracle regret?
#
# Q3. Does increasing information improve the coverage of the
#     realised-policy prediction interval?
#
# Q4. Under what h2 x information conditions does RAM-N fail
#     its local-normal diagnostic?
#
# Q5. How much does ignoring PEC change uncertainty:
#
#       V_diag / V_full ?
#
# Q6. Is there selection-induced optimism:
#
#       E[ predicted gain - true realised gain ] ?
#
#
# IMPORTANT:
# The selected top-5 set is DATA-DEPENDENT in each outer replicate,
# but once selected, it is held fixed inside RAM-N/bootstrap/Bayesian.
#
# Thus this simulation evaluates uncertainty for the REALISED
# EBV-selected policy, not the population-level performance of the
# selection rule itself.
#
# ============================================================


# ======================================================================
# FINAL BAYESIAN METHOD FOR SIMULATION II
# h2 ~ Beta(1,1), log(sigma_P^2) ~ N(0,1)
#
# This definition intentionally overrides the earlier validation
# bayesian_animal_validation() function below.
# ======================================================================

final_logsumexp2 <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

final_h2_from_eta <- function(eta_A, eta_e) {
  xi <- eta_A - eta_e
  ans <- numeric(length(xi))
  pos <- xi >= 0
  ans[pos] <- 1 / (1 + exp(-xi[pos]))
  ex <- exp(xi[!pos])
  ans[!pos] <- ex / (1 + ex)
  pmin(pmax(ans, 1e-12), 1 - 1e-12)
}

final_log_prior_beta11_eta <- function(
    eta_A, eta_e,
    logP_mean = 0,
    logP_sd = 1) {

  h2 <- final_h2_from_eta(eta_A, eta_e)
  logP <- final_logsumexp2(eta_A, eta_e)

  dbeta(h2, 1, 1, log = TRUE) +
    dnorm(logP, logP_mean, logP_sd, log = TRUE) +
    log(h2) + log1p(-h2)
}

build_final_bayes_grid_animal <- function(
    y, prep, fit = NULL,
    ncoarse = 31,
    ngrid = 41,
    log_drop = 14,
    logP_mean = 0,
    logP_sd = 1) {

  if (is.null(fit)) fit <- fit_reml_animal(y, prep)

  g1c <- seq(fit$lower[1], fit$upper[1], length.out = ncoarse)
  g2c <- seq(fit$lower[2], fit$upper[2], length.out = ncoarse)

  coarse <- expand.grid(eta_A = g1c, eta_e = g2c)

  coarse$logL <- vapply(
    seq_len(nrow(coarse)),
    function(i) {
      log_reml_eta(
        as.numeric(coarse[i, c("eta_A", "eta_e")]),
        y, prep
      )
    },
    numeric(1)
  )

  coarse$logPrior <- final_log_prior_beta11_eta(
    coarse$eta_A, coarse$eta_e,
    logP_mean = logP_mean,
    logP_sd = logP_sd
  )

  coarse$logPost <- coarse$logL + coarse$logPrior
  mx <- max(coarse$logPost)
  keep <- coarse$logPost >= (mx - log_drop)

  if (!any(keep)) {
    stop("Final Bayesian coarse grid retained no points.")
  }

  d1 <- diff(g1c)[1]
  d2 <- diff(g2c)[1]

  lo1 <- max(fit$lower[1], min(coarse$eta_A[keep]) - 2*d1)
  hi1 <- min(fit$upper[1], max(coarse$eta_A[keep]) + 2*d1)
  lo2 <- max(fit$lower[2], min(coarse$eta_e[keep]) - 2*d2)
  hi2 <- min(fit$upper[2], max(coarse$eta_e[keep]) + 2*d2)

  g1 <- seq(lo1, hi1, length.out = ngrid)
  g2 <- seq(lo2, hi2, length.out = ngrid)

  grid <- expand.grid(eta_A = g1, eta_e = g2)

  grid$logL <- vapply(
    seq_len(nrow(grid)),
    function(i) {
      log_reml_eta(
        as.numeric(grid[i, c("eta_A", "eta_e")]),
        y, prep
      )
    },
    numeric(1)
  )

  grid$logPrior <- final_log_prior_beta11_eta(
    grid$eta_A, grid$eta_e,
    logP_mean = logP_mean,
    logP_sd = logP_sd
  )

  grid$logPost <- grid$logL + grid$logPrior
  w <- exp(grid$logPost - max(grid$logPost))
  w <- w / sum(w)

  grid$h2 <- final_h2_from_eta(grid$eta_A, grid$eta_e)

  i1 <- rep(seq_along(g1), times = length(g2))
  i2 <- rep(seq_along(g2), each = length(g1))
  edge <- (
    i1 == 1 | i1 == length(g1) |
    i2 == 1 | i2 == length(g2)
  )

  list(
    grid = grid,
    weight = w,
    fit = fit,
    edge_mass = sum(w[edge])
  )
}

bayesian_animal_validation <- function(
    y, prep, cstar,
    M = 2000,
    ngrid = 41,
    target = 0,
    seed_inner = NULL,
    logP_mean = 0,
    logP_sd = 1,
    ...) {

  if (!is.null(seed_inner)) set.seed(seed_inner)

  fit <- fit_reml_animal(y, prep)

  bg <- build_final_bayes_grid_animal(
    y = y,
    prep = prep,
    fit = fit,
    ngrid = ngrid,
    logP_mean = logP_mean,
    logP_sd = logP_sd
  )

  id <- sample(
    seq_len(nrow(bg$grid)),
    size = M,
    replace = TRUE,
    prob = bg$weight
  )

  uid <- unique(id)
  cache_mu <- numeric(length(uid))
  cache_v <- numeric(length(uid))
  names(cache_mu) <- as.character(uid)
  names(cache_v) <- as.character(uid)

  for (j in seq_along(uid)) {
    ii <- uid[j]
    eta <- as.numeric(bg$grid[ii, c("eta_A", "eta_e")])
    gs <- policy_stats_animal(eta, y, prep, cstar)
    cache_mu[j] <- gs$mu
    cache_v[j] <- gs$v
  }

  mu <- cache_mu[as.character(id)]
  v <- cache_v[as.character(id)]

  delta <- mu + sqrt(pmax(v, 0)) * rnorm(M)
  h2_draw <- bg$grid$h2[id]

  list(
    valid = TRUE,
    mean = mean(delta),
    sd = sd(delta),
    interval = safe_quantile(delta, c(0.025, 0.975)),
    p = mean(delta > target),
    delta = delta,
    mean_h2_post = mean(h2_draw),
    median_h2_post = median(h2_draw),
    edge_mass = bg$edge_mass,
    fit = fit
  )
}



# ======================================================================
# SIMULATION I: BALANCED SIRE MODEL, FIXED POLICY
# ======================================================================

SIMI_ETA_LOWER <- c(-12, -12)
SIMI_ETA_UPPER <- c(10, 10)

simI_simulate <- function(
    m = 24, k = 8,
    sigma2_A = 30, sigma2_e = 70,
    mu = 100,
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  u <- rnorm(m, 0, sqrt(sigma2_A))
  e <- matrix(
    rnorm(m*k, 0, sqrt(sigma2_e)),
    nrow = m, ncol = k
  )

  list(
    y = mu + 0.5*u + e,
    u = u
  )
}

simI_reml_parts <- function(y) {
  m <- nrow(y)
  k <- ncol(y)
  ym <- rowMeans(y)
  ybar <- mean(ym)
  SSW <- sum(sweep(y, 1, ym, "-")^2)
  SSB <- sum((ym - ybar)^2)

  list(
    m = m, k = k,
    ym = ym, ybar = ybar,
    SSW = SSW, SSB = SSB
  )
}

simI_log_reml <- function(eta, y) {
  z <- simI_reml_parts(y)
  sigma2_A <- exp(eta[1])
  sigma2_e <- exp(eta[2])

  q <- 0.25*sigma2_A + sigma2_e/z$k
  dfW <- z$m*(z$k - 1)
  dfB <- z$m - 1

  -0.5*(
    dfW*log(sigma2_e) + z$SSW/sigma2_e +
    dfB*log(q) + z$SSB/q
  )
}

simI_fit_reml <- function(y) {
  z <- simI_reml_parts(y)

  sigma2_e0 <- max(
    z$SSW/(z$m*(z$k - 1)),
    exp(SIMI_ETA_LOWER[2])
  )

  q0 <- z$SSB/(z$m - 1)

  sigma2_A0 <- max(
    4*(q0 - sigma2_e0/z$k),
    exp(SIMI_ETA_LOWER[1])
  )

  f <- function(eta) -simI_log_reml(eta, y)

  opt <- optim(
    log(c(sigma2_A0, sigma2_e0)),
    f,
    method = "L-BFGS-B",
    lower = SIMI_ETA_LOWER,
    upper = SIMI_ETA_UPPER
  )

  H <- optimHess(opt$par, f)
  Hs <- (H + t(H))/2
  ee <- eigen(Hs, symmetric = TRUE)
  lam <- pmax(ee$values, 1e-8)

  Veta <- (
    ee$vectors %*%
    diag(1/lam) %*%
    t(ee$vectors)
  )

  list(
    eta = opt$par,
    theta = exp(opt$par),
    Veta = Veta,
    H = Hs,
    raw_eigen = ee$values,
    convergence = opt$convergence,
    lower = SIMI_ETA_LOWER,
    upper = SIMI_ETA_UPPER
  )
}

simI_fit_reml_fast <- function(y) {
  z <- simI_reml_parts(y)
  dfW <- z$m*(z$k - 1)
  dfB <- z$m - 1

  sigma2_e <- z$SSW/dfW
  q <- z$SSB/dfB
  sigma2_A <- 4*(q - sigma2_e/z$k)

  min_A <- exp(SIMI_ETA_LOWER[1])
  eta_e <- log(sigma2_e)

  if (
    is.finite(sigma2_A) &&
    sigma2_A > min_A &&
    eta_e > SIMI_ETA_LOWER[2] &&
    eta_e < SIMI_ETA_UPPER[2]
  ) {
    return(log(c(sigma2_A, sigma2_e)))
  }

  oe <- optimize(
    function(le) {
      -simI_log_reml(
        c(SIMI_ETA_LOWER[1], le),
        y
      )
    },
    interval = c(
      SIMI_ETA_LOWER[2],
      SIMI_ETA_UPPER[2]
    )
  )

  c(SIMI_ETA_LOWER[1], oe$minimum)
}

simI_policy_c <- function(
    m,
    selected = 1:3,
    gamma = 0.5) {

  w <- rep(0, m)
  w[selected] <- 1/length(selected)
  b <- rep(1/m, m)

  gamma*(w - b)
}

simI_gain_stats <- function(
    eta, y, cstar) {

  z <- simI_reml_parts(y)

  sigma2_A <- exp(eta[1])
  sigma2_e <- exp(eta[2])

  q <- 0.25*sigma2_A + sigma2_e/z$k
  Bcoef <- 0.5*sigma2_A/q

  uhat <- Bcoef*(z$ym - z$ybar)
  mu <- sum(cstar*uhat)

  PEV_i <- sigma2_A - (0.5*sigma2_A)^2/q

  # sum(cstar)=0, so the common intercept-related covariance component
  # cancels for this fixed contrast.
  v <- PEV_i*sum(cstar^2)

  list(
    mu = mu,
    v = max(v, 0),
    uhat = uhat
  )
}

simI_rmvn_rows <- function(n, mu, Sigma) {
  R <- chol((Sigma + t(Sigma))/2)
  Z <- matrix(
    rnorm(n*length(mu)),
    nrow = n,
    ncol = length(mu)
  )
  sweep(Z %*% R, 2, mu, "+")
}

simI_RAM_N <- function(
    y, cstar,
    M = 2000,
    target = 0,
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  fit <- simI_fit_reml(y)
  se <- sqrt(diag(fit$Veta))
  reasons <- character()

  if (
    fit$convergence != 0 ||
    any(!is.finite(fit$eta))
  ) {
    reasons <- c(reasons, "REML_fit")
  }

  if (
    any(!is.finite(fit$raw_eigen)) ||
    min(fit$raw_eigen) <= 1e-6
  ) {
    reasons <- c(reasons, "Hessian")
  }

  if (
    fit$eta[1] <= SIMI_ETA_LOWER[1] + 0.5
  ) {
    reasons <- c(reasons, "variance_boundary")
  }

  if (
    any(!is.finite(se)) ||
    any(se > 4)
  ) {
    reasons <- c(reasons, "large_eta_SE")
  }

  if (length(reasons) > 0) {
    return(
      list(
        valid = FALSE,
        reasons = unique(reasons),
        fit = fit
      )
    )
  }

  eta_draw <- simI_rmvn_rows(
    M, fit$eta, fit$Veta
  )

  mu <- numeric(M)
  v <- numeric(M)

  for (i in seq_len(M)) {
    g <- simI_gain_stats(
      eta_draw[i, ],
      y,
      cstar
    )
    mu[i] <- g$mu
    v[i] <- g$v
  }

  delta <- mu + sqrt(pmax(v, 0))*rnorm(M)

  list(
    valid = TRUE,
    mean = mean(delta),
    sd = sd(delta),
    interval = safe_quantile(
      delta,
      c(0.025, 0.975)
    ),
    p = mean(delta > target),
    delta = delta,
    fit = fit,
    reasons = character()
  )
}

simI_parametric_bootstrap <- function(
    y, cstar,
    B = 500,
    target = 0,
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  fit <- simI_fit_reml(y)
  sigma2_A <- fit$theta[1]
  sigma2_e <- fit$theta[2]
  mu0 <- mean(y)

  current <- simI_gain_stats(
    fit$eta, y, cstar
  )$mu

  err <- numeric(B)

  for (b in seq_len(B)) {
    d <- simI_simulate(
      m = nrow(y),
      k = ncol(y),
      sigma2_A = sigma2_A,
      sigma2_e = sigma2_e,
      mu = mu0
    )

    eta_b <- simI_fit_reml_fast(d$y)

    dhat_b <- simI_gain_stats(
      eta_b, d$y, cstar
    )$mu

    truth_b <- sum(cstar*d$u)
    err[b] <- dhat_b - truth_b
  }

  q <- safe_quantile(
    err,
    c(0.025, 0.975)
  )

  pseudo_delta <- current - err

  list(
    valid = TRUE,
    mean = current - mean(err),
    sd = sd(err),
    interval = c(
      current - q[2],
      current - q[1]
    ),
    p = mean(pseudo_delta > target),
    error = err,
    fit = fit
  )
}

simI_log_prior_beta11 <- function(
    eta_A, eta_e,
    total_variance_center = 100,
    log_total_sd = 1) {

  ratio_A <- final_h2_from_eta(
    eta_A, eta_e
  )

  log_total <- final_logsumexp2(
    eta_A, eta_e
  )

  # The prior is on the variance ratio sigmaA2/(sigmaA2+sigmae2).
  # Because the observation model contains 0.5*u_i, this is not called
  # the record-scale heritability here.
  dbeta(ratio_A, 1, 1, log = TRUE) +
    dnorm(
      log_total,
      mean = log(total_variance_center),
      sd = log_total_sd,
      log = TRUE
    ) +
    log(ratio_A) +
    log1p(-ratio_A)
}

simI_bayes_beta11 <- function(
    y, cstar,
    M = 2000,
    ngrid = 41,
    target = 0,
    total_variance_center = 100,
    log_total_sd = 1,
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  fit <- simI_fit_reml(y)

  ncoarse <- 31
  log_drop <- 14

  g1c <- seq(
    SIMI_ETA_LOWER[1],
    SIMI_ETA_UPPER[1],
    length.out = ncoarse
  )

  g2c <- seq(
    SIMI_ETA_LOWER[2],
    SIMI_ETA_UPPER[2],
    length.out = ncoarse
  )

  coarse <- expand.grid(
    eta_A = g1c,
    eta_e = g2c
  )

  coarse$logL <- vapply(
    seq_len(nrow(coarse)),
    function(i) {
      simI_log_reml(
        as.numeric(
          coarse[i, c("eta_A", "eta_e")]
        ),
        y
      )
    },
    numeric(1)
  )

  coarse$logPrior <- simI_log_prior_beta11(
    coarse$eta_A,
    coarse$eta_e,
    total_variance_center = total_variance_center,
    log_total_sd = log_total_sd
  )

  coarse$logPost <- coarse$logL + coarse$logPrior

  mx <- max(coarse$logPost)
  keep <- coarse$logPost >= (mx - log_drop)

  if (!any(keep)) {
    stop("Simulation I Bayesian coarse grid retained no points.")
  }

  d1 <- diff(g1c)[1]
  d2 <- diff(g2c)[1]

  lo1 <- max(
    SIMI_ETA_LOWER[1],
    min(coarse$eta_A[keep]) - 2*d1
  )
  hi1 <- min(
    SIMI_ETA_UPPER[1],
    max(coarse$eta_A[keep]) + 2*d1
  )
  lo2 <- max(
    SIMI_ETA_LOWER[2],
    min(coarse$eta_e[keep]) - 2*d2
  )
  hi2 <- min(
    SIMI_ETA_UPPER[2],
    max(coarse$eta_e[keep]) + 2*d2
  )

  g1 <- seq(lo1, hi1, length.out = ngrid)
  g2 <- seq(lo2, hi2, length.out = ngrid)

  grid <- expand.grid(
    eta_A = g1,
    eta_e = g2
  )

  grid$logL <- vapply(
    seq_len(nrow(grid)),
    function(i) {
      simI_log_reml(
        as.numeric(
          grid[i, c("eta_A", "eta_e")]
        ),
        y
      )
    },
    numeric(1)
  )

  grid$logPrior <- simI_log_prior_beta11(
    grid$eta_A,
    grid$eta_e,
    total_variance_center = total_variance_center,
    log_total_sd = log_total_sd
  )

  logPost <- grid$logL + grid$logPrior
  w <- exp(logPost - max(logPost))
  w <- w/sum(w)

  id <- sample(
    seq_len(nrow(grid)),
    M,
    replace = TRUE,
    prob = w
  )

  uid <- unique(id)
  cache_mu <- numeric(length(uid))
  cache_v <- numeric(length(uid))
  names(cache_mu) <- as.character(uid)
  names(cache_v) <- as.character(uid)

  for (j in seq_along(uid)) {
    ii <- uid[j]
    eta <- as.numeric(
      grid[ii, c("eta_A", "eta_e")]
    )
    g <- simI_gain_stats(
      eta, y, cstar
    )
    cache_mu[j] <- g$mu
    cache_v[j] <- g$v
  }

  mu <- cache_mu[as.character(id)]
  v <- cache_v[as.character(id)]
  delta <- mu + sqrt(pmax(v, 0))*rnorm(M)

  ratio_draw <- final_h2_from_eta(
    grid$eta_A[id],
    grid$eta_e[id]
  )

  list(
    valid = TRUE,
    mean = mean(delta),
    sd = sd(delta),
    interval = safe_quantile(
      delta,
      c(0.025, 0.975)
    ),
    p = mean(delta > target),
    delta = delta,
    mean_variance_ratio_A = mean(ratio_draw),
    fit = fit
  )
}

simI_pack <- function(x) {
  if (
    is.null(x$valid) ||
    !isTRUE(x$valid)
  ) {
    return(NULL)
  }

  c(
    mean = x$mean,
    sd = x$sd,
    lo = unname(x$interval[1]),
    hi = unname(x$interval[2]),
    p = x$p
  )
}

simI_one_outer <- function(
    m, k,
    sigma2_A, sigma2_e,
    selected = 1:3,
    gamma = 0.5,
    target = 0,
    M = 2000,
    B = 500,
    ngrid = 41,
    seed_data,
    seed_RAMN,
    seed_boot,
    seed_bayes) {

  d <- simI_simulate(
    m = m,
    k = k,
    sigma2_A = sigma2_A,
    sigma2_e = sigma2_e,
    seed = seed_data
  )

  cstar <- simI_policy_c(
    m,
    selected = selected,
    gamma = gamma
  )

  truth <- sum(cstar*d$u)

  fit <- simI_fit_reml(d$y)
  g <- simI_gain_stats(
    fit$eta, d$y, cstar
  )

  sd0 <- sqrt(g$v)

  conditional <- c(
    mean = g$mu,
    sd = sd0,
    lo = g$mu - 1.96*sd0,
    hi = g$mu + 1.96*sd0,
    p = 1 - pnorm(
      target,
      mean = g$mu,
      sd = max(sd0, 1e-12)
    )
  )

  rn <- simI_RAM_N(
    d$y, cstar,
    M = M,
    target = target,
    seed = seed_RAMN
  )

  pb <- simI_parametric_bootstrap(
    d$y, cstar,
    B = B,
    target = target,
    seed = seed_boot
  )

  by <- simI_bayes_beta11(
    d$y, cstar,
    M = M,
    ngrid = ngrid,
    target = target,
    total_variance_center = sigma2_A + sigma2_e,
    seed = seed_bayes
  )

  ratio_hat <- fit$theta[1]/sum(fit$theta)

  list(
    truth = truth,
    conditional = conditional,
    RAM_N = simI_pack(rn),
    bootstrap = simI_pack(pb),
    Bayesian_Beta11 = simI_pack(by),
    RAM_N_reasons = if (isTRUE(rn$valid)) character() else rn$reasons,
    reml_diag = c(
      variance_ratio_hat = ratio_hat,
      sigma2A_hat = fit$theta[1],
      sigma2e_hat = fit$theta[2],
      hessian_min = min(fit$raw_eigen),
      etaA_SE = sqrt(fit$Veta[1, 1]),
      near_boundary = as.numeric(
        fit$eta[1] <= SIMI_ETA_LOWER[1] + 0.5
      )
    ),
    bayes_ratio_mean = by$mean_variance_ratio_A
  )
}

simI_summarize_method <- function(
    res, name, target = 0) {

  x <- lapply(res, `[[`, name)

  ok <- !vapply(
    x,
    is.null,
    logical(1)
  )

  x <- x[ok]

  if (length(x) == 0) {
    return(
      c(
        valid = 0,
        bias = NA,
        RMSE = NA,
        mean_SE = NA,
        MSE_ratio = NA,
        inclusion = NA,
        width = NA,
        Brier = NA,
        mean_p = NA,
        event_rate = NA
      )
    )
  }

  truth <- vapply(
    res[ok],
    `[[`,
    numeric(1),
    "truth"
  )

  mat <- do.call(rbind, x)
  err <- mat[, "mean"] - truth
  mc_mse <- mean(err^2)
  mean_est_mse <- mean(mat[, "sd"]^2)
  z <- as.numeric(truth > target)

  c(
    valid = length(x),
    bias = mean(err),
    RMSE = sqrt(mc_mse),
    mean_SE = mean(mat[, "sd"]),
    MSE_ratio = mean_est_mse/mc_mse,
    inclusion = mean(
      truth >= mat[, "lo"] &
      truth <= mat[, "hi"]
    ),
    width = mean(
      mat[, "hi"] - mat[, "lo"]
    ),
    Brier = mean((mat[, "p"] - z)^2),
    mean_p = mean(mat[, "p"]),
    event_rate = mean(z)
  )
}

simI_summarize_scenario <- function(
    res, target = 0) {

  methods <- c(
    "conditional",
    "RAM_N",
    "bootstrap",
    "Bayesian_Beta11"
  )

  tab <- t(
    vapply(
      methods,
      function(m) {
        simI_summarize_method(
          res, m, target
        )
      },
      numeric(10)
    )
  )

  rownames(tab) <- methods

  cbind(
    tab,
    valid_rate = tab[, "valid"]/length(res)
  )
}

simI_summarize_diagnostics <- function(res) {
  D <- do.call(
    rbind,
    lapply(
      res,
      `[[`,
      "reml_diag"
    )
  )

  ramn_valid <- !vapply(
    lapply(
      res,
      `[[`,
      "RAM_N"
    ),
    is.null,
    logical(1)
  )

  c(
    mean_variance_ratio_hat = mean(D[, "variance_ratio_hat"]),
    median_variance_ratio_hat = median(D[, "variance_ratio_hat"]),
    mean_sigma2A_hat = mean(D[, "sigma2A_hat"]),
    mean_sigma2e_hat = mean(D[, "sigma2e_hat"]),
    mean_hessian_min = mean(D[, "hessian_min"]),
    median_etaA_SE = median(D[, "etaA_SE"], na.rm = TRUE),
    boundary_rate = mean(D[, "near_boundary"]),
    RAM_N_valid_rate = mean(ramn_valid),
    mean_Bayes_variance_ratio = mean(
      vapply(
        res,
        `[[`,
        numeric(1),
        "bayes_ratio_mean"
      )
    )
  )
}



# ======================================================================
# CHECKPOINT / OUTPUT HELPERS
# ======================================================================

final_dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(
      path,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  invisible(path)
}

final_safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

final_add_coverage_mcse <- function(tab) {
  if (
    !"inclusion" %in% names(tab) ||
    !"valid" %in% names(tab)
  ) {
    return(tab)
  }

  tab$coverage_MCSE <- sqrt(
    pmax(
      tab$inclusion*(1 - tab$inclusion)/
        pmax(tab$valid, 1),
      0
    )
  )

  tab
}

final_write_csv <- function(x, path) {
  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  invisible(path)
}


# ======================================================================
# FINAL RUNNER: SIMULATION I
# ======================================================================

run_simulation_I_final <- function(
    S = 1000,
    M = 2000,
    B = 500,
    ngrid = 41,
    target = 0,
    checkpoint_every = 25,
    output_dir = file.path(
      getwd(),
      "simulation_I_II_final_output",
      "Simulation_I"
    ),
    resume = TRUE) {

  final_dir_create(output_dir)

  scenarios <- list(
    Scenario_A = list(
      m = 24,
      k = 8,
      sigma2_A = 30,
      sigma2_e = 70,
      seed = 20260801
    ),
    Scenario_B = list(
      m = 16,
      k = 5,
      sigma2_A = 10,
      sigma2_e = 90,
      seed = 20260802
    )
  )

  final_summary <- list()

  for (sc_name in names(scenarios)) {
    sc <- scenarios[[sc_name]]

    cat(
      "\n====================================================\n",
      "Simulation I: ", sc_name, "\n",
      "m=", sc$m,
      ", k=", sc$k,
      ", sigma2_A=", sc$sigma2_A,
      ", sigma2_e=", sc$sigma2_e,
      "\n====================================================\n",
      sep = ""
    )

    checkpoint_file <- file.path(
      output_dir,
      paste0(final_safe_name(sc_name), "_checkpoint.rds")
    )

    final_file <- file.path(
      output_dir,
      paste0(final_safe_name(sc_name), "_final.rds")
    )

    if (
      resume &&
      file.exists(final_file)
    ) {
      cat("Loading completed scenario from checkpoint.\n")
      done <- readRDS(final_file)
      final_summary[[sc_name]] <- done$summary
      next
    }

    set.seed(sc$seed)

    seed_matrix <- matrix(
      sample.int(
        .Machine$integer.max,
        S*4
      ),
      nrow = S,
      ncol = 4
    )

    res <- vector("list", S)
    start_at <- 1L

    if (
      resume &&
      file.exists(checkpoint_file)
    ) {
      cp <- readRDS(checkpoint_file)

      if (
        identical(cp$S, S) &&
        identical(cp$scenario, sc_name)
      ) {
        res <- cp$res

        missing_idx <- which(
          vapply(res, is.null, logical(1))
        )

        if (length(missing_idx) == 0) {
          start_at <- S + 1L
        } else {
          start_at <- min(missing_idx)
        }

        cat(
          "Resuming at outer replicate ",
          start_at,
          ".\n",
          sep = ""
        )
      }
    }

    if (start_at <= S) {
      for (s in seq.int(start_at, S)) {

        res[[s]] <- simI_one_outer(
          m = sc$m,
          k = sc$k,
          sigma2_A = sc$sigma2_A,
          sigma2_e = sc$sigma2_e,
          target = target,
          M = M,
          B = B,
          ngrid = ngrid,
          seed_data = seed_matrix[s, 1],
          seed_RAMN = seed_matrix[s, 2],
          seed_boot = seed_matrix[s, 3],
          seed_bayes = seed_matrix[s, 4]
        )

        if (
          s %% checkpoint_every == 0 ||
          s == S
        ) {
          saveRDS(
            list(
              S = S,
              scenario = sc_name,
              settings = sc,
              res = res
            ),
            checkpoint_file
          )

          cat(
            "Completed ",
            s,
            " / ",
            S,
            "\n",
            sep = ""
          )
        }
      }
    }

    method_tab <- simI_summarize_scenario(
      res,
      target = target
    )

    method_df <- data.frame(
      scenario = sc_name,
      method = rownames(method_tab),
      method_tab,
      row.names = NULL,
      check.names = FALSE
    )

    method_df <- final_add_coverage_mcse(method_df)

    diag <- simI_summarize_diagnostics(res)

    diag_df <- data.frame(
      scenario = sc_name,
      t(diag),
      row.names = NULL,
      check.names = FALSE
    )

    one_summary <- list(
      scenario = sc,
      methods = method_df,
      diagnostics = diag_df
    )

    saveRDS(
      list(
        S = S,
        scenario = sc_name,
        settings = sc,
        summary = one_summary,
        raw = res
      ),
      final_file
    )

    final_summary[[sc_name]] <- one_summary
  }

  method_all <- do.call(
    rbind,
    lapply(final_summary, `[[`, "methods")
  )

  diag_all <- do.call(
    rbind,
    lapply(final_summary, `[[`, "diagnostics")
  )

  final_write_csv(
    method_all,
    file.path(
      output_dir,
      "Simulation_I_method_performance.csv"
    )
  )

  final_write_csv(
    diag_all,
    file.path(
      output_dir,
      "Simulation_I_REML_RAMN_diagnostics.csv"
    )
  )

  out <- list(
    settings = list(
      S = S,
      M = M,
      B = B,
      ngrid = ngrid,
      target = target
    ),
    scenarios = final_summary,
    method_table = method_all,
    diagnostic_table = diag_all
  )

  saveRDS(
    out,
    file.path(output_dir, "Simulation_I_ALL.rds")
  )

  cat("\nSimulation I completed.\n")
  invisible(out)
}


# ======================================================================
# FINAL RUNNER: SIMULATION II
# ======================================================================

run_simulation_II_final <- function(
    S = 1000,
    h2_values = c(0.05, 0.20, 0.40),
    info_fractions = c(
      low = 0.40,
      moderate = 0.70,
      high = 1.00
    ),
    n_select = 5,
    sigma2_P = 1,
    beta = 0,
    target = 0,
    M = 2000,
    B = 500,
    ngrid = 41,
    seed_info = 20260850,
    seed_outer = 20260851,
    checkpoint_every = 25,
    output_dir = file.path(
      getwd(),
      "simulation_I_II_final_output",
      "Simulation_II"
    ),
    resume = TRUE,
    keep_raw_in_master = FALSE) {

  final_dir_create(output_dir)

  base_design <- build_simII_design()

  info_bundle <- build_information_level_designs(
    base_design = base_design,
    info_fractions = info_fractions,
    seed_info = seed_info
  )

  cat(
    "\n====================================================\n",
    "Simulation II: realised EBV-selection policy\n",
    "====================================================\n",
    sep = ""
  )

  cat("\nCandidate relationship summary:\n")
  print(
    relationship_summary(
      base_design$A,
      base_design$candidates
    )
  )

  print_information_design(info_bundle)

  set.seed(seed_outer)

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
      S*n_h*n_i*3
    ),
    dim = c(S, n_h, n_i, 3)
  )

  scenario_results <- list()

  for (ih in seq_along(h2_values)) {
    h2 <- h2_values[ih]

    for (ii in seq_along(info_bundle$info_levels)) {
      info <- info_bundle$info_levels[[ii]]

      key <- paste0(
        "h2_",
        h2,
        "__",
        info$name
      )

      cat(
        "\n====================================================\n",
        "Simulation II: ",
        key,
        " | n_pheno=",
        info$n_pheno,
        "\n====================================================\n",
        sep = ""
      )

      checkpoint_file <- file.path(
        output_dir,
        paste0(final_safe_name(key), "_checkpoint.rds")
      )

      final_file <- file.path(
        output_dir,
        paste0(final_safe_name(key), "_final.rds")
      )

      if (
        resume &&
        file.exists(final_file)
      ) {
        cat("Loading completed scenario from checkpoint.\n")
        done <- readRDS(final_file)
        scenario_results[[key]] <- done$summary

        if (
          keep_raw_in_master &&
          !is.null(done$raw)
        ) {
          scenario_results[[key]]$raw <- done$raw
        }

        next
      }

      res <- vector("list", S)
      start_at <- 1L

      if (
        resume &&
        file.exists(checkpoint_file)
      ) {
        cp <- readRDS(checkpoint_file)

        if (
          identical(cp$S, S) &&
          identical(cp$key, key)
        ) {
          res <- cp$res

          missing_idx <- which(
            vapply(res, is.null, logical(1))
          )

          if (length(missing_idx) == 0) {
            start_at <- S + 1L
          } else {
            start_at <- min(missing_idx)
          }

          cat(
            "Resuming at outer replicate ",
            start_at,
            ".\n",
            sep = ""
          )
        }
      }

      if (start_at <= S) {
        for (s in seq.int(start_at, S)) {

          res[[s]] <- one_realised_policy_scenario(
            latent = latent_list[[s]],
            info_bundle = info_bundle,
            info_level = info,
            h2 = h2,
            sigma2_P = sigma2_P,
            beta = beta,
            target = target,
            n_select = n_select,
            M = M,
            B = B,
            ngrid = ngrid,
            seed_RAMN = seed_methods[s, ih, ii, 1],
            seed_boot = seed_methods[s, ih, ii, 2],
            seed_bayes = seed_methods[s, ih, ii, 3]
          )

          if (
            s %% checkpoint_every == 0 ||
            s == S
          ) {
            saveRDS(
              list(
                S = S,
                key = key,
                h2 = h2,
                information = info$name,
                res = res
              ),
              checkpoint_file
            )

            cat(
              "Completed ",
              s,
              " / ",
              S,
              "\n",
              sep = ""
            )
          }
        }
      }

      summary_one <- list(
        h2 = h2,
        information = info$name,
        n_pheno = info$n_pheno,
        methods = summarize_realised_policy_methods(
          res,
          target = target
        ),
        selection = summarize_realised_policy_selection(res),
        diagnostics = summarize_realised_policy_reml(res)
      )

      saveRDS(
        list(
          S = S,
          key = key,
          summary = summary_one,
          raw = res
        ),
        final_file
      )

      scenario_results[[key]] <- summary_one

      if (keep_raw_in_master) {
        scenario_results[[key]]$raw <- res
      }
    }
  }

  out <- list(
    design = base_design,
    information_design = info_bundle,
    result = scenario_results,
    settings = list(
      S = S,
      h2_values = h2_values,
      info_fractions = info_fractions,
      n_select = n_select,
      sigma2_P = sigma2_P,
      beta = beta,
      target = target,
      M = M,
      B = B,
      ngrid = ngrid,
      Bayesian_prior =
        "h2~Beta(1,1); log(sigmaP2)~N(0,1)",
      seed_info = seed_info,
      seed_outer = seed_outer
    )
  )

  method_table <- make_deep_method_table(out)
  method_table <- final_add_coverage_mcse(method_table)

  selection_table <- make_deep_selection_table(out)
  diagnostic_table <- make_deep_diagnostic_table(out)
  pec_table <- make_PEC_comparison_table(out)
  info_gain <- make_information_gain_table(out)

  final_write_csv(
    method_table,
    file.path(
      output_dir,
      "Simulation_II_method_performance.csv"
    )
  )

  final_write_csv(
    selection_table,
    file.path(
      output_dir,
      "Simulation_II_selection_performance.csv"
    )
  )

  final_write_csv(
    diagnostic_table,
    file.path(
      output_dir,
      "Simulation_II_REML_RAMN_diagnostics.csv"
    )
  )

  final_write_csv(
    pec_table,
    file.path(
      output_dir,
      "Simulation_II_PEC_comparison.csv"
    )
  )

  final_write_csv(
    info_gain,
    file.path(
      output_dir,
      "Simulation_II_information_gain.csv"
    )
  )

  out$method_table <- method_table
  out$selection_table <- selection_table
  out$diagnostic_table <- diagnostic_table
  out$pec_table <- pec_table
  out$information_gain_table <- info_gain

  saveRDS(
    out,
    file.path(output_dir, "Simulation_II_ALL.rds")
  )

  cat("\nSimulation II completed.\n")
  invisible(out)
}


# ======================================================================
# RUN BOTH SIMULATIONS FROM ONE COMMAND
# ======================================================================

run_all_final <- function(
    S_I = 1000,
    S_II = 1000,
    M = 2000,
    B = 500,
    ngrid = 41,
    checkpoint_every = 25,
    root_output_dir = file.path(
      getwd(),
      "simulation_I_II_final_output"
    ),
    resume = TRUE) {

  final_dir_create(root_output_dir)

  cat(
    "\n####################################################\n",
    "FINAL SIMULATION I + II\n",
    "####################################################\n",
    sep = ""
  )

  simI <- run_simulation_I_final(
    S = S_I,
    M = M,
    B = B,
    ngrid = ngrid,
    checkpoint_every = checkpoint_every,
    output_dir = file.path(
      root_output_dir,
      "Simulation_I"
    ),
    resume = resume
  )

  simII <- run_simulation_II_final(
    S = S_II,
    M = M,
    B = B,
    ngrid = ngrid,
    checkpoint_every = checkpoint_every,
    output_dir = file.path(
      root_output_dir,
      "Simulation_II"
    ),
    resume = resume
  )

  out <- list(
    Simulation_I = simI,
    Simulation_II = simII,
    settings = list(
      S_I = S_I,
      S_II = S_II,
      M = M,
      B = B,
      ngrid = ngrid,
      checkpoint_every = checkpoint_every,
      root_output_dir = normalizePath(
        root_output_dir,
        winslash = "/",
        mustWork = FALSE
      )
    )
  )

  saveRDS(
    out,
    file.path(
      root_output_dir,
      "Simulation_I_II_FINAL_ALL.rds"
    )
  )

  cat(
    "\n####################################################\n",
    "ALL FINAL SIMULATIONS COMPLETED\n",
    "Output directory:\n",
    normalizePath(
      root_output_dir,
      winslash = "/",
      mustWork = FALSE
    ),
    "\n####################################################\n",
    sep = ""
  )

  invisible(out)
}


# ======================================================================
# SMALL INTEGRATED SELF-TEST
# ======================================================================

self_test_final_simulations <- function() {

  cat("\nRunning integrated self-test...\n")

  tmp <- file.path(
    tempdir(),
    "genetic_gain_final_selftest"
  )

  unlink(
    tmp,
    recursive = TRUE,
    force = TRUE
  )

  i <- run_simulation_I_final(
    S = 2,
    M = 100,
    B = 20,
    ngrid = 21,
    checkpoint_every = 1,
    output_dir = file.path(tmp, "I"),
    resume = FALSE
  )

  ii <- run_simulation_II_final(
    S = 2,
    h2_values = c(0.20),
    info_fractions = c(
      low = 0.50,
      high = 1.00
    ),
    M = 100,
    B = 20,
    ngrid = 21,
    checkpoint_every = 1,
    output_dir = file.path(tmp, "II"),
    resume = FALSE
  )

  stopifnot(
    nrow(i$method_table) == 8,
    all(
      c(
        "conditional",
        "RAM_N",
        "bootstrap",
        "Bayesian_Beta11"
      ) %in% i$method_table$method
    ),
    all(
      c(
        "conditional",
        "conditional_diag",
        "RAM_N_policy",
        "bootstrap_policy",
        "Bayesian_policy"
      ) %in% ii$method_table$method
    )
  )

  cat(
    "\nIntegrated self-test completed successfully.\n"
  )

  invisible(
    list(
      Simulation_I = i,
      Simulation_II = ii
    )
  )
}


# ======================================================================
# USAGE
# ======================================================================
#
#
# Optional one-time structural test:
#   self_test_final_simulations()
#
# Run both final simulations:
#   final <- run_all_final()
#
# For 2000 outer replicates:
#   final <- run_all_final(
#     S_I = 2000,
#     S_II = 2000
#   )
#
# Results are written automatically under:
#   ./simulation_I_II_final_output/
#
# Resume is TRUE by default.
# ======================================================================



# ======================================================================
# SIMULATION II EXTENSION
# BOOTSTRAP GENERATING-MODEL DIAGNOSTIC
# ======================================================================
#
# PURPOSE
# -------
# The previous true-VC oracle showed that, with the realised data-selected
# policy c(Y) held fixed, using the TRUE variance components restores
# approximately nominal coverage for the conditional interval.
#
# The next question is therefore:
#
#   Why did the realised-policy parametric bootstrap still undercover?
#
# This extension isolates ONE feature of the bootstrap:
# the variance components used to GENERATE the inner bootstrap world.
#
# In each OUTER replicate:
#   1. generate the same observed dataset as the final Simulation II,
#   2. REML -> EBLUP -> select top 5 -> define c(Y),
#   3. KEEP THE SAME c(Y) FIXED in all inner calculations,
#   4. compare two bootstrap generators:
#
#      (A) plug-in generator:
#          y*, G* are generated with theta_hat from the outer dataset
#
#      (B) true-VC generator (oracle diagnostic):
#          y*, G* are generated with theta_true
#
#   5. in BOTH bootstrap worlds, each inner replicate re-estimates variance
#      components by REML and computes EBLUP using the same fixed c(Y).
#
# IMPORTANT
# ---------
# - There is NO reselection in the inner bootstrap.
# - The OUTER estimand is unchanged: G = c(Y)'u.
# - The estimator applied to the observed data is unchanged: REML-EBLUP.
# - The basic/error bootstrap interval is unchanged.
# - For the cleanest one-factor comparison, the same fitted fixed-effect
#   mean X*beta_hat from the outer plug-in fit is used in BOTH generators.
#   Thus the main intended change is the covariance/variance-component
#   generator only.
# - The SAME standard-normal inner random numbers are used for plug-in and
#   true-VC generators. This makes the generator comparison paired and
#   reduces inner Monte Carlo noise.
#
# MAIN SCIENTIFIC INTERPRETATION
# ------------------------------
# If true-VC-generator bootstrap coverage is close to 0.95 while plug-in
# bootstrap remains low, the remaining bootstrap undercoverage is largely
# attributable to the plug-in generating model for variance components.
#
# If true-VC-generator bootstrap still undercovers materially, then the
# next diagnostic should focus on interval construction / finite-sample
# bootstrap behavior rather than attributing the deficit to selection.
#
# DEFAULT DESIGN
# --------------
# Same as the previous production Simulation II:
#   h2 = 0.05, 0.20, 0.40
#   information = low, moderate, high (~60, 90, 120 phenotypes)
#   20 candidate males; top 5 selected by REML-EBLUP
#   outer replicates = 1000
#   bootstrap replicates = 500
#   seed_info  = 20260850
#   seed_outer = 20260851
#
# ONE-FILE USE
# ------------
# This file is intended to be appended to the complete production file.
# After sourcing the combined file, run e.g.
#
#   boot20 <- run_bootstrap_generator_diagnostic_final(
#     S = 20,
#     B = 100,
#     checkpoint_every = 5,
#     resume = FALSE
#   )
#
# Then the production run:
#
#   boot1000 <- run_bootstrap_generator_diagnostic_final(
#     S = 1000,
#     B = 500,
#     checkpoint_every = 25,
#     resume = TRUE
#   )
#
# Main tables:
#   boot1000$method_table
#   boot1000$comparison_table
#   boot1000$diagnostic_table
#
# ======================================================================


# ----------------------------------------------------------------------
# 1. Fixed-eta conditional interval, used only as an oracle reference
# ----------------------------------------------------------------------

bootdiag_fixed_eta_interval <- function(
    eta,
    y,
    prep,
    cstar,
    target = 0) {

  g <- policy_stats_animal(
    eta = eta,
    y = y,
    prep = prep,
    cstar = cstar,
    calc_diag = FALSE
  )

  sd_use <- max(g$sd, 1e-12)

  list(
    valid = TRUE,
    mean = g$mu,
    sd = g$sd,
    interval = g$mu + c(-1, 1) * 1.96 * g$sd,
    p = 1 - pnorm(
      target,
      mean = g$mu,
      sd = sd_use
    )
  )
}


# ----------------------------------------------------------------------
# 2. GLS beta under supplied covariance
# ----------------------------------------------------------------------

bootdiag_gls_beta <- function(
    y,
    prep,
    sigma2_A,
    sigma2_e) {

  X <- prep$X

  ViX <- Vinv_apply(
    prep,
    sigma2_A,
    sigma2_e,
    X
  )

  Viy <- Vinv_apply(
    prep,
    sigma2_A,
    sigma2_e,
    y
  )

  drop(
    solve(
      crossprod(X, ViX),
      crossprod(X, Viy)
    )
  )
}


# ----------------------------------------------------------------------
# 3. Build joint covariance of (y*, G*) for a fixed policy
# ----------------------------------------------------------------------

bootdiag_joint_generator <- function(
    prep,
    cstar,
    sigma2_A,
    sigma2_e,
    beta_generator) {

  if (!is.finite(sigma2_A) || sigma2_A <= 0) {
    stop("sigma2_A must be positive and finite.")
  }

  if (!is.finite(sigma2_e) || sigma2_e <= 0) {
    stop("sigma2_e must be positive and finite.")
  }

  X <- prep$X
  A <- prep$A
  ids <- prep$pheno_ids

  beta_generator <- as.numeric(beta_generator)

  if (length(beta_generator) != ncol(X)) {
    stop("beta_generator has incompatible length.")
  }

  cAc <- drop(
    crossprod(
      cstar,
      A %*% cstar
    )
  )

  a_obs_c <- drop(
    A[ids, , drop = FALSE] %*% cstar
  )

  V <- (
    sigma2_A * prep$K +
      sigma2_e * diag(prep$nobs)
  )

  h <- sigma2_A * a_obs_c
  var_delta <- sigma2_A * cAc

  Sjoint <- rbind(
    cbind(V, h),
    c(h, var_delta)
  )

  Sjoint <- (Sjoint + t(Sjoint)) / 2

  # Small numerical stabilizer only.
  diag(Sjoint) <- diag(Sjoint) + 1e-12

  mean_joint <- c(
    drop(X %*% beta_generator),
    0
  )

  list(
    mean = mean_joint,
    Sigma = Sjoint,
    chol = chol(Sjoint)
  )
}


# ----------------------------------------------------------------------
# 4. Bootstrap with a user-specified generating covariance
# ----------------------------------------------------------------------
#
# The observed-data estimator remains REML-EBLUP from the actual y.
# Each inner replicate also re-estimates VC by REML.
# Only the inner DATA GENERATOR can be switched between theta_hat and
# theta_true.
#
# z_joint can be supplied so that competing generators use exactly the
# same standard-normal random numbers.
# ----------------------------------------------------------------------

parametric_bootstrap_fixed_policy_generator <- function(
    y,
    prep,
    cstar,
    sigma2_A_gen,
    sigma2_e_gen,
    B = 500,
    target = 0,
    current_fit = NULL,
    beta_generator = NULL,
    z_joint = NULL,
    return_inner = FALSE) {

  if (is.null(current_fit)) {
    current_fit <- fit_reml_animal(
      y,
      prep
    )
  }

  gs0 <- policy_stats_animal(
    current_fit$eta,
    y,
    prep,
    cstar
  )

  current <- gs0$mu

  # Default generator mean uses the observed outer plug-in GLS beta.
  # This keeps the fixed-effect generator identical between plug-in and
  # true-VC covariance generators in the main diagnostic.
  if (is.null(beta_generator)) {
    beta_generator <- bootdiag_gls_beta(
      y = y,
      prep = prep,
      sigma2_A = current_fit$theta[1],
      sigma2_e = current_fit$theta[2]
    )
  }

  joint <- bootdiag_joint_generator(
    prep = prep,
    cstar = cstar,
    sigma2_A = sigma2_A_gen,
    sigma2_e = sigma2_e_gen,
    beta_generator = beta_generator
  )

  d <- length(joint$mean)

  if (is.null(z_joint)) {
    z_joint <- matrix(
      rnorm(B * d),
      nrow = B,
      ncol = d
    )
  }

  if (!is.matrix(z_joint) ||
      nrow(z_joint) != B ||
      ncol(z_joint) != d) {
    stop("z_joint must be a B x (nobs+1) matrix.")
  }

  draws <- sweep(
    z_joint %*% joint$chol,
    2,
    joint$mean,
    "+"
  )

  err <- rep(NA_real_, B)
  se_b <- rep(NA_real_, B)
  h2_b <- rep(NA_real_, B)
  conv_b <- rep(NA_real_, B)
  boundary_b <- rep(NA_real_, B)

  for (b in seq_len(B)) {

    yb <- draws[
      b,
      seq_len(prep$nobs)
    ]

    truth_b <- draws[
      b,
      prep$nobs + 1L
    ]

    fit_b <- fit_reml_animal(
      yb,
      prep
    )

    gs_b <- policy_stats_animal(
      fit_b$eta,
      yb,
      prep,
      cstar
    )

    err[b] <- gs_b$mu - truth_b
    se_b[b] <- gs_b$sd

    h2_b[b] <- (
      fit_b$theta[1] /
        sum(fit_b$theta)
    )

    conv_b[b] <- as.numeric(
      fit_b$convergence == 0
    )

    boundary_b[b] <- as.numeric(
      fit_b$eta[1] - fit_b$lower[1] < 0.05
    )
  }

  ok <- is.finite(err)

  if (sum(ok) < max(20, ceiling(0.80 * B))) {
    return(
      list(
        valid = FALSE,
        reason = "too_few_valid_inner_bootstrap_replicates",
        n_valid = sum(ok),
        B = B
      )
    )
  }

  err_ok <- err[ok]
  se_ok <- se_b[ok]

  q <- safe_quantile(
    err_ok,
    c(0.025, 0.975)
  )

  interval <- c(
    current - q[2],
    current - q[1]
  )

  pseudo_delta <- current - err_ok

  out <- list(
    valid = TRUE,
    mean = current - mean(err_ok),
    sd = sd(err_ok),
    interval = interval,
    p = mean(pseudo_delta > target),
    current = current,
    current_se = gs0$sd,
    error_mean = mean(err_ok),
    error_sd = sd(err_ok),
    error_q025 = q[1],
    error_q975 = q[2],
    n_valid = length(err_ok),
    valid_fraction = length(err_ok) / B,
    inner_REML_convergence_rate = mean(conv_b[ok], na.rm = TRUE),
    inner_boundary_rate = mean(boundary_b[ok], na.rm = TRUE),
    inner_mean_h2_hat = mean(h2_b[ok], na.rm = TRUE),
    inner_median_h2_hat = median(h2_b[ok], na.rm = TRUE),
    generator_sigma2_A = sigma2_A_gen,
    generator_sigma2_e = sigma2_e_gen,
    generator_h2 = sigma2_A_gen / (sigma2_A_gen + sigma2_e_gen)
  )

  if (return_inner) {
    out$error <- err
    out$inner_se <- se_b
    out$inner_h2_hat <- h2_b
    out$inner_converged <- conv_b
    out$inner_boundary <- boundary_b
  }

  out
}


# ----------------------------------------------------------------------
# 5. Pack a method to the standard five-number interface
# ----------------------------------------------------------------------

bootdiag_pack_method <- function(x) {

  if (is.null(x) ||
      is.null(x$valid) ||
      !isTRUE(x$valid)) {
    return(NULL)
  }

  c(
    mean = x$mean,
    sd = x$sd,
    lo = unname(x$interval[1]),
    hi = unname(x$interval[2]),
    p = x$p
  )
}


# ----------------------------------------------------------------------
# 6. One OUTER replicate for one h2 x information scenario
# ----------------------------------------------------------------------

one_bootstrap_generator_diagnostic_scenario <- function(
    latent,
    info_bundle,
    info_level,
    h2,
    sigma2_P = 1,
    beta = 0,
    target = 0,
    n_select = 5,
    B = 500,
    seed_boot,
    return_inner = FALSE) {

  base <- info_bundle$base_design

  dat <- make_dataset_from_latent(
    latent = latent,
    info_level = info_level,
    sigma2_P = sigma2_P,
    h2 = h2,
    beta = beta
  )

  # ------------------------------------------------------------
  # Observed-data REML + EBLUP selection
  # ------------------------------------------------------------

  fit_sel <- fit_reml_animal(
    dat$y,
    info_level$prep
  )

  pred <- predict_ebv_validation(
    fit_sel$eta,
    dat$y,
    info_level$prep
  )

  candidates <- base$candidates

  ord <- order(
    pred$uhat[candidates],
    decreasing = TRUE
  )

  selected <- candidates[
    ord[seq_len(n_select)]
  ]

  cstar <- make_policy(
    n_animals = nrow(base$ped),
    candidates = candidates,
    selected = selected
  )

  truth <- sum(
    cstar * dat$u
  )

  # ------------------------------------------------------------
  # Observed-data plug-in conditional reference
  # ------------------------------------------------------------

  gs_plugin <- policy_stats_animal(
    fit_sel$eta,
    dat$y,
    info_level$prep,
    cstar,
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

  # ------------------------------------------------------------
  # True-VC conditional oracle reference: same c(Y), no reselection
  # ------------------------------------------------------------

  eta_true <- log(
    c(
      dat$sigma2_A,
      dat$sigma2_e
    )
  )

  conditional_trueVC_oracle <- bootdiag_fixed_eta_interval(
    eta = eta_true,
    y = dat$y,
    prep = info_level$prep,
    cstar = cstar,
    target = target
  )

  # ------------------------------------------------------------
  # The fixed-effect generator is intentionally common to both
  # bootstrap worlds. Only the covariance generator changes.
  # ------------------------------------------------------------

  beta_hat_plugin <- bootdiag_gls_beta(
    y = dat$y,
    prep = info_level$prep,
    sigma2_A = fit_sel$theta[1],
    sigma2_e = fit_sel$theta[2]
  )

  # Same standardized inner draws for both generators.
  set.seed(seed_boot)

  z_joint <- matrix(
    rnorm(B * (info_level$prep$nobs + 1L)),
    nrow = B,
    ncol = info_level$prep$nobs + 1L
  )

  # ------------------------------------------------------------
  # Bootstrap A: conventional plug-in generator
  # ------------------------------------------------------------

  boot_plugin <- parametric_bootstrap_fixed_policy_generator(
    y = dat$y,
    prep = info_level$prep,
    cstar = cstar,
    sigma2_A_gen = fit_sel$theta[1],
    sigma2_e_gen = fit_sel$theta[2],
    B = B,
    target = target,
    current_fit = fit_sel,
    beta_generator = beta_hat_plugin,
    z_joint = z_joint,
    return_inner = return_inner
  )

  # ------------------------------------------------------------
  # Bootstrap B: TRUE-VC generator (oracle diagnostic)
  # Inner estimation is still REML, exactly as above.
  # ------------------------------------------------------------

  boot_trueVC_gen <- parametric_bootstrap_fixed_policy_generator(
    y = dat$y,
    prep = info_level$prep,
    cstar = cstar,
    sigma2_A_gen = dat$sigma2_A,
    sigma2_e_gen = dat$sigma2_e,
    B = B,
    target = target,
    current_fit = fit_sel,
    beta_generator = beta_hat_plugin,
    z_joint = z_joint,
    return_inner = return_inner
  )

  h2_hat_outer <- (
    fit_sel$theta[1] /
      sum(fit_sel$theta)
  )

  out <- list(
    truth = truth,
    h2 = h2,
    information = info_level$name,
    n_pheno = info_level$n_pheno,
    selected = selected,
    cstar = cstar,
    sigma2_A_true = dat$sigma2_A,
    sigma2_e_true = dat$sigma2_e,
    sigma2_A_hat = fit_sel$theta[1],
    sigma2_e_hat = fit_sel$theta[2],
    h2_hat = h2_hat_outer,
    lower_boundary_outer = as.numeric(
      fit_sel$eta[1] - fit_sel$lower[1] < 0.05
    ),
    conditional_plugin =
      bootdiag_pack_method(
        conditional_plugin
      ),
    conditional_trueVC_oracle =
      bootdiag_pack_method(
        conditional_trueVC_oracle
      ),
    bootstrap_plugin_generator =
      bootdiag_pack_method(
        boot_plugin
      ),
    bootstrap_trueVC_generator =
      bootdiag_pack_method(
        boot_trueVC_gen
      ),
    bootstrap_plugin_diagnostics = boot_plugin,
    bootstrap_trueVC_diagnostics = boot_trueVC_gen
  )

  out
}


# ----------------------------------------------------------------------
# 7. Standard method-performance summary
# ----------------------------------------------------------------------

summarize_bootstrap_generator_method <- function(
    res,
    method,
    target = 0) {

  x <- lapply(
    res,
    `[[`,
    method
  )

  ok <- !vapply(
    x,
    is.null,
    logical(1)
  )

  x_ok <- x[ok]

  truth <- vapply(
    res,
    `[[`,
    numeric(1),
    "truth"
  )

  if (length(x_ok) == 0L) {
    return(
      c(
        valid = 0,
        bias = NA,
        RMSE = NA,
        mean_SE = NA,
        MSE_ratio = NA,
        inclusion = NA,
        width = NA,
        Brier = NA,
        mean_p = NA,
        event_rate = mean(truth > target)
      )
    )
  }

  truth_ok <- truth[ok]

  get_component <- function(nm) {
    vapply(
      x_ok,
      function(z) z[[nm]],
      numeric(1)
    )
  }

  mu <- get_component("mean")
  se <- get_component("sd")
  lo <- get_component("lo")
  hi <- get_component("hi")
  pp <- get_component("p")

  err <- mu - truth_ok
  mse <- mean(err^2)

  c(
    valid = length(x_ok),
    bias = mean(err),
    RMSE = sqrt(mse),
    mean_SE = mean(se),
    MSE_ratio = if (mse > 0) mean(se^2) / mse else NA_real_,
    inclusion = mean(
      truth_ok >= lo &
        truth_ok <= hi
    ),
    width = mean(hi - lo),
    Brier = mean(
      (pp - as.numeric(truth_ok > target))^2
    ),
    mean_p = mean(pp),
    event_rate = mean(truth_ok > target)
  )
}


summarize_bootstrap_generator_methods <- function(
    res,
    target = 0) {

  methods <- c(
    "conditional_plugin",
    "conditional_trueVC_oracle",
    "bootstrap_plugin_generator",
    "bootstrap_trueVC_generator"
  )

  tab <- t(
    vapply(
      methods,
      function(m) {
        summarize_bootstrap_generator_method(
          res,
          m,
          target = target
        )
      },
      numeric(10)
    )
  )

  data.frame(
    method = methods,
    tab,
    row.names = NULL,
    check.names = FALSE
  )
}


# ----------------------------------------------------------------------
# 8. Bootstrap-mechanism diagnostics
# ----------------------------------------------------------------------

summarize_bootstrap_generator_diagnostics <- function(res) {

  get_outer <- function(nm) {
    vapply(
      res,
      `[[`,
      numeric(1),
      nm
    )
  }

  get_boot_diag <- function(which_boot, nm) {
    vapply(
      res,
      function(z) {
        x <- z[[which_boot]]
        if (is.null(x) ||
            is.null(x$valid) ||
            !isTRUE(x$valid) ||
            is.null(x[[nm]])) {
          return(NA_real_)
        }
        as.numeric(x[[nm]])
      },
      numeric(1)
    )
  }

  h2_true <- get_outer("h2")
  h2_hat <- get_outer("h2_hat")

  c(
    true_h2 = mean(h2_true),
    mean_outer_h2_hat = mean(h2_hat),
    median_outer_h2_hat = median(h2_hat),
    outer_h2_bias = mean(h2_hat - h2_true),
    outer_lower_boundary_rate = mean(
      get_outer("lower_boundary_outer")
    ),
    plugin_inner_mean_h2_hat = mean(
      get_boot_diag(
        "bootstrap_plugin_diagnostics",
        "inner_mean_h2_hat"
      ),
      na.rm = TRUE
    ),
    trueVC_inner_mean_h2_hat = mean(
      get_boot_diag(
        "bootstrap_trueVC_diagnostics",
        "inner_mean_h2_hat"
      ),
      na.rm = TRUE
    ),
    plugin_inner_boundary_rate = mean(
      get_boot_diag(
        "bootstrap_plugin_diagnostics",
        "inner_boundary_rate"
      ),
      na.rm = TRUE
    ),
    trueVC_inner_boundary_rate = mean(
      get_boot_diag(
        "bootstrap_trueVC_diagnostics",
        "inner_boundary_rate"
      ),
      na.rm = TRUE
    ),
    plugin_inner_REML_convergence = mean(
      get_boot_diag(
        "bootstrap_plugin_diagnostics",
        "inner_REML_convergence_rate"
      ),
      na.rm = TRUE
    ),
    trueVC_inner_REML_convergence = mean(
      get_boot_diag(
        "bootstrap_trueVC_diagnostics",
        "inner_REML_convergence_rate"
      ),
      na.rm = TRUE
    ),
    plugin_mean_boot_error_bias = mean(
      get_boot_diag(
        "bootstrap_plugin_diagnostics",
        "error_mean"
      ),
      na.rm = TRUE
    ),
    trueVC_mean_boot_error_bias = mean(
      get_boot_diag(
        "bootstrap_trueVC_diagnostics",
        "error_mean"
      ),
      na.rm = TRUE
    )
  )
}


# ----------------------------------------------------------------------
# 9. Compact comparison row for one scenario
# ----------------------------------------------------------------------

make_bootstrap_generator_comparison <- function(methods) {

  get_row <- function(method) {
    methods[
      methods$method == method,
      ,
      drop = FALSE
    ]
  }

  cond <- get_row("conditional_plugin")
  oracle <- get_row("conditional_trueVC_oracle")
  bp <- get_row("bootstrap_plugin_generator")
  bt <- get_row("bootstrap_trueVC_generator")

  data.frame(
    conditional_coverage = cond$inclusion,
    trueVC_oracle_coverage = oracle$inclusion,
    plugin_bootstrap_coverage = bp$inclusion,
    trueVC_bootstrap_coverage = bt$inclusion,
    trueVC_minus_plugin_bootstrap_coverage =
      bt$inclusion - bp$inclusion,
    plugin_bootstrap_MSE_ratio = bp$MSE_ratio,
    trueVC_bootstrap_MSE_ratio = bt$MSE_ratio,
    plugin_bootstrap_bias = bp$bias,
    trueVC_bootstrap_bias = bt$bias,
    plugin_bootstrap_RMSE = bp$RMSE,
    trueVC_bootstrap_RMSE = bt$RMSE,
    plugin_bootstrap_mean_SE = bp$mean_SE,
    trueVC_bootstrap_mean_SE = bt$mean_SE,
    plugin_bootstrap_width = bp$width,
    trueVC_bootstrap_width = bt$width,
    row.names = NULL
  )
}


# ----------------------------------------------------------------------
# 10. Main production runner
# ----------------------------------------------------------------------

run_bootstrap_generator_diagnostic_final <- function(
    S = 1000,
    h2_values = c(0.05, 0.20, 0.40),
    info_fractions = c(
      low = 0.40,
      moderate = 0.70,
      high = 1.00
    ),
    n_select = 5,
    sigma2_P = 1,
    beta = 0,
    target = 0,
    B = 500,
    seed_info = 20260850,
    seed_outer = 20260851,
    checkpoint_every = 25,
    output_dir = file.path(
      getwd(),
      "simulation_I_II_final_output",
      paste0(
        "Simulation_II_bootstrap_generator_diagnostic_S",
        S,
        "_B",
        B
      )
    ),
    resume = TRUE,
    keep_raw_in_master = FALSE,
    return_inner = FALSE) {

  final_dir_create(output_dir)

  base_design <- build_simII_design()

  info_bundle <- build_information_level_designs(
    base_design = base_design,
    info_fractions = info_fractions,
    seed_info = seed_info
  )

  cat(
    "\n====================================================\n",
    "Simulation II: bootstrap generating-model diagnostic\n",
    "====================================================\n",
    sep = ""
  )

  print_information_design(info_bundle)

  # Recreate the SAME outer latent seeds and method-seed array structure
  # as the previous final Simulation II. This permits direct verification
  # against previous scenario RDS files when the default design is used.
  set.seed(seed_outer)

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

  scenario_results <- list()

  for (ih in seq_along(h2_values)) {

    h2 <- h2_values[ih]

    for (ii in seq_along(info_bundle$info_levels)) {

      info <- info_bundle$info_levels[[ii]]

      key <- paste0(
        "h2_",
        h2,
        "__",
        info$name
      )

      cat(
        "\n====================================================\n",
        "Bootstrap diagnostic: ",
        key,
        " | n_pheno=",
        info$n_pheno,
        "\n====================================================\n",
        sep = ""
      )

      checkpoint_file <- file.path(
        output_dir,
        paste0(
          final_safe_name(key),
          "_bootstrap_generator_checkpoint.rds"
        )
      )

      final_file <- file.path(
        output_dir,
        paste0(
          final_safe_name(key),
          "_bootstrap_generator_final.rds"
        )
      )

      if (resume && file.exists(final_file)) {

        done <- readRDS(final_file)

        final_compatible <- (
          !is.null(done$S) &&
          !is.null(done$B) &&
          !is.null(done$key) &&
          identical(done$S, S) &&
          identical(done$B, B) &&
          identical(done$key, key)
        )

        if (isTRUE(final_compatible)) {
          cat(
            "Loading completed scenario from compatible final RDS ",
            "(S=", S, ", B=", B, ").\n",
            sep = ""
          )

          scenario_results[[key]] <- done$summary

          if (keep_raw_in_master &&
              !is.null(done$raw)) {
            scenario_results[[key]]$raw <- done$raw
          }

          next
        } else {
          cat(
            "Existing final RDS is incompatible with this run ",
            "(requested S=", S, ", B=", B, "); ignoring it.\n",
            sep = ""
          )
        }
      }

      res <- vector(
        "list",
        S
      )

      start_at <- 1L

      if (resume && file.exists(checkpoint_file)) {

        cp <- readRDS(checkpoint_file)

        if (identical(cp$S, S) &&
            identical(cp$key, key) &&
            identical(cp$B, B)) {

          res <- cp$res

          missing_idx <- which(
            vapply(
              res,
              is.null,
              logical(1)
            )
          )

          if (length(missing_idx) == 0L) {
            start_at <- S + 1L
          } else {
            start_at <- min(missing_idx)
          }

          cat(
            "Resuming at outer replicate ",
            start_at,
            ".\n",
            sep = ""
          )
        }
      }

      if (start_at <= S) {

        for (s in seq.int(start_at, S)) {

          res[[s]] <-
            one_bootstrap_generator_diagnostic_scenario(
              latent = latent_list[[s]],
              info_bundle = info_bundle,
              info_level = info,
              h2 = h2,
              sigma2_P = sigma2_P,
              beta = beta,
              target = target,
              n_select = n_select,
              B = B,
              # Use the same bootstrap seed slot as the previous final
              # Simulation II (method index 2 = bootstrap).
              seed_boot = seed_methods[s, ih, ii, 2],
              return_inner = return_inner
            )

          if (
            s %% checkpoint_every == 0 ||
            s == S
          ) {

            saveRDS(
              list(
                S = S,
                B = B,
                key = key,
                h2 = h2,
                information = info$name,
                res = res
              ),
              checkpoint_file
            )

            cat(
              "Completed ",
              s,
              " / ",
              S,
              "\n",
              sep = ""
            )
          }
        }
      }

      method_summary <-
        summarize_bootstrap_generator_methods(
          res,
          target = target
        )

      method_summary$h2 <- h2
      method_summary$information <- info$name
      method_summary$n_pheno <- info$n_pheno

      method_summary <- method_summary[
        , c(
          "h2",
          "information",
          "n_pheno",
          "method",
          setdiff(
            names(method_summary),
            c(
              "h2",
              "information",
              "n_pheno",
              "method"
            )
          )
        ),
        drop = FALSE
      ]

      comparison <-
        make_bootstrap_generator_comparison(
          method_summary
        )

      comparison$h2 <- h2
      comparison$information <- info$name
      comparison$n_pheno <- info$n_pheno

      comparison <- comparison[
        , c(
          "h2",
          "information",
          "n_pheno",
          setdiff(
            names(comparison),
            c(
              "h2",
              "information",
              "n_pheno"
            )
          )
        ),
        drop = FALSE
      ]

      diagnostic <-
        as.data.frame(
          t(
            summarize_bootstrap_generator_diagnostics(
              res
            )
          ),
          check.names = FALSE
        )

      diagnostic$h2 <- h2
      diagnostic$information <- info$name
      diagnostic$n_pheno <- info$n_pheno

      diagnostic <- diagnostic[
        , c(
          "h2",
          "information",
          "n_pheno",
          setdiff(
            names(diagnostic),
            c(
              "h2",
              "information",
              "n_pheno"
            )
          )
        ),
        drop = FALSE
      ]

      summary_one <- list(
        h2 = h2,
        information = info$name,
        n_pheno = info$n_pheno,
        methods = method_summary,
        comparison = comparison,
        diagnostics = diagnostic
      )

      saveRDS(
        list(
          S = S,
          B = B,
          key = key,
          summary = summary_one,
          raw = res
        ),
        final_file
      )

      scenario_results[[key]] <- summary_one

      if (keep_raw_in_master) {
        scenario_results[[key]]$raw <- res
      }
    }
  }

  method_table <- do.call(
    rbind,
    lapply(
      scenario_results,
      `[[`,
      "methods"
    )
  )

  rownames(method_table) <- NULL
  method_table <- final_add_coverage_mcse(method_table)

  comparison_table <- do.call(
    rbind,
    lapply(
      scenario_results,
      `[[`,
      "comparison"
    )
  )

  rownames(comparison_table) <- NULL

  diagnostic_table <- do.call(
    rbind,
    lapply(
      scenario_results,
      `[[`,
      "diagnostics"
    )
  )

  rownames(diagnostic_table) <- NULL

  final_write_csv(
    method_table,
    file.path(
      output_dir,
      "Simulation_II_bootstrap_generator_method_performance.csv"
    )
  )

  final_write_csv(
    comparison_table,
    file.path(
      output_dir,
      "Simulation_II_bootstrap_generator_comparison.csv"
    )
  )

  final_write_csv(
    diagnostic_table,
    file.path(
      output_dir,
      "Simulation_II_bootstrap_generator_diagnostics.csv"
    )
  )

  out <- list(
    design = base_design,
    information_design = info_bundle,
    result = scenario_results,
    method_table = method_table,
    comparison_table = comparison_table,
    diagnostic_table = diagnostic_table,
    settings = list(
      S = S,
      B = B,
      h2_values = h2_values,
      info_fractions = info_fractions,
      n_select = n_select,
      sigma2_P = sigma2_P,
      beta = beta,
      target = target,
      seed_info = seed_info,
      seed_outer = seed_outer,
      output_dir = normalizePath(
        output_dir,
        winslash = "/",
        mustWork = FALSE
      ),
      estimand =
        "true expected genetic superiority c(Y)'u for the REML-EBLUP-selected realised policy; c(Y) fixed within bootstrap",
      main_comparison =
        "plug-in VC generator vs true-VC generator; same outer estimator, same c(Y), same beta generator, same inner standard-normal draws, REML re-estimated inside both bootstraps"
    )
  )

  saveRDS(
    out,
    file.path(
      output_dir,
      "Simulation_II_bootstrap_generator_ALL.rds"
    )
  )

  cat("\nBootstrap generating-model diagnostic completed.\n")
  cat("\nMethod performance:\n")
  print(method_table)
  cat("\nGenerator comparison:\n")
  print(comparison_table)
  cat("\nMechanism diagnostics:\n")
  print(diagnostic_table)

  invisible(out)
}


# ----------------------------------------------------------------------
# 11. Verify exact reproduction of the previous plug-in bootstrap
# ----------------------------------------------------------------------
#
# IMPORTANT: Exact bootstrap reproduction is meaningful only when the diagnostic
# run itself uses the same production S=1000 and B=500 seed schedule. A pilot
# run with S=20 and/or B=100 is expected to differ in packed bootstrap output,
# even though selected sets and truth can match exactly. For a cheap exact check,
# use verify_plugin_bootstrap_reproduction_quick() above.
#
# This check compares selected sets, truth and packed bootstrap output.
# ----------------------------------------------------------------------

verify_bootstrap_generator_against_previous <- function(
    diagnostic_output,
    previous_dir = file.path(
      getwd(),
      "simulation_I_II_final_output",
      "Simulation_II"
    ),
    diagnostic_dir = diagnostic_output$settings$output_dir,
    tolerance = 1e-10) {

  if (!dir.exists(previous_dir)) {
    stop("previous_dir does not exist: ", previous_dir)
  }

  if (is.null(diagnostic_dir) ||
      !dir.exists(diagnostic_dir)) {
    stop("diagnostic_dir could not be resolved.")
  }

  rows <- list()
  k <- 1L

  for (key in names(diagnostic_output$result)) {

    old_file <- file.path(
      previous_dir,
      paste0(
        final_safe_name(key),
        "_final.rds"
      )
    )

    new_file <- file.path(
      diagnostic_dir,
      paste0(
        final_safe_name(key),
        "_bootstrap_generator_final.rds"
      )
    )

    if (!file.exists(old_file) ||
        !file.exists(new_file)) {

      rows[[k]] <- data.frame(
        scenario = key,
        previous_file_found = file.exists(old_file),
        diagnostic_file_found = file.exists(new_file),
        selected_match_rate = NA_real_,
        max_abs_truth_difference = NA_real_,
        max_abs_bootstrap_packed_difference = NA_real_,
        exact_within_tolerance = NA,
        stringsAsFactors = FALSE
      )

      k <- k + 1L
      next
    }

    old_raw <- readRDS(old_file)$raw
    new_raw <- readRDS(new_file)$raw

    if (is.null(old_raw) || is.null(new_raw)) {

      rows[[k]] <- data.frame(
        scenario = key,
        previous_file_found = TRUE,
        diagnostic_file_found = TRUE,
        selected_match_rate = NA_real_,
        max_abs_truth_difference = NA_real_,
        max_abs_bootstrap_packed_difference = NA_real_,
        exact_within_tolerance = NA,
        stringsAsFactors = FALSE
      )

      k <- k + 1L
      next
    }

    n <- min(
      length(old_raw),
      length(new_raw)
    )

    selected_match <- vapply(
      seq_len(n),
      function(i) {
        identical(
          as.integer(old_raw[[i]]$selected),
          as.integer(new_raw[[i]]$selected)
        )
      },
      logical(1)
    )

    truth_diff <- vapply(
      seq_len(n),
      function(i) {
        abs(
          old_raw[[i]]$truth -
            new_raw[[i]]$truth
        )
      },
      numeric(1)
    )

    boot_diff <- vapply(
      seq_len(n),
      function(i) {

        a <- old_raw[[i]]$bootstrap_policy
        b <- new_raw[[i]]$bootstrap_plugin_generator

        if (is.null(a) || is.null(b)) {
          return(NA_real_)
        }

        max(
          abs(
            as.numeric(a) -
              as.numeric(b)
          )
        )
      },
      numeric(1)
    )

    mx_boot <- if (all(is.na(boot_diff))) {
      NA_real_
    } else {
      max(boot_diff, na.rm = TRUE)
    }

    exact <- (
      all(truth_diff <= tolerance) &&
        (is.na(mx_boot) || mx_boot <= tolerance) &&
        all(selected_match)
    )

    rows[[k]] <- data.frame(
      scenario = key,
      previous_file_found = TRUE,
      diagnostic_file_found = TRUE,
      selected_match_rate = mean(selected_match),
      max_abs_truth_difference = max(truth_diff),
      max_abs_bootstrap_packed_difference = mx_boot,
      exact_within_tolerance = exact,
      stringsAsFactors = FALSE
    )

    k <- k + 1L
  }

  do.call(
    rbind,
    rows
  )
}



# ----------------------------------------------------------------------
# 11B. QUICK EXACT REPRODUCTION CHECK AGAINST THE ORIGINAL S=1000, B=500 RUN
# ----------------------------------------------------------------------
#
# IMPORTANT:
# The method-seed array in the original production Simulation II was created
# AFTER drawing S outer latent seeds. Therefore a pilot run with S=20 does NOT
# have the same method seeds as the original S=1000 production run.
# Also, packed bootstrap summaries cannot be exactly equal when B differs.
#
# This helper reconstructs the ORIGINAL production seed schedule with
# S_reference=1000 and B_reference=500, including the otherwise easy-to-miss
# step of generating all S_reference latent replicates BEFORE method seeds.
# It evaluates only a small number of outer replicates after reconstructing
# that RNG state. It computes ONLY the conventional plug-in bootstrap branch,
# so exact reproduction can be checked cheaply before the new diagnostic is run.
#
# Recommended:
#   quick_check <- verify_plugin_bootstrap_reproduction_quick(
#     outer_replicates = 1:2
#   )
#   print(quick_check)
#
# Expected for every scenario:
#   selected_match = TRUE
#   truth_difference = 0 (within numerical tolerance)
#   bootstrap_packed_difference = 0 (within numerical tolerance)
#   exact_within_tolerance = TRUE
# ----------------------------------------------------------------------

verify_plugin_bootstrap_reproduction_quick <- function(
    outer_replicates = 1L,
    S_reference = 1000L,
    B_reference = 500L,
    h2_values = c(0.05, 0.20, 0.40),
    info_fractions = c(
      low = 0.40,
      moderate = 0.70,
      high = 1.00
    ),
    n_select = 5,
    sigma2_P = 1,
    beta = 0,
    target = 0,
    seed_info = 20260850,
    seed_outer = 20260851,
    previous_dir = file.path(
      getwd(),
      "simulation_I_II_final_output",
      "Simulation_II"
    ),
    tolerance = 1e-10) {

  outer_replicates <- as.integer(outer_replicates)

  if (length(outer_replicates) == 0L ||
      any(outer_replicates < 1L) ||
      any(outer_replicates > S_reference)) {
    stop("outer_replicates must be between 1 and S_reference.")
  }

  if (!dir.exists(previous_dir)) {
    stop("previous_dir does not exist: ", previous_dir)
  }

  base_design <- build_simII_design()

  info_bundle <- build_information_level_designs(
    base_design = base_design,
    info_fractions = info_fractions,
    seed_info = seed_info
  )

  n_h <- length(h2_values)
  n_i <- length(info_bundle$info_levels)

  # Reconstruct EXACT original production random-number schedule.
  set.seed(seed_outer)

  seed_latent_reference <- sample.int(
    .Machine$integer.max,
    S_reference
  )

  # IMPORTANT: the original production runner generated ALL latent
  # replicates before sampling method seeds. Because
  # simulate_latent_replicate() calls set.seed() and consumes random
  # numbers, this changes the global RNG state. We must reproduce that
  # step exactly before constructing seed_methods_reference.
  latent_reference <- lapply(
    seed_latent_reference,
    function(ss) {
      simulate_latent_replicate(
        info_bundle,
        seed = ss
      )
    }
  )

  seed_methods_reference <- array(
    sample.int(
      .Machine$integer.max,
      S_reference * n_h * n_i * 3L
    ),
    dim = c(
      S_reference,
      n_h,
      n_i,
      3L
    )
  )

  rows <- list()
  k <- 1L

  for (ih in seq_along(h2_values)) {

    h2 <- h2_values[ih]

    for (ii in seq_along(info_bundle$info_levels)) {

      info <- info_bundle$info_levels[[ii]]

      key <- paste0(
        "h2_",
        h2,
        "__",
        info$name
      )

      old_file <- file.path(
        previous_dir,
        paste0(
          final_safe_name(key),
          "_final.rds"
        )
      )

      if (!file.exists(old_file)) {
        stop("Previous production scenario file not found: ", old_file)
      }

      old_raw <- readRDS(old_file)$raw

      if (is.null(old_raw)) {
        stop("Previous production RDS has no raw results: ", old_file)
      }

      for (s in outer_replicates) {

        latent <- latent_reference[[s]]

        dat <- make_dataset_from_latent(
          latent = latent,
          info_level = info,
          sigma2_P = sigma2_P,
          h2 = h2,
          beta = beta
        )

        fit_sel <- fit_reml_animal(
          dat$y,
          info$prep
        )

        pred <- predict_ebv_validation(
          fit_sel$eta,
          dat$y,
          info$prep
        )

        candidates <- base_design$candidates

        ord <- order(
          pred$uhat[candidates],
          decreasing = TRUE
        )

        selected <- candidates[
          ord[seq_len(n_select)]
        ]

        cstar <- make_policy(
          n_animals = nrow(base_design$ped),
          candidates = candidates,
          selected = selected
        )

        truth <- sum(cstar * dat$u)

        beta_hat_plugin <- bootdiag_gls_beta(
          y = dat$y,
          prep = info$prep,
          sigma2_A = fit_sel$theta[1],
          sigma2_e = fit_sel$theta[2]
        )

        set.seed(
          seed_methods_reference[
            s,
            ih,
            ii,
            2L
          ]
        )

        z_joint <- matrix(
          rnorm(
            B_reference *
              (info$prep$nobs + 1L)
          ),
          nrow = B_reference,
          ncol = info$prep$nobs + 1L
        )

        pb <- parametric_bootstrap_fixed_policy_generator(
          y = dat$y,
          prep = info$prep,
          cstar = cstar,
          sigma2_A_gen = fit_sel$theta[1],
          sigma2_e_gen = fit_sel$theta[2],
          B = B_reference,
          target = target,
          current_fit = fit_sel,
          beta_generator = beta_hat_plugin,
          z_joint = z_joint,
          return_inner = FALSE
        )

        packed_new <- bootdiag_pack_method(pb)
        packed_old <- old_raw[[s]]$bootstrap_policy

        selected_match <- identical(
          as.integer(old_raw[[s]]$selected),
          as.integer(selected)
        )

        truth_difference <- abs(
          old_raw[[s]]$truth - truth
        )

        bootstrap_difference <- if (
          is.null(packed_old) ||
          is.null(packed_new)
        ) {
          NA_real_
        } else {
          max(
            abs(
              as.numeric(packed_old) -
                as.numeric(packed_new)
            )
          )
        }

        exact <- (
          selected_match &&
          truth_difference <= tolerance &&
          !is.na(bootstrap_difference) &&
          bootstrap_difference <= tolerance
        )

        rows[[k]] <- data.frame(
          scenario = key,
          outer_replicate = s,
          selected_match = selected_match,
          truth_difference = truth_difference,
          bootstrap_packed_difference = bootstrap_difference,
          exact_within_tolerance = exact,
          stringsAsFactors = FALSE
        )

        k <- k + 1L
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


# ----------------------------------------------------------------------
# 12. Small structural self-test
# ----------------------------------------------------------------------

self_test_bootstrap_generator_diagnostic <- function() {

  cat("\nRunning bootstrap-generator diagnostic self-test...\n")

  tmp <- file.path(
    tempdir(),
    "bootstrap_generator_diagnostic_selftest"
  )

  unlink(
    tmp,
    recursive = TRUE,
    force = TRUE
  )

  z <- run_bootstrap_generator_diagnostic_final(
    S = 2,
    h2_values = c(0.20),
    info_fractions = c(
      low = 0.50,
      high = 1.00
    ),
    B = 20,
    seed_info = 123,
    seed_outer = 456,
    checkpoint_every = 1,
    output_dir = tmp,
    resume = FALSE,
    keep_raw_in_master = TRUE,
    return_inner = FALSE
  )

  stopifnot(
    all(
      c(
        "conditional_plugin",
        "conditional_trueVC_oracle",
        "bootstrap_plugin_generator",
        "bootstrap_trueVC_generator"
      ) %in% z$method_table$method
    ),
    nrow(z$comparison_table) == 2,
    nrow(z$diagnostic_table) == 2
  )

  cat("\nBootstrap-generator diagnostic self-test completed.\n")

  invisible(z)
}


# ======================================================================
# END BOOTSTRAP GENERATING-MODEL DIAGNOSTIC EXTENSION
# ======================================================================


# ======================================================================
# EXTENSION: OUTER vs TRUE-VC INNER BOOTSTRAP THREE-COMPONENT
#            PREDICTION-ERROR DECOMPOSITION
# ======================================================================
#
# PURPOSE
# -------
# Diagnose WHY the true-VC generating bootstrap still under-covers.
#
# OUTER, for the same realised policy c(Y):
#
#   e_total   = c'uhat(hat(theta)) - c'u
#   e_knownVC = c'uhat(theta_true) - c'u
#   e_VC      = c'{uhat(hat(theta)) - uhat(theta_true)}
#
# so that
#
#   e_total = e_knownVC + e_VC.
#
# TRUE-VC INNER BOOTSTRAP, holding the OUTER c(Y) fixed:
#
#   e_total*   = c'uhat*(hat(theta)*) - G*
#   e_knownVC* = c'uhat*(theta_true)  - G*
#   e_VC*      = c'{uhat*(hat(theta)*) - uhat*(theta_true)}
#
# so that
#
#   e_total* = e_knownVC* + e_VC*.
#
# The central question is whether the true-VC bootstrap reproduces the
# OUTER second moments of each component:
#
#   E[e_knownVC^2], E[e_VC^2], 2E[e_knownVC e_VC].
#
# IMPORTANT
# ---------
# * No reselection is performed inside the bootstrap.
# * The same outer selected set c(Y) is used throughout each inner world.
# * The inner generator is exactly the TRUE-VC generator used in the
#   previous bootstrap-generator diagnostic, including the same plug-in
#   GLS beta generator and the same bootstrap seed schedule.
# * Therefore e_total* can be checked against the previously saved
#   true-VC bootstrap error mean, SD, and quantiles replicate by replicate.
# * This runner re-runs ONLY the true-VC inner bootstrap branch because
#   the previous production run did not retain all inner vectors.
#
# ======================================================================

ovib_dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

ovib_safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

ovib_write_csv <- function(x, path) {
  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  invisible(path)
}

ovib_pop_var <- function(x) {
  m <- mean(x)
  mean((x - m)^2)
}

ovib_pop_cov <- function(x, y) {
  mx <- mean(x)
  my <- mean(y)
  mean((x - mx) * (y - my))
}

ovib_safe_ratio <- function(num, den, tol = 1e-15) {
  if (!is.finite(num) || !is.finite(den) || abs(den) <= tol) {
    return(NA_real_)
  }
  num / den
}

ovib_safe_cor <- function(x, y) {
  if (length(x) < 2L ||
      !all(is.finite(x)) ||
      !all(is.finite(y)) ||
      sd(x) == 0 || sd(y) == 0) {
    return(NA_real_)
  }
  cor(x, y)
}


# ----------------------------------------------------------------------
# Known-VC linear prediction weights
# ----------------------------------------------------------------------
# For fixed theta and fixed c,
#
#   mu_true(y) = c' uhat(theta_true) = h' P y = (P h)' y.
#
# This avoids repeating a full known-VC mixed-model calculation for every
# inner bootstrap sample. It is algebraically identical to
# policy_stats_animal(eta_true, y, prep, cstar)$mu.
# ----------------------------------------------------------------------

ovib_knownVC_weight <- function(
    prep,
    cstar,
    sigma2_A,
    sigma2_e) {

  X <- prep$X
  A <- prep$A
  ids <- prep$pheno_ids

  a_obs_c <- drop(
    A[ids, , drop = FALSE] %*% cstar
  )

  h <- sigma2_A * a_obs_c

  ViX <- Vinv_apply(
    prep,
    sigma2_A,
    sigma2_e,
    X
  )

  Vih <- Vinv_apply(
    prep,
    sigma2_A,
    sigma2_e,
    h
  )

  XtViX <- crossprod(X, ViX)

  Ph <- Vih - ViX %*% solve(
    XtViX,
    crossprod(X, Vih)
  )

  drop(Ph)
}


# ----------------------------------------------------------------------
# One outer replicate: regenerate the TRUE-VC inner bootstrap and split
# its error into known-VC and VC-induced components.
# ----------------------------------------------------------------------

ovib_one_outer <- function(
    outer_raw,
    dat,
    prep,
    cstar,
    B,
    seed_boot,
    tolerance = 1e-10,
    keep_inner_vectors = FALSE) {

  if (is.null(outer_raw$bootstrap_trueVC_diagnostics) ||
      !isTRUE(outer_raw$bootstrap_trueVC_diagnostics$valid)) {
    stop("Previous true-VC bootstrap diagnostics are missing or invalid.")
  }

  sigma2_A_true <- as.numeric(dat$sigma2_A)
  sigma2_e_true <- as.numeric(dat$sigma2_e)

  sigma2_A_hat <- as.numeric(outer_raw$sigma2_A_hat)
  sigma2_e_hat <- as.numeric(outer_raw$sigma2_e_hat)

  # Same beta generator as the previous true-VC bootstrap diagnostic.
  beta_hat_plugin <- bootdiag_gls_beta(
    y = dat$y,
    prep = prep,
    sigma2_A = sigma2_A_hat,
    sigma2_e = sigma2_e_hat
  )

  joint <- bootdiag_joint_generator(
    prep = prep,
    cstar = cstar,
    sigma2_A = sigma2_A_true,
    sigma2_e = sigma2_e_true,
    beta_generator = beta_hat_plugin
  )

  # Same seed and standard-normal draw construction as V3.
  set.seed(seed_boot)

  z_joint <- matrix(
    rnorm(B * (prep$nobs + 1L)),
    nrow = B,
    ncol = prep$nobs + 1L
  )

  draws <- sweep(
    z_joint %*% joint$chol,
    2,
    joint$mean,
    "+"
  )

  # Known-VC predictor is linear in y for fixed theta and c.
  w_true <- ovib_knownVC_weight(
    prep = prep,
    cstar = cstar,
    sigma2_A = sigma2_A_true,
    sigma2_e = sigma2_e_true
  )

  e_total <- rep(NA_real_, B)
  e_knownVC <- rep(NA_real_, B)
  e_VC <- rep(NA_real_, B)

  h2_hat <- rep(NA_real_, B)
  converged <- rep(NA_real_, B)
  boundary <- rep(NA_real_, B)

  for (b in seq_len(B)) {

    yb <- draws[
      b,
      seq_len(prep$nobs)
    ]

    truth_b <- draws[
      b,
      prep$nobs + 1L
    ]

    fit_b <- fit_reml_animal(
      yb,
      prep
    )

    gs_hat <- policy_stats_animal(
      fit_b$eta,
      yb,
      prep,
      cstar,
      calc_diag = FALSE
    )

    mu_hat_b <- as.numeric(gs_hat$mu)
    mu_true_b <- as.numeric(sum(w_true * yb))

    e_total[b] <- mu_hat_b - truth_b
    e_knownVC[b] <- mu_true_b - truth_b
    e_VC[b] <- mu_hat_b - mu_true_b

    h2_hat[b] <- fit_b$theta[1] / sum(fit_b$theta)

    converged[b] <- as.numeric(
      fit_b$convergence == 0
    )

    boundary[b] <- as.numeric(
      fit_b$eta[1] - fit_b$lower[1] < 0.05
    )
  }

  ok <- (
    is.finite(e_total) &
      is.finite(e_knownVC) &
      is.finite(e_VC)
  )

  if (sum(ok) < max(20L, ceiling(0.80 * B))) {
    stop(
      "Too few valid inner replicates: ",
      sum(ok), " / ", B
    )
  }

  et <- e_total[ok]
  ek <- e_knownVC[ok]
  ev <- e_VC[ok]

  identity_error <- max(
    abs(et - (ek + ev))
  )

  if (identity_error > tolerance) {
    stop(
      "Inner error identity failed: max residual = ",
      format(identity_error, digits = 16)
    )
  }

  # ------------------------------------------------------------
  # Exact reproduction check against previous true-VC bootstrap
  # ------------------------------------------------------------

  prev <- outer_raw$bootstrap_trueVC_diagnostics

  q_total <- safe_quantile(
    et,
    c(0.025, 0.975)
  )

  reproduction_diffs <- c(
    error_mean = mean(et) - prev$error_mean,
    error_sd = sd(et) - prev$error_sd,
    error_q025 = q_total[1] - prev$error_q025,
    error_q975 = q_total[2] - prev$error_q975,
    mean_h2 = mean(h2_hat[ok], na.rm = TRUE) - prev$inner_mean_h2_hat,
    boundary_rate = mean(boundary[ok], na.rm = TRUE) - prev$inner_boundary_rate,
    convergence_rate = mean(converged[ok], na.rm = TRUE) - prev$inner_REML_convergence_rate
  )

  max_reproduction_difference <- max(
    abs(reproduction_diffs),
    na.rm = TRUE
  )

  if (max_reproduction_difference > tolerance) {
    stop(
      "Previous true-VC bootstrap was not exactly reproduced. ",
      "max absolute difference = ",
      format(max_reproduction_difference, digits = 16)
    )
  }

  # ------------------------------------------------------------
  # Inner conditional moments for this outer dataset
  # ------------------------------------------------------------

  out <- list(
    n_valid = length(et),
    valid_fraction = length(et) / B,

    bias_total = mean(et),
    bias_knownVC = mean(ek),
    bias_VC = mean(ev),

    E2_total = mean(et^2),
    E2_knownVC = mean(ek^2),
    E2_VC = mean(ev^2),
    E_cross = 2 * mean(ek * ev),

    # Sample variance/covariance: directly aligned with sd(error)^2 used
    # by the bootstrap method summary.
    Var_total_sample = var(et),
    Var_knownVC_sample = var(ek),
    Var_VC_sample = var(ev),
    Cov_cross_sample = 2 * cov(ek, ev),

    # Population versions obey the exact second-moment decomposition with
    # the conditional bootstrap bias terms.
    Var_total_pop = ovib_pop_var(et),
    Var_knownVC_pop = ovib_pop_var(ek),
    Var_VC_pop = ovib_pop_var(ev),
    Cov_cross_pop = 2 * ovib_pop_cov(ek, ev),

    corr_knownVC_VC = ovib_safe_cor(ek, ev),

    mean_h2_hat = mean(h2_hat[ok], na.rm = TRUE),
    boundary_rate = mean(boundary[ok], na.rm = TRUE),
    convergence_rate = mean(converged[ok], na.rm = TRUE),

    max_inner_identity_error = identity_error,
    max_previous_bootstrap_reproduction_difference =
      max_reproduction_difference
  )

  # Exact conditional identities.
  out$E2_identity_error <- (
    out$E2_total -
      (
        out$E2_knownVC +
          out$E2_VC +
          out$E_cross
      )
  )

  out$Var_pop_identity_error <- (
    out$Var_total_pop -
      (
        out$Var_knownVC_pop +
          out$Var_VC_pop +
          out$Cov_cross_pop
      )
  )

  if (keep_inner_vectors) {
    out$e_total <- et
    out$e_knownVC <- ek
    out$e_VC <- ev
    out$h2_hat <- h2_hat[ok]
    out$boundary <- boundary[ok]
    out$converged <- converged[ok]
  }

  out
}


# ----------------------------------------------------------------------
# Summarise OUTER errors for one scenario.
# ----------------------------------------------------------------------

ovib_outer_summary <- function(res) {

  truth <- vapply(
    res,
    function(z) as.numeric(z$truth),
    numeric(1)
  )

  plugin <- vapply(
    res,
    function(z) as.numeric(z$conditional_plugin["mean"]),
    numeric(1)
  )

  knownVC <- vapply(
    res,
    function(z) as.numeric(z$conditional_trueVC_oracle["mean"]),
    numeric(1)
  )

  et <- plugin - truth
  ek <- knownVC - truth
  ev <- plugin - knownVC

  iderr <- max(abs(et - (ek + ev)))

  list(
    bias_total = mean(et),
    bias_knownVC = mean(ek),
    bias_VC = mean(ev),

    MSE_total = mean(et^2),
    MSE_knownVC = mean(ek^2),
    MSE_VC = mean(ev^2),
    MSE_cross = 2 * mean(ek * ev),

    Var_total_pop = ovib_pop_var(et),
    Var_knownVC_pop = ovib_pop_var(ek),
    Var_VC_pop = ovib_pop_var(ev),
    Var_cross_pop = 2 * ovib_pop_cov(ek, ev),

    corr_knownVC_VC = ovib_safe_cor(ek, ev),
    max_identity_error = iderr
  )
}


# ----------------------------------------------------------------------
# Summarise TRUE-VC inner-bootstrap decomposition across outer datasets.
# We average conditional moments over the S outer replicates.
# ----------------------------------------------------------------------

ovib_inner_summary <- function(inner_res) {

  grab <- function(name) {
    vapply(
      inner_res,
      function(z) as.numeric(z[[name]]),
      numeric(1)
    )
  }

  list(
    mean_bias_total = mean(grab("bias_total")),
    mean_bias_knownVC = mean(grab("bias_knownVC")),
    mean_bias_VC = mean(grab("bias_VC")),

    mean_E2_total = mean(grab("E2_total")),
    mean_E2_knownVC = mean(grab("E2_knownVC")),
    mean_E2_VC = mean(grab("E2_VC")),
    mean_E_cross = mean(grab("E_cross")),

    mean_Var_total_sample = mean(grab("Var_total_sample")),
    mean_Var_knownVC_sample = mean(grab("Var_knownVC_sample")),
    mean_Var_VC_sample = mean(grab("Var_VC_sample")),
    mean_Cov_cross_sample = mean(grab("Cov_cross_sample")),

    mean_Var_total_pop = mean(grab("Var_total_pop")),
    mean_Var_knownVC_pop = mean(grab("Var_knownVC_pop")),
    mean_Var_VC_pop = mean(grab("Var_VC_pop")),
    mean_Cov_cross_pop = mean(grab("Cov_cross_pop")),

    mean_corr_knownVC_VC = mean(
      grab("corr_knownVC_VC"),
      na.rm = TRUE
    ),

    mean_h2_hat = mean(grab("mean_h2_hat"), na.rm = TRUE),
    mean_boundary_rate = mean(grab("boundary_rate"), na.rm = TRUE),
    mean_convergence_rate = mean(grab("convergence_rate"), na.rm = TRUE),

    max_inner_identity_error = max(
      grab("max_inner_identity_error"),
      na.rm = TRUE
    ),

    max_E2_identity_error = max(
      abs(grab("E2_identity_error")),
      na.rm = TRUE
    ),

    max_Var_pop_identity_error = max(
      abs(grab("Var_pop_identity_error")),
      na.rm = TRUE
    ),

    max_previous_bootstrap_reproduction_difference = max(
      grab("max_previous_bootstrap_reproduction_difference"),
      na.rm = TRUE
    )
  )
}


# ----------------------------------------------------------------------
# Build one scenario comparison row.
# ----------------------------------------------------------------------

ovib_compare_one <- function(
    h2,
    information,
    n_pheno,
    outer,
    inner,
    existing_trueVC_bootstrap_MSE_ratio = NA_real_) {

  data.frame(
    h2 = h2,
    information = information,
    n_pheno = n_pheno,

    outer_bias_total = outer$bias_total,
    outer_bias_knownVC = outer$bias_knownVC,
    outer_bias_VC = outer$bias_VC,

    inner_mean_bias_total = inner$mean_bias_total,
    inner_mean_bias_knownVC = inner$mean_bias_knownVC,
    inner_mean_bias_VC = inner$mean_bias_VC,

    outer_MSE_total = outer$MSE_total,
    outer_MSE_knownVC = outer$MSE_knownVC,
    outer_MSE_VC = outer$MSE_VC,
    outer_MSE_cross = outer$MSE_cross,

    inner_mean_E2_total = inner$mean_E2_total,
    inner_mean_E2_knownVC = inner$mean_E2_knownVC,
    inner_mean_E2_VC = inner$mean_E2_VC,
    inner_mean_E_cross = inner$mean_E_cross,

    # Direct component reproduction ratios based on second moments.
    E2_reproduction_total = ovib_safe_ratio(
      inner$mean_E2_total,
      outer$MSE_total
    ),
    E2_reproduction_knownVC = ovib_safe_ratio(
      inner$mean_E2_knownVC,
      outer$MSE_knownVC
    ),
    E2_reproduction_VC = ovib_safe_ratio(
      inner$mean_E2_VC,
      outer$MSE_VC
    ),

    RMS_reproduction_total = sqrt(
      pmax(
        ovib_safe_ratio(inner$mean_E2_total, outer$MSE_total),
        0
      )
    ),
    RMS_reproduction_knownVC = sqrt(
      pmax(
        ovib_safe_ratio(inner$mean_E2_knownVC, outer$MSE_knownVC),
        0
      )
    ),
    RMS_reproduction_VC = sqrt(
      pmax(
        ovib_safe_ratio(inner$mean_E2_VC, outer$MSE_VC),
        0
      )
    ),

    # This is the quantity most directly comparable with the existing
    # bootstrap method's MSE_ratio because its reported SE is sd(error*).
    inner_mean_sampleVar_total_over_outer_MSE = ovib_safe_ratio(
      inner$mean_Var_total_sample,
      outer$MSE_total
    ),
    existing_trueVC_bootstrap_MSE_ratio =
      existing_trueVC_bootstrap_MSE_ratio,

    inner_mean_sampleVar_knownVC_over_outer_MSE_knownVC = ovib_safe_ratio(
      inner$mean_Var_knownVC_sample,
      outer$MSE_knownVC
    ),
    inner_mean_sampleVar_VC_over_outer_MSE_VC = ovib_safe_ratio(
      inner$mean_Var_VC_sample,
      outer$MSE_VC
    ),

    outer_corr_knownVC_VC = outer$corr_knownVC_VC,
    inner_mean_corr_knownVC_VC = inner$mean_corr_knownVC_VC,

    inner_mean_h2_hat = inner$mean_h2_hat,
    inner_mean_boundary_rate = inner$mean_boundary_rate,
    inner_mean_convergence_rate = inner$mean_convergence_rate,

    row.names = NULL,
    check.names = FALSE
  )
}


# ----------------------------------------------------------------------
# Long component table, useful for figures and manuscript tables.
# ----------------------------------------------------------------------

ovib_component_table_one <- function(
    h2,
    information,
    n_pheno,
    outer,
    inner) {

  components <- c(
    "total",
    "knownVC",
    "VC",
    "cross"
  )

  outer_value <- c(
    outer$MSE_total,
    outer$MSE_knownVC,
    outer$MSE_VC,
    outer$MSE_cross
  )

  inner_value <- c(
    inner$mean_E2_total,
    inner$mean_E2_knownVC,
    inner$mean_E2_VC,
    inner$mean_E_cross
  )

  ratio <- vapply(
    seq_along(components),
    function(i) {
      ovib_safe_ratio(
        inner_value[i],
        outer_value[i]
      )
    },
    numeric(1)
  )

  data.frame(
    h2 = h2,
    information = information,
    n_pheno = n_pheno,
    component = components,
    outer_second_moment = outer_value,
    inner_mean_second_moment = inner_value,
    inner_over_outer = ratio,
    inner_minus_outer = inner_value - outer_value,
    row.names = NULL,
    check.names = FALSE
  )
}


# ----------------------------------------------------------------------
# Main runner
# ----------------------------------------------------------------------

run_outer_vs_trueVC_inner_decomposition <- function(
    boot1000 = NULL,
    prior_diagnostic_dir = NULL,
    output_dir = NULL,
    checkpoint_every = 10,
    resume = TRUE,
    tolerance = 1e-10,
    keep_per_outer_in_master = FALSE,
    keep_inner_vectors = FALSE) {

  # ------------------------------------------------------------
  # Resolve settings from the completed V3 object when available.
  # ------------------------------------------------------------

  if (!is.null(boot1000)) {

    if (is.null(boot1000$settings)) {
      stop("boot1000 has no $settings element.")
    }

    st <- boot1000$settings

    S <- st$S
    B <- st$B
    h2_values <- st$h2_values
    info_fractions <- st$info_fractions
    n_select <- st$n_select
    sigma2_P <- st$sigma2_P
    beta <- st$beta
    target <- st$target
    seed_info <- st$seed_info
    seed_outer <- st$seed_outer

    if (is.null(prior_diagnostic_dir)) {
      prior_diagnostic_dir <- st$output_dir
    }

  } else {

    S <- 1000L
    B <- 500L
    h2_values <- c(0.05, 0.20, 0.40)
    info_fractions <- c(
      low = 0.40,
      moderate = 0.70,
      high = 1.00
    )
    n_select <- 5L
    sigma2_P <- 1
    beta <- 0
    target <- 0
    seed_info <- 20260850
    seed_outer <- 20260851

    if (is.null(prior_diagnostic_dir)) {
      prior_diagnostic_dir <- file.path(
        getwd(),
        "simulation_I_II_final_output",
        paste0(
          "Simulation_II_bootstrap_generator_diagnostic_S",
          S,
          "_B",
          B
        )
      )
    }
  }

  if (is.null(output_dir)) {
    output_dir <- file.path(
      getwd(),
      "simulation_I_II_final_output",
      paste0(
        "Simulation_II_outer_vs_trueVC_inner_decomposition_S",
        S,
        "_B",
        B
      )
    )
  }

  ovib_dir_create(output_dir)

  if (!dir.exists(prior_diagnostic_dir)) {
    stop(
      "Previous bootstrap-generator diagnostic directory not found: ",
      prior_diagnostic_dir
    )
  }

  # ------------------------------------------------------------
  # Rebuild exact design and exact RNG schedule used in V3.
  # Note that latent_list MUST be generated before seed_methods.
  # ------------------------------------------------------------

  base_design <- build_simII_design()

  info_bundle <- build_information_level_designs(
    base_design = base_design,
    info_fractions = info_fractions,
    seed_info = seed_info
  )

  cat(
    "\n====================================================\n",
    "Outer vs TRUE-VC inner bootstrap error decomposition\n",
    "====================================================\n",
    sep = ""
  )

  print_information_design(info_bundle)

  set.seed(seed_outer)

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

  comparison_rows <- list()
  component_rows <- list()
  scenario_objects <- list()

  row_id <- 1L

  for (ih in seq_along(h2_values)) {

    h2 <- h2_values[ih]

    for (ii in seq_along(info_bundle$info_levels)) {

      info <- info_bundle$info_levels[[ii]]

      key <- paste0(
        "h2_",
        h2,
        "__",
        info$name
      )

      cat(
        "\n====================================================\n",
        "Three-component diagnostic: ",
        key,
        " | n_pheno=",
        info$n_pheno,
        "\n====================================================\n",
        sep = ""
      )

      prior_file <- file.path(
        prior_diagnostic_dir,
        paste0(
          final_safe_name(key),
          "_bootstrap_generator_final.rds"
        )
      )

      if (!file.exists(prior_file)) {
        stop(
          "Previous scenario final RDS not found: ",
          prior_file
        )
      }

      prior <- readRDS(prior_file)

      compatible <- (
        !is.null(prior$S) &&
          !is.null(prior$B) &&
          !is.null(prior$key) &&
          identical(prior$S, S) &&
          identical(prior$B, B) &&
          identical(prior$key, key) &&
          !is.null(prior$raw) &&
          length(prior$raw) == S
      )

      if (!isTRUE(compatible)) {
        stop(
          "Previous scenario RDS is incompatible with requested run: ",
          key
        )
      }

      outer_raw <- prior$raw

      checkpoint_file <- file.path(
        output_dir,
        paste0(
          ovib_safe_name(key),
          "_outer_vs_inner_checkpoint.rds"
        )
      )

      final_file <- file.path(
        output_dir,
        paste0(
          ovib_safe_name(key),
          "_outer_vs_inner_final.rds"
        )
      )

      if (resume && file.exists(final_file)) {

        done <- readRDS(final_file)

        compatible_final <- (
          !is.null(done$S) &&
            !is.null(done$B) &&
            !is.null(done$key) &&
            identical(done$S, S) &&
            identical(done$B, B) &&
            identical(done$key, key)
        )

        if (isTRUE(compatible_final)) {
          cat("Loading completed compatible decomposition.\n")

          comparison_rows[[row_id]] <- done$comparison
          component_rows[[row_id]] <- done$component_table
          scenario_objects[[key]] <- done$summary

          if (keep_per_outer_in_master && !is.null(done$per_outer)) {
            scenario_objects[[key]]$per_outer <- done$per_outer
          }

          row_id <- row_id + 1L
          next
        }
      }

      inner_res <- vector("list", S)
      start_at <- 1L

      if (resume && file.exists(checkpoint_file)) {

        cp <- readRDS(checkpoint_file)

        compatible_cp <- (
          !is.null(cp$S) &&
            !is.null(cp$B) &&
            !is.null(cp$key) &&
            identical(cp$S, S) &&
            identical(cp$B, B) &&
            identical(cp$key, key)
        )

        if (isTRUE(compatible_cp)) {

          inner_res <- cp$inner_res

          missing_idx <- which(
            vapply(
              inner_res,
              is.null,
              logical(1)
            )
          )

          if (length(missing_idx) == 0L) {
            start_at <- S + 1L
          } else {
            start_at <- min(missing_idx)
          }

          cat(
            "Resuming at outer replicate ",
            start_at,
            ".\n",
            sep = ""
          )
        }
      }

      if (start_at <= S) {

        for (s in seq.int(start_at, S)) {

          dat <- make_dataset_from_latent(
            latent = latent_list[[s]],
            info_level = info,
            sigma2_P = sigma2_P,
            h2 = h2,
            beta = beta
          )

          zouter <- outer_raw[[s]]

          # ------------------------------------------------------
          # Strong outer identity checks before any inner work.
          # ------------------------------------------------------

          cstar <- as.numeric(zouter$cstar)

          truth_rebuilt <- sum(
            cstar * dat$u
          )

          truth_diff <- abs(
            truth_rebuilt - as.numeric(zouter$truth)
          )

          if (truth_diff > tolerance) {
            stop(
              "Outer truth mismatch at scenario ", key,
              ", replicate ", s,
              ": ", format(truth_diff, digits = 16)
            )
          }

          if (length(cstar) != nrow(base_design$ped)) {
            stop(
              "cstar length mismatch at scenario ", key,
              ", replicate ", s
            )
          }

          # Recompute the two outer point predictions to verify the data
          # reconstruction and saved variance components.
          eta_hat_outer <- log(
            c(
              zouter$sigma2_A_hat,
              zouter$sigma2_e_hat
            )
          )

          eta_true_outer <- log(
            c(
              dat$sigma2_A,
              dat$sigma2_e
            )
          )

          mu_hat_check <- policy_stats_animal(
            eta_hat_outer,
            dat$y,
            info$prep,
            cstar,
            calc_diag = FALSE
          )$mu

          mu_true_check <- policy_stats_animal(
            eta_true_outer,
            dat$y,
            info$prep,
            cstar,
            calc_diag = FALSE
          )$mu

          mu_hat_diff <- abs(
            as.numeric(mu_hat_check) -
              as.numeric(zouter$conditional_plugin["mean"])
          )

          mu_true_diff <- abs(
            as.numeric(mu_true_check) -
              as.numeric(zouter$conditional_trueVC_oracle["mean"])
          )

          if (max(mu_hat_diff, mu_true_diff) > tolerance) {
            stop(
              "Outer point-prediction reproduction failed at scenario ",
              key, ", replicate ", s,
              ". max diff = ",
              format(max(mu_hat_diff, mu_true_diff), digits = 16)
            )
          }

          inner_res[[s]] <- ovib_one_outer(
            outer_raw = zouter,
            dat = dat,
            prep = info$prep,
            cstar = cstar,
            B = B,
            seed_boot = seed_methods[s, ih, ii, 2],
            tolerance = tolerance,
            keep_inner_vectors = keep_inner_vectors
          )

          if (
            s %% checkpoint_every == 0L ||
            s == S
          ) {

            saveRDS(
              list(
                S = S,
                B = B,
                key = key,
                h2 = h2,
                information = info$name,
                inner_res = inner_res
              ),
              checkpoint_file
            )

            cat(
              "Completed ",
              s,
              " / ",
              S,
              "\n",
              sep = ""
            )
          }
        }
      }

      # ----------------------------------------------------------
      # Scenario summaries
      # ----------------------------------------------------------

      outer_summary <- ovib_outer_summary(
        outer_raw
      )

      inner_summary <- ovib_inner_summary(
        inner_res
      )

      existing_mse_ratio <- NA_real_

      if (!is.null(prior$summary$methods)) {
        md <- prior$summary$methods
        zz <- md[
          md$method == "bootstrap_trueVC_generator",
          ,
          drop = FALSE
        ]
        if (nrow(zz) == 1L && "MSE_ratio" %in% names(zz)) {
          existing_mse_ratio <- as.numeric(zz$MSE_ratio[1])
        }
      }

      comparison <- ovib_compare_one(
        h2 = h2,
        information = info$name,
        n_pheno = info$n_pheno,
        outer = outer_summary,
        inner = inner_summary,
        existing_trueVC_bootstrap_MSE_ratio = existing_mse_ratio
      )

      component_table <- ovib_component_table_one(
        h2 = h2,
        information = info$name,
        n_pheno = info$n_pheno,
        outer = outer_summary,
        inner = inner_summary
      )

      # ----------------------------------------------------------
      # Scenario-level validation checks
      # ----------------------------------------------------------

      # IMPORTANT: the existing bootstrap MSE_ratio does NOT use the
      # uncorrected outer plug-in predictor in its denominator.
      # summarize_bootstrap_generator_method() defines the bootstrap point
      # estimate as
      #
      #   current - mean(error*)
      #
      # for each outer data set.  Therefore, exact reproduction of the
      # previously reported bootstrap MSE_ratio must use the MSE of that
      # bias-corrected bootstrap point estimate.
      #
      # By contrast,
      #   mean_Var_total_sample / outer_summary$MSE_total
      # is a NEW diagnostic quantity: it asks how much of the actual outer
      # plug-in prediction MSE is reproduced by the inner bootstrap variance.
      # It is intentionally retained in comparison_table, but it must NOT be
      # forced to equal the previously reported bootstrap MSE_ratio.

      outer_plugin_error_vec <- vapply(
        outer_raw,
        function(z) {
          as.numeric(z$conditional_plugin["mean"]) -
            as.numeric(z$truth)
        },
        numeric(1)
      )

      inner_error_mean_vec <- vapply(
        inner_res,
        function(z) as.numeric(z$bias_total),
        numeric(1)
      )

      bootstrap_point_error_vec <- (
        outer_plugin_error_vec - inner_error_mean_vec
      )

      existing_denominator_recomputed <- mean(
        bootstrap_point_error_vec^2
      )

      existing_mse_ratio_recomputed <- (
        inner_summary$mean_Var_total_sample /
          existing_denominator_recomputed
      )

      mse_ratio_diff <- abs(
        existing_mse_ratio_recomputed - existing_mse_ratio
      )

      validation <- data.frame(
        check = c(
          "outer_error_identity",
          "inner_error_identity",
          "inner_E2_identity",
          "inner_population_variance_identity",
          "previous_trueVC_bootstrap_exact_reproduction",
          "existing_trueVC_MSE_ratio_reproduction"
        ),
        max_abs_error = c(
          outer_summary$max_identity_error,
          inner_summary$max_inner_identity_error,
          inner_summary$max_E2_identity_error,
          inner_summary$max_Var_pop_identity_error,
          inner_summary$max_previous_bootstrap_reproduction_difference,
          mse_ratio_diff
        ),
        within_tolerance = c(
          outer_summary$max_identity_error <= tolerance,
          inner_summary$max_inner_identity_error <= tolerance,
          inner_summary$max_E2_identity_error <= tolerance,
          inner_summary$max_Var_pop_identity_error <= tolerance,
          inner_summary$max_previous_bootstrap_reproduction_difference <= tolerance,
          mse_ratio_diff <= tolerance
        ),
        row.names = NULL
      )

      if (!all(validation$within_tolerance)) {
        # Save and print the scenario-level validation table BEFORE stopping,
        # so a failed check can be diagnosed without repeating the inner run.
        validation_debug_file <- file.path(
          output_dir,
          paste0(
            ovib_safe_name(key),
            "_outer_vs_inner_validation_FAILED.csv"
          )
        )

        ovib_write_csv(validation, validation_debug_file)
        print(validation)

        stop(
          "Scenario validation failed for ",
          key,
          ". Validation table saved to: ",
          validation_debug_file
        )
      }

      per_outer <- data.frame(
        outer_replicate = seq_len(S),
        bias_total_star = vapply(inner_res, `[[`, numeric(1), "bias_total"),
        bias_knownVC_star = vapply(inner_res, `[[`, numeric(1), "bias_knownVC"),
        bias_VC_star = vapply(inner_res, `[[`, numeric(1), "bias_VC"),
        E2_total_star = vapply(inner_res, `[[`, numeric(1), "E2_total"),
        E2_knownVC_star = vapply(inner_res, `[[`, numeric(1), "E2_knownVC"),
        E2_VC_star = vapply(inner_res, `[[`, numeric(1), "E2_VC"),
        E_cross_star = vapply(inner_res, `[[`, numeric(1), "E_cross"),
        Var_total_star = vapply(inner_res, `[[`, numeric(1), "Var_total_sample"),
        Var_knownVC_star = vapply(inner_res, `[[`, numeric(1), "Var_knownVC_sample"),
        Var_VC_star = vapply(inner_res, `[[`, numeric(1), "Var_VC_sample"),
        CrossVar_star = vapply(inner_res, `[[`, numeric(1), "Cov_cross_sample"),
        corr_knownVC_VC_star = vapply(inner_res, `[[`, numeric(1), "corr_knownVC_VC"),
        mean_h2_hat_star = vapply(inner_res, `[[`, numeric(1), "mean_h2_hat"),
        boundary_rate_star = vapply(inner_res, `[[`, numeric(1), "boundary_rate"),
        convergence_rate_star = vapply(inner_res, `[[`, numeric(1), "convergence_rate"),
        row.names = NULL,
        check.names = FALSE
      )

      summary_one <- list(
        h2 = h2,
        information = info$name,
        n_pheno = info$n_pheno,
        outer = outer_summary,
        inner = inner_summary,
        comparison = comparison,
        component_table = component_table,
        validation = validation
      )

      saveRDS(
        list(
          S = S,
          B = B,
          key = key,
          summary = summary_one,
          comparison = comparison,
          component_table = component_table,
          validation = validation,
          per_outer = per_outer,
          inner_res = if (keep_inner_vectors) inner_res else NULL
        ),
        final_file
      )

      comparison_rows[[row_id]] <- comparison
      component_rows[[row_id]] <- component_table
      scenario_objects[[key]] <- summary_one

      if (keep_per_outer_in_master) {
        scenario_objects[[key]]$per_outer <- per_outer
      }

      row_id <- row_id + 1L
    }
  }

  comparison_table <- do.call(
    rbind,
    comparison_rows
  )
  rownames(comparison_table) <- NULL

  component_table <- do.call(
    rbind,
    component_rows
  )
  rownames(component_table) <- NULL

  # ------------------------------------------------------------
  # Compact table: the three key reproduction ratios.
  # ------------------------------------------------------------

  compact_table <- comparison_table[
    , c(
      "h2",
      "information",
      "n_pheno",
      "outer_MSE_total",
      "outer_MSE_knownVC",
      "outer_MSE_VC",
      "inner_mean_E2_total",
      "inner_mean_E2_knownVC",
      "inner_mean_E2_VC",
      "E2_reproduction_total",
      "E2_reproduction_knownVC",
      "E2_reproduction_VC",
      "inner_mean_sampleVar_total_over_outer_MSE",
      "existing_trueVC_bootstrap_MSE_ratio",
      "outer_corr_knownVC_VC",
      "inner_mean_corr_knownVC_VC"
    ),
    drop = FALSE
  ]

  # ------------------------------------------------------------
  # Global validation table.
  # ------------------------------------------------------------

  all_validation <- do.call(
    rbind,
    lapply(
      names(scenario_objects),
      function(k) {
        vv <- scenario_objects[[k]]$validation
        vv$scenario <- k
        vv[, c("scenario", "check", "max_abs_error", "within_tolerance")]
      }
    )
  )

  rownames(all_validation) <- NULL

  global_check <- aggregate(
    max_abs_error ~ check,
    data = all_validation,
    FUN = max
  )

  global_check$within_tolerance <- (
    global_check$max_abs_error <= tolerance
  )

  # ------------------------------------------------------------
  # Write outputs.
  # ------------------------------------------------------------

  ovib_write_csv(
    comparison_table,
    file.path(
      output_dir,
      "outer_vs_trueVC_inner_comparison.csv"
    )
  )

  ovib_write_csv(
    compact_table,
    file.path(
      output_dir,
      "outer_vs_trueVC_inner_compact.csv"
    )
  )

  ovib_write_csv(
    component_table,
    file.path(
      output_dir,
      "outer_vs_trueVC_inner_component_table.csv"
    )
  )

  ovib_write_csv(
    all_validation,
    file.path(
      output_dir,
      "outer_vs_trueVC_inner_validation_all.csv"
    )
  )

  ovib_write_csv(
    global_check,
    file.path(
      output_dir,
      "outer_vs_trueVC_inner_global_check.csv"
    )
  )

  out <- list(
    settings = list(
      S = S,
      B = B,
      h2_values = h2_values,
      info_fractions = info_fractions,
      n_select = n_select,
      sigma2_P = sigma2_P,
      beta = beta,
      target = target,
      seed_info = seed_info,
      seed_outer = seed_outer,
      prior_diagnostic_dir = normalizePath(
        prior_diagnostic_dir,
        winslash = "/",
        mustWork = FALSE
      ),
      output_dir = normalizePath(
        output_dir,
        winslash = "/",
        mustWork = FALSE
      ),
      estimand = paste0(
        "same realised outer policy c(Y); true-VC inner generator; ",
        "no inner reselection"
      )
    ),
    scenarios = scenario_objects,
    comparison_table = comparison_table,
    compact_table = compact_table,
    component_table = component_table,
    validation_table = all_validation,
    global_check = global_check
  )

  saveRDS(
    out,
    file.path(
      output_dir,
      "Simulation_II_outer_vs_trueVC_inner_decomposition_ALL.rds"
    )
  )

  cat(
    "\nOuter vs TRUE-VC inner decomposition completed.\n"
  )

  cat("\nGlobal checks:\n")
  print(global_check)

  cat("\nCompact comparison:\n")
  print(compact_table)

  invisible(out)
}


# ----------------------------------------------------------------------
# Optional small end-to-end self-test.
# This creates a temporary V3 diagnostic first, then decomposes it.
# ----------------------------------------------------------------------

self_test_outer_vs_trueVC_inner_decomposition <- function() {

  tmp <- file.path(
    tempdir(),
    "ovib_selftest"
  )

  unlink(
    tmp,
    recursive = TRUE,
    force = TRUE
  )

  prior_dir <- file.path(tmp, "prior")
  dec_dir <- file.path(tmp, "decomposition")

  z <- run_bootstrap_generator_diagnostic_final(
    S = 2,
    h2_values = c(0.20),
    info_fractions = c(low = 0.50),
    B = 20,
    checkpoint_every = 1,
    output_dir = prior_dir,
    resume = FALSE,
    keep_raw_in_master = FALSE,
    return_inner = FALSE
  )

  d <- run_outer_vs_trueVC_inner_decomposition(
    boot1000 = z,
    prior_diagnostic_dir = prior_dir,
    output_dir = dec_dir,
    checkpoint_every = 1,
    resume = FALSE,
    tolerance = 1e-9
  )

  stopifnot(
    all(d$global_check$within_tolerance),
    nrow(d$compact_table) == 1L
  )

  cat(
    "\nOuter-vs-inner decomposition self-test completed successfully.\n"
  )

  invisible(d)
}


# ======================================================================
# RECOMMENDED PRODUCTION USAGE
# ======================================================================
#
# After the completed V3 object boot1000 exists in the workspace:
#
#   dec_inner1000 <- run_outer_vs_trueVC_inner_decomposition(
#     boot1000 = boot1000,
#     checkpoint_every = 10,
#     resume = TRUE
#   )
#
# Then inspect:
#
#   print(dec_inner1000$global_check)
#   print(dec_inner1000$compact_table)
#   print(dec_inner1000$comparison_table)
#   print(dec_inner1000$component_table)
#
# The most important columns are:
#
#   E2_reproduction_knownVC
#   E2_reproduction_VC
#   inner_mean_sampleVar_total_over_outer_MSE
#   existing_trueVC_bootstrap_MSE_ratio
#
# Interpretation:
#
# * E2_reproduction_knownVC ~ 1, but E2_reproduction_VC << 1
#     -> bootstrap mainly misses the VC-induced prediction-error component.
#
# * E2_reproduction_knownVC << 1 as well
#     -> holding c(Y) fixed also fails to reproduce part of the known-VC
#        post-selection error distribution; selection dependence itself is
#        implicated in addition to VC-induced error.
#
# * Both component ratios ~ 1 but coverage remains low
#     -> second moments are reproduced; investigate interval shape / tails /
#        studentization rather than variance-component magnitude.
#
# ======================================================================

# ======================================================================
# NEXT MECHANISM DIAGNOSTIC
# Simulation II: INNER RESELECTION and VC-induced prediction error
# ======================================================================
#
# PURPOSE
# -------
# The previous outer-vs-inner decomposition showed that the TRUE-VC
# fixed-realised-policy bootstrap reproduces the known-VC component well,
# but reproduces only about 11-12% of the outer VC-induced second moment.
#
# This diagnostic asks whether the missing VC-induced error is generated by
# the interaction
#
#   variance-component estimation -> EBLUP perturbation -> policy construction.
#
# IMPORTANT: this is a MECHANISM DIAGNOSTIC, not a replacement CI.
# Reselecting inside bootstrap changes the inferential object from the observed
# realised policy toward selection-rule behaviour.  We use reselection only to
# diagnose why the fixed-c bootstrap misses VC-induced prediction error.
#
# For each TRUE-VC inner bootstrap dataset y* we calculate three versions of
# the VC-induced shift delta_u* = uhat*(hat(theta)*) - uhat*(theta0):
#
#   1. FIXED OUTER POLICY
#      eVC_fixed* = c_outer' delta_u*
#
#   2. PLUGIN-RESELECTED POLICY  [main mechanism diagnostic]
#      select top 5 using uhat*(hat(theta)*)
#      construct c_plugin*
#      eVC_plugin_selected* = c_plugin*' delta_u*
#
#   3. TRUE-VC-RESELECTED POLICY [negative-control / alignment diagnostic]
#      select top 5 using uhat*(theta0)
#      construct c_trueVC*
#      eVC_trueVC_selected* = c_trueVC*' delta_u*
#
# The main comparison is
#
#   R_VC_fixed
#     = mean_inner[(eVC_fixed*)^2] / mean_outer[eVC^2]
#
#   R_VC_plugin_selected
#     = mean_inner[(eVC_plugin_selected*)^2] / mean_outer[eVC^2].
#
# If R_VC_fixed remains ~0.11-0.12 but R_VC_plugin_selected rises strongly
# toward 1, this directly supports VC estimation x data-dependent policy
# construction as the mechanism missed by the fixed-realised-policy bootstrap.
#
# The code deliberately regenerates the SAME TRUE-VC inner y* draws used in
# the previous diagnostic (same V3 RNG schedule).  Therefore the fixed-policy
# branch should reproduce the previous E2_reproduction_VC when the full S and
# B are used.
# ======================================================================

rsel_safe_ratio <- function(num, den) {
  if (!is.finite(num) || !is.finite(den) || den == 0) return(NA_real_)
  num / den
}

rsel_safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3L || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  suppressWarnings(cor(x, y, method = method))
}

rsel_top_selected <- function(uhat, candidates, n_select) {
  ord <- order(
    uhat[candidates],
    decreasing = TRUE
  )
  candidates[ord[seq_len(n_select)]]
}

rsel_policy_from_selected <- function(
    n_animals,
    candidates,
    selected) {

  make_policy(
    n_animals = n_animals,
    candidates = candidates,
    selected = selected
  )
}

rsel_overlap_rate <- function(a, b) {
  length(intersect(as.integer(a), as.integer(b))) / length(a)
}

rsel_exact_set_match <- function(a, b) {
  identical(sort(as.integer(a)), sort(as.integer(b)))
}

# ----------------------------------------------------------------------
# One outer replicate
# ----------------------------------------------------------------------

rsel_one_outer <- function(
    outer_raw,
    dat,
    prep,
    candidates,
    n_select,
    B,
    B_reference = B,
    seed_boot,
    tolerance = 1e-10,
    keep_inner_vectors = FALSE) {

  if (is.null(outer_raw$bootstrap_trueVC_diagnostics) ||
      !isTRUE(outer_raw$bootstrap_trueVC_diagnostics$valid)) {
    stop("Previous true-VC bootstrap diagnostics are missing or invalid.")
  }

  c_outer <- as.numeric(outer_raw$cstar)
  n_animals <- length(c_outer)

  sigma2_A_true <- as.numeric(dat$sigma2_A)
  sigma2_e_true <- as.numeric(dat$sigma2_e)

  sigma2_A_hat_outer <- as.numeric(outer_raw$sigma2_A_hat)
  sigma2_e_hat_outer <- as.numeric(outer_raw$sigma2_e_hat)

  # Same generator mean as the previous TRUE-VC generator diagnostic.
  beta_hat_plugin <- bootdiag_gls_beta(
    y = dat$y,
    prep = prep,
    sigma2_A = sigma2_A_hat_outer,
    sigma2_e = sigma2_e_hat_outer
  )

  # Reproduce the same joint generator used previously.  Only the y* part is
  # needed for the reselection diagnostic, but using the identical joint
  # construction preserves the previous RNG mapping exactly.
  joint <- bootdiag_joint_generator(
    prep = prep,
    cstar = c_outer,
    sigma2_A = sigma2_A_true,
    sigma2_e = sigma2_e_true,
    beta_generator = beta_hat_plugin
  )

  if (B < 1L || B_reference < B) {
    stop("Require 1 <= B <= B_reference.")
  }

  set.seed(seed_boot)

  # Generate the full reference RNG matrix first, then keep the first B rows.
  # Thus a B=100 pilot uses an exact subset of the B=500 production random
  # numbers rather than a different column-major matrix mapping.
  z_joint_reference <- matrix(
    rnorm(B_reference * (prep$nobs + 1L)),
    nrow = B_reference,
    ncol = prep$nobs + 1L
  )

  z_joint <- z_joint_reference[
    seq_len(B),
    ,
    drop = FALSE
  ]

  draws <- sweep(
    z_joint %*% joint$chol,
    2,
    joint$mean,
    "+"
  )

  eta_true <- log(c(sigma2_A_true, sigma2_e_true))

  # Exact old fixed-c known-VC weight.  This is used only to validate that the
  # fixed branch is identical to the previous decomposition.
  w_true_outer <- ovib_knownVC_weight(
    prep = prep,
    cstar = c_outer,
    sigma2_A = sigma2_A_true,
    sigma2_e = sigma2_e_true
  )

  eVC_fixed <- rep(NA_real_, B)
  eVC_fixed_vector <- rep(NA_real_, B)
  eVC_plugin_selected <- rep(NA_real_, B)
  eVC_trueVC_selected <- rep(NA_real_, B)

  overlap_plugin_true <- rep(NA_real_, B)
  exact_plugin_true <- rep(NA_real_, B)
  rank_cor_plugin_true <- rep(NA_real_, B)

  h2_hat <- rep(NA_real_, B)
  converged <- rep(NA_real_, B)
  boundary <- rep(NA_real_, B)

  for (b in seq_len(B)) {

    yb <- draws[b, seq_len(prep$nobs)]

    fit_b <- fit_reml_animal(
      yb,
      prep
    )

    # EBLUP under REML-estimated VC.
    pred_hat <- predict_ebv_validation(
      fit_b$eta,
      yb,
      prep
    )$uhat

    # BLUP under the true simulation VC, using the same y*.
    pred_true <- predict_ebv_validation(
      eta_true,
      yb,
      prep
    )$uhat

    delta_u <- pred_hat - pred_true

    # ----------------------------------------------------------
    # A. Fixed outer realised policy: exact reproduction branch
    # ----------------------------------------------------------

    gs_hat_fixed <- policy_stats_animal(
      fit_b$eta,
      yb,
      prep,
      c_outer,
      calc_diag = FALSE
    )

    mu_hat_fixed <- as.numeric(gs_hat_fixed$mu)
    mu_true_fixed <- as.numeric(sum(w_true_outer * yb))

    eVC_fixed[b] <- mu_hat_fixed - mu_true_fixed
    eVC_fixed_vector[b] <- sum(c_outer * delta_u)

    # ----------------------------------------------------------
    # B. Reselect using the plugin REML-EBLUP: main diagnostic
    # ----------------------------------------------------------

    selected_plugin <- rsel_top_selected(
      pred_hat,
      candidates,
      n_select
    )

    c_plugin <- rsel_policy_from_selected(
      n_animals = n_animals,
      candidates = candidates,
      selected = selected_plugin
    )

    eVC_plugin_selected[b] <- sum(c_plugin * delta_u)

    # ----------------------------------------------------------
    # C. Reselect using true-VC BLUP: negative-control branch
    # ----------------------------------------------------------

    selected_true <- rsel_top_selected(
      pred_true,
      candidates,
      n_select
    )

    c_true <- rsel_policy_from_selected(
      n_animals = n_animals,
      candidates = candidates,
      selected = selected_true
    )

    eVC_trueVC_selected[b] <- sum(c_true * delta_u)

    overlap_plugin_true[b] <- rsel_overlap_rate(
      selected_plugin,
      selected_true
    )

    exact_plugin_true[b] <- as.numeric(
      rsel_exact_set_match(selected_plugin, selected_true)
    )

    rank_cor_plugin_true[b] <- rsel_safe_cor(
      pred_hat[candidates],
      pred_true[candidates],
      method = "spearman"
    )

    h2_hat[b] <- fit_b$theta[1] / sum(fit_b$theta)
    converged[b] <- as.numeric(fit_b$convergence == 0)
    boundary[b] <- as.numeric(
      fit_b$eta[1] - fit_b$lower[1] < 0.05
    )
  }

  fixed_vector_diff <- abs(eVC_fixed - eVC_fixed_vector)

  ok <- (
    is.finite(eVC_fixed) &
      is.finite(eVC_fixed_vector) &
      is.finite(eVC_plugin_selected) &
      is.finite(eVC_trueVC_selected)
  )

  if (sum(ok) < max(20L, ceiling(0.80 * B))) {
    stop(
      "Too few valid inner reselection replicates: ",
      sum(ok), " / ", B
    )
  }

  ef <- eVC_fixed[ok]
  ep <- eVC_plugin_selected[ok]
  et <- eVC_trueVC_selected[ok]

  max_fixed_vector_identity_error <- max(
    fixed_vector_diff[ok],
    na.rm = TRUE
  )

  if (max_fixed_vector_identity_error > tolerance) {
    stop(
      "Fixed-c VC shift identity failed: max residual = ",
      format(max_fixed_vector_identity_error, digits = 16)
    )
  }

  out <- list(
    n_valid = length(ef),
    valid_fraction = length(ef) / B,

    bias_fixed = mean(ef),
    bias_plugin_selected = mean(ep),
    bias_trueVC_selected = mean(et),

    E2_fixed = mean(ef^2),
    E2_plugin_selected = mean(ep^2),
    E2_trueVC_selected = mean(et^2),

    Var_fixed = var(ef),
    Var_plugin_selected = var(ep),
    Var_trueVC_selected = var(et),

    mean_abs_fixed = mean(abs(ef)),
    mean_abs_plugin_selected = mean(abs(ep)),
    mean_abs_trueVC_selected = mean(abs(et)),

    mean_overlap_plugin_true = mean(
      overlap_plugin_true[ok],
      na.rm = TRUE
    ),
    exact_set_match_rate_plugin_true = mean(
      exact_plugin_true[ok],
      na.rm = TRUE
    ),
    mean_rank_cor_plugin_true = mean(
      rank_cor_plugin_true[ok],
      na.rm = TRUE
    ),

    mean_h2_hat = mean(h2_hat[ok], na.rm = TRUE),
    boundary_rate = mean(boundary[ok], na.rm = TRUE),
    convergence_rate = mean(converged[ok], na.rm = TRUE),

    max_fixed_vector_identity_error = max_fixed_vector_identity_error
  )

  if (keep_inner_vectors) {
    out$eVC_fixed <- ef
    out$eVC_plugin_selected <- ep
    out$eVC_trueVC_selected <- et
    out$overlap_plugin_true <- overlap_plugin_true[ok]
    out$exact_plugin_true <- exact_plugin_true[ok]
    out$rank_cor_plugin_true <- rank_cor_plugin_true[ok]
  }

  out
}

# ----------------------------------------------------------------------
# Aggregate one scenario
# ----------------------------------------------------------------------

rsel_summarize_scenario <- function(
    inner_res,
    outer_raw_subset,
    h2,
    information,
    n_pheno,
    reference_fixed_ratio = NA_real_) {

  grab <- function(name) {
    vapply(
      inner_res,
      function(z) as.numeric(z[[name]]),
      numeric(1)
    )
  }

  outer <- ovib_outer_summary(outer_raw_subset)

  mean_E2_fixed <- mean(grab("E2_fixed"))
  mean_E2_plugin <- mean(grab("E2_plugin_selected"))
  mean_E2_true <- mean(grab("E2_trueVC_selected"))

  R_fixed <- rsel_safe_ratio(
    mean_E2_fixed,
    outer$MSE_VC
  )

  R_plugin <- rsel_safe_ratio(
    mean_E2_plugin,
    outer$MSE_VC
  )

  R_true <- rsel_safe_ratio(
    mean_E2_true,
    outer$MSE_VC
  )

  data.frame(
    h2 = h2,
    information = information,
    n_pheno = n_pheno,
    n_outer_inner = length(inner_res),

    outer_bias_VC = outer$bias_VC,
    outer_MSE_VC = outer$MSE_VC,

    inner_bias_fixed = mean(grab("bias_fixed")),
    inner_bias_plugin_selected = mean(grab("bias_plugin_selected")),
    inner_bias_trueVC_selected = mean(grab("bias_trueVC_selected")),

    inner_mean_E2_fixed = mean_E2_fixed,
    inner_mean_E2_plugin_selected = mean_E2_plugin,
    inner_mean_E2_trueVC_selected = mean_E2_true,

    R_VC_fixed = R_fixed,
    R_VC_plugin_selected = R_plugin,
    R_VC_trueVC_selected = R_true,

    reselection_amplification = rsel_safe_ratio(
      mean_E2_plugin,
      mean_E2_fixed
    ),

    plugin_vs_true_selection_E2_ratio = rsel_safe_ratio(
      mean_E2_plugin,
      mean_E2_true
    ),

    mean_overlap_plugin_true = mean(
      grab("mean_overlap_plugin_true"),
      na.rm = TRUE
    ),
    exact_set_match_rate_plugin_true = mean(
      grab("exact_set_match_rate_plugin_true"),
      na.rm = TRUE
    ),
    mean_rank_cor_plugin_true = mean(
      grab("mean_rank_cor_plugin_true"),
      na.rm = TRUE
    ),

    mean_inner_h2_hat = mean(
      grab("mean_h2_hat"),
      na.rm = TRUE
    ),
    mean_inner_boundary_rate = mean(
      grab("boundary_rate"),
      na.rm = TRUE
    ),
    mean_inner_convergence_rate = mean(
      grab("convergence_rate"),
      na.rm = TRUE
    ),

    max_fixed_vector_identity_error = max(
      grab("max_fixed_vector_identity_error"),
      na.rm = TRUE
    ),

    reference_previous_R_VC_fixed = reference_fixed_ratio,
    reference_fixed_ratio_difference = if (
      is.finite(reference_fixed_ratio)
    ) {
      abs(R_fixed - reference_fixed_ratio)
    } else {
      NA_real_
    },

    row.names = NULL,
    check.names = FALSE
  )
}

# ----------------------------------------------------------------------
# Main runner
# ----------------------------------------------------------------------

run_inner_reselection_VC_diagnostic <- function(
    boot1000,
    dec_inner1000 = NULL,
    outer_indices = NULL,
    B_inner = NULL,
    prior_diagnostic_dir = NULL,
    output_dir = NULL,
    checkpoint_every = 10,
    resume = TRUE,
    tolerance = 1e-10,
    keep_inner_vectors = FALSE,
    keep_per_outer_in_master = FALSE) {

  if (is.null(boot1000) || is.null(boot1000$settings)) {
    stop("boot1000 with a $settings element is required.")
  }

  st <- boot1000$settings

  S <- as.integer(st$S)
  B_reference <- as.integer(st$B)

  if (is.null(B_inner)) {
    B <- B_reference
  } else {
    B <- as.integer(B_inner)
  }

  if (!is.finite(B) || B < 20L || B > B_reference) {
    stop(
      "B_inner must be between 20 and the original boot1000 B=",
      B_reference,
      "."
    )
  }

  h2_values <- st$h2_values
  info_fractions <- st$info_fractions
  n_select <- as.integer(st$n_select)
  sigma2_P <- st$sigma2_P
  beta <- st$beta
  seed_info <- st$seed_info
  seed_outer <- st$seed_outer

  if (is.null(outer_indices)) {
    outer_indices <- seq_len(S)
  }

  outer_indices <- sort(unique(as.integer(outer_indices)))

  if (length(outer_indices) == 0L ||
      any(outer_indices < 1L) ||
      any(outer_indices > S)) {
    stop("outer_indices must be within 1:S.")
  }

  full_outer_run <- identical(
    outer_indices,
    seq_len(S)
  )

  if (is.null(prior_diagnostic_dir)) {
    prior_diagnostic_dir <- st$output_dir
  }

  if (is.null(prior_diagnostic_dir) ||
      !dir.exists(prior_diagnostic_dir)) {
    stop(
      "Previous bootstrap-generator diagnostic directory not found: ",
      prior_diagnostic_dir
    )
  }

  if (is.null(output_dir)) {
    output_dir <- file.path(
      getwd(),
      "simulation_I_II_final_output",
      paste0(
        "Simulation_II_inner_reselection_VC_diagnostic_Sinner",
        length(outer_indices),
        "_B",
        B
      )
    )
  }

  ovib_dir_create(output_dir)

  # ------------------------------------------------------------
  # Rebuild the exact V3 design and RNG schedule.
  # ------------------------------------------------------------

  base_design <- build_simII_design()

  info_bundle <- build_information_level_designs(
    base_design = base_design,
    info_fractions = info_fractions,
    seed_info = seed_info
  )

  candidates <- base_design$candidates

  cat(
    "\n====================================================\n",
    "INNER RESELECTION VC-ERROR DIAGNOSTIC\n",
    "====================================================\n",
    "Inner outer-replicate count: ", length(outer_indices), " / ", S, "\n",
    "B per outer replicate: ", B, "\n",
    "This is a mechanism diagnostic; it is NOT a replacement CI.\n",
    sep = ""
  )

  print_information_design(info_bundle)

  set.seed(seed_outer)

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

  scenario_rows <- list()
  scenario_objects <- list()
  row_id <- 1L

  for (ih in seq_along(h2_values)) {

    h2 <- h2_values[ih]

    for (ii in seq_along(info_bundle$info_levels)) {

      info <- info_bundle$info_levels[[ii]]

      key <- paste0(
        "h2_",
        h2,
        "__",
        info$name
      )

      cat(
        "\n====================================================\n",
        "Reselection diagnostic: ", key,
        " | n_pheno=", info$n_pheno,
        "\n====================================================\n",
        sep = ""
      )

      prior_file <- file.path(
        prior_diagnostic_dir,
        paste0(
          final_safe_name(key),
          "_bootstrap_generator_final.rds"
        )
      )

      if (!file.exists(prior_file)) {
        stop("Previous scenario final RDS not found: ", prior_file)
      }

      prior <- readRDS(prior_file)

      compatible <- (
        !is.null(prior$S) &&
          !is.null(prior$B) &&
          !is.null(prior$key) &&
          identical(as.integer(prior$S), S) &&
          identical(as.integer(prior$B), B_reference) &&
          identical(prior$key, key) &&
          !is.null(prior$raw) &&
          length(prior$raw) == S
      )

      if (!isTRUE(compatible)) {
        stop(
          "Previous scenario RDS is incompatible with this run: ",
          key
        )
      }

      outer_raw_full <- prior$raw
      outer_raw_subset <- outer_raw_full[outer_indices]

      scenario_tag <- paste0(
        ovib_safe_name(key),
        "__nOuter",
        length(outer_indices)
      )

      checkpoint_file <- file.path(
        output_dir,
        paste0(
          scenario_tag,
          "_reselection_checkpoint.rds"
        )
      )

      final_file <- file.path(
        output_dir,
        paste0(
          scenario_tag,
          "_reselection_final.rds"
        )
      )

      if (resume && file.exists(final_file)) {

        done <- readRDS(final_file)

        compatible_final <- (
          !is.null(done$S) &&
            !is.null(done$B) &&
            !is.null(done$key) &&
            !is.null(done$outer_indices) &&
            identical(as.integer(done$S), S) &&
            identical(as.integer(done$B), B) &&
            identical(done$key, key) &&
            identical(as.integer(done$outer_indices), outer_indices)
        )

        if (isTRUE(compatible_final)) {
          cat("Loading completed compatible reselection scenario.\n")
          scenario_rows[[row_id]] <- done$summary_table
          scenario_objects[[key]] <- done$summary_object
          row_id <- row_id + 1L
          next
        }
      }

      inner_res <- vector("list", length(outer_indices))
      start_pos <- 1L

      if (resume && file.exists(checkpoint_file)) {

        cp <- readRDS(checkpoint_file)

        compatible_cp <- (
          !is.null(cp$S) &&
            !is.null(cp$B) &&
            !is.null(cp$key) &&
            !is.null(cp$outer_indices) &&
            identical(as.integer(cp$S), S) &&
            identical(as.integer(cp$B), B) &&
            identical(cp$key, key) &&
            identical(as.integer(cp$outer_indices), outer_indices)
        )

        if (isTRUE(compatible_cp)) {

          inner_res <- cp$inner_res

          missing_pos <- which(
            vapply(inner_res, is.null, logical(1))
          )

          if (length(missing_pos) == 0L) {
            start_pos <- length(outer_indices) + 1L
          } else {
            start_pos <- min(missing_pos)
          }

          cat(
            "Resuming at subset position ", start_pos,
            " (outer replicate ",
            if (start_pos <= length(outer_indices)) outer_indices[start_pos] else NA,
            ").\n",
            sep = ""
          )
        }
      }

      if (start_pos <= length(outer_indices)) {

        for (pos in seq.int(start_pos, length(outer_indices))) {

          s <- outer_indices[pos]

          dat <- make_dataset_from_latent(
            latent = latent_list[[s]],
            info_level = info,
            sigma2_P = sigma2_P,
            h2 = h2,
            beta = beta
          )

          zouter <- outer_raw_full[[s]]

          # Strong reconstruction check.
          c_outer <- as.numeric(zouter$cstar)
          truth_rebuilt <- sum(c_outer * dat$u)
          truth_diff <- abs(truth_rebuilt - as.numeric(zouter$truth))

          if (truth_diff > tolerance) {
            stop(
              "Outer truth mismatch at ", key,
              ", replicate ", s,
              ": ", format(truth_diff, digits = 16)
            )
          }

          inner_res[[pos]] <- rsel_one_outer(
            outer_raw = zouter,
            dat = dat,
            prep = info$prep,
            candidates = candidates,
            n_select = n_select,
            B = B,
            B_reference = B_reference,
            # IMPORTANT: index 2 is the same TRUE-VC bootstrap seed used by
            # the previous decomposition and generator diagnostic.
            seed_boot = seed_methods[s, ih, ii, 2],
            tolerance = tolerance,
            keep_inner_vectors = keep_inner_vectors
          )

          if (
            pos %% checkpoint_every == 0L ||
            pos == length(outer_indices)
          ) {

            saveRDS(
              list(
                S = S,
                B = B,
                key = key,
                outer_indices = outer_indices,
                inner_res = inner_res
              ),
              checkpoint_file
            )

            cat(
              "Completed subset position ", pos,
              " / ", length(outer_indices),
              " (outer replicate ", s, ")\n",
              sep = ""
            )
          }
        }
      }

      if (any(vapply(inner_res, is.null, logical(1)))) {
        stop("Scenario contains incomplete inner results: ", key)
      }

      # ----------------------------------------------------------
      # Optional exact reference to the previous decomposition.
      # This comparison is meaningful only for the full S run.
      # ----------------------------------------------------------

      reference_fixed_ratio <- NA_real_

      if (
        full_outer_run &&
        identical(B, B_reference) &&
        !is.null(dec_inner1000) &&
        !is.null(dec_inner1000$compact_table)
      ) {

        dd <- dec_inner1000$compact_table

        idx <- which(
          abs(dd$h2 - h2) < 1e-12 &
            as.character(dd$information) == info$name
        )

        if (length(idx) == 1L &&
            "E2_reproduction_VC" %in% names(dd)) {
          reference_fixed_ratio <- as.numeric(
            dd$E2_reproduction_VC[idx]
          )
        }
      }

      summary_table <- rsel_summarize_scenario(
        inner_res = inner_res,
        outer_raw_subset = outer_raw_subset,
        h2 = h2,
        information = info$name,
        n_pheno = info$n_pheno,
        reference_fixed_ratio = reference_fixed_ratio
      )

      validation <- data.frame(
        check = c(
          "fixed_vector_identity",
          "previous_fixed_R_VC_reproduction"
        ),
        max_abs_error = c(
          summary_table$max_fixed_vector_identity_error,
          summary_table$reference_fixed_ratio_difference
        ),
        within_tolerance = c(
          summary_table$max_fixed_vector_identity_error <= tolerance,
          if (is.finite(summary_table$reference_fixed_ratio_difference)) {
            summary_table$reference_fixed_ratio_difference <= max(tolerance, 1e-8)
          } else {
            NA
          }
        ),
        row.names = NULL
      )

      # Do not discard a long run merely because an optional reference check
      # is unavailable.  Structural fixed-vector identity remains mandatory.
      if (!isTRUE(validation$within_tolerance[1])) {
        stop(
          "Structural reselection validation failed for ", key,
          ". max fixed-vector identity error = ",
          format(validation$max_abs_error[1], digits = 16)
        )
      }

      if (
        isFALSE(validation$within_tolerance[2])
      ) {
        warning(
          "Previous fixed R_VC was not reproduced exactly for ", key,
          ". Inspect validation before interpretation."
        )
      }

      per_outer <- data.frame(
        outer_replicate = outer_indices,
        bias_fixed = vapply(inner_res, `[[`, numeric(1), "bias_fixed"),
        bias_plugin_selected = vapply(inner_res, `[[`, numeric(1), "bias_plugin_selected"),
        bias_trueVC_selected = vapply(inner_res, `[[`, numeric(1), "bias_trueVC_selected"),
        E2_fixed = vapply(inner_res, `[[`, numeric(1), "E2_fixed"),
        E2_plugin_selected = vapply(inner_res, `[[`, numeric(1), "E2_plugin_selected"),
        E2_trueVC_selected = vapply(inner_res, `[[`, numeric(1), "E2_trueVC_selected"),
        mean_overlap_plugin_true = vapply(inner_res, `[[`, numeric(1), "mean_overlap_plugin_true"),
        exact_set_match_rate_plugin_true = vapply(inner_res, `[[`, numeric(1), "exact_set_match_rate_plugin_true"),
        mean_rank_cor_plugin_true = vapply(inner_res, `[[`, numeric(1), "mean_rank_cor_plugin_true"),
        mean_h2_hat = vapply(inner_res, `[[`, numeric(1), "mean_h2_hat"),
        boundary_rate = vapply(inner_res, `[[`, numeric(1), "boundary_rate"),
        convergence_rate = vapply(inner_res, `[[`, numeric(1), "convergence_rate"),
        row.names = NULL,
        check.names = FALSE
      )

      summary_object <- list(
        h2 = h2,
        information = info$name,
        n_pheno = info$n_pheno,
        outer_indices = outer_indices,
        summary_table = summary_table,
        validation = validation
      )

      saveRDS(
        list(
          S = S,
          B = B,
          key = key,
          outer_indices = outer_indices,
          summary_table = summary_table,
          summary_object = summary_object,
          validation = validation,
          per_outer = per_outer,
          inner_res = if (keep_inner_vectors) inner_res else NULL
        ),
        final_file
      )

      write.csv(
        summary_table,
        file.path(
          output_dir,
          paste0(scenario_tag, "_summary.csv")
        ),
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )

      write.csv(
        validation,
        file.path(
          output_dir,
          paste0(scenario_tag, "_validation.csv")
        ),
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )

      write.csv(
        per_outer,
        file.path(
          output_dir,
          paste0(scenario_tag, "_per_outer.csv")
        ),
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )

      scenario_rows[[row_id]] <- summary_table
      scenario_objects[[key]] <- summary_object

      if (keep_per_outer_in_master) {
        scenario_objects[[key]]$per_outer <- per_outer
      }

      row_id <- row_id + 1L
    }
  }

  compact_table <- do.call(rbind, scenario_rows)
  rownames(compact_table) <- NULL

  validation_all <- do.call(
    rbind,
    lapply(
      names(scenario_objects),
      function(k) {
        vv <- scenario_objects[[k]]$validation
        vv$scenario <- k
        vv[, c("scenario", "check", "max_abs_error", "within_tolerance")]
      }
    )
  )

  rownames(validation_all) <- NULL

  # A compact scientific table with only the variables needed for the main
  # interpretation.
  key_table <- compact_table[
    , c(
      "h2",
      "information",
      "n_pheno",
      "outer_MSE_VC",
      "R_VC_fixed",
      "R_VC_plugin_selected",
      "R_VC_trueVC_selected",
      "reselection_amplification",
      "plugin_vs_true_selection_E2_ratio",
      "mean_overlap_plugin_true",
      "exact_set_match_rate_plugin_true",
      "mean_rank_cor_plugin_true",
      "outer_bias_VC",
      "inner_bias_fixed",
      "inner_bias_plugin_selected",
      "inner_bias_trueVC_selected"
    ),
    drop = FALSE
  ]

  write.csv(
    compact_table,
    file.path(output_dir, "inner_reselection_VC_compact_full.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    key_table,
    file.path(output_dir, "inner_reselection_VC_key_table.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    validation_all,
    file.path(output_dir, "inner_reselection_VC_validation_all.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  out <- list(
    settings = list(
      S_reference = S,
      B_reference = B_reference,
      B_inner = B,
      outer_indices = outer_indices,
      full_outer_run = full_outer_run,
      h2_values = h2_values,
      info_fractions = info_fractions,
      n_select = n_select,
      sigma2_P = sigma2_P,
      beta = beta,
      seed_info = seed_info,
      seed_outer = seed_outer,
      prior_diagnostic_dir = normalizePath(
        prior_diagnostic_dir,
        winslash = "/",
        mustWork = FALSE
      ),
      output_dir = normalizePath(
        output_dir,
        winslash = "/",
        mustWork = FALSE
      ),
      diagnostic_target = paste0(
        "mechanism diagnostic: VC-induced prediction error under fixed outer c, ",
        "plugin-EBLUP inner reselection, and true-VC-BLUP inner reselection"
      )
    ),
    compact_table = compact_table,
    key_table = key_table,
    validation = validation_all,
    result = scenario_objects
  )

  saveRDS(
    out,
    file.path(output_dir, "inner_reselection_VC_diagnostic_master.rds")
  )

  cat(
    "\n====================================================\n",
    "Inner reselection VC diagnostic completed.\n",
    "====================================================\n",
    sep = ""
  )

  cat("\nKey table:\n")
  print(key_table)

  cat("\nValidation:\n")
  print(validation_all)

  invisible(out)
}

# ----------------------------------------------------------------------
# Convenience interpretation printer
# ----------------------------------------------------------------------

print_inner_reselection_VC_diagnostic <- function(x, digits = 4) {

  if (is.null(x$key_table)) {
    stop("Object has no $key_table.")
  }

  z <- x$key_table

  cat("\nMain VC reproduction ratios\n")
  print(
    round(
      z[, c(
        "h2",
        "n_pheno",
        "R_VC_fixed",
        "R_VC_plugin_selected",
        "R_VC_trueVC_selected",
        "reselection_amplification"
      )],
      digits
    )
  )

  cat("\nSelection sensitivity diagnostics\n")
  print(
    round(
      z[, c(
        "h2",
        "n_pheno",
        "mean_overlap_plugin_true",
        "exact_set_match_rate_plugin_true",
        "mean_rank_cor_plugin_true"
      )],
      digits
    )
  )

  invisible(z)
}

# ======================================================================
# RECOMMENDED USE
# ======================================================================
#
# 1) QUICK STRUCTURAL CHECK (B=100, only 2 outer replicates/scenario):
#
#   rsel_check <- run_inner_reselection_VC_diagnostic(
#     boot1000 = boot1000,
#     dec_inner1000 = dec_inner1000,
#     outer_indices = 1:2,
#     B_inner = 100,
#     checkpoint_every = 1,
#     resume = FALSE
#   )
#   print(rsel_check$key_table)
#   print(rsel_check$validation)
#
# 2) OPTIONAL PILOT (20 outer replicates/scenario, B=100):
#
#   rsel20 <- run_inner_reselection_VC_diagnostic(
#     boot1000 = boot1000,
#     dec_inner1000 = dec_inner1000,
#     outer_indices = 1:20,
#     B_inner = 100,
#     checkpoint_every = 5,
#     resume = TRUE
#   )
#   print_inner_reselection_VC_diagnostic(rsel20)
#
# 3) FULL PRODUCTION RUN (1000 outer x 500 inner x 9 scenarios):
#
#   rsel1000 <- run_inner_reselection_VC_diagnostic(
#     boot1000 = boot1000,
#     dec_inner1000 = dec_inner1000,
#     checkpoint_every = 10,
#     resume = TRUE
#   )
#
#   print(rsel1000$validation)
#   print(rsel1000$key_table)
#   print_inner_reselection_VC_diagnostic(rsel1000)
#
# INTERPRETATION
# --------------
# - R_VC_fixed should reproduce the previous ~0.11-0.12 in the full run.
# - If R_VC_plugin_selected increases toward 1, the missing VC error is
#   largely attributable to VC-estimation x policy-construction interaction.
# - R_VC_trueVC_selected is a useful negative-control: if plugin-based
#   selection amplifies eVC much more than true-VC-based selection, this is
#   direct evidence that selection is aligning the policy with the VC-induced
#   EBLUP perturbation itself.
# - Do NOT reinterpret the reselection branch as a CI for the originally
#   observed realised policy.  It is only a mechanism diagnostic.
# ======================================================================


# ======================================================================
# EXTENSION: SAME-SAMPLE vs CROSS-SAMPLE POLICY CONSTRUCTION
#            FOR VC-INDUCED PREDICTION ERROR
# ======================================================================
#
# SCIENTIFIC QUESTION
# -------------------
# The previous inner-reselection diagnostic showed:
#
#   R_VC_fixed            ~ 0.11-0.12
#   R_VC_plugin_selected  ~ 1
#   R_VC_trueVC_selected  ~ 0.8-1
#
# This proves that rebuilding a data-dependent policy inside the bootstrap
# restores most of the VC-induced second moment.  However, it does NOT yet
# distinguish two mechanisms:
#
#   (A) marginal policy variation:
#       c* varies from one inner replicate to another;
#
#   (B) same-sample dependence / alignment:
#       c*(y*) is built from the SAME y* that also generates
#       delta_u*(y*) = uhat*(hat theta*) - uhat*(theta0).
#
# This extension separates these mechanisms WITHOUT changing the marginal
# distribution of the selected policy.
#
# For inner replicate b, define
#
#   delta_u_b = uhat_b(hat theta_b) - uhat_b(theta0).
#
# Same-sample plugin selection:
#
#   eVC_same_plugin,b = c_plugin,b' delta_u_b,
#
# where c_plugin,b is the top-k policy from uhat_b(hat theta_b).
#
# Cross-sample plugin selection:
#
#   eVC_cross_plugin,b = c_plugin,pi(b)' delta_u_b,
#
# where pi(b) is a non-zero cyclic shift.  Thus c_plugin,pi(b) has the SAME
# marginal distribution as c_plugin,b, but comes from an independent inner
# replicate and is therefore independent of delta_u_b conditional on the
# outer dataset and generating model.
#
# Several cyclic shifts are averaged.  This uses NO additional REML fits once
# the B inner fits have been computed for an outer dataset.
#
# The same construction is repeated for true-VC-based selection as a control.
#
# MAIN INTERPRETATION
# -------------------
# If
#
#   R_VC_same_plugin ~ 1
#   but
#   R_VC_cross_plugin collapses toward R_VC_fixed (~0.11),
#
# then the large reselection effect is NOT caused simply by c* varying.  It is
# caused by SAME-SAMPLE dependence/alignment between the selected policy and
# the VC-induced EBLUP perturbation.
#
# If R_VC_cross_plugin remains high, then marginal variation of the policy is
# itself sufficient to explain much of the amplification.
#
# IMPORTANT
# ---------
# This remains a mechanism diagnostic.  Cross-sample and same-sample
# reselection branches are NOT proposed confidence intervals for the observed
# realised policy c(Y_obs).
# ======================================================================

crosssel_safe_ratio <- function(num, den) {
  if (!is.finite(num) || !is.finite(den) || den == 0) return(NA_real_)
  num / den
}

crosssel_choose_shifts <- function(B, n_shifts = 5L) {
  B <- as.integer(B)
  n_shifts <- as.integer(n_shifts)

  if (B < 4L) stop("B must be at least 4.")
  if (n_shifts < 1L) stop("n_shifts must be >= 1.")

  # Spread shifts across the full range.  All are non-zero, so the policy
  # always comes from a different inner replicate.
  cand <- unique(
    pmax(
      1L,
      pmin(
        B - 1L,
        as.integer(round(seq(1, B - 1L, length.out = n_shifts)))
      )
    )
  )

  cand[cand %% B != 0L]
}

crosssel_shift_index <- function(B, shift) {
  # For evaluation replicate b, use policy from replicate pi(b)=b+shift mod B.
  ((seq_len(B) - 1L + as.integer(shift)) %% B) + 1L
}

crosssel_policy_contrast <- function(delta_candidate, selected_ids, candidates) {
  # c'delta_u for the policy used throughout this study:
  # mean(selected delta_u) - mean(all-candidate delta_u).
  idx <- match(as.integer(selected_ids), as.integer(candidates))
  if (anyNA(idx)) stop("Selected ID not found in candidate set.")
  mean(delta_candidate[idx]) - mean(delta_candidate)
}

crosssel_one_outer <- function(
    outer_raw,
    dat,
    prep,
    candidates,
    n_select,
    B,
    B_reference = B,
    seed_boot,
    n_cross_shifts = 5L,
    tolerance = 1e-10,
    keep_inner_vectors = FALSE) {

  if (is.null(outer_raw$bootstrap_trueVC_diagnostics) ||
      !isTRUE(outer_raw$bootstrap_trueVC_diagnostics$valid)) {
    stop("Previous true-VC bootstrap diagnostics are missing or invalid.")
  }

  c_outer <- as.numeric(outer_raw$cstar)
  n_animals <- length(c_outer)

  sigma2_A_true <- as.numeric(dat$sigma2_A)
  sigma2_e_true <- as.numeric(dat$sigma2_e)

  sigma2_A_hat_outer <- as.numeric(outer_raw$sigma2_A_hat)
  sigma2_e_hat_outer <- as.numeric(outer_raw$sigma2_e_hat)

  beta_hat_plugin <- bootdiag_gls_beta(
    y = dat$y,
    prep = prep,
    sigma2_A = sigma2_A_hat_outer,
    sigma2_e = sigma2_e_hat_outer
  )

  joint <- bootdiag_joint_generator(
    prep = prep,
    cstar = c_outer,
    sigma2_A = sigma2_A_true,
    sigma2_e = sigma2_e_true,
    beta_generator = beta_hat_plugin
  )

  if (B < 20L || B_reference < B) {
    stop("Require 20 <= B <= B_reference.")
  }

  set.seed(seed_boot)

  # EXACT same RNG mapping as the preceding diagnostics.
  z_joint_reference <- matrix(
    rnorm(B_reference * (prep$nobs + 1L)),
    nrow = B_reference,
    ncol = prep$nobs + 1L
  )

  z_joint <- z_joint_reference[
    seq_len(B),
    ,
    drop = FALSE
  ]

  draws <- sweep(
    z_joint %*% joint$chol,
    2,
    joint$mean,
    "+"
  )

  eta_true <- log(c(sigma2_A_true, sigma2_e_true))

  # Candidate-only storage is sufficient because c is zero outside candidates.
  n_cand <- length(candidates)
  delta_candidate <- matrix(
    NA_real_,
    nrow = B,
    ncol = n_cand
  )

  selected_plugin <- matrix(
    NA_integer_,
    nrow = B,
    ncol = n_select
  )

  selected_true <- matrix(
    NA_integer_,
    nrow = B,
    ncol = n_select
  )

  e_fixed <- rep(NA_real_, B)
  e_same_plugin <- rep(NA_real_, B)
  e_same_true <- rep(NA_real_, B)

  h2_hat <- rep(NA_real_, B)
  converged <- rep(NA_real_, B)
  boundary <- rep(NA_real_, B)

  outer_selected <- candidates[c_outer[candidates] > 0]
  if (length(outer_selected) != n_select) {
    stop("Could not recover the outer selected set from c_outer.")
  }

  for (b in seq_len(B)) {

    yb <- draws[b, seq_len(prep$nobs)]

    fit_b <- fit_reml_animal(
      yb,
      prep
    )

    pred_hat <- predict_ebv_validation(
      fit_b$eta,
      yb,
      prep
    )$uhat

    pred_true <- predict_ebv_validation(
      eta_true,
      yb,
      prep
    )$uhat

    dc <- pred_hat[candidates] - pred_true[candidates]
    delta_candidate[b, ] <- dc

    sp <- rsel_top_selected(
      pred_hat,
      candidates,
      n_select
    )

    st <- rsel_top_selected(
      pred_true,
      candidates,
      n_select
    )

    selected_plugin[b, ] <- sp
    selected_true[b, ] <- st

    e_fixed[b] <- crosssel_policy_contrast(
      dc,
      outer_selected,
      candidates
    )

    e_same_plugin[b] <- crosssel_policy_contrast(
      dc,
      sp,
      candidates
    )

    e_same_true[b] <- crosssel_policy_contrast(
      dc,
      st,
      candidates
    )

    h2_hat[b] <- fit_b$theta[1] / sum(fit_b$theta)
    converged[b] <- as.numeric(fit_b$convergence == 0)
    boundary[b] <- as.numeric(
      fit_b$eta[1] - fit_b$lower[1] < 0.05
    )
  }

  ok <- (
    is.finite(e_fixed) &
      is.finite(e_same_plugin) &
      is.finite(e_same_true) &
      apply(is.finite(delta_candidate), 1L, all)
  )

  if (sum(ok) < max(20L, ceiling(0.80 * B))) {
    stop(
      "Too few valid inner replicates: ",
      sum(ok), " / ", B
    )
  }

  # All standard runs should be fully valid; retaining the filtering makes the
  # code robust if an unusual numerical failure occurs.
  valid_idx <- which(ok)
  dc_mat <- delta_candidate[valid_idx, , drop = FALSE]
  sp_mat <- selected_plugin[valid_idx, , drop = FALSE]
  st_mat <- selected_true[valid_idx, , drop = FALSE]

  ef <- e_fixed[valid_idx]
  esp <- e_same_plugin[valid_idx]
  est <- e_same_true[valid_idx]

  Bv <- length(valid_idx)
  shifts <- crosssel_choose_shifts(Bv, n_shifts = n_cross_shifts)

  # Cross-sample policy contrasts.  Each shift is a derangement.  We pool all
  # Bv x n_shift contrasts when calculating moments.
  cross_plugin <- numeric(Bv * length(shifts))
  cross_true <- numeric(Bv * length(shifts))

  pos <- 1L
  for (sh in shifts) {
    jj <- crosssel_shift_index(Bv, sh)

    for (b in seq_len(Bv)) {
      rng <- pos + b - 1L

      cross_plugin[rng] <- crosssel_policy_contrast(
        dc_mat[b, ],
        sp_mat[jj[b], ],
        candidates
      )

      cross_true[rng] <- crosssel_policy_contrast(
        dc_mat[b, ],
        st_mat[jj[b], ],
        candidates
      )
    }

    pos <- pos + Bv
  }

  # A direct same-sample identity check against the vector c construction.
  same_plugin_vector_check <- numeric(Bv)
  for (b in seq_len(Bv)) {
    ctmp <- make_policy(
      n_animals = n_animals,
      candidates = candidates,
      selected = sp_mat[b, ]
    )

    # Reconstruct candidate-only delta into a full vector only for validation.
    full_delta <- numeric(n_animals)
    full_delta[candidates] <- dc_mat[b, ]
    same_plugin_vector_check[b] <- sum(ctmp * full_delta)
  }

  max_same_plugin_identity_error <- max(
    abs(esp - same_plugin_vector_check),
    na.rm = TRUE
  )

  if (max_same_plugin_identity_error > tolerance) {
    stop(
      "Same-sample plugin contrast identity failed: ",
      format(max_same_plugin_identity_error, digits = 16)
    )
  }

  out <- list(
    n_valid = Bv,
    valid_fraction = Bv / B,
    shifts = shifts,

    bias_fixed = mean(ef),
    bias_same_plugin = mean(esp),
    bias_cross_plugin = mean(cross_plugin),
    bias_same_true = mean(est),
    bias_cross_true = mean(cross_true),

    E2_fixed = mean(ef^2),
    E2_same_plugin = mean(esp^2),
    E2_cross_plugin = mean(cross_plugin^2),
    E2_same_true = mean(est^2),
    E2_cross_true = mean(cross_true^2),

    Var_fixed = var(ef),
    Var_same_plugin = var(esp),
    Var_cross_plugin = var(cross_plugin),
    Var_same_true = var(est),
    Var_cross_true = var(cross_true),

    mean_h2_hat = mean(h2_hat[valid_idx], na.rm = TRUE),
    boundary_rate = mean(boundary[valid_idx], na.rm = TRUE),
    convergence_rate = mean(converged[valid_idx], na.rm = TRUE),

    max_same_plugin_identity_error = max_same_plugin_identity_error
  )

  if (keep_inner_vectors) {
    out$e_fixed <- ef
    out$e_same_plugin <- esp
    out$e_cross_plugin <- cross_plugin
    out$e_same_true <- est
    out$e_cross_true <- cross_true
    out$selected_plugin <- sp_mat
    out$selected_true <- st_mat
    out$delta_candidate <- dc_mat
  }

  out
}

crosssel_summarize_scenario <- function(
    inner_res,
    outer_raw_subset,
    h2,
    information,
    n_pheno,
    reference_row = NULL) {

  grab <- function(name) {
    vapply(
      inner_res,
      function(z) as.numeric(z[[name]]),
      numeric(1)
    )
  }

  outer <- ovib_outer_summary(outer_raw_subset)
  den <- outer$MSE_VC

  mean_E2_fixed <- mean(grab("E2_fixed"))
  mean_E2_same_plugin <- mean(grab("E2_same_plugin"))
  mean_E2_cross_plugin <- mean(grab("E2_cross_plugin"))
  mean_E2_same_true <- mean(grab("E2_same_true"))
  mean_E2_cross_true <- mean(grab("E2_cross_true"))

  R_fixed <- crosssel_safe_ratio(mean_E2_fixed, den)
  R_same_plugin <- crosssel_safe_ratio(mean_E2_same_plugin, den)
  R_cross_plugin <- crosssel_safe_ratio(mean_E2_cross_plugin, den)
  R_same_true <- crosssel_safe_ratio(mean_E2_same_true, den)
  R_cross_true <- crosssel_safe_ratio(mean_E2_cross_true, den)

  ref_fixed <- NA_real_
  ref_same_plugin <- NA_real_
  ref_same_true <- NA_real_

  if (!is.null(reference_row) && nrow(reference_row) == 1L) {
    if ("R_VC_fixed" %in% names(reference_row)) {
      ref_fixed <- as.numeric(reference_row$R_VC_fixed[1])
    }
    if ("R_VC_plugin_selected" %in% names(reference_row)) {
      ref_same_plugin <- as.numeric(reference_row$R_VC_plugin_selected[1])
    }
    if ("R_VC_trueVC_selected" %in% names(reference_row)) {
      ref_same_true <- as.numeric(reference_row$R_VC_trueVC_selected[1])
    }
  }

  data.frame(
    h2 = h2,
    information = information,
    n_pheno = n_pheno,
    n_outer = length(inner_res),
    outer_MSE_VC = den,

    R_VC_fixed = R_fixed,
    R_VC_same_plugin = R_same_plugin,
    R_VC_cross_plugin = R_cross_plugin,
    R_VC_same_trueVC = R_same_true,
    R_VC_cross_trueVC = R_cross_true,

    same_over_cross_plugin = crosssel_safe_ratio(
      mean_E2_same_plugin,
      mean_E2_cross_plugin
    ),
    same_over_cross_trueVC = crosssel_safe_ratio(
      mean_E2_same_true,
      mean_E2_cross_true
    ),
    cross_over_fixed_plugin = crosssel_safe_ratio(
      mean_E2_cross_plugin,
      mean_E2_fixed
    ),
    cross_over_fixed_trueVC = crosssel_safe_ratio(
      mean_E2_cross_true,
      mean_E2_fixed
    ),

    inner_bias_fixed = mean(grab("bias_fixed")),
    inner_bias_same_plugin = mean(grab("bias_same_plugin")),
    inner_bias_cross_plugin = mean(grab("bias_cross_plugin")),
    inner_bias_same_trueVC = mean(grab("bias_same_true")),
    inner_bias_cross_trueVC = mean(grab("bias_cross_true")),

    mean_inner_h2_hat = mean(grab("mean_h2_hat"), na.rm = TRUE),
    mean_boundary_rate = mean(grab("boundary_rate"), na.rm = TRUE),
    mean_convergence_rate = mean(grab("convergence_rate"), na.rm = TRUE),

    max_same_plugin_identity_error = max(
      grab("max_same_plugin_identity_error"),
      na.rm = TRUE
    ),

    reference_R_fixed = ref_fixed,
    reference_R_same_plugin = ref_same_plugin,
    reference_R_same_trueVC = ref_same_true,

    diff_reference_R_fixed = if (is.finite(ref_fixed)) abs(R_fixed - ref_fixed) else NA_real_,
    diff_reference_R_same_plugin = if (is.finite(ref_same_plugin)) abs(R_same_plugin - ref_same_plugin) else NA_real_,
    diff_reference_R_same_trueVC = if (is.finite(ref_same_true)) abs(R_same_true - ref_same_true) else NA_real_,

    row.names = NULL,
    check.names = FALSE
  )
}

run_same_vs_cross_selection_VC_diagnostic <- function(
    boot1000,
    rsel1000 = NULL,
    outer_indices = NULL,
    B_inner = NULL,
    prior_diagnostic_dir = NULL,
    output_dir = NULL,
    checkpoint_every = 10,
    resume = TRUE,
    tolerance = 1e-10,
    n_cross_shifts = 5L,
    keep_inner_vectors = FALSE,
    keep_per_outer_in_master = FALSE) {

  if (is.null(boot1000) || is.null(boot1000$settings)) {
    stop("boot1000 with a $settings element is required.")
  }

  st <- boot1000$settings

  S <- as.integer(st$S)
  B_reference <- as.integer(st$B)

  if (is.null(B_inner)) {
    B <- B_reference
  } else {
    B <- as.integer(B_inner)
  }

  if (!is.finite(B) || B < 20L || B > B_reference) {
    stop(
      "B_inner must be between 20 and the original boot1000 B=",
      B_reference,
      "."
    )
  }

  h2_values <- st$h2_values
  info_fractions <- st$info_fractions
  n_select <- as.integer(st$n_select)
  sigma2_P <- st$sigma2_P
  beta <- st$beta
  seed_info <- st$seed_info
  seed_outer <- st$seed_outer

  if (is.null(outer_indices)) {
    outer_indices <- seq_len(S)
  }

  outer_indices <- sort(unique(as.integer(outer_indices)))

  if (length(outer_indices) == 0L ||
      any(outer_indices < 1L) ||
      any(outer_indices > S)) {
    stop("outer_indices must be within 1:S.")
  }

  full_outer_run <- identical(outer_indices, seq_len(S))

  if (is.null(prior_diagnostic_dir)) {
    prior_diagnostic_dir <- st$output_dir
  }

  if (is.null(prior_diagnostic_dir) ||
      !dir.exists(prior_diagnostic_dir)) {
    stop(
      "Previous bootstrap-generator diagnostic directory not found: ",
      prior_diagnostic_dir
    )
  }

  if (is.null(output_dir)) {
    output_dir <- file.path(
      getwd(),
      "simulation_I_II_final_output",
      paste0(
        "Simulation_II_same_vs_cross_selection_VC_Sinner",
        length(outer_indices),
        "_B",
        B,
        "_K",
        n_cross_shifts
      )
    )
  }

  ovib_dir_create(output_dir)

  base_design <- build_simII_design()

  info_bundle <- build_information_level_designs(
    base_design = base_design,
    info_fractions = info_fractions,
    seed_info = seed_info
  )

  candidates <- base_design$candidates

  cat(
    "\n====================================================\n",
    "SAME-SAMPLE vs CROSS-SAMPLE SELECTION VC DIAGNOSTIC\n",
    "====================================================\n",
    "Outer-replicate count: ", length(outer_indices), " / ", S, "\n",
    "B per outer replicate: ", B, "\n",
    "Cross-policy shifts per outer replicate: ", n_cross_shifts, "\n",
    "Mechanism diagnostic only; NOT a replacement CI.\n",
    sep = ""
  )

  print_information_design(info_bundle)

  set.seed(seed_outer)

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

  scenario_rows <- list()
  scenario_objects <- list()
  row_id <- 1L

  for (ih in seq_along(h2_values)) {

    h2 <- h2_values[ih]

    for (ii in seq_along(info_bundle$info_levels)) {

      info <- info_bundle$info_levels[[ii]]

      key <- paste0(
        "h2_",
        h2,
        "__",
        info$name
      )

      cat(
        "\n====================================================\n",
        "Same-vs-cross diagnostic: ", key,
        " | n_pheno=", info$n_pheno,
        "\n====================================================\n",
        sep = ""
      )

      prior_file <- file.path(
        prior_diagnostic_dir,
        paste0(
          final_safe_name(key),
          "_bootstrap_generator_final.rds"
        )
      )

      if (!file.exists(prior_file)) {
        stop("Previous scenario final RDS not found: ", prior_file)
      }

      prior <- readRDS(prior_file)

      compatible <- (
        !is.null(prior$S) &&
          !is.null(prior$B) &&
          !is.null(prior$key) &&
          identical(as.integer(prior$S), S) &&
          identical(as.integer(prior$B), B_reference) &&
          identical(prior$key, key) &&
          !is.null(prior$raw) &&
          length(prior$raw) == S
      )

      if (!isTRUE(compatible)) {
        stop(
          "Previous scenario RDS is incompatible with this run: ",
          key
        )
      }

      outer_raw_full <- prior$raw
      outer_raw_subset <- outer_raw_full[outer_indices]

      scenario_tag <- paste0(
        ovib_safe_name(key),
        "__nOuter",
        length(outer_indices)
      )

      checkpoint_file <- file.path(
        output_dir,
        paste0(scenario_tag, "_crosssel_checkpoint.rds")
      )

      final_file <- file.path(
        output_dir,
        paste0(scenario_tag, "_crosssel_final.rds")
      )

      if (resume && file.exists(final_file)) {

        done <- readRDS(final_file)

        compatible_final <- (
          !is.null(done$S) &&
            !is.null(done$B) &&
            !is.null(done$key) &&
            !is.null(done$outer_indices) &&
            !is.null(done$n_cross_shifts) &&
            identical(as.integer(done$S), S) &&
            identical(as.integer(done$B), B) &&
            identical(done$key, key) &&
            identical(as.integer(done$outer_indices), outer_indices) &&
            identical(as.integer(done$n_cross_shifts), as.integer(n_cross_shifts))
        )

        if (isTRUE(compatible_final)) {
          cat("Loading completed compatible same-vs-cross scenario.\n")
          scenario_rows[[row_id]] <- done$summary_table
          scenario_objects[[key]] <- done$summary_object
          row_id <- row_id + 1L
          next
        }
      }

      inner_res <- vector("list", length(outer_indices))
      start_pos <- 1L

      if (resume && file.exists(checkpoint_file)) {

        cp <- readRDS(checkpoint_file)

        compatible_cp <- (
          !is.null(cp$S) &&
            !is.null(cp$B) &&
            !is.null(cp$key) &&
            !is.null(cp$outer_indices) &&
            !is.null(cp$n_cross_shifts) &&
            identical(as.integer(cp$S), S) &&
            identical(as.integer(cp$B), B) &&
            identical(cp$key, key) &&
            identical(as.integer(cp$outer_indices), outer_indices) &&
            identical(as.integer(cp$n_cross_shifts), as.integer(n_cross_shifts))
        )

        if (isTRUE(compatible_cp)) {
          inner_res <- cp$inner_res
          missing_pos <- which(vapply(inner_res, is.null, logical(1)))

          if (length(missing_pos) == 0L) {
            start_pos <- length(outer_indices) + 1L
          } else {
            start_pos <- min(missing_pos)
          }

          cat(
            "Resuming at subset position ", start_pos, ".\n",
            sep = ""
          )
        }
      }

      if (start_pos <= length(outer_indices)) {

        for (pos in seq.int(start_pos, length(outer_indices))) {

          s <- outer_indices[pos]

          dat <- make_dataset_from_latent(
            latent = latent_list[[s]],
            info_level = info,
            sigma2_P = sigma2_P,
            h2 = h2,
            beta = beta
          )

          zouter <- outer_raw_full[[s]]

          c_outer <- as.numeric(zouter$cstar)
          truth_rebuilt <- sum(c_outer * dat$u)
          truth_diff <- abs(
            truth_rebuilt - as.numeric(zouter$truth)
          )

          if (truth_diff > tolerance) {
            stop(
              "Outer truth mismatch at ", key,
              ", replicate ", s,
              ": ", format(truth_diff, digits = 16)
            )
          }

          inner_res[[pos]] <- crosssel_one_outer(
            outer_raw = zouter,
            dat = dat,
            prep = info$prep,
            candidates = candidates,
            n_select = n_select,
            B = B,
            B_reference = B_reference,
            seed_boot = seed_methods[s, ih, ii, 2],
            n_cross_shifts = n_cross_shifts,
            tolerance = tolerance,
            keep_inner_vectors = keep_inner_vectors
          )

          if (
            pos %% checkpoint_every == 0L ||
            pos == length(outer_indices)
          ) {
            saveRDS(
              list(
                S = S,
                B = B,
                key = key,
                outer_indices = outer_indices,
                n_cross_shifts = n_cross_shifts,
                inner_res = inner_res
              ),
              checkpoint_file
            )

            cat(
              "Completed subset position ", pos,
              " / ", length(outer_indices),
              " (outer replicate ", s, ")\n",
              sep = ""
            )
          }
        }
      }

      if (any(vapply(inner_res, is.null, logical(1)))) {
        stop("Scenario contains incomplete inner results: ", key)
      }

      reference_row <- NULL

      if (
        full_outer_run &&
        identical(B, B_reference) &&
        !is.null(rsel1000) &&
        !is.null(rsel1000$key_table)
      ) {
        rr <- rsel1000$key_table
        idx <- which(
          abs(rr$h2 - h2) < 1e-12 &
            as.character(rr$information) == info$name
        )
        if (length(idx) == 1L) {
          reference_row <- rr[idx, , drop = FALSE]
        }
      }

      summary_table <- crosssel_summarize_scenario(
        inner_res = inner_res,
        outer_raw_subset = outer_raw_subset,
        h2 = h2,
        information = info$name,
        n_pheno = info$n_pheno,
        reference_row = reference_row
      )

      validation <- data.frame(
        check = c(
          "same_plugin_vector_identity",
          "previous_R_fixed_reproduction",
          "previous_R_same_plugin_reproduction",
          "previous_R_same_trueVC_reproduction"
        ),
        max_abs_error = c(
          summary_table$max_same_plugin_identity_error,
          summary_table$diff_reference_R_fixed,
          summary_table$diff_reference_R_same_plugin,
          summary_table$diff_reference_R_same_trueVC
        ),
        within_tolerance = c(
          summary_table$max_same_plugin_identity_error <= tolerance,
          if (is.finite(summary_table$diff_reference_R_fixed)) {
            summary_table$diff_reference_R_fixed <= max(tolerance, 1e-8)
          } else NA,
          if (is.finite(summary_table$diff_reference_R_same_plugin)) {
            summary_table$diff_reference_R_same_plugin <= max(tolerance, 1e-8)
          } else NA,
          if (is.finite(summary_table$diff_reference_R_same_trueVC)) {
            summary_table$diff_reference_R_same_trueVC <= max(tolerance, 1e-8)
          } else NA
        ),
        row.names = NULL
      )

      if (!isTRUE(validation$within_tolerance[1])) {
        stop(
          "Structural same-vs-cross validation failed for ", key,
          "."
        )
      }

      # Full production runs supplied with rsel1000 should reproduce all three
      # previous ratios exactly apart from floating-point rounding.
      ref_rows <- which(!is.na(validation$within_tolerance[-1L])) + 1L
      if (length(ref_rows) > 0L &&
          any(!validation$within_tolerance[ref_rows])) {
        warning(
          "One or more previous reselection ratios were not reproduced for ",
          key,
          ". Inspect validation before interpretation."
        )
      }

      per_outer <- data.frame(
        outer_replicate = outer_indices,
        E2_fixed = vapply(inner_res, `[[`, numeric(1), "E2_fixed"),
        E2_same_plugin = vapply(inner_res, `[[`, numeric(1), "E2_same_plugin"),
        E2_cross_plugin = vapply(inner_res, `[[`, numeric(1), "E2_cross_plugin"),
        E2_same_true = vapply(inner_res, `[[`, numeric(1), "E2_same_true"),
        E2_cross_true = vapply(inner_res, `[[`, numeric(1), "E2_cross_true"),
        bias_fixed = vapply(inner_res, `[[`, numeric(1), "bias_fixed"),
        bias_same_plugin = vapply(inner_res, `[[`, numeric(1), "bias_same_plugin"),
        bias_cross_plugin = vapply(inner_res, `[[`, numeric(1), "bias_cross_plugin"),
        bias_same_true = vapply(inner_res, `[[`, numeric(1), "bias_same_true"),
        bias_cross_true = vapply(inner_res, `[[`, numeric(1), "bias_cross_true"),
        row.names = NULL,
        check.names = FALSE
      )

      summary_object <- list(
        h2 = h2,
        information = info$name,
        n_pheno = info$n_pheno,
        outer_indices = outer_indices,
        summary_table = summary_table,
        validation = validation
      )

      saveRDS(
        list(
          S = S,
          B = B,
          key = key,
          outer_indices = outer_indices,
          n_cross_shifts = n_cross_shifts,
          summary_table = summary_table,
          summary_object = summary_object,
          validation = validation,
          per_outer = per_outer,
          inner_res = if (keep_inner_vectors) inner_res else NULL
        ),
        final_file
      )

      write.csv(
        summary_table,
        file.path(output_dir, paste0(scenario_tag, "_summary.csv")),
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

      scenario_rows[[row_id]] <- summary_table
      scenario_objects[[key]] <- summary_object

      if (keep_per_outer_in_master) {
        scenario_objects[[key]]$per_outer <- per_outer
      }

      row_id <- row_id + 1L
    }
  }

  key_table <- do.call(rbind, scenario_rows)
  rownames(key_table) <- NULL

  validation_all <- do.call(
    rbind,
    lapply(
      names(scenario_objects),
      function(k) {
        vv <- scenario_objects[[k]]$validation
        vv$scenario <- k
        vv[, c("scenario", "check", "max_abs_error", "within_tolerance")]
      }
    )
  )
  rownames(validation_all) <- NULL

  out <- list(
    settings = list(
      S_reference = S,
      B_reference = B_reference,
      outer_indices = outer_indices,
      B_inner = B,
      n_cross_shifts = n_cross_shifts,
      h2_values = h2_values,
      info_fractions = info_fractions,
      n_select = n_select,
      estimand = paste0(
        "mechanism diagnostic: same-sample versus independent-cross-sample ",
        "policy construction for VC-induced EBLUP perturbation"
      ),
      output_dir = normalizePath(
        output_dir,
        winslash = "/",
        mustWork = FALSE
      )
    ),
    key_table = key_table,
    validation = validation_all,
    result = scenario_objects
  )

  write.csv(
    key_table,
    file.path(output_dir, "same_vs_cross_selection_VC_key_table.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    validation_all,
    file.path(output_dir, "same_vs_cross_selection_VC_validation.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  saveRDS(
    out,
    file.path(output_dir, "same_vs_cross_selection_VC_master.rds")
  )

  cat(
    "\n====================================================\n",
    "Same-vs-cross selection VC diagnostic completed.\n",
    "====================================================\n",
    sep = ""
  )

  cat("\nKey table:\n")
  print(key_table)

  cat("\nValidation:\n")
  print(validation_all)

  invisible(out)
}

print_same_vs_cross_selection_VC_diagnostic <- function(x, digits = 4) {

  if (is.null(x$key_table)) stop("Object has no $key_table.")

  z <- x$key_table

  cat("\nVC second-moment reproduction: fixed vs same-sample vs cross-sample\n")
  print(
    round(
      z[, c(
        "h2",
        "n_pheno",
        "R_VC_fixed",
        "R_VC_same_plugin",
        "R_VC_cross_plugin",
        "R_VC_same_trueVC",
        "R_VC_cross_trueVC",
        "same_over_cross_plugin"
      )],
      digits
    )
  )

  cat("\nBias diagnostic\n")
  print(
    round(
      z[, c(
        "h2",
        "n_pheno",
        "inner_bias_fixed",
        "inner_bias_same_plugin",
        "inner_bias_cross_plugin",
        "inner_bias_same_trueVC",
        "inner_bias_cross_trueVC"
      )],
      digits
    )
  )

  invisible(z)
}

# ======================================================================
# RECOMMENDED USE
# ======================================================================
#
# 1) PILOT: 20 outer x 100 inner, five independent cyclic cross-pairings
#
#   cross20 <- run_same_vs_cross_selection_VC_diagnostic(
#     boot1000 = boot1000,
#     rsel1000 = rsel1000,
#     outer_indices = 1:20,
#     B_inner = 100,
#     n_cross_shifts = 5,
#     checkpoint_every = 5,
#     resume = TRUE
#   )
#
#   print_same_vs_cross_selection_VC_diagnostic(cross20)
#
# 2) FULL PRODUCTION: 1000 outer x 500 inner
#
#   cross1000 <- run_same_vs_cross_selection_VC_diagnostic(
#     boot1000 = boot1000,
#     rsel1000 = rsel1000,
#     n_cross_shifts = 5,
#     checkpoint_every = 10,
#     resume = TRUE
#   )
#
#   print(cross1000$validation)
#   print(cross1000$key_table)
#   print_same_vs_cross_selection_VC_diagnostic(cross1000)
#
# PRIMARY INTERPRETATION
# ----------------------
# - R_VC_same_plugin should reproduce the previous ~0.9-1.05 result.
# - If R_VC_cross_plugin falls back near R_VC_fixed (~0.11-0.12), then
#   same-sample selection/error dependence is the key source of amplification.
# - If R_VC_cross_plugin remains near R_VC_same_plugin, then marginal variation
#   in selected policies, rather than same-sample alignment, is sufficient.
# - At low h2, compare inner_bias_same_plugin with inner_bias_cross_plugin.
#   Collapse of the positive bias under cross-pairing would directly show that
#   the bias is selection-induced alignment rather than a marginal VC effect.
# ======================================================================


# ======================================================================
# STEP 5
# TWO-LAYER SELECTION-VC ALIGNMENT DECOMPOSITION
# ======================================================================
#
# PURPOSE
# -------
# Decompose the SAME-SAMPLE plugin-selected VC-induced error
#
#   e_plugin = c_plugin' d
#
# where
#
#   d = uhat(hat(theta)) - uhat(theta_true),
#
# into two layers:
#
#   e_plugin
#     = e_common + e_switch
#
#   e_common = c_trueVC' d
#
#   e_switch = (c_plugin - c_trueVC)' d
#            = e_plugin - e_common
#            >= 0  (for exact top-k selection, apart from numerical tolerance).
#
# Interpretation:
#
#   e_common:
#     common-data alignment that remains even when the selected set is formed
#     from TRUE-VC EBVs.  c_trueVC and d still come from the SAME inner data.
#
#   e_switch:
#     additional increment caused by changing the selected top-k set from the
#     TRUE-VC EBV ranking to the plugin REML-EBLUP ranking.
#
# SECOND-MOMENT DECOMPOSITION
# ---------------------------
#
# Exact within each scenario:
#
#   E(e_plugin^2)
#     = E(e_common^2)
#       + E(e_switch^2)
#       + 2 E(e_common * e_switch)
#
# A second, hierarchical decomposition uses an INDEPENDENT cross-sample
# TRUE-VC-selected policy as the baseline:
#
#   M_plugin
#     = M_cross_true
#       + (M_common - M_cross_true)
#       + (M_plugin - M_common)
#
# where
#
#   baseline independent-policy component:
#       M_cross_true
#
#   common-data alignment excess:
#       M_common - M_cross_true
#
#   policy-switch second-moment increment:
#       M_plugin - M_common
#
# These three terms sum EXACTLY to M_plugin.
#
# IMPORTANT
# ---------
# This is a mechanism diagnostic, not a replacement interval method.
# It preserves the same TRUE-VC inner generator and the exact same RNG mapping
# used by the preceding diagnostics so that completed production results can be
# reproduced numerically.
#
# ONE-FILE USE
# ------------
# This file contains the complete preceding same-vs-cross diagnostic code plus
# this new STEP 5 extension.  No other R file needs to be sourced.
#
# PILOT:
#
#
#   layer20 <- run_two_layer_alignment_decomposition(
#     boot1000 = boot1000,
#     cross_reference = cross1000,
#     outer_indices = 1:20,
#     B_inner = 100,
#     checkpoint_every = 5,
#     resume = FALSE
#   )
#
#   print_two_layer_alignment_decomposition(layer20)
#
# PRODUCTION:
#
#   layer1000 <- run_two_layer_alignment_decomposition(
#     boot1000 = boot1000,
#     cross_reference = cross1000,
#     checkpoint_every = 10,
#     resume = TRUE
#   )
#
#   print(layer1000$validation)
#   print(layer1000$key_table)
#   print_two_layer_alignment_decomposition(layer1000)
#
# ======================================================================


twolayer_safe_ratio <- function(a, b) {
  ifelse(
    is.finite(a) & is.finite(b) & b != 0,
    a / b,
    NA_real_
  )
}


twolayer_exact_set_match_rows <- function(a, b) {

  if (!all(dim(a) == dim(b))) {
    stop("Selection matrices must have identical dimensions.")
  }

  vapply(
    seq_len(nrow(a)),
    function(i) {
      identical(
        sort(as.integer(a[i, ])),
        sort(as.integer(b[i, ]))
      )
    },
    logical(1)
  )
}


twolayer_overlap_rows <- function(a, b) {

  if (!all(dim(a) == dim(b))) {
    stop("Selection matrices must have identical dimensions.")
  }

  k <- ncol(a)

  vapply(
    seq_len(nrow(a)),
    function(i) {
      length(
        intersect(
          as.integer(a[i, ]),
          as.integer(b[i, ])
        )
      ) / k
    },
    numeric(1)
  )
}


twolayer_summarize_inner <- function(
    z,
    tolerance = 1e-10) {

  required <- c(
    "e_same_plugin",
    "e_same_true",
    "e_cross_true",
    "e_fixed",
    "selected_plugin",
    "selected_true",
    "delta_candidate"
  )

  miss <- setdiff(required, names(z))
  if (length(miss) > 0L) {
    stop(
      "Inner vectors are missing from crosssel_one_outer result: ",
      paste(miss, collapse = ", ")
    )
  }

  e_plugin <- as.numeric(z$e_same_plugin)
  e_common <- as.numeric(z$e_same_true)
  e_cross_true <- as.numeric(z$e_cross_true)
  e_fixed <- as.numeric(z$e_fixed)

  if (length(e_plugin) != length(e_common)) {
    stop("e_plugin and e_common lengths differ.")
  }

  e_switch <- e_plugin - e_common

  sp <- z$selected_plugin
  st <- z$selected_true

  exact_match <- twolayer_exact_set_match_rows(sp, st)
  overlap <- twolayer_overlap_rows(sp, st)

  # ------------------------------------------------------------
  # Algebraic identities
  # ------------------------------------------------------------

  max_error_identity <- max(
    abs(e_plugin - (e_common + e_switch)),
    na.rm = TRUE
  )

  M_plugin <- mean(e_plugin^2)
  M_common <- mean(e_common^2)
  M_switch_sq <- mean(e_switch^2)
  two_cross <- 2 * mean(e_common * e_switch)

  second_moment_identity_error <- abs(
    M_plugin -
      (M_common + M_switch_sq + two_cross)
  )

  # Deterministic top-k result:
  # (c_plugin - c_trueVC)'d >= 0.
  min_switch <- min(e_switch, na.rm = TRUE)
  max_negative_switch <- max(
    pmax(-e_switch, 0),
    na.rm = TRUE
  )

  nonnegative_switch_rate <- mean(
    e_switch >= -tolerance,
    na.rm = TRUE
  )

  # When the selected sets are exactly the same, switch must be zero.
  max_abs_switch_when_same <- if (any(exact_match)) {
    max(abs(e_switch[exact_match]), na.rm = TRUE)
  } else {
    NA_real_
  }

  changed <- !exact_match

  mean_switch_when_changed <- if (any(changed)) {
    mean(e_switch[changed], na.rm = TRUE)
  } else {
    0
  }

  median_switch_when_changed <- if (any(changed)) {
    median(e_switch[changed], na.rm = TRUE)
  } else {
    0
  }

  q95_switch_when_changed <- if (any(changed)) {
    unname(
      quantile(
        e_switch[changed],
        probs = 0.95,
        names = FALSE,
        type = 8,
        na.rm = TRUE
      )
    )
  } else {
    0
  }

  M_cross_true <- mean(e_cross_true^2)
  M_fixed <- mean(e_fixed^2)

  bias_plugin <- mean(e_plugin)
  bias_common <- mean(e_common)
  bias_switch <- mean(e_switch)
  bias_cross_true <- mean(e_cross_true)
  bias_fixed <- mean(e_fixed)

  # ------------------------------------------------------------
  # Hierarchical second-moment decomposition:
  #
  # M_plugin
  #   = M_cross_true
  #   + (M_common - M_cross_true)
  #   + (M_plugin - M_common)
  # ------------------------------------------------------------

  M_baseline <- M_cross_true
  M_common_alignment_excess <- M_common - M_cross_true
  M_policy_switch_increment <- M_plugin - M_common

  hierarchy_identity_error <- abs(
    M_plugin -
      (
        M_baseline +
          M_common_alignment_excess +
          M_policy_switch_increment
      )
  )

  # Bias analogue:
  #
  # E(e_plugin)
  #   = E(e_cross_true)
  #   + [E(e_common)-E(e_cross_true)]
  #   + E(e_switch)
  bias_common_alignment <- bias_common - bias_cross_true

  bias_hierarchy_identity_error <- abs(
    bias_plugin -
      (
        bias_cross_true +
          bias_common_alignment +
          bias_switch
      )
  )

  out <- list(
    n_valid = length(e_plugin),

    bias_plugin = bias_plugin,
    bias_common = bias_common,
    bias_switch = bias_switch,
    bias_cross_true = bias_cross_true,
    bias_fixed = bias_fixed,
    bias_common_alignment = bias_common_alignment,

    M_plugin = M_plugin,
    M_common = M_common,
    M_switch_sq = M_switch_sq,
    two_common_switch = two_cross,
    M_cross_true = M_cross_true,
    M_fixed = M_fixed,

    M_baseline = M_baseline,
    M_common_alignment_excess = M_common_alignment_excess,
    M_policy_switch_increment = M_policy_switch_increment,

    baseline_fraction_of_plugin =
      twolayer_safe_ratio(M_baseline, M_plugin),

    common_alignment_fraction_of_plugin =
      twolayer_safe_ratio(
        M_common_alignment_excess,
        M_plugin
      ),

    policy_switch_fraction_of_plugin =
      twolayer_safe_ratio(
        M_policy_switch_increment,
        M_plugin
      ),

    common_raw_fraction_of_plugin =
      twolayer_safe_ratio(M_common, M_plugin),

    switch_sq_fraction_of_plugin =
      twolayer_safe_ratio(M_switch_sq, M_plugin),

    two_common_switch_fraction_of_plugin =
      twolayer_safe_ratio(two_cross, M_plugin),

    cross_true_over_fixed =
      twolayer_safe_ratio(M_cross_true, M_fixed),

    exact_set_match_rate = mean(exact_match),
    mean_overlap = mean(overlap),
    mean_number_replaced =
      ncol(sp) * (1 - mean(overlap)),

    mean_switch_when_changed = mean_switch_when_changed,
    median_switch_when_changed = median_switch_when_changed,
    q95_switch_when_changed = q95_switch_when_changed,

    min_switch = min_switch,
    max_negative_switch = max_negative_switch,
    nonnegative_switch_rate = nonnegative_switch_rate,

    max_abs_switch_when_same = max_abs_switch_when_same,
    max_error_identity = max_error_identity,
    second_moment_identity_error = second_moment_identity_error,
    hierarchy_identity_error = hierarchy_identity_error,
    bias_hierarchy_identity_error = bias_hierarchy_identity_error
  )

  out
}


twolayer_summarize_scenario <- function(
    per_outer,
    outer_MSE_VC,
    h2,
    information,
    n_pheno,
    reference_row = NULL,
    tolerance = 1e-10) {

  grab <- function(name) {
    vapply(
      per_outer,
      function(z) as.numeric(z[[name]]),
      numeric(1)
    )
  }

  mean_M_plugin <- mean(grab("M_plugin"))
  mean_M_common <- mean(grab("M_common"))
  mean_M_switch_sq <- mean(grab("M_switch_sq"))
  mean_two_cross <- mean(grab("two_common_switch"))
  mean_M_cross_true <- mean(grab("M_cross_true"))
  mean_M_fixed <- mean(grab("M_fixed"))

  M_baseline <- mean_M_cross_true
  M_common_alignment_excess <-
    mean_M_common - mean_M_cross_true
  M_policy_switch_increment <-
    mean_M_plugin - mean_M_common

  hierarchy_sum <- (
    M_baseline +
      M_common_alignment_excess +
      M_policy_switch_increment
  )

  bias_plugin <- mean(grab("bias_plugin"))
  bias_common <- mean(grab("bias_common"))
  bias_switch <- mean(grab("bias_switch"))
  bias_cross_true <- mean(grab("bias_cross_true"))
  bias_common_alignment <- (
    bias_common - bias_cross_true
  )

  # Exact second-moment expansion on the SAME-sample scale.
  second_moment_expansion_error <- abs(
    mean_M_plugin -
      (
        mean_M_common +
          mean_M_switch_sq +
          mean_two_cross
      )
  )

  hierarchy_error <- abs(
    mean_M_plugin - hierarchy_sum
  )

  bias_hierarchy_error <- abs(
    bias_plugin -
      (
        bias_cross_true +
          bias_common_alignment +
          bias_switch
      )
  )

  ref_same_plugin <- NA_real_
  ref_same_true <- NA_real_
  ref_cross_true <- NA_real_
  ref_fixed <- NA_real_

  if (!is.null(reference_row) && nrow(reference_row) == 1L) {
    if ("R_VC_same_plugin" %in% names(reference_row)) {
      ref_same_plugin <- as.numeric(
        reference_row$R_VC_same_plugin[1]
      )
    }
    if ("R_VC_same_trueVC" %in% names(reference_row)) {
      ref_same_true <- as.numeric(
        reference_row$R_VC_same_trueVC[1]
      )
    }
    if ("R_VC_cross_trueVC" %in% names(reference_row)) {
      ref_cross_true <- as.numeric(
        reference_row$R_VC_cross_trueVC[1]
      )
    }
    if ("R_VC_fixed" %in% names(reference_row)) {
      ref_fixed <- as.numeric(
        reference_row$R_VC_fixed[1]
      )
    }
  }

  R_plugin <- twolayer_safe_ratio(
    mean_M_plugin,
    outer_MSE_VC
  )

  R_common <- twolayer_safe_ratio(
    mean_M_common,
    outer_MSE_VC
  )

  R_cross_true <- twolayer_safe_ratio(
    mean_M_cross_true,
    outer_MSE_VC
  )

  R_fixed <- twolayer_safe_ratio(
    mean_M_fixed,
    outer_MSE_VC
  )

  out <- data.frame(
    h2 = h2,
    information = information,
    n_pheno = n_pheno,
    n_outer = length(per_outer),
    outer_MSE_VC = outer_MSE_VC,

    # Core moments
    M_plugin = mean_M_plugin,
    M_common = mean_M_common,
    M_switch_sq = mean_M_switch_sq,
    two_common_switch = mean_two_cross,
    M_cross_true = mean_M_cross_true,
    M_fixed = mean_M_fixed,

    # Exact hierarchical decomposition of M_plugin
    M_baseline_independent = M_baseline,
    M_common_alignment_excess =
      M_common_alignment_excess,
    M_policy_switch_increment =
      M_policy_switch_increment,

    baseline_fraction_of_plugin =
      twolayer_safe_ratio(
        M_baseline,
        mean_M_plugin
      ),

    common_alignment_fraction_of_plugin =
      twolayer_safe_ratio(
        M_common_alignment_excess,
        mean_M_plugin
      ),

    policy_switch_fraction_of_plugin =
      twolayer_safe_ratio(
        M_policy_switch_increment,
        mean_M_plugin
      ),

    # Exact same-sample square expansion
    common_raw_fraction_of_plugin =
      twolayer_safe_ratio(
        mean_M_common,
        mean_M_plugin
      ),

    switch_sq_fraction_of_plugin =
      twolayer_safe_ratio(
        mean_M_switch_sq,
        mean_M_plugin
      ),

    two_common_switch_fraction_of_plugin =
      twolayer_safe_ratio(
        mean_two_cross,
        mean_M_plugin
      ),

    # Reproduction ratios
    R_VC_plugin = R_plugin,
    R_VC_common_trueVC_selected = R_common,
    R_VC_cross_trueVC = R_cross_true,
    R_VC_fixed = R_fixed,

    # Bias decomposition
    bias_plugin = bias_plugin,
    bias_common = bias_common,
    bias_switch = bias_switch,
    bias_cross_true = bias_cross_true,
    bias_common_alignment =
      bias_common_alignment,

    # Selection diagnostics
    exact_set_match_rate =
      mean(grab("exact_set_match_rate")),
    mean_overlap =
      mean(grab("mean_overlap")),
    mean_number_replaced =
      mean(grab("mean_number_replaced")),

    mean_switch_when_changed =
      mean(grab("mean_switch_when_changed")),
    median_switch_when_changed =
      mean(grab("median_switch_when_changed")),
    q95_switch_when_changed =
      mean(grab("q95_switch_when_changed")),

    # Deterministic/nonnegative switch check
    min_switch =
      min(grab("min_switch")),
    max_negative_switch =
      max(grab("max_negative_switch")),
    nonnegative_switch_rate =
      mean(grab("nonnegative_switch_rate")),

    max_abs_switch_when_same =
      max(
        grab("max_abs_switch_when_same"),
        na.rm = TRUE
      ),

    # Algebraic checks
    max_inner_error_identity =
      max(grab("max_error_identity")),
    max_inner_second_moment_identity_error =
      max(grab("second_moment_identity_error")),
    max_inner_hierarchy_identity_error =
      max(grab("hierarchy_identity_error")),
    max_inner_bias_hierarchy_identity_error =
      max(grab("bias_hierarchy_identity_error")),

    scenario_second_moment_expansion_error =
      second_moment_expansion_error,
    scenario_hierarchy_identity_error =
      hierarchy_error,
    scenario_bias_hierarchy_identity_error =
      bias_hierarchy_error,

    # References to preceding full same-vs-cross run
    reference_R_VC_same_plugin =
      ref_same_plugin,
    reference_R_VC_same_trueVC =
      ref_same_true,
    reference_R_VC_cross_trueVC =
      ref_cross_true,
    reference_R_VC_fixed =
      ref_fixed,

    diff_reference_R_VC_same_plugin =
      if (is.finite(ref_same_plugin)) {
        abs(R_plugin - ref_same_plugin)
      } else {
        NA_real_
      },

    diff_reference_R_VC_same_trueVC =
      if (is.finite(ref_same_true)) {
        abs(R_common - ref_same_true)
      } else {
        NA_real_
      },

    diff_reference_R_VC_cross_trueVC =
      if (is.finite(ref_cross_true)) {
        abs(R_cross_true - ref_cross_true)
      } else {
        NA_real_
      },

    diff_reference_R_VC_fixed =
      if (is.finite(ref_fixed)) {
        abs(R_fixed - ref_fixed)
      } else {
        NA_real_
      },

    row.names = NULL,
    check.names = FALSE
  )

  out
}


run_two_layer_alignment_decomposition <- function(
    boot1000,
    cross_reference = NULL,
    outer_indices = NULL,
    B_inner = NULL,
    prior_diagnostic_dir = NULL,
    output_dir = NULL,
    checkpoint_every = 10,
    resume = TRUE,
    tolerance = 1e-10,
    n_cross_shifts = 5L,
    keep_per_outer_in_master = FALSE) {

  if (is.null(boot1000) || is.null(boot1000$settings)) {
    stop("boot1000 with a $settings element is required.")
  }

  st <- boot1000$settings

  S <- as.integer(st$S)
  B_reference <- as.integer(st$B)

  if (is.null(B_inner)) {
    B <- B_reference
  } else {
    B <- as.integer(B_inner)
  }

  if (!is.finite(B) || B < 20L || B > B_reference) {
    stop(
      "B_inner must be between 20 and original boot1000 B=",
      B_reference,
      "."
    )
  }

  h2_values <- st$h2_values
  info_fractions <- st$info_fractions
  n_select <- as.integer(st$n_select)
  sigma2_P <- st$sigma2_P
  beta <- st$beta
  seed_info <- st$seed_info
  seed_outer <- st$seed_outer

  if (is.null(outer_indices)) {
    outer_indices <- seq_len(S)
  }

  outer_indices <- sort(unique(as.integer(outer_indices)))

  if (
    length(outer_indices) == 0L ||
      any(outer_indices < 1L) ||
      any(outer_indices > S)
  ) {
    stop("outer_indices must be within 1:S.")
  }

  full_reference_run <- (
    identical(outer_indices, seq_len(S)) &&
      identical(B, B_reference)
  )

  if (is.null(prior_diagnostic_dir)) {
    prior_diagnostic_dir <- st$output_dir
  }

  if (
    is.null(prior_diagnostic_dir) ||
      !dir.exists(prior_diagnostic_dir)
  ) {
    stop(
      "Previous bootstrap-generator diagnostic directory not found: ",
      prior_diagnostic_dir
    )
  }

  if (is.null(output_dir)) {
    output_dir <- file.path(
      getwd(),
      "simulation_I_II_final_output",
      paste0(
        "Simulation_II_two_layer_alignment_Sinner",
        length(outer_indices),
        "_B",
        B,
        "_K",
        n_select
      )
    )
  }

  ovib_dir_create(output_dir)

  base_design <- build_simII_design()

  info_bundle <- build_information_level_designs(
    base_design = base_design,
    info_fractions = info_fractions,
    seed_info = seed_info
  )

  candidates <- base_design$candidates

  cat(
    "\n====================================================\n",
    "TWO-LAYER SELECTION-VC ALIGNMENT DECOMPOSITION\n",
    "====================================================\n",
    "Outer-replicate count: ",
    length(outer_indices),
    " / ",
    S,
    "\n",
    "B per outer replicate: ",
    B,
    "\n",
    "Selected candidates: ",
    n_select,
    "\n",
    "Mechanism diagnostic only; NOT a replacement CI.\n",
    sep = ""
  )

  print_information_design(info_bundle)

  # ------------------------------------------------------------
  # Rebuild EXACT same outer latent draws and method seeds as the
  # preceding production simulations.
  # ------------------------------------------------------------

  set.seed(seed_outer)

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

  scenario_rows <- list()
  scenario_objects <- list()
  validation_rows <- list()
  row_id <- 1L

  for (ih in seq_along(h2_values)) {

    h2 <- h2_values[ih]

    for (ii in seq_along(info_bundle$info_levels)) {

      info <- info_bundle$info_levels[[ii]]

      key <- paste0(
        "h2_",
        h2,
        "__",
        info$name
      )

      cat(
        "\n====================================================\n",
        "Two-layer decomposition: ",
        key,
        " | n_pheno=",
        info$n_pheno,
        "\n====================================================\n",
        sep = ""
      )

      prior_file <- file.path(
        prior_diagnostic_dir,
        paste0(
          final_safe_name(key),
          "_bootstrap_generator_final.rds"
        )
      )

      if (!file.exists(prior_file)) {
        stop(
          "Previous scenario final RDS not found: ",
          prior_file
        )
      }

      prior <- readRDS(prior_file)

      compatible <- (
        !is.null(prior$S) &&
          !is.null(prior$B) &&
          !is.null(prior$key) &&
          identical(as.integer(prior$S), S) &&
          identical(as.integer(prior$B), B_reference) &&
          identical(prior$key, key) &&
          !is.null(prior$raw) &&
          length(prior$raw) == S
      )

      if (!isTRUE(compatible)) {
        stop(
          "Previous scenario RDS is incompatible: ",
          key
        )
      }

      outer_raw_full <- prior$raw
      outer_raw_subset <- outer_raw_full[outer_indices]

      # Outer denominator from the exact previous outer results.
      outer_summary <- ovib_outer_summary(
        outer_raw_subset
      )

      outer_MSE_VC <- as.numeric(
        outer_summary$MSE_VC
      )

      scenario_tag <- paste0(
        ovib_safe_name(key),
        "__nOuter",
        length(outer_indices)
      )

      checkpoint_file <- file.path(
        output_dir,
        paste0(
          scenario_tag,
          "_twolayer_checkpoint.rds"
        )
      )

      final_file <- file.path(
        output_dir,
        paste0(
          scenario_tag,
          "_twolayer_final.rds"
        )
      )

      if (resume && file.exists(final_file)) {

        done <- readRDS(final_file)

        compatible_final <- (
          !is.null(done$S) &&
            !is.null(done$B) &&
            !is.null(done$key) &&
            !is.null(done$outer_indices) &&
            identical(as.integer(done$S), S) &&
            identical(as.integer(done$B), B) &&
            identical(done$key, key) &&
            identical(
              as.integer(done$outer_indices),
              outer_indices
            )
        )

        if (isTRUE(compatible_final)) {

          cat(
            "Loading completed compatible two-layer scenario.\n"
          )

          scenario_rows[[row_id]] <-
            done$summary_table

          scenario_objects[[key]] <-
            done$summary_object

          validation_rows[[row_id]] <-
            done$validation

          row_id <- row_id + 1L
          next
        }
      }

      per_outer <- vector(
        "list",
        length(outer_indices)
      )

      start_pos <- 1L

      if (
        resume &&
          file.exists(checkpoint_file)
      ) {

        cp <- readRDS(checkpoint_file)

        compatible_cp <- (
          !is.null(cp$S) &&
            !is.null(cp$B) &&
            !is.null(cp$key) &&
            !is.null(cp$outer_indices) &&
            identical(as.integer(cp$S), S) &&
            identical(as.integer(cp$B), B) &&
            identical(cp$key, key) &&
            identical(
              as.integer(cp$outer_indices),
              outer_indices
            )
        )

        if (isTRUE(compatible_cp)) {

          per_outer <- cp$per_outer

          missing_pos <- which(
            vapply(
              per_outer,
              is.null,
              logical(1)
            )
          )

          if (length(missing_pos) == 0L) {
            start_pos <-
              length(outer_indices) + 1L
          } else {
            start_pos <- min(missing_pos)
          }

          cat(
            "Resuming at subset position ",
            start_pos,
            ".\n",
            sep = ""
          )
        }
      }

      if (start_pos <= length(outer_indices)) {

        for (
          pos in seq.int(
            start_pos,
            length(outer_indices)
          )
        ) {

          s <- outer_indices[pos]

          dat <- make_dataset_from_latent(
            latent = latent_list[[s]],
            info_level = info,
            sigma2_P = sigma2_P,
            h2 = h2,
            beta = beta
          )

          zouter <- outer_raw_full[[s]]

          # Exact outer data reconstruction check.
          c_outer <- as.numeric(
            zouter$cstar
          )

          truth_rebuilt <- sum(
            c_outer * dat$u
          )

          truth_diff <- abs(
            truth_rebuilt -
              as.numeric(zouter$truth)
          )

          if (truth_diff > tolerance) {
            stop(
              "Outer truth mismatch at ",
              key,
              ", replicate ",
              s,
              ": ",
              format(
                truth_diff,
                digits = 16
              )
            )
          }

          # Reuse the SAME true-VC inner generator and the SAME seed
          # as the prior same-vs-cross/reselection diagnostics.
          zin <- crosssel_one_outer(
            outer_raw = zouter,
            dat = dat,
            prep = info$prep,
            candidates = candidates,
            n_select = n_select,
            B = B,
            B_reference = B_reference,
            seed_boot =
              seed_methods[
                s,
                ih,
                ii,
                2
              ],
            n_cross_shifts =
              n_cross_shifts,
            tolerance = tolerance,
            keep_inner_vectors = TRUE
          )

          per_outer[[pos]] <-
            twolayer_summarize_inner(
              zin,
              tolerance = tolerance
            )

          if (
            pos %% checkpoint_every == 0L ||
              pos == length(outer_indices)
          ) {

            saveRDS(
              list(
                S = S,
                B = B,
                key = key,
                outer_indices =
                  outer_indices,
                per_outer =
                  per_outer
              ),
              checkpoint_file
            )

            cat(
              "Completed subset position ",
              pos,
              " / ",
              length(outer_indices),
              " (outer replicate ",
              s,
              ")\n",
              sep = ""
            )
          }
        }
      }

      if (
        any(
          vapply(
            per_outer,
            is.null,
            logical(1)
          )
        )
      ) {
        stop(
          "Scenario contains incomplete results: ",
          key
        )
      }

      # ----------------------------------------------------------
      # Reference row is used only for the exact full 1000 x 500
      # production reproduction check.
      # ----------------------------------------------------------

      reference_row <- NULL

      if (
        full_reference_run &&
          !is.null(cross_reference) &&
          !is.null(cross_reference$key_table)
      ) {

        rr <- cross_reference$key_table[
          abs(
            cross_reference$key_table$h2 -
              h2
          ) < 1e-12 &
            cross_reference$key_table$information ==
              info$name,
          ,
          drop = FALSE
        ]

        if (nrow(rr) == 1L) {
          reference_row <- rr
        }
      }

      summary_table <-
        twolayer_summarize_scenario(
          per_outer = per_outer,
          outer_MSE_VC =
            outer_MSE_VC,
          h2 = h2,
          information =
            info$name,
          n_pheno =
            info$n_pheno,
          reference_row =
            reference_row,
          tolerance = tolerance
        )

      # ----------------------------------------------------------
      # Validation
      # ----------------------------------------------------------

      validation <- data.frame(
        scenario = key,
        check = c(
          "switch_nonnegative",
          "switch_zero_when_same_set",
          "inner_error_identity",
          "inner_second_moment_identity",
          "inner_hierarchy_identity",
          "scenario_second_moment_identity",
          "scenario_hierarchy_identity",
          "scenario_bias_hierarchy_identity",
          "previous_same_plugin_R_reproduction",
          "previous_same_trueVC_R_reproduction",
          "previous_cross_trueVC_R_reproduction",
          "previous_fixed_R_reproduction"
        ),
        max_abs_error = c(
          summary_table$max_negative_switch,
          summary_table$max_abs_switch_when_same,
          summary_table$max_inner_error_identity,
          summary_table$max_inner_second_moment_identity_error,
          summary_table$max_inner_hierarchy_identity_error,
          summary_table$scenario_second_moment_expansion_error,
          summary_table$scenario_hierarchy_identity_error,
          summary_table$scenario_bias_hierarchy_identity_error,
          summary_table$diff_reference_R_VC_same_plugin,
          summary_table$diff_reference_R_VC_same_trueVC,
          summary_table$diff_reference_R_VC_cross_trueVC,
          summary_table$diff_reference_R_VC_fixed
        ),
        row.names = NULL,
        check.names = FALSE
      )

      validation$within_tolerance <- (
        is.na(validation$max_abs_error) |
          validation$max_abs_error <=
            tolerance
      )

      structural_checks <- validation$check %in% c(
        "switch_nonnegative",
        "switch_zero_when_same_set",
        "inner_error_identity",
        "inner_second_moment_identity",
        "inner_hierarchy_identity",
        "scenario_second_moment_identity",
        "scenario_hierarchy_identity",
        "scenario_bias_hierarchy_identity"
      )

      if (
        !all(
          validation$within_tolerance[
            structural_checks
          ]
        )
      ) {
        print(validation)
        stop(
          "Structural validation failed for ",
          key
        )
      }

      # For a full production run, prior reproduction must also hold.
      if (
        full_reference_run &&
          !is.null(reference_row)
      ) {

        prior_checks <- grepl(
          "^previous_",
          validation$check
        )

        if (
          !all(
            validation$within_tolerance[
              prior_checks
            ]
          )
        ) {
          print(validation)
          stop(
            "Previous-result reproduction failed for ",
            key
          )
        }
      }

      summary_object <- list(
        key = key,
        h2 = h2,
        information = info$name,
        n_pheno = info$n_pheno,
        outer_indices = outer_indices,
        summary_table = summary_table,
        validation = validation
      )

      if (keep_per_outer_in_master) {
        summary_object$per_outer <-
          per_outer
      }

      saveRDS(
        list(
          S = S,
          B = B,
          key = key,
          outer_indices =
            outer_indices,
          summary_table =
            summary_table,
          summary_object =
            summary_object,
          validation =
            validation,
          per_outer =
            per_outer
        ),
        final_file
      )

      write.csv(
        summary_table,
        file.path(
          output_dir,
          paste0(
            scenario_tag,
            "_twolayer_summary.csv"
          )
        ),
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )

      per_outer_df <- do.call(
        rbind,
        lapply(
          seq_along(per_outer),
          function(j) {
            z <- per_outer[[j]]
            data.frame(
              outer_replicate =
                outer_indices[j],
              as.data.frame(
                z,
                check.names = FALSE
              ),
              row.names = NULL,
              check.names = FALSE
            )
          }
        )
      )

      write.csv(
        per_outer_df,
        file.path(
          output_dir,
          paste0(
            scenario_tag,
            "_twolayer_per_outer.csv"
          )
        ),
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )

      write.csv(
        validation,
        file.path(
          output_dir,
          paste0(
            scenario_tag,
            "_twolayer_validation.csv"
          )
        ),
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )

      scenario_rows[[row_id]] <-
        summary_table

      scenario_objects[[key]] <-
        summary_object

      validation_rows[[row_id]] <-
        validation

      row_id <- row_id + 1L
    }
  }

  key_table <- do.call(
    rbind,
    scenario_rows
  )

  rownames(key_table) <- NULL

  validation_all <- do.call(
    rbind,
    validation_rows
  )

  rownames(validation_all) <- NULL

  out <- list(
    settings = list(
      S_reference = S,
      B_reference = B_reference,
      outer_indices = outer_indices,
      B_inner = B,
      h2_values = h2_values,
      info_fractions =
        info_fractions,
      n_select = n_select,
      n_cross_shifts =
        n_cross_shifts,
      estimand = paste0(
        "mechanism diagnostic: two-layer decomposition of ",
        "same-sample selection-VC alignment"
      ),
      decomposition = paste0(
        "e_plugin = e_common + e_switch; ",
        "M_plugin = M_cross_true + ",
        "(M_common-M_cross_true) + ",
        "(M_plugin-M_common)"
      ),
      output_dir = normalizePath(
        output_dir,
        winslash = "/",
        mustWork = FALSE
      )
    ),
    key_table = key_table,
    validation = validation_all,
    result = scenario_objects
  )

  write.csv(
    key_table,
    file.path(
      output_dir,
      "two_layer_alignment_key_table.csv"
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    validation_all,
    file.path(
      output_dir,
      "two_layer_alignment_validation.csv"
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  saveRDS(
    out,
    file.path(
      output_dir,
      "two_layer_alignment_master.rds"
    )
  )

  cat(
    "\n====================================================\n",
    "Two-layer alignment decomposition completed.\n",
    "====================================================\n",
    sep = ""
  )

  cat("\nKey table:\n")
  print(key_table)

  cat("\nValidation:\n")
  print(validation_all)

  cat(
    "\nOutput directory:\n",
    normalizePath(
      output_dir,
      winslash = "/",
      mustWork = FALSE
    ),
    "\n",
    sep = ""
  )

  invisible(out)
}


print_two_layer_alignment_decomposition <- function(
    x,
    digits = 4) {

  if (
    is.null(x) ||
      is.null(x$key_table)
  ) {
    stop(
      "Object has no $key_table."
    )
  }

  z <- x$key_table

  cat(
    "\n====================================================\n",
    "TWO-LAYER SELECTION-VC ALIGNMENT DECOMPOSITION\n",
    "====================================================\n",
    sep = ""
  )

  cat(
    "\nHierarchical second-moment decomposition\n",
    "M_plugin = baseline + common-data alignment + policy-switch increment\n",
    sep = ""
  )

  show1 <- data.frame(
    h2 = z$h2,
    n_pheno = z$n_pheno,

    R_plugin = round(
      z$R_VC_plugin,
      digits
    ),

    baseline_fraction = round(
      z$baseline_fraction_of_plugin,
      digits
    ),

    common_alignment_fraction = round(
      z$common_alignment_fraction_of_plugin,
      digits
    ),

    policy_switch_fraction = round(
      z$policy_switch_fraction_of_plugin,
      digits
    ),

    common_plus_switch = round(
      z$common_alignment_fraction_of_plugin +
        z$policy_switch_fraction_of_plugin,
      digits
    ),

    cross_true_over_fixed = round(
      z$M_cross_true / z$M_fixed,
      digits
    ),

    row.names = NULL,
    check.names = FALSE
  )

  print(show1)

  cat(
    "\nExact SAME-sample square expansion\n",
    "E(e_plugin^2) = E(e_common^2) + E(e_switch^2) + 2E(e_common*e_switch)\n",
    sep = ""
  )

  show2 <- data.frame(
    h2 = z$h2,
    n_pheno = z$n_pheno,

    common_raw_fraction = round(
      z$common_raw_fraction_of_plugin,
      digits
    ),

    switch_sq_fraction = round(
      z$switch_sq_fraction_of_plugin,
      digits
    ),

    cross_term_fraction = round(
      z$two_common_switch_fraction_of_plugin,
      digits
    ),

    sum = round(
      z$common_raw_fraction_of_plugin +
        z$switch_sq_fraction_of_plugin +
        z$two_common_switch_fraction_of_plugin,
      digits
    ),

    row.names = NULL,
    check.names = FALSE
  )

  print(show2)

  cat(
    "\nBias decomposition\n",
    "E(e_plugin) = E(e_cross,true) + common-data alignment bias + switch bias\n",
    sep = ""
  )

  show3 <- data.frame(
    h2 = z$h2,
    n_pheno = z$n_pheno,

    bias_plugin = round(
      z$bias_plugin,
      digits
    ),

    bias_cross_true = round(
      z$bias_cross_true,
      digits
    ),

    common_alignment_bias = round(
      z$bias_common_alignment,
      digits
    ),

    switch_bias = round(
      z$bias_switch,
      digits
    ),

    row.names = NULL,
    check.names = FALSE
  )

  print(show3)

  cat(
    "\nSelection / deterministic-switch diagnostics\n",
    sep = ""
  )

  show4 <- data.frame(
    h2 = z$h2,
    n_pheno = z$n_pheno,

    exact_set_match_rate = round(
      z$exact_set_match_rate,
      digits
    ),

    mean_overlap = round(
      z$mean_overlap,
      digits
    ),

    mean_number_replaced = round(
      z$mean_number_replaced,
      digits
    ),

    mean_switch_when_changed = round(
      z$mean_switch_when_changed,
      digits
    ),

    min_switch = signif(
      z$min_switch,
      4
    ),

    nonnegative_switch_rate = round(
      z$nonnegative_switch_rate,
      digits
    ),

    row.names = NULL,
    check.names = FALSE
  )

  print(show4)

  cat(
    "\nInterpretation reminder:\n",
    "- baseline_fraction: independent cross-sample TRUE-VC-selected policy baseline.\n",
    "- common_alignment_fraction: excess caused by using a TRUE-VC-selected policy\n",
    "  from the SAME inner data as d.\n",
    "- policy_switch_fraction: additional second-moment increment from replacing\n",
    "  the TRUE-VC-selected top-k set by the plugin REML-EBLUP top-k set.\n",
    "- These three fractions sum exactly to 1 (up to numerical rounding).\n",
    "- e_switch should be nonnegative for every inner replicate under exact top-k\n",
    "  selection; this is a deterministic optimality result, not a stochastic claim.\n",
    sep = ""
  )

  invisible(
    list(
      hierarchy = show1,
      square_expansion = show2,
      bias = show3,
      selection = show4
    )
  )
}


# ======================================================================
# END STEP 5
# ======================================================================


# ======================================================================
# STEP 7 FIXED v2
# SELECTION-INTENSITY SENSITIVITY ANALYSIS
# k = 2, 5, 10 selected from 20 candidates
# FIX: exact stored outer REML / joint-generator reproduction for k=5
# ======================================================================
#
# PURPOSE
# -------
# Test whether the selection-VC alignment mechanism found at k=5 persists
# when selection intensity changes.
#
# IMPORTANT COMPUTATIONAL DESIGN
# ------------------------------
# For each OUTER replicate and each h2 x information scenario:
#
#   1. generate the OUTER data once;
#   2. fit OUTER REML once;
#   3. generate each TRUE-VC inner data set once;
#   4. fit inner REML once per inner replicate;
#   5. compute plugin and TRUE-VC EBVs once;
#   6. evaluate ALL requested k values from those SAME EBVs.
#
# Therefore k=2,5,10 does NOT triple the expensive REML fitting work.
# The three k values differ only in their top-k selection contrasts after
# the common predictions have been obtained.
#
# PRIMARY DIAGNOSTICS FOR EACH k
# ------------------------------
#
# Let
#
#   d = uhat(hat(theta)) - uhat(theta_true).
#
# SAME-sample plugin:
#
#   e_same,plugin(k) = c_plugin,k' d
#
# SAME-sample TRUE-VC-selected policy:
#
#   e_same,true(k) = c_true,k' d
#
# CROSS-sample plugin / TRUE-VC policy:
#
#   e_cross,plugin(k)
#   e_cross,true(k)
#
# OUTER VC-induced error:
#
#   e_outer,VC(k)
#     = c_plugin,k(Y)' [
#         uhat(hat(theta);Y) - uhat(theta_true;Y)
#       ].
#
# Main ratios:
#
#   R_same,plugin(k)
#     = E_inner[e_same,plugin(k)^2] / E_outer[e_outer,VC(k)^2]
#
#   R_cross,plugin(k)
#     = E_inner[e_cross,plugin(k)^2] / E_outer[e_outer,VC(k)^2]
#
#   P_alignment(k)
#     = 1 - M_cross,plugin(k) / M_same,plugin(k)
#
# Two-layer decomposition:
#
#   M_same,plugin(k)
#     = M_cross,true(k)
#       + {M_same,true(k)-M_cross,true(k)}
#       + {M_same,plugin(k)-M_same,true(k)}
#
# Components divided by M_same,plugin(k):
#
#   baseline_fraction(k)
#   common_alignment_fraction(k)
#   policy_switch_fraction(k)
#
# and
#
#   common_share_of_alignment(k)
#     = common_alignment /
#       (common_alignment + policy_switch).
#
# The policy-switch ERROR itself satisfies
#
#   e_switch(k)
#     = e_same,plugin(k) - e_same,true(k) >= 0
#
# for exact top-k selection.
#
# TARGET OF THIS STEP
# -------------------
# This is a sensitivity / generalisation diagnostic.  It is not a new
# inferential interval method.
#
# PILOT
# -----
#
#
#   ksens20 <- run_selection_intensity_sensitivity(
#     boot1000 = boot1000,
#     layer_reference = layer1000,
#     k_values = c(2,5,10),
#     outer_indices = 1:20,
#     B_inner = 100,
#     checkpoint_every = 5,
#     resume = FALSE,
#     MC_R = 2000
#   )
#
#   print_selection_intensity_sensitivity(ksens20)
#   print(ksens20$validation)
#
# PRODUCTION
# ----------
#
#   ksens1000 <- run_selection_intensity_sensitivity(
#     boot1000 = boot1000,
#     layer_reference = layer1000,
#     k_values = c(2,5,10),
#     checkpoint_every = 10,
#     resume = TRUE,
#     MC_R = 5000
#   )
#
#   print_selection_intensity_sensitivity(ksens1000)
#   print(ksens1000$key_table)
#   print(ksens1000$mc_ci_table)
#   print(ksens1000$validation)
#
# ======================================================================


# ----------------------------------------------------------------------
# Small utilities
# ----------------------------------------------------------------------

ksens_safe_ratio <- function(a, b) {
  ifelse(
    is.finite(a) & is.finite(b) & b != 0,
    a / b,
    NA_real_
  )
}


ksens_q <- function(x, probs = c(0.025, 0.975)) {
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


ksens_choose_shifts <- function(B, n_shifts = 5L) {

  if (B < 2L) {
    stop("B must be >= 2 for cross-sample pairing.")
  }

  max_shifts <- min(
    as.integer(n_shifts),
    B - 1L
  )

  # Deterministic non-zero cyclic shifts, spread across the available range.
  cand <- unique(
    pmax(
      1L,
      pmin(
        B - 1L,
        as.integer(
          round(
            seq(
              from = 1,
              to = B - 1L,
              length.out = max_shifts
            )
          )
        )
      )
    )
  )

  cand[cand != 0L]
}


ksens_shift_index <- function(B, shift) {
  ((seq_len(B) - 1L + as.integer(shift)) %% B) + 1L
}


ksens_top_selected <- function(
    uhat,
    candidates,
    k) {

  candidates[
    order(
      uhat[candidates],
      decreasing = TRUE
    )[seq_len(k)]
  ]
}


ksens_policy_contrast <- function(
    delta_candidate,
    selected_ids,
    candidates) {

  idx <- match(
    as.integer(selected_ids),
    as.integer(candidates)
  )

  if (anyNA(idx)) {
    stop("Selected ID not found in candidate set.")
  }

  mean(delta_candidate[idx]) -
    mean(delta_candidate)
}


ksens_exact_set_match <- function(a, b) {
  identical(
    sort(as.integer(a)),
    sort(as.integer(b))
  )
}


ksens_overlap <- function(a, b) {
  length(
    intersect(
      as.integer(a),
      as.integer(b)
    )
  ) / length(a)
}


# ----------------------------------------------------------------------
# One OUTER replicate, ALL k values simultaneously.
# ----------------------------------------------------------------------

ksens_one_outer_all_k <- function(
    outer_raw,
    dat,
    prep,
    candidates,
    k_values,
    B,
    B_reference,
    seed_boot,
    n_cross_shifts = 5L,
    tolerance = 1e-10) {

  n_animals <- nrow(prep$A)
  n_cand <- length(candidates)

  if (any(k_values < 1L) ||
      any(k_values >= n_cand)) {
    stop(
      "Each k must satisfy 1 <= k < number of candidates."
    )
  }

  if (B < 20L || B > B_reference) {
    stop("Require 20 <= B <= B_reference.")
  }

  # ------------------------------------------------------------
  # OUTER plugin and TRUE-VC EBLUP
  #
  # IMPORTANT FIX:
  # Use the STORED outer REML estimates and stored realised k=5 policy
  # from the completed bootstrap-generator diagnostic rather than
  # refitting REML here.  This is required for exact k=5 reproduction.
  # ------------------------------------------------------------

  if (is.null(outer_raw)) {
    stop("outer_raw is required for exact reference reproduction.")
  }

  c_outer_ref <- as.numeric(outer_raw$cstar)

  if (length(c_outer_ref) != n_animals) {
    stop("outer_raw$cstar has incompatible length.")
  }

  sigma2_A_hat_outer <- as.numeric(
    outer_raw$sigma2_A_hat
  )

  sigma2_e_hat_outer <- as.numeric(
    outer_raw$sigma2_e_hat
  )

  if (
    !is.finite(sigma2_A_hat_outer) ||
      sigma2_A_hat_outer <= 0 ||
      !is.finite(sigma2_e_hat_outer) ||
      sigma2_e_hat_outer <= 0
  ) {
    stop("Stored outer variance-component estimates are invalid.")
  }

  eta_outer_plugin <- log(
    c(
      sigma2_A_hat_outer,
      sigma2_e_hat_outer
    )
  )

  pred_outer_plugin <- predict_ebv_validation(
    eta_outer_plugin,
    dat$y,
    prep
  )$uhat

  eta_true <- log(
    c(
      dat$sigma2_A,
      dat$sigma2_e
    )
  )

  pred_outer_true <- predict_ebv_validation(
    eta_true,
    dat$y,
    prep
  )$uhat

  delta_outer_candidate <- (
    pred_outer_plugin[candidates] -
      pred_outer_true[candidates]
  )

  # EXACT same beta generator as the previous TRUE-VC diagnostics.
  beta_hat_plugin <- bootdiag_gls_beta(
    y = dat$y,
    prep = prep,
    sigma2_A = sigma2_A_hat_outer,
    sigma2_e = sigma2_e_hat_outer
  )

  # EXACT same joint generator as crosssel_one_outer().
  # Although only y* is used below, retaining the (y*, G*) joint
  # covariance and the same Cholesky construction guarantees the same
  # numerical RNG mapping as the completed reference diagnostic.
  joint <- bootdiag_joint_generator(
    prep = prep,
    cstar = c_outer_ref,
    sigma2_A = dat$sigma2_A,
    sigma2_e = dat$sigma2_e,
    beta_generator = beta_hat_plugin
  )

  set.seed(seed_boot)

  z_joint_reference <- matrix(
    rnorm(
      B_reference * (prep$nobs + 1L)
    ),
    nrow = B_reference,
    ncol = prep$nobs + 1L
  )

  z_joint <- z_joint_reference[
    seq_len(B),
    ,
    drop = FALSE
  ]

  draws <- sweep(
    z_joint %*% joint$chol,
    2,
    joint$mean,
    "+"
  )

  y_draws <- draws[
    ,
    seq_len(prep$nobs),
    drop = FALSE
  ]

  # Stored realised k=5 selected set from the completed outer analysis.
  outer_selected_ref <- candidates[
    c_outer_ref[candidates] > 0
  ]

  if (
    5L %in% k_values &&
      length(outer_selected_ref) != 5L
  ) {
    stop(
      "Could not recover the stored outer k=5 selected set from cstar."
    )
  }

  # ------------------------------------------------------------
  # Common inner predictions.  These expensive fits are shared by ALL k.
  # ------------------------------------------------------------

  delta_candidate <- matrix(
    NA_real_,
    nrow = B,
    ncol = n_cand
  )

  pred_plugin_candidate <- matrix(
    NA_real_,
    nrow = B,
    ncol = n_cand
  )

  pred_true_candidate <- matrix(
    NA_real_,
    nrow = B,
    ncol = n_cand
  )

  h2_hat <- rep(NA_real_, B)
  converged <- rep(NA_real_, B)
  boundary <- rep(NA_real_, B)

  for (b in seq_len(B)) {

    yb <- y_draws[b, ]

    fit_b <- fit_reml_animal(
      yb,
      prep
    )

    pred_hat <- predict_ebv_validation(
      fit_b$eta,
      yb,
      prep
    )$uhat

    pred_true <- predict_ebv_validation(
      eta_true,
      yb,
      prep
    )$uhat

    pred_plugin_candidate[b, ] <-
      pred_hat[candidates]

    pred_true_candidate[b, ] <-
      pred_true[candidates]

    delta_candidate[b, ] <- (
      pred_hat[candidates] -
        pred_true[candidates]
    )

    h2_hat[b] <- (
      fit_b$theta[1] /
        sum(fit_b$theta)
    )

    converged[b] <- as.numeric(
      fit_b$convergence == 0
    )

    boundary[b] <- as.numeric(
      fit_b$eta[1] -
        fit_b$lower[1] < 0.05
    )
  }

  ok <- (
    apply(
      is.finite(delta_candidate),
      1L,
      all
    ) &
      apply(
        is.finite(pred_plugin_candidate),
        1L,
        all
      ) &
      apply(
        is.finite(pred_true_candidate),
        1L,
        all
      )
  )

  if (
    sum(ok) <
      max(
        20L,
        ceiling(0.80 * B)
      )
  ) {
    stop(
      "Too few valid inner replicates: ",
      sum(ok),
      " / ",
      B
    )
  }

  valid_idx <- which(ok)

  dc <- delta_candidate[
    valid_idx,
    ,
    drop = FALSE
  ]

  pp <- pred_plugin_candidate[
    valid_idx,
    ,
    drop = FALSE
  ]

  pt <- pred_true_candidate[
    valid_idx,
    ,
    drop = FALSE
  ]

  Bv <- nrow(dc)

  shifts <- ksens_choose_shifts(
    Bv,
    n_shifts = n_cross_shifts
  )

  # ------------------------------------------------------------
  # Evaluate every k using the SAME prediction matrices.
  # ------------------------------------------------------------

  by_k <- vector(
    "list",
    length(k_values)
  )

  names(by_k) <- paste0(
    "k",
    k_values
  )

  for (ik in seq_along(k_values)) {

    k <- as.integer(k_values[ik])

    # Outer realised plugin policy for this k.
    outer_selected_plugin_calc <- ksens_top_selected(
      pred_outer_plugin,
      candidates,
      k
    )

    # For k=5, use the stored realised outer policy exactly.
    # The calculated ranking is checked against it.
    if (k == 5L) {

      outer_selected_plugin <- outer_selected_ref

      outer_k5_set_match <- as.numeric(
        ksens_exact_set_match(
          outer_selected_plugin_calc,
          outer_selected_ref
        )
      )

      if (outer_k5_set_match != 1) {
        stop(
          "Stored and recalculated outer k=5 selected sets differ."
        )
      }

    } else {

      outer_selected_plugin <- outer_selected_plugin_calc
      outer_k5_set_match <- NA_real_
    }

    outer_eVC <- ksens_policy_contrast(
      delta_outer_candidate,
      outer_selected_plugin,
      candidates
    )

    # Direct k=5 outer-error reproduction check against the stored
    # conditional plugin and TRUE-VC oracle means.
    outer_k5_eVC_reference_difference <- NA_real_

    if (
      k == 5L &&
        !is.null(outer_raw$conditional_plugin) &&
        !is.null(outer_raw$conditional_trueVC_oracle)
    ) {

      ref_outer_eVC <- (
        as.numeric(
          outer_raw$conditional_plugin["mean"]
        ) -
          as.numeric(
            outer_raw$conditional_trueVC_oracle["mean"]
          )
      )

      outer_k5_eVC_reference_difference <- abs(
        outer_eVC - ref_outer_eVC
      )

      if (
        is.finite(outer_k5_eVC_reference_difference) &&
          outer_k5_eVC_reference_difference > tolerance
      ) {
        stop(
          "Stored k=5 outer VC error was not reproduced: ",
          format(
            outer_k5_eVC_reference_difference,
            digits = 16
          )
        )
      }
    }

    selected_plugin <- matrix(
      NA_integer_,
      nrow = Bv,
      ncol = k
    )

    selected_true <- matrix(
      NA_integer_,
      nrow = Bv,
      ncol = k
    )

    e_fixed <- rep(NA_real_, Bv)
    e_same_plugin <- rep(NA_real_, Bv)
    e_same_true <- rep(NA_real_, Bv)

    exact_match <- rep(NA_real_, Bv)
    overlap <- rep(NA_real_, Bv)

    for (b in seq_len(Bv)) {

      # Candidate-only orderings.
      ord_plugin <- order(
        pp[b, ],
        decreasing = TRUE
      )

      ord_true <- order(
        pt[b, ],
        decreasing = TRUE
      )

      sp <- candidates[
        ord_plugin[seq_len(k)]
      ]

      st <- candidates[
        ord_true[seq_len(k)]
      ]

      selected_plugin[b, ] <- sp
      selected_true[b, ] <- st

      e_fixed[b] <- ksens_policy_contrast(
        dc[b, ],
        outer_selected_plugin,
        candidates
      )

      e_same_plugin[b] <- ksens_policy_contrast(
        dc[b, ],
        sp,
        candidates
      )

      e_same_true[b] <- ksens_policy_contrast(
        dc[b, ],
        st,
        candidates
      )

      exact_match[b] <- as.numeric(
        ksens_exact_set_match(
          sp,
          st
        )
      )

      overlap[b] <- ksens_overlap(
        sp,
        st
      )
    }

    e_switch <- (
      e_same_plugin -
        e_same_true
    )

    # Deterministic top-k inequality.
    max_negative_switch <- max(
      pmax(
        -e_switch,
        0
      ),
      na.rm = TRUE
    )

    nonnegative_switch_rate <- mean(
      e_switch >= -tolerance
    )

    if (
      max_negative_switch >
        tolerance
    ) {
      stop(
        "Top-k switch nonnegativity failed for k=",
        k,
        ": ",
        format(
          max_negative_switch,
          digits = 16
        )
      )
    }

    # Independent cross-sample policies.
    e_cross_plugin <- numeric(
      Bv * length(shifts)
    )

    e_cross_true <- numeric(
      Bv * length(shifts)
    )

    pos <- 1L

    for (sh in shifts) {

      jj <- ksens_shift_index(
        Bv,
        sh
      )

      for (b in seq_len(Bv)) {

        rng <- pos + b - 1L

        e_cross_plugin[rng] <-
          ksens_policy_contrast(
            dc[b, ],
            selected_plugin[jj[b], ],
            candidates
          )

        e_cross_true[rng] <-
          ksens_policy_contrast(
            dc[b, ],
            selected_true[jj[b], ],
            candidates
          )
      }

      pos <- pos + Bv
    }

    M_same_plugin <- mean(
      e_same_plugin^2
    )

    M_same_true <- mean(
      e_same_true^2
    )

    M_cross_plugin <- mean(
      e_cross_plugin^2
    )

    M_cross_true <- mean(
      e_cross_true^2
    )

    M_fixed <- mean(
      e_fixed^2
    )

    M_switch_sq <- mean(
      e_switch^2
    )

    two_common_switch <- (
      2 *
        mean(
          e_same_true *
            e_switch
        )
    )

    second_moment_identity_error <- abs(
      M_same_plugin -
        (
          M_same_true +
            M_switch_sq +
            two_common_switch
        )
    )

    baseline_fraction <- ksens_safe_ratio(
      M_cross_true,
      M_same_plugin
    )

    common_alignment_fraction <- ksens_safe_ratio(
      M_same_true -
        M_cross_true,
      M_same_plugin
    )

    policy_switch_fraction <- ksens_safe_ratio(
      M_same_plugin -
        M_same_true,
      M_same_plugin
    )

    alignment_total_fraction <- (
      common_alignment_fraction +
        policy_switch_fraction
    )

    hierarchy_identity_error <- abs(
      (
        baseline_fraction +
          common_alignment_fraction +
          policy_switch_fraction
      ) - 1
    )

    changed <- exact_match < 0.5

    mean_switch_when_changed <- if (any(changed)) {
      mean(
        e_switch[changed]
      )
    } else {
      0
    }

    by_k[[ik]] <- list(
      k = k,
      selected_fraction = k / n_cand,
      n_valid = Bv,

      # OUTER error for denominator
      outer_eVC = outer_eVC,
      outer_k5_set_match = outer_k5_set_match,
      outer_k5_eVC_reference_difference =
        outer_k5_eVC_reference_difference,

      # moments
      M_same_plugin = M_same_plugin,
      M_cross_plugin = M_cross_plugin,
      M_same_true = M_same_true,
      M_cross_true = M_cross_true,
      M_fixed = M_fixed,
      M_switch_sq = M_switch_sq,
      two_common_switch = two_common_switch,

      # first moments
      bias_same_plugin =
        mean(e_same_plugin),
      bias_cross_plugin =
        mean(e_cross_plugin),
      bias_same_true =
        mean(e_same_true),
      bias_cross_true =
        mean(e_cross_true),
      bias_switch =
        mean(e_switch),

      # alignment / decomposition
      P_alignment_plugin =
        1 -
          ksens_safe_ratio(
            M_cross_plugin,
            M_same_plugin
          ),

      baseline_fraction =
        baseline_fraction,

      common_alignment_fraction =
        common_alignment_fraction,

      policy_switch_fraction =
        policy_switch_fraction,

      alignment_total_fraction =
        alignment_total_fraction,

      common_share_of_alignment =
        ksens_safe_ratio(
          common_alignment_fraction,
          alignment_total_fraction
        ),

      switch_share_of_alignment =
        ksens_safe_ratio(
          policy_switch_fraction,
          alignment_total_fraction
        ),

      cross_plugin_over_fixed =
        ksens_safe_ratio(
          M_cross_plugin,
          M_fixed
        ),

      cross_true_over_fixed =
        ksens_safe_ratio(
          M_cross_true,
          M_fixed
        ),

      # selection
      exact_set_match_rate =
        mean(exact_match),
      mean_overlap =
        mean(overlap),
      mean_number_replaced =
        k * (1 - mean(overlap)),
      mean_switch_when_changed =
        mean_switch_when_changed,

      # structural validation
      max_negative_switch =
        max_negative_switch,
      nonnegative_switch_rate =
        nonnegative_switch_rate,
      second_moment_identity_error =
        second_moment_identity_error,
      hierarchy_identity_error =
        hierarchy_identity_error
    )
  }

  list(
    by_k = by_k,
    mean_h2_hat =
      mean(
        h2_hat[valid_idx],
        na.rm = TRUE
      ),
    boundary_rate =
      mean(
        boundary[valid_idx],
        na.rm = TRUE
      ),
    convergence_rate =
      mean(
        converged[valid_idx],
        na.rm = TRUE
      )
  )
}


# ----------------------------------------------------------------------
# Summarise one h2 x information scenario across OUTER replicates.
# ----------------------------------------------------------------------

ksens_summarize_scenario <- function(
    per_outer,
    k_values,
    h2,
    information,
    n_pheno) {

  rows <- vector(
    "list",
    length(k_values)
  )

  for (ik in seq_along(k_values)) {

    k <- as.integer(k_values[ik])
    kn <- paste0("k", k)

    grab <- function(name) {
      vapply(
        per_outer,
        function(z) {
          as.numeric(
            z$by_k[[kn]][[name]]
          )
        },
        numeric(1)
      )
    }

    outer_eVC <- grab(
      "outer_eVC"
    )

    outer_MSE_VC <- mean(
      outer_eVC^2
    )

    M_same_plugin <- mean(
      grab("M_same_plugin")
    )

    M_cross_plugin <- mean(
      grab("M_cross_plugin")
    )

    M_same_true <- mean(
      grab("M_same_true")
    )

    M_cross_true <- mean(
      grab("M_cross_true")
    )

    M_fixed <- mean(
      grab("M_fixed")
    )

    M_switch_sq <- mean(
      grab("M_switch_sq")
    )

    two_common_switch <- mean(
      grab("two_common_switch")
    )

    baseline_fraction <- ksens_safe_ratio(
      M_cross_true,
      M_same_plugin
    )

    common_alignment_fraction <- ksens_safe_ratio(
      M_same_true -
        M_cross_true,
      M_same_plugin
    )

    policy_switch_fraction <- ksens_safe_ratio(
      M_same_plugin -
        M_same_true,
      M_same_plugin
    )

    alignment_total_fraction <- (
      common_alignment_fraction +
        policy_switch_fraction
    )

    rows[[ik]] <- data.frame(
      h2 = h2,
      information = information,
      n_pheno = n_pheno,
      k = k,
      selected_fraction =
        k / 20,
      n_outer =
        length(per_outer),

      outer_MSE_VC =
        outer_MSE_VC,

      M_same_plugin =
        M_same_plugin,
      M_cross_plugin =
        M_cross_plugin,
      M_same_true =
        M_same_true,
      M_cross_true =
        M_cross_true,
      M_fixed =
        M_fixed,

      R_same_plugin =
        ksens_safe_ratio(
          M_same_plugin,
          outer_MSE_VC
        ),

      R_cross_plugin =
        ksens_safe_ratio(
          M_cross_plugin,
          outer_MSE_VC
        ),

      R_same_true =
        ksens_safe_ratio(
          M_same_true,
          outer_MSE_VC
        ),

      R_cross_true =
        ksens_safe_ratio(
          M_cross_true,
          outer_MSE_VC
        ),

      R_fixed =
        ksens_safe_ratio(
          M_fixed,
          outer_MSE_VC
        ),

      P_alignment_plugin =
        1 -
          ksens_safe_ratio(
            M_cross_plugin,
            M_same_plugin
          ),

      baseline_fraction =
        baseline_fraction,

      common_alignment_fraction =
        common_alignment_fraction,

      policy_switch_fraction =
        policy_switch_fraction,

      alignment_total_fraction =
        alignment_total_fraction,

      common_share_of_alignment =
        ksens_safe_ratio(
          common_alignment_fraction,
          alignment_total_fraction
        ),

      switch_share_of_alignment =
        ksens_safe_ratio(
          policy_switch_fraction,
          alignment_total_fraction
        ),

      common_raw_fraction =
        ksens_safe_ratio(
          M_same_true,
          M_same_plugin
        ),

      switch_sq_fraction =
        ksens_safe_ratio(
          M_switch_sq,
          M_same_plugin
        ),

      cross_term_fraction =
        ksens_safe_ratio(
          two_common_switch,
          M_same_plugin
        ),

      cross_plugin_over_fixed =
        ksens_safe_ratio(
          M_cross_plugin,
          M_fixed
        ),

      cross_true_over_fixed =
        ksens_safe_ratio(
          M_cross_true,
          M_fixed
        ),

      bias_same_plugin =
        mean(
          grab(
            "bias_same_plugin"
          )
        ),

      bias_cross_plugin =
        mean(
          grab(
            "bias_cross_plugin"
          )
        ),

      bias_same_true =
        mean(
          grab(
            "bias_same_true"
          )
        ),

      bias_cross_true =
        mean(
          grab(
            "bias_cross_true"
          )
        ),

      bias_switch =
        mean(
          grab(
            "bias_switch"
          )
        ),

      exact_set_match_rate =
        mean(
          grab(
            "exact_set_match_rate"
          )
        ),

      mean_overlap =
        mean(
          grab(
            "mean_overlap"
          )
        ),

      mean_number_replaced =
        mean(
          grab(
            "mean_number_replaced"
          )
        ),

      mean_switch_when_changed =
        mean(
          grab(
            "mean_switch_when_changed"
          )
        ),

      max_negative_switch =
        max(
          grab(
            "max_negative_switch"
          )
        ),

      nonnegative_switch_rate =
        mean(
          grab(
            "nonnegative_switch_rate"
          )
        ),

      max_second_moment_identity_error =
        max(
          grab(
            "second_moment_identity_error"
          )
        ),

      max_hierarchy_identity_error =
        max(
          grab(
            "hierarchy_identity_error"
          )
        ),

      min_outer_k5_set_match =
        if (k == 5L) {
          min(
            grab("outer_k5_set_match"),
            na.rm = TRUE
          )
        } else {
          NA_real_
        },

      max_outer_k5_eVC_reference_difference =
        if (k == 5L) {
          max(
            grab("outer_k5_eVC_reference_difference"),
            na.rm = TRUE
          )
        } else {
          NA_real_
        },

      row.names = NULL,
      check.names = FALSE
    )
  }

  do.call(
    rbind,
    rows
  )
}


# ----------------------------------------------------------------------
# Monte Carlo intervals from OUTER resampling.
# ----------------------------------------------------------------------

ksens_MC_CI <- function(
    per_outer,
    k_values,
    h2,
    information,
    n_pheno,
    R = 5000L,
    seed = 20260824L) {

  S <- length(per_outer)

  set.seed(seed)

  idx <- matrix(
    sample.int(
      S,
      size = S * R,
      replace = TRUE
    ),
    nrow = S,
    ncol = R
  )

  rows <- list()
  rr <- 1L

  quantities <- c(
    "P_alignment_plugin",
    "baseline_fraction",
    "common_alignment_fraction",
    "policy_switch_fraction",
    "alignment_total_fraction",
    "common_share_of_alignment",
    "R_same_plugin",
    "R_cross_plugin",
    "cross_plugin_over_fixed",
    "exact_set_match_rate",
    "mean_overlap"
  )

  for (k in k_values) {

    kn <- paste0("k", k)

    grab <- function(name) {
      vapply(
        per_outer,
        function(z) {
          as.numeric(
            z$by_k[[kn]][[name]]
          )
        },
        numeric(1)
      )
    }

    # Base per-outer quantities/moments.
    outer_eVC <- grab("outer_eVC")

    Mp <- grab("M_same_plugin")
    Mcp <- grab("M_cross_plugin")
    Mt <- grab("M_same_true")
    Mct <- grab("M_cross_true")
    Mf <- grab("M_fixed")

    exm <- grab("exact_set_match_rate")
    ov <- grab("mean_overlap")

    mean_idx <- function(x) {
      colMeans(
        matrix(
          x[idx],
          nrow = S,
          ncol = R
        )
      )
    }

    outer_MSE <- mean_idx(
      outer_eVC^2
    )

    b_Mp <- mean_idx(Mp)
    b_Mcp <- mean_idx(Mcp)
    b_Mt <- mean_idx(Mt)
    b_Mct <- mean_idx(Mct)
    b_Mf <- mean_idx(Mf)

    b_R_same <- ksens_safe_ratio(
      b_Mp,
      outer_MSE
    )

    b_R_cross <- ksens_safe_ratio(
      b_Mcp,
      outer_MSE
    )

    b_P <- 1 -
      ksens_safe_ratio(
        b_Mcp,
        b_Mp
      )

    b_baseline <- ksens_safe_ratio(
      b_Mct,
      b_Mp
    )

    b_common <- ksens_safe_ratio(
      b_Mt - b_Mct,
      b_Mp
    )

    b_switch <- ksens_safe_ratio(
      b_Mp - b_Mt,
      b_Mp
    )

    b_total <- b_common + b_switch

    b_common_share <- ksens_safe_ratio(
      b_common,
      b_total
    )

    b_cross_fixed <- ksens_safe_ratio(
      b_Mcp,
      b_Mf
    )

    b_exact <- mean_idx(exm)
    b_overlap <- mean_idx(ov)

    boot <- list(
      P_alignment_plugin =
        b_P,
      baseline_fraction =
        b_baseline,
      common_alignment_fraction =
        b_common,
      policy_switch_fraction =
        b_switch,
      alignment_total_fraction =
        b_total,
      common_share_of_alignment =
        b_common_share,
      R_same_plugin =
        b_R_same,
      R_cross_plugin =
        b_R_cross,
      cross_plugin_over_fixed =
        b_cross_fixed,
      exact_set_match_rate =
        b_exact,
      mean_overlap =
        b_overlap
    )

    # Point estimate from the original OUTER sample.
    pt <- ksens_summarize_scenario(
      per_outer = per_outer,
      k_values = k,
      h2 = h2,
      information = information,
      n_pheno = n_pheno
    )

    for (nm in quantities) {

      qq <- ksens_q(
        boot[[nm]]
      )

      rows[[rr]] <- data.frame(
        h2 = h2,
        information =
          information,
        n_pheno =
          n_pheno,
        k = k,
        quantity = nm,
        estimate =
          as.numeric(pt[[nm]][1]),
        MC_boot_SE =
          sd(
            boot[[nm]],
            na.rm = TRUE
          ),
        CI95_lo = qq[1],
        CI95_hi = qq[2],
        row.names = NULL,
        check.names = FALSE
      )

      rr <- rr + 1L
    }
  }

  do.call(
    rbind,
    rows
  )
}


# ----------------------------------------------------------------------
# Main runner
# ----------------------------------------------------------------------

run_selection_intensity_sensitivity <- function(
    boot1000,
    layer_reference = NULL,
    k_values = c(2L, 5L, 10L),
    outer_indices = NULL,
    B_inner = NULL,
    prior_diagnostic_dir = NULL,
    output_dir = NULL,
    checkpoint_every = 10,
    resume = TRUE,
    tolerance = 1e-10,
    n_cross_shifts = 5L,
    MC_R = 5000L,
    MC_seed = 20260824L,
    keep_per_outer_in_master = FALSE) {

  if (
    is.null(boot1000) ||
      is.null(boot1000$settings)
  ) {
    stop(
      "boot1000 with $settings is required."
    )
  }

  st <- boot1000$settings

  S <- as.integer(st$S)
  B_reference <- as.integer(st$B)

  if (is.null(B_inner)) {
    B <- B_reference
  } else {
    B <- as.integer(B_inner)
  }

  k_values <- sort(
    unique(
      as.integer(k_values)
    )
  )

  if (
    length(k_values) < 2L
  ) {
    warning(
      "Only one k supplied; sensitivity comparison will be limited."
    )
  }

  if (
    any(k_values < 1L) ||
      any(k_values >= 20L)
  ) {
    stop(
      "For the default 20 candidates, require 1 <= k < 20."
    )
  }

  if (
    B < 20L ||
      B > B_reference
  ) {
    stop(
      "B_inner must satisfy 20 <= B_inner <= ",
      B_reference
    )
  }

  if (is.null(outer_indices)) {
    outer_indices <- seq_len(S)
  }

  outer_indices <- sort(
    unique(
      as.integer(outer_indices)
    )
  )

  if (
    length(outer_indices) == 0L ||
      any(outer_indices < 1L) ||
      any(outer_indices > S)
  ) {
    stop(
      "outer_indices must lie within 1:S."
    )
  }

  h2_values <- st$h2_values
  info_fractions <- st$info_fractions
  sigma2_P <- st$sigma2_P
  beta <- st$beta
  seed_info <- st$seed_info
  seed_outer <- st$seed_outer

  if (is.null(prior_diagnostic_dir)) {
    prior_diagnostic_dir <- st$output_dir
  }

  if (
    is.null(prior_diagnostic_dir) ||
      !dir.exists(prior_diagnostic_dir)
  ) {
    stop(
      "Previous bootstrap-generator diagnostic directory not found: ",
      prior_diagnostic_dir
    )
  }

  if (is.null(output_dir)) {
    output_dir <- file.path(
      getwd(),
      "simulation_I_II_final_output",
      paste0(
        "Simulation_II_selection_intensity_FIXED_v2_S",
        length(outer_indices),
        "_B",
        B,
        "_K",
        paste(
          k_values,
          collapse = "_"
        )
      )
    )
  }

  ovib_dir_create(
    output_dir
  )

  base_design <- build_simII_design()

  info_bundle <- build_information_level_designs(
    base_design = base_design,
    info_fractions = info_fractions,
    seed_info = seed_info
  )

  candidates <- base_design$candidates

  if (length(candidates) != 20L) {
    stop(
      "This STEP 7 implementation expects the established 20-candidate design."
    )
  }

  cat(
    "\n====================================================\n",
    "SELECTION-INTENSITY SENSITIVITY\n",
    "====================================================\n",
    "k values: ",
    paste(
      k_values,
      collapse = ", "
    ),
    "\n",
    "Outer replicates: ",
    length(outer_indices),
    " / ",
    S,
    "\n",
    "Inner B: ",
    B,
    "\n",
    "All k share the SAME inner REML fits.\n",
    sep = ""
  )

  print_information_design(
    info_bundle
  )

  # ------------------------------------------------------------
  # Reconstruct established outer RNG schedule.
  # ------------------------------------------------------------

  set.seed(seed_outer)

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
  n_i <- length(
    info_bundle$info_levels
  )

  seed_methods <- array(
    sample.int(
      .Machine$integer.max,
      S * n_h * n_i * 3
    ),
    dim = c(
      S,
      n_h,
      n_i,
      3
    )
  )

  key_rows <- list()
  ci_rows <- list()
  validation_rows <- list()
  scenario_objects <- list()

  row_id <- 1L
  ci_id <- 1L

  full_reference_run <- (
    identical(
      outer_indices,
      seq_len(S)
    ) &&
      identical(
        B,
        B_reference
      )
  )

  for (ih in seq_along(h2_values)) {

    h2 <- h2_values[ih]

    for (
      ii in seq_along(
        info_bundle$info_levels
      )
    ) {

      info <- info_bundle$info_levels[[ii]]

      key <- paste0(
        "h2_",
        h2,
        "__",
        info$name
      )

      cat(
        "\n====================================================\n",
        "k sensitivity: ",
        key,
        " | n_pheno=",
        info$n_pheno,
        "\n====================================================\n",
        sep = ""
      )

      # --------------------------------------------------------
      # Load the exact completed bootstrap-generator scenario.
      # Its stored outer REML estimates and cstar define the reference
      # used by the preceding same-vs-cross / two-layer diagnostics.
      # --------------------------------------------------------

      prior_file <- file.path(
        prior_diagnostic_dir,
        paste0(
          final_safe_name(key),
          "_bootstrap_generator_final.rds"
        )
      )

      if (!file.exists(prior_file)) {
        stop(
          "Previous scenario final RDS not found: ",
          prior_file
        )
      }

      prior <- readRDS(prior_file)

      prior_compatible <- (
        !is.null(prior$S) &&
          !is.null(prior$B) &&
          !is.null(prior$key) &&
          !is.null(prior$raw) &&
          identical(as.integer(prior$S), S) &&
          identical(as.integer(prior$B), B_reference) &&
          identical(prior$key, key) &&
          length(prior$raw) == S
      )

      if (!isTRUE(prior_compatible)) {
        stop(
          "Previous scenario RDS is incompatible: ",
          key
        )
      }

      outer_raw_full <- prior$raw

      scenario_tag <- paste0(
        ovib_safe_name(key),
        "__K",
        paste(
          k_values,
          collapse = "_"
        ),
        "__nOuter",
        length(outer_indices)
      )

      checkpoint_file <- file.path(
        output_dir,
        paste0(
          scenario_tag,
          "_ksens_checkpoint.rds"
        )
      )

      final_file <- file.path(
        output_dir,
        paste0(
          scenario_tag,
          "_ksens_final.rds"
        )
      )

      if (
        resume &&
          file.exists(final_file)
      ) {

        done <- readRDS(
          final_file
        )

        compatible <- (
          !is.null(done$S) &&
            !is.null(done$B) &&
            !is.null(done$key) &&
            !is.null(done$k_values) &&
            !is.null(done$outer_indices) &&
            identical(
              as.integer(done$S),
              S
            ) &&
            identical(
              as.integer(done$B),
              B
            ) &&
            identical(
              done$key,
              key
            ) &&
            identical(
              as.integer(done$k_values),
              k_values
            ) &&
            identical(
              as.integer(done$outer_indices),
              outer_indices
            )
        )

        if (isTRUE(compatible)) {

          cat(
            "Loading completed compatible scenario.\n"
          )

          key_rows[[row_id]] <-
            done$summary_table

          if (
            !is.null(done$mc_ci_table)
          ) {
            ci_rows[[ci_id]] <-
              done$mc_ci_table
            ci_id <- ci_id + 1L
          }

          validation_rows[[row_id]] <-
            done$validation

          scenario_objects[[key]] <-
            done$summary_object

          row_id <- row_id + 1L
          next
        }
      }

      per_outer <- vector(
        "list",
        length(outer_indices)
      )

      start_pos <- 1L

      if (
        resume &&
          file.exists(
            checkpoint_file
          )
      ) {

        cp <- readRDS(
          checkpoint_file
        )

        compatible_cp <- (
          !is.null(cp$S) &&
            !is.null(cp$B) &&
            !is.null(cp$key) &&
            !is.null(cp$k_values) &&
            !is.null(cp$outer_indices) &&
            identical(
              as.integer(cp$S),
              S
            ) &&
            identical(
              as.integer(cp$B),
              B
            ) &&
            identical(
              cp$key,
              key
            ) &&
            identical(
              as.integer(cp$k_values),
              k_values
            ) &&
            identical(
              as.integer(cp$outer_indices),
              outer_indices
            )
        )

        if (isTRUE(compatible_cp)) {

          per_outer <- cp$per_outer

          missing_pos <- which(
            vapply(
              per_outer,
              is.null,
              logical(1)
            )
          )

          if (
            length(missing_pos) == 0L
          ) {
            start_pos <-
              length(outer_indices) + 1L
          } else {
            start_pos <-
              min(missing_pos)
          }

          cat(
            "Resuming at subset position ",
            start_pos,
            ".\n",
            sep = ""
          )
        }
      }

      if (
        start_pos <=
          length(outer_indices)
      ) {

        for (
          pos in seq.int(
            start_pos,
            length(outer_indices)
          )
        ) {

          s <- outer_indices[pos]

          dat <- make_dataset_from_latent(
            latent =
              latent_list[[s]],
            info_level =
              info,
            sigma2_P =
              sigma2_P,
            h2 =
              h2,
            beta =
              beta
          )

          per_outer[[pos]] <-
            ksens_one_outer_all_k(
              outer_raw = outer_raw_full[[s]],
              dat = dat,
              prep = info$prep,
              candidates =
                candidates,
              k_values =
                k_values,
              B = B,
              B_reference =
                B_reference,
              seed_boot =
                seed_methods[
                  s,
                  ih,
                  ii,
                  2
                ],
              n_cross_shifts =
                n_cross_shifts,
              tolerance =
                tolerance
            )

          if (
            pos %%
              checkpoint_every ==
              0L ||
              pos ==
                length(
                  outer_indices
                )
          ) {

            saveRDS(
              list(
                S = S,
                B = B,
                key = key,
                k_values =
                  k_values,
                outer_indices =
                  outer_indices,
                per_outer =
                  per_outer
              ),
              checkpoint_file
            )

            cat(
              "Completed subset position ",
              pos,
              " / ",
              length(outer_indices),
              " (outer replicate ",
              s,
              ")\n",
              sep = ""
            )
          }
        }
      }

      if (
        any(
          vapply(
            per_outer,
            is.null,
            logical(1)
          )
        )
      ) {
        stop(
          "Incomplete scenario: ",
          key
        )
      }

      summary_table <-
        ksens_summarize_scenario(
          per_outer =
            per_outer,
          k_values =
            k_values,
          h2 =
            h2,
          information =
            info$name,
          n_pheno =
            info$n_pheno
        )

      mc_ci_table <- NULL

      if (
        !is.null(MC_R) &&
          is.finite(MC_R) &&
          MC_R > 0
      ) {

        mc_ci_table <-
          ksens_MC_CI(
            per_outer =
              per_outer,
            k_values =
              k_values,
            h2 =
              h2,
            information =
              info$name,
            n_pheno =
              info$n_pheno,
            R =
              as.integer(MC_R),
            seed =
              as.integer(
                MC_seed +
                  100L * ih +
                  ii
              )
          )
      }

      # ----------------------------------------------------------
      # Structural validation + k=5 exact reproduction against
      # the completed layer1000 reference when full production
      # settings are used.
      # ----------------------------------------------------------

      validation <- data.frame(
        scenario = character(),
        k = integer(),
        check = character(),
        max_abs_error = numeric(),
        within_tolerance = logical(),
        stringsAsFactors = FALSE
      )

      add_check <- function(
          k,
          check,
          err) {

        data.frame(
          scenario = key,
          k = k,
          check = check,
          max_abs_error =
            as.numeric(err),
          within_tolerance =
            is.na(err) ||
              as.numeric(err) <=
                tolerance,
          stringsAsFactors =
            FALSE
        )
      }

      vv <- list()
      vi <- 1L

      for (
        r in seq_len(
          nrow(summary_table)
        )
      ) {

        k <- summary_table$k[r]

        vv[[vi]] <- add_check(
          k,
          "switch_nonnegative",
          summary_table$max_negative_switch[r]
        )
        vi <- vi + 1L

        vv[[vi]] <- add_check(
          k,
          "second_moment_identity",
          summary_table$max_second_moment_identity_error[r]
        )
        vi <- vi + 1L

        vv[[vi]] <- add_check(
          k,
          "hierarchy_identity",
          summary_table$max_hierarchy_identity_error[r]
        )
        vi <- vi + 1L

        if (k == 5L) {

          vv[[vi]] <- add_check(
            k,
            "outer_k5_selected_set_reproduction",
            abs(
              summary_table$min_outer_k5_set_match[r] - 1
            )
          )
          vi <- vi + 1L

          vv[[vi]] <- add_check(
            k,
            "outer_k5_VC_error_reproduction",
            summary_table$max_outer_k5_eVC_reference_difference[r]
          )
          vi <- vi + 1L
        }
      }

      # k=5 should reproduce the completed full reference exactly
      # (within floating-point tolerance) under the full 1000 x 500 run.
      if (
        full_reference_run &&
          5L %in% k_values &&
          !is.null(layer_reference) &&
          !is.null(
            layer_reference$key_table
          )
      ) {

        ref <- layer_reference$key_table[
          abs(
            layer_reference$key_table$h2 -
              h2
          ) < 1e-12 &
            layer_reference$key_table$information ==
              info$name,
          ,
          drop = FALSE
        ]

        cur <- summary_table[
          summary_table$k == 5L,
          ,
          drop = FALSE
        ]

        if (
          nrow(ref) == 1L &&
            nrow(cur) == 1L
        ) {

          checks <- list(
            R_same_plugin =
              abs(
                cur$R_same_plugin -
                  ref$R_VC_plugin
              ),

            R_same_true =
              abs(
                cur$R_same_true -
                  ref$R_VC_common_trueVC_selected
              ),

            R_cross_true =
              abs(
                cur$R_cross_true -
                  ref$R_VC_cross_trueVC
              ),

            baseline_fraction =
              abs(
                cur$baseline_fraction -
                  ref$baseline_fraction_of_plugin
              ),

            common_alignment_fraction =
              abs(
                cur$common_alignment_fraction -
                  ref$common_alignment_fraction_of_plugin
              ),

            policy_switch_fraction =
              abs(
                cur$policy_switch_fraction -
                  ref$policy_switch_fraction_of_plugin
              ),

            exact_set_match_rate =
              abs(
                cur$exact_set_match_rate -
                  ref$exact_set_match_rate
              )
          )

          for (
            nm in names(checks)
          ) {
            vv[[vi]] <- add_check(
              5L,
              paste0(
                "k5_reference_",
                nm
              ),
              checks[[nm]]
            )
            vi <- vi + 1L
          }
        }
      }

      validation <- do.call(
        rbind,
        vv
      )

      structural <- !grepl(
        "^k5_reference_",
        validation$check
      )

      if (
        !all(
          validation$within_tolerance[
            structural
          ]
        )
      ) {
        print(validation)
        stop(
          "Structural validation failed for ",
          key
        )
      }

      if (
        full_reference_run &&
          5L %in% k_values &&
          any(
            grepl(
              "^k5_reference_",
              validation$check
            )
          )
      ) {

        ref_idx <- grepl(
          "^k5_reference_",
          validation$check
        )

        if (
          !all(
            validation$within_tolerance[
              ref_idx
            ]
          )
        ) {
          print(validation)
          stop(
            "k=5 reference reproduction failed for ",
            key
          )
        }
      }

      summary_object <- list(
        key = key,
        h2 = h2,
        information =
          info$name,
        n_pheno =
          info$n_pheno,
        k_values =
          k_values,
        outer_indices =
          outer_indices,
        summary_table =
          summary_table,
        mc_ci_table =
          mc_ci_table,
        validation =
          validation
      )

      if (
        keep_per_outer_in_master
      ) {
        summary_object$per_outer <-
          per_outer
      }

      saveRDS(
        list(
          S = S,
          B = B,
          key = key,
          k_values =
            k_values,
          outer_indices =
            outer_indices,
          summary_table =
            summary_table,
          mc_ci_table =
            mc_ci_table,
          validation =
            validation,
          per_outer =
            per_outer,
          summary_object =
            summary_object
        ),
        final_file
      )

      write.csv(
        summary_table,
        file.path(
          output_dir,
          paste0(
            scenario_tag,
            "_ksens_summary.csv"
          )
        ),
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )

      if (
        !is.null(
          mc_ci_table
        )
      ) {
        write.csv(
          mc_ci_table,
          file.path(
            output_dir,
            paste0(
              scenario_tag,
              "_ksens_MC_CI.csv"
            )
          ),
          row.names = FALSE,
          fileEncoding = "UTF-8"
        )
      }

      write.csv(
        validation,
        file.path(
          output_dir,
          paste0(
            scenario_tag,
            "_ksens_validation.csv"
          )
        ),
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )

      key_rows[[row_id]] <-
        summary_table

      validation_rows[[row_id]] <-
        validation

      scenario_objects[[key]] <-
        summary_object

      if (
        !is.null(
          mc_ci_table
        )
      ) {
        ci_rows[[ci_id]] <-
          mc_ci_table
        ci_id <- ci_id + 1L
      }

      row_id <- row_id + 1L
    }
  }

  key_table <- do.call(
    rbind,
    key_rows
  )

  rownames(key_table) <- NULL

  validation_all <- do.call(
    rbind,
    validation_rows
  )

  rownames(validation_all) <- NULL

  mc_ci_all <- if (
    length(ci_rows) > 0L
  ) {
    z <- do.call(
      rbind,
      ci_rows
    )
    rownames(z) <- NULL
    z
  } else {
    NULL
  }

  # ------------------------------------------------------------
  # Compact k-comparison table.
  # ------------------------------------------------------------

  compact_table <- key_table[
    ,
    c(
      "h2",
      "information",
      "n_pheno",
      "k",
      "selected_fraction",
      "R_same_plugin",
      "R_cross_plugin",
      "P_alignment_plugin",
      "baseline_fraction",
      "common_alignment_fraction",
      "policy_switch_fraction",
      "common_share_of_alignment",
      "cross_plugin_over_fixed",
      "exact_set_match_rate",
      "mean_overlap",
      "mean_number_replaced"
    )
  ]

  out <- list(
    settings = list(
      S_reference =
        S,
      B_reference =
        B_reference,
      outer_indices =
        outer_indices,
      B_inner =
        B,
      k_values =
        k_values,
      n_cross_shifts =
        n_cross_shifts,
      MC_R =
        MC_R,
      estimand = paste0(
        "selection-intensity sensitivity of ",
        "same-sample selection-VC alignment"
      ),
      output_dir = normalizePath(
        output_dir,
        winslash = "/",
        mustWork = FALSE
      )
    ),
    key_table =
      key_table,
    compact_table =
      compact_table,
    mc_ci_table =
      mc_ci_all,
    validation =
      validation_all,
    result =
      scenario_objects
  )

  write.csv(
    key_table,
    file.path(
      output_dir,
      "selection_intensity_key_table.csv"
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    compact_table,
    file.path(
      output_dir,
      "selection_intensity_compact_table.csv"
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  if (!is.null(mc_ci_all)) {
    write.csv(
      mc_ci_all,
      file.path(
        output_dir,
        "selection_intensity_MC_CI_long.csv"
      ),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
  }

  write.csv(
    validation_all,
    file.path(
      output_dir,
      "selection_intensity_validation.csv"
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  saveRDS(
    out,
    file.path(
      output_dir,
      "selection_intensity_master.rds"
    )
  )

  cat(
    "\n====================================================\n",
    "Selection-intensity sensitivity completed.\n",
    "====================================================\n",
    sep = ""
  )

  cat("\nCompact table:\n")
  print(compact_table)

  cat("\nValidation:\n")
  print(validation_all)

  cat(
    "\nOutput directory:\n",
    normalizePath(
      output_dir,
      winslash = "/",
      mustWork = FALSE
    ),
    "\n",
    sep = ""
  )

  invisible(out)
}


# ----------------------------------------------------------------------
# Compact printer
# ----------------------------------------------------------------------

print_selection_intensity_sensitivity <- function(
    x,
    digits = 4) {

  if (
    is.null(x) ||
      is.null(x$compact_table)
  ) {
    stop(
      "Object has no $compact_table."
    )
  }

  z <- x$compact_table

  cat(
    "\n====================================================\n",
    "SELECTION-INTENSITY SENSITIVITY: k = ",
    paste(
      sort(unique(z$k)),
      collapse = ", "
    ),
    "\n====================================================\n",
    sep = ""
  )

  show <- data.frame(
    h2 = z$h2,
    n_pheno =
      z$n_pheno,
    k = z$k,

    R_same = round(
      z$R_same_plugin,
      digits
    ),

    R_cross = round(
      z$R_cross_plugin,
      digits
    ),

    P_alignment = round(
      z$P_alignment_plugin,
      digits
    ),

    baseline = round(
      z$baseline_fraction,
      digits
    ),

    common_alignment = round(
      z$common_alignment_fraction,
      digits
    ),

    policy_switch = round(
      z$policy_switch_fraction,
      digits
    ),

    common_share_alignment = round(
      z$common_share_of_alignment,
      digits
    ),

    cross_over_fixed = round(
      z$cross_plugin_over_fixed,
      digits
    ),

    exact_match = round(
      z$exact_set_match_rate,
      digits
    ),

    mean_overlap = round(
      z$mean_overlap,
      digits
    ),

    row.names = NULL,
    check.names = FALSE
  )

  print(show)

  cat(
    "\nWhat to look for:\n",
    "1. P_alignment should remain large across k if the mechanism is not specific\n",
    "   to selecting 5 of 20 candidates.\n",
    "2. cross_over_fixed should remain near 1 if marginal policy randomness alone\n",
    "   is still insufficient to produce the same-sample amplification.\n",
    "3. common_share_alignment shows whether common-data dependence remains the\n",
    "   dominant component as selection intensity changes.\n",
    "4. k=5 in the full S=1000, B=500 run is checked against the completed\n",
    "   two-layer reference and should reproduce it numerically.\n",
    sep = ""
  )

  invisible(show)
}


# ======================================================================
# END STEP 7
# ======================================================================

# ----------------------------------------------------------------------------
# Optional helper: translate legacy internal method labels to manuscript labels
# ----------------------------------------------------------------------------
paper_method_label <- function(x) {
  y <- as.character(x)
  y[y %in% c("Bayesian_Beta11", "Bayesian_policy")] <- "RL-UP"
  y[y %in% c("RAM_N", "RAM_N_policy")] <- "RAM-N"
  y[y %in% c("bootstrap", "bootstrap_policy")] <- "Bootstrap"
  y[y %in% c("conditional", "conditional_plugin")] <- "Conditional"
  y
}

# ======================================================================
# STANDALONE TRUE-VC REFERENCE RUNNER
# Preserved verbatim from the attached final production code.
# It keeps the realised REML-EBLUP c(Y) fixed and replaces only the
# variance components used for prediction/PEC; it does not reselect.
# ======================================================================

fixed_eta_policy_interval <- function(
    eta,
    y,
    prep,
    cstar,
    target = 0) {

  g <- policy_stats_animal(
    eta = eta,
    y = y,
    prep = prep,
    cstar = cstar,
    calc_diag = FALSE
  )

  sd_use <- max(g$sd, 1e-12)

  list(
    valid = TRUE,
    mean = as.numeric(g$mu),
    sd = as.numeric(g$sd),
    interval = as.numeric(
      g$mu + c(-1, 1) * 1.96 * g$sd
    ),
    p = as.numeric(
      1 - pnorm(
        target,
        mean = g$mu,
        sd = sd_use
      )
    )
  )
}

one_true_VC_oracle_scenario <- function(
    latent,
    info_bundle,
    info_level,
    h2,
    sigma2_P = 1,
    beta = 0,
    target = 0,
    n_select = 5) {

  base <- info_bundle$base_design

  # Reconstruct the same dataset as the production Simulation II.
  dat <- make_dataset_from_latent(
    latent = latent,
    info_level = info_level,
    sigma2_P = sigma2_P,
    h2 = h2,
    beta = beta
  )

  # ----------------------------------------------------------
  # Step 1: REML-EBLUP selection exactly as in Simulation II
  # ----------------------------------------------------------

  fit_sel <- fit_reml_animal(
    y = dat$y,
    prep = info_level$prep
  )

  pred_sel <- predict_ebv_validation(
    eta = fit_sel$eta,
    y = dat$y,
    prep = info_level$prep
  )

  candidates <- base$candidates

  ord <- order(
    pred_sel$uhat[candidates],
    decreasing = TRUE
  )

  selected <- candidates[
    ord[seq_len(n_select)]
  ]

  cstar <- make_policy(
    n_animals = nrow(base$ped),
    candidates = candidates,
    selected = selected
  )

  # Same estimand as the production realised-policy analysis:
  # true expected genetic superiority of this data-selected policy.
  truth <- sum(cstar * dat$u)

  # ----------------------------------------------------------
  # Step 2: REML plug-in conditional reference
  # ----------------------------------------------------------
  # Use the SAME REML fit that generated the selected set.

  cond <- fixed_eta_policy_interval(
    eta = fit_sel$eta,
    y = dat$y,
    prep = info_level$prep,
    cstar = cstar,
    target = target
  )

  # ----------------------------------------------------------
  # Step 3: TRUE-VC oracle
  # ----------------------------------------------------------
  # Crucial: cstar is NOT changed and candidates are NOT reselected.

  eta_true <- log(
    c(
      dat$sigma2_A,
      dat$sigma2_e
    )
  )

  oracle <- fixed_eta_policy_interval(
    eta = eta_true,
    y = dat$y,
    prep = info_level$prep,
    cstar = cstar,
    target = target
  )

  h2_hat <- fit_sel$theta[1] / sum(fit_sel$theta)

  list(
    truth = as.numeric(truth),
    h2 = h2,
    information = info_level$name,
    n_pheno = info_level$n_pheno,

    selected = selected,
    cstar = cstar,

    true_sigma2_A = dat$sigma2_A,
    true_sigma2_e = dat$sigma2_e,
    h2_hat = as.numeric(h2_hat),
    sigma2A_hat = as.numeric(fit_sel$theta[1]),
    sigma2e_hat = as.numeric(fit_sel$theta[2]),

    conditional = pack_validation_method(cond),
    true_VC_oracle = pack_validation_method(oracle),

    conditional_error = as.numeric(cond$mean - truth),
    oracle_error = as.numeric(oracle$mean - truth),
    point_shift_oracle_minus_plugin = as.numeric(
      oracle$mean - cond$mean
    ),
    SE_ratio_oracle_over_plugin = if (
      is.finite(cond$sd) && cond$sd > 0
    ) {
      as.numeric(oracle$sd / cond$sd)
    } else {
      NA_real_
    },

    reml_diag = extract_reml_diagnostics(
      fit_sel,
      true_h2 = h2
    )
  )
}

summarize_true_VC_oracle_scenario <- function(
    res,
    target = 0) {

  cond <- summarize_simII(
    res = res,
    name = "conditional",
    target = target
  )

  oracle <- summarize_simII(
    res = res,
    name = "true_VC_oracle",
    target = target
  )

  method_table <- rbind(
    conditional = cond,
    true_VC_oracle = oracle
  )

  method_table <- as.data.frame(
    method_table,
    stringsAsFactors = FALSE
  )

  method_table$method <- rownames(method_table)
  rownames(method_table) <- NULL

  # Put method first.
  method_table <- method_table[
    , c(
      "method",
      setdiff(names(method_table), "method")
    ),
    drop = FALSE
  ]

  get_num <- function(name) {
    vapply(
      res,
      `[[`,
      numeric(1),
      name
    )
  }

  D <- do.call(
    rbind,
    lapply(
      res,
      `[[`,
      "reml_diag"
    )
  )

  comparison <- data.frame(
    conditional_coverage = unname(cond["inclusion"]),
    oracle_coverage = unname(oracle["inclusion"]),
    coverage_improvement = unname(
      oracle["inclusion"] - cond["inclusion"]
    ),

    conditional_MSE_ratio = unname(cond["MSE_ratio"]),
    oracle_MSE_ratio = unname(oracle["MSE_ratio"]),

    conditional_RMSE = unname(cond["RMSE"]),
    oracle_RMSE = unname(oracle["RMSE"]),

    conditional_mean_SE = unname(cond["mean_SE"]),
    oracle_mean_SE = unname(oracle["mean_SE"]),

    mean_point_shift_oracle_minus_plugin = mean(
      get_num("point_shift_oracle_minus_plugin")
    ),
    RMSE_point_shift_oracle_minus_plugin = sqrt(
      mean(
        get_num("point_shift_oracle_minus_plugin")^2
      )
    ),
    mean_SE_ratio_oracle_over_plugin = mean(
      get_num("SE_ratio_oracle_over_plugin"),
      na.rm = TRUE
    ),

    mean_h2_hat = mean(D[, "h2_hat"]),
    median_h2_hat = median(D[, "h2_hat"]),
    lower_boundary_rate = mean(D[, "near_lower_A"]),
    stringsAsFactors = FALSE
  )

  list(
    methods = method_table,
    comparison = comparison
  )
}

run_true_VC_oracle_final <- function(
    S = 1000,
    h2_values = c(0.05, 0.20, 0.40),
    info_fractions = c(
      low = 0.40,
      moderate = 0.70,
      high = 1.00
    ),
    n_select = 5,
    sigma2_P = 1,
    beta = 0,
    target = 0,
    seed_info = 20260850,
    seed_outer = 20260851,
    checkpoint_every = 25,
    output_dir = file.path(
      getwd(),
      "simulation_true_VC_oracle_output"
    ),
    resume = TRUE,
    keep_raw_in_master = FALSE) {

  final_dir_create(output_dir)

  base_design <- build_simII_design()

  info_bundle <- build_information_level_designs(
    base_design = base_design,
    info_fractions = info_fractions,
    seed_info = seed_info
  )

  cat(
    "\n====================================================\n",
    "Simulation II: TRUE-VC ORACLE diagnostic\n",
    "Same c(Y), no reselection under true VC\n",
    "====================================================\n",
    sep = ""
  )

  cat("\nCandidate relationship summary:\n")
  print(
    relationship_summary(
      base_design$A,
      base_design$candidates
    )
  )

  print_information_design(info_bundle)

  # IMPORTANT: same outer seed construction as production Simulation II.
  set.seed(seed_outer)

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

  scenario_results <- list()

  for (ih in seq_along(h2_values)) {

    h2 <- h2_values[ih]

    for (ii in seq_along(info_bundle$info_levels)) {

      info <- info_bundle$info_levels[[ii]]

      key <- paste0(
        "h2_",
        h2,
        "__",
        info$name
      )

      cat(
        "\n====================================================\n",
        "TRUE-VC oracle: ",
        key,
        " | n_pheno=",
        info$n_pheno,
        "\n====================================================\n",
        sep = ""
      )

      checkpoint_file <- file.path(
        output_dir,
        paste0(
          final_safe_name(key),
          "_oracle_checkpoint.rds"
        )
      )

      final_file <- file.path(
        output_dir,
        paste0(
          final_safe_name(key),
          "_oracle_final.rds"
        )
      )

      if (
        resume &&
        file.exists(final_file)
      ) {
        cat("Loading completed oracle scenario.\n")

        done <- readRDS(final_file)
        scenario_results[[key]] <- done$summary

        if (
          keep_raw_in_master &&
          !is.null(done$raw)
        ) {
          scenario_results[[key]]$raw <- done$raw
        }

        next
      }

      res <- vector("list", S)
      start_at <- 1L

      if (
        resume &&
        file.exists(checkpoint_file)
      ) {

        cp <- readRDS(checkpoint_file)

        compatible <-
          identical(cp$S, S) &&
          identical(cp$key, key) &&
          identical(cp$seed_info, seed_info) &&
          identical(cp$seed_outer, seed_outer)

        if (compatible) {

          res <- cp$res

          missing_idx <- which(
            vapply(
              res,
              is.null,
              logical(1)
            )
          )

          if (length(missing_idx) == 0L) {
            start_at <- S + 1L
          } else {
            start_at <- min(missing_idx)
          }

          cat(
            "Resuming at outer replicate ",
            start_at,
            ".\n",
            sep = ""
          )
        }
      }

      if (start_at <= S) {

        for (s in seq.int(start_at, S)) {

          res[[s]] <- one_true_VC_oracle_scenario(
            latent = latent_list[[s]],
            info_bundle = info_bundle,
            info_level = info,
            h2 = h2,
            sigma2_P = sigma2_P,
            beta = beta,
            target = target,
            n_select = n_select
          )

          if (
            s %% checkpoint_every == 0 ||
            s == S
          ) {

            saveRDS(
              list(
                S = S,
                key = key,
                h2 = h2,
                information = info$name,
                seed_info = seed_info,
                seed_outer = seed_outer,
                res = res
              ),
              checkpoint_file
            )

            cat(
              "Completed ",
              s,
              " / ",
              S,
              "\n",
              sep = ""
            )
          }
        }
      }

      sm <- summarize_true_VC_oracle_scenario(
        res = res,
        target = target
      )

      sm$methods$h2 <- h2
      sm$methods$information <- info$name
      sm$methods$n_pheno <- info$n_pheno

      sm$methods <- sm$methods[
        , c(
          "h2",
          "information",
          "n_pheno",
          "method",
          setdiff(
            names(sm$methods),
            c(
              "h2",
              "information",
              "n_pheno",
              "method"
            )
          )
        ),
        drop = FALSE
      ]

      sm$comparison$h2 <- h2
      sm$comparison$information <- info$name
      sm$comparison$n_pheno <- info$n_pheno

      sm$comparison <- sm$comparison[
        , c(
          "h2",
          "information",
          "n_pheno",
          setdiff(
            names(sm$comparison),
            c(
              "h2",
              "information",
              "n_pheno"
            )
          )
        ),
        drop = FALSE
      ]

      summary_one <- list(
        h2 = h2,
        information = info$name,
        n_pheno = info$n_pheno,
        methods = sm$methods,
        comparison = sm$comparison
      )

      saveRDS(
        list(
          S = S,
          key = key,
          summary = summary_one,
          raw = res
        ),
        final_file
      )

      scenario_results[[key]] <- summary_one

      if (keep_raw_in_master) {
        scenario_results[[key]]$raw <- res
      }
    }
  }

  method_table <- do.call(
    rbind,
    lapply(
      scenario_results,
      `[[`,
      "methods"
    )
  )

  rownames(method_table) <- NULL
  method_table <- final_add_coverage_mcse(method_table)

  comparison_table <- do.call(
    rbind,
    lapply(
      scenario_results,
      `[[`,
      "comparison"
    )
  )

  rownames(comparison_table) <- NULL

  final_write_csv(
    method_table,
    file.path(
      output_dir,
      "Simulation_II_true_VC_oracle_method_performance.csv"
    )
  )

  final_write_csv(
    comparison_table,
    file.path(
      output_dir,
      "Simulation_II_true_VC_oracle_comparison.csv"
    )
  )

  out <- list(
    design = base_design,
    information_design = info_bundle,
    result = scenario_results,
    method_table = method_table,
    comparison_table = comparison_table,
    settings = list(
      S = S,
      h2_values = h2_values,
      info_fractions = info_fractions,
      n_select = n_select,
      sigma2_P = sigma2_P,
      beta = beta,
      target = target,
      seed_info = seed_info,
      seed_outer = seed_outer,
      output_dir = normalizePath(
        output_dir,
        winslash = "/",
        mustWork = FALSE
      ),
      estimand =
        "true expected genetic superiority c(Y)'u for the REML-EBLUP-selected policy; c(Y) fixed for oracle calculation",
      oracle =
        "same c(Y), true variance components, no reselection"
    )
  )

  saveRDS(
    out,
    file.path(
      output_dir,
      "Simulation_II_true_VC_oracle_ALL.rds"
    )
  )

  cat("\nTRUE-VC oracle simulation completed.\n")
  cat("\nMethod performance:\n")
  print(method_table)
  cat("\nOracle comparison:\n")
  print(comparison_table)

  invisible(out)
}

verify_true_VC_oracle_against_previous <- function(
    oracle_output,
    previous_dir,
    oracle_dir = oracle_output$settings$output_dir,
    tolerance = 1e-10) {

  if (!dir.exists(previous_dir)) {
    stop("previous_dir does not exist: ", previous_dir)
  }

  if (is.null(oracle_dir) || !dir.exists(oracle_dir)) {
    stop(
      "oracle_dir could not be resolved. Supply oracle_dir explicitly."
    )
  }

  rows <- list()
  k <- 1L

  for (key in names(oracle_output$result)) {

    old_file <- file.path(
      previous_dir,
      paste0(
        final_safe_name(key),
        "_final.rds"
      )
    )

    oracle_file <- file.path(
      oracle_dir,
      paste0(
        final_safe_name(key),
        "_oracle_final.rds"
      )
    )

    # Prefer raw attached to the master output when available; otherwise
    # load the scenario-level oracle RDS.
    new_raw <- oracle_output$result[[key]]$raw

    if (is.null(new_raw) && file.exists(oracle_file)) {
      new_raw <- readRDS(oracle_file)$raw
    }

    if (is.null(new_raw)) {
      rows[[k]] <- data.frame(
        scenario = key,
        previous_file_found = file.exists(old_file),
        oracle_raw_found = FALSE,
        selected_match_rate = NA_real_,
        max_abs_truth_difference = NA_real_,
        truth_match_within_tolerance = NA,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
      next
    }

    if (!file.exists(old_file)) {
      rows[[k]] <- data.frame(
        scenario = key,
        previous_file_found = FALSE,
        oracle_raw_found = TRUE,
        selected_match_rate = NA_real_,
        max_abs_truth_difference = NA_real_,
        truth_match_within_tolerance = NA,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
      next
    }

    old_raw <- readRDS(old_file)$raw

    if (is.null(old_raw)) {
      rows[[k]] <- data.frame(
        scenario = key,
        previous_file_found = TRUE,
        oracle_raw_found = TRUE,
        selected_match_rate = NA_real_,
        max_abs_truth_difference = NA_real_,
        truth_match_within_tolerance = NA,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
      next
    }

    n <- min(length(old_raw), length(new_raw))

    selected_match <- vapply(
      seq_len(n),
      function(i) {
        identical(
          as.integer(old_raw[[i]]$selected),
          as.integer(new_raw[[i]]$selected)
        )
      },
      logical(1)
    )

    truth_diff <- vapply(
      seq_len(n),
      function(i) {
        abs(
          old_raw[[i]]$truth -
          new_raw[[i]]$truth
        )
      },
      numeric(1)
    )

    rows[[k]] <- data.frame(
      scenario = key,
      previous_file_found = TRUE,
      oracle_raw_found = TRUE,
      selected_match_rate = mean(selected_match),
      max_abs_truth_difference = max(truth_diff),
      truth_match_within_tolerance = all(
        truth_diff <= tolerance
      ),
      stringsAsFactors = FALSE
    )

    k <- k + 1L
  }

  do.call(rbind, rows)
}

# ======================================================================
# END STANDALONE TRUE-VC REFERENCE RUNNER
# ======================================================================

