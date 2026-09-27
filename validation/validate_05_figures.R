# ============================================================================
# Smoke test for 05_GSE_make_figures.R using synthetic analysis-shaped CSVs.
# This checks file/column plumbing only; it does not validate manuscript values.
# ============================================================================

.this <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.vdir <- if (!is.null(.this)) dirname(normalizePath(.this, winslash = "/", mustWork = FALSE)) else file.path(getwd(), "validation")
.root <- normalizePath(file.path(.vdir, ".."), winslash = "/", mustWork = FALSE)
source(file.path(.root, "05_GSE_make_figures.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")

run_GSE_figure_smoke_test <- function() {
  td <- file.path(tempdir(), "gse_figure_smoke")
  unlink(td, recursive = TRUE, force = TRUE)
  dir.create(td, recursive = TRUE, showWarnings = FALSE)

  scen <- data.frame(
    h2 = rep(c(0.05, 0.20, 0.40), each = 3),
    n_pheno = rep(c(60, 90, 120), times = 3)
  )
  meth <- c("conditional", "RAM_N_policy", "Bayesian_policy", "bootstrap_policy")
  main <- do.call(rbind, lapply(seq_len(nrow(scen)), function(i) {
    data.frame(h2 = scen$h2[i], n_pheno = scen$n_pheno[i], method = meth,
               inclusion = pmin(0.99, c(0.5, 0.9, 0.94, 0.8) + i / 1000))
  }))
  tv <- data.frame(h2 = scen$h2, n_pheno = scen$n_pheno,
                   method = "true_VC_oracle", inclusion = 0.95)
  sc <- data.frame(h2 = scen$h2, n_pheno = scen$n_pheno,
                   R_VC_fixed = 0.12, R_VC_same_plugin = 0.98, R_VC_cross_plugin = 0.12)
  st <- data.frame(h2 = scen$h2, n_pheno = scen$n_pheno,
                   baseline_fraction_of_plugin = 0.12,
                   common_alignment_fraction_of_plugin = 0.81,
                   policy_switch_fraction_of_plugin = 0.07)

  p_main <- file.path(td, "main.csv"); p_tv <- file.path(td, "truevc.csv")
  p_sc <- file.path(td, "samecross.csv"); p_st <- file.path(td, "staged.csv")
  write.csv(main, p_main, row.names = FALSE); write.csv(tv, p_tv, row.names = FALSE)
  write.csv(sc, p_sc, row.names = FALSE); write.csv(st, p_st, row.names = FALSE)

  a <- make_GSE_coverage_figure(p_main, p_tv,
        file.path(td, "coverage.pdf"), file.path(td, "coverage.png"))
  b <- make_GSE_same_cross_figure(p_sc,
        file.path(td, "samecross.pdf"), file.path(td, "samecross.png"))
  c <- make_GSE_staged_decomposition_figure(p_st,
        file.path(td, "staged.pdf"), file.path(td, "staged.png"))
  d <- make_GSE_staged_average_waterfall_figure(p_st,
        file.path(td, "staged_average.pdf"), file.path(td, "staged_average.png"))

  files <- c(a$pdf, a$png, b$pdf, b$png, c$pdf, c$png, d$pdf, d$png)
  stopifnot(all(file.exists(files)), all(file.info(files)$size > 0))
  cat("Figure-generation smoke test: PASS\n")
  invisible(files)
}

if (sys.nframe() == 0L) run_GSE_figure_smoke_test()
