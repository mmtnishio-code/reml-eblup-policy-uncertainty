# ============================================================================
# GSE: FIGURE REPRODUCTION FROM ANALYSIS OUTPUTS
# Base R only.  No manuscript values are hard-coded into plotting functions.
# ============================================================================

.gse_assert_columns <- function(x, required, label) {
  miss <- setdiff(required, names(x))
  if (length(miss) > 0L) {
    stop(label, " is missing columns: ", paste(miss, collapse = ", "))
  }
  invisible(TRUE)
}

.gse_read_csv <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path)
  z <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(z) == 0L) stop(label, " is empty: ", path)
  z
}

.gse_scenario_label <- function(h2, n_pheno) {
  paste0("h2=", format(h2, trim = TRUE), "\nn=", n_pheno)
}

make_GSE_coverage_figure <- function(
    simulation_II_csv = file.path(
      getwd(), "simulation_I_II_final_output", "Simulation_II",
      "Simulation_II_method_performance.csv"
    ),
    trueVC_csv = file.path(
      getwd(), "bootstrap_diagnostic",
      "Simulation_II_bootstrap_generator_method_performance.csv"
    ),
    output_pdf = file.path(getwd(), "Figure_coverage_comparison.pdf"),
    output_png = file.path(getwd(), "Figure_coverage_comparison.png"),
    width = 9,
    height = 5.5) {

  main <- .gse_read_csv(simulation_II_csv, "Simulation-II method table")
  tv <- .gse_read_csv(trueVC_csv, "True-VC/bootstrap diagnostic method table")

  .gse_assert_columns(main, c("h2", "n_pheno", "method", "inclusion"), "Simulation-II method table")
  .gse_assert_columns(tv, c("h2", "n_pheno", "method", "inclusion"), "True-VC table")

  keep <- c("conditional", "RAM_N_policy", "Bayesian_policy", "bootstrap_policy")
  z <- main[main$method %in% keep, c("h2", "n_pheno", "method", "inclusion"), drop = FALSE]
  truevc_method <- if ("conditional_trueVC_oracle" %in% tv$method) {
    "conditional_trueVC_oracle"
  } else if ("true_VC_oracle" %in% tv$method) {
    "true_VC_oracle"
  } else {
    stop("True-VC input has neither conditional_trueVC_oracle nor true_VC_oracle method rows.")
  }
  tvz <- tv[tv$method == truevc_method, c("h2", "n_pheno", "method", "inclusion"), drop = FALSE]
  tvz$method <- "TrueVC_reference"
  z <- rbind(z, tvz)

  map <- c(
    conditional = "Conditional",
    RAM_N_policy = "RAM-N",
    Bayesian_policy = "RL-UP",
    bootstrap_policy = "Fixed-policy bootstrap",
    TrueVC_reference = "True-VC reference"
  )
  z$label <- unname(map[z$method])
  if (any(is.na(z$label))) stop("Unexpected method name in coverage inputs.")

  scen <- unique(z[, c("h2", "n_pheno")])
  scen <- scen[order(scen$h2, scen$n_pheno), , drop = FALSE]
  scen$key <- paste(scen$h2, scen$n_pheno, sep = "__")
  z$key <- paste(z$h2, z$n_pheno, sep = "__")
  z$scenario_index <- match(z$key, scen$key)

  method_order <- unname(map[c("conditional", "RAM_N_policy", "Bayesian_policy", "bootstrap_policy", "TrueVC_reference")])
  mat <- matrix(NA_real_, nrow = nrow(scen), ncol = length(method_order),
                dimnames = list(NULL, method_order))
  for (j in seq_along(method_order)) {
    q <- z[z$label == method_order[j], , drop = FALSE]
    mat[q$scenario_index, j] <- q$inclusion
  }
  if (any(!is.finite(mat))) stop("Coverage figure matrix contains missing/non-finite values.")

  draw <- function() {
    x <- seq_len(nrow(scen))
    plot(x, mat[, 1], type = "n", ylim = range(c(mat, 0.95), finite = TRUE),
         xaxt = "n", xlab = "Simulation-II scenario", ylab = "95% interval coverage")
    axis(1, at = x, labels = .gse_scenario_label(scen$h2, scen$n_pheno), las = 2, cex.axis = 0.8)
    abline(h = 0.95, lty = 2)
    for (j in seq_len(ncol(mat))) {
      lines(x, mat[, j], type = "b", pch = j, lty = j)
    }
    legend("bottomright", legend = colnames(mat), pch = seq_len(ncol(mat)),
           lty = seq_len(ncol(mat)), bty = "n", cex = 0.85)
  }

  pdf(output_pdf, width = width, height = height)
  draw()
  dev.off()
  png(output_png, width = width, height = height, units = "in", res = 180)
  draw()
  dev.off()

  invisible(list(data = z, matrix = mat, pdf = output_pdf, png = output_png))
}

make_GSE_same_cross_figure <- function(
    same_cross_csv = file.path(getwd(), "same_cross", "same_vs_cross_selection_VC_key_table.csv"),
    output_pdf = file.path(getwd(), "Figure_same_cross.pdf"),
    output_png = file.path(getwd(), "Figure_same_cross.png"),
    width = 9,
    height = 5.5) {

  z <- .gse_read_csv(same_cross_csv, "same/cross key table")
  .gse_assert_columns(
    z,
    c("h2", "n_pheno", "R_VC_fixed", "R_VC_same_plugin", "R_VC_cross_plugin"),
    "same/cross key table"
  )
  z <- z[order(z$h2, z$n_pheno), , drop = FALSE]
  mat <- cbind(
    `Fixed policy` = z$R_VC_fixed,
    `Same sample` = z$R_VC_same_plugin,
    `Cross sample` = z$R_VC_cross_plugin
  )
  if (any(!is.finite(mat))) stop("same/cross figure input contains non-finite values.")

  draw <- function() {
    x <- seq_len(nrow(z))
    plot(x, mat[, 1], type = "n", ylim = range(c(0, mat), finite = TRUE), xaxt = "n",
         xlab = "Simulation-II scenario",
         ylab = "Ratio to outer VC-error second moment")
    axis(1, at = x, labels = .gse_scenario_label(z$h2, z$n_pheno), las = 2, cex.axis = 0.8)
    abline(h = 1, lty = 2)
    for (j in seq_len(ncol(mat))) lines(x, mat[, j], type = "b", pch = j, lty = j)
    legend("topright", legend = colnames(mat), pch = seq_len(ncol(mat)),
           lty = seq_len(ncol(mat)), bty = "n")
  }

  pdf(output_pdf, width = width, height = height); draw(); dev.off()
  png(output_png, width = width, height = height, units = "in", res = 180); draw(); dev.off()
  invisible(list(data = z, matrix = mat, pdf = output_pdf, png = output_png))
}

make_GSE_staged_decomposition_figure <- function(
    staged_csv = file.path(getwd(), "staged", "two_layer_alignment_key_table.csv"),
    output_pdf = file.path(getwd(), "Figure_staged_decomposition.pdf"),
    output_png = file.path(getwd(), "Figure_staged_decomposition.png"),
    width = 9,
    height = 5.5) {

  z <- .gse_read_csv(staged_csv, "staged-decomposition key table")
  .gse_assert_columns(
    z,
    c("h2", "n_pheno", "baseline_fraction_of_plugin",
      "common_alignment_fraction_of_plugin", "policy_switch_fraction_of_plugin"),
    "staged-decomposition key table"
  )
  z <- z[order(z$h2, z$n_pheno), , drop = FALSE]
  mat <- rbind(
    `Fixed/independent baseline` = z$baseline_fraction_of_plugin,
    `Common-data increment` = z$common_alignment_fraction_of_plugin,
    `Policy-switch increment` = z$policy_switch_fraction_of_plugin
  )
  if (any(!is.finite(mat))) stop("Staged figure input contains non-finite values.")

  # These are sequential net increments, not orthogonal variance components.
  draw <- function() {
    barplot(mat, beside = FALSE,
            names.arg = .gse_scenario_label(z$h2, z$n_pheno), las = 2,
            ylab = "Sequential second-moment fraction", xlab = "Simulation-II scenario",
            ylim = c(0, max(1, colSums(mat))))
    abline(h = 1, lty = 2)
    legend("topright", legend = rownames(mat), fill = seq_len(nrow(mat)), bty = "n", cex = 0.85)
  }

  pdf(output_pdf, width = width, height = height); draw(); dev.off()
  png(output_png, width = width, height = height, units = "in", res = 180); draw(); dev.off()
  invisible(list(data = z, matrix = mat, pdf = output_pdf, png = output_png))
}


make_GSE_staged_average_waterfall_figure <- function(
    staged_csv = file.path(getwd(), "staged", "two_layer_alignment_key_table.csv"),
    output_pdf = file.path(getwd(), "Figure_3_staged_average_waterfall.pdf"),
    output_png = file.path(getwd(), "Figure_3_staged_average_waterfall.png"),
    width = 7.5,
    height = 5.5) {

  z <- .gse_read_csv(staged_csv, "staged-decomposition key table")
  required <- c(
    "baseline_fraction_of_plugin",
    "common_alignment_fraction_of_plugin",
    "policy_switch_fraction_of_plugin"
  )
  .gse_assert_columns(z, required, "staged-decomposition key table")

  x <- z[, required, drop = FALSE]
  if (any(!is.finite(as.matrix(x)))) {
    stop("Staged-average figure input contains missing/non-finite values.")
  }

  # Manuscript Figure 3 reports the arithmetic mean across the available
  # Simulation-II scenarios. These are sequential net increments from the
  # staged replacement diagnostic; they are not orthogonal variance components.
  values <- colMeans(x)
  names(values) <- c(
    "Fixed-policy stage",
    "Common-data increment",
    "Policy-switch increment"
  )
  cumulative <- cumsum(values)
  starts <- c(0, head(cumulative, -1L))
  ends <- cumulative

  summary_table <- data.frame(
    component = names(values),
    value = unname(values),
    cumulative = unname(cumulative),
    stringsAsFactors = FALSE
  )

  draw <- function() {
    ylim <- range(c(0, 1, starts, ends), finite = TRUE)
    pad <- max(0.04, 0.08 * diff(ylim))
    ylim <- c(ylim[1] - pad, ylim[2] + pad)

    plot(NA_real_, NA_real_,
         xlim = c(0.45, 3.55), ylim = ylim, xaxt = "n",
         xlab = "Sequential diagnostic stage",
         ylab = "Cumulative fraction of same-sample second moment")
    axis(1, at = 1:3, labels = c(
      "Fixed-policy\nstage",
      "+ Common-data\nincrement",
      "+ Policy-switch\nincrement"
    ))
    abline(h = 1, lty = 2)

    fill <- gray(c(0.85, 0.65, 0.45))
    for (j in 1:3) {
      rect(j - 0.32, starts[j], j + 0.32, ends[j],
           col = fill[j], border = "black")
      if (j < 3L) {
        segments(j + 0.32, ends[j], j + 1 - 0.32, ends[j], lty = 3)
      }
      text(j, (starts[j] + ends[j]) / 2,
           labels = sprintf("%.3f", values[j]))
    }
    text(3.45, ends[3], labels = sprintf("Total = %.3f", ends[3]),
         adj = c(1, -0.25))
  }

  pdf(output_pdf, width = width, height = height)
  draw()
  dev.off()
  png(output_png, width = width, height = height, units = "in", res = 180)
  draw()
  dev.off()

  invisible(list(
    data = z,
    summary = summary_table,
    values = values,
    cumulative = cumulative,
    pdf = output_pdf,
    png = output_png
  ))
}

# Example after production analyses:
# make_GSE_coverage_figure(
#   simulation_II_csv = "simulation_I_II_final_output/Simulation_II/Simulation_II_method_performance.csv",
#   trueVC_csv = "bootstrap_diagnostic/Simulation_II_bootstrap_generator_method_performance.csv"
# )
# make_GSE_same_cross_figure("same_cross/same_vs_cross_selection_VC_key_table.csv")
# make_GSE_staged_decomposition_figure("staged/two_layer_alignment_key_table.csv")
