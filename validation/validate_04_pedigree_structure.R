# ============================================================================
# Validation for 04_GSE_pedigree_structure_sensitivity.R
# Base R only.
# ============================================================================

.this <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.vdir <- if (!is.null(.this)) dirname(normalizePath(.this, winslash = "/", mustWork = FALSE)) else file.path(getwd(), "validation")
.root <- normalizePath(file.path(.vdir, ".."), winslash = "/", mustWork = FALSE)

run_GSE_pedigree_design_validation <- function() {
  GSE_PEDIGREE_RUN_MODE <<- "design_only"
  source(file.path(.root, "04_GSE_pedigree_structure_sensitivity.R"),
         echo = FALSE, chdir = FALSE, encoding = "UTF-8", local = .GlobalEnv)

  exp <- pedigree_structure_expected_results()
  fields <- c(
    "candidate_mean_relationship",
    "candidate_q90_relationship",
    "candidate_max_relationship",
    "candidate_mean_inbreeding",
    "mean_candidate_to_record_relationship"
  )
  obs <- relationship_table[match(exp$structure, relationship_table$structure), c("structure", fields), drop = FALSE]
  diffs <- abs(as.matrix(obs[, fields, drop = FALSE]) - as.matrix(exp[, fields, drop = FALSE]))

  # The archived design targets above are stored to six decimal places.
  # A value exactly half a unit in the sixth decimal place (for example,
  # 0.0296875 versus archived 0.029688) can be represented in binary
  # floating point as infinitesimally larger than 5e-7.  Allow only that
  # machine-roundoff margin; this changes validation tolerance only and does
  # not alter pedigree generation or any scientific calculation.
  archived_rounding_tolerance <- 0.5e-6 + 100 * .Machine$double.eps

  if (any(!is.finite(diffs)) || any(diffs > archived_rounding_tolerance)) {
    print(obs, row.names = FALSE)
    cat("Maximum absolute design-summary difference: ",
        format(max(diffs, na.rm = TRUE), digits = 17), "\n", sep = "")
    stop("Pedigree design does not reproduce the archived relationship summaries.")
  }
  if (!all(info_check$n_pheno == 90L) || !all(info_check$candidate_records == 20L)) {
    stop("Pedigree information-level design check failed.")
  }
  cat("Pedigree design-only validation: PASS\n")
  invisible(obs)
}

validate_GSE_pedigree_output <- function(path) {
  GSE_PEDIGREE_RUN_MODE <<- "design_only"
  source(file.path(.root, "04_GSE_pedigree_structure_sensitivity.R"),
         echo = FALSE, chdir = FALSE, encoding = "UTF-8", local = .GlobalEnv)
  validate_pedigree_structure_production(path, stop_on_failure = TRUE)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L || identical(args[[1]], "design")) {
    run_GSE_pedigree_design_validation()
  } else if (identical(args[[1]], "validate")) {
    if (length(args) < 2L) stop("validate mode requires pedigree_structure_sensitivity_FINAL_comparison.csv")
    z <- validate_GSE_pedigree_output(args[[2]])
    print(z, row.names = FALSE)
  } else if (identical(args[[1]], "pilot")) {
    options(gse.pedigree.pilot.S = 2L, gse.pedigree.pilot.B = 20L)
    options(gse.output.root = file.path(tempdir(), "gse_pedigree_pilot"))
    GSE_PEDIGREE_RUN_MODE <- "pilot"
    source(file.path(.root, "04_GSE_pedigree_structure_sensitivity.R"),
           echo = FALSE, chdir = FALSE, encoding = "UTF-8")
    cat("Pedigree S=2, B=20 pilot completed. Archived-value validation is intentionally not applied to a pilot.\n")
  } else {
    stop("Use design, validate, or pilot.")
  }
}

# Optional cross-analysis check: baseline pedigree sensitivity should reproduce
# the h2=0.20, n=90 baseline main/True-VC summaries when all analyses use the
# archived settings and RNG streams.
validate_GSE_pedigree_baseline_against_main <- function(
    pedigree_csv,
    main_method_csv,
    trueVC_method_csv,
    tolerance = 1e-8,
    stop_on_failure = TRUE) {

  ped <- read.csv(pedigree_csv, stringsAsFactors = FALSE, check.names = FALSE)
  main <- read.csv(main_method_csv, stringsAsFactors = FALSE, check.names = FALSE)
  tv <- read.csv(trueVC_method_csv, stringsAsFactors = FALSE, check.names = FALSE)

  p <- ped[ped$structure == "baseline", , drop = FALSE]
  m <- main[main$h2 == 0.20 & main$n_pheno == 90 & main$method == "conditional", , drop = FALSE]
  t <- tv[tv$h2 == 0.20 & tv$n_pheno == 90 & tv$method %in% c("true_VC_oracle", "conditional_trueVC_oracle"), , drop = FALSE]

  if (nrow(p) != 1L || nrow(m) != 1L || nrow(t) != 1L) {
    stop("Could not isolate unique baseline/main/True-VC rows for h2=0.20, n=90.")
  }

  d <- c(
    conditional_coverage = abs(p$conditional_coverage - m$inclusion),
    conditional_MSE_ratio = abs(p$conditional_MSE_ratio - m$MSE_ratio),
    trueVC_coverage = abs(p$trueVC_coverage - t$inclusion),
    trueVC_MSE_ratio = abs(p$trueVC_MSE_ratio - t$MSE_ratio)
  )
  ok <- all(is.finite(d)) && all(d <= tolerance)

  if (isTRUE(stop_on_failure) && !ok) {
    print(d)
    stop("Baseline pedigree sensitivity does not exactly match the corresponding main/True-VC production rows.")
  }

  cat("Baseline pedigree vs main/True-VC cross-analysis check: ", if (ok) "PASS" else "CHECK REQUIRED", "\n", sep = "")
  invisible(d)
}
