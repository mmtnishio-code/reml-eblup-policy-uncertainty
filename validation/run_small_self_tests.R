# ============================================================================
# End-to-end small self-tests for the GSE public release.
# Base R only.
#
# Note: the validated mechanism-diagnostic functions require B_inner >= 20, so
# this release uses B=20 for smoke tests even though smaller B would be faster.
# ============================================================================

.this <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.vdir <- if (!is.null(.this)) dirname(normalizePath(.this, winslash = "/", mustWork = FALSE)) else file.path(getwd(), "validation")
.root <- normalizePath(file.path(.vdir, ".."), winslash = "/", mustWork = FALSE)

source(file.path(.root, "01_GSE_main_and_diagnostics.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")
source(file.path(.root, "06_GSE_procedure_level_reselection_bootstrap.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")

cat("\n[1/8] Integrated Simulation-I/II self-test\n")
self_test_final_simulations()

cat("\n[2/8] Bootstrap-generator self-test\n")
self_test_bootstrap_generator_diagnostic()

cat("\n[3/8] True-VC inner-decomposition self-test\n")
self_test_outer_vs_trueVC_inner_decomposition()

cat("\n[4/8] same/cross -> staged -> k sensitivity smoke test\n")
td <- file.path(tempdir(), "gse_mechanism_smoke")
unlink(td, recursive = TRUE, force = TRUE)
boot <- run_bootstrap_generator_diagnostic_final(
  S = 2,
  h2_values = 0.20,
  info_fractions = c(moderate = 0.70),
  B = 20,
  seed_info = 20260850,
  seed_outer = 20260851,
  checkpoint_every = 1,
  output_dir = file.path(td, "boot"),
  resume = FALSE,
  keep_raw_in_master = FALSE,
  return_inner = FALSE
)
rsel <- run_inner_reselection_VC_diagnostic(
  boot1000 = boot,
  outer_indices = 1:2,
  B_inner = 20,
  prior_diagnostic_dir = file.path(td, "boot"),
  output_dir = file.path(td, "reselection"),
  checkpoint_every = 1,
  resume = FALSE,
  tolerance = 1e-8,
  keep_inner_vectors = TRUE
)
cross <- run_same_vs_cross_selection_VC_diagnostic(
  boot1000 = boot,
  rsel1000 = rsel,
  outer_indices = 1:2,
  B_inner = 20,
  prior_diagnostic_dir = file.path(td, "boot"),
  output_dir = file.path(td, "same_cross"),
  checkpoint_every = 1,
  resume = FALSE,
  tolerance = 1e-8,
  n_cross_shifts = 2L,
  keep_inner_vectors = TRUE
)
layer <- run_two_layer_alignment_decomposition(
  boot1000 = boot,
  cross_reference = cross,
  outer_indices = 1:2,
  B_inner = 20,
  prior_diagnostic_dir = file.path(td, "boot"),
  output_dir = file.path(td, "staged"),
  checkpoint_every = 1,
  resume = FALSE,
  tolerance = 1e-8,
  n_cross_shifts = 2L
)
ksens <- run_selection_intensity_sensitivity(
  boot1000 = boot,
  layer_reference = layer,
  k_values = c(2L, 5L, 10L),
  outer_indices = 1:2,
  B_inner = 20,
  prior_diagnostic_dir = file.path(td, "boot"),
  output_dir = file.path(td, "selection_intensity"),
  checkpoint_every = 1,
  resume = FALSE,
  tolerance = 1e-8,
  n_cross_shifts = 2L,
  MC_R = 0L,
  MC_seed = 20260824L
)
stopifnot(
  nrow(cross$key_table) == 1L,
  nrow(layer$key_table) == 1L,
  all(c(2L, 5L, 10L) %in% ksens$key_table$k)
)
cat("Mechanism-diagnostic smoke test: PASS\n")

cat("\n[5/8] Procedure-level reselection-bootstrap inner smoke test\n")
self_test_reselection_bootstrap_master_only(B = 20)

cat("\n[6/8] RL-UP RNG-fixed prior-sensitivity self-test\n")
source(file.path(.root, "02_GSE_RLUP_prior_sensitivity_RNG_FIXED.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")
run_RLUP_release_self_test()

cat("\n[7/8] RAM-N matched-replicate synthetic self-test\n")
source(file.path(.root, "03_GSE_RAMN_matched_replicates.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")
self_test_RAMN_matched_replicate_comparison()

cat("\n[8/8] Pedigree design + figure plumbing\n")
source(file.path(.vdir, "validate_04_pedigree_structure.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")
run_GSE_pedigree_design_validation()
source(file.path(.vdir, "validate_05_figures.R"), echo = FALSE, chdir = FALSE, encoding = "UTF-8")
run_GSE_figure_smoke_test()

cat("\n============================================================\n")
cat("ALL SMALL GSE PUBLIC-RELEASE SELF-TESTS PASSED\n")
cat("============================================================\n")
