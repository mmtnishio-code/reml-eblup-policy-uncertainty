# ============================================================================
# Parse/static validation for the GSE public release
# Base R only.
# ============================================================================

.args <- commandArgs(trailingOnly = TRUE)
.this <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.validation_dir <- if (!is.null(.this)) {
  dirname(normalizePath(.this, winslash = "/", mustWork = FALSE))
} else {
  file.path(getwd(), "validation")
}
.root <- normalizePath(file.path(.validation_dir, ".."), winslash = "/", mustWork = FALSE)

required_root_files <- c(
  "01_GSE_main_and_diagnostics.R",
  "02_GSE_RLUP_prior_sensitivity_RNG_FIXED.R",
  "03_GSE_RAMN_matched_replicates.R",
  "04_GSE_pedigree_structure_sensitivity.R",
  "05_GSE_make_figures.R",
  "06_GSE_procedure_level_reselection_bootstrap.R",
  "README.md"
)

missing <- required_root_files[!file.exists(file.path(.root, required_root_files))]
if (length(missing) > 0L) {
  stop("Missing required release files: ", paste(missing, collapse = ", "))
}

r_files <- c(
  list.files(.root, pattern = "\\.R$", full.names = TRUE),
  list.files(.validation_dir, pattern = "\\.R$", full.names = TRUE)
)
r_files <- sort(unique(normalizePath(r_files, winslash = "/", mustWork = TRUE)))

cat("Parsing ", length(r_files), " R files...\n", sep = "")
for (f in r_files) {
  parse(file = f, keep.source = FALSE)
  cat("  PASS: ", basename(f), "\n", sep = "")
}

# Check user-specific absolute paths in public source files.  This intentionally
# ignores ordinary relative paths and manuscript prose.
text_files <- c(
  list.files(.root, pattern = "\\.R$", full.names = TRUE),
  list.files(.validation_dir, pattern = "\\.R$", full.names = TRUE)
)
for (f in unique(text_files)) {
  if (identical(normalizePath(f, winslash = "/", mustWork = FALSE),
                normalizePath(file.path(.validation_dir, "validation_parse_all.R"), winslash = "/", mustWork = FALSE))) {
    next
  }
  x <- readLines(f, warn = FALSE, encoding = "UTF-8")
  bad <- grepl(
    "([A-Za-z]:[/\\\\]Users[/\\\\]|/Users/|/home/[^/]+/|OneDrive)",
    x,
    perl = TRUE
  )
  if (any(bad)) {
    stop(
      "User-specific absolute path found in ", basename(f),
      " at line(s): ", paste(which(bad), collapse = ", ")
    )
  }
}

# Detect duplicate direct top-level function assignments in each R file.
# Historical helper names inside different function bodies are not counted.
get_top_function_names <- function(path) {
  ex <- parse(file = path, keep.source = FALSE)
  out <- character()
  for (e in ex) {
    if (
      is.call(e) && length(e) >= 3L &&
      identical(e[[1]], as.name("<-")) &&
      is.symbol(e[[2]]) &&
      is.call(e[[3]]) &&
      identical(e[[3]][[1]], as.name("function"))
    ) {
      out <- c(out, as.character(e[[2]]))
    }
  }
  out
}

for (f in r_files) {
  nm <- get_top_function_names(f)
  dup <- unique(nm[duplicated(nm)])
  if (length(dup) > 0L) {
    stop(
      "Duplicate top-level function definition(s) in ", basename(f), ": ",
      paste(dup, collapse = ", ")
    )
  }
}

# Source the main core and verify the functions required by the paper.
source(file.path(.root, "01_GSE_main_and_diagnostics.R"),
       echo = FALSE, chdir = FALSE, encoding = "UTF-8")
source(file.path(.root, "06_GSE_procedure_level_reselection_bootstrap.R"),
       echo = FALSE, chdir = FALSE, encoding = "UTF-8")

required_functions <- c(
  "make_moderate_pedigree",
  "make_A",
  "build_simII_design",
  "fit_reml_animal",
  "predict_ebv_validation",
  "make_policy",
  "policy_stats_animal",
  "run_simulation_II_final",
  "run_all_final",
  "run_true_VC_oracle_final",
  "run_bootstrap_generator_diagnostic_final",
  "run_outer_vs_trueVC_inner_decomposition",
  "run_inner_reselection_VC_diagnostic",
  "run_same_vs_cross_selection_VC_diagnostic",
  "run_two_layer_alignment_decomposition",
  "run_selection_intensity_sensitivity",
  "validate_boot1000_outer_reconstruction",
  "parametric_bootstrap_reselection_from_master_one",
  "run_reselection_bootstrap_from_master",
  "self_test_reselection_bootstrap_master_only",
  "rbootm_make_supplementary_S10",
  "rbootm_make_supplementary_S11"
)

missing_fun <- required_functions[
  !vapply(required_functions, exists, logical(1), mode = "function")
]
if (length(missing_fun) > 0L) {
  stop("Missing required functions after sourcing 01: ", paste(missing_fun, collapse = ", "))
}

cat("\nAll R files parsed successfully.\n")
cat("No user-specific absolute paths were detected.\n")
cat("No duplicate direct top-level function definitions were detected.\n")
cat("All required main-analysis and procedure-level reselection functions exist.\n")
