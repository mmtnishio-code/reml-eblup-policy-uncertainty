# ============================================================================
# GSE PUBLIC RELEASE -- PEDIGREE-STRUCTURE SENSITIVITY (Supplementary S9)
#
# This release preserves the attached, validated analysis sequence and the
# actual output object/column names.  Only execution plumbing was changed:
#   * no workspace clearing;
#   * no setwd() or personal absolute path;
#   * the cumulative public core is sourced from 01_GSE_main_and_diagnostics.R;
#   * output is written under a relative release output directory.
# ============================================================================

# ============================================================
# GSE
# PEDIGREE-STRUCTURE SENSITIVITY ANALYSIS
#
# Primary sensitivity setting:
#   h2 = 0.20
#   information = moderate
#   n_pheno = 90
#   k = 5
#
# Pedigree structures:
#   weak     : broad uniform parent use
#   baseline : original manuscript pedigree
#   strong   : concentrated weighted parent use
#
# IMPORTANT:
#   First run with:
#       RUN_MODE <- "design_only"
#
#   After checking pedigree separation:
#       RUN_MODE <- "pilot"
#
#   Final:
#       RUN_MODE <- "production"
# ============================================================


# Public release: do not clear the caller's workspace.

options(
  warn = 1,
  stringsAsFactors = FALSE
)


# ============================================================
# 0. SETTINGS
# ============================================================

if (!exists("GSE_PEDIGREE_RUN_MODE", inherits = FALSE)) {
  GSE_PEDIGREE_RUN_MODE <- "design_only"
}
RUN_MODE <- GSE_PEDIGREE_RUN_MODE
# Allowed: "design_only", "pilot", "production"


.gse_ped_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.gse_ped_code_dir <- if (!is.null(.gse_ped_script)) {
  dirname(normalizePath(.gse_ped_script, winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

# 01 is the validated cumulative public core.  When this script is sourced
# indirectly from validation/run_small_self_tests.R, sys.frame(1)$ofile can
# refer to the validation wrapper rather than to this file.  Therefore: (i) do
# not re-source 01 if the required core functions are already loaded, and
# (ii) if they are not loaded, search this directory, its parent, and getwd().
.gse_ped_required_core <- c(
  "build_simII_design",
  "relationship_summary",
  "run_bootstrap_generator_diagnostic_final",
  "run_inner_reselection_VC_diagnostic",
  "run_same_vs_cross_selection_VC_diagnostic",
  "run_two_layer_alignment_decomposition"
)

.gse_ped_core_loaded <- all(vapply(
  .gse_ped_required_core,
  exists,
  logical(1),
  mode = "function",
  inherits = TRUE
))

if (!.gse_ped_core_loaded) {
  .gse_ped_core_candidates <- unique(c(
    file.path(.gse_ped_code_dir, "01_GSE_main_and_diagnostics.R"),
    file.path(.gse_ped_code_dir, "..", "01_GSE_main_and_diagnostics.R"),
    file.path(getwd(), "01_GSE_main_and_diagnostics.R")
  ))
  .gse_ped_core_candidates <- normalizePath(
    .gse_ped_core_candidates,
    winslash = "/",
    mustWork = FALSE
  )
  .gse_ped_hits <- .gse_ped_core_candidates[file.exists(.gse_ped_core_candidates)]

  if (length(.gse_ped_hits) == 0L) {
    stop(
      "Cannot locate 01_GSE_main_and_diagnostics.R. ",
      "Place it in the release root (beside 04_GSE_pedigree_structure_sensitivity.R) ",
      "or source 01_GSE_main_and_diagnostics.R before this script."
    )
  }

  .gse_ped_core <- .gse_ped_hits[[1L]]
  source(
    .gse_ped_core,
    echo = FALSE,
    chdir = FALSE,
    encoding = "UTF-8"
  )
}


stopifnot(
  exists("build_simII_design", mode = "function"),
  exists("relationship_summary", mode = "function"),
  exists("run_bootstrap_generator_diagnostic_final", mode = "function"),
  exists("run_inner_reselection_VC_diagnostic", mode = "function"),
  exists("run_same_vs_cross_selection_VC_diagnostic", mode = "function"),
  exists("run_two_layer_alignment_decomposition", mode = "function")
)


# Keep the validated original builder.
build_simII_design_original <- build_simII_design


# ============================================================
# 1. DEFINE THREE PEDIGREE STRUCTURES
# ============================================================

pedigree_spec <- data.frame(

  structure = c(
    "weak",
    "baseline",
    "strong"
  ),

  sire_fraction = c(
    1.00,
    0.75,
    0.50
  ),

  dam_fraction = c(
    1.00,
    0.95,
    0.80
  ),

  parent_use = c(
    "uniform",
    "uniform",
    "weighted"
  ),

  stringsAsFactors = FALSE
)


# Same number of founders, generations, animals/generation,
# candidate number, and pedigree seed.
#
# Thus only the concentration of parent use is changed.

make_structure_design <- function(structure_name) {

  z <- pedigree_spec[
    pedigree_spec$structure == structure_name,
    ,
    drop = FALSE
  ]

  stopifnot(nrow(z) == 1L)

  build_simII_design_original(
    n_founders = 24,
    n_gen = 3,
    n_per_gen = 40,

    sire_fraction =
      z$sire_fraction,

    dam_fraction =
      z$dam_fraction,

    parent_use =
      z$parent_use,

    n_select = 5,

    seed_pedigree =
      20260816
  )
}


designs <- lapply(
  pedigree_spec$structure,
  make_structure_design
)

names(designs) <-
  pedigree_spec$structure


# ============================================================
# 2. PEDIGREE-STRUCTURE DIAGNOSTICS
# ============================================================

safe_q90 <- function(x) {

  if (length(x) == 0L) {
    return(NA_real_)
  }

  unname(
    quantile(
      x,
      0.90,
      names = FALSE
    )
  )
}


pedigree_structure_summary <- function(
    design,
    structure_name) {

  A <- design$A
  ped <- design$ped
  candidates <- design$candidates
  pheno_ids <- design$pheno_ids

  # ----------------------------------------------------------
  # Candidate-candidate relationships
  # ----------------------------------------------------------

  rel_cand <-
    relationship_summary(
      A,
      candidates
    )


  # ----------------------------------------------------------
  # Candidate to other phenotyped animals
  # ----------------------------------------------------------

  other_records <-
    setdiff(
      pheno_ids,
      candidates
    )

  cross_rel <-
    as.numeric(
      A[
        candidates,
        other_records,
        drop = FALSE
      ]
    )


  # ----------------------------------------------------------
  # Parent concentration in latest generation
  # ----------------------------------------------------------

  gmax <-
    max(ped$generation)

  last_gen <-
    ped[
      ped$generation == gmax,
      ,
      drop = FALSE
    ]

  sire_tab <-
    table(last_gen$sire)

  dam_tab <-
    table(last_gen$dam)


  data.frame(

    structure =
      structure_name,

    n_animals =
      nrow(ped),

    n_pheno_full =
      length(pheno_ids),

    n_candidates =
      length(candidates),

    candidate_mean_relationship =
      unname(
        rel_cand["mean_offdiag"]
      ),

    candidate_median_relationship =
      unname(
        rel_cand["median_offdiag"]
      ),

    candidate_q90_relationship =
      unname(
        rel_cand["q90_offdiag"]
      ),

    candidate_max_relationship =
      unname(
        rel_cand["max_offdiag"]
      ),

    candidate_mean_inbreeding =
      unname(
        rel_cand["mean_inbreeding"]
      ),

    candidate_max_inbreeding =
      unname(
        rel_cand["max_inbreeding"]
      ),

    mean_candidate_to_record_relationship =
      mean(cross_rel),

    q90_candidate_to_record_relationship =
      safe_q90(cross_rel),

    n_unique_sires_last_generation =
      length(sire_tab),

    n_unique_dams_last_generation =
      length(dam_tab),

    max_sire_offspring_share =
      max(sire_tab) /
      nrow(last_gen),

    max_dam_offspring_share =
      max(dam_tab) /
      nrow(last_gen),

    row.names = NULL
  )
}


relationship_table <-
  do.call(
    rbind,
    lapply(
      names(designs),
      function(nm) {
        pedigree_structure_summary(
          designs[[nm]],
          nm
        )
      }
    )
  )


cat(
  "\n============================================================\n",
  "PEDIGREE STRUCTURE DIAGNOSTIC\n",
  "============================================================\n",
  sep = ""
)

print(
  relationship_table,
  row.names = FALSE
)


# ------------------------------------------------------------
# Confirm that moderate information remains n = 90
# ------------------------------------------------------------

info_check <- do.call(
  rbind,
  lapply(
    names(designs),
    function(nm) {

      ib <-
        build_information_level_designs(
          base_design =
            designs[[nm]],

          info_fractions =
            c(
              moderate = 0.70
            ),

          seed_info =
            20260850
        )

      x <- ib$info_levels[[1]]

      data.frame(
        structure = nm,
        n_pheno = x$n_pheno,
        candidate_records =
          x$n_candidate_records,
        noncandidate_records =
          x$n_noncandidate_records,
        row.names = NULL
      )
    }
  )
)


cat(
  "\n============================================================\n",
  "INFORMATION DESIGN CHECK\n",
  "============================================================\n",
  sep = ""
)

print(
  info_check,
  row.names = FALSE
)


stopifnot(
  all(
    info_check$n_pheno == 90L
  ),
  all(
    info_check$candidate_records == 20L
  )
)


# ============================================================
# 3. STOP HERE FOR THE FIRST RUN
# ============================================================

if (RUN_MODE == "design_only") {

  cat(
    "\n============================================================\n",
    "DESIGN-ONLY CHECK COMPLETED\n",
    "============================================================\n",
    "No Monte Carlo analysis has been run.\n",
    "Check that candidate relationships increase clearly:\n",
    "    weak < baseline < strong\n",
    "before running the pilot or production analysis.\n",
    "============================================================\n",
    sep = ""
  )

} else {


# ============================================================
# 4. PILOT / PRODUCTION SETTINGS
# ============================================================

  if (RUN_MODE == "pilot") {

    S <- as.integer(getOption("gse.pedigree.pilot.S", 50L))
    B <- as.integer(getOption("gse.pedigree.pilot.B", 50L))
    if (S < 2L) stop("Pilot S must be >= 2.")
    if (B < 20L) stop("Pilot B must be >= 20 because the validated inner diagnostics require B_inner >= 20.")

    run_label <-
      "pilot_S50_B50"

  } else if (RUN_MODE == "production") {

    S <- 1000L
    B <- 500L

    run_label <-
      "production_S1000_B500"

  } else {

    stop(
      "RUN_MODE must be design_only, pilot, or production."
    )
  }


  h2_value <- 0.20

  info_fractions <-
    c(
      moderate = 0.70
    )

  n_select <- 5L


  output_base <- getOption(
    "gse.output.root",
    file.path(getwd(), "GSE_output")
  )

  output_root <- file.path(
    output_base,
    paste0(
      "GSE_pedigree_structure_sensitivity_",
      run_label
    )
  )


  dir.create(
    output_root,
    recursive = TRUE,
    showWarnings = FALSE
  )


# ============================================================
# 5. VALIDATION HELPER
# ============================================================

  all_validation_true <- function(v) {

    if (
      is.null(v) ||
      !"within_tolerance" %in%
        names(v)
    ) {
      return(FALSE)
    }

    w <-
      v$within_tolerance[
        !is.na(
          v$within_tolerance
        )
      ]

    length(w) > 0L &&
      all(w)
  }


# ============================================================
# 6. RUN ONE PEDIGREE STRUCTURE
# ============================================================

  run_one_pedigree_structure <- function(
      structure_name,
      design) {

    cat(
      "\n\n############################################################\n",
      "PEDIGREE STRUCTURE: ",
      toupper(structure_name),
      "\n############################################################\n",
      sep = ""
    )


    structure_dir <-
      file.path(
        output_root,
        structure_name
      )

    boot_dir <-
      file.path(
        structure_dir,
        "03_bootstrap"
      )

    rsel_dir <-
      file.path(
        structure_dir,
        "05_reselection"
      )

    cross_dir <-
      file.path(
        structure_dir,
        "06_same_cross"
      )

    layer_dir <-
      file.path(
        structure_dir,
        "07_staged"
      )


    dir.create(
      structure_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )


    # --------------------------------------------------------
    # IMPORTANT:
    # Temporarily replace build_simII_design().
    #
    # Existing validated analysis functions therefore use
    # exactly this pedigree, while all inferential machinery
    # remains unchanged.
    # --------------------------------------------------------

    original_builder <-
      get(
        "build_simII_design",
        envir = .GlobalEnv
      )


    assign(
      "build_simII_design",
      function(...) {
        design
      },
      envir = .GlobalEnv
    )


    on.exit(
      assign(
        "build_simII_design",
        original_builder,
        envir = .GlobalEnv
      ),
      add = TRUE
    )


    # --------------------------------------------------------
    # Check the installed design
    # --------------------------------------------------------

    test_design <-
      build_simII_design()

    stopifnot(
      identical(
        test_design$ped,
        design$ped
      ),
      max(
        abs(
          test_design$A -
            design$A
        )
      ) < 1e-12
    )


    cat(
      "\nCandidate relationship summary:\n"
    )

    print(
      relationship_summary(
        design$A,
        design$candidates
      )
    )


# ============================================================
# STEP A
# Bootstrap-generator diagnostic
#
# Provides:
#   Conditional plug-in
#   True-VC reference
#   outer data required by subsequent diagnostics
# ============================================================

    cat(
      "\n------------------------------------------------------------\n",
      "STEP A: BOOTSTRAP / TRUE-VC REFERENCE\n",
      "------------------------------------------------------------\n",
      sep = ""
    )


    boot <-
      run_bootstrap_generator_diagnostic_final(

        S = S,

        h2_values =
          h2_value,

        info_fractions =
          info_fractions,

        n_select =
          n_select,

        sigma2_P = 1,

        beta = 0,

        target = 0,

        B = B,

        seed_info =
          20260850,

        seed_outer =
          20260851,

        checkpoint_every =
          25,

        output_dir =
          boot_dir,

        resume =
          TRUE,

        keep_raw_in_master =
          FALSE,

        return_inner =
          FALSE
      )


    stopifnot(
      nrow(
        boot$method_table
      ) == 4L
    )


# ============================================================
# STEP B
# Inner reselection
# ============================================================

    cat(
      "\n------------------------------------------------------------\n",
      "STEP B: INNER RESELECTION\n",
      "------------------------------------------------------------\n",
      sep = ""
    )


    rsel <-
      run_inner_reselection_VC_diagnostic(

        boot1000 =
          boot,

        dec_inner1000 =
          NULL,

        outer_indices =
          seq_len(S),

        B_inner =
          B,

        prior_diagnostic_dir =
          boot_dir,

        output_dir =
          rsel_dir,

        checkpoint_every =
          25,

        resume =
          TRUE,

        tolerance =
          1e-9,

        keep_inner_vectors =
          TRUE,

        keep_per_outer_in_master =
          FALSE
      )


    stopifnot(
      all_validation_true(
        rsel$validation
      )
    )


# ============================================================
# STEP C
# SAME vs CROSS
# ============================================================

    cat(
      "\n------------------------------------------------------------\n",
      "STEP C: SAME vs CROSS\n",
      "------------------------------------------------------------\n",
      sep = ""
    )


    cross <-
      run_same_vs_cross_selection_VC_diagnostic(

        boot1000 =
          boot,

        rsel1000 =
          rsel,

        outer_indices =
          seq_len(S),

        B_inner =
          B,

        prior_diagnostic_dir =
          boot_dir,

        output_dir =
          cross_dir,

        checkpoint_every =
          25,

        resume =
          TRUE,

        tolerance =
          1e-9,

        n_cross_shifts =
          5L,

        keep_inner_vectors =
          TRUE,

        keep_per_outer_in_master =
          FALSE
      )


    stopifnot(
      all_validation_true(
        cross$validation
      )
    )


# ============================================================
# STEP D
# Staged / two-layer decomposition
# ============================================================

    cat(
      "\n------------------------------------------------------------\n",
      "STEP D: STAGED DECOMPOSITION\n",
      "------------------------------------------------------------\n",
      sep = ""
    )


    layer <-
      run_two_layer_alignment_decomposition(

        boot1000 =
          boot,

        cross_reference =
          cross,

        outer_indices =
          seq_len(S),

        B_inner =
          B,

        prior_diagnostic_dir =
          boot_dir,

        output_dir =
          layer_dir,

        checkpoint_every =
          25,

        resume =
          TRUE,

        tolerance =
          1e-9,

        n_cross_shifts =
          5L,

        keep_per_outer_in_master =
          FALSE
      )


    stopifnot(
      all_validation_true(
        layer$validation
      )
    )


# ============================================================
# 7. EXTRACT CORE RESULTS
# ============================================================

    mt <-
      boot$method_table


    cond <-
      mt[
        mt$method ==
          "conditional_plugin",
        ,
        drop = FALSE
      ]


    trueVC <-
      mt[
        mt$method ==
          "conditional_trueVC_oracle",
        ,
        drop = FALSE
      ]


    stopifnot(
      nrow(cond) == 1L,
      nrow(trueVC) == 1L
    )


    cr <-
      cross$key_table

    ly <-
      layer$key_table


    stopifnot(
      nrow(cr) == 1L,
      nrow(ly) == 1L
    )


    # Exact staged identity
    fraction_sum <-
      ly$baseline_fraction_of_plugin +
      ly$common_alignment_fraction_of_plugin +
      ly$policy_switch_fraction_of_plugin


    stopifnot(
      abs(
        fraction_sum - 1
      ) < 1e-9
    )


    P_alignment <-
      1 -
      cr$R_VC_cross_plugin /
      cr$R_VC_same_plugin


    rel <-
      pedigree_structure_summary(
        design,
        structure_name
      )


    result_row <-
      data.frame(

        structure =
          structure_name,

        h2 =
          h2_value,

        n_pheno =
          cond$n_pheno,

        k =
          n_select,

        # --------------------------------------
        # Pedigree structure
        # --------------------------------------

        candidate_mean_relationship =
          rel$candidate_mean_relationship,

        candidate_q90_relationship =
          rel$candidate_q90_relationship,

        candidate_max_relationship =
          rel$candidate_max_relationship,

        candidate_mean_inbreeding =
          rel$candidate_mean_inbreeding,

        mean_candidate_to_record_relationship =
          rel$mean_candidate_to_record_relationship,


        # --------------------------------------
        # Conditional vs True-VC
        # --------------------------------------

        conditional_coverage =
          cond$inclusion,

        conditional_MSE_ratio =
          cond$MSE_ratio,

        trueVC_coverage =
          trueVC$inclusion,

        trueVC_MSE_ratio =
          trueVC$MSE_ratio,


        # --------------------------------------
        # Same / cross mechanism
        # --------------------------------------

        R_fixed =
          cr$R_VC_fixed,

        R_same =
          cr$R_VC_same_plugin,

        R_cross =
          cr$R_VC_cross_plugin,

        same_over_cross =
          cr$R_VC_same_plugin /
          cr$R_VC_cross_plugin,

        P_alignment =
          P_alignment,


        # --------------------------------------
        # Staged decomposition
        # --------------------------------------

        baseline_fraction =
          ly$baseline_fraction_of_plugin,

        common_alignment_fraction =
          ly$common_alignment_fraction_of_plugin,

        policy_switch_fraction =
          ly$policy_switch_fraction_of_plugin,


        # --------------------------------------
        # Selection stability
        # --------------------------------------

        exact_set_match_rate =
          ly$exact_set_match_rate,

        mean_overlap =
          ly$mean_overlap,

        mean_number_replaced =
          ly$mean_number_replaced,

        stringsAsFactors = FALSE
      )


    write.csv(
      result_row,
      file.path(
        structure_dir,
        paste0(
          "pedigree_structure_",
          structure_name,
          "_core_results.csv"
        )
      ),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )


    saveRDS(
      list(
        structure = structure_name,
        design = design,
        relationship_summary = rel,
        boot = boot,
        rsel = rsel,
        cross = cross,
        layer = layer,
        core_result = result_row
      ),
      file.path(
        structure_dir,
        paste0(
          "pedigree_structure_",
          structure_name,
          "_complete.rds"
        )
      )
    )


    cat(
      "\nCore result:\n"
    )

    print(
      result_row,
      row.names = FALSE
    )


    result_row
  }


# ============================================================
# 8. RUN ALL THREE STRUCTURES
# ============================================================

  results <-
    lapply(
      names(designs),
      function(nm) {

        run_one_pedigree_structure(
          structure_name = nm,
          design = designs[[nm]]
        )
      }
    )


  final_table <-
    do.call(
      rbind,
      results
    )


  rownames(final_table) <-
    NULL


# ============================================================
# 9. FINAL COMPARISON
# ============================================================

  cat(
    "\n\n============================================================\n",
    "PEDIGREE STRUCTURE SENSITIVITY: FINAL COMPARISON\n",
    "============================================================\n",
    sep = ""
  )


  print(
    final_table,
    row.names = FALSE
  )


  write.csv(
    final_table,
    file.path(
      output_root,
      "pedigree_structure_sensitivity_FINAL_comparison.csv"
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )


  write.csv(
    relationship_table,
    file.path(
      output_root,
      "pedigree_structure_relationship_summary.csv"
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )


  saveRDS(
    list(
      settings = list(
        RUN_MODE = RUN_MODE,
        S = S,
        B = B,
        h2 = h2_value,
        information_fraction = 0.70,
        k = n_select,
        seed_pedigree = 20260816,
        seed_info = 20260850,
        seed_outer = 20260851
      ),
      pedigree_spec =
        pedigree_spec,
      relationship_table =
        relationship_table,
      final_table =
        final_table
    ),
    file.path(
      output_root,
      "pedigree_structure_sensitivity_FINAL.rds"
    )
  )


  cat(
    "\n============================================================\n",
    "STATUS: PASSED\n",
    "Output directory:\n",
    normalizePath(
      output_root,
      winslash = "/",
      mustWork = FALSE
    ),
    "\n============================================================\n",
    sep = ""
  )

}

# ============================================================
# 10. ARCHIVED PRODUCTION-RESULT VALIDATION
# ============================================================

pedigree_structure_expected_results <- function() {
  data.frame(
    structure = c("weak", "baseline", "strong"),
    candidate_mean_relationship = c(0.078454, 0.089145, 0.208470),
    candidate_q90_relationship = c(0.140625, 0.171875, 0.375000),
    candidate_max_relationship = c(0.406250, 0.531250, 0.687500),
    candidate_mean_inbreeding = c(0.029688, 0.031250, 0.087500),
    mean_candidate_to_record_relationship = c(0.075570, 0.086523, 0.154672),
    conditional_coverage = c(0.626, 0.646, 0.635),
    conditional_MSE_ratio = c(0.304926, 0.304527, 0.342428),
    trueVC_coverage = c(0.938, 0.942, 0.961),
    trueVC_MSE_ratio = c(1.007175, 0.968949, 1.050055),
    R_fixed = c(0.106999, 0.112134, 0.134213),
    R_same = c(0.961151, 0.946948, 1.029450),
    R_cross = c(0.107236, 0.112826, 0.135382),
    P_alignment = c(0.888430, 0.880853, 0.868491),
    baseline_fraction = c(0.111404, 0.118783, 0.130986),
    common_alignment_fraction = c(0.823538, 0.808662, 0.771699),
    policy_switch_fraction = c(0.065057, 0.072555, 0.097315),
    exact_set_match_rate = c(0.673050, 0.652632, 0.607570),
    mean_overlap = c(0.929577, 0.924080, 0.911843),
    mean_number_replaced = c(0.352114, 0.379600, 0.440786),
    stringsAsFactors = FALSE
  )
}

validate_pedigree_structure_production <- function(
    x,
    rounded_tolerance = 0.0006,
    stop_on_failure = TRUE) {

  obs <- if (is.character(x) && length(x) == 1L) {
    read.csv(x, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.data.frame(x)) {
    x
  } else if (is.list(x) && !is.null(x$final_table)) {
    x$final_table
  } else {
    stop("x must be the final comparison data.frame, final RDS object, or CSV path.")
  }

  exp <- pedigree_structure_expected_results()
  fields <- setdiff(names(exp), "structure")
  miss <- setdiff(c("structure", fields), names(obs))
  if (length(miss) > 0L) stop("Missing pedigree-sensitivity columns: ", paste(miss, collapse = ", "))

  obs2 <- obs[match(exp$structure, obs$structure), c("structure", fields), drop = FALSE]
  if (any(is.na(match(exp$structure, obs$structure)))) stop("One or more expected structures are missing.")

  checks <- data.frame(structure = exp$structure, stringsAsFactors = FALSE)
  all_ok <- TRUE
  for (nm in fields) {
    d <- abs(as.numeric(obs2[[nm]]) - as.numeric(exp[[nm]]))
    checks[[paste0(nm, "_abs_diff")]] <- d
    ok <- is.finite(d) & d <= rounded_tolerance
    checks[[paste0(nm, "_ok")]] <- ok
    all_ok <- all_ok && all(ok)
  }

  # Mechanism-level check central to Supplementary S9.
  mechanism_ok <- all(
    obs2$R_same > 5 * obs2$R_cross &
      abs(obs2$R_cross - obs2$R_fixed) < 0.01
  )
  all_ok <- all_ok && mechanism_ok

  if (isTRUE(stop_on_failure) && !all_ok) {
    print(checks, row.names = FALSE)
    stop("Pedigree-structure production validation failed. Investigate provenance/settings rather than forcing values.")
  }

  cat("Pedigree-structure production validation: ", if (all_ok) "PASS" else "CHECK REQUIRED", "\n", sep = "")
  cat("Mechanism R_same >> R_cross ~= R_fixed: ", if (mechanism_ok) "PASS" else "FAIL", "\n", sep = "")
  invisible(checks)
}
