# GSE public reproducibility code

Public/reproducibility code for the manuscript on finite-sample prediction uncertainty of a **realised REML-EBLUP selection policy**.

The inferential target is

```text
G = c(Y)^T u
```

where `c(Y)` is reconstructed from REML-EBLUP in every **outer** Monte Carlo replicate. In a sire-only biological interpretation the expected next-generation response may be `G/2`, but `G/2` is **not** the simulation estimand in this repository. Some historical internal variable names contain `DeltaG`; these names are retained only where changing them would risk altering validated code, and they refer to `G` in the analyses released here.

## Release principles

The priority of this repository is agreement with the code that generated the manuscript results. The validated numerical algorithms were therefore not refactored merely for style. Public-release changes are limited to path handling, removal of one demonstrably overwritten legacy function definition, the explicitly required RL-UP RNG-order correction, lightweight safety/validation wrappers, figure-generation code, documentation, and the separately validated procedure-level reselection-bootstrap add-on used in the final manuscript.

If an archived production-value validator fails, **do not change the code to force the manuscript value**. First investigate the input RDS/CSV, run settings, RNG order, checkpoint provenance, R version, and BLAS/LAPACK environment.

## Files

### Analysis scripts

- `01_GSE_main_and_diagnostics.R`  
  Cumulative validated core for the paper. Includes pedigree generation, numerator relationship matrix `A`, REML, EBLUP/full PEC, realised top-k policy construction, Conditional, RAM-N, the historically named `Bayesian_policy` / `Bayesian_Beta11` implementation corresponding to manuscript **RL-UP**, fixed-policy parametric bootstrap, True-VC reference, VC-error diagnostics, inner reselection, same/cross diagnostics, staged decomposition, and selection-intensity sensitivity (`k = 2, 5, 10`).

- `02_GSE_RLUP_prior_sensitivity_RNG_FIXED.R`  
  Final RL-UP prior-sensitivity analysis. Compares `h2 ~ Beta(1,1)`, `Beta(2,2)`, and `Beta(0.5,0.5)` with `log(sigma_P^2) ~ N(0,1^2)`. Uses the **final RNG order** described below.

- `03_GSE_RAMN_matched_replicates.R`  
  Post-processes existing Simulation-II production RDS files and compares Conditional / RL-UP / RAM-N on exactly the outer replicates for which RAM-N passed the production diagnostics. No new outer simulation is generated.

- `04_GSE_pedigree_structure_sensitivity.R`  
  Supplementary pedigree-structure sensitivity analysis for weak / baseline / strong structures at `h2=0.20`, `n_pheno=90`, `k=5`. The extraction uses the actual validated object/column names (`R_VC_fixed`, `R_VC_same_plugin`, `R_VC_cross_plugin`, staged-fraction fields, overlap fields), not provisional manuscript aliases.

- `05_GSE_make_figures.R`  
  Regenerates the 9-scenario coverage figure, same/cross figure, and staged-decomposition figure from analysis CSV output. Manuscript plotting values are not manually embedded in the plotting functions.

- `06_GSE_procedure_level_reselection_bootstrap.R`  
  Final procedure-level parametric bootstrap with inner reselection. It reconstructs the original Simulation-II outer samples from the saved bootstrap master object, reproduces the existing fixed-policy bootstrap, samples full breeding-value vectors conditional on the same inner phenotype samples, repeats REML -> EBLUP -> top-k selection within every inner replicate, and evaluates basic-bootstrap coverage/calibration. It also retains the paired fixed-policy branch and per-outer error/boundary diagnostics used for Main Table 9 and Supplementary Tables S10-S11. **This is distinct from the inner-reselection mechanism diagnostic already present in `01_GSE_main_and_diagnostics.R`.**

### Validation scripts

- `validation/validation_parse_all.R` — R `parse()` check for every public `.R` file; required-function check; direct top-level duplicate-function check; personal absolute-path scan.
- `validation/run_small_self_tests.R` — small end-to-end tests, including main Simulation I/II, bootstrap/True-VC, same/cross, staged decomposition, `k=2,5,10`, the formal procedure-level reselection bootstrap, RL-UP prior sensitivity, RAM-N matched processing, pedigree design, and figure plumbing.
- `validation/validate_01_main_production.R` — archived-value check for main Simulation-II coverage/MSE-ratio tables and RAM-N applicability.
- `validation/validation_trueVC_reference.R` — True-VC smoke test and archived-value check.
- `validation/validate_02_RLUP_prior_sensitivity.R` — RL-UP self-test, archived S=1000 check, or full production driver.
- `validation/validate_03_RAMN_matched.R` — synthetic self-test and Supplementary S6 archived-value check.
- `validation/validate_04_pedigree_structure.R` — design-only pedigree check, optional S=2/B=20 pilot, and Supplementary S9 archived-value check.
- `validation/validate_05_figures.R` — synthetic CSV-to-PDF/PNG smoke test.
- `validation/validate_06_reselection_bootstrap.R` — smoke test and archived-value validation for the final procedure-level reselection bootstrap; can also run the exact outer-reconstruction check from the saved bootstrap master object.
- `validation/validate_all_production_outputs.R` — one command to validate any completed production outputs supplied to it.
- `validation/VALIDATION_REPORT.md` — build-time audit record and remaining R-runtime step.

## Software requirements

The scientific analysis code uses **base R only**; no contributed R packages are required.

Use an R 4.x release. The exact minimum R version has not been formally established, so the final repository/Zenodo archive should include the `sessionInfo()` produced by the machine used for the definitive S=1000 runs. Small platform-level differences can arise from numerical optimisation and BLAS/LAPACK implementations.

## Production settings

### Main Simulation II

- pedigree seed: `20260816`
- information-design seed: `20260850`
- outer seed: `20260851`
- `h2 = 0.05, 0.20, 0.40`
- phenotype-information levels: approximately `n = 60, 90, 120`
- candidates: 20 latest-generation males
- main realised policy: top 5 by REML-EBLUP, equal use among selected candidates relative to equal use of all 20 candidates
- outer replicates: `S = 1000`
- RAM-N draws: `M = 2000`
- fixed-policy parametric bootstrap: `B = 500`
- RL-UP final grid: `41 x 41` after adaptive localisation
- main checkpoint interval: 25 outer replicates
- procedure-level reselection bootstrap: `S = 1000`, `B = 500`, same 9 Simulation-II scenarios and the exact original plug-in inner-phenotype RNG stream

### Pedigree-structure sensitivity

Fixed setting: `h2=0.20`, `n_pheno=90`, `k=5`, `S=1000`, `B=500`.

All structures use 24 founders, 3 generations, 40 animals/generation, and pedigree seed `20260816`.

| structure | sire_fraction | dam_fraction | parent_use |
|---|---:|---:|---|
| weak | 1.00 | 1.00 | uniform |
| baseline | 0.75 | 0.95 | uniform |
| strong | 0.50 | 0.80 | weighted |

## First: validate the release itself

From the repository root:

```bash
Rscript validation/validation_parse_all.R
Rscript validation/run_small_self_tests.R
```

The mechanism-diagnostic functions in the validated source require `B_inner >= 20`; therefore the smoke tests use `S=2, B=20` where an inner bootstrap is required. This is a test setting only and is not the manuscript production setting.

## Main analysis: recommended production sequence

Start from a clean R session in the repository directory.

```r
source("01_GSE_main_and_diagnostics.R")

root <- "paper_output"
dir.create(root, recursive = TRUE, showWarnings = FALSE)

# Main Simulation I + II
final <- run_all_final(
  S_I = 1000,
  S_II = 1000,
  M = 2000,
  B = 500,
  ngrid = 41,
  checkpoint_every = 25,
  root_output_dir = file.path(root, "main"),
  resume = TRUE
)

# Standalone True-VC reference: same realised c(Y), no reselection under true VC
oracle1000 <- run_true_VC_oracle_final(
  S = 1000,
  checkpoint_every = 25,
  output_dir = file.path(root, "trueVC"),
  resume = TRUE
)

# Bootstrap generating-model diagnostic and source object for later mechanism diagnostics
boot1000 <- run_bootstrap_generator_diagnostic_final(
  S = 1000,
  B = 500,
  checkpoint_every = 25,
  output_dir = file.path(root, "bootstrap_diagnostic"),
  resume = TRUE
)

# Outer vs inner VC-error decomposition
dec1000 <- run_outer_vs_trueVC_inner_decomposition(
  boot1000 = boot1000,
  prior_diagnostic_dir = file.path(root, "bootstrap_diagnostic"),
  output_dir = file.path(root, "outer_inner_decomposition"),
  checkpoint_every = 10,
  resume = TRUE
)

# Inner reselection mechanism diagnostic
rsel1000 <- run_inner_reselection_VC_diagnostic(
  boot1000 = boot1000,
  dec_inner1000 = dec1000,
  B_inner = 500,
  prior_diagnostic_dir = file.path(root, "bootstrap_diagnostic"),
  output_dir = file.path(root, "reselection"),
  checkpoint_every = 10,
  resume = TRUE
)

# Same-sample vs cross-sample diagnostic
cross1000 <- run_same_vs_cross_selection_VC_diagnostic(
  boot1000 = boot1000,
  rsel1000 = rsel1000,
  B_inner = 500,
  prior_diagnostic_dir = file.path(root, "bootstrap_diagnostic"),
  output_dir = file.path(root, "same_cross"),
  n_cross_shifts = 5,
  checkpoint_every = 10,
  resume = TRUE
)

# Sequential/staged decomposition
layer1000 <- run_two_layer_alignment_decomposition(
  boot1000 = boot1000,
  cross_reference = cross1000,
  B_inner = 500,
  prior_diagnostic_dir = file.path(root, "bootstrap_diagnostic"),
  output_dir = file.path(root, "staged"),
  n_cross_shifts = 5,
  checkpoint_every = 10,
  resume = TRUE
)

# Selection-intensity sensitivity: k = 2, 5, 10
ksens1000 <- run_selection_intensity_sensitivity(
  boot1000 = boot1000,
  layer_reference = layer1000,
  k_values = c(2, 5, 10),
  B_inner = 500,
  prior_diagnostic_dir = file.path(root, "bootstrap_diagnostic"),
  output_dir = file.path(root, "selection_intensity"),
  checkpoint_every = 10,
  resume = TRUE,
  MC_R = 5000
)
```

These calculations are computationally intensive. Resume only from checkpoints generated under the same script version, settings, and RNG scheme.

## Procedure-level reselection bootstrap: final manuscript add-on

The mechanism diagnostic called "inner reselection" inside `01_GSE_main_and_diagnostics.R` asks why the fixed-policy bootstrap misses variance-component-induced EBLUP variation.  The final manuscript additionally evaluates a **procedure-level basic bootstrap interval** in which REML, EBLUP, and top-k policy formation are all repeated inside each bootstrap sample.  These are related but different analyses.

The procedure-level analysis uses the saved master object produced by `run_bootstrap_generator_diagnostic_final()` and deliberately leaves the validated `01` core unchanged.

```r
source("01_GSE_main_and_diagnostics.R", encoding = "UTF-8")
source("06_GSE_procedure_level_reselection_bootstrap.R", encoding = "UTF-8")

# Use the master RDS written by run_bootstrap_generator_diagnostic_final().
boot1000 <- readRDS(
  "paper_output/bootstrap_diagnostic/Simulation_II_bootstrap_generator_ALL.rds"
)

# Cheap new-function test first.
self_test_reselection_bootstrap_master_only(B = 20)

# Strongly recommended before the expensive bootstrap: reconstruct all original
# outer samples and verify Conditional / True-VC summaries against the immutable
# values stored in boot1000.  This fits outer REML models but does no inner B=500 run.
outer_check <- validate_boot1000_outer_reconstruction(
  boot1000 = boot1000,
  tolerance = 1e-8,
  stop_on_failure = TRUE
)
stopifnot(isTRUE(outer_check$all_ok))

# Optional pilot: do not compare pilot values with full S=1000/B=500 manuscript
# values as though they were the same Monte Carlo run.
rboot20 <- run_reselection_bootstrap_from_master(
  boot1000 = boot1000,
  outer_indices = 1:20,
  B_inner = 100,
  output_dir = "paper_output/reselection_bootstrap_pilot",
  checkpoint_every = 5,
  resume = FALSE
)
print_reselection_bootstrap_from_master(rboot20)

# Final production run used for the manuscript.
rboot1000 <- run_reselection_bootstrap_from_master(
  boot1000 = boot1000,
  B_inner = 500,
  output_dir = "paper_output/reselection_bootstrap",
  checkpoint_every = 10,
  resume = TRUE,
  tolerance = 1e-8,
  stop_on_full_validation_failure = TRUE
)
print_reselection_bootstrap_from_master(rboot1000)
stopifnot(all(rboot1000$validation_table$within_tolerance %in% TRUE))

# Recreate the two final supplementary diagnostic tables from saved per-outer
# diagnostics; no REML/EBLUP/bootstrap rerun is required for this step.
rbootm_write_manuscript_diagnostic_tables(
  rboot1000,
  output_dir = "paper_output/reselection_bootstrap"
)
```

The full run writes, among other files:

- `reselection_bootstrap_main_comparison.csv` — fixed vs paired-fixed vs reselection coverage/calibration; Main Table 9;
- `reselection_bootstrap_validation.csv` — exact reproduction checks against immutable stored results;
- `reselection_bootstrap_error_diagnostics_per_outer.csv` — per-outer error, interval, policy, and boundary diagnostics;
- `reselection_bootstrap_master_only_ALL.rds` — master object for the new analysis;
- `Supplementary_Table_S10_reselection_diagnostics.csv` and `Supplementary_Table_S11_boundary_strata.csv` — final post-processed diagnostics.

At full `S=1000, B=500`, do not interpret the new reselection result unless the stored Conditional, True-VC, and fixed-policy bootstrap reproduction checks all pass.  The new code does not overwrite the historical manuscript results stored in `boot1000`.

## RL-UP prior sensitivity: final RNG-fixed run

```r
source("01_GSE_main_and_diagnostics.R")
source("02_GSE_RLUP_prior_sensitivity_RNG_FIXED.R")

# Fast structural/numerical test
run_RLUP_release_self_test()

# Manuscript production settings are frozen inside this wrapper
sens1000 <- run_RLUP_prior_sensitivity_production(
  output_dir = "RLUP_prior_sensitivity_S1000_RNG_FIXED",
  keep_raw = TRUE,
  resume = TRUE
)
```

The wrapper fixes `S=1000`, `M=2000`, `31x31 -> 41x41` adaptive grids, the three primary priors, and the manuscript seeds.

### Why the RNG order was corrected

The final main Simulation-II RNG sequence is:

```r
set.seed(seed_outer)
seed_latent <- sample.int(...)
latent_list <- lapply(seed_latent, ...)
seed_methods <- array(sample.int(...))
```

An earlier prior-sensitivity extension created `seed_methods` before `latent_list`. `simulate_latent_replicate()` changes R's RNG state, so changing this order changes the later method-seed stream even with the same `seed_outer`. The final sensitivity script therefore fixes the order to **latent-list construction before method-seed construction** and stores an RNG-order tag in checkpoint signatures to prevent silent reuse of old-order checkpoints.

### Why the sensitivity baseline can differ slightly from the main RL-UP row

The prior-sensitivity script deliberately constructs a **common fine grid from the union of the high-weight regions across all compared priors** and then uses common random numbers across priors. The main Simulation-II RL-UP calculation instead constructs its adaptive grid for the single baseline prior. Therefore the Beta(1,1) sensitivity row is the baseline **within the prior-sensitivity experiment** and need not be bit-for-bit identical to the standalone main-analysis RL-UP row. The archived values for each analysis are validated separately.

Archived RNG-fixed prior-sensitivity values include, for `h2=0.05`:

- Beta(1,1): coverage `0.905 / 0.908 / 0.930`; MSE ratio `0.691 / 0.809 / 0.898`
- Beta(2,2): coverage `0.829 / 0.855 / 0.885`; MSE ratio `0.483 / 0.568 / 0.644`
- Beta(0.5,0.5): coverage `0.922 / 0.927 / 0.936`; MSE ratio `0.928 / 1.089 / 1.139`

The observed maximum edge mass in the archived run was approximately `0.000978`, well below the diagnostic limit `0.10`.

## RAM-N matched-replicate comparison

This script does **not** simulate new outer data. It reads the scenario-level `*_final.rds` files written by the main Simulation-II production runner and restricts all three methods to exactly the replicates satisfying

```r
!is.null(replicate$RAM_N_policy)
```

Run:

```r
source("03_GSE_RAMN_matched_replicates.R")

matched <- run_RAMN_matched_replicate_comparison(
  input_dir = "paper_output/main/Simulation_II",
  output_dir = "paper_output/RAMN_matched",
  R_boot = 5000
)
```

The RAM-N applicability criteria are implementation diagnostics for the local normal approximation. They are **not universal theoretical thresholds**.

## Pedigree-structure sensitivity

Run the design check first:

```r
GSE_PEDIGREE_RUN_MODE <- "design_only"
source("04_GSE_pedigree_structure_sensitivity.R")
```

Optional reduced pilot:

```r
options(gse.pedigree.pilot.S = 2L, gse.pedigree.pilot.B = 20L)
options(gse.output.root = "paper_output")
GSE_PEDIGREE_RUN_MODE <- "pilot"
source("04_GSE_pedigree_structure_sensitivity.R")
```

Final run:

```r
options(gse.output.root = "paper_output")
GSE_PEDIGREE_RUN_MODE <- "production"
source("04_GSE_pedigree_structure_sensitivity.R")
```

The central mechanism check is whether `R_same >> R_cross ~= R_fixed` persists as pedigree relatedness changes.

## Figures

After production outputs exist:

```r
source("05_GSE_make_figures.R")

make_GSE_coverage_figure(
  simulation_II_csv = "paper_output/main/Simulation_II/Simulation_II_method_performance.csv",
  trueVC_csv = "paper_output/trueVC/Simulation_II_true_VC_oracle_method_performance.csv",
  output_pdf = "paper_output/Figure_coverage_comparison.pdf",
  output_png = "paper_output/Figure_coverage_comparison.png"
)

make_GSE_same_cross_figure(
  same_cross_csv = "paper_output/same_cross/same_vs_cross_selection_VC_key_table.csv",
  output_pdf = "paper_output/Figure_same_cross.pdf",
  output_png = "paper_output/Figure_same_cross.png"
)

make_GSE_staged_average_waterfall_figure(
  staged_csv = "paper_output/staged/two_layer_alignment_key_table.csv",
  output_pdf = "paper_output/Figure_3_staged_average_waterfall.pdf",
  output_png = "paper_output/Figure_3_staged_average_waterfall.png"
)
```

## Scientific meaning of the diagnostic analyses

### True-VC reference

The True-VC reference **retains the realised policy `c(Y)`** formed from ordinary REML-EBLUP. It does not select a new policy using the true variance components. Only the variance components used in prediction/PEC are replaced by the known simulation values. It is therefore a simulation-based mechanism/reference analysis, not a method available with real data.

### Fixed-policy parametric bootstrap

For each outer dataset, the observed `c_obs` is held fixed throughout the inner bootstrap. This targets uncertainty for the already realised policy. It does not automatically reproduce the outer mechanism in which policy formation and the variance-component-induced EBLUP perturbation are generated by the same dataset.

### Same versus cross

- **same**: policy and VC-induced EBV perturbation are constructed from the same inner sample.
- **cross**: the policy is constructed from an independent inner sample while preserving the same marginal policy-generation rule.

`R_same` and `R_cross` are **ratios**: the corresponding inner second moments divided by the outer VC-error second moment. They are not raw second moments.

The diagnostic signature of interest is

```text
R_same >> R_cross ~= R_fixed
```

### Staged decomposition

Values such as `0.120 / 0.812 / 0.068` are **not orthogonal variance components**. They are sequential net changes in the second moment as conditions are replaced step by step: fixed-policy baseline -> common-data dependence -> policy-switch contribution. The final increment can contain cross terms, so the values must not be described as independent variance shares.

### Procedure-level reselection bootstrap

For each inner phenotype sample, the final add-on repeats REML, EBLUP, and top-k selection and evaluates

```text
e*_reselect = c(Y*)^T uhat* - c(Y*)^T u*.
```

The full latent breeding-value vector `u*` is sampled conditional on the same `y*` by Matheron's Gaussian conditional-simulation rule under the outer REML plug-in variance components. The same `u*` is also used in a paired fixed-policy branch, isolating the effect of rebuilding the policy. This procedure-level bootstrap is used to evaluate a basic interval; it does not replace the same/cross or staged analyses, which remain mechanism diagnostics.

### RL-UP

RL-UP is **not a full Bayesian posterior analysis**. It is restricted-likelihood uncertainty propagation using a proper prior weight and low-dimensional weighted integration. The historical internal object name `Bayesian_policy` is retained for computational provenance.

Ahlinder & Waldmann's Bayesian OCS and Fogg et al.'s robust OCS concern **which policy to choose under uncertainty**. This repository concerns whether uncertainty for an **already realised REML-EBLUP policy** is calibrated.

## Output-to-manuscript map

| Analysis output | Manuscript role |
|---|---|
| `Simulation_II_method_performance.csv` | Main Tables 2-3 and 9-condition coverage figure |
| `Simulation_II_REML_RAMN_diagnostics.csv` | Main Table 4 / RAM-N diagnostic description |
| `Simulation_II_true_VC_oracle_method_performance.csv` | Main Table 5 True-VC reference |
| `Simulation_II_bootstrap_generator_method_performance.csv` | plug-in vs True-VC bootstrap generating-model diagnostic |
| `Simulation_II_bootstrap_generator_diagnostics.csv` | Supplementary bootstrap-generator diagnostics |
| `outer_vs_trueVC_inner_component_table.csv` | Main Table 6 mechanism underlying component reproduction |
| `inner_reselection_VC_key_table.csv` | intermediate reselection mechanism diagnostic |
| `same_vs_cross_selection_VC_key_table.csv` | Main Table 7 |
| `two_layer_alignment_key_table.csv` | Main Table 8 and staged mechanism figures |
| `selection_intensity_key_table.csv`, `selection_intensity_compact_table.csv`, `selection_intensity_MC_CI_long.csv` | Supplementary S1-S4 / selection-intensity and staged-detail material |
| `RLUP_prior_sensitivity_method_table.csv` and paired-difference files | Supplementary S5 |
| `RAMN_matched_method_performance.csv`, `RAMN_matched_applicability.csv` | Supplementary S6 |
| True-VC generating-bootstrap outputs | Supplementary S7 |
| main method table + Wilson calculation | Supplementary S8 |
| `pedigree_structure_sensitivity_FINAL_comparison.csv` | Supplementary S9 |
| `reselection_bootstrap_main_comparison.csv` | Main Table 9: fixed-policy vs procedure-level reselection bootstrap |
| `reselection_bootstrap_error_diagnostics_per_outer.csv` | Source diagnostics for Supplementary S10-S11 |
| `Supplementary_Table_S10_reselection_diagnostics.csv` | Supplementary S10 |
| `Supplementary_Table_S11_boundary_strata.csv` | Supplementary S11 |

## Production-output validation

Individual examples:

```bash
Rscript validation/validate_01_main_production.R \
  paper_output/main/Simulation_II/Simulation_II_method_performance.csv

Rscript validation/validation_trueVC_reference.R \
  paper_output/trueVC/Simulation_II_true_VC_oracle_method_performance.csv

Rscript validation/validate_02_RLUP_prior_sensitivity.R validate \
  RLUP_prior_sensitivity_S1000_RNG_FIXED/RLUP_prior_sensitivity_MASTER.rds

Rscript validation/validate_03_RAMN_matched.R validate \
  paper_output/RAMN_matched

Rscript validation/validate_04_pedigree_structure.R validate \
  paper_output/GSE_pedigree_structure_sensitivity_production_S1000_B500/pedigree_structure_sensitivity_FINAL_comparison.csv

Rscript validation/validate_06_reselection_bootstrap.R validate \
  paper_output/reselection_bootstrap/reselection_bootstrap_master_only_ALL.rds
```

Or validate all available production outputs in one call:

```bash
Rscript validation/validate_all_production_outputs.R \
  --main paper_output/main/Simulation_II/Simulation_II_method_performance.csv \
  --truevc paper_output/trueVC/Simulation_II_true_VC_oracle_method_performance.csv \
  --rlup RLUP_prior_sensitivity_S1000_RNG_FIXED/RLUP_prior_sensitivity_MASTER.rds \
  --matched paper_output/RAMN_matched \
  --pedigree paper_output/GSE_pedigree_structure_sensitivity_production_S1000_B500/pedigree_structure_sensitivity_FINAL_comparison.csv \
  --reselection paper_output/reselection_bootstrap/reselection_bootstrap_master_only_ALL.rds
```

## Checkpoints and RNG safety

For exact reproduction:

1. start from a clean R session;
2. use the documented production settings;
3. do not change the order of random-number generation;
4. do not resume from a checkpoint created by a different code version or RNG scheme;
5. preserve the scenario-level final RDS files used for matched-replicate analyses;
6. archive `sessionInfo()` with the definitive run.

The RL-UP sensitivity runner has an explicit RNG-order signature. For other long-running analyses, use separate output directories whenever settings or code versions change.

## GitHub / Zenodo release checklist

Before public deposit:

- run `validation/validation_parse_all.R`;
- run `validation/run_small_self_tests.R`;
- validate completed S=1000 outputs, including the procedure-level reselection-bootstrap master object;
- save the final `sessionInfo()`;
- add the final manuscript citation/DOI when available;
- choose and add an explicit software license if desired;
- create a versioned GitHub release/tag;
- archive that release in Zenodo and place the Zenodo DOI in the manuscript's code-availability statement.
