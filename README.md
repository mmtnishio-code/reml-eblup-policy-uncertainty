# REML-EBLUP policy uncertainty

Reproducibility code for a simulation study of finite-sample prediction uncertainty for a realised selection policy formed from REML-EBLUP.

## Scientific target

The main inferential target is

```text
G = c(Y)^T u
```

where `u` is the vector of true breeding values and `c(Y)` is the realised policy contrast constructed from REML-EBLUP using the observed data `Y`.

The repository evaluates uncertainty for the realised policy itself. In a sire-only interpretation, the expected next-generation response may be `G/2`, but `G/2` is not the simulation estimand used here.

All analyses are simulation-based; no empirical animal-level or personal data are included.

## Repository contents

- `01_GSE_main_and_diagnostics.R` — main simulations and mechanism diagnostics, including REML, EBLUP, Conditional PEV, RAM-N, RL-UP, fixed-policy bootstrap, True-VC reference, same/cross diagnostics, staged decomposition, and selection-intensity sensitivity.
- `02_GSE_RLUP_prior_sensitivity_RNG_FIXED.R` — RL-UP prior-sensitivity analysis.
- `03_GSE_RAMN_matched_replicates.R` — matched-replicate comparison of Conditional, RL-UP, and RAM-N.
- `04_GSE_pedigree_structure_sensitivity.R` — pedigree-structure sensitivity analysis.
- `05_GSE_make_figures.R` — figure-generation code from analysis CSV outputs.
- `06_GSE_procedure_level_reselection_bootstrap.R` — procedure-level bootstrap that repeats REML, EBLUP, and top-k selection within each bootstrap sample.
- `validation/` — parse checks, small end-to-end tests, analysis-specific archived-value checks, and production-output validation.

## Software requirements

The scientific analysis code uses base R only; no contributed R packages are required.

Use an R 4.x release. Small numerical differences can occur across R, BLAS/LAPACK, and operating-system versions, so the final archived release should retain the exact `sessionInfo()` from the definitive production environment.

## Main simulation settings

### Simulation II

- pedigree seed: `20260816`
- information-design seed: `20260850`
- outer seed: `20260851`
- heritability: `h2 = 0.05, 0.20, 0.40`
- phenotype-information levels: approximately `n = 60, 90, 120`
- candidates: 20 latest-generation males
- realised policy: top 5 candidates by REML-EBLUP, with equal use among selected candidates
- reference policy: equal use of all 20 candidates
- outer replicates: `S = 1000`
- RAM-N draws: `M = 2000`
- fixed-policy parametric bootstrap: `B = 500`
- RL-UP adaptive grid: `31 x 31` coarse grid followed by a `41 x 41` fine grid
- procedure-level reselection bootstrap: `S = 1000`, `B = 500`

### Pedigree-structure sensitivity

The sensitivity analysis fixes `h2 = 0.20`, approximately `n = 90`, `k = 5`, `S = 1000`, and `B = 500`, and compares weak, baseline, and strong pedigree structures.

## Quick validation

From the repository root:

```bash
Rscript validation/validation_parse_all.R
Rscript validation/run_small_self_tests.R
```

Run these checks before long production analyses.

## Main production workflow

Start from a clean R session in the repository root.

```r
source("01_GSE_main_and_diagnostics.R")

root <- "paper_output"
dir.create(root, recursive = TRUE, showWarnings = FALSE)

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

oracle1000 <- run_true_VC_oracle_final(
  S = 1000,
  checkpoint_every = 25,
  output_dir = file.path(root, "trueVC"),
  resume = TRUE
)

boot1000 <- run_bootstrap_generator_diagnostic_final(
  S = 1000,
  B = 500,
  checkpoint_every = 25,
  output_dir = file.path(root, "bootstrap_diagnostic"),
  resume = TRUE
)

dec1000 <- run_outer_vs_trueVC_inner_decomposition(
  boot1000 = boot1000,
  prior_diagnostic_dir = file.path(root, "bootstrap_diagnostic"),
  output_dir = file.path(root, "outer_inner_decomposition"),
  checkpoint_every = 10,
  resume = TRUE
)

rsel1000 <- run_inner_reselection_VC_diagnostic(
  boot1000 = boot1000,
  dec_inner1000 = dec1000,
  B_inner = 500,
  prior_diagnostic_dir = file.path(root, "bootstrap_diagnostic"),
  output_dir = file.path(root, "reselection"),
  checkpoint_every = 10,
  resume = TRUE
)

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

The production calculations are computationally intensive. Resume only from checkpoints generated with the same code version, settings, and RNG scheme.

## Procedure-level reselection bootstrap

The formal procedure-level bootstrap in `06_GSE_procedure_level_reselection_bootstrap.R` is distinct from the inner-reselection mechanism diagnostic in `01_GSE_main_and_diagnostics.R`.

```r
source("01_GSE_main_and_diagnostics.R", encoding = "UTF-8")
source("06_GSE_procedure_level_reselection_bootstrap.R", encoding = "UTF-8")

boot1000 <- readRDS(
  "paper_output/bootstrap_diagnostic/Simulation_II_bootstrap_generator_ALL.rds"
)

outer_check <- validate_boot1000_outer_reconstruction(
  boot1000 = boot1000,
  tolerance = 1e-8,
  stop_on_failure = TRUE
)
stopifnot(isTRUE(outer_check$all_ok))

rboot1000 <- run_reselection_bootstrap_from_master(
  boot1000 = boot1000,
  B_inner = 500,
  output_dir = "paper_output/reselection_bootstrap",
  checkpoint_every = 10,
  resume = TRUE,
  tolerance = 1e-8,
  keep_inner_vectors = FALSE,
  keep_per_outer_in_master = TRUE,
  stop_on_full_validation_failure = TRUE
)

print_reselection_bootstrap_from_master(rboot1000)
stopifnot(all(rboot1000$validation_table$within_tolerance %in% TRUE))
```

## Additional analyses

RL-UP prior sensitivity:

```r
source("01_GSE_main_and_diagnostics.R")
source("02_GSE_RLUP_prior_sensitivity_RNG_FIXED.R")

sens1000 <- run_RLUP_prior_sensitivity_production(
  output_dir = "RLUP_prior_sensitivity_S1000_RNG_FIXED",
  keep_raw = TRUE,
  resume = TRUE
)
```

RAM-N matched-replicate comparison:

```r
source("03_GSE_RAMN_matched_replicates.R")
matched <- run_RAMN_matched_replicate_comparison(
  input_dir = "paper_output/main/Simulation_II",
  output_dir = "paper_output/RAMN_matched",
  R_boot = 5000
)
```

Pedigree-structure sensitivity:

```r
options(gse.output.root = "paper_output")
GSE_PEDIGREE_RUN_MODE <- "production"
source("04_GSE_pedigree_structure_sensitivity.R")
```

Figures:

```r
source("05_GSE_make_figures.R")
```

The figure functions read analysis CSV outputs rather than hard-coded manuscript values.

## Main output-to-manuscript map

| Output | Manuscript role |
|---|---|
| `Simulation_II_method_performance.csv` | Main method-performance tables |
| `Simulation_II_REML_RAMN_diagnostics.csv` | REML/RAM-N diagnostics |
| `Simulation_II_true_VC_oracle_method_performance.csv` | True-VC reference |
| `same_vs_cross_selection_VC_key_table.csv` | Same/cross diagnostic |
| `two_layer_alignment_key_table.csv` | Staged decomposition |
| `reselection_bootstrap_main_comparison.csv` | Fixed-policy vs procedure-level reselection bootstrap |
| `Supplementary_Table_S10_reselection_diagnostics.csv` | Reselection error diagnostics |
| `Supplementary_Table_S11_boundary_strata.csv` | Boundary-stratified diagnostics |

## Production-output validation

Use the analysis-specific validators in `validation/` or the combined validator after production outputs are available.

```bash
Rscript validation/validate_all_production_outputs.R \
  --main paper_output/main/Simulation_II/Simulation_II_method_performance.csv \
  --truevc paper_output/trueVC/Simulation_II_true_VC_oracle_method_performance.csv \
  --rlup RLUP_prior_sensitivity_S1000_RNG_FIXED/RLUP_prior_sensitivity_MASTER.rds \
  --matched paper_output/RAMN_matched \
  --pedigree paper_output/GSE_pedigree_structure_sensitivity_production_S1000_B500/pedigree_structure_sensitivity_FINAL_comparison.csv \
  --reselection paper_output/reselection_bootstrap/reselection_bootstrap_master_only_ALL.rds
```

If a validator fails, check the input files, run settings, RNG order, checkpoint provenance, R version, and numerical environment before changing scientific code.

## Reproducibility notes

- Start production analyses from a clean R session.
- Do not change the documented RNG order.
- Do not resume from checkpoints created with a different code version or RNG scheme.
- Keep large generated outputs and checkpoints outside the version-controlled source tree unless intentionally archived.
- Archive the definitive `sessionInfo()` with the versioned release.

## Code availability

GitHub repository:

https://github.com/mmtnishio-code/reml-eblup-policy-uncertainty

A versioned archival copy and DOI will be added after the final release is deposited in Zenodo.

## Citation

Please cite the associated manuscript and the archived Zenodo release once their bibliographic information is available.

## License

No software license is specified yet. Add the intended license before the archival release if reuse under an explicit open-source license is desired.
