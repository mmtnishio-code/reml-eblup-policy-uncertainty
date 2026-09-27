# Validation report for this assembled public release

## Checks executed in the current build environment

The attached source files were inspected directly before assembly. The release was then checked by a lexical/static audit that ignores comments and quoted strings while balancing `()`, `[]`, and `{}`.

Executed checks:

- all release `.R` files passed the delimiter/string-balance audit;
- no duplicate top-level function definitions remain in the six public analysis R files;
- the only confirmed overwritten duplicate in the integrated public core was the earlier `bayesian_animal_validation()` definition; the later production definition was retained;
- `01_GSE_main_and_diagnostics.R` contains the required production functions for pedigree/A, REML, EBLUP/policy, RAM-N, RL-UP internal implementation, bootstrap, True-VC diagnostics, same/cross, staged decomposition, and `k=2,5,10` sensitivity;
- `02_GSE_RLUP_prior_sensitivity_RNG_FIXED.R` was inspected to confirm the order `set.seed -> seed_latent -> latent_list -> seed_methods` inside `run_RLUP_prior_sensitivity_final()`;
- an RNG-order version tag was added to RL-UP checkpoint signatures to prevent accidental reuse of old-order checkpoints;
- personal `C:/Users/...` and OneDrive paths from the RAM-N example and pedigree-sensitivity runner were removed; production paths are relative or function arguments/options;
- the pedigree-sensitivity release uses the actual attached output fields (`R_VC_fixed`, `R_VC_same_plugin`, `R_VC_cross_plugin`, `baseline_fraction_of_plugin`, `common_alignment_fraction_of_plugin`, `policy_switch_fraction_of_plugin`, `exact_set_match_rate`, `mean_overlap`, `mean_number_replaced`) rather than guessed aliases;
- production-output validators were added for the archived RL-UP prior-sensitivity, True-VC reference, and weak/baseline/strong pedigree-sensitivity summaries;
- figure code reads CSV outputs and does not hard-code manuscript plotting values;
- the final procedure-level reselection bootstrap is released separately as `06_GSE_procedure_level_reselection_bootstrap.R`, preserving the validated 2026-09-25 production logic;
- `06` reconstructs the immutable outer Simulation-II samples, reproduces the historical fixed-policy bootstrap at full S/B before interpreting the new result, and retains per-outer diagnostics used for Supplementary S10-S11;
- public-release post-processing helpers for S10-S11 were added without changing REML, EBLUP, bootstrap generation, selection, or interval calculations;
- `validation/validate_06_reselection_bootstrap.R` was added to check the final S=1000/B=500 add-on object and archived manuscript summaries.

## R runtime limitation in this build environment

`R`/`Rscript` is not installed in the current artifact-build container. Therefore **R's own `parse()` and numerical self-tests could not be executed here**. They are not claimed as having passed.

The release includes the exact commands needed to complete those checks in an R environment:

```bash
Rscript validation/validation_parse_all.R
Rscript validation/run_small_self_tests.R
Rscript validation/validate_06_reselection_bootstrap.R smoke
```

After the S=1000 outputs are available, also run the output validators described in `README.md`. A mismatch must be investigated; do not modify code to force the archived rounded values.

## Why this limitation does not justify changing the scientific code

The user's priority was agreement with the validated manuscript-generating source. The release therefore makes only narrow, auditable changes: public path/source plumbing, the explicitly requested RL-UP RNG-order correction, removal of one function definition that was demonstrably overwritten later during `source()`, and validation/figure wrappers. No attempt was made to reimplement or "clean up" numerical algorithms in another language simply to manufacture a local test run.
