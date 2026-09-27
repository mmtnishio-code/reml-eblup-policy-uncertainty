# Source provenance

This public release is assembled from the validated manuscript-generating R scripts rather than reimplemented from the manuscript text.

- `01_GSE_main_and_diagnostics.R`: validated cumulative core for the main simulations and mechanism diagnostics.
- `02_GSE_RLUP_prior_sensitivity_RNG_FIXED.R`: final RL-UP prior-sensitivity implementation with the validated RNG order.
- `03_GSE_RAMN_matched_replicates.R`: validated matched-replicate post-processing for RAM-N comparisons.
- `04_GSE_pedigree_structure_sensitivity.R`: validated pedigree-structure sensitivity analysis. The 2026-09-27 public copy contains only a path-resolution hotfix for nested `source()` calls under Windows/RStudio; the scientific calculations are unchanged.
- `05_GSE_make_figures.R`: public figure-generation code based on analysis CSV outputs.
- `06_GSE_procedure_level_reselection_bootstrap.R`: final master-RDS-only, diagnostic-preserving procedure-level reselection bootstrap. It reconstructs the original outer samples, reproduces the archived fixed-policy bootstrap before interpreting reselection results, samples full breeding-value vectors conditional on the same inner phenotype draws, and retains the per-outer diagnostics used in Supplementary Tables S10–S11.

The validation folder contains static checks, small self-tests, archived-value checks for analyses 01–06, and the combined production-output validator.
