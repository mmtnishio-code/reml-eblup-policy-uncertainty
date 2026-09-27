# Release notes — 2026-09-27 release candidate

This release candidate supersedes the earlier five-script public package.

Changes in this candidate:
- adds `06_GSE_procedure_level_reselection_bootstrap.R`, the final procedure-level parametric bootstrap with inner REML–EBLUP reselection used for Main Table 9 and Supplementary Tables S10–S11;
- adds `validation/validate_06_reselection_bootstrap.R` and integrates script 06 into the public self-tests and production-output validation workflow;
- updates `README.md` and the public execution guide;
- applies the Windows/RStudio nested-`source()` path hotfix to `04_GSE_pedigree_structure_sensitivity.R` without changing its scientific calculations;
- retains the validated scientific logic in scripts 01–05.

Before creating GitHub/Zenodo v1.0.0, run the full validation sequence described in README.md and save the definitive `sessionInfo()`.
