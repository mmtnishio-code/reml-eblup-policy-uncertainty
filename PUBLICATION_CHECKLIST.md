# GitHub / Zenodo publication checklist

Release target: `2026-09-08-final-r4`

## Already completed

- Release-wide parse/path/duplicate-function validation: PASS.
- Seven-stage small end-to-end self-test suite after the r4 repair: PASS.
- Main Simulation-II S=1000 archived-value validation: PASS.
- True-VC S=1000 archived-value validation: PASS.
- RL-UP RNG-fixed prior-sensitivity S=1000 validation: PASS.
- RAM-N matched-replicate validation: PASS.
- Pedigree-structure sensitivity validation and main/True-VC baseline cross-check: PASS.
- Combined production-output validator: PASS.
- Mechanism diagnostics and selection-intensity sensitivity checks: completed and internally validated.
- Definitive `sessionInfo()` archived as `sessionInfo_final.txt`.
- SHA256 manifest generated for the release tree after environment archival.

## Before uploading the release

1. Decide the software license. Do not imply a license until a `LICENSE` file has actually been added.

2. Create the GitHub repository/release. Suggested tag: `v2026.09.08-r4` (or use the project's preferred semantic versioning scheme).

3. Upload the exact final ZIP and/or tag the exact repository commit corresponding to this archive.

4. Connect the GitHub repository to Zenodo, create an archived release, and record the version DOI and concept DOI.

5. Replace the placeholders in `CODE_AVAILABILITY_TEMPLATE.md` with the real GitHub URL and Zenodo DOI. Add the final manuscript citation/DOI when assigned.

6. Preserve the final ZIP SHA256 outside the archive (for example in release notes or an internal submission record).

## Do not change after freeze without revalidation

Do not modify REML, EBLUP, policy construction, RAM-N, RL-UP integration, bootstrap, True-VC, same/cross, staged-decomposition, selection-intensity, or pedigree scientific calculations after this freeze. If any R code changes, rerun at least the parse check and the relevant small tests; if a scientific calculation changes, treat it as a new analysis version and repeat the corresponding production validation.
