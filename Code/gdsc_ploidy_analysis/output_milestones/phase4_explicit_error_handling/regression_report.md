# Regression Report

- Baseline: `Code/gdsc_ploidy_analysis/output_milestones/phase3_structured_pubchem`
- Candidate: `Code/gdsc_ploidy_analysis/output_milestones/phase4_explicit_error_handling`
- Final result changed: `no scientific-result change; regenerated binary artifacts differ`
- Difference count: `1`

## Difference Summary

- `checksums.tsv`: changed (md5 2be450f745dd5598261ce88b615b2fcb -> 7836a5f6b7a8e14c6397c4c2eb2c401f)

## Attribution

The first and only tracked difference is `checksums.tsv`. `drugsVsPloidyCorr.RData`, `coxIn.RData`, exported correlations, saved annotations, group counts, and workbook sheet TSV exports match Phase 3.

The changed hashes are limited to regenerated `drugsVsPloidyCorr.xlsx`, `drugsVsPloidyCorr.pdf`, and `ploidyVsDrugSensitivity.pdf`. Replacing silent `try()` calls with explicit helpers did not change successful computations.

Failure-path tests passed with `Rscript Code/gdsc_ploidy_analysis/tests/test_error_handling.R`, confirming enrichment failures include cancer, direction, drug count, group count, and the original error, and plotting failures include cancer and plotted drug count.
