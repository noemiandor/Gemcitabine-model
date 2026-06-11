# Regression Report

- Baseline: `Code/gdsc_ploidy_analysis/output_milestones/phase1_cached_annotations`
- Candidate: `Code/gdsc_ploidy_analysis/output_milestones/phase2_refresh_entrypoint`
- Final result changed: `no scientific-result change; regenerated binary artifacts differ`
- Difference count: `1`

## Difference Summary

- `checksums.tsv`: changed (md5 c1700678a6bd33f20be7ed6add5487e8 -> 38fe01ed833fba60b4cc34712014d38f)

## Attribution

The first and only tracked difference is `checksums.tsv`. `drugsVsPloidyCorr.RData`, `coxIn.RData`, exported correlations, saved annotations, group counts, and workbook sheet TSV exports match Phase 1.

The changed hashes are limited to regenerated `drugsVsPloidyCorr.xlsx`, `drugsVsPloidyCorr.pdf`, and `ploidyVsDrugSensitivity.pdf`. This milestone added a separate refresh entrypoint and did not change the main analysis code path, so there is no scientific-result change.

The refresh entrypoint was tested in reuse-only mode with `--limit=3`; the refresh log showed all selected drugs were reused from the existing cache.
