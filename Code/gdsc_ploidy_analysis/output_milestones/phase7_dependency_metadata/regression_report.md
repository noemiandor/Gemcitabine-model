# Regression Report

- Baseline: `Code/gdsc_ploidy_analysis/output_milestones/phase6_duplicate_resolution`
- Candidate: `Code/gdsc_ploidy_analysis/output_milestones/phase7_dependency_metadata`
- Final result changed: `no scientific-result change; metadata outputs added`
- Difference count: `1`

## Difference Summary

- `checksums.tsv`: changed (md5 b9fa5ce4c14dfaf7b1e0a64d495b5ce5 -> 8ea1b98df67e3b1ea7e8212a4bffad3f)

## Attribution

The first and only tracked difference is `checksums.tsv`. `drugsVsPloidyCorr.RData`, `coxIn.RData`, exported correlations, saved annotations, group counts, duplicate QC files, and workbook sheet TSV exports match Phase 6.

Expected differences:

- `input_manifest.tsv` was added and records MD5 checksums for the GDSC dose-response input, ploidy input, small molecule metadata, manual overrides, and cached annotations.
- `session_info.txt` was added and records R, platform, locale, and attached package versions.
- Excel/PDF hashes changed because binary outputs were regenerated.

Startup dependency checks passed for `EnrichIntersect`, `xlsx`, and `plyr`, with local `.Rlibs` prepended when present.
