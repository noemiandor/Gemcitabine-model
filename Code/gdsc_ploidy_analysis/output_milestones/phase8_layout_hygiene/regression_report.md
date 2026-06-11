# Regression Report

- Baseline: `Code/gdsc_ploidy_analysis/output_milestones/phase7_dependency_metadata`
- Candidate: `Code/gdsc_ploidy_analysis/output_milestones/phase8_layout_hygiene`
- Final result changed: `no scientific-result change; layout and generated-binary changes only`
- Difference count: `1`

## Difference Summary

- `checksums.tsv`: changed (md5 962397586d2945fb892545251c40bc6f -> 41a620f4d519ba8f2b1a5b4609d9e61c)

## Attribution

The first and only tracked difference is `checksums.tsv`. `drugsVsPloidyCorr.RData`, `coxIn.RData`, exported correlations, saved annotations, group counts, duplicate QC files, and workbook sheet TSV exports match Phase 7.

Expected differences:

- `input_manifest.tsv` changed because input files moved from the flat layout into `data/raw/`, `data/manual/`, and `data/derived/`.
- Excel/PDF hashes changed because binary outputs were regenerated.
- Maintained code moved under `src/`, historical scripts moved under `references/`, and root wrapper scripts preserved the existing `Rscript Code/gdsc_ploidy_analysis/run_gdsc_ploidy_analysis.R` and `Rscript Code/gdsc_ploidy_analysis/refresh_pubchem_annotations.R` commands.
- Generated `output/`, generated milestone internals, and local `.Rlibs/` are now ignored by Git. Previously tracked generated outputs and `.Rlibs` were removed from the Git index with files left on disk.
