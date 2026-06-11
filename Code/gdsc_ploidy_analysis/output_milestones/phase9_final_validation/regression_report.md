# Regression Report

- Baseline: `Code/gdsc_ploidy_analysis/output_milestones/phase8_layout_hygiene`
- Candidate: `Code/gdsc_ploidy_analysis/output_milestones/phase9_final_validation`
- Final result changed: `no scientific-result change; regenerated binary artifacts differ`
- Difference count: `1`

## Difference Summary

- `checksums.tsv`: changed (md5 41a620f4d519ba8f2b1a5b4609d9e61c -> 0639f13fe82bef052f2505afd772d253)

## Attribution

The first and only tracked difference is `checksums.tsv`. `drugsVsPloidyCorr.RData`, `coxIn.RData`, exported correlations, saved annotations, group counts, duplicate QC files, input manifest, session info, and workbook sheet TSV exports match Phase 8.

The changed hashes are limited to regenerated `drugsVsPloidyCorr.xlsx`, `drugsVsPloidyCorr.pdf`, and `ploidyVsDrugSensitivity.pdf`. These binary artifacts are expected to change across reruns even when the exported scientific intermediates are identical.

Final validation commands run:

- `Rscript Code/gdsc_ploidy_analysis/tests/smoke_test.R`
- `Rscript Code/gdsc_ploidy_analysis/run_gdsc_ploidy_analysis.R --output-dir=Code/gdsc_ploidy_analysis/output_milestones/phase9_final_validation/outputs`
- `Rscript Code/gdsc_ploidy_analysis/tests/create_analysis_baseline.R --baseline-dir=Code/gdsc_ploidy_analysis/output_milestones/phase9_final_validation --output-dir=Code/gdsc_ploidy_analysis/output_milestones/phase9_final_validation/outputs`
- `Rscript Code/gdsc_ploidy_analysis/tests/compare_analysis_outputs.R --baseline-dir=Code/gdsc_ploidy_analysis/output_milestones/phase8_layout_hygiene --candidate-dir=Code/gdsc_ploidy_analysis/output_milestones/phase9_final_validation --allow-differences=true`
