# Regression Report

- Baseline: `Code/gdsc_ploidy_analysis/output_milestones/phase4_explicit_error_handling`
- Candidate: `Code/gdsc_ploidy_analysis/output_milestones/phase5_metric_labels_required_groups`
- Final result changed: `expected presentation/QC change only; numeric scientific intermediates unchanged`
- Difference count: `1`

## Difference Summary

- `checksums.tsv`: changed (md5 7836a5f6b7a8e14c6397c4c2eb2c401f -> ba9f196f83916f0378c4d434252321f0)

## Attribution

The first and only tracked difference is `checksums.tsv`. `drugsVsPloidyCorr.RData`, `coxIn.RData`, exported correlations, saved annotations, group counts, and workbook sheet TSV exports match Phase 4.

Expected differences:

- `drugsVsPloidyCorr.pdf` changed because scatter plots now label the x-axis as `GDSC Z-score` instead of relying on an unlabeled metric axis.
- `ploidyVsDrugSensitivity.pdf` changed because the barplot x-axis now states `Pearson r between ploidy and drug sensitivity (GDSC Z-score)` instead of incorrectly saying `IC50`.
- `drugsVsPloidyCorr.xlsx` was regenerated, but the exported workbook sheet TSVs match Phase 4.
- New output metadata/QC files were added: `run_config.tsv`, `drug_category_counts_before_filter.tsv`, and `drug_category_counts_after_filter.tsv`.

Required-group validation passed because both `SIGNALING` and `CYTOTOXIC` are present after category filtering.
