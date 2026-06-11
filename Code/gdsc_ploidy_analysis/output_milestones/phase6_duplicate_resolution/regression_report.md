# Regression Report

- Baseline: `Code/gdsc_ploidy_analysis/output_milestones/phase5_metric_labels_required_groups`
- Candidate: `Code/gdsc_ploidy_analysis/output_milestones/phase6_duplicate_resolution`
- Final result changed: `yes; intentional duplicate-resolution scientific-result change`
- Difference count: `6`

## Difference Summary

- `checksums.tsv`: changed (md5 ba9f196f83916f0378c4d434252321f0 -> b9fa5ce4c14dfaf7b1e0a64d495b5ce5)
- `intermediate/coxIn_saved.tsv`: changed (md5 ed2799c1cb3dc755b705ac66150c4eeb -> 138d86e2940a9d6ec2fa4d282f79b046)
- `intermediate/drug_ploidy_correlations.tsv`: changed (md5 6ffc78a0971bd5f55e63982147550732 -> 788019033579a35fc296ebc0e721be27)
- `intermediate/workbook_highpIsSens.tsv`: changed (md5 27207692372cb6d346f7bfb605af5f94 -> b390eee98800aefd0ac2a8bc6985a061)
- `intermediate/workbook_lowpIsSens.tsv`: changed (md5 04371dc1bd1c2f200fa90e65d320b836 -> 8621551558e44021d0a06630f7a6bfc6)
- `result_summary.tsv`: changed (md5 674078a31d3d0e2bccd00bd1bce3b7d3 -> 45c9b729ca161865aa32db2e6e42dbe3)

## Attribution

First meaningful divergence: `intermediate/drug_ploidy_correlations.tsv`.

Implementation cause: duplicate `DRUG_NAME`/`CELL_LINE_NAME` rows are now resolved once, before correlation analysis, using `lowest_rmse` with stable identifier tie-breakers. Previously the script kept the first row encountered inside each drug/cancer subset, making results dependent on input row order.

Duplicate audit:

- Input rows: 242,036
- Duplicate rows: 13,232
- Duplicate drug/cell-line keys: 6,611
- Removed rows after deterministic resolution: 6,621
- Strategy: `lowest_rmse`

Observed result impact:

- Correlation rows compared: 6,161
- Correlations changed: 569
- Largest observed correlation change: `LIHC` / `Oxaliplatin`, baseline `0.10842299`, candidate `-0.64512225`, delta `-0.7535452`
- `lowpIsSens` workbook changed cells: 217
- `highpIsSens` workbook changed cells: 218

`coxIn_saved.tsv` changed by row order because the correlation result ordering changed; the drug set stayed fixed at 286 rows with zero added or removed drugs.

This is an accepted scientific-result change because `lowest_rmse` is a deterministic quality-based rule and removes the previous input-order dependency.
