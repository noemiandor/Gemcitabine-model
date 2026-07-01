# GDSC Ploidy Analysis

This module generates the GDSC drug-response versus ploidy analysis used for
Figure 1 and related public-data support panels.

## Inputs

The active workflow uses three local inputs:

- `data/raw/GDSC2_fitted_dose_response_24Jul22.txt`
- `data/raw/ploidyAcrossCellLines_V1.txt`
- `data/manual/drug_class_final_used_with_primary_secondary_corrected.xlsx`

The reviewed workbook is the only drug-class assignment source. Its
`primary_anticancer_class` column is the enrichment grouping variable.
`secondary_suggested_class` and other evidence columns are carried only as
metadata. PubChem caches, generated curation TSVs, legacy class schemas, and
manual fallback tables are no longer consumed by the main analysis.

## Outputs

Running the analysis writes a reproducible result directory containing:

- `drugsVsPloidyCorr_primary_secondary_Z_SCORE.xlsx`
- `drugsVsPloidyCorr.xlsx`, a compatibility copy of the active workbook
- `ploidy_enrichment_panels_primary_secondary_ABC.png` and `.pdf`
- `ploidy_enrichment_clustered_primary_secondary_lowpIsSens_clustermap.png`
  and `.pdf`
- `ploidy_enrichment_clustered_primary_secondary_highpIsSens_clustermap.png`
  and `.pdf`
- `ploidy_enrichment_clustered_primary_secondary_shared_order.png` and `.pdf`
- `ploidyVsDrugSensitivity.pdf`
- `metadata/run_config.tsv`
- `metadata/run_parameters.tsv`
- `metadata/category_mode_artifacts.tsv`
- `metadata/input_manifest.tsv`
- `metadata/session_info.txt`
- drug-ploidy correlation tables under `tables/`
- primary/secondary drug-class audit, count, validation, and review tables
  under `tables/`
- duplicate-resolution and cell-line matching QC tables under `qc/`

The canonical drug-response metric is `Z_SCORE`, labeled in figures as
`GDSC Z-score`. Positive correlations indicate higher ploidy is associated with
higher values of this metric; negative correlations indicate higher ploidy is
associated with lower values of this metric.

## R Dependencies

- `EnrichIntersect`
- `openxlsx`
- `data.table`

If `Code/gdsc_ploidy_analysis/.Rlibs` exists locally, it is prepended to
`.libPaths()` before package checks. `.Rlibs` is ignored by Git.

## Commands

From the repository root:

```sh
Rscript Code/gdsc_ploidy_analysis/tests/smoke_test.R
Rscript Code/gdsc_ploidy_analysis/run_gdsc_ploidy_analysis.R
```

The main run accepts:

- `--analysis-mode=dev|manuscript`
- `--drug-class-workbook=<path>`
- `--enrichment-permute-n=<n>`, default `300`
- `--output-dir=<path>`

`--category-mode=primary_secondary` is accepted for transitional compatibility.
Other category-mode values are intentionally unsupported.

For the root workflow:

```sh
bash Manager.sh --mode check-only --modules gdsc
bash Manager.sh --mode standard --modules gdsc --overwrite
```

## Notes

- Drug classes are reviewed anticancer mechanism classes based on the workbook,
  not raw PubChem categories.
- The workflow fails fast if the workbook has duplicate normalized drug keys,
  missing primary classes, mismatched class counts, or incomplete coverage of
  the correlation-eligible drug universe.
- Classes with five or fewer drugs are retained but flagged as exploratory in
  `tables/drug_class_primary_secondary_counts.tsv` and
  `tables/drug_class_primary_secondary_review_summary.tsv`.
- Duplicate `DRUG_NAME`/`CELL_LINE_NAME` records are resolved with
  `lowest_rmse`, which requires an `RMSE` column and uses stable identifier
  tie-breakers. The duplicate audit, selected-row audit, and resolution summary
  are written under `qc/`.
