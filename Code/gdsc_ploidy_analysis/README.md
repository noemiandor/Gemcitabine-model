# GDSC Ploidy Analysis

This folder packages the GDSC-versus-ploidy analysis used in the gemcitabine manuscript into a self-contained module.

## Contents

- `run_gdsc_ploidy_analysis.R`: backward-compatible wrapper for the main entrypoint
- `refresh_pubchem_annotations.R`: backward-compatible wrapper for the network-dependent annotation refresh entrypoint
- `src/`: maintained analysis, annotation, dependency, PubChem, plotting, and validation code
- `data/raw/`: raw local input files required by the analysis
- `data/manual/`: manual category override inputs
- `data/derived/`: reviewed derived inputs, including the cached PubChem annotation table
- `references/`: historical scripts retained for provenance, not maintained entrypoints
- `output/`: generated results, ignored by Git

## Inputs

The analysis uses:

- `data/raw/GDSC2_fitted_dose_response_24Jul22.txt`
- `data/raw/ploidyAcrossCellLines_V1.txt`
- `data/raw/small_molecule_20200407234909.csv`
- `data/manual/custom_set_candidate.tsv`
- `data/manual/drug_class_category_schema.tsv`
- `data/manual/drug_class_final_curated.tsv`
- `data/derived/pubchem_drug_annotations.tsv`
- `tests/fixtures/wrong_assignment_examples.tsv`

## Outputs

Running the analysis writes:

- `output/drugsVsPloidyCorr.pdf`
- `output/drugsVsPloidyCorr.RData`
- `output/coxIn.RData`
- `output/drugsVsPloidyCorr.xlsx`
- `output/drugsVsPloidyCorr_curated_Z_SCORE.xlsx`
- `output/ploidy_enrichment_panels_curated_ABC.png`
- `output/ploidy_enrichment_panels_curated_ABC.pdf`
- `output/ploidy_enrichment_clustered_curated_lowpIsSens_clustermap.png`
- `output/ploidy_enrichment_clustered_curated_lowpIsSens_clustermap.pdf`
- `output/ploidy_enrichment_clustered_curated_highpIsSens_clustermap.png`
- `output/ploidy_enrichment_clustered_curated_highpIsSens_clustermap.pdf`
- `output/ploidy_enrichment_clustered_curated_shared_order.png`
- `output/ploidy_enrichment_clustered_curated_shared_order.pdf`
- `output/ploidyVsDrugSensitivity.pdf`
- `output/metadata/run_config.tsv`
- `output/metadata/run_parameters.tsv`
- `output/metadata/category_mode_artifacts.tsv`
- `output/metadata/input_manifest.tsv`
- `output/metadata/session_info.txt`
- `output/tables/cell_line_key_collisions_gdsc.tsv`
- `output/tables/cell_line_key_collisions_ploidy.tsv`
- `output/tables/cell_line_matching_delta_raw_vs_normalized.tsv`
- `output/tables/drug_ploidy_correlations_by_cancer_Z_SCORE.tsv`
- `output/tables/drug_ploidy_correlations_by_cancer_Z_SCORE.xlsx`
- `output/tables/drug_ploidy_correlations_by_cancer_all_metrics.tsv`
- `output/tables/drug_ploidy_correlations_by_cancer_all_metrics.xlsx`
- `output/tables/drug_ploidy_correlations_by_cancer_Z_SCORE_normalized_key_SUPPLEMENTAL.tsv`
- `output/tables/drug_ploidy_correlations_by_cancer_Z_SCORE_normalized_key_SUPPLEMENTAL.xlsx`
- `output/tables/correlation_delta_raw_vs_normalized_Z_SCORE.tsv`
- `output/tables/drug_annotation_audit.tsv`
- `output/tables/drug_class_final_curated_TEMPLATE.tsv`
- `output/tables/drug_class_correlation_eligible_drugs.tsv`
- `output/tables/drug_class_curation_evidence_by_drug_id.tsv`
- `output/tables/drug_class_curation_input.tsv`
- `output/tables/drug_class_final_used.tsv`
- `output/tables/drug_class_assignment_diff.tsv`
- `output/tables/drug_class_unassigned_failures.tsv`
- `output/tables/drug_class_category_counts.tsv`
- `output/tables/drug_class_reviewed_exclusions.tsv`
- `output/tables/drug_class_low_count_actions.tsv`
- `output/tables/class_enrichment_curated_Z_SCORE.tsv`
- `output/tables/class_enrichment_curated_metadata.tsv`
- `output/tables/class_enrichment_selected_drugs_curated_Z_SCORE.tsv`
- `output/tables/drugsVsPloidyCorr_curated_Z_SCORE.xlsx`
- `output/tables/ploidyVsDrugSensitivity_plot_values_Z_SCORE.tsv`
- `output/tables/ploidyVsDrugSensitivity_pages/ploidyVsDrugSensitivity_page_<NN>_<cancer_type>.tsv`
- `output/tables/gemcitabine_rank_summary_Z_SCORE.tsv`
- `output/tables/gemcitabine_rank_summary_all_metrics.tsv`
- `output/tables/all_cancers_tissue_adjusted_ploidy_models_Z_SCORE.tsv`
- `output/tables/all_cancers_tissue_adjusted_ploidy_models_all_metrics.tsv`
- `output/qc/duplicate_drug_cell_line_records.tsv`
- `output/qc/duplicate_drug_cell_line_selection.tsv`
- `output/qc/duplicate_resolution_summary.tsv`

The default sensitivity metric is `Z_SCORE`, labeled in figures as `GDSC Z-score`. Positive correlations indicate higher ploidy is associated with higher values of this metric; negative correlations indicate higher ploidy is associated with lower values of this metric.

## R package dependencies

- `EnrichIntersect`
- `xlsx`
- `plyr`
- `openxlsx`
- `data.table`
- `ChemmineR`
- `jsonlite`
- `httr2`

The main analysis requires `EnrichIntersect`, `xlsx`, `plyr`, `openxlsx`, and `data.table`. Annotation refresh additionally requires `ChemmineR`, `jsonlite`, and `httr2`. If `Code/gdsc_ploidy_analysis/.Rlibs` exists locally, it is prepended to `.libPaths()` before package checks. `.Rlibs` is ignored by Git.

## Commands

From the repository root:

```sh
Rscript Code/gdsc_ploidy_analysis/tests/smoke_test.R
Rscript Code/gdsc_ploidy_analysis/run_gdsc_ploidy_analysis.R
```

The main run accepts `--analysis-mode=dev` or `--analysis-mode=manuscript`. `dev` is the default and uses 300 enrichment permutations; `manuscript` uses 10000 unless overridden with `--enrichment-permute-n=<n>`. Category assignment defaults to `--category-mode=curated`, which requires `data/manual/drug_class_category_schema.tsv` and `data/manual/drug_class_final_curated.tsv`. Manuscript mode refuses `legacy` or `proposal` category modes. The selected mode, category inputs, and checksums are recorded in `output/metadata/run_parameters.tsv` and `output/metadata/category_mode_artifacts.tsv`.

To run an isolated full validation:

```sh
Rscript Code/gdsc_ploidy_analysis/run_gdsc_ploidy_analysis.R \
  --output-dir=Code/gdsc_ploidy_analysis/output_milestones/manual_validation/outputs

Rscript Code/gdsc_ploidy_analysis/tests/create_analysis_baseline.R \
  --baseline-dir=Code/gdsc_ploidy_analysis/output_milestones/manual_validation \
  --output-dir=Code/gdsc_ploidy_analysis/output_milestones/manual_validation/outputs

Rscript Code/gdsc_ploidy_analysis/tests/compare_analysis_outputs.R \
  --baseline-dir=Code/gdsc_ploidy_analysis/output_milestones/phase8_layout_hygiene \
  --candidate-dir=Code/gdsc_ploidy_analysis/output_milestones/manual_validation \
  --allow-differences=true
```

The comparison report identifies whether final outputs changed and where the first meaningful divergence occurred. Regenerated Excel/PDF hashes can differ even when RData and exported TSV intermediates match.

## Notes

- `data/manual/custom_set_candidate.tsv` is a reconstructed manual override table inferred from surviving outputs and rerun behavior.
- The original script relied on an undefined in-memory object called `custom.set`; this package makes that dependency explicit.
- The main analysis uses `data/derived/pubchem_drug_annotations.tsv` as raw evidence by default and does not contact PubChem during normal runs. Final manuscript categories come from `data/manual/drug_class_final_curated.tsv`, validated against `data/manual/drug_class_category_schema.tsv`.
- To refresh PubChem annotations, run `Rscript refresh_pubchem_annotations.R --output-cache=data/derived/pubchem_drug_annotations.refresh.tsv` from `Code/gdsc_ploidy_analysis/` and review the TSV diff before replacing the canonical cache. When a drug subset is requested, the refresh preserves unrequested cached rows by default; use `--subset-output` only when a subset-only cache is intentional.
- PubChem categories are raw evidence only. The curated workflow does not infer final biology from the last semicolon-separated PubChem category, does not use broad administrative labels such as `Antineoplastic Agents` or `Signal Transduction Inhibitors` as manuscript final classes, and orders enrichment result matrices from the reviewed category schema.
- Duplicate `DRUG_NAME`/`CELL_LINE_NAME` records are resolved with `lowest_rmse` by default, which requires an `RMSE` column and uses stable identifier tie-breakers. The duplicate audit, selected-row audit, and resolution summary are written under `output/qc/`.
- `references/annotate_from_pubchem_legacy.R` is retained as a historical near-verbatim port of the original script. The refresh entrypoint uses structured PubChem JSON parsing in `src/pubchem_client.R`.
- `matlab` is no longer required by this module because `isempty()` is implemented locally.
