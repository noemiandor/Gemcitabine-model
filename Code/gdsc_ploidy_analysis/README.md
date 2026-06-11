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
- `data/derived/pubchem_drug_annotations.tsv`

## Outputs

Running the analysis writes:

- `output/drugsVsPloidyCorr.pdf`
- `output/drugsVsPloidyCorr.RData`
- `output/coxIn.RData`
- `output/drugsVsPloidyCorr.xlsx`
- `output/ploidyVsDrugSensitivity.pdf`
- `output/metadata/run_config.tsv`
- `output/metadata/input_manifest.tsv`
- `output/metadata/session_info.txt`
- `output/tables/drug_category_counts_before_filter.tsv`
- `output/tables/drug_category_counts_after_filter.tsv`
- `output/qc/duplicate_drug_cell_line_records.tsv`
- `output/qc/duplicate_resolution_summary.tsv`

The default sensitivity metric is `Z_SCORE`, labeled in figures as `GDSC Z-score`. Positive correlations indicate higher ploidy is associated with higher values of this metric; negative correlations indicate higher ploidy is associated with lower values of this metric.

## R package dependencies

- `EnrichIntersect`
- `xlsx`
- `plyr`
- `ChemmineR`
- `jsonlite`
- `httr2`

The main analysis requires `EnrichIntersect`, `xlsx`, and `plyr`. Annotation refresh additionally requires `ChemmineR`, `jsonlite`, and `httr2`. If `Code/gdsc_ploidy_analysis/.Rlibs` exists locally, it is prepended to `.libPaths()` before package checks. `.Rlibs` is ignored by Git.

## Commands

From the repository root:

```sh
Rscript Code/gdsc_ploidy_analysis/tests/smoke_test.R
Rscript Code/gdsc_ploidy_analysis/run_gdsc_ploidy_analysis.R
```

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
- The main analysis uses `data/derived/pubchem_drug_annotations.tsv` by default and does not contact PubChem during normal runs.
- To refresh PubChem annotations, run `Rscript refresh_pubchem_annotations.R --output-cache=data/derived/pubchem_drug_annotations.refresh.tsv` from `Code/gdsc_ploidy_analysis/` and review the TSV diff before replacing the canonical cache.
- `SIGNALING` and `CYTOTOXIC` are required after category normalization and filtering because they are used to order enrichment result matrices.
- Duplicate `DRUG_NAME`/`CELL_LINE_NAME` records are resolved with `lowest_rmse` by default, with stable identifier tie-breakers. The duplicate audit and resolution summary are written under `output/qc/`.
- `references/annotate_from_pubchem_legacy.R` is retained as a historical near-verbatim port of the original script. The refresh entrypoint uses structured PubChem JSON parsing in `src/pubchem_client.R`.
- `matlab` is no longer required by this module because `isempty()` is implemented locally.
