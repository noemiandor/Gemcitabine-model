# CCLE Ploidy Analysis

This folder packages the CCLE breast-cell-line drug sensitivity versus ploidy barplot analysis into a self-contained module.

The default entrypoint reproduces the historical `BREAST_CCLE_PathwaysPloidy_depracated.R` numeric result, but labels it by the actual response column: `BreastCancerDrugSensitivity Z Score`. Lower Z Score values are treated as more sensitive.

## Contents

- `run_ccle_ploidy_analysis.R`: command-line entrypoint
- `src/`: maintained analysis, plotting, metadata, and helper code
- `data/raw/`: local input files required by the analysis
- `data/manual/`: manually curated breast CCLE ploidy table from the historical scripts
- `data/derived/`: reduced derived inputs needed for filtering, including CCLE expression column names
- `references/`: historical scripts retained for provenance
- `tests/`: smoke and baseline comparison scripts
- `baseline/`: reference outputs for the default Z Score analysis
- `output/`: generated results, ignored by Git if the parent repository ignore rules include it

## Inputs

The default Z Score analysis uses:

- `data/raw/Cell_app_export.txt`
- `data/raw/DrugAliases.txt`
- `data/raw/breast_cancer_drug_sensitivity/*.tsv`
- `data/manual/breast_ccle_ploidy.tsv`
- `data/derived/ccle_expression_columns.tsv`

`data/derived/ccle_expression_columns.tsv` stores only the sample names from the original CCLE expression matrix because the drug-ploidy plot only needs expression availability for cell-line filtering, not expression values.

The optional grbrowser analysis also uses:

- `data/raw/grbrowser/*.tsv`

## Outputs

Running the default analysis writes:

- `output/ccle_drug_ploidy_correlations_z_score.pdf`
- `output/ccle_drug_ploidy_correlations_z_score.tsv`
- `output/drug_ploidy_correlations_z_score.RData`
- `output/result_summary.tsv`
- `output/metadata/run_config.tsv`
- `output/metadata/run_parameters.tsv`
- `output/metadata/input_manifest.tsv`
- `output/metadata/session_info.txt`
- `output/tables/drug_coverage_z_score.tsv`
- `output/tables/drug_ploidy_correlations_all_z_score.tsv`
- `output/tables/drug_ploidy_correlations_plotted_z_score.tsv`

## R Dependencies

The maintained analysis uses base R only. If `.Rlibs/` exists locally, it is prepended to `.libPaths()` for consistency with neighboring packaged analyses.

## Commands

From the repository root:

```sh
Rscript Code/ccle_ploidy_analysis/tests/smoke_test.R
Rscript Code/ccle_ploidy_analysis/run_ccle_ploidy_analysis.R
```

To write outputs to a custom directory:

```sh
Rscript Code/ccle_ploidy_analysis/run_ccle_ploidy_analysis.R \
  --output-dir=Code/ccle_ploidy_analysis/output
```

To run the optional grbrowser metric path:

```sh
Rscript Code/ccle_ploidy_analysis/run_ccle_ploidy_analysis.R \
  --metric=GR_AOC \
  --metric-source=grbrowser
```

To run true grbrowser IC50 values, use:

```sh
Rscript Code/ccle_ploidy_analysis/run_ccle_ploidy_analysis.R \
  --metric=IC50 \
  --metric-source=grbrowser
```

To compare current outputs with the checked baseline:

```sh
Rscript Code/ccle_ploidy_analysis/tests/compare_analysis_outputs.R
```

## Notes

- The default Z Score analysis intentionally follows `BREAST_CCLE_PathwaysPloidy_depracated.R`, including the `0.8 * max(drug coverage)` drug filter and `abs(correlation) >= 0.2` plotting threshold.
- Positive correlations in the default Z Score output are colored as low-ploidy sensitive because lower Z Score values are documented as more sensitive.
- `--metric=IC50 --metric-source=legacy` is accepted as a deprecated alias for the default Z Score analysis and emits a warning. Canonical default outputs use `z_score` filenames.
- `BREAST_CCLE_PathwaysPloidy.R` is retained as provenance for the later MEP-LINCS `GR_AOC` plot, while `BREAST_CCLE_PathwaysPloidy_depracated.R` is retained as provenance for the default Z Score plot.
