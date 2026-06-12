# CCLE Ploidy Analysis

This folder packages the CCLE breast-cell-line drug sensitivity versus ploidy barplot analysis into a self-contained module.

The default entrypoint reproduces the historical IC50-labeled plot from `BREAST_CCLE_PathwaysPloidy_depracated.R`. That historical plot reads `BreastCancerDrugSensitivity` files and correlates ploidy with their `Z Score` column, while labeling the x-axis as `Pearson (ploidy, IC50)`.

## Contents

- `run_ccle_ploidy_analysis.R`: command-line entrypoint
- `src/`: maintained analysis, plotting, metadata, and helper code
- `data/raw/`: local input files required by the analysis
- `data/manual/`: manually curated breast CCLE ploidy table from the historical scripts
- `data/derived/`: reduced derived inputs needed for filtering, including CCLE expression column names
- `references/`: historical scripts retained for provenance
- `tests/`: smoke and baseline comparison scripts
- `baseline/`: reference outputs for the default IC50-labeled analysis
- `output/`: generated results, ignored by Git if the parent repository ignore rules include it

## Inputs

The default IC50-labeled analysis uses:

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

- `output/ccle_drug_ploidy_correlations_ic50.pdf`
- `output/ccle_drug_ploidy_correlations_ic50.tsv`
- `output/drug_ploidy_correlations_ic50.RData`
- `output/result_summary.tsv`
- `output/metadata/run_config.tsv`
- `output/metadata/run_parameters.tsv`
- `output/metadata/input_manifest.tsv`
- `output/metadata/session_info.txt`
- `output/tables/drug_coverage_ic50.tsv`
- `output/tables/drug_ploidy_correlations_all_ic50.tsv`
- `output/tables/drug_ploidy_correlations_plotted_ic50.tsv`

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

To compare current outputs with the checked baseline:

```sh
Rscript Code/ccle_ploidy_analysis/tests/compare_analysis_outputs.R
```

## Notes

- The default IC50-labeled analysis intentionally follows `BREAST_CCLE_PathwaysPloidy_depracated.R`, including the `0.8 * max(drug coverage)` drug filter and `abs(correlation) >= 0.2` plotting threshold.
- Positive correlations in the default IC50-labeled output are colored as low-ploidy sensitive, matching the historical legend.
- `BREAST_CCLE_PathwaysPloidy.R` is retained as provenance for the later MEP-LINCS `GR_AOC` plot, while `BREAST_CCLE_PathwaysPloidy_depracated.R` is retained as provenance for the default IC50-labeled plot.
