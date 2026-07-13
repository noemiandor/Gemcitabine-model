# Essential pseudotime-TGI figure, data, and statistics index

This directory is a tracked snapshot of the essential CellCycle pseudotime-TGI results.

## Source

The snapshot is generated from the standalone workflow:

```text
Code/in-vivo/04h_pseudotime_TGI_essential.R
Code/in-vivo/04h_pseudotime_TGI_essential_util.R
```

The corresponding full result directory is:

```text
/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04h_pseudotime_TGI_essential
```

The full analysis can start from explicitly supplied CellCycle and NonCellCycle cell-level CSV files. If either path is omitted, the entrypoint reconstructs the missing table from these repository-local inputs:

```text
Data/in-vivo/scvelo_cell_metrics.csv
Data/in-vivo/all_ploidy.tsv
Data/in-vivo/sample_info.xlsx
Data/in-vivo/dt_Gem_VT_20241223_v4.xlsx
```

By default, reconstructed cell-level tables are written under `Data/in-vivo`. A different destination can be selected with `--generated_input_root=<directory>`. The NonCellCycle input is used only when deriving sample mean end-timepoint ploidy; all figures and reported associations use CellCycle cells.

The raw-input paths can be overridden individually with `--scvelo_metrics_input`, `--cell_ploidy_input`, `--sample_info_input`, and `--growth_curve_input`. Supplying both `--cellcycle_input` and `--noncellcycle_input` skips reconstruction and uses those finished tables directly.

## Methods

The same six figures, plotting-data tables, and statistical outputs are provided for four peer-level grouping methods:

| Method directory | Grouping definition |
|---|---|
| `initial_ploidy` | Original sample initial ploidy group (`2N` or `4N`) |
| `ETP_fixed_threshold_2_25` | End-timepoint mean ploidy threshold of 2.25 |
| `ETP_boundary_stress_threshold_2_375` | End-timepoint mean ploidy boundary-stress threshold of 2.375 |
| `ETP_reference_balanced_threshold_2_24` | Untreated-reference-balanced end-timepoint mean ploidy threshold of 2.24 |

## Directory structure

```text
pseudotime_TGI/
├── Figures/
│   └── <method>/
├── stats/
│   └── <method>/
├── plot_data/
│   └── <method>/
└── README.md
```

Each method contains six PDF figures, thirteen statistical CSV files, and six plotting-data CSV files.

## Figure-to-table mapping

| Figure | Plotting-data table | Primary statistical file(s) | Supporting statistical file(s) |
|---|---|---|---|
| `CellCycle_direct_group_ecdf_comparisons.pdf` | `CellCycle_direct_group_ecdf_comparisons_plot_data.csv` | `CellCycle_direct_group_ecdf_comparisons_9panel_tests.csv` | None |
| `CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref.pdf` | `CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref_plot_data.csv` | `CellCycle_TGI_associations_ecdf_rmse.csv` | `CellCycle_primary_TGI_robustness_summary.csv`; `CellCycle_primary_TGI_leave_one_out.csv`; `CellCycle_primary_TGI_bootstrap.csv` |
| `CellCycle_TGI_AUC_vs_mean_ETP.pdf` | `CellCycle_TGI_AUC_vs_mean_ETP_plot_data.csv` | `CellCycle_TGI_associations_mean_ETP.csv` | `CellCycle_mean_ETP_TGI_robustness_summary.csv`; `CellCycle_mean_ETP_TGI_leave_one_out.csv`; `CellCycle_mean_ETP_TGI_bootstrap.csv` |
| `CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose.pdf` | `CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose_plot_data.csv` | `CellCycle_AUC_TGI_shift_ploidy_dose_models.csv`; `CellCycle_AUC_TGI_shift_ploidy_dose_model_summaries.csv` | `CellCycle_TGI_associations_ecdf_rmse.csv` supplies the displayed primary correlation and permutation P value. |
| `CellCycle_TGI_association_within_dose_centered.pdf` | `CellCycle_TGI_association_within_dose_centered_plot_data.csv` | `residualized_TGI_associations.csv` | None |
| `CellCycle_ecdf_rmse_vs_ploidy.pdf` | `CellCycle_ecdf_rmse_vs_ploidy_plot_data.csv` | `ploidy_confounding_tests.csv` | None |

Each plotting-data CSV contains the exact rows used by its PDF. Scatter-plot tables also retain the Pearson and Spearman statistics, asymptotic and permutation P values, permutation mode, permutation count, sample count, and annotation text printed in the figure.

## Primary-row filters

For `CellCycle_TGI_associations_ecdf_rmse.csv`, the primary row displayed in the equal-sample-reference association figures is identified by:

```text
compartment = CellCycle
sample_set = treated
reference_type = primary_equal_sample_reference
shift_metric = ecdf_rmse
tgi_measure = TGI_percent_auc
pre_specified_primary = TRUE
```

For `CellCycle_TGI_associations_mean_ETP.csv`, the AUC-based mean-ETP row is identified by:

```text
compartment = CellCycle
sample_set = treated
predictor_label = sample_mean_ETP
tgi_measure = TGI_percent_auc
```

For `residualized_TGI_associations.csv`, the displayed within-dose-centered result is identified by:

```text
compartment = CellCycle
analysis = within_dose_centered
method = pearson
```

For `ploidy_confounding_tests.csv`, the displayed all-sample ploidy association is identified by:

```text
compartment = CellCycle
sample_set = all
shift_metric = ecdf_rmse
ploidy_measure = mean_cell_ploidy
```

## Figures-only regeneration

The figures can be regenerated directly from the tracked plotting-data and statistical tables without the cell-level inputs and without rerunning statistical tests:

```bash
Rscript Code/in-vivo/04h_pseudotime_TGI_essential.R \
  --figures_only=TRUE \
  --tables_root='Figs/pseudotime_TGI' \
  --output_root='Figs/pseudotime_TGI' \
  --methods=all \
  --workers=4 \
  --overwrite=TRUE
```

In figures-only mode, the workflow reads only `plot_data/` and the required files in `stats/`. It replaces only the selected `Figures/<method>/` directories and does not modify `stats/`, `plot_data/`, or this README.

## Inventory

- Four method directories
- Six PDF figures per method
- Thirteen statistical CSV files per method
- Six plotting-data CSV files per method
- Twenty-four PDF files in total
- Fifty-two statistical CSV files in total
- Twenty-four plotting-data CSV files in total
- No PNG, PPT, or PPTX files
