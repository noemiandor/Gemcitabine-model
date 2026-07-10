# Pseudotime-TGI figure and statistics index

This directory is a curated snapshot of selected CellCycle pseudotime-TGI results.

## Source

All files were copied from:

```text
/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04h_pseudotime_TGI
```

The result directory above remains the source of truth. Files in this directory are organized copies for figure review and downstream use.

## Methods

The same six figures and their associated statistical outputs are provided for four peer-level grouping methods:

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
│   ├── initial_ploidy/
│   ├── ETP_fixed_threshold_2_25/
│   ├── ETP_boundary_stress_threshold_2_375/
│   └── ETP_reference_balanced_threshold_2_24/
├── stats/
│   ├── initial_ploidy/
│   ├── ETP_fixed_threshold_2_25/
│   ├── ETP_boundary_stress_threshold_2_375/
│   └── ETP_reference_balanced_threshold_2_24/
└── README.md
```

Each method contains six PDF figures and thirteen CSV statistical files.

## Figure-to-statistics mapping

| Figure | Primary statistical file(s) | Supporting statistical file(s) | Interpretation of the matching result |
|---|---|---|---|
| `CellCycle_direct_group_ecdf_comparisons.pdf` | `CellCycle_direct_group_ecdf_comparisons_9panel_tests.csv` | None | The CSV contains the permutation tests and effect sizes for all nine ECDF comparison panels. |
| `CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref.pdf` | `CellCycle_TGI_associations_ecdf_rmse.csv` | `CellCycle_primary_TGI_robustness_summary.csv`; `CellCycle_primary_TGI_leave_one_out.csv`; `CellCycle_primary_TGI_bootstrap.csv` | The plotted primary association uses treated CellCycle samples, the equal-sample reference, ECDF RMSE, and AUC-based TGI. The three supporting files report robustness, leave-one-out influence, and bootstrap uncertainty. |
| `CellCycle_TGI_AUC_vs_mean_ETP.pdf` | `CellCycle_TGI_associations_mean_ETP.csv` | `CellCycle_mean_ETP_TGI_robustness_summary.csv`; `CellCycle_mean_ETP_TGI_leave_one_out.csv`; `CellCycle_mean_ETP_TGI_bootstrap.csv` | The plotted association uses treated CellCycle samples, sample mean ETP as the predictor, and AUC-based TGI as the response. The supporting files report the corresponding robustness analyses. |
| `CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose.pdf` | `CellCycle_AUC_TGI_shift_ploidy_dose_models.csv`; `CellCycle_AUC_TGI_shift_ploidy_dose_model_summaries.csv` | None | These files contain the shift-only model, the ploidy- and dose-adjusted model, and the shift-by-ploidy interaction model shown by the grouped scatter plot context. |
| `CellCycle_TGI_association_within_dose_centered.pdf` | `residualized_TGI_associations.csv` | None | Use the row with `compartment = CellCycle` and `analysis = within_dose_centered` for the correlation displayed in the figure. |
| `CellCycle_ecdf_rmse_vs_ploidy.pdf` | `ploidy_confounding_tests.csv` | None | Use rows with `compartment = CellCycle`, `sample_set = all`, `shift_metric = ecdf_rmse`, and `ploidy_measure = mean_cell_ploidy`; Pearson and Spearman results are both reported. |

## Primary-row filters

For `CellCycle_TGI_associations_ecdf_rmse.csv`, the primary row plotted in `CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref.pdf` is identified by:

```text
compartment = CellCycle
sample_set = treated
reference_type = primary_equal_sample_reference
shift_metric = ecdf_rmse
tgi_measure = TGI_percent_auc
pre_specified_primary = TRUE
```

For `CellCycle_TGI_associations_mean_ETP.csv`, the AUC-based mean-ETP row plotted in `CellCycle_TGI_AUC_vs_mean_ETP.pdf` is identified by:

```text
compartment = CellCycle
sample_set = treated
predictor_label = sample_mean_ETP
tgi_measure = TGI_percent_auc
```

## Normalized association filename

The ECDF RMSE association table had a historical engine-specific filename in the source results:

| Method | Source filename | Curated filename |
|---|---|---|
| `initial_ploidy` | `CellCycle_TGI_associations_v5.csv` | `CellCycle_TGI_associations_ecdf_rmse.csv` |
| All three ETP methods | `CellCycle_TGI_associations_v4.csv` | `CellCycle_TGI_associations_ecdf_rmse.csv` |

Only the copied filename was normalized. The CSV contents and source files were not modified.

## Inventory

- Four method directories
- Six PDF figures per method
- Thirteen CSV statistical files per method
- Twenty-four PDF files in total
- Fifty-two CSV files in total
- No PNG, PPT, or PPTX files
