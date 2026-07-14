# Essential pseudotime-TGI analysis

This directory was generated directly from the two 04h cell-level input CSV files.

The standalone workflow uses base R plus the `ggplot2` and `ggrepel` packages.

## Selected TGI outcome

- Measure: `TGI_percent_Day_17`.
- Outcome: `day`.
- Matched untreated-control summary: `mean`.
- TGI day: `17`.

TGI is calculated as `100 * (1 - mouse tumor-growth delta / matched-control reference delta)`. Controls are matched by initial ploidy. `--control-summary` chooses whether the matched untreated-control deltas are summarized by their mean, median, or maximum; it does not summarize the treated mice.

CLI options: `--tgi-outcome=auc|day`, `--control-summary=mean|median|max`, and (for the day outcome) `--tgi-day=<available day>`.

The established output filenames retain `AUC` where present for backward compatibility. The selected outcome is recorded in the statistical `tgi_measure` fields, in scatter-plot metadata, and in plot titles and axes.

## Contents

- `Figures/<method>/`: seven CellCycle PDF figures per method.
- `stats/<method>/`: thirteen statistical CSV files per method.
- `plot_data/<method>/`: seven plotting-data CSV files per method, one for each figure.

## Methods

- `initial_ploidy`: original 2N/4N sample grouping.
- `ETP_fixed_threshold_2_25`: sample mean ETP threshold 2.25.
- `ETP_boundary_stress_threshold_2_375`: sample mean ETP threshold 2.375.
- `ETP_reference_balanced_threshold_2_24`: sample mean ETP threshold 2.24.

## Figure data

Each plotting-data CSV contains the exact rows used by its PDF. Scatter-plot tables also contain the Pearson and Spearman statistics, permutation P values, permutation mode, and the annotation text printed in the figure.

`CellCycle_direct_group_ecdf_comparisons_selected_3panel.pdf` retains original grid positions `(1,1)`, `(3,2)`, and `(3,3)` (panels 1, 8, and 9). Its plotting-data CSV is a subset of the full direct-ECDF plotting table, and both figures reuse `CellCycle_direct_group_ecdf_comparisons_11panel_tests.csv`; no duplicate statistics file is produced.

`CellCycle_ecdf_rmse_vs_ploidy.pdf` uses treated mice only for the correlation and one threshold-independent reference formed by averaging the ECDF of each of the eight 0 mg/kg mice with equal sample weight. Initial ploidy, rather than the threshold-defined ETP group, supplies the point shape.

Figures can be regenerated without the cell-level inputs and without rerunning statistical tests by using `--figures_only=TRUE --tables_root=<existing-output-root>`. In this mode the workflow reads only `plot_data/` and the required files in `stats/`, and writes only `Figures/`.

The NonCellCycle input is used only when deriving sample mean end-timepoint ploidy. All figures and reported associations use CellCycle cells.
