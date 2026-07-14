# Essential pseudotime-TGI analysis

This directory was generated from the CellCycle and NonCellCycle 04h cell-level input tables.

When either cell-level input path is omitted, the standalone entrypoint reconstructs the missing table from `scvelo_cell_metrics.csv`, `all_ploidy.tsv`, `sample_info.xlsx`, and `dt_Gem_VT_20241223_v4.xlsx` under `Data/in-vivo`.

The standalone workflow uses base R plus the `ggplot2`, `ggrepel`, `readr`, and `readxl` packages.

## Selected TGI outcome

The CLI default is Day-17 TGI using the mean of initial-ploidy-matched untreated controls. AUC and the other control summaries remain available through explicit CLI options.

- Measure: `TGI_percent_Day_17`.
- Outcome: `day`.
- Matched untreated-control summary: `mean`.
- TGI day: `17`.

TGI is calculated as `100 * (1 - mouse tumor-growth delta / matched-control reference delta)`. Controls are matched by initial ploidy. `--control-summary` chooses whether the matched untreated-control deltas are summarized by their mean, median, or maximum; it does not summarize the treated mice.

CLI options: `--tgi-outcome=auc|day`, `--control-summary=mean|median|max`, and (for the day outcome) `--tgi-day=<available day>`. With no TGI options, these resolve to `--tgi-outcome=day --tgi-day=17 --control-summary=mean`.

The established output filenames retain `AUC` where present for backward compatibility. The selected outcome is recorded in the statistical `tgi_measure` fields, in scatter-plot metadata, and in plot titles and axes.

## Contents

- `Figures/<method>/`: eight CellCycle PDF figures per method.
- `stats/<method>/`: fourteen statistical CSV files per method.
- `plot_data/<method>/`: eight plotting-data CSV files per method, one for each figure.
- `Figures/TGI_calculation/`: two shared PDFs explaining the selected TGI calculation.
- `stats/TGI_calculation/`: the per-mouse calculation components and run metadata.
- `plot_data/TGI_calculation/`: the exact data behind both shared calculation figures.

## Methods

- `initial_ploidy`: original 2N/4N sample grouping.
- `ETP_fixed_threshold_2_25`: sample mean ETP threshold 2.25.
- `ETP_boundary_stress_threshold_2_375`: sample mean ETP threshold 2.375.
- `ETP_reference_balanced_threshold_2_24`: sample mean ETP threshold 2.24.

## Figure data

Each plotting-data CSV contains the exact rows used by its PDF. Scatter-plot tables also contain the Pearson and Spearman statistics, permutation P values, permutation mode, and the annotation text printed in the figure.

`CellCycle_direct_group_ecdf_comparisons_selected_3panel.pdf` retains original grid positions `(1,1)`, `(3,2)`, and `(3,3)` (panels 1, 8, and 9). Its plotting-data CSV is a subset of the full direct-ECDF plotting table, and both figures reuse `CellCycle_direct_group_ecdf_comparisons_11panel_tests.csv`; no duplicate statistics file is produced.

`CellCycle_TGI_group_boxplot.pdf` compares the selected CLI TGI measure between 2N and 4N treated tumors for the initial-ploidy method, or between ETP-lower and ETP-higher treated tumors for an ETP method. The boxes are adjacent independent groups, not one-to-one mouse pairs. The displayed primary P value comes from a group-label permutation stratified by dose, and the reported effect is adjusted for dose.

`TGI_calculation_growth_trajectories.pdf` shows every mouse's baseline-adjusted tumor-growth trajectory, facets mice by initial ploidy, overlays the selected summary of matched untreated controls, and marks Day 17. `TGI_calculation_components_by_treated_mouse.pdf` shows the treated-mouse growth delta and its matched-control reference for every treated mouse, with TGI recalculated from the plotted values.

`CellCycle_ecdf_rmse_vs_ploidy.pdf` uses treated mice only for the correlation and one threshold-independent reference formed by averaging the ECDF of each of the eight 0 mg/kg mice with equal sample weight. Initial ploidy, rather than the threshold-defined ETP group, supplies the point shape.

Figures can be regenerated without the cell-level inputs and without rerunning statistical tests by using `--figures_only=TRUE --tables_root=<existing-output-root>`. In this mode the workflow reads only `plot_data/` and the required files in `stats/`, and writes only `Figures/`.

The NonCellCycle input is used only when deriving sample mean end-timepoint ploidy. All figures and reported associations use CellCycle cells.
