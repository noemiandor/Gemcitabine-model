# TGI Calculation Process

This document describes how TGI was calculated for the exported data used with:

`fig4_pseudotime_distribution_per_sample.pdf`

Output CSV:

`fig4_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv`

## Input Data

The tumor growth curve data were read from:

`/Users/4482173/Documents/GitHub/Gemcitabine-model/Data/in-vivo/dt_Gem_VT_20241223_v4.xlsx`

The analysis used `Sheet1`.

The key columns are:

- `Sequencing IDs`: sample identifier used to match the scRNA-seq sample.
- `harvest`: source label containing initial ploidy and Gemcitabine dose.
- `Day_0`: baseline tumor volume.
- `Day_38`: final tumor volume in this dataset.

Rows without a valid `Sequencing IDs` value or without finite tumor-volume measurements were excluded from TGI calculation.

## Sample ID Matching

The scVelo/04f sample IDs do not include the trailing `-HM` suffix, while some growth-curve `Sequencing IDs` do. Therefore, sample IDs were normalized before joining:

- `2N-A1-0-HM` becomes `2N-A1-0`
- `2N-A1-LR-HM` becomes `2N-A1-LR`
- `4N-A5-0` remains `4N-A5-0`
- `A5-4N-L` remains `A5-4N-L`

This normalized ID was then joined to the 04f `sampleID`.

## Initial Ploidy and Dose Assignment

Initial ploidy was inferred from `Sequencing IDs` and `harvest`:

- Samples containing `2N` were assigned `2N`.
- Samples containing `4N`, including IDs such as `A5-4N-L`, were assigned `4N`.

Gemcitabine dose was inferred from the `harvest` label:

`SUM159-<initial ploidy>-<dose>-<sample position>_harvest`

For example:

- `SUM159-2N-0-O_harvest` -> `0mg/kg`
- `SUM159-2N-30-L_harvest` -> `30mg/kg`
- `SUM159-4N-120-RR_harvest` -> `120mg/kg`

## Tumor Growth Delta

For each sample:

```text
delta = final tumor volume - baseline tumor volume
```

In this dataset:

```text
delta = Day_38 - Day_0
```

The helper code supports using the last finite `Day_*` column as the final volume when `final_day` is not explicitly supplied. For the current file, the last finite day for all exported samples is `Day_38`.

## Control Matching

Controls were matched separately by initial ploidy:

- 2N treated samples were compared against the mean delta of 2N `0mg/kg` controls.
- 4N treated samples were compared against the mean delta of 4N `0mg/kg` controls.

The matched control means used in this run were:

| initial ploidy | control dose | control n | mean control delta |
|---|---:|---:|---:|
| 2N | 0mg/kg | 4 | 2107.105 |
| 4N | 0mg/kg | 4 | 3315.015 |

## TGI Formula

TGI was calculated as:

```text
TGI_percent = 100 * (1 - delta_treated / mean_delta_control)
```

where `mean_delta_control` is the ploidy-matched `0mg/kg` control mean.

This formula was applied to every sample, including the control samples themselves. Therefore, individual control samples can have positive or negative TGI depending on whether their own delta is below or above the matched control-group mean.

## Exported CSV Structure

The exported CSV is cell-level because it is intended to reproduce the per-sample pseudotime distribution in fig4. Each row is one target Tumor cell from clusters `4c`, `6`, or `10`.

Key columns:

- `cell_id`: cell barcode/ID.
- `sample_id`: normalized sample ID.
- `cluster`: target cluster (`4c`, `6`, or `10`).
- `initial_ploidy`: sample initial ploidy (`2N` or `4N`).
- `gemcitabine_dose`: Gemcitabine dose label.
- `gemcitabine_dose_mg_per_kg`: numeric dose.
- `pseudotime`: scVelo pseudotime used in fig4.
- `cell_ploidy`: cell-level ploidy value.
- `average_ploidy`: sample-level mean `cell_ploidy` across the exported target cells.
- `tumor_volume_baseline_day`: baseline day used for delta.
- `tumor_volume_baseline`: baseline tumor volume.
- `tumor_volume_final_day`: final day used for delta.
- `tumor_volume_final`: final tumor volume.
- `tumor_volume_delta`: final minus baseline tumor volume.
- `matched_control_mean_delta`: initial-ploidy-matched control mean delta.
- `matched_control_n`: number of control samples in the matched control group.
- `TGI_percent`: calculated tumor growth inhibition.

## Code Location

The reusable helper function is:

`calculate_in_vivo_tgi()` in `/Users/4482173/Documents/GitHub/Gemcitabine-model/Code/in-vivo/Utils.R`

The 04f export logic is in:

`/Users/4482173/Documents/GitHub/Gemcitabine-model/Code/in-vivo/04f_cell_cycle_results.R`
