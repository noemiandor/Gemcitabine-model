# Audit: 04h Pseudotime/TGI Analysis

Date: 2026-07-05

## Files Found

- Current analysis script: `Code/in-vivo/04h_pseudotime_TGI_analysis.R`
- Manuscript/result text: `dose_pseudotime_conclusion.tex`
- Cell-cycle input CSV: `Data/in-vivo/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv`
- Non-cell-cycle input CSV: `Data/in-vivo/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv`
- Tumor growth workbook used upstream for TGI export: `Data/in-vivo/dt_Gem_VT_20241223_v4.xlsx`
- Existing result directory: `/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04h_pseudotime_TGI_analysis`
- Existing result archive: `/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04h_pseudotime_TGI_analysis.zip`
- External expression/pseudobulk-adjacent files found under `/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results`, including Seurat RDS objects, scVelo/PAGA h5ad outputs, and pseudobulk tables from 04f/04g.

## Required Inputs

The current 04h script requires two cell-level CSV exports. Both contain:

- `sample_id`, `cluster`, `initial_ploidy`, `gemcitabine_dose`, `gemcitabine_dose_mg_per_kg`
- cell-level `pseudotime` and `cell_ploidy`
- sample-level average/median ploidy and pseudotime summaries
- endpoint TGI columns for every available `Day_*`
- `TGI_percent_auc`
- tumor-volume columns needed to reconstruct endpoint and AUC growth summaries

## Existing Outputs

The existing 04h output directory already contains:

- `CellCycle/` and `NonCellCycle/` data, stats, figures, and summaries
- `Comparison/` paired comparison outputs
- AUC/growth-curve plots
- compatibility outputs for the earlier `pseudotimeAssociations.R` style
- `reported_dose_group_tests.csv`
- `reported_tgi_association_stats.csv`
- `tex_ready_summary.md`

The current output is useful, but it should not be overwritten by the revised analysis. The revised analysis should use a new `04h_pseudotime_TGI_analysis_v2` result directory.

## End-to-End Reproducibility

The current script is runnable from the repository if the two exported CSV files exist in `Data/in-vivo`. It uses `ggplot2` and otherwise mostly base R. However, it has a machine-specific default result root:

```r
"/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results"
```

The script can be redirected with `RESULTS_ROOT`, but command-line arguments are not the primary interface. A revised script should accept `--input_root`, `--results_root`, `--seed`, `--n_perm`, `--n_boot`, `--n_downsample`, and `--run_optional_methods`.

## Package Dependencies

The current script uses:

- required: `ggplot2`
- base R statistics: `stats`

The revised analysis can additionally use installed optional packages where available, but must fail gracefully when optional packages or expression inputs are absent.

## Inference Unit

The current ECDF dose-group tests permute sample-level labels, which is appropriate for biological replicate-level inference. The existing TGI association is also computed at the sample level among treated tumors.

The main weakness is the shift-metric reference:

```r
reference_pseudotime <- df$pseudotime[df$sample_id %in% reference_samples]
```

This pools all cells from ploidy-matched control samples. It can overweight controls with more recovered cells and therefore makes the reference ECDF partially cell-count weighted. The revised primary metric should instead build a sample-level ECDF for each control and average those ECDFs with equal sample weight. The pooled-cell reference can remain as a sensitivity analysis only.

## Cell Pooling Risks

- Density plots and raw visual summaries pool cells and should be interpreted descriptively.
- The current shift metric uses pooled control-cell pseudotime for the reference, which can bias per-sample distances when control sample cell counts differ.
- The dose ECDF test itself is sample-aware, because it computes sample ECDFs and permutes sample labels.

## Current Manuscript Claims

The current `dose_pseudotime_conclusion.tex` claims are mostly directionally supported, with wording caveats:

- Supported: CellCycle pseudotime distributions differ between untreated and treated samples.
- Supported if reproduced in v2: the CellCycle treatment effect remains after initial-ploidy-stratified sample-level permutation.
- Supported: 30 mg/kg and 120 mg/kg separately differ from untreated, while 30 vs 120 is not clearly different; therefore the result should be described as treatment-associated, not dose-dependent.
- Supported but should stay cautious: treated-only CellCycle shift magnitude has a positive nominal association with AUC-based TGI.
- Not strong enough for definitive wording: a predictive biomarker claim, global ploidy independence, or a statistically proven CellCycle-specific TGI association.
- Needs formal support: CellCycle vs NonCellCycle specificity should include paired compartment tests and a correlation-difference permutation rather than relying only on one significant and one non-significant test.

## Revised Analysis Target

`Code/in-vivo/04h_pseudotime_TGI_analysis_v2.R` should preserve the original biological framing while adding:

- primary sample-equal control-reference ECDF shift metrics;
- pooled-cell reference sensitivity;
- sample-aware and ploidy-stratified ECDF dose tests;
- exact or exhaustive TGI-label permutation for small treated-sample association tests;
- leave-one-out, bootstrap, cell-count downsampling, ploidy/dose confounding checks;
- formal paired CellCycle vs NonCellCycle comparison;
- dose-specific TGI exploration for 30 vs 120 mg/kg, including AUC TGI, every available endpoint TGI column, within-dose shift-TGI correlations, and an underpowered interaction sensitivity model;
- composition and cluster-specific analyses;
- optional pseudobulk/gene-program/neighborhood/growth-model analyses when inputs are available;
- revised manuscript text and final report.
