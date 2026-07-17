# Figure 7 reproducibility module

This module generates the six approved Figure 7 source panels as matched PDF and
300-DPI PNG files. It does not assemble the final A-F manuscript composite,
matching the Figure 1-6 workflow.

## Frozen routine analysis

- Panels 7A-7E are recomputed from the two tracked cell-level tables under `Data/in-vivo/figure7/processed/`.
- TGI is the Day-17 endpoint statistic, recalculated for each treated mouse using the mean Day-17 growth delta of untreated controls matched by initial ploidy.
- Panels 7C-7E contain exactly eight treated mice at 30 or 120 mg/kg. Untreated mice contribute references only.
- Panel 7D uses the reference-balanced ETP threshold 2.24 and equal-mouse untreated ECDF references.
- Panel 7F is rendered from immutable compact 04i tables in `Data/in-vivo/figure7/saved_state_pathway/taoli_04i_etp2_24_day17_v1/`.

The canonical panel-7F tables were not available when this implementation was written. Standard mode therefore fails before creating analysis output until all eight files and their reviewed SHA-256 values replace the `REQUIRED_CANONICAL_SHA256` markers in `figure7_config.yaml`. An embedded report raster is not accepted as plotting data.

## Commands

Routine manager execution:

```bash
bash Manager.sh --mode standard --modules in_vivo_figure7 --run-id <run_id>
```

Until the canonical panel-7F tables are available, explicitly generate and materialize only 7A-7E:

```bash
bash Manager.sh --mode standard --modules in_vivo_figure7 \
  --figure7-panels-ae-only --run-id <run_id>
```

This mode records a five-panel contract and does not read, validate, render, or
materialize panel 7F. Each included panel is written in both PDF and PNG format.

Standalone rendering from an immutable completed run (does not rerun statistics):

```bash
Rscript Code/in-vivo/figure7/run_figure7.R \
  --mode=render-only \
  --config=Code/in-vivo/figure7/figure7_config.yaml \
  --source-run-dir=Results/in-vivo/figure7/runs/<source>_figure7 \
  --output-dir=/new/empty/output
```

The manager's `panels-only` mode does not invoke this R script; it materializes six existing PDFs from an explicit `--source-run-id`.

Full panel-F recomputation requires both an explicit Seurat RDS and a pinned local gene-set artifact. It never queries live `msigdbr`. The path is currently guarded because the canonical compact reference and the mixed human/mouse feature policy have not been approved; it fails rather than silently changing panel 7F.

## Output contract

Each successful run has `figures/`, `tables/`, `metadata/`, and `logs/`. The
default six-panel contract contains these PDF/PNG pairs:

1. `panel_7A_day17_tgi_calculation.{pdf,png}`
2. `panel_7B_cellcycle_selected_ecdf_comparisons.{pdf,png}`
3. `panel_7C_day17_tgi_by_initial_ploidy.{pdf,png}`
4. `panel_7D_day17_tgi_vs_centered_ecdf_shift.{pdf,png}`
5. `panel_7E_day17_tgi_vs_mean_etp.{pdf,png}`
6. `panel_7F_pseudotime_state_pathway_activity.{pdf,png}`

An explicit `--panel-set=a-e`/`--figure7-panels-ae-only` run instead contains
exactly the first five pairs, records `panel_set=a-e`, and excludes all panel-F
inputs and outputs.

Plotting data, exact-permutation tests, the complete compact state-pathway audit chain, frozen-reference comparison, run settings, panel contract, and session information are retained alongside the PDFs.

## Tests

```bash
Rscript Code/in-vivo/figure7/tests/testthat.R
```

The tests parse all module files, reproduce the frozen A-E numerical results, enforce treated-only outcomes and selected ECDF IDs 1/8/9, exercise the strict panel-F contract with generated non-scientific fixtures, and verify fail-fast output behavior.
