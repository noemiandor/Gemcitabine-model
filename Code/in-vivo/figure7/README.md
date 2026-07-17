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

The eight canonical panel-7F tables are a read-only export from the exact
`ETP_reference_balanced_threshold_2_24` analysis used by
`04i_pseudotime_state_pathways_report.html`, with accumulated pseudotime interval
0.30-0.49. Their reviewed SHA-256 values are pinned in `figure7_config.yaml`;
an embedded report raster is not accepted as plotting data.

## Commands

Routine manager execution:

```bash
bash Manager.sh --mode standard --modules in_vivo_figure7 --run-id <run_id>
```

To export the canonical panel-7F reference from a completed 04i result tree and
then generate and publish Figure 7 in one Manager run:

```bash
bash Manager.sh \
  --mode standard \
  --modules in_vivo_figure7 \
  --run-id <run_id> \
  --figure7-state-pathway-results-root /path/to/04i_pseudotime_state_pathways
```

The supplied results root is normalized and recorded in the Manager export
metadata, module-run notes, Figure 7 `run_config.tsv`, and copied panel-7F
provenance. The eight exported TSVs are retained under that Manager run's
`artifacts/figure7_state_pathway_reference/` directory; the tracked frozen
reference is not overwritten.

To explicitly generate and materialize only 7A-7E:

```bash
bash Manager.sh --mode standard --modules in_vivo_figure7 \
  --figure7-panels-ae-only --run-id <run_id>
```

This mode records a five-panel contract and does not read, validate, render, or
materialize panel 7F. Each included panel is written in both PDF and PNG format.

To run the standalone Figure 7 workflow directly on the HPC and generate all
six source panels:

```bash
module load R/4.4.2-gfbf-2024a

cd /share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels_figures

figure7_output_dir="/share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels_figures/Results/in-vivo/figure7/runs/manual_$(date +%Y%m%d_%H%M%S)_figure7"

Rscript Code/in-vivo/figure7/run_figure7.R \
  --mode=standard \
  --panel-set=a-f \
  --config=Code/in-vivo/figure7/figure7_config.yaml \
  --cellcycle-input=Data/in-vivo/figure7/processed/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv \
  --non-cellcycle-input=Data/in-vivo/figure7/processed/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv \
  --saved-state-pathway-dir=Data/in-vivo/figure7/saved_state_pathway/taoli_04i_etp2_24_day17_v1 \
  --state-pathway-results-root=/share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels/Results/04i_pseudotime_state_pathways \
  --output-dir="${figure7_output_dir}"

echo "Figure 7 results: ${figure7_output_dir}"
```

Run this command in the HPC shell rather than at an interactive R prompt. The
timestamp creates a new output directory for every run. The `standard` mode
recomputes panels 7A-7E and renders panel 7F from the pinned canonical 04i
tables without refitting the 04i model.

Standalone rendering from an immutable completed run (does not rerun statistics):

```bash
Rscript Code/in-vivo/figure7/run_figure7.R \
  --mode=render-only \
  --config=Code/in-vivo/figure7/figure7_config.yaml \
  --source-run-dir=Results/in-vivo/figure7/runs/<source>_figure7 \
  --output-dir=/new/empty/output
```

The manager's `panels-only` mode does not invoke this R script; it materializes six existing PDFs from an explicit `--source-run-id`.

The frozen reference can be regenerated from the unchanged completed 04i result
tree with:

```bash
Rscript Code/in-vivo/figure7/export_04i_state_pathway_reference.R \
  --results-root=/path/to/04i_pseudotime_state_pathways
```

On the HPC, run the exporter from the repository root with explicit source and
output paths:

```bash
cd /share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels_figures

Rscript Code/in-vivo/figure7/export_04i_state_pathway_reference.R \
  --results-root=/share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels/Results/04i_pseudotime_state_pathways \
  --report-html=/share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels/Results/04i_pseudotime_state_pathways/report/04i_pseudotime_state_pathways_report.html \
  --output-dir=/share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels_figures/Results/in-vivo/figure7/hpc_export_recheck_20260717/taoli_04i_etp2_24_day17_v1
```

Run this command in the HPC shell rather than from an interactive R prompt. The
final output-directory basename must remain `taoli_04i_etp2_24_day17_v1`, and
the parent directory must be new because the exporter refuses to overwrite an
existing canonical export.

The exporter verifies the report, source-input, interval-config, and source-table
SHA-256 values before reading results, then writes the eight TSVs atomically. It
does not refit the model or query gene sets. The source analysis did not record a
separate MSigDB release identifier, so provenance retains
`gene_set_release=not_recorded_in_04i_manifest` and the recorded `msigdbr`
package version instead.

Full panel-F recomputation remains a separate guarded path: it requires both an
explicit Seurat RDS and a pinned local gene-set artifact, and it never queries
live `msigdbr`. This prevents a recomputation from silently changing panel 7F.

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

The tests parse all module files, reproduce the frozen A-E numerical results,
enforce treated-only outcomes and selected ECDF IDs 1/8/9, validate the tracked
canonical 04i reference and its lineage, exercise the strict panel-F contract
with generated non-scientific fixtures, and verify fail-fast output behavior.
