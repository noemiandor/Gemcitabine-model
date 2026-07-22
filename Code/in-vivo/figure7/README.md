# Figure 7 reproducibility module

This module generates the six approved Figure 7 source panels as matched PDF and
300-DPI PNG files. It does not assemble the final A-F manuscript composite,
matching the Figure 1-6 workflow.

## Frozen routine analysis

- Panels 7A-7E are recomputed from the two tracked cell-level tables under `Data/in-vivo/figure7/processed/`.
- TGI is the Day-17 endpoint statistic, recalculated for each treated mouse using the mean Day-17 growth delta of untreated controls matched by initial ploidy.
- Panels 7C-7E contain exactly eight treated mice at 30 or 120 mg/kg. Untreated mice contribute references only.
- Panel 7D uses the reference-balanced ETP threshold 2.24 and equal-mouse untreated ECDF references.
- Panel 7F is rendered from immutable compact state-pathway tables in `Data/in-vivo/figure7/saved_state_pathway/taoli_state_pathway_etp2_24_day17_v1/`.

The eight canonical panel-7F tables are a read-only export from the exact
`ETP_reference_balanced_threshold_2_24` analysis used by
`pseudotime_state_pathways_report.html`, with accumulated pseudotime interval
0.30-0.49. Their reviewed SHA-256 values are pinned in `figure7_config.yaml`;
an embedded report raster is not accepted as plotting data.

## Commands

Routine manager execution:

```bash
bash Manager.sh --mode standard --run-id <run_id>
```

Figure 7 is part of the default manuscript module set. To run only Figure 7,
add `--modules in_vivo_figure7`.

The canonical manuscript run remains Day-17 TGI materialized to
`figures/Figure7/`. To generate a Day-24 supplementary variant without
overwriting the canonical assets:

```bash
bash Manager.sh \
  --mode standard \
  --modules in_vivo_figure7 \
  --run-id <run_id> \
  --figure7-tgi-day 24 \
  --figure7-figure-name Figure7_Supplement
```

Likewise, a Day-31 variant can be kept separately:

```bash
bash Manager.sh \
  --mode standard \
  --modules in_vivo_figure7 \
  --run-id <run_id> \
  --figure7-tgi-day 31 \
  --figure7-figure-name Figure7_Supplement2
```

The selected endpoint day is used consistently for the matched-control TGI
calculation, plot labels, statistical tables, run metadata, panel contract, and
day-bearing filenames. Panel 7B and panel 7F are scientifically independent of
the TGI endpoint and are regenerated unchanged into the selected destination.

To export the canonical panel-7F reference from a completed state-pathway result tree and
then generate and publish Figure 7 in one Manager run:

```bash
bash Manager.sh \
  --mode standard \
  --modules in_vivo_figure7 \
  --run-id <run_id> \
  --figure7-state-pathway-results-root /path/to/pseudotime_state_pathways
```

The supplied results root is normalized and recorded in the Manager export
metadata, module-run notes, and Figure 7 `run_config.tsv`. Canonical panel-7F
provenance is deliberately location-independent: it records stable report/source
identifiers and checksums, never the runtime filesystem location. The eight
exported TSVs are retained under that Manager run's
`artifacts/figure7_state_pathway_reference/` directory. After Figure 7 and its
manifests complete successfully, Manager refreshes the same eight TSVs using
atomic per-file replacement under
`Data/in-vivo/figure7/saved_state_pathway/taoli_state_pathway_etp2_24_day17_v1/` and
records that publication in `metadata/figure7_state_pathway_materialization.tsv`.
Failed Figure 7 runs do not refresh the tracked canonical Data reference.

To explicitly generate and materialize only 7A-7E:

```bash
bash Manager.sh --mode standard --modules in_vivo_figure7 \
  --figure7-panels-ae-only --run-id <run_id>
```

This mode records a five-panel contract and does not read, validate, render, or
materialize panel 7F. Each included panel is written in both PDF and PNG format.

## End-to-end workflow from loom and Seurat

`run_figure7.R --mode=full-workflow` can start from raw loom files, resume from
`scvelo_cell_metrics.csv`, or reuse an existing complete CellCycle/NonCellCycle
cell-table pair. The preflight order is:

1. If both cell-level tables exist, skip scVelo and cell-level generation.
2. If neither cell-level table exists but scVelo metrics exist, generate only
   the two cell-level tables.
3. If none of those intermediates exist, generate scVelo metrics from loom and
   Seurat first, then generate the two cell-level tables.
4. If only one cell-level table exists, fail rather than mix partial outputs.

If scVelo or panel 7F requires raw inputs and no explicit local paths were
provided, the workflow uses the open dataset at DOI
[`10.5281/zenodo.21463392`](https://doi.org/10.5281/zenodo.21463392). The pinned
contract contains 18 loom files plus the final Seurat RDS (11,395,115,098 bytes,
approximately 10.61 GiB). Files are cached under
`Data/in-vivo/figure7/raw/zenodo_21463392/`, downloaded through resumable `.part`
files, and promoted only after the Zenodo size/MD5 checks pass. The Seurat RDS
also must match the frozen SHA-256 in `figure7_config.yaml`. A valid cache is
verified and reused without network transfer.

When `aria2c` is available, the downloader runs multiple files concurrently and
uses multiple HTTP range connections per file. The defaults are four concurrent
files and two connections per file. Configure them with `--download-workers`
and `--download-connections-per-file` (each 1-16). The fallback order is
external wget, external curl, then R libcurl; resumable `.part` files are retained
after an interrupted download.

After the cell-level pair is ready, the workflow generates the state-pathway
result tree, exports the compact panel-7F reference, and renders Figure 7A-F.
For files generated inside `full-workflow`, validation uses the complete table
schema, frozen scientific parameters, row/grid counts, pathway-selection audit,
and recorded observed SHA-256 values. Byte-for-byte hashes from the reviewed
snapshot remain mandatory in `standard` mode, but are not used to reject
numerically equivalent newly serialized CSV/TSV files.
The run metadata records the actual scVelo (when present), CellCycle, and
NonCellCycle SHA-256 values used by that invocation. When the scVelo stage reads
the Seurat RDS, it also exports every `obj@meta.data` column, with cell barcodes
in the first `cell` column and `UMAP_1`/`UMAP_2` from the Seurat `umap`
reduction, to `Data/in-vivo/seurat_metadata.csv` by default.
Use `--seurat-metadata-output` to select another persistent path.

Automatic download from the pinned Zenodo record:

```bash
Rscript Code/in-vivo/figure7/run_figure7.R \
  --mode=full-workflow \
  --python=/path/to/scvelo/python \
  --cell-ploidy-input=/path/to/all_ploidy.tsv \
  --sample-info-input=/path/to/sample_info.xlsx \
  --growth-curve-input=/path/to/dt_Gem_VT_20241223_v4.xlsx \
  --download-workers=4 \
  --download-connections-per-file=2 \
  --output-dir=/path/to/new_figure7_run
```

Use `--raw-data-dir=/path/to/cache` to move the cache, or
`--download-missing-raw=false` to prohibit network downloads and fail if the
cache is missing or corrupt.

To use existing raw files explicitly and skip Zenodo materialization:

```bash
Rscript Code/in-vivo/figure7/run_figure7.R \
  --mode=full-workflow \
  --loom-root=/path/to/velocyto_loom \
  --seurat-rds=/path/to/integrated_sct_cca_seurat_final_reclustered.rds \
  --seurat-metadata-output=Data/in-vivo/seurat_metadata.csv \
  --python=/path/to/scvelo/python \
  --cell-ploidy-input=/path/to/all_ploidy.tsv \
  --sample-info-input=/path/to/sample_info.xlsx \
  --growth-curve-input=/path/to/dt_Gem_VT_20241223_v4.xlsx \
  --intermediate-dir=/path/to/figure7_intermediates \
  --config=Code/in-vivo/figure7/figure7_config.yaml \
  --output-dir=/path/to/new_figure7_run
```

To inspect which stage would run without creating outputs:

```bash
Rscript Code/in-vivo/figure7/run_figure7.R \
  --mode=full-workflow \
  --preflight-only=true \
  --loom-root=/path/to/velocyto_loom \
  --seurat-rds=/path/to/seurat.rds \
  --cell-ploidy-input=/path/to/all_ploidy.tsv \
  --sample-info-input=/path/to/sample_info.xlsx \
  --growth-curve-input=/path/to/growth_curve.xlsx \
  --output-dir=/path/to/planned_figure7_run
```

Intermediate files are reused by default. `--overwrite-intermediates=true` is
required to replace a partial pair or an existing state-pathway run directory.
The final Figure 7 output directory must still be new and empty.

To run the standalone Figure 7 workflow directly on the HPC and generate all
six source panels:

```bash
module load Python/3.12.3-GCCcore-13.3.0
module load R/4.4.2-gfbf-2024a

cd /share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels_figures

figure7_output_dir="/share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels_figures/Results/in-vivo/figure7/runs/manual_$(date +%Y%m%d_%H%M%S)_figure7"

Rscript Code/in-vivo/figure7/run_figure7.R \
  --mode=standard \
  --panel-set=a-f \
  --config=Code/in-vivo/figure7/figure7_config.yaml \
  --cellcycle-input=Data/in-vivo/figure7/processed/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv \
  --non-cellcycle-input=Data/in-vivo/figure7/processed/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv \
  --saved-state-pathway-dir=Data/in-vivo/figure7/saved_state_pathway/taoli_state_pathway_etp2_24_day17_v1 \
  --state-pathway-results-root=/share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels/Results/pseudotime_state_pathways \
  --output-dir="${figure7_output_dir}"

echo "Figure 7 results: ${figure7_output_dir}"
```

Run this command in the HPC shell rather than at an interactive R prompt. The
timestamp creates a new output directory for every run. The `standard` mode
recomputes panels 7A-7E and renders panel 7F from the pinned canonical state-pathway
tables without refitting the state-pathway model.

Standalone rendering from an immutable completed run (does not rerun statistics):

```bash
Rscript Code/in-vivo/figure7/run_figure7.R \
  --mode=render-only \
  --config=Code/in-vivo/figure7/figure7_config.yaml \
  --source-run-dir=Results/in-vivo/figure7/runs/<source>_figure7 \
  --output-dir=/new/empty/output
```

The manager's `panels-only` mode does not invoke this R script; it materializes six existing PDFs from an explicit `--source-run-id`.

The frozen reference can be regenerated from the unchanged completed state-pathway result
tree with:

```bash
Rscript Code/in-vivo/figure7/export_state_pathway_reference.R \
  --results-root=/path/to/pseudotime_state_pathways \
  --output-dir=Results/in-vivo/figure7/reference_exports/<run_id>/taoli_state_pathway_etp2_24_day17_v1
```

The exporter audits the source CSV/manifest tables from the state-pathway result tree. It
does not read or checksum the rendered HTML report.

On the HPC, run the exporter from the repository root with explicit source and
output paths:

```bash
cd /share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels_figures

Rscript Code/in-vivo/figure7/export_state_pathway_reference.R \
  --results-root=/share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels/Results/pseudotime_state_pathways \
  --output-dir=/share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels_figures/Results/in-vivo/figure7/hpc_export_recheck_20260717/taoli_state_pathway_etp2_24_day17_v1
```

Run this command in the HPC shell rather than from an interactive R prompt. The
final output-directory basename must remain `taoli_state_pathway_etp2_24_day17_v1`, and
the parent directory must be new because the exporter refuses to overwrite an
existing canonical export.

By default, the exporter verifies the source input checksum records and the
exact SHA-256 values of all eight consumed scientific source tables before
parsing those tables, then writes the eight TSVs atomically. The report hash is
recorded provenance and the rendered HTML itself is not read in strict export
mode. The full workflow sets `--verify-source-checksums=FALSE` for its newly
generated result tree, applies the same structural and scientific audits, and
records the observed source-table hashes plus the support-script SHA-256. If
that support workflow did not create an HTML report, provenance records
`not_generated_by_support_workflow` rather than attributing the reviewed report
hash to the new run. The export step itself does not refit the model or query
gene sets. Runtime source paths remain in Manager/run
metadata and are excluded from the immutable canonical provenance. The source
analysis did not record a separate MSigDB release identifier, so provenance retains
`gene_set_release=not_recorded_in_source_manifest` and the recorded `msigdbr`
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
canonical state-pathway reference and its lineage, exercise the strict panel-F contract
with generated non-scientific fixtures, and verify fail-fast output behavior.
