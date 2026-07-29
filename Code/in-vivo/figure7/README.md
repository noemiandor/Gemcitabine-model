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

## Raw-data fallback and intermediate reuse

`full-refit` is the end-to-end fallback. It is intentionally separate from the
lightweight routine mode:

```bash
bash Manager.sh \
  --mode full-refit \
  --modules in_vivo_figure7,si_figures \
  --run-id <run_id>
```

With those two modules selected, Manager does not run Figures 1-6. It first
validates the downstream frozen/generated caches. A complete valid cache avoids
all raw-data access. When a missing downstream stage needs the final Seurat
object, there are two supported source boundaries:

- pass `--figure7-cellranger-root /path/to/cellranger` to reconstruct the final
  object from the 18 `filtered_feature_bc_matrix.h5` inputs and reuse the five
  fingerprinted Seurat stages under the shared
  `--figure7-seurat-upstream-dir`; or
- omit that option to reuse or download the deposited final Seurat RDS pinned
  by `zenodo_required_files.tsv`.

Figure 7A-7E additionally reuse or download the 18 deposited loom files. The
complete Zenodo fallback is about 10.61 GiB. Manager then runs only the missing
figure-facing stages:

1. reconstruct or validate the shared final Seurat object when required;
2. calculate scVelo pseudotime and derive the CellCycle and NonCellCycle
   Day-17 TGI tables used by 7A-7E;
3. fit the state-pathway model and export a compact generated reference for 7F;
4. build the 11 plot-facing Supplementary Figure tables and render the four
   composites.

The upstream Seurat reconstruction is the exact narrow sequence needed from
Tao's `01_data.R`, `01a_cell_cycle.R`, `02b_cluster_refine.R`,
`02d_manual_cluster_merge.R`, and `03_final_cluster.R`: sample QC/integration,
cell-cycle annotation, reviewed UMAP cluster refinement, fixed manual merges,
and final cluster filtering/reduction. Marker surveys, exploratory plots, and
other downstream analyses are not ported. `--jobs` is recorded and forwarded,
but these reviewed Seurat stages always use one scientific worker.

No FASTQ-to-Cell-Ranger invocation or FASTQ collection was available in the
source work, so the 18 H5 matrices are the earliest executable expression-data
boundary. Likewise, no complete Numbat/karyotyping workflow producing
`all_ploidy.tsv` was available; that checksum-pinned table is therefore a
versioned source artifact rather than a generated cache. The deposited final
RDS remains the verified fallback when the external H5 boundary is unavailable.

The endpoint-ploidy table, sample workbook, and growth-curve workbook are
versioned source artifacts. Their revision and SHA-256 values are pinned in
`figure7_config.yaml`. Raw downloads and generated intermediates stay below
`Results/`; they are not publication inputs and are not committed.

Every reusable stage has a dependency and output fingerprint. A subsequent run
reuses a stage only when its source hashes, relevant configuration, code,
runtime package contract, and output hashes validate. Stale automatically
managed caches are preserved with a `.stale.<timestamp>` suffix before
regeneration; explicit external paths are never modified. To require an
already-populated raw cache, add `--figure7-no-download-missing-raw`.

The raw-generated 7F reference deliberately remains separate from the reviewed
canonical v1 reference. Until the deferred GRCh/GRCm correction is implemented,
the fallback preserves Tao's mixed-species behavior and records
`canonical_publication_allowed=false`. Routine `standard` runs continue to use
the reviewed frozen 7F input.

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

Routine and raw-fallback runs do not refresh the tracked reviewed v1 panel-7F
reference. Full-workflow instead writes a separately identified
`runtime_state_pathway_legacy_mixed_v1` generated reference below the run
intermediates and marks it noncanonical.

To strictly re-export the reviewed canonical panel-7F reference from the
completed 04i result tree, then generate and publish Figure 7 in one Manager
run:

```bash
bash Manager.sh \
  --mode standard \
  --modules in_vivo_figure7 \
  --run-id <run_id> \
  --figure7-state-pathway-results-root /path/to/04i_pseudotime_state_pathways
```

This explicit option runs
`Code/in-vivo/figure7/export_04i_state_pathway_reference.R`. The exporter
requires the reviewed report hash, source revision, source input/config
checksums, and exact hashes of all eight scientific source tables. It records
the normalized runtime source root only in run metadata; canonical provenance
uses stable identifiers and checksums and therefore remains portable.

The verified export is retained under the Manager run's
`artifacts/figure7_state_pathway_reference/` directory. Only after Figure 7 and
its manifests finish successfully does Manager atomically materialize the eight
TSVs under
`Data/in-vivo/figure7/saved_state_pathway/taoli_04i_etp2_24_day17_v1/` and
write `metadata/figure7_state_pathway_materialization.tsv`. A failed Figure 7
run cannot refresh the tracked reference.

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
  --saved-state-pathway-dir=Data/in-vivo/figure7/saved_state_pathway/taoli_04i_etp2_24_day17_v1 \
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

The older `full-analysis` entrypoint remains available for comparison against a
separately supplied gene-set artifact. Manager's `full-refit` path instead uses
the raw-data fallback above, pins `msigdbr` and the MSigDB release, and keeps its
generated legacy-mixed reference noncanonical.

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
canonical 04i reference and its lineage, exercise the strict panel-F contract,
verify cached-stage fingerprints and tamper rejection, validate the complete
Zenodo manifest, and confirm that missing raw inputs fail before output is
created when downloading is disabled.
