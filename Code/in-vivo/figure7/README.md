# Figure 7 reproducibility module

This module generates six reviewed scientific source panels as matched PDF and
300-DPI PNG files and assembles the manuscript-facing A-K composite
`Figure7_reviewed_GRCh.png`. The composite reuses Supplementary Figure 4A-C/E
and Supplementary Figure 7B through the shared production plotting helper;
their supplementary copies are retained. Routine A-F source runs are
publication eligible and use the byte-pinned human-only panel-7F v2 reference.
The historical mixed-feature v1 reference remains available only for audit,
while newly generated raw-refit references remain noncanonical until
separately reviewed.

## Frozen routine analysis

- Panels 7A-7E are recomputed from the two tracked cell-level tables under `Data/in-vivo/figure7/processed/`.
- TGI is the Day-17 endpoint statistic, recalculated for each treated mouse using the mean Day-17 growth delta of untreated controls matched by initial ploidy.
- Panels 7C-7E contain exactly eight treated mice at 30 or 120 mg/kg. Untreated mice contribute references only.
- Panel 7D uses the reference-balanced ETP threshold 2.24 and equal-mouse untreated ECDF references.
- Panel 7F is rendered from immutable compact tables in `Data/in-vivo/figure7/saved_state_pathway/state_pathway_grch_human_only_etp2_24_day17_v2/`.
- The A-K manuscript composite binds the exact reviewed 11-table SI cache and
  displays, in first-citation order: source 7A, source 7C, SI4A-C, SI4E, SI7B,
  source 7B, source 7F, and source 7D-E.

The eight reviewed panel-7F files retain exact `GRCh38-` features before
expression filtering, symbol resolution, model fitting, and Homo sapiens GSEA.
They use the `ETP_reference_balanced_threshold_2_24` analysis and accumulated
pseudotime interval 0.30-0.49. The displayed selection contains 21
FDR-significant pathways (Hallmark 5, Reactome 8, GO biological process 8):
collection-wide BH-adjusted P <= 0.05, up to four per sign and collection, with
no nonsignificant backfill. All eight SHA-256 values and the exact approved
retry7 lineage are pinned in `figure7_config.yaml` and the reviewed provenance.
An embedded report raster is not accepted as plotting data.

## Raw-data fallback and intermediate reuse

`full-refit` is the end-to-end fallback. It is intentionally separate from the
lightweight routine mode:

```bash
bash Manager.sh \
  --mode full-refit \
  --modules in_vivo_figure7,si_figures \
  --run-id <run_id>
```

With those two modules selected, Manager does not run Figures 1-6. When the
full A-K figure is requested, Manager automatically schedules `si_figures`
before `in_vivo_figure7`, even if only the latter was listed. It first validates
corrected generated caches; the historical panel-7F reference and reviewed SI
plot-only cache do not satisfy an explicit full refit. A complete,
lineage-valid generated cache avoids raw-data access. When a missing downstream
stage needs the final Seurat object, there are two supported source boundaries:

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
2. build the run-scoped, human-only 11-table supplementary cache and render
   Supplementary Figures 4-7;
3. calculate scVelo pseudotime and derive the CellCycle and NonCellCycle
   Day-17 TGI tables used by 7A-7E;
4. fit the state-pathway model and export a compact generated reference for 7F;
5. pass the generated supplementary cache and its complete raw-input lineage
   into Figure 7 and assemble the A-K review candidate.

The generated-cache handoff is required in `full-refit`: Figure 7 will not
silently fall back to the reviewed supplementary cache. The resulting
`Figure7_generated_GRCh_candidate.png` is explicitly noncanonical. Only the
exact reviewed cache may produce `Figure7_reviewed_GRCh.png` in routine mode.

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

The byte-pinned v1 panel-7F reference is retained only for historical audit
because it was fitted from mixed human/mouse features. It can still be rendered
for comparison by explicitly supplying its directory, but that A-F run is
never publication eligible. Routine `standard` uses the reviewed human-only v2
reference. The corrected full-workflow path retains exact `GRCh38-` count rows
before expression filtering, symbol resolution, modeling, and GSEA. It writes
a separate generated human-only reference with
`canonical_publication_allowed=false`; raw reruns do not inherit the approval
of the exact frozen v2 bytes. GSEA starts with the configured simple-permutation budget,
retries only unresolved pathways at increasing pinned budgets, recomputes BH
adjustment across each complete collection, and fails closed if any pathway
still lacks finite statistics at the configured cap. The generated human-only
panel displays only pathways with collection-wide BH-adjusted P <= 0.05, then
takes up to four pathways in each direction and collection. It never backfills
a direction with nonsignificant pathways, so generated collection and panel row
counts may be smaller than the historical frozen 8/24-pathway layout. A-E runs
remain publication eligible.

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

Routine and raw-fallback runs do not refresh either tracked frozen reference.
Full-workflow instead writes a separately identified
`runtime_state_pathway_grch_human_only_v2` generated reference below the run
intermediates and marks it noncanonical.

To strictly re-export the historical mixed-feature panel-7F reference from the
completed 04i result tree for audit and comparison:

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
the normalized runtime source root only in run metadata; the portable artifact
provenance uses stable identifiers and checksums.

The verified export is retained only under the Manager run's
`artifacts/figure7_state_pathway_reference/` directory. Manager labels the A-F
result noncanonical, does not publish it to `figures/Figure7/`, and never copies
the eight TSVs into tracked `Data/`. This option preserves the original v1
audit chain; it is not a substitute for a reviewed human-only v2 reference.

To explicitly generate and materialize only 7A-7E:

```bash
bash Manager.sh --mode standard --modules in_vivo_figure7 \
  --figure7-panels-ae-only --run-id <run_id>
```

This mode records a five-panel contract and does not read, validate, render, or
materialize panel 7F. Each included panel is written in both PDF and PNG format.

To run the standalone canonical Figure 7 workflow directly on the HPC and
render all six panels:

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
  --saved-state-pathway-dir=Data/in-vivo/figure7/saved_state_pathway/state_pathway_grch_human_only_etp2_24_day17_v2 \
  --si-table-cache-dir=Data/in-vivo/SIfigures \
  --output-dir="${figure7_output_dir}"

echo "Figure 7 results: ${figure7_output_dir}"
```

Run this command in the HPC shell rather than at an interactive R prompt. The
timestamp creates a new output directory for every run. The `standard` mode
recomputes source panels 7A-7E, renders source panel 7F from the pinned reviewed
human-only tables without refitting the model, and assembles the A-K main
composite from those sources and the reviewed SI cache; the resulting run is
canonical.

Standalone rendering from an immutable completed run (does not rerun statistics):

```bash
Rscript Code/in-vivo/figure7/run_figure7.R \
  --mode=render-only \
  --config=Code/in-vivo/figure7/figure7_config.yaml \
  --source-run-dir=Results/in-vivo/figure7/runs/<source>_figure7 \
  --output-dir=/new/empty/output
```

The manager's `panels-only` mode does not invoke this R script. For Figure 7 it
materializes the exact recorded A-K composite/source-panel or explicit A-E
contract from an explicit `--source-run-id`. Historical and generated/
noncanonical runs are rejected.

The older artifact-driven `full-analysis` entrypoint remains only as an
explicit guard that directs callers to the corrected workflow. Manager's
`full-refit` path uses the raw-data fallback above, pins `msigdbr` and the
MSigDB release, applies the exact human-only feature policy, and keeps its
generated reference noncanonical until scientific review.

## Output contract

Each successful run has `figures/`, `tables/`, `metadata/`, and `logs/`. The
default six-panel source contract contains these PDF/PNG pairs:

1. `panel_7A_day17_tgi_calculation.{pdf,png}`
2. `panel_7B_cellcycle_selected_ecdf_comparisons.{pdf,png}`
3. `panel_7C_day17_tgi_by_initial_ploidy.{pdf,png}`
4. `panel_7D_day17_tgi_vs_centered_ecdf_shift.{pdf,png}`
5. `panel_7E_day17_tgi_vs_mean_etp.{pdf,png}`
6. `panel_7F_pseudotime_state_pathway_activity.{pdf,png}`

It additionally contains `Figure7_reviewed_GRCh.png`, whose panel contract is:
A=7A, B=7C, C=SI4A, D=SI4B, E=SI4C, F=SI4E, G=SI7B, H=7B, I=7F, J=7D,
and K=7E. The reviewed SI manifest hash and this ordered mapping are recorded in
run metadata and enforced during materialization.

A full-refit run has the same scientific panel mapping but writes
`Figure7_generated_GRCh_candidate.png`. Its generated supplementary-cache
manifest, analysis-input manifest, run configuration, and provenance are all
hash-bound in Figure 7 metadata; `canonical_publication_allowed=false` prevents
the candidate from being materialized as a manuscript asset.

An explicit `--panel-set=a-e`/`--figure7-panels-ae-only` run instead contains
exactly the first five pairs, records `panel_set=a-e`, and excludes all panel-F
inputs and outputs. It also excludes the A-K composite because the main figure
requires the reviewed state-pathway and SI-cache contracts together.

Plotting data, exact-permutation tests, the complete compact state-pathway audit chain, frozen-reference comparison, run settings, panel contract, and session information are retained alongside the PDFs.

## Tests

```bash
Rscript Code/in-vivo/figure7/tests/testthat.R
```

The tests parse all module files, reproduce the frozen A-E numerical results,
enforce treated-only outcomes and selected ECDF IDs 1/8/9, validate the tracked
historical 04i reference and its lineage, exercise the strict panel-F
reviewed-v2 publication guard and generated human-only contract, verify cached-stage
fingerprints and tamper rejection, validate the complete Zenodo manifest, and
confirm that missing raw inputs fail before output is created when downloading
is disabled.
