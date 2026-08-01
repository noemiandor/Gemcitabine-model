# Figure 7 reproducibility module

This module generates six scientific source panels as matched PDF and 300-DPI
PNG files and can assemble the manuscript-facing A-K composite
`Figure7_reviewed_GRCh.png`. The composite reuses Supplementary Figure 4A-C/E
and Supplementary Figure 7B through the shared production plotting helper;
their supplementary copies are retained. The former human-only panel-7F v2
reference is retained for audit but is superseded because its model adjusted
for a run-confounded endpoint-CN-score group. The reviewed v3 reference instead
adjusts for injected initial ploidy and is the canonical panel-7F source. Panel
K shows the unadjusted mouse-level Pearson association between the
checksum-pinned mean endpoint tumor-cell ploidy and TGI, with unrestricted
exact enumeration of all 8! TGI-label permutations; endpoint ploidy is not
used as a nuisance covariate in the state-pathway model.

## Frozen routine analysis

- Panels 7A-7D are recomputed from the two tracked plot-facing cell-level
  tables under `Data/in-vivo/figure7/processed/`. Source panel 7E (main panel
  K) additionally reads the exact six-column, checksum-pinned
  `Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv`: 14,125 cells across 16 CBS
  files. This remains the complete immutable CBS inventory. Every scored cell
  must also occur, with the identical file, barcode, and value, in the exact
  9,832-cell QC-passed CellCycle + NonCellCycle union represented in the
  plot-facing tables. The other 4,293 inventory cells remain provenance-only
  and are excluded from scoring. Panel K uses the resulting 5,335 treated
  cells.
  That newer artifact has its own source revision (`dcdb62f2252...`) in the
  config; the surrounding legacy source revision does not claim it existed in
  the earlier snapshot.
- TGI is the Day-17 endpoint statistic, recalculated for each treated mouse using the mean Day-17 growth delta of untreated controls matched by initial ploidy.
- Panels 7C-7E contain exactly eight treated mice at 30 or 120 mg/kg. Untreated mice contribute references only.
- Panel 7D (main panel J) uses equal-mouse untreated ECDF references matched by
  injected initial ploidy. The treated-mouse association is Pearson r =
  0.7399455 with exact within-dose permutation P = 0.0173611 (576 labelings).
- Panel 7F (main panel I) uses model
  `initial_ploidy_adjusted_grch_human_only_v3`; its nuisance terms are dose and
  injected initial ploidy, never endpoint CN score or an endpoint-derived
  threshold group.
- Panel 7E (main panel K) reports the descriptive mouse-level association
  between mean endpoint tumor-cell ploidy and Day-17 TGI (Pearson r =
  -0.6984010; asymptotic P = 0.0540070; exact unrestricted permutation P =
  0.0591270; 40,320 assignments). It is an unadjusted association, not
  evidence for an independent or causal terminal-ploidy effect.
- The A-K manuscript composite binds the exact reviewed 11-table SI cache and
  displays, in first-citation order: source 7A, source 7C, SI4A-C, SI4E, SI7B,
  source 7B, source 7F, and source 7D-E.

The eight reviewed v3 panel-7F files retain exact `GRCh38-` features before
expression filtering, symbol resolution, model fitting, and Homo sapiens GSEA.
They use the accumulated pseudotime interval 0.30-0.49. Displayed pathways must
have collection-wide BH-adjusted P <= 0.05; the selector takes up to four per
sign and collection and never backfills with nonsignificant pathways. The
reviewed provenance authenticates all eight compact tables, the exact input
hashes, model/design audit, feature-species boundary, MSigDB release, and
selection rule. An embedded report raster is not accepted as plotting data.

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

Figure 7A-7D additionally reuse or download the 18 deposited loom files. Source
panel 7E/main panel K binds the complete combined CBS inventory; in
`full-refit`, Manager regenerates that run-scoped table from the manifest-pinned
16 CBS matrices and requires it to reproduce the canonical checksum. It then
restricts scoring to exact file+barcode keys retained in the final Seurat tumor
universe and represented by the two processed tables. The
complete Zenodo fallback is about 10.61 GiB. Manager then runs only the missing
figure-facing stages:

1. reconstruct or validate the shared final Seurat object when required;
2. build the run-scoped, human-only 11-table supplementary cache and render
   Supplementary Figures 4-7;
3. calculate scVelo pseudotime and derive the CellCycle and NonCellCycle
   Day-17 TGI tables used by 7A-7D and for sample/TGI metadata in 7E;
4. regenerate and checksum-validate the complete six-column, 14,125-cell CBS
   inventory, then calculate 7E's per-mouse scores from the exact 9,832-cell
   QC-passed union (5,335 cells across the eight treated tumors);
5. fit the state-pathway model and export a compact generated reference for 7F;
6. pass the generated supplementary cache and its complete raw-input lineage
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
boundary. The 16 tracked downstream NUMBAT-derived CBS matrices reproduce the
SI6E chromosome-state view and recompute every value in the canonical
`scRNAseq_Numbat/all_ploidy.csv` (and its reduced `all_ploidy.tsv` projection).
Two checksum-pinned, project-designated lineage-matched karyotype proxies
provide the 2N-A7M and 4N-A5M reference distributions used by SI6F; they are
not the same-passage A6M/A4M inocula. The reference calculation multiplies the
autosomal length-weighted estimate by one plus the `chr999` unassigned-extra-
DNA fraction, matching the recorded source-workflow policy. It gives a 4N
proxy mean of 4.98623 and an eight-mouse terminal 4N-origin mean of 2.32156
(53.44% lower); the corresponding 2N values are 2.29335 and 2.13534. This is a
descriptive cross-assay comparison, not a formal test or evidence about when
the reduction occurred. These downstream
artifacts do not establish a complete, versioned upstream NUMBAT inference in
this repository. External A03_Numbat storage contains additional upstream
artifacts, but their completeness and exact correspondence to all 16 canonical
CBS exports have not been established; no executable upstream NUMBAT workflow
is retained here.
The CBS filenames, byte sizes, and SHA-256 values are pinned in
`scRNAseq_Numbat/cbs_manifest.tsv` and
`scRNAseq_Numbat/injected_reference/reference_manifest.tsv`. Those matrices
and the canonical combined table form the executable downstream boundary, not
a reproducible raw NUMBAT run. The
deposited final RDS remains the verified fallback when the external H5 boundary
is unavailable.

The endpoint-ploidy table, sample workbook, and growth-curve workbook are
versioned source artifacts. Their revision and SHA-256 values are pinned in
`figure7_config.yaml`. Raw downloads and generated intermediates stay below
`Results/`; they are not publication inputs and are not committed.
If the reduced endpoint-ploidy TSV is missing while the 16 reviewed CBS files
are present, Manager regenerates the exact checksum-pinned table in the current
run's artifact directory and reuses it without modifying tracked inputs.

Every reusable stage has a dependency and output fingerprint. A subsequent run
reuses a stage only when its source hashes, relevant configuration, code,
runtime package contract, and output hashes validate. Stale automatically
managed caches are preserved with a `.stale.<timestamp>` suffix before
regeneration; explicit external paths are never modified. To require an
already-populated raw cache, add `--figure7-no-download-missing-raw`.

The byte-pinned v1 panel-7F reference is retained only for historical audit
because it was fitted from mixed human/mouse features. It can still be rendered
for comparison by explicitly supplying its directory, but that A-F run is
never publication eligible. The previously approved human-only v2 reference is
also superseded for inference because its nuisance term was derived from the
run-confounded endpoint CN score. The corrected full-workflow path retains
exact `GRCh38-` count rows before expression filtering, symbol resolution,
modeling, and GSEA. It writes a separate generated human-only,
initial-ploidy-adjusted reference with
`canonical_publication_allowed=false`; raw reruns do not inherit the approval
of the exact reviewed v3 bytes. GSEA starts with the configured
simple-permutation budget,
retries only unresolved pathways at increasing pinned budgets, recomputes BH
adjustment across each complete collection, and fails closed if any pathway
still lacks finite statistics at the configured cap. The generated human-only
panel displays only pathways with collection-wide BH-adjusted P <= 0.05, then
takes up to four pathways in each direction and collection. It never backfills
a direction with nonsignificant pathways, so generated collection and panel row
counts may be smaller than the historical frozen 8/24-pathway layout. Panel K
must use only the checksum-pinned endpoint table and the recorded raw
mouse-level Pearson/unrestricted-permutation contract; temporary historical
reconstructions are not publication inputs.

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

After scientific review, the exact five consumed files per endpoint (7A plot,
7C plot/test, and 7E plot/test) plus the two run configs are frozen in
`Data/in-vivo/figure7/saved_tgi_sensitivity/tgi_day24_day31_curated_cbs_v3_raw_pearson/`.
Canonical SI8 assembly reads that compact tracked bundle by default; the full
40-file result runs are not publication dependencies.

The earlier `tgi_day24_day31_all_cbs_v1` bundle was removed because it
summarized cells outside the final QC tumor universe. The assembler rejects
that obsolete policy and accepts only the curated raw-Pearson v3 publication
source.

Assemble the reviewed composite with:

```bash
Rscript Code/in-vivo/figure7/assemble_tgi_sensitivity.R \
  --output-dir=figures/Figure7_Supplement
```

The assembler validates the reviewed raw mouse-level panel-K contract in both
runs, then creates a five-panel Day-24/Day-31 composite plus a
portable checksum provenance table. It also writes the separate standard-schema
`si8_manifest.tsv` for the final SI8 PDF/PNG. The pre-existing 13-row
`manifest.tsv` remains the Day-24 Figure 7 source-panel contract; mixing SI8
rows into it would conflate two products and Manager would overwrite them. The
assembler does not refit either endpoint.
After a provenance-only code correction, pass `--provenance-only=true` to
rebuild the provenance and SI8 manifest without rewriting the reviewed PNG or
PDF bytes.

The selected endpoint day is used consistently for the matched-control TGI
calculation, plot labels, statistical tables, run metadata, panel contract, and
day-bearing filenames. Panel 7B and panel 7F are scientifically independent of
the TGI endpoint and are regenerated unchanged into the selected destination.

Routine and raw-fallback runs do not refresh either tracked frozen reference.
Full-workflow instead writes a separately identified
`state_pathway_grch_human_only_initial_ploidy_day17_v3_candidate` generated reference below the run
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
audit chain; it is not a substitute for the reviewed human-only v3 reference.

To explicitly generate and materialize only 7A-7E:

```bash
bash Manager.sh --mode standard --modules in_vivo_figure7 \
  --figure7-panels-ae-only --run-id <run_id>
```

This mode records a five-panel contract and does not read, validate, render, or
materialize panel 7F. Each included panel is written in both PDF and PNG format.

The superseded v2 directory remains byte-pinned and is exercised by its
dedicated audit validator and tests. Routine rendering intentionally refuses
to substitute it for reviewed v3, preventing an endpoint-CN-score-adjusted
panel from re-entering a canonical composite.

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
