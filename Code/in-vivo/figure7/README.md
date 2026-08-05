# Figure 7 reproducibility module

This module generates six scientific source panels as matched vector PDF and
300-DPI PNG files and assembles the manuscript-facing A-L composite as
`Figure7_reviewed_GRCh.{pdf,png}`. The composite is built natively at the
7.1 x 10.645 inch full-page target rather than by shrinking a large-format
canvas. Its fixed six-row layout is A/B; C/D/E; F/G; H; I/J; K/L. The composite
reuses Supplementary Figure 4A-C/E and Supplementary Figure 7B through the
shared production plotting helper; their supplementary copies are retained.
The former human-only panel-7F v2
reference is retained for audit but is superseded because its model adjusted
for a run-confounded endpoint-CN-score group. The reviewed pointwise-v4 reference
adjusts for injected initial ploidy and directly uses the computed
0.296--0.486 density-support interval as the canonical panel-7F source. Panel
L shows the unadjusted mouse-level Pearson association between the
checksum-pinned mean endpoint tumor-cell ploidy and TGI, with unrestricted
exact enumeration of all 8! TGI-label permutations; endpoint ploidy is not
used as a nuisance covariate in the state-pathway model.

## Frozen routine analysis

- Panels 7A-7D are recomputed from the two tracked plot-facing cell-level
  tables under `Data/in-vivo/figure7/processed/`. Source panel 7E (main panel
  L) additionally reads the exact six-column, checksum-pinned
  `Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv`: 14,125 cells across 16 CBS
  files. This remains the complete immutable CBS inventory. Every scored cell
  must also occur, with the identical file, barcode, and value, in the exact
  9,832-cell QC-passed CellCycle + NonCellCycle union represented in the
  plot-facing tables. The other 4,293 inventory cells remain provenance-only
  and are excluded from scoring. Panel L uses the resulting 5,335 treated
  cells.
  That newer artifact has its own source revision (`dcdb62f2252...`) in the
  config; the surrounding legacy source revision does not claim it existed in
  the earlier snapshot.
- TGI is the Day-17 endpoint statistic, recalculated for each treated mouse using the mean Day-17 growth delta of untreated controls matched by initial ploidy.
- Panels 7C-7E contain exactly eight treated mice at 30 or 120 mg/kg. Untreated mice contribute references only.
- Source panel 7B computes the three equal-mouse ECDF comparisons used in main
  panel H and the equal-mouse Gaussian-density contrast displayed as SI4I. A common pooled,
  label-invariant bandwidth is pinned in `density_localization_config.yaml`.
  Exact treatment-label enumeration within injected-origin strata gives 4,900
  assignments: pointwise positive support spans 0.296--0.486, while
  studentized max-absolute-T family-wise support spans 0.414--0.426 (global
  exact P = 0.0473469). The raw density excess peaks at pseudotime 0.452;
  standardized max-absolute-T evidence is strongest at 0.420.
- Panel 7D (main panel K) uses equal-mouse untreated ECDF references matched by
  injected initial ploidy. The treated-mouse association is Pearson r =
  0.7399455 with exact within-dose permutation P = 0.0173611 (576 labelings).
- Panel 7F (main panel I) uses model
  `initial_ploidy_adjusted_grch_human_only_pointwise_interval_v4`; its nuisance terms are dose and
  injected initial ploidy, never endpoint CN score or an endpoint-derived
  threshold group. Its primary interval is computed directly as the unique
  connected positive pointwise-supported density region (0.296--0.486), and
  its two equal-width flanks are then derived as [0.106,0.296) and
  (0.486,0.676]. Static interval fields are validation expectations and cannot
  select a different modeling contrast.
- Panel 7E (main panel L) reports the descriptive mouse-level association
  between mean endpoint tumor-cell ploidy and Day-17 TGI (Pearson r =
  -0.6984010; asymptotic P = 0.0540070; exact unrestricted permutation P =
  0.0591270; 40,320 assignments). It is an unadjusted association, not
  evidence for an independent or causal terminal-ploidy effect.
- Main panel J is the live-grob cell-by-chromosome view of the exact 9,832-cell
  final-QC NUMBAT universe, with injected-origin, gemcitabine-dose, and mouse
  annotation bars and a compact key defining every annotation color. Its
  continuous copy-number palette uses a gold anchor at copy number 3 so that
  the common one-copy-gain state remains visible on a white page.
- The A-L manuscript composite binds the exact reviewed 11-table SI cache and
  displays in manuscript reading order: source 7A, source 7C, SI4A-C, SI4E, SI7B,
  source 7B, source 7F, the exact QC-filtered NUMBAT heatmap, and source 7D-E. Publication styling changes only
  layout, typography, reader-facing labels, and legend placement; it does not
  change the panel identities, source tables, fitted models, tests, or values.

The nine reviewed pointwise-v4 panel-7F files retain exact `GRCh38-` features before
expression filtering, symbol resolution, model fitting, and Homo sapiens GSEA.
They include the computed interval definition and use pseudotime 0.296--0.486.
Displayed pathways must
have collection-wide BH-adjusted P <= 0.05; the selector takes up to four per
sign and collection and never backfills with nonsignificant pathways. The
reviewed provenance authenticates all nine compact reference files,
including the interval definition, together with the exact input hashes, model/design audit,
feature-species boundary, MSigDB release, and
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
full A-L figure is requested, Manager automatically schedules `si_figures`
before `in_vivo_figure7`, even if only the latter was listed. It first validates
corrected generated caches; the reviewed frozen panel-7F reference and reviewed
SI plot-only cache do not satisfy an explicit full refit. A complete,
lineage-valid generated cache avoids raw-data access. When a missing downstream
stage needs the final Seurat object, select the source boundary explicitly:

- `--figure7-scrna-source rds` validates or downloads the deposited final RDS
  and skips `Code/in-vivo/scRNA_Seq_analysis`;
- `--figure7-scrna-source h5` requires the reviewed 18-sample Cell Ranger
  inventory, runs the cluster-SIF/full-SIF standalone workflow, and accepts its
  final RDS only after the semantic audit against the deposited Zenodo RDS
  passes. Metadata, cluster assignments, identities, count matrices, and graphs
  must match exactly; floating assay and PCA values use the recorded numerical
  tolerances, and UMAP uses same-axis correlation plus per-cell displacement
  limits. Seurat command timestamps are non-gating, while command parameters,
  calls, assays, and seeds must match exactly. Downstream reuse additionally
  validates the standalone `run_manifest.tsv`, `final_artifact_runtime.tsv`,
  `PIPELINE_COMPLETE.txt`, and `rds_semantic_audit/AUDIT_COMPLETE.txt` chain.
  A passed H5 audit makes that generated RDS the mandatory input to both
  Supplementary Figures and Figure 7; failure stops the run without falling
  back to the deposited RDS.

The complete deposited bundle includes 18 loom files, one final Seurat RDS,
18 Cell Ranger H5 files, and 11 support/provenance files. All 48 manifest rows
are size/MD5 validated when selected. For every Zenodo-hosted object, including
support/provenance files, the Zenodo manifest/API size and MD5 are the
authoritative integrity contract; narrative values inside a support document
are not used as replacement checksums. The canonical H5 layout is
`Data/in-vivo/figure7/raw/zenodo_21463392/SUM-159/A02_cellRanger/*-Count-HM/outs/*-Count-HM_filtered_feature_bc_matrix.h5`.
Zenodo publishes the H5 files with flat basenames; the downloader materializes
each one under its sample-specific canonical directory above.

Figure 7A-7D additionally reuse or download the 18 deposited loom files. Source
panel 7E/main panel L binds the complete combined CBS inventory; in
`full-refit`, Manager regenerates that run-scoped table from the manifest-pinned
16 CBS matrices and requires it to reproduce the canonical checksum. It then
restricts scoring to exact file+barcode keys retained in the final Seurat tumor
universe and represented by the two processed tables. The
complete Zenodo deposit is about 11.56 GiB. The RDS workflow selects the loom,
RDS, and support roles without downloading the alternative H5 inputs. Manager
then runs only the missing
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
   into Figure 7 and assemble the A-L review candidate.

The generated-cache handoff is required in `full-refit`: Figure 7 will not
silently fall back to the reviewed supplementary cache. The resulting
`Figure7_generated_GRCh_candidate.{pdf,png}` pair is explicitly noncanonical.
Only the exact reviewed cache may produce `Figure7_reviewed_GRCh.{pdf,png}` in
routine mode.

To copy a completed H5 full-refit candidate into a reviewable figure tree,
use `--publish-generated-candidate` together with an explicit isolated
`--figure-root`, for example `figures/H5_fullrefit_<run_id>`. This opt-in never
writes the canonical `figures/` tree, never updates canonical `latest`
pointers, and records `canonical_publication_allowed=false` in the published
manifest rows. Candidate publication revalidates the complete semantic-audit
report set and the bound generated/Zenodo RDS identities.

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
main Figure 7J chromosome-state view and recompute every value in the canonical
`scRNAseq_Numbat/all_ploidy.csv` (and its reduced `all_ploidy.tsv` projection).
Two checksum-pinned, project-designated lineage-matched karyotype proxies
provide the 2N-A7M and 4N-A5M reference distributions used by SI6E; they are
not the same-passage A6M/A4M inocula. The reference calculation adds the
`chr999` value, expressed in haploid-genome-equivalent units of unassigned DNA,
to the autosomal length-weighted estimate. It gives a 4N proxy mean of 3.94651
and an eight-mouse terminal 4N-origin mean of 2.31855 (41.25% lower); the
corresponding 2N values are 2.15124 and 2.13551. This is a
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
`figure7_config.yaml`. Deposited raw inputs stay below
`Data/in-vivo/figure7/raw/zenodo_21463392`; generated intermediates stay below
`Results/`. Neither is a publication input or committed.
If the reduced endpoint-ploidy TSV is missing while the 16 reviewed CBS files
are present, Manager regenerates the exact checksum-pinned table in the current
run's artifact directory and reuses it without modifying tracked inputs.

Every reusable stage has a dependency and output fingerprint. A subsequent run
reuses a stage only when its source hashes, relevant configuration, code,
runtime package contract, and output hashes validate. Stale automatically
managed caches are preserved with a `.stale.<timestamp>` suffix before
regeneration; explicit external paths are never modified. To require an
already-populated raw cache, add `--figure7-no-download-missing-raw`.

The previously approved human-only v2 reference is superseded for inference
because its nuisance term was derived from the run-confounded endpoint CN
score. The corrected full-workflow path retains
exact `GRCh38-` count rows before expression filtering, symbol resolution,
modeling, and GSEA. It writes a separate generated human-only,
initial-ploidy-adjusted reference with
`canonical_publication_allowed=false`; raw reruns do not inherit the approval
of the exact reviewed pointwise-v4 bytes. GSEA starts with the configured
simple-permutation budget,
retries only unresolved pathways at increasing pinned budgets, recomputes BH
adjustment across each complete collection, and fails closed if any pathway
still lacks finite statistics at the configured cap. The generated human-only
panel displays only pathways with collection-wide BH-adjusted P <= 0.05, then
takes up to four pathways in each direction and collection. It never backfills
a direction with nonsignificant pathways, so generated collection and panel row
counts therefore vary with the significant results. Main panel L
must use only the checksum-pinned endpoint table and the recorded raw
mouse-level Pearson/unrestricted-permutation contract; temporary unreviewed
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

The assembler validates the reviewed raw mouse-level endpoint-ploidy contract in both
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

Routine and raw-fallback runs do not refresh tracked reviewed or audit reference directories.
Full-workflow instead writes a separately identified
`state_pathway_grch_human_only_initial_ploidy_day17_pointwise_v4_candidate` generated reference below the run
intermediates and marks it noncanonical.

To explicitly generate and materialize only 7A-7E:

```bash
bash Manager.sh --mode standard --modules in_vivo_figure7 \
  --figure7-panels-ae-only --run-id <run_id>
```

This mode records a five-panel contract and does not read, validate, render, or
materialize panel 7F. Each included panel is written in both PDF and PNG format.

The superseded v2 directory remains byte-pinned and is exercised by its
dedicated audit validator and tests. Routine rendering intentionally refuses
to substitute it for reviewed pointwise-v4, preventing an endpoint-CN-score-adjusted
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
materializes the exact recorded A-L composite/source-panel or explicit A-E
contract from an explicit `--source-run-id`. Noncanonical runs are rejected.

The older artifact-driven `full-analysis` entrypoint remains only as an
explicit guard that directs callers to the corrected workflow. Manager's
`full-refit` path uses the raw-data fallback above, pins `msigdbr` and the
MSigDB release, applies the exact human-only feature policy, and keeps its
generated reference noncanonical until scientific review.

## Publication-scale composite and visual-QC package

The normal Figure 7 renderer is the authoritative compositor. It rebuilds the
A-L figure directly from live ggplot and heatmap grob objects at 7.1 x 10.645
inches, writing a 300-DPI PNG and a vector PDF. It never assembles the final
figure from exported panel rasters. The fixed six rows allocate extra display
area to the two heatmaps and the two mouse-level association panels while
preserving the enforced panel identity and manuscript reading order:

1. A/B
2. C/D/E
3. F/G
4. H
5. I/J
6. K/L

Panel H retains only the three equal-mouse ECDF comparisons on their shared,
unzoomed 0--1.05 scale. The density-localization view is reproduced as the
new full-width SI4I. Main panels I and J share the penultimate row, with J
showing the exact QC-filtered NUMBAT copy-number matrix.

The presentation-only audit package under `figures/Figure7/polishing/` can be
rebuilt independently with:

```bash
scripts/agentRrunner.sh \
  figures/Figure7/polishing/scripts/polish_figures.R --phase all
```

This command invokes the same scientific panel builders and frozen inputs as
the normal renderer; it does not refit models or define an alternative
scientific analysis. The package records the explicit A-L identity map, target
dimensions, adopted layout, optimizer diagnostic, rebuild command, input and
output hashes, byte-identity report, and print-size visual-QC assessment.
Audit subpanel PNGs are inspection artifacts only. Canonical manuscript
materialization remains the responsibility of `Manager.sh` and its recorded
run contract.

## Output contract

Each successful run has `figures/`, `tables/`, `metadata/`, and `logs/`. The
default six-panel source contract contains these PDF/PNG pairs:

1. `panel_7A_day17_tgi_calculation.{pdf,png}`
2. `panel_7B_cellcycle_selected_ecdf_comparisons.{pdf,png}`
3. `panel_7C_day17_tgi_by_initial_ploidy.{pdf,png}`
4. `panel_7D_day17_tgi_vs_centered_ecdf_shift.{pdf,png}`
5. `panel_7E_day17_tgi_vs_mean_etp.{pdf,png}`
6. `panel_7F_pseudotime_state_pathway_activity.{pdf,png}`

It additionally contains the publication-scale
`Figure7_reviewed_GRCh.{pdf,png}` pair, whose panel contract is: A=7A, B=7C,
C=SI4A, D=SI4B, E=SI4C, F=SI4E, G=SI7B, H=7B, I=7F, J=7J,
K=7D, and L=7E. The copy-number heatmap is published only as main J. The
fixed injected-origin/dose/mouse block order is preserved in J, while cells are
hierarchically clustered separately within each mouse using Euclidean distance
on the displayed chromosome profiles and Ward.D2 linkage. The reviewed SI
manifest hash, target dimensions, ordered mapping, and row-ordering policy are
recorded in run metadata and enforced during materialization. The exact
9,832-row display order is exported as
`tables/panel_7J_copy_number_row_order.tsv`.

A full-refit run has the same scientific panel mapping but writes
`Figure7_generated_GRCh_candidate.{pdf,png}`. Its generated supplementary-cache
manifest, analysis-input manifest, run configuration, and provenance are all
hash-bound in Figure 7 metadata; `canonical_publication_allowed=false` prevents
the candidate from being materialized as a canonical manuscript asset. It may
be copied only as an explicitly labeled generated candidate under an isolated
figure root.

An explicit `--panel-set=a-e`/`--figure7-panels-ae-only` run instead contains
exactly the first five pairs, records `panel_set=a-e`, and excludes all panel-F
inputs and outputs. It also excludes the A-L composite because the main figure
requires the reviewed state-pathway and SI-cache contracts together.

Plotting data, exact-permutation tests, the complete compact state-pathway audit chain, frozen-reference comparison, run settings, panel contract, and session information are retained alongside the PDFs.
Source panel 7B, whose ECDF component is main H and whose localization is SI4I,
additionally writes
`panel_7B_density_localization_grid.tsv`,
`panel_7B_density_localization_intervals.tsv`, and
`panel_7B_density_localization_test.tsv`; render-only mode requires and
hash-binds the separate density-localization configuration.

## Tests

```bash
Rscript Code/in-vivo/figure7/tests/testthat.R
```

The tests parse all module files, reproduce the frozen A-E numerical results,
enforce treated-only outcomes and selected ECDF IDs 1/8/9, reproduce the
equal-mouse density localization under all 4,900 origin-stratified assignments,
pin both support intervals and the simultaneous critical value, validate the
reviewed pointwise-v4 reference and its lineage, exercise the strict panel-F
reviewed-v2 publication guard and generated human-only contract, verify cached-stage
fingerprints and tamper rejection, validate the complete Zenodo manifest, and
confirm that missing raw inputs fail before output is created when downloading
is disabled.
