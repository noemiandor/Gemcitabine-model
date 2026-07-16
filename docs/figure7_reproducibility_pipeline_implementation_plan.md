# Figure 7 Reproducibility Pipeline Implementation Plan

## Status

Drafted on 2026-07-16 for promotion to `main` after implementation and validation.

This plan defines a narrow manuscript workflow for the six approved Figure 7 source panels. It follows the `main`-branch Figure 1–6 pattern: a manager-controlled module writes an immutable result run, input/output manifests record provenance, selected source panels are materialized under `figures/Figure7/`, and the figure manifest is validated. It does not automate final composite assembly because the existing Figure 1–6 workflow does not assemble the Overleaf composites either.

## Objective

Promote only the code, compact inputs, configuration, and tests needed to reproduce the approved Figure 7 panels A–F. Integrate them into the root `Manager.sh` workflow without merging the divergent TaoLi/in-vivo branch history or importing the large historical output trees.

The routine workflow must:

1. generate exactly six Figure 7 source-panel PDFs;
2. preserve the approved Day-17 TGI, matched-control, ploidy, ETP-threshold, and pseudotime-state choices;
3. record all inputs and outputs with checksums;
4. support rematerialization from an immutable result run with `--mode panels-only`;
5. fail on a missing, duplicated, or unexpected Figure 7 panel;
6. leave Figure 1–6 materialization unchanged.

Statistical tables, panel plotting-data tables, logs, and metadata are required provenance and are not counted as extra panels. No unrequested PDF, PNG, or other figure may be emitted by the Figure 7 module.

## Frozen Figure 7 panel contract

The screenshot `Screenshot 2026-07-16 at 5.18.38 PM.png` is the visual reference. The canonical source panels and settings are:

| Panel | Approved content | Current source provenance | Frozen settings | Canonical module output |
|---|---|---|---|---|
| 7A | Tumor-growth trajectories explaining the Day-17 TGI calculation | `04h_tgi_calculation_plots_util.R`; current `TGI_calculation_growth_trajectories.pdf` | endpoint Day 17; mean of initial-ploidy-matched untreated controls | `figures/panel_7A_day17_tgi_calculation.pdf` |
| 7B | Selected CellCycle mean-ECDF comparisons | initial-ploidy `CellCycle_direct_group_ecdf_comparisons_selected_3panel.pdf` | original panel IDs 1, 8, and 9, corresponding to grid positions `(1,1)`, `(3,2)`, and `(3,3)` | `figures/panel_7B_cellcycle_selected_ecdf_comparisons.pdf` |
| 7C | Day-17 TGI in initial 2N versus 4N treated tumors | initial-ploidy `CellCycle_TGI_group_boxplot.pdf` | independent tumors; dose-stratified group-label permutation; not a paired-mouse test | `figures/panel_7C_day17_tgi_by_initial_ploidy.pdf` |
| 7D | Within-dose-centered TGI/ECDF-shift association | ETP reference-balanced `CellCycle_TGI_association_within_dose_centered.pdf` | ETP threshold 2.24; Day-17 mean-control TGI | `figures/panel_7D_day17_tgi_vs_centered_ecdf_shift.pdf` |
| 7E | Day-17 TGI versus sample mean ETP | ETP reference-balanced `CellCycle_TGI_AUC_vs_mean_ETP.pdf` | ETP threshold 2.24; Day-17 mean-control TGI | `figures/panel_7E_day17_tgi_vs_mean_etp.pdf` |
| 7F | Pathway activity across the accumulated CellCycle pseudotime state | TaoLi `04i` reference-balanced ETP 2.24 `primary_state_pathway_activity_heatmap.pdf` | accumulated interval `[0.30, 0.49]`; ETP group threshold 2.24; top positive and negative pathways per collection using the approved activity table | `figures/panel_7F_pseudotime_state_pathway_activity.pdf` |

The misleading legacy `AUC` token in the source filename for panel 7E must not appear in the manuscript-facing asset name or caption. The plotted endpoint is Day-17 TGI, not AUC TGI.

## Current implementation sources and promotion policy

The relevant implementations are split across histories that diverged from `main` before the current manuscript manager and output-contract work:

- A–E: local `codex/pseudotime-tgi-state-analysis` at `77caec93`;
- F: `origin/TaoLi` at `ea9e669b`;
- Figure 1–6 manager contract: `main` at `52bca332` when this plan was drafted.

Implementation must start on a fresh branch from the then-current `main`. Do not merge either feature branch and do not cherry-pick their broad commits. Those histories include unrelated in-vivo scripts, large generated `Figs/` trees, cached files, and changes that predate the current manager architecture.

Transplant or refactor only the narrow call graphs needed by panels 7A–F. In particular:

- do not import `Rplots.pdf`, `__pycache__`, report HTML, existing `Figs/pseudotime_TGI/**`, or `Figs/pseudotime_state_pathways/**` as analysis code;
- do not replace the `main` version of `Code/in-vivo/Utils.R` with the divergent TaoLi version;
- move the two `04i` dependencies `get_assay_matrix()` and `clean_gene_symbols()` into module-local helpers;
- preserve attribution in comments for functions extracted from the `04h` and `04i` implementations.

## Target architecture

### One manager module with two internal analysis stages

Add one manuscript module named `in_vivo_figure7` rather than overloading the existing legacy `in_vivo` module. The entrypoint should be:

```text
Code/in-vivo/figure7/run_figure7.R
```

Recommended module layout:

```text
Code/in-vivo/figure7/
├── run_figure7.R
├── figure7_config.yaml
├── src/
│   ├── common_io.R
│   ├── tgi_data.R
│   ├── tgi_statistics.R
│   ├── tgi_panels.R
│   ├── state_pathway_analysis.R
│   └── state_pathway_panel.R
└── tests/
    ├── test_figure7_contract.R
    ├── test_figure7_statistics.R
    ├── test_figure7_panels_only.R
    └── fixtures/
```

The entrypoint should orchestrate two internal stages but expose one output contract:

- the TGI/pseudotime stage computes panels 7A–7E;
- the state-pathway stage supplies panel 7F.

The two stages must not invoke the current broad `04h` or `04i` entrypoints because those write many unrequested figures. Refactor the underlying pure data/statistics/plot functions and call only the six approved panel builders.

### Analysis, rendering, and materialization modes

The module needs two manager-driven analysis modes and one standalone rendering mode:

1. `standard`:
   - recompute panels 7A–7E from the two tracked cell-level analysis tables;
   - render panel 7F from a frozen, tracked pathway-activity plotting table produced by the approved full `04i` run;
   - verify that all frozen settings and table checksums match the config;
   - this is the routine manuscript mode and must not require the large external Seurat object.
2. `full-analysis`:
   - recompute panels 7A–7E as in `standard`;
   - recompute the panel 7F gene model, GSEA, fitted trajectory, pathway selection, and activity table from raw RNA counts in an explicitly supplied Seurat RDS;
   - compare the newly generated activity table and pathway ordering with the frozen reference before replacing any canonical result;
   - record the RDS checksum, package versions, gene-set release, feature/species policy, and all model parameters.
3. standalone module `render-only`:
   - read plotting-data and statistics tables from an existing immutable Figure 7 run;
   - regenerate the same six PDFs into a new, explicitly supplied output directory without rerunning statistical tests or pathway modeling;
   - never modify the source run and never use a mutable `latest.txt` pointer.

Preserve the root manager's existing semantics: `Manager.sh --mode panels-only` must not invoke module analysis or rendering. It should locate a concrete existing Figure 7 result run and copy its six existing PDFs into `figures/Figure7/` through the materializer.

Add a generic `--source-run-id ID` option for `Manager.sh --mode panels-only`. `--run-id` identifies the new manager/materialization operation, while `--source-run-id` identifies the immutable module runs to materialize. This avoids reusing or overwriting the original manager run directory. Reject `--source-run-id` outside panels-only mode and reject `latest.txt` or arbitrary directory inference. Generated figure-manifest rows must record the source analysis run ID in their `run_id` field; the fresh materialization-operation ID belongs in manager metadata or a separate notes/operation field, never in place of the source run ID.

Use one explicit interface for optional full pathway recomputation:

```text
--mode full-refit --modules in_vivo_figure7 \
  --figure7-full-analysis --figure7-seurat-rds /absolute/path/object.rds
```

Require both Figure 7 flags together. `full-refit` alone must not make the default multi-module manager unexpectedly require a Seurat RDS. Document the interface in `Manager.sh --help`; never infer an RDS from a developer-specific `/Volumes` or HPC path.

Define Figure 7 behavior for every existing global manager mode:

| Manager mode | Figure 7 behavior |
|---|---|
| `check-only` | Validate routine inputs; validate the RDS and pinned gene sets as well when both full-analysis flags are supplied. Do not run analysis. |
| `saved-fit` | Same Figure 7 behavior as `standard`: recompute 7A–7E and render 7F from the frozen compact state-pathway analysis. |
| `standard` | Recompute 7A–7E and render 7F from the frozen compact state-pathway analysis. |
| `full-refit` without both Figure 7 full-analysis flags | Same Figure 7 behavior as `standard`; do not require the RDS. |
| `full-refit` with both Figure 7 full-analysis flags | Recompute 7F from the RDS, write the full compact audit chain, and compare against the frozen reference. |
| `panels-only` | Run no Figure 7 code; materialize the six existing PDFs from `--source-run-id`. |

### Canonical run layout

Each run should be rooted at:

```text
Results/in-vivo/figure7/runs/<run_id>_figure7/
```

Required contents:

```text
<run_id>_figure7/
├── figures/
│   ├── panel_7A_day17_tgi_calculation.pdf
│   ├── panel_7B_cellcycle_selected_ecdf_comparisons.pdf
│   ├── panel_7C_day17_tgi_by_initial_ploidy.pdf
│   ├── panel_7D_day17_tgi_vs_centered_ecdf_shift.pdf
│   ├── panel_7E_day17_tgi_vs_mean_etp.pdf
│   └── panel_7F_pseudotime_state_pathway_activity.pdf
├── tables/
│   ├── panel_7A_plot_data.tsv
│   ├── panel_7B_plot_data.tsv
│   ├── panel_7B_tests.tsv
│   ├── panel_7C_plot_data.tsv
│   ├── panel_7C_test.tsv
│   ├── panel_7D_plot_data.tsv
│   ├── panel_7D_test.tsv
│   ├── panel_7E_plot_data.tsv
│   ├── panel_7E_test.tsv
│   ├── panel_7F_plot_data.tsv
│   ├── panel_7F_selected_pathway_gsea.tsv
│   ├── panel_7F_leading_edge_genes.tsv
│   ├── state_pathway_gene_ranking_complete.tsv
│   ├── state_pathway_gsea_complete.tsv
│   ├── state_pathway_sample_bin_coverage.tsv
│   ├── state_pathway_design_qc.tsv
│   └── state_pathway_frozen_reference_comparison.tsv
├── metadata/
│   ├── input_manifest.tsv
│   ├── output_manifest.tsv
│   ├── run_config.tsv
│   ├── panel_contract.tsv
│   ├── session_info.txt
│   └── state_pathway_provenance.tsv
└── logs/
    ├── stdout.log
    └── stderr.log
```

Full-analysis mode must add the complete compact audit chain under `tables/state_pathway_full_analysis/`: the modeled-gene ranking, complete GSEA results for every configured collection, sample/bin coverage, design/model QC, selected-pathway leading edges and activity, and frozen-reference comparison. It must still emit only the six PDFs listed above. Standard mode must copy the corresponding frozen compact reference tables into its run so its panel-F selection remains auditable without the external RDS. Standard and render-only modes must never copy or regenerate the other `04h`/`04i` figures. Manager panels-only may only materialize those six existing PDFs.

## Input contract

### Routine tracked inputs

Promote these compact analysis inputs to `main` if they are not already there:

```text
Data/in-vivo/figure7/processed/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv
Data/in-vivo/figure7/processed/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv
Data/in-vivo/figure7/saved_state_pathway/<reference_id>/panel_7F_pathway_activity_plot_data.tsv
Data/in-vivo/figure7/saved_state_pathway/<reference_id>/panel_7F_selected_pathway_gsea.tsv
Data/in-vivo/figure7/saved_state_pathway/<reference_id>/panel_7F_leading_edge_genes.tsv
Data/in-vivo/figure7/saved_state_pathway/<reference_id>/state_pathway_gene_ranking_complete.tsv
Data/in-vivo/figure7/saved_state_pathway/<reference_id>/state_pathway_gsea_complete.tsv
Data/in-vivo/figure7/saved_state_pathway/<reference_id>/state_pathway_sample_bin_coverage.tsv
Data/in-vivo/figure7/saved_state_pathway/<reference_id>/state_pathway_design_qc.tsv
Data/in-vivo/figure7/saved_state_pathway/<reference_id>/state_pathway_provenance.tsv
```

The two cell-level tables are approximately 2.7 MB and 6.4 MB in the current branch and are sufficient for routine reconstruction of panels 7A–7E. Promote reviewed copies into the Figure 7 `processed/` namespace rather than making the module depend on feature-branch locations. The full pathway-activity table for panel 7F was written by the `04i` workflow but was not committed in `ea9e669b`; it must be exported from the canonical TaoLi/HPC result before implementation can be considered complete.

Use an immutable, versioned `<reference_id>` directory for the saved state-pathway analysis, analogous to the saved-fit directories used by the PKPD module. Freeze that ID in `figure7_config.yaml`; do not select the reference through a mutable `latest` pointer.

Create `panel_7F_pathway_activity_plot_data.tsv` as a panel-ready frozen subset, not as an unqualified copy of the full activity table. It must contain only the displayed pathways and must include stable `collection_id`, `collection_display_order`, `pathway_id`, `pathway_display_order`, `selected_direction`, and `selected_rank_within_direction` columns. The routine renderer must use those explicit order columns and must not rerank or break ties at render time. Retain the selected GSEA rows and leading-edge membership in the adjacent tracked audit tables.

The complete compact gene ranking and all candidate GSEA rows are mandatory frozen-analysis inputs, not optional archival detail. They make it possible to verify that the displayed pathways truly satisfy the top-four-per-sign selector. The sample/bin coverage and design QC tables document the experimental-unit support and fitted model behind that ranking. The external Seurat object and large pseudobulk count matrix remain external; these compact derived tables do not.

`state_pathway_provenance.tsv` must identify at least:

- the full-analysis run directory and report identifier;
- the exact `04i` code revision;
- the Seurat RDS SHA-256 and assay/layer;
- CellCycle and NonCellCycle metadata SHA-256 values;
- frozen interval configuration SHA-256;
- ETP method and threshold;
- spline degrees of freedom, number of bins, minimum cells per bin, grid size, and seed;
- gene-set source, release/version, species, collections, and size filters;
- feature/species handling policy;
- the GSEA ranking and pathway-selection rule;
- the activity-table SHA-256.

`Manager.sh::input_paths_for_module()` must enumerate `figure7_config.yaml`, both cell-level CSVs, every tracked panel-F table above, and `state_pathway_provenance.tsv`. Full-analysis mode must additionally enumerate the external RDS and pinned gene-set artifact. These are the exact files passed to the shared input-manifest writer; no scientifically relevant input may exist only in a command comment or README.

The module must fail if the observed inputs disagree with `figure7_config.yaml` or the provenance table.

### Optional raw-input refresh for A–E

The source branch can regenerate the two cell-level tables from:

```text
Data/in-vivo/scvelo_cell_metrics.csv
Data/in-vivo/all_ploidy.tsv
Data/in-vivo/sample_info.xlsx
Data/in-vivo/dt_Gem_VT_20241223_v4.xlsx
```

Keep this as an explicit input-refresh workflow, not a hidden side effect of the routine Figure 7 run. A refresh must write to:

```text
Results/in-vivo/figure7_input_refresh/runs/<run_id>/{tables,metadata,logs,qc}/
```

It must compare schemas/checksums and sample-level values with the frozen tables and require review before promotion. Reviewed refresh outputs may be deliberately copied into `Data/in-vivo/figure7/processed/`; the refresh command must never overwrite reusable `Data/` inputs directly.

### Optional full-analysis input for F

Full-analysis mode additionally requires:

```text
--figure7-seurat-rds /absolute/path/to/integrated_...rds
```

The RDS must be treated as an external, checksummed input. It must never be copied into the repository or result run. A missing RDS is a hard failure in full-analysis mode and irrelevant in standard, render-only, and manager panels-only modes.

For long-term reproducibility, freeze the MSigDB gene-set tables or record a versioned local gene-set artifact. A live `msigdbr` query without an audited release is not sufficient for canonical full recomputation.

## Frozen statistical and display specifications

### Shared A–E TGI specification

Freeze these values in `figure7_config.yaml` and validate them at runtime:

```yaml
tgi:
  outcome: day
  day: 17
  matched_control_summary: mean
  matched_control_group: initial_ploidy
  treated_doses_mg_per_kg: [30, 120]
statistics:
  seed: 1
  permutations: 10000
  exact_enumeration_when_feasible: true
etp:
  method: reference_balanced
  threshold: 2.24
```

Panel-specific rules:

- 7A: show individual mouse baseline-adjusted growth trajectories, facet by initial ploidy, overlay the mean matched untreated-control trajectory, and highlight Day 17.
- 7B: compute the full internal direct-comparison object if required by the test code, but add a numeric `comparison_id` and retain only IDs 1, 8, and 9 with labels `1. 0 vs treated`, `8. 4N: 0 vs treated`, and `9. 2N: 0 vs treated` in that order. Do not write the full 11-panel PDF.
- 7C: compare treated initial-2N and initial-4N tumors as independent groups. Use a dose-stratified label permutation. Do not describe or implement this as a paired boxplot.
- 7D: use the ETP reference-balanced 2.24 grouping to build equal-sample untreated references and correlate dose-centered ECDF RMSE with dose-centered Day-17 TGI.
- 7E: correlate sample mean endpoint ploidy with Day-17 TGI and retain the 2.24 ETP grouping only for shape/display metadata.

Exact enumeration supersedes the nominal 10,000 Monte Carlo permutations when the finite assignment space can be enumerated. The approved panels currently use 36 arrangements for 7C, 576 within-stratum assignments for 7D, and 40,320 TGI-label assignments for 7E.

Do not trust the embedded `TGI_percent_Day_17` column without verification. At runtime, recompute each treated mouse's value as:

```text
100 * (1 - tumor_volume_delta_Day_17 /
  mean(tumor_volume_delta_Day_17 of initial-ploidy-matched untreated mice))
```

Fail if the recomputed and embedded values differ beyond the configured numerical tolerance. Assert that panels 7C–7E contain exactly the eight treated mice at 30 or 120 mg/kg; untreated mice may define references and appear in panels 7A/7B but must not enter the panel 7C–7E outcome rows.

### Panel F specification

Freeze these values:

```yaml
state_pathways:
  model: ETP_reference_balanced_threshold_2_24
  accumulated_interval: {start: 0.30, end: 0.49, include_start: true, include_end: true}
  left_neighbor: {start: 0.11, end: 0.30, include_start: true, include_end: false}
  right_neighbor: {start: 0.49, end: 0.68, include_start: false, include_end: true}
  assay: RNA
  counts_layer: counts
  pseudotime_bins: 20
  minimum_cells_per_sample_bin: 10
  expression_filter: "CPM > 1 in at least 6 retained sample-bin observations"
  normalization: TMM
  observation_model: voom
  mouse_block: duplicateCorrelation
  spline_df: 5
  nuisance_terms: [dose_factor, ETP_reference_balanced_threshold_2_24_group]
  treatment_by_pseudotime_interaction: false
  empirical_bayes: robust
  grid_size: 501
  contrast: "mean(primary grid) - 0.5 * mean(left grid) - 0.5 * mean(right grid)"
  gsea_rank_statistic: moderated_t
  pathway_activity: "row-standardize fitted expression per leading-edge gene, then average genes"
  collections: [H, C2:CP:REACTOME, C5:GO:BP]
  top_positive_per_collection: 4
  top_negative_per_collection: 4
  pathway_selector: "finite adjusted P; no FDR cutoff; top four per sign and collection"
```

The standard renderer must use the pathway and collection order stored in the approved panel-ready table; it must not rerank pathways at render time. The dotted vertical lines must be placed at 0.30 and 0.49. Full analysis uses TMM-normalized sample-by-bin pseudobulk counts, a common pseudotime spline plus dose and ETP-group nuisance terms, mouse-blocked correlation, robust empirical Bayes, and no treatment-by-pseudotime interaction. The primary contrast is the equal-grid mean inside the accumulated interval minus equally weighted left and right neighboring means.

Full-analysis pathway reselection must reproduce the original selector exactly within each collection:

1. retain positive rows with `NES > 0` and finite adjusted P, with no FDR-significance cutoff;
2. order positives by adjusted P ascending, NES descending, then pathway ID ascending, and retain four;
3. retain negative rows with `NES < 0` and finite adjusted P, again with no FDR-significance cutoff;
4. order negatives by adjusted P ascending, NES ascending, then pathway ID ascending, and retain four;
5. concatenate the positive and negative selections using the frozen collection/pathway display-order convention.

The lack of an FDR cutoff is intentional: the approved heatmap includes the highest-ranked positive Hallmark pathways even though no positive Hallmark pathway passes FDR 0.05.

Prespecify the full-analysis comparison with the frozen reference:

- exact match for collection IDs, pathway IDs, display order, selected direction/rank, leading-edge gene IDs, and the 501-point pseudotime grid;
- `abs(new - reference) <= 1e-10 + 1e-6 * abs(reference)` for NES and adjusted P values;
- maximum absolute standardized-activity difference `<= 1e-6` and activity RMSE `<= 1e-7` over keyed pathway/grid rows;
- a machine-readable comparison table containing every keyed difference and threshold;
- a nonzero exit plus a review artifact when any threshold fails.

Keep these tolerances in the frozen config. They may be changed only through an explicitly reviewed reference update, not during a failing run.

Do not silently change the current feature/species policy during promotion. Before approving a new full-analysis reference, explicitly resolve the mixed human/mouse feature issue identified in review. If a corrected tumor-species filter or ortholog mapping changes panel 7F, treat that as a scientific figure revision requiring approval, not as an implementation-only change.

## Root manager integration

Modify `Manager.sh` on the fresh main-based branch as follows:

1. add `in_vivo_figure7` to the accepted module names;
2. add `Results/in-vivo/figure7/runs/<run_id>_figure7` to `module_run_dir()`;
3. add the routine tracked inputs to `input_paths_for_module()`;
4. conditionally add the explicit Seurat RDS in full-analysis mode;
5. add `command_for_module()` arguments for the config, mode, inputs, output directory, and optional RDS;
6. include `in_vivo_figure7` in the `panels-only` branch without invoking its entrypoint;
7. expose and validate any Figure 7-specific CLI flags in `usage()` and argument parsing;
8. add generic `--source-run-id` resolution for panels-only materialization while keeping the current manager operation under a distinct `--run-id`;
9. use the existing `run_module()` path so logs, input/output manifests, immutable run handling, and `latest.txt` metadata follow the existing contract.

Once all routine inputs are tracked and standard mode passes on a clean checkout, append `in_vivo_figure7` to the manager's default module list. Until then it may be opt-in, but Figure 7 promotion is not complete until the default manuscript command includes it or the repository explicitly documents why Figure 7 remains optional.

Add a row to `docs/manuscript_figure_module_registry.tsv` with:

- module: `in_vivo_figure7`;
- domain: `in-vivo`;
- entrypoint: `Code/in-vivo/figure7/run_figure7.R`;
- output root: `Results/in-vivo/figure7`;
- routine and optional full-analysis inputs;
- expected output: six source-panel PDFs plus tables and metadata;
- default mode: `standard` using the frozen pathway-activity table.

Do not reuse the current legacy `in_vivo` module name. Its `main`-branch command runs `pseudotimeAssociations.R` and has a different, pending output contract.

## Manuscript asset materialization

Add six `PANEL_SPECS` entries in `Code/tools/materialize_figure_assets.py`. Keep each module-run-relative source distinct from its manuscript-facing destination:

```text
source: figures/panel_7A_day17_tgi_calculation.pdf
  -> asset: figures/Figure7/panel_7A_day17_tgi_calculation.pdf
source: figures/panel_7B_cellcycle_selected_ecdf_comparisons.pdf
  -> asset: figures/Figure7/panel_7B_cellcycle_selected_ecdf_comparisons.pdf
source: figures/panel_7C_day17_tgi_by_initial_ploidy.pdf
  -> asset: figures/Figure7/panel_7C_day17_tgi_by_initial_ploidy.pdf
source: figures/panel_7D_day17_tgi_vs_centered_ecdf_shift.pdf
  -> asset: figures/Figure7/panel_7D_day17_tgi_vs_centered_ecdf_shift.pdf
source: figures/panel_7E_day17_tgi_vs_mean_etp.pdf
  -> asset: figures/Figure7/panel_7E_day17_tgi_vs_mean_etp.pdf
source: figures/panel_7F_pseudotime_state_pathway_activity.pdf
  -> asset: figures/Figure7/panel_7F_pseudotime_state_pathway_activity.pdf
```

The materializer should create `figures/Figure7/manifest.tsv` with exactly six generated-panel rows. For both standard and panels-only materialization, each row's `run_id`, `source_file`, and `result_run_dir` must consistently identify the same source analysis run.

While adding Figure 7, close two generic contract gaps:

1. require the complete expected panel set for every touched figure and reject missing or duplicate panel IDs, so a partial Figure 7 cannot validate silently;
2. pass the module run root into `generated_row()` explicitly instead of deriving `result_run_dir` from a source file's parent directory.

Completeness must be conditional on the module runs supplied to the materializer and each panel spec's `optional` flag. For example, Figure 3I remains optional when the `lci_overlays` module is not supplied. For the single `in_vivo_figure7` module, all six specs are mandatory whenever that module run is supplied. Define how manual/external rows participate in duplicate checks without requiring an optional generated source that was not selected.

Extend `validate_figure_manifest()` or add an equivalent generic expected-panel-set validation. Do not hard-code a one-off check that only works for Figure 7 if the generic manifest contract can express it cleanly.

Also make shared manifest validation portable across clean checkouts. When both paths are present, prefer an existing `repo_relative_path` resolved against the current checkout; fall back to `absolute_path` only for genuinely external inputs. Add regression tests demonstrating that existing Figure 1–6 manifests still validate after the repository is moved to a different absolute path.

## Documentation updates

Update these main-branch documents during implementation:

- `docs/FigureCodeMap.md`: add Figure 7 panels A–F, source functions, inputs, canonical result run, and manuscript-facing assets;
- `docs/manuscript_figure_module_registry.tsv`: add `in_vivo_figure7`;
- `docs/manuscript_figure_output_standardization_plan.md`: record Figure 7 as implemented and note the saved-analysis/full-analysis distinction;
- `figures/Figure7/manifest.tsv`: generated by the materializer, never handwritten;
- a module README under `Code/in-vivo/figure7/README.md`: exact commands, modes, dependencies, output inventory, and interpretation caveats.

The documentation must state that the pipeline materializes six source panels and does not assemble the final A–F composite, matching the Figure 1–6 behavior.

## Validation and regression tests

### Unit and static validation

Require:

- successful parsing of every new R file;
- `bash -n Manager.sh`;
- Python compilation/tests for changed materializer and manifest helpers;
- schema tests for every plotting/statistics table;
- tests that unknown modes, missing inputs, mismatched settings, and nonempty output directories fail clearly;
- a test that the Figure 7 module emits exactly six PDFs and no PNG/SVG/TIFF or extra PDF;
- a test that materialization fails for any missing, duplicated, or unexpected Figure 7 panel.

### Frozen numerical regression checks

At minimum, assert the currently approved values within explicit numerical tolerances:

- panel 7B numeric `comparison_id` values are exactly `1`, `8`, and `9`, with the approved string labels, in that display order;
- panel 7C adjusted difference `(4N - 2N)` is approximately `-39.5232168` percentage points and exact dose-stratified permutation `P = 0.0555556`;
- panel 7D Pearson `r = 0.8105173` and exact permutation `P = 0.0104167`;
- panel 7E Pearson `r = -0.6984010` and exact permutation `P = 0.0591270`;
- all A–E metadata report `TGI_percent_Day_17`, outcome `day`, matched-control summary `mean`, and Day 17;
- panel 7F contains the approved collection labels and pathway order, exactly eight selected pathways per collection unless the approved frozen table documents a different count;
- panel 7F vertical boundaries are exactly 0.30 and 0.49.

If a code cleanup changes these values, stop and determine whether it is a bug fix or a scientific revision. Do not update reference values merely to make tests pass.

### Visual regression checks

PDF byte hashes are not an appropriate sole regression target because embedded creation metadata can vary. Instead:

1. verify each PDF is nonempty, single-page, and parseable;
2. extract text and assert the expected title, axes, panel labels, and annotations;
3. rasterize at a fixed DPI and compare with approved reference rasters using a documented perceptual tolerance;
4. compare plotting-data and statistics tables deterministically across two clean runs.

### Manager integration checks

Run, with a fresh run ID:

```bash
bash Manager.sh \
  --mode check-only \
  --modules in_vivo_figure7 \
  --run-id figure7_check
```

Run the canonical standard analysis and materialization:

```bash
bash Manager.sh \
  --mode standard \
  --modules in_vivo_figure7 \
  --run-id <fresh_run_id>
```

Then rematerialize from that exact run with `--mode panels-only`. Verify:

```bash
bash Manager.sh \
  --mode panels-only \
  --modules in_vivo_figure7 \
  --source-run-id <canonical_analysis_run_id> \
  --run-id <fresh_materialization_run_id>
```

This command must copy existing PDFs from the immutable source run; it must not rerun R code or overwrite the source manager run. Verify:

- both module input/output manifests validate;
- `figures/Figure7/manifest.tsv` validates and contains exactly six generated rows;
- every manifest `result_run_dir` points to the Figure 7 module root;
- every manifest `run_id` equals `<canonical_analysis_run_id>`, not `<fresh_materialization_run_id>`;
- no source path comes from `Figs/`, `Downloads`, a report HTML, or `latest.txt`;
- the six rematerialized panels match the standard-run plotting data and approved visual references;
- Figure 1–6 materialization and manifests are unchanged.

For full-analysis validation, run:

```bash
bash Manager.sh \
  --mode full-refit \
  --modules in_vivo_figure7 \
  --figure7-full-analysis \
  --figure7-seurat-rds /absolute/path/to/object.rds \
  --run-id <fresh_full_analysis_run_id>
```

Confirm the new panel 7F plotting table agrees with the approved frozen table under the prespecified comparison rules. A difference must produce a review artifact and a nonzero exit unless replacement is explicitly authorized.

## Implementation sequence

### Phase 0: recover and freeze missing panel-F provenance

1. Export `primary_state_pathway_activity_heatmap_plot_data.csv` from the exact TaoLi/HPC run used for the screenshot.
2. Export the associated interval config, complete modeled-gene ranking, complete GSEA tables for all configured collections, model parameters, sample/bin coverage, design QC, package versions, gene-set version, cell-ID/count audit, and input checksums.
3. Freeze a panel-ready subset with explicit collection/pathway display order, selected direction/rank, selected GSEA rows, and leading-edge membership.
4. Confirm the frozen table regenerates the screenshot's panel 7F without rerunning the model.
5. Resolve or explicitly freeze the feature/species policy before approving a future full-analysis reference.

This phase is a hard prerequisite. The PDF alone is insufficient as the canonical panel-F analysis input.

### Phase 1: create the narrow module on current `main`

1. Create a fresh feature branch from current `main`.
2. Add `Code/in-vivo/figure7/` and the frozen config.
3. Extract/refactor the minimal 04h call graph for A–E.
4. Add the frozen-table renderer for F.
5. Add exact output-inventory and numerical regression tests.

### Phase 2: add optional full pathway recomputation

1. Extract the necessary 04i model/GSEA/activity functions into module-local files.
2. remove developer-specific paths and the dependency on divergent `Utils.R`;
3. require an explicit external RDS and pinned gene sets;
4. write the mandatory complete compact gene-ranking, GSEA, coverage, design-QC, provenance, and comparison tables;
5. ensure only panel 7F is rendered from this stage.

### Phase 3: connect the manuscript manager

1. Register `in_vivo_figure7` in `Manager.sh` and the module registry.
2. Add the six materializer specifications.
3. harden panel-set completeness and module-root provenance validation.
4. add `--source-run-id` panels-only support without changing the no-analysis semantics.
5. make repo-relative manifest validation portable across checkout locations.
6. update FigureCodeMap and module documentation.

### Phase 4: canonical run and promotion

1. Run check-only, unit, standard, panels-only, and full-analysis validations.
2. Review the six materialized assets against the approved screenshot.
3. Inspect `git diff --check`, manifests, output inventory, and worktree scope.
4. Commit one reviewed canonical standard Figure 7 result run under `Results/in-vivo/figure7/runs/`, including its six source PDFs, tables, metadata/manifests, and logs required by the output contract.
5. Materialize and commit `figures/Figure7/` plus its generated six-row manifest from that concrete committed run.
6. Keep code/config/tests and canonical data/result artifacts in separate commits where practical, but all are required for promotion.
7. Promote the reviewed commits to `main`; do not merge the broad source branches.

## Definition of done

Figure 7 is ready for `main` only when all of the following are true:

- a clean checkout can run the routine Figure 7 module without developer-specific paths;
- the module writes exactly six source-panel PDFs with the frozen settings;
- panel 7F has a tracked plotting table and complete full-analysis provenance, not only a PDF;
- panel 7F's complete compact gene ranking, candidate GSEA results, coverage, and design QC are retained so its displayed selection can be audited;
- full-analysis mode can recompute panel 7F when the external RDS is supplied;
- input/output manifests and the six-row Figure 7 manifest validate;
- one canonical standard Figure 7 result run is committed so the committed figure manifest validates on a clean checkout;
- manager panels-only rematerialization works from that concrete immutable run via `--source-run-id` without rerunning analysis;
- standalone render-only regeneration writes to a new output location and leaves the source run unchanged;
- numerical and visual regressions match the approved Figure 7 reference;
- Figure 1–6 behavior remains unchanged;
- no bulk TaoLi/04h/04i output trees or unrelated in-vivo code are promoted;
- the final diff has been independently reviewed for scientific settings, software architecture, and reproducibility.
