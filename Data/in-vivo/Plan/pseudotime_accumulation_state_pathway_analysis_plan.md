# Pseudotime Accumulation-State Pathway Analysis: Implementation Plan

## Objective

Identify the biological state in which gemcitabine-treated tumor cells accumulate along the existing CellCycle pseudotime trajectory.

The analysis will use the treatment-versus-control cell-density comparison only once, to locate and freeze the accumulated pseudotime region. Gene-expression modeling will then annotate that region relative to neighboring pseudotime states. It will not test treated versus control expression at a fixed pseudotime.

The intended conclusion has the following form:

> Gemcitabine-treated tumors contain an increased proportion of cells occupying a pseudotime state characterized by pathways X, Y, and Z.

This analysis will not, by itself, claim that treatment speeds up or slows down a pathway. Accumulation is compatible with a bottleneck or delayed exit, but also with altered entry, survival, or trajectory branching.

## Non-goals

The primary analysis will not:

- perform another treated-versus-control gene-expression contrast;
- optimize the interval endpoints to obtain stronger gene or pathway results;
- reroot, reverse, or rescale the established pseudotime;
- relate pathway scores to TGI;
- infer physical transition rates from cross-sectional cell abundance.

## Analysis estimand

The primary estimand is the transcriptional program associated with the already identified accumulation state:

```text
mean fitted expression within the frozen pseudotime interval
minus
mean fitted expression in the immediately adjacent pseudotime intervals
```

The fitted pseudotime effect will be common across treatment groups. Dose and initial ploidy may enter as additive nuisance terms, but there will be no treatment-by-pseudotime, dose-by-pseudotime, or ploidy-by-pseudotime interaction in the primary model.

Consequently:

- A positive gene statistic means that the gene characterizes the accumulated state more strongly than its neighboring states.
- A positive pathway normalized enrichment score (NES) means that the pathway characterizes the accumulated state.
- Neither result means that treatment directly upregulated the gene or pathway.

## Pre-specified pseudotime regions

The region must be frozen before inspecting expression or pathway results. The current equal-mouse density analysis supports the following pre-specification:

| Role | Pseudotime interval | Purpose |
|---|---:|---|
| Primary accumulated-state interval | `[0.30, 0.49]` | Connected region of pointwise treated-cell excess, rounded from approximately `0.296-0.486` |
| Left neighboring interval | `[0.11, 0.30)` | Equal-width state immediately before the primary interval |
| Right neighboring interval | `(0.49, 0.68]` | Equal-width state immediately after the primary interval |
| Broad sensitivity interval | `[0.11, 0.61]` | Broader region of numerical treated-cell excess |
| Strict-support sensitivity interval | `[0.414, 0.426]` | Narrow region supported by the simultaneous density-difference band; use only if cell coverage is adequate |
| Depleted-region negative control | `[0.80, 0.98]` | Region in which treated cells were depleted rather than enriched |

The primary adjacent-state contrast will give the two flanks equal weight, irrespective of their cell counts:

```text
primary-window mean - 0.5 * (left-flank mean + right-flank mean)
```

The means above refer to values evaluated on an equally spaced pseudotime grid, not cell-weighted means. This prevents the treatment-associated density shift from defining the expression contrast a second time.

Before implementation, record these endpoints in a machine-readable file such as:

```text
Code/in-vivo/04i_pseudotime_state_pathways_config.yaml
```

The record should also include:

- the pseudotime scenario and root orientation;
- the density method and common bandwidth;
- the density-analysis input checksum;
- the code revision used to derive the interval;
- the date on which the interval was frozen;
- an explicit statement that expression data were not inspected when selecting the interval.

If a frozen interval has insufficient expression coverage, do not move its endpoints after viewing pathway results. Report the coverage failure and use only the pre-specified sensitivity interval that passes the coverage rule.

## Existing repository components

### Cell-level pseudotime metadata

Use:

```text
Data/in-vivo/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv
```

The required columns already include:

```text
cell_id
sample_id
cluster
initial_ploidy
gemcitabine_dose
gemcitabine_dose_mg_per_kg
pseudotime
```

This file contains metadata and pseudotime but no gene-expression matrix. Only its CellCycle cells should enter this analysis. The NonCellCycle 04h input is not needed.

### Preferred expression source

Use raw RNA counts from the `RNA/counts` layer of the same Seurat object used upstream for the trajectory, resolving the first available object in the same order used by the trajectory scripts:

```text
<Results_root>/03c_cluster_annotation/04_objects/integrated_sct_cca_seurat_cluster_final_annotation.rds
<Results_root>/03b_manual_cluster_merge/objects/integrated_sct_cca_seurat_final_manual_merge.rds
<Results_root>/03_final_cluster/03_objects/integrated_sct_cca_seurat_final_reclustered.rds
```

The implementation should allow an explicit `--seurat-rds` override and should use `get_assay_matrix()` from `Code/in-vivo/Utils.R` for Seurat-version compatibility.

Do not use the existing scVelo H5AD `X` matrix as though it were raw counts. `Code/in-vivo/04f_cell_cycle_results.R` already checks this matrix and can find that it is non-integer-like, in which case it deliberately skips count-based pseudobulk differential expression. An H5AD source is acceptable only if a documented raw-count layer exists, is non-negative and integer-like, and maps to the same cells and genes.

### Reusable enrichment machinery

Reuse or extract the following behavior from `Code/in-vivo/03b_cluster_annotation_and_GSEA.R`:

- human Hallmark gene sets from `msigdbr`;
- gene-symbol cleaning;
- `fgseaMultilevel()` with a documented fallback;
- leading-edge gene export;
- NES/FDR tables and enrichment plots.

The focused signatures and H5AD sparse-reading helpers in `Code/in-vivo/04f_cell_cycle_results.R` are useful for secondary visualization, but its current treatment/dose-by-bin GSEA is not the primary analysis requested here.

## Proposed implementation

Add a standalone workflow rather than expanding the 04h TGI workflow:

```text
Code/in-vivo/04i_pseudotime_state_pathways.R
Code/in-vivo/04i_pseudotime_state_pathways_util.R
Code/in-vivo/04i_pseudotime_state_pathways_config.yaml
```

The separation is intentional: 04h establishes where cells accumulate, whereas 04i annotates what that state is.

### Command-line interface

Support at least:

```text
--cell-metadata=<path>
--seurat-rds=<path>
--config=<path>
--output-root=<path>
--assay=RNA
--counts-layer=counts
--n-pseudotime-bins=20
--min-cells-per-sample-bin=10
--spline-df=5
--gene-set-collections=H
--seed=1
--workers=4
--overwrite=FALSE
```

The default full output can live under the configured in-vivo results root, with a compact tracked snapshot under:

```text
Figs/pseudotime_state_pathways/
```

Expected R dependencies are `Seurat`, `Matrix`, `yaml`, `digest`, `edgeR`, `limma`, `splines`, `fgsea`, `msigdbr`, `dplyr`, `readr`, `tidyr`, `ggplot2`, and `patchwork`. Dependency checks should occur before any outputs are written.

### Phase 1: Input and identity audit

1. Read the frozen interval configuration before loading expression.
2. Read the CellCycle metadata and retain finite pseudotime values on the expected scale.
3. Require unique `cell_id` values and a single `sample_id`, dose, and initial-ploidy assignment per cell.
4. Read raw `RNA/counts` and check that values are sparse, non-negative, and integer-like.
5. Match metadata `cell_id` values to expression column names exactly.
6. Write all unmatched metadata and expression identifiers to audit tables.
7. Fail if the metadata-to-expression match rate is below 99%; the expected outcome is 100% for the cells used by 04h.
8. Confirm that every retained sample belongs to only one dose and initial-ploidy group.
9. Record package versions, input paths, file checksums, pseudotime range, cell counts, gene counts, and sample counts.

Do not silently repair barcodes by stripping sample prefixes or `-1` suffixes. If normalization is necessary, make it an explicit, one-to-one mapping and export the crosswalk.

Require pseudotime to remain on the established `[0, 1]` orientation. Values may be checked with a small numerical tolerance, but the workflow must not renormalize them, because doing so would invalidate the frozen interval endpoints.

### Phase 2: Sample-aware pseudobulk construction

Construct fixed, equally spaced pseudotime bins across `[0, 1]`. The default is 20 bins, subject to a preflight coverage report.

For each mouse and pseudotime bin:

1. Sum raw counts across cells.
2. Record the number of contributing cells.
3. Retain sample, dose, initial ploidy, bin boundaries, and bin midpoint.
4. Exclude sample-bin observations with fewer than 10 cells from model fitting, but retain them in the coverage table.
5. Remove genes with insufficient expression using one fixed rule, for example CPM greater than 1 in at least as many sample-bin observations as the number of mice contributing to the primary interval.
6. Apply TMM library normalization before `voom` transformation.

Pseudobulking ensures that cells are not treated as independent biological replicates. Multiple bins from the same mouse remain repeated observations and must be modeled as such.

Write a preflight table showing, for every mouse and every frozen region:

- number of cells;
- number of retained bins;
- library size;
- whether the mouse contributes to the primary contrast;
- the reason for any exclusion.

The primary contrast should proceed only if at least six mice contribute usable observations to the primary interval and both flanks. Coverage must include both untreated and treated mice and both initial-ploidy groups so that nuisance terms are estimable. This is a coverage criterion, not a treatment-expression comparison.

### Phase 3: Continuous common-pseudotime model

Use a scalable sample-aware model, preferably `limma-voom` with splines and within-mouse correlation:

```text
expression ~ spline(pseudotime, df = 5) + dose + initial_ploidy
block = sample_id
```

Implementation details:

1. Evaluate a natural-spline basis at each retained bin midpoint.
2. Include dose as an additive factor (`0`, `30`, and `120 mg/kg`) and initial ploidy as an additive factor (`2N` and `4N`).
3. Estimate repeated-observation correlation with `duplicateCorrelation(sample_id)` and refit with that consensus correlation.
4. Use robust empirical-Bayes moderation.
5. Do not include cluster as a covariate because cluster/state differences along pseudotime are part of the biology to be annotated.
6. Do not include treatment-by-pseudotime, dose-by-pseudotime, or ploidy-by-pseudotime interactions.

The additive dose term prevents a pseudotime-independent treatment expression shift from being mislabeled as the state program. It is a nuisance adjustment, not a treated-versus-control contrast. As a sensitivity analysis, repeat the common-smooth model without the dose term and report rank and NES concordance.

The model should generate fitted expression on a dense, equally spaced pseudotime grid for every retained gene. The grid must not be weighted by the observed number of cells at each pseudotime.

`tradeSeq` may be used as a secondary curve-shape check, but it should not be the primary inferential model here because a standard cell-level `tradeSeq` fit would treat cells as independent and would not preserve mouse-level replication. The primary `limma-voom`/blocking design is chosen specifically to retain continuous pseudotime while keeping the mouse as the experimental unit.

### Phase 4: Gene-ranking statistics

Generate two pre-specified rankings.

#### Primary ranking: accumulated state versus neighboring states

Build a linear contrast of the spline coefficients:

```text
mean fitted expression over [0.30, 0.49]
- 0.5 * mean fitted expression over [0.11, 0.30)
- 0.5 * mean fitted expression over (0.49, 0.68]
```

Rank all tested genes by the moderated t statistic for this contrast. Retain the fitted log-expression difference, standard error, t statistic, nominal P value, and BH FDR for reporting, but do not filter genes by significance before GSEA.

#### Secondary ranking: state specificity over the full trajectory

Calculate:

```text
mean fitted expression over [0.30, 0.49]
minus
mean fitted expression outside [0.30, 0.49]
```

Also record whether each gene's fitted maximum lies inside the primary interval. This ranking tests whether a program is specific to the state rather than merely changing monotonically through it.

For duplicate gene symbols, retain the row with the largest absolute primary t statistic and export the duplicate-resolution table.

### Phase 5: Preranked pathway analysis

Run preranked GSEA on all modeled genes.

Primary collection:

- MSigDB Hallmark (`H`), human.

Optional secondary collections, declared before inspecting results:

- Reactome (`C2:CP:REACTOME`);
- Gene Ontology Biological Process (`C5:GO:BP`).

Use `fgseaMultilevel()` with fixed `minSize`, `maxSize`, seed, and `eps = 0`. Apply BH correction within each gene-set collection and ranking. Export the complete table, not only significant pathways.

For every pathway, report:

- NES and FDR;
- ranking definition;
- gene-set size;
- leading-edge genes;
- mean fitted pathway activity across pseudotime;
- whether the pathway peaks inside the primary interval;
- sensitivity results for the broad and strict-support intervals;
- leave-one-mouse-out direction and rank stability.

Pathway activity curves for visualization should be computed from standardized fitted expression of the pathway's leading-edge genes. These curves annotate the trajectory and must not be split into treated and control curves in the primary figures.

### Phase 6: Robustness analyses

Run the following pre-specified checks without redefining the primary window:

1. Broad interval `[0.11, 0.61]`.
2. Strict-support interval `[0.414, 0.426]`, only if it passes the frozen coverage rule.
3. Primary model without the additive dose nuisance term.
4. Primary model with 15 and 25 pseudotime bins.
5. Primary model with spline degrees of freedom 4 and 6.
6. Leave one mouse out, refit, rerank, and rerun Hallmark GSEA.
7. Fit the same common-pseudotime model separately within untreated and treated cells as a diagnostic only; do not test the two curves against one another. A pathway should be labeled a robust trajectory-state marker when its direction is concordant in both subsets and in the combined primary model.
8. Use `[0.80, 0.98]` as a negative-control state annotation. It should produce a biologically distinct program and must not be described as a treated-cell accumulation state.

Summarize robustness with rank correlation for genes, NES correlation for pathways, sign consistency, and the proportion of leave-one-mouse-out runs in which each top pathway remains in the top 20.

## Required outputs

Use a deterministic directory structure:

```text
04i_pseudotime_state_pathways/
├── 00_manifest/
├── 01_qc/
├── 02_pseudobulk/
├── 03_gene_models/
├── 04_gsea/
├── 05_figures/
├── 06_sensitivity/
└── README.md
```

### Tables

At minimum, write:

```text
00_manifest/analysis_parameters.csv
00_manifest/input_checksums.csv
00_manifest/frozen_interval_definition.csv
01_qc/cell_id_join_audit.csv
01_qc/sample_region_coverage.csv
01_qc/expression_source_audit.csv
02_pseudobulk/sample_bin_metadata.csv
03_gene_models/gene_primary_adjacent_state_contrast.csv
03_gene_models/gene_full_trajectory_state_specificity.csv
03_gene_models/gene_symbol_resolution.csv
04_gsea/hallmark_primary_adjacent_state_gsea.csv
04_gsea/hallmark_full_trajectory_state_specificity_gsea.csv
04_gsea/hallmark_leading_edge_genes.csv
04_gsea/pathway_activity_over_pseudotime.csv
06_sensitivity/pathway_robustness_summary.csv
06_sensitivity/leave_one_mouse_out_hallmark_gsea.csv
```

### Figures

At minimum, generate:

1. `pseudotime_density_with_frozen_state.pdf`: the equal-mouse density difference with the frozen primary interval and sensitivity intervals marked.
2. `primary_state_hallmark_gsea.pdf`: NES and FDR for the top positive and negative Hallmark pathways.
3. `primary_state_pathway_activity_heatmap.pdf`: treatment-agnostic pathway activity over pseudotime with the frozen interval outlined.
4. `primary_state_top_pathway_curves.pdf`: common fitted curves for top pathways, without treatment-specific curves.
5. `primary_state_leading_edge_gene_curves.pdf`: common fitted curves for leading-edge genes.
6. `pathway_robustness_summary.pdf`: NES comparison across frozen sensitivity analyses and leave-one-mouse-out stability.
7. `sample_region_coverage.pdf`: cells and usable pseudobulk bins per mouse in the primary interval and flanks.

Every plot must have a corresponding plotting-data CSV.

The tracked snapshot should contain the interval specification, manifests, QC summaries, gene rankings, GSEA results, plotting-data tables, PDFs, and README. It should not duplicate the raw Seurat object or large pseudobulk count matrices; those remain in the full results directory and are referenced by checksum.

## Interpretation rules

Use the following language consistently:

| Result | Allowed interpretation | Interpretation to avoid |
|---|---|---|
| Positive primary NES | The pathway characterizes the accumulated pseudotime state relative to adjacent states | Treatment upregulates or accelerates the pathway |
| Negative primary NES | The pathway is relatively less active in the accumulated state than in adjacent states | Treatment directly suppresses the pathway |
| Pathway peak inside the interval | The pathway marks the state in which treated cells accumulate | The pathway caused the accumulation |
| Robust treated-cell excess in the interval | Treatment changes the abundance of cells occupying this state | Treatment definitively slows progression through this state |
| State enrichment plus independent kinetic evidence | A delay or bottleneck may be supported | Make no speed claim without the independent evidence |

The primary report should therefore say, for example:

> Gemcitabine-treated tumors contain an increased proportion of cells in pseudotime 0.30-0.49. This state is characterized by positive enrichment of pathways X, Y, and Z relative to the immediately preceding and following trajectory states.

If justified, a separate sentence may state that the accumulation is consistent with a bottleneck at that state, while explicitly noting that cross-sectional pseudotime does not establish transition speed.

## Validation and acceptance criteria

The implementation is complete when all of the following hold:

- The frozen interval file is read-only input to the expression analysis and its checksum is recorded.
- No expression-derived quantity is used to move or resize the primary interval.
- At least 99% of CellCycle metadata cell IDs match the raw-count matrix, or the workflow fails with an audit table.
- The primary expression source is verified raw, non-negative, integer-like RNA counts.
- Every fitted observation is a mouse-by-pseudotime-bin pseudobulk; no single cell is treated as an independent replicate.
- Repeated bins from the same mouse are blocked in the model.
- The primary model contains a common pseudotime smooth and no treatment-by-pseudotime interaction.
- The adjacent-state contrast gives the left and right flanks equal weight and uses an equal-grid rather than cell-density-weighted average.
- GSEA uses a full ranked gene list and reports BH-adjusted results.
- All figures have exact plotting-data tables.
- Results are reproducible from a clean output directory with a fixed seed.
- A smoke test using a small gene and cell subset exercises input validation, pseudobulk construction, spline-contrast construction, and GSEA-table creation.
- A regression test confirms that permuting treatment labels while leaving pseudotime and expression unchanged does not change the unadjusted common-smooth state contrast; for the dose-adjusted model, changes must be limited to the nuisance adjustment and documented.
- The README and final summary use state-enrichment language rather than causal speed claims.

## Suggested implementation order

1. Add and freeze the YAML interval specification.
2. Implement raw-count loading and the exact cell-ID crosswalk.
3. Implement coverage reporting and mouse-by-bin pseudobulk matrices.
4. Implement the common spline model and verify the primary contrast algebra on simulated curves.
5. Add primary and secondary gene rankings.
6. Extract/reuse Hallmark `fgsea` helpers from 03b and add complete result exports.
7. Add pathway and leading-edge fitted-curve summaries.
8. Add figures and their plotting-data tables.
9. Add sensitivity and leave-one-mouse-out runs.
10. Add smoke/regression tests and a figures-only regeneration mode.
11. Run the full workflow, review QC before interpreting pathways, and materialize the compact tracked snapshot.

## Key design decision

The treatment comparison answers **where cells accumulate**. The common pseudotime expression model answers **what biological state is located there**. Keeping these as separate stages prevents a second treatment contrast from changing the biological question and supports the intended conclusion without overstating what cross-sectional pseudotime can establish.
