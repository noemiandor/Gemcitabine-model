# Figure 7 Reproducibility Pipeline

## Status and authority

Figure 7 is implemented as the scoped `in_vivo_figure7` manuscript module.
The executable contract is defined by `Manager.sh`,
`Code/in-vivo/figure7/figure7_config.yaml`, the module README, and the test
suites. This document summarizes the current architecture and promotion rules.

## Scientific and display contract

The module produces six scientific source panels, 7A--7F, as matched PDF and
PNG files. It also assembles the manuscript-facing A--L composite using the
configured first-citation mapping:

```text
A=7A, B=7C, C=SI4A, D=SI4B, E=SI4C, F=SI4E,
G=SI7B, H=7B-ECDF, I=7F, J=QC-filtered copy-number heatmap,
K=7D, L=7E
```

The shared supplementary plotting implementation remains the single source of
truth for SI4A--C/E and SI7B; their copies in the supplement and main Figure 7
must be regenerated from the same reviewed tables.

The source-panel analyses are:

- 7A: tumor-growth trajectories and the Day-17 TGI calculation;
- 7B: selected CellCycle pseudotime-distribution comparisons;
- 7C: Day-17 TGI by injected initial ploidy;
- 7D: within-dose-centered TGI versus the injected-origin-matched ECDF shift;
- 7E: Day-17 TGI versus mean endpoint tumor-cell ploidy, using the raw
  mouse-level Pearson correlation and exact unrestricted enumeration of all
  8! TGI-label permutations;
- 7F: pathway activity across the directly computed 0.296--0.486 CellCycle
  pseudotime interval, using the reviewed GRCh-only,
  injected-initial-ploidy-adjusted pointwise-v4 reference.

Panel 7F retains only pathways with collection-wide BH-adjusted P <= 0.05,
selects up to four pathways per direction and collection, and never fills a
direction with nonsignificant pathways. Endpoint tumor-cell ploidy is not a
nuisance covariate in that model.

## Analysis modes

### Routine mode

Routine execution recomputes 7A--7E from the tracked, plot-facing cell tables
and renders 7F from the checksum-pinned reviewed pointwise-v4 compact reference. It also
validates the reviewed shared SI cache and the exact 9,832-cell copy-number
inputs before assembling the A--L composite. The density-localization component
of source 7B is regenerated as Supplementary Figure 4I.
This mode is intentionally lightweight and does not require raw loom files or
a large Seurat object.

### Full refit

`full-refit` is the cache-first raw-data fallback. When Figure 7 requires the
shared SI panels, Manager schedules `si_figures` before `in_vivo_figure7`.
Valid downstream intermediates are reused; missing stages are reconstructed
from the earliest available reviewed boundary:

- a deposited final Seurat RDS or 18 Cell Ranger filtered-feature H5 files;
- 18 loom files for pseudotime/velocity-derived cell metrics;
- the pinned human gene-set release;
- the sample-information and tumor-growth workbooks; and
- the endpoint-ploidy inventory and its 16 reviewed CBS matrices.

Before expression filtering, symbol resolution, model fitting, or GSEA, the
state-pathway workflow retains exact `GRCh38-` count rows, excludes exact
`GRCm39-` rows, and rejects unrecognized prefixes. It computes the unique
connected positive pointwise-supported interval before loading expression
counts and uses that exact interval and its derived equal-width flanks for the
model and GSEA. A raw refit writes a versioned generated pointwise-v4 candidate and marks it noncanonical until its exact
scientific outputs are reviewed and promoted.

### Render-only and panels-only

The standalone `render-only` mode regenerates panels from a concrete immutable
Figure 7 result run without rerunning statistics. Manager `panels-only` does
not invoke analysis or rendering; it materializes an explicit source run
selected by `--source-run-id`. Neither mode may infer a mutable `latest`
pointer or publish a noncanonical reference.

## Input and cell-universe contract

Routine 7A--7E inputs are:

```text
Data/in-vivo/figure7/processed/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv
Data/in-vivo/figure7/processed/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv
Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv
```

The CellCycle and NonCellCycle tables must represent the exact reviewed final
Seurat universe: 35,513 cells in nine retained clusters, including 9,832 tumor
cells and exactly 5,335 treated tumor cells. Clusters excluded by the reviewed
QC policy may not re-enter any Figure 7 or SI4--7 analysis through an upstream
fallback or cache.

The endpoint-ploidy file contains the complete 14,125-cell, 16-file CBS
inventory. Panel 7E scores the exact QC-passed plot-table tumor-cell union and
requires every file/barcode/value tuple to match the complete inventory. The
eight treated tumors contribute 5,335 scored cells.

The reviewed panel-7F reference consists of these nine compact files under a
versioned reference directory:

```text
state_pathway_interval_definition.tsv
panel_7F_pathway_activity_plot_data.tsv
panel_7F_selected_pathway_gsea.tsv
panel_7F_leading_edge_genes.tsv
state_pathway_gene_ranking_complete.tsv
state_pathway_gsea_complete.tsv
state_pathway_sample_bin_coverage.tsv
state_pathway_design_qc.tsv
state_pathway_provenance.tsv
```

The renderer uses the stored collection and pathway display order and does not
rerank at render time. The complete ranking, GSEA, coverage, design, and
provenance tables remain mandatory so the displayed selection can be audited
without the raw Seurat object.

## Output and provenance contract

Each analysis run has `figures/`, `tables/`, `metadata/`, and `logs/`
directories. Required metadata include portable input/output manifests, the
run configuration, exact panel contract, session information, reference
identity, source hashes, and shared-SI cache lineage.

Publication materialization requires all of the following:

1. the exact reviewed pointwise-v4 reference identity and checksums;
2. canonical-publication approval in both run metadata and reference
   provenance;
3. the reviewed shared-SI cache identity and manifest;
4. the exact A--L panel mapping and output inventory;
5. the reviewed Figure 7L numerical/statistical contract;
6. the exact 9,832-cell Figure 7J copy-number contract; and
7. complete, relocatable input and output manifests.

Any missing, duplicated, unexpected, noncanonical, or checksum-mismatched
asset causes materialization to fail before manuscript-facing files are
written.

## Manager usage

Routine Figure 7 only:

```bash
bash Manager.sh \
  --mode standard \
  --modules in_vivo_figure7 \
  --run-id <run_id>
```

Cache-first end-to-end Figure 7 and SI4--7 regeneration:

```bash
bash Manager.sh \
  --mode full-refit \
  --modules in_vivo_figure7,si_figures \
  --run-id <run_id>
```

Materialization from an immutable reviewed run:

```bash
bash Manager.sh \
  --mode panels-only \
  --modules in_vivo_figure7 \
  --source-run-id <analysis_run_id> \
  --run-id <materialization_operation_id>
```

## Validation and promotion

Changes are ready for promotion only when:

- the R contract and workflow tests pass;
- the Manager CLI and materializer tests pass;
- standard execution reproduces the reviewed source panels and A--L
  composite;
- full-refit either reuses a valid lineage or reconstructs every missing stage
  from pinned inputs;
- the exact reviewed cell universe and human-feature boundary remain enforced;
- all manifests validate in a relocated checkout; and
- visual review confirms readable, correctly ordered, unclipped panels at the
  target manuscript dimensions.

Only the scoped Figure 7/SI4--7 code, compact scientific inputs, reusable
intermediate contracts, targeted tests, and reviewed manuscript assets belong
in this workflow.
