# Supplementary Figures 4-7

`generate_supplementary_figures.R` generates Supplementary Figures 4-7. By
default it first validates and uses the complete 32-file table cache under
`Data/in-vivo/SIfigures/`; in that plot-only mode it does not open the raw
Seurat RDS and does not rerun differential expression, ORA, or GSEA.

When the cache is absent or `--force-reanalysis` is supplied, it regenerates
the tables from the Figure 7 input bundle:

- `Data/in-vivo/seurat_metadata.csv`
- `Data/in-vivo/scvelo_cell_metrics.csv`
- `Data/in-vivo/all_ploidy.tsv`
- the raw `integrated_sct_cca_seurat_final_reclustered.rds`
- optional `Data/in-vivo/seurat_metadata_provenance.tsv`
- `Code/in-vivo/figure7/figure7_config.yaml`

The canonical Manager module is `si_figures`. It writes an immutable source run
to `Results/in-vivo/SI_figures/runs/<run_id>_si_figures/` after Figure 7
prerequisites are available. Manager validates the exact 32-table inventory,
publishes the generated PDF/PNG panels to `figures/Supplementary/`, and
atomically publishes those tables plus `manifest.tsv` to
`Data/in-vivo/SIfigures/`. An existing canonical cache is backed up into the
Manager run before replacement.

The source run contains:

- Supplementary Figure 4: panels A-H and composite, comparing Tumor and CellLine cells. The first row is A-D, followed by E-F and G-H.
- Supplementary Figure 5: panels A-I and composite, describing initial ploidy, dose, mouse, and cluster composition.
- Supplementary Figure 6: panels A-D and composite, describing endpoint ploidy from `all_ploidy.tsv`.
- Supplementary Figure 7: panels A-B and composite, reproducing the cluster Hallmark ORA and GSEA heatmaps from a standalone raw-Seurat-RDS workflow.

The immutable source run also includes all composition counts/proportions, the endpoint-ploidy join audit, per-cluster DEG/ORA/GSEA tables, heatmap matrices, exact panel contract, input hashes, session information, and run provenance.

SI Figure 7 differential expression uses a resource-aware two-level future
layout. It detects the CPUs available to the current process, reserves two CPUs
when possible, and also detects the effective Slurm/cgroup memory limit. The
loaded Seurat object size, a measured 1.5-fold per-worker allowance, and 90% of
the memory allocation define a conservative worker budget, preventing CPU-rich but memory-limited jobs from launching one
large Seurat task per cluster. Up to one outer worker per cluster is launched,
and the safe worker budget is distributed across active clusters by quotient
and remainder, with at least one worker per active cluster. Linux multicore
workers share the read-only Seurat object through fork semantics; unsupported
environments fall back to sequential execution. CPU, memory, object-size,
outer-cluster, per-cluster, and globals-limit values are recorded in
`metadata/run_config.tsv`.

The strict materializer reconciles the exact 27-logical-panel/54-file inventory
and manifests before publication. Tables and run metadata remain under the
immutable Results run as well as the validated table cache; tables are not
copied into `figures/`.

The Seurat RDS and loom inputs are themselves processed artifacts with a series of preceding data-generation and analysis steps. Their reviewed workflow description, Seurat provenance, loom-generation parameters, runtime environments, package inventories, session information, and checksums are archived at [Zenodo DOI 10.5281/zenodo.21463392](https://doi.org/10.5281/zenodo.21463392). After every successful Manager `si_figures` run, a uniform offline index is published under `metadata/upstream_analysis/` as `README.md`, `zenodo_document_manifest.tsv`, and `provenance.tsv`; the index is also included in the standard input/output manifests.
