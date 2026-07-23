# Supplementary Figure 4

`generate_supplementary_figure4.R` draws the in-vivo tumor-cell UMAP and cluster-composition panels from the Figure 7 scVelo-stage input bundle and its shared configuration:

- `Data/in-vivo/seurat_metadata.csv`
- `Data/in-vivo/scvelo_cell_metrics.csv`
- `Data/in-vivo/seurat_metadata_provenance.tsv`
- `Code/in-vivo/figure7/figure7_config.yaml`

The canonical Manager module is `si_figure4`. It writes an immutable source run to `Results/in-vivo/SI_figure4/runs/<run_id>_si_figure4/` before the `in_vivo_figure7` module. Input resolution is bundle-based: Manager first uses the three published `Data/in-vivo` files, then the complete three-file bundle in the configured Figure 7 intermediate directory, and otherwise runs the Figure 7 scVelo input-preparation stage. Partial bundles are rejected. All three files are schema-validated, staged, checksum-verified, and atomically published before plotting.

The source run contains seven PDF/PNG pairs: panels `SuppFig4A` through `SuppFig4F` plus the full composite. Their canonical filenames begin with `panel_SuppFig4`; Manager publishes them without renaming to `figures/Supplementary/` and records them in `figures/Supplementary/manifest.tsv`.

The immutable source run also contains the formal reproducibility artifacts:

- `tables/si_figure4_cell_metadata.csv`: one row per source-object cell with cell ID, UMAP coordinates, sample, cluster ID/annotation/order/color, initial ploidy, treatment/dose, CellCycle classification, and explicit inclusion/exclusion fields.
- `tables/si_figure4_cluster_key.tsv`: cluster annotation, classification, continuous display order, color, and all/included cell counts.
- `tables/cluster_composition_by_mouse.csv`: cluster counts, mouse denominators, proportions, and denominator definition.
- `tables/cluster_composition_by_ploidy_dose.csv`: mouse-weighted group summaries plus contributing mouse and cell counts.
- `metadata/si_figure4_provenance.tsv`: source Seurat RDS checksum, reductions, clustering resolution/fields, QC and inclusion policy, artifact checksums, and row counts.

The strict materializer reconciles these tables against one another and against their manifests before it publishes any SI Figure4 panel. Tables and run metadata remain under the immutable Results run and are not copied into `figures/`.

The Seurat RDS and loom inputs are themselves processed artifacts with a series of preceding data-generation and analysis steps. Their reviewed workflow description, Seurat provenance, loom-generation parameters, runtime environments, package inventories, session information, and checksums are archived at [Zenodo DOI 10.5281/zenodo.21463392](https://doi.org/10.5281/zenodo.21463392). After every successful Manager `si_figure4` run, a uniform offline index is published under `metadata/upstream_analysis/` as `README.md`, `zenodo_document_manifest.tsv`, and `provenance.tsv`; the index is also included in the standard input/output manifests.
