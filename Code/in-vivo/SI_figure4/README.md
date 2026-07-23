# Supplementary Figure 4

`generate_supplementary_figure4.R` draws the in-vivo tumor-cell UMAP and cluster-composition panels from the paired Figure 7 scVelo-stage outputs:

- `Data/in-vivo/seurat_metadata.csv`
- `Data/in-vivo/scvelo_cell_metrics.csv`

The canonical Manager module is `si_figure4`. It writes an immutable source run to `Results/in-vivo/SI_figure4/runs/<run_id>_si_figure4/` before the `in_vivo_figure7` module. Input resolution is pairwise: Manager first uses the two published `Data/in-vivo` files, then the two files in the configured Figure 7 intermediate directory, and otherwise runs the Figure 7 scVelo input-preparation stage. Both intermediate files are validated, staged, and checksum-verified before publication and plotting.

The source run contains seven PDF/PNG pairs: panels `SuppFig4A` through `SuppFig4F` plus the full composite. Their canonical filenames begin with `panel_SuppFig4`; Manager publishes them without renaming to `figures/Supplementary/` and records them in `figures/Supplementary/manifest.tsv`.

Tables and run metadata remain under the immutable Results run and are not copied into `figures/`.
