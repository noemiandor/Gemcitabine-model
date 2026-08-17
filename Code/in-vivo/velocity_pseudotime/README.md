# Supplementary RNA-velocity/pseudotime panel

This module produces the standalone Supplementary Figure 10 panel needed to
support the manuscript statement that 2,881 xenograft-derived cells in
clusters 4c, 6, and 10 were ordered by stochastic scVelo velocity pseudotime
with cluster 6 as the root.

The panel contains all four parts of that claim:

- the reviewed Seurat UMAP positions for the exact 2,881-cell tumor
  `CellCycle` subset;
- a binned RNA-velocity vector field derived from
  `scvelo.tl.velocity_embedding(..., basis="umap")` after the stochastic
  velocity graph is calculated;
- cell color on a fixed 0--1 velocity-pseudotime scale; and
- labels for clusters 4c, 6, and 10, with every cluster-6 cell declared as a
  root cell and cluster 6 visibly labeled and outlined as the root set.

Cluster 6 is a pre-specified trajectory-root choice inherited from the
reviewed analysis; it is not selected or optimized by this plotting module.
The red dashed outline is a display device only and uses the 90% radial core
of cluster 6. It does not redefine which cells are roots: the table contract
requires `is_root = TRUE` for every and only cluster-6 cell.

## Cache-first execution

Once reviewed, the small plot-facing table and checksum-bound provenance live
at:

```text
Data/in-vivo/figure7/processed/cellcycle_velocity_pseudotime_umap.tsv
Data/in-vivo/figure7/processed/cellcycle_velocity_pseudotime_umap.provenance.tsv
```

A routine run reads only that bundle:

```bash
bash Manager.sh \
  --modules in_vivo_velocity_pseudotime \
  --run-id velocity_panel
```

The module is intentionally opt-in until that frozen table has been generated
from the raw source boundary and scientifically reviewed. After review, it can
be added to Manager's default module list without introducing any raw scVelo
dependency into ordinary figure regeneration.

Full-workflow mode first reuses a valid frozen bundle. It may reuse a
run-scoped generated bundle only after matching its recorded hashes to the
current CellCycle/SI inputs, builder, plotting helper, and scVelo generator.
The corresponding embedding cache also requires a sidecar that binds the
embedding hash to the selected final Seurat RDS, an exact 18-file loom
inventory (bytes and SHA-256 for every file), the Figure 7 config, and the
environment lock. An embedding TSV without those sidecars is never reused.
If any dependency differs, the affected stage is regenerated:

```bash
bash Manager.sh \
  --mode full-refit \
  --modules in_vivo_velocity_pseudotime \
  --run-id velocity_panel_raw \
  --figure7-python /path/to/python-with-pinned-scvelo \
  --figure7-raw-seurat-rds /path/to/integrated_sct_cca_seurat_final_reclustered.rds \
  --figure7-loom-root /path/to/velocyto_loom
```

The source stage computes velocity and pseudotime on the complete reviewed
35,513-cell Seurat universe using the deposited loom layers, reviewed Seurat
PCA/UMAP, 30 neighbors, 2,000 genes, and stochastic scVelo. It exports cell
IDs, UMAP positions, UMAP velocity components, velocity pseudotime, cluster,
context, sample, and root flags. The table builder then matches the exact
2,881-cell tumor `CellCycle` membership against the reviewed processed table
and SI metadata; mismatched pseudotime, UMAP, cluster, sample, context, or root
assignments are fatal.

Raw-generated output is marked noncanonical and is not materialized into the
manuscript-facing asset tree until its compact table has been reviewed and
frozen. After that promotion, Manager materializes only the final PDF and PNG
as:

```text
figures/Supplementary/panel_SuppFig10_velocity_pseudotime.pdf
figures/Supplementary/panel_SuppFig10_velocity_pseudotime.png
```

Run metadata records repository-relative locators, never absolute HPC or
workstation paths. Both Manager and the asset materializer independently
require the exact canonical table/provenance locators, root/count/vector
contract, checksums, and byte-identical run plot table. A generated, external,
or alternate `--frozen-table` bundle remains useful for review but cannot be
published by changing a flag in its run metadata.

No container, vendored package archive, or package snapshot is introduced.

Run the focused tests with:

```bash
Rscript Code/in-vivo/velocity_pseudotime/tests/test_velocity_pseudotime_panel.R
python3 -m unittest Code.tools.tests.test_manager_velocity_pseudotime_cli
```
