# Supplementary Figures 4-7

This directory contains the plot-only generator for Supplementary Figures 4-7.
It uses the 11 reviewed, plot-facing tables under
`Data/in-vivo/SIfigures/` and writes four final composite figures as PDF/PNG
pairs.

Run it directly with:

```bash
Rscript Code/in-vivo/SI_figures/generate_supplementary_figures.R \
  --output-dir Results/in-vivo/SI_figures/runs/example_si_figures
```

Or run the same frozen-input workflow through Manager:

```bash
bash Manager.sh \
  --modules si_figures \
  --run-id example
```

The generator intentionally cannot download raw data, load the Seurat RDS,
perform differential expression, or rerun ORA/GSEA. Its only analysis
dependencies are `ggplot2`, `patchwork`, `pheatmap`, and `yaml`.

The cache contains:

- canonical cell metadata and cluster keys used by Figures 4-6;
- six plot-facing composition/audit tables;
- the two 20-by-9 Hallmark matrices plotted in Figure 7.

`Code/tools/validate_si_figures_table_cache.py` enforces the exact 11-file
inventory, schemas, cell/cluster reconciliation, endpoint-ploidy contract,
finite matrix values, and portable SHA-256 manifest.

## SI Figure 7 species policy

The cells in this analysis are human tumor/cell-line cells aligned to a
combined human/mouse reference. The reviewed SI Figure 7 matrices therefore
use a human-only feature policy: retain `GRCh*` features and remove `GRCm39-*`
features before symbol cleanup, duplicate resolution, ORA, and GSEA. The
frozen matrices were recalculated from Tao's cluster DEG cache using MSigDB
2026.1.Hs Hallmark gene sets.

This explicit filter replaces the previous mixed behavior, where species
prefixes were stripped but symbol case was preserved before querying human
gene sets. In ORA, uppercased duplicate keys allowed human/mouse symbol pairs
to compete for marker selection. In GSEA, case-distinct mouse rows generally
could not match the human Hallmark symbols but still occupied positions in the
ranked vector and changed the enrichment statistic.

Only the final composite files are materialized under
`figures/Supplementary/`. Individual subpanels are constructed in memory as
part of each composite and are not published as duplicate derivatives.
