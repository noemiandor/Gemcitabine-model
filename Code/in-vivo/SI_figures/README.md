# Supplementary Figures 4-7

This directory contains the Supplementary Figures 4-7 orchestrator, a
plot-only renderer, and the narrowly scoped raw-table builder. Routine mode
uses the 11 reviewed plot-facing tables under `Data/in-vivo/SIfigures/` and
writes four final composite figures as PDF/PNG pairs.

Run it directly with:

```bash
Rscript Code/in-vivo/SI_figures/run_supplementary_figures.R \
  --mode=plot-only \
  --output-dir Results/in-vivo/SI_figures/runs/example_si_figures
```

Or run the same frozen-input workflow through Manager:

```bash
bash Manager.sh \
  --modules si_figures \
  --run-id example
```

The renderer intentionally cannot download raw data, load the Seurat RDS,
perform differential expression, or rerun ORA/GSEA. Its only analysis
dependencies are `ggplot2`, `patchwork`, `pheatmap`, and `yaml`.

To rebuild Supplementary Figures 4-7 from the shared Seurat source boundary
without running Figures 1-6:

```bash
bash Manager.sh \
  --mode full-refit \
  --modules si_figures \
  --run-id example_raw_si
```

`run_supplementary_figures.R --mode=full-workflow` validates the 11-table cache
final Seurat object plus the versioned endpoint-ploidy table. The Seurat object
can be reused/reconstructed from the same five-stage Cell Ranger H5 cache as
Figure 7 by passing `--figure7-cellranger-root` through Manager; both modules
use `--figure7-seurat-upstream-dir`. If that boundary is unavailable, the
checksum-pinned deposited final RDS is the fallback. SI Figures 4-7 do not
require scVelo; an explicitly supplied scVelo table is audit-only.

A previously generated 11-table cache is reused before selecting or opening a
Seurat source. Its manifest binds the exact source RDS, transitive upstream
dependencies, cumulative scientific code contracts, scoped SI configuration,
and stage-specific runtime contract. Missing archived ancestors are
manifest-attested; any still-present final, partial-stage, or H5 source must
match. If more than one cache is lineage-compatible, the run stops unless one
is selected explicitly with `--generated-cache-dir`.

The builder writes exactly the same 11 plot-facing table names to a run-scoped
cache. Cluster differential-expression work is separately resumable: each
completed cluster has a table and dependency sidecar under the cache's stable
`work/` directory, and is reused only after schema, hash, and fingerprint
validation. ORA and GSEA working data remain below that directory. Only the
four composite PDF/PNG pairs are materialized under
`figures/Supplementary/`.

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

The raw fallback currently preserves that previous Tao behavior because the
requested GRCh/GRCm repair has been explicitly deferred until the entire raw
pipeline is in place. Raw-rebuilt SI7 matrices and composites are therefore
marked `canonical_publication_allowed=false`; they cannot overwrite or validate
as the canonical human-only cache.

Individual subpanels are constructed in memory as part of each composite and
are not published as duplicate derivatives.
