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

Composition panels SI4E/G and SI5F-I share
`normalized_composition.R`. Cells are first converted to within-sample
cluster proportions; samples are then averaged with equal weight inside each
displayed group, so sequencing depth cannot determine a sample's influence.
Positive enrichment is assessed on those independent-sample proportions by
an exact one-group-versus-rest label-permutation test inside the relevant
nuisance strata: initial ploidy for SI4E and SI5G/H, context for SI4G, and dose
for SI5I. Benjamini-Hochberg correction covers every group-by-cluster contrast
in a panel, and asterisks mark only positive enrichments at FDR <= 0.05. SI5F
has one biological sample per displayed mouse; it is therefore explicitly
descriptive and emits no inferential stars rather than treating cells as
replicates.

Main Figure 7J is reproduced from the 16 tracked downstream
`Data/in-vivo/scRNAseq_Numbat/*.sps.cbs` matrices. The matrices contain 14,125
cells and use two different segment schemas. Their exact filenames, byte sizes,
and SHA-256 values are pinned in `scRNAseq_Numbat/cbs_manifest.tsv` and checked
before rendering. The complete matrix collection is a validated source, not
the analysis universe: the frozen endpoint-ploidy audit restricts Figure 7J and SI6E to
the exact 9,832 tumor cells retained after final Seurat QC (including 5,335
treated cells), so discarded clusters 3, 4, 9, and 9c cannot re-enter through
the downstream CBS files. Because those schemas cannot be
assumed to share a coordinate build, the renderer does not align breakpoints
or project them to common loci. Instead, it computes each cell's
length-weighted mean across the available CBS segments of each autosome and
plots the resulting 9,832-by-22 matrix in chromosome order. Rows are grouped
by injected 2N/4N origin, dose, and mouse, then ordered by the post-processed
copy-number score; neither rows nor columns are clustered. Mouse, gemcitabine
dose, and injected origin are displayed as row annotations. The renderer exports the plotted
matrix, cell order, and per-file/per-schema represented-base-pair audit.

SI6E compares the project-designated injected-cell karyotype references (20 2N
A7M and 16 4N A5M metaphases) with one terminal NUMBAT-derived mean per mouse;
dose is encoded by color. For the karyotype references, the autosomal
length-weighted estimate is added to the `chr999` value, which is expressed in
haploid-genome-equivalent units of unassigned DNA. This converts the intermediate
assigned-autosomal means of 2.00997 and 3.51561 to final 2N- and 4N-reference
means of 2.15124 and 3.94651. The terminal mouse-balanced means are 2.13551 and
2.31855, respectively: descriptive changes of -0.01573 (-0.73%) for 2N and
-1.62796 (-41.25%) for 4N. The 4N-minus-2N separation contracts from 1.79527 in
the references to 0.18304 at endpoint (89.80%). Every terminal 4N-origin cell
estimate is below the minimum 4N-reference metaphase. The A7M and A5M matrices
are project-designated lineage-matched proxies, not the same-passage A6M and A4M
inocula. The cross-assay comparison is therefore descriptive and has no P value:
each proxy represents one culture-level biological unit, and the independently
processed endpoint runs use different schemas/calibration. Metadata retain the
exact reference-cell values, mouse means, ranges, changes, source commit,
confirmed `chr999` unit interpretation, input hashes, and this inference boundary.

The tracked CBS matrices reproduce this downstream panel and recompute and
validate every score/coverage value in the canonical
`Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv` (and its reduced
`Data/in-vivo/all_ploidy.tsv` projection). The two injected-reference matrices
and their checksums are under `scRNAseq_Numbat/injected_reference/`. These files do
not constitute a complete NUMBAT run: allele-count inputs, clone posteriors,
consensus segment outputs preceding these matrices, phylogeny, configuration,
logs, and an executable upstream inference workflow remain unavailable.

SI4A-C/E and SI7A/B are constructed by `shared_context_panels.R`. SI4I
recomputes the shared density-localization analysis from the exact tracked
2,881-cell pseudotime table. Main Figure 7 calls the same helper for its copies
of SI4A-C/E and SI7B, with display tags assigned by the A-L compositor; fixed shuffle keys keep the UMAP point order
identical between main and supplementary copies.

To rebuild Supplementary Figures 4-7 from the shared Seurat source boundary
without running Figures 1-6:

```bash
bash Manager.sh \
  --mode full-refit \
  --modules si_figures \
  --run-id example_raw_si
```

For a full-refit A-L Figure 7, Manager automatically runs this module before
`in_vivo_figure7` (including when only `in_vivo_figure7` was requested). It
passes the run-scoped table cache, analysis-input manifest, run configuration,
and SI provenance to Figure 7. That generated-cache composite is a
noncanonical review candidate; canonical routine Figure 7 remains bound to the
exact reviewed cache below `Data/in-vivo/SIfigures/`.

`run_supplementary_figures.R --mode=full-workflow` deliberately bypasses the
reviewed plot-only cache: it first reuses a lineage-valid generated human-only
cache when one exists, otherwise it rebuilds the tables from the final Seurat
object plus the versioned endpoint-ploidy table. The Seurat object can be
reused/reconstructed from the same five-stage
Cell Ranger H5 cache as Figure 7 by passing `--figure7-cellranger-root` through
Manager; both modules use `--figure7-seurat-upstream-dir`. If that boundary is
unavailable, the checksum-pinned deposited final RDS is the fallback. SI
Figures 4-7 do not require scVelo; an explicitly supplied scVelo table is
audit-only.

A previously generated 11-table cache is reused before selecting or opening a
Seurat source. Its manifest binds the exact source RDS, transitive upstream
dependencies, cumulative scientific code contracts, scoped SI configuration,
and stage-specific runtime contract. Missing archived ancestors are
manifest-attested; any still-present final, partial-stage, or H5 source must
match. If more than one cache is lineage-compatible, the run stops unless one
is selected explicitly with `--generated-cache-dir`.

The builder writes exactly the same 11 frozen support-table names to a run-scoped
cache. Cluster differential-expression work is separately resumable: each
completed cluster has a table and dependency sidecar under the cache's stable
`work/` directory, and is reused only after schema, hash, and fingerprint
validation. ORA and GSEA working data remain below that directory. Only the
four composite PDF/PNG pairs are materialized under
`figures/Supplementary/`.

The cache contains:

- canonical cell metadata and cluster keys used by Figures 4-6;
- six composition/audit tables retained by the established cache contract;
- the two 20-by-9 Hallmark matrices plotted in Figure 7.

The normalized SI4E/G and SI5F-I estimates and tests are derived directly from
canonical cell metadata. The SI4 context/ploidy tables still supply the raw
count panels SI4F/H; the four historical SI5 aggregate tables are retained and
validated for backward-compatible audit only and are not plotting inputs.

`Code/tools/validate_si_figures_table_cache.py` enforces the exact 11-file
inventory, schemas, cell/cluster reconciliation, endpoint-ploidy contract,
finite matrix values, and portable SHA-256 manifest.

The renderer writes the exact composition estimates and test results to
`metadata/normalized_composition_plot_data.tsv` and
`metadata/normalized_composition_enrichment_tests.tsv`. Run the targeted
statistical and SI7 heatmap contract tests with:

```bash
Rscript Code/in-vivo/SI_figures/tests/test_normalized_composition.R
Rscript Code/in-vivo/SI_figures/tests/test_si7_heatmap_clustering.R
Rscript Code/in-vivo/SI_figures/tests/test_shared_context_panels.R
```

## SI Figure 7 species policy

The cells in this analysis are human tumor/cell-line cells aligned to a
combined human/mouse reference. The reviewed SI Figure 7 matrices therefore
use a human-only feature policy: retain exact `GRCh38-` features, exclude
`GRCm39-` features, and reject unclassified features before fresh RNA
normalization, differential expression, symbol cleanup, duplicate resolution,
ORA, and GSEA. The frozen matrices are the exact scientifically reviewed
outputs from raw-refit run
`grch_human_only_v2_20260729_raw_refit_retry3_si_figures`, using MSigDB
2026.1.Hs Hallmark gene sets. The other nine SI4-7 frozen tables remain
unchanged.

Both SI7 heatmaps hierarchically cluster the pathway rows and cluster columns;
the rendered panels therefore include dendrograms on both axes. Clustering
changes display order only and does not modify the reviewed matrix values.

This explicit filter replaces the previous mixed behavior, where species
prefixes were stripped but symbol case was preserved before querying human
gene sets. In ORA, uppercased duplicate keys allowed human/mouse symbol pairs
to compete for marker selection. In GSEA, case-distinct mouse rows generally
could not match the human Hallmark symbols but still occupied positions in the
ranked vector and changed the enrichment statistic.

The raw fallback now builds a new SI7-only Seurat object from exact `GRCh38-`
RNA counts and runs fresh `LogNormalize` before differential expression, ORA,
and GSEA. It does not reuse the mixed-species normalized data layer. Species
counts, the policy/helper hashes, and the human-only normalization contract are
recorded in the run metadata. The exact approved retry3 SI7 matrices are now
the reviewed routine cache. Future raw-rebuilt matrices and composites remain
marked `canonical_publication_allowed=false`; they cannot overwrite or
impersonate the reviewed cache without a new explicit promotion.

Individual subpanels are constructed in memory as part of each composite and
are not published as duplicate derivatives.
