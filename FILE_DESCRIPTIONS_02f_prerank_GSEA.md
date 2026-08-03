# 02f Pre-ranked GSEA: File and Figure Guide

This document explains the contents of the `02f_prerank_GSEA` result directory produced by `Code/in-vivo/02f_prerank_GSEA.R`.

The goal of this step is to:
- reuse per-cluster DEG tables from `02e_cluster_annotation`,
- build a ranked gene list for each `cluster_final`,
- run pre-ranked Hallmark GSEA,
- run pre-ranked KEGG GSEA,
- export per-cluster result tables,
- generate summary tables and cross-cluster heatmaps.

The workflow prefers previously generated DEG files from:

- `02e_cluster_annotation/01_cluster_final_vs_rest_DEG`

If a marker file is missing, the script recomputes only the missing `cluster_final` DEG table with `Seurat::FindMarkers`.

## Directory Structure

- `00_summary/`
  Combined result tables, per-cluster summary table, ranking tables, and the plain-text run summary.
- `01_hallmark_prerank_GSEA/`
  Per-cluster Hallmark GSEA tables and top-pathway plots.
- `02_kegg_prerank_GSEA/`
  Per-cluster KEGG GSEA tables and top-pathway plots.
- `03_plots/`
  Cross-cluster NES heatmaps for Hallmark and KEGG results.

## Gene Set Sources

### Hallmark

Hallmark gene sets are retrieved through `msigdbr` from:

- `MSigDB H`

### KEGG

KEGG-like pathway gene sets are retrieved through `msigdbr` from the first available source among:

- `MSigDB C2:CP:KEGG_MEDICUS`
- `MSigDB C2:CP:KEGG_LEGACY`
- `MSigDB C2:CP:KEGG`

The exact collection actually used in a given run is recorded in:

- `00_summary/run_summary.txt`
- the `collection_source` column of the output result tables

## Ranking Strategy and Key Calculations

### Input DEG Tables

The script reads per-cluster marker tables from `02e_cluster_annotation/01_cluster_final_vs_rest_DEG/cluster_<ID>/cluster_<ID>_vs_rest_markers.csv`.

The expected DEG columns are:

- a recognized logFC column such as `avg_log2FC`,
- `p_val` if available, otherwise `p_val_adj`,
- `gene` and optionally `gene_symbol`.

### Gene Normalization and Deduplication

Before ranking:

- prefixes such as `GRCh38-`, `GRCm39-`, and `hg38-` are removed,
- trailing version suffixes are removed,
- symbols are uppercased for matching,
- duplicated genes are collapsed by retaining the row with the largest absolute ranking statistic,
- ties are further broken by larger absolute logFC and then smaller p value.

### Ranking Metric

The pre-ranked statistic is:

`rank_stat = sign(logFC) * -log10(p_value)`

with a tiny logFC-based tie-break term added internally to keep stable ordering when the main statistic ties.

Interpretation:

- positive values indicate genes enriched in the focal cluster,
- negative values indicate genes depleted in the focal cluster relative to the rest,
- larger absolute values indicate stronger evidence.

### Pre-ranked GSEA

The script runs `fgsea::fgseaMultilevel` using:

- the ranked vector defined above,
- minimum gene-set size: 15,
- maximum gene-set size: 500,
- exact epsilon setting: `eps = 0`.

The current summary threshold for significance is:

- `padj < 0.05`

## 00_summary

### `Hallmark_prerank_GSEA_all.csv`

Combined Hallmark GSEA results across all clusters.

Columns:

- `cluster`: `cluster_final` label
- `collection`: collection label, here `Hallmark`
- `collection_source`: the exact gene-set source used in the run
- `pathway`: original pathway ID from MSigDB
- `pathway_label`: cleaned human-readable pathway name used in figures
- `size`: number of ranked genes from this pathway that were present in the analysis
- `ES`: raw enrichment score from GSEA
- `NES`: normalized enrichment score
- `pval`: raw fgsea p value
- `padj`: Benjamini-Hochberg adjusted p value
- `log2err`: fgsea uncertainty estimate when provided by the installed fgsea version
- `direction`: `up` if `NES >= 0`, `down` if `NES < 0`
- `leadingEdge_genes`: semicolon-separated list of leading-edge genes reported by fgsea

### `KEGG_prerank_GSEA_all.csv`

Combined KEGG GSEA results across all clusters.

Column definitions are the same as in `Hallmark_prerank_GSEA_all.csv`.

### `cluster_prerank_GSEA_summary.csv`

One-row-per-cluster summary of ranked-gene counts and leading significant pathways.

Columns:

- `cluster`: `cluster_final` label
- `n_ranked_genes`: number of genes retained after cleaning and deduplication
- `n_hallmark_sig`: number of Hallmark pathways with `padj < 0.05`
- `top_hallmark_up`: top positively enriched Hallmark pathway, chosen by smallest `padj`, then largest `NES`
- `top_hallmark_up_nes`: `NES` of the top positive Hallmark pathway
- `top_hallmark_down`: top negatively enriched Hallmark pathway, chosen by smallest `padj`, then most negative `NES`
- `top_hallmark_down_nes`: `NES` of the top negative Hallmark pathway
- `n_kegg_sig`: number of KEGG pathways with `padj < 0.05`
- `top_kegg_up`: top positively enriched KEGG pathway
- `top_kegg_up_nes`: `NES` of the top positive KEGG pathway
- `top_kegg_down`: top negatively enriched KEGG pathway
- `top_kegg_down_nes`: `NES` of the top negative KEGG pathway

### `rankings/cluster_<ID>_prerank_stats.csv`

The ranked-gene table actually supplied to fgsea for one cluster.

Columns:

- `gene_key`: cleaned uppercase gene symbol used internally for matching
- `gene_symbol`: cleaned display gene symbol
- `rank_stat`: the pre-ranked statistic `sign(logFC) * -log10(p_value)`
- `lfc_value`: numeric log fold change used to determine direction
- `p_value_num`: numeric p value used in ranking
- `p_val_adj`: adjusted p value from the DEG input table, if available
- `pct.1`: fraction of expressing cells in the focal cluster, if available
- `pct.2`: fraction of expressing cells outside the focal cluster, if available

### `run_summary.txt`

Plain-text execution summary describing:

- input DEG directory,
- fallback Seurat object,
- cluster labels analyzed,
- Hallmark and KEGG collection sources,
- ranking metric,
- significance cutoff,
- gene-set size limits,
- whether any DEG files had to be recomputed.

## 01_hallmark_prerank_GSEA

Each `cluster_<ID>/` directory contains Hallmark GSEA outputs for one `cluster_final`.

### `cluster_<ID>_Hallmark_prerank_GSEA.csv`

Per-cluster Hallmark GSEA table.

Column definitions are identical to those described for `00_summary/Hallmark_prerank_GSEA_all.csv`, but restricted to one cluster.

Rows are ordered by:

- ascending adjusted p value,
- descending absolute normalized enrichment score,
- pathway name.

### `cluster_<ID>_Hallmark_prerank_GSEA_top.pdf`, `.png`, `.tiff`

Top-pathway bar plot for Hallmark GSEA in one cluster.

Data source:

- `cluster_<ID>_Hallmark_prerank_GSEA.csv`

How terms are selected:

- if significant pathways exist, the plot uses only pathways with `padj < 0.05`,
- the top positive pathways are selected by smallest `padj`, then largest `NES`,
- the top negative pathways are selected by smallest `padj`, then most negative `NES`,
- up to 10 pathways per direction are shown.

Axes:

- x axis: pathway labels
- y axis: `NES`

Legend:

- fill color encodes `direction`
- red bars represent `up` enrichment (`NES >= 0`)
- blue bars represent `down` enrichment (`NES < 0`)

## 02_kegg_prerank_GSEA

Each `cluster_<ID>/` directory contains KEGG GSEA outputs for one `cluster_final`.

### `cluster_<ID>_KEGG_prerank_GSEA.csv`

Per-cluster KEGG GSEA table.

Column definitions are identical to those described for `00_summary/KEGG_prerank_GSEA_all.csv`, but restricted to one cluster.

### `cluster_<ID>_KEGG_prerank_GSEA_top.pdf`, `.png`, `.tiff`

Top-pathway bar plot for KEGG GSEA in one cluster.

Data source:

- `cluster_<ID>_KEGG_prerank_GSEA.csv`

Axes:

- x axis: pathway labels
- y axis: `NES`

Legend:

- fill color encodes `direction`
- red bars represent `up` enrichment
- blue bars represent `down` enrichment

Selection rules are the same as for the Hallmark top plot.

## 03_plots

### `Hallmark_prerank_GSEA_NES_heatmap.pdf`, `.png`, `.tiff`

Cross-cluster heatmap of significant Hallmark GSEA results.

Data source:

- `00_summary/Hallmark_prerank_GSEA_all.csv`

Matrix definition:

- rows: the top 20 Hallmark pathways selected across clusters
- columns: `cluster_final` labels
- eligible rows: only pathways with `padj < 0.05`
- row selection across clusters: by smallest `padj`, then largest absolute `NES`
- cell value: `NES`
- if the same cluster and pathway appear more than once, the row with the smallest `padj` is used
- clusters with no significant entry for a displayed pathway receive `0`

Legend:

- the heatmap color bar represents `NES`
- positive values correspond to enrichment among genes ranked toward the top of the cluster-specific list
- negative values correspond to enrichment among genes ranked toward the bottom of the list

### `KEGG_prerank_GSEA_NES_heatmap.pdf`, `.png`, `.tiff`

Cross-cluster heatmap of significant KEGG GSEA results.

Data source:

- `00_summary/KEGG_prerank_GSEA_all.csv`

Matrix definition and legend are the same as for the Hallmark heatmap, except that the pathways come from the KEGG collection used in the run.

## Notes on Interpretation

- `ES` and `NES` are directional. A positive `NES` means the pathway tends to accumulate near the top of the ranked gene list, which corresponds to genes higher in the focal cluster than in the rest.
- A negative `NES` means the pathway tends to accumulate near the bottom of the ranked list, which corresponds to relative depletion in the focal cluster.
- `leadingEdge_genes` are not all genes in the pathway. They are the subset driving the enrichment signal according to fgsea.
- The exact KEGG pathway universe depends on which MSigDB KEGG subcollection is available in the installed `msigdbr` version. Always consult `collection_source` when comparing results across runs or environments.
