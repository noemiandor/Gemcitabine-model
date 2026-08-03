# 02e Cluster Annotation: File and Figure Guide

This document explains the contents of the `02e_cluster_annotation` result directory produced by `Code/in-vivo/02e_cluster_annotation.R`.

The goal of this step is to:
- merge selected refined Seurat clusters into `cluster_final`,
- compute `cluster_final` versus rest differential expression,
- select an up-regulated marker set for Hallmark over-representation analysis (ORA),
- assign Hallmark-based cluster annotations,
- write annotation labels back to the Seurat object,
- export summary tables and diagnostic plots.

Unless otherwise stated, all cluster IDs in this directory refer to `cluster_final`, not the original `seurat_cluster_refine` labels.

## Directory Structure

- `00_summary/`
  Central summary tables describing cluster mapping, cell counts, differential expression summaries, Hallmark ORA summaries, and the final annotation calls.
- `01_cluster_final_vs_rest_DEG/`
  Per-cluster marker tables used as input to Hallmark ORA.
- `02_hallmark_ORA/`
  Per-cluster Hallmark ORA result tables and top-term bar plots.
- `03_plots/`
  Cross-cluster summary plots, including composition plots, heatmaps, and UMAP visualizations.
- `04_objects/`
  The annotated Seurat object with new cluster-level annotation metadata columns.

## Cluster Definition Used in This Step

The script creates `cluster_final` by applying the following merges to `seurat_cluster_refine`:

- `1 -> 0`
- `7 -> 0`
- `11 -> 10`
- `12 -> 10`

All downstream DEG, ORA, summaries, and plots in this directory are based on `cluster_final`.

## Selection Rules and Key Calculations

### Differential Expression

Per cluster, differential expression is computed with `Seurat::FindMarkers` using:

- comparison: `cluster_final X` versus all remaining cells,
- assay: `RNA`,
- slot: `data`,
- test: Wilcoxon,
- `min.pct = 0.10`,
- `logfc.threshold = 0`.

The DEG table is sorted by:

- ascending adjusted p value (`p_val_adj`),
- then descending absolute log fold change.

### Gene Set Used for Hallmark ORA

The Hallmark ORA input gene list is derived from the DEG table as follows:

- gene symbols are cleaned to remove prefixes such as `GRCh38-`,
- duplicate genes are collapsed by keeping the row with the largest absolute log fold change, then the smallest adjusted p value,
- genes must satisfy all of the following:
  - `p_val_adj < 0.05`
  - `abs(logFC) >= 0.25`
  - `abs(pct.1 - pct.2) >= 0.05`
  - positive logFC only
- the top 100 up-regulated genes are kept after ranking by descending logFC and ascending adjusted p value.

### Hallmark ORA

Hallmark gene sets come from MSigDB Hallmark (`H`) through `msigdbr`.

The enrichment test is a hypergeometric over-representation test using:

- query genes: the top 100 up-regulated genes described above,
- universe genes: all non-missing cleaned gene symbols present in the RNA assay,
- minimum pathway size: 15,
- maximum pathway size: 500,
- minimum overlap: 3 genes,
- multiple-testing correction: Benjamini-Hochberg FDR (`p_adj`).

### Annotation Score

After Hallmark ORA, the script computes an integrated annotation score for significant Hallmark terms only (`p_adj < 0.05`).

For each Hallmark term:

- `overlap_query_fraction = overlap / query_size`
- `overlap_set_fraction = overlap / set_size`
- `overlap_score_raw = mean(overlap_query_fraction, overlap_set_fraction)`
- `marker_mean_abs_logfc = mean(abs(logFC))` across overlapping genes
- `marker_mean_pct1 = mean(pct.1)` across overlapping genes
- `marker_mean_abs_delta_pct = mean(abs(pct.1 - pct.2))` across overlapping genes

Three raw components are min-max scaled within the cluster:

- `expression_score_scaled` from `marker_mean_abs_logfc`
- `detection_score_scaled` from `marker_mean_pct1`
- `overlap_score_scaled` from `overlap_score_raw`

The final integrated score is:

`annotation_score = (expression + detection + overlap) / 3`

because all three weights are set to one third in the current script.

The primary, secondary, and tertiary annotation labels are the top Hallmark terms after sorting by:

- descending `annotation_score`,
- descending `expression_score_raw`,
- descending `detection_score_raw`,
- descending `overlap_score_raw`,
- then pathway name.

## 00_summary

### `cluster_merge_mapping_for_annotation.csv`

Maps original refined cluster labels to `cluster_final`.

Columns:

- `original_cluster`: label in `seurat_cluster_refine`
- `analysis_cluster`: resulting `cluster_final` label after hard-coded merges
- `retained_cluster_label`: original label if retained unchanged, blank if absorbed into another cluster
- `merge_status`: `retained` or `absorbed`

### `cluster_final_cell_counts.csv`

Number of cells assigned to each `cluster_final`.

Columns:

- `cluster`: `cluster_final` label
- `n_cells`: number of cells in that cluster

### `hallmark_ora_universe_genes.csv`

The Hallmark ORA background universe.

Columns:

- `gene_symbol`: cleaned gene symbol present in the RNA assay and eligible for ORA background counting

### `cluster_final_deg_summary.csv`

One row per `cluster_final`, summarizing DEG availability and the ORA input set.

Columns:

- `cluster`: `cluster_final` label
- `n_cells`: number of cells in the cluster
- `tested_genes`: number of genes returned in the `FindMarkers` table
- `sig_genes`: number of DEG rows with `p_val_adj < 0.05`
- `deg_for_ora`: number of genes retained in the top-100 up-regulated Hallmark ORA input set
- `up_genes_for_ora`: count of ORA-input genes labeled `direction = "up"`
- `down_genes_for_ora`: always `0` in the current version because ORA uses only up-regulated genes

### `cluster_final_Hallmark_ORA_all.csv`

Combined Hallmark ORA table across all clusters.

Columns:

- `cluster`: `cluster_final` label
- `annotation_rank`: rank after integrated-score sorting within the cluster
- `pathway`: original Hallmark pathway ID from MSigDB
- `hallmark_label`: human-readable Hallmark label
- `set_size`: number of background genes in the Hallmark set after intersecting with the universe
- `query_size`: number of ORA-input genes submitted for that cluster
- `overlap`: number of ORA-input genes found in the Hallmark set
- `overlap_query_fraction`: `overlap / query_size`
- `overlap_set_fraction`: `overlap / set_size`
- `overlap_score_raw`: mean of the two overlap fractions above
- `marker_mean_abs_logfc`: mean absolute log fold change of overlapping genes
- `marker_mean_pct1`: mean `pct.1` of overlapping genes
- `marker_mean_abs_delta_pct`: mean absolute `pct.1 - pct.2` of overlapping genes
- `n_overlap_up`: number of overlapping genes with positive direction in the ORA input table
- `n_overlap_down`: number of overlapping genes with negative direction in the ORA input table; expected to be `0` in the current workflow
- `expression_score_raw`: identical to `marker_mean_abs_logfc`
- `detection_score_raw`: identical to `marker_mean_pct1`
- `expression_score_scaled`: min-max scaled expression component within the cluster
- `detection_score_scaled`: min-max scaled detection component within the cluster
- `overlap_score_scaled`: min-max scaled overlap component within the cluster
- `annotation_score`: mean of the three scaled components
- `gene_ratio`: formatted overlap fraction as `overlap/query_size`
- `bg_ratio`: formatted background fraction as `set_size/universe_size`
- `odds_ratio`: enrichment odds ratio from the hypergeometric table
- `p_value`: raw hypergeometric p value
- `p_adj`: Benjamini-Hochberg adjusted p value
- `overlap_genes`: semicolon-separated overlapping gene symbols

### `cluster_final_annotation_summary.csv`

The final per-cluster annotation call table.

Columns:

- `cluster`: `cluster_final` label
- `n_cells`: number of cells in the cluster
- `n_deg_for_ora`: number of genes in the ORA input table for that cluster
- `n_up_deg_for_ora`: count of ORA-input up-regulated genes
- `n_down_deg_for_ora`: count of ORA-input down-regulated genes
- `n_significant_hallmarks`: number of Hallmark ORA rows with `p_adj < 0.05`
- `annotation_primary`: top Hallmark label by integrated annotation score
- `annotation_secondary`: second-ranked Hallmark label
- `annotation_tertiary`: third-ranked Hallmark label
- `annotation_multi`: semicolon-separated concatenation of the top ranked labels
- `top_hallmark_1_annotation_score`: integrated score of the primary annotation
- `top_hallmark_2_annotation_score`: integrated score of the secondary annotation
- `top_hallmark_3_annotation_score`: integrated score of the tertiary annotation
- `top_hallmark_1_fdr`: adjusted p value of the primary annotation
- `top_hallmark_2_fdr`: adjusted p value of the secondary annotation
- `top_hallmark_3_fdr`: adjusted p value of the tertiary annotation
- `note`: interpretation note, including whether significant Hallmark evidence was found

### `run_summary.txt`

Plain-text execution summary describing:

- input object,
- output directory,
- cluster labels,
- merge rules,
- ORA gene-selection rule,
- annotation-score rule,
- major output locations.

## 01_cluster_final_vs_rest_DEG

Each `cluster_<ID>/` directory contains DEG files for one `cluster_final`.

### `cluster_<ID>_vs_rest_markers.csv`

The full differential expression table for `cluster_final <ID>` versus all other cells.

Common columns:

- `p_val`: raw p value from `FindMarkers`
- `avg_log2FC` or another recognized logFC column: effect size used by the script
- `pct.1`: fraction of cells expressing the gene in the focal cluster
- `pct.2`: fraction of cells expressing the gene in all remaining cells
- `p_val_adj`: adjusted p value
- `gene`: gene identifier, often with a genome-prefix such as `GRCh38-`
- `gene_symbol`: cleaned gene symbol when available

### `cluster_<ID>_vs_rest_top100_up_for_Hallmark_ORA.csv`

The actual Hallmark ORA input gene table used by the current version of the script.

This file contains the top 100 up-regulated genes that pass the filtering rules described above.

Additional useful columns may include:

- `gene_key`: uppercase normalized gene symbol used for uniqueness matching
- `lfc_value`: numeric log fold change used in filtering
- `delta_pct`: `pct.1 - pct.2`
- `abs_logfc`: absolute log fold change
- `abs_delta_pct`: absolute detection difference
- `direction`: always `up` in the current workflow

### `cluster_<ID>_vs_rest_top50_each_for_Hallmark_ORA.csv`

This file is a legacy artifact from an older version of the workflow that used top up-regulated and down-regulated genes separately for ORA.

It is present in the directory but is not used by the current `02e_cluster_annotation.R` script.

## 02_hallmark_ORA

Each `cluster_<ID>/` directory contains one Hallmark ORA result table and one bar plot.

### `cluster_<ID>_Hallmark_ORA.csv`

Per-cluster Hallmark ORA table with the same column definitions as `00_summary/cluster_final_Hallmark_ORA_all.csv`, but restricted to a single cluster.

Rows are ordered by integrated annotation rank.

### `cluster_<ID>_Hallmark_ORA_top.pdf`, `.png`, `.tiff`

Bar plot of the top Hallmark terms for one cluster.

Data source:

- `cluster_<ID>_Hallmark_ORA.csv`

How bars are chosen:

- if `annotation_score` is available, the plot uses the integrated-score ranking,
- otherwise it falls back to significant ORA terms ordered by FDR.

Axes:

- y axis: either `Integrated annotation score` or `-log10(FDR)`, depending on available columns
- x axis: Hallmark term labels

Legend:

- no categorical legend is used in this plot
- bar height itself is the encoded quantity

## 03_plots

### `sample_type_cluster_final_proportion.csv`

Composition table used for the stacked bar plot of cluster distribution within each sample type.

Columns:

- `sample_type`: group label used on the x axis
- `cluster_final`: stacked category
- `count`: number of cells in that `sample_type` and `cluster_final`
- `proportion`: `count / sum(count within the same sample_type)`

### `cluster_final_sample_type_proportion.csv`

Composition table used for the stacked bar plot of sample-type distribution within each cluster.

Columns:

- `cluster_final`: group label used on the x axis
- `sample_type`: stacked category
- `count`: number of cells in that `cluster_final` and `sample_type`
- `proportion`: `count / sum(count within the same cluster_final)`

### `stack_sample_type_by_cluster_final.pdf`, `.png`, `.tiff`

Stacked bar plot showing the proportion of each `cluster_final` within each `sample_type`.

Data source:

- `sample_type_cluster_final_proportion.csv`

Axes:

- x axis: `sample_type`
- y axis: within-sample-type proportion

Legend:

- fill color encodes `cluster_final`

### `stack_cluster_final_by_sample_type.pdf`, `.png`, `.tiff`

Stacked bar plot showing the proportion of each `sample_type` within each `cluster_final`.

Data source:

- `cluster_final_sample_type_proportion.csv`

Axes:

- x axis: `cluster_final`
- y axis: within-cluster proportion

Legend:

- fill color encodes `sample_type`

### `Hallmark_annotation_heatmap.pdf`, `.tiff`

Cross-cluster heatmap of Hallmark annotation signals.

Data source:

- `00_summary/cluster_final_Hallmark_ORA_all.csv`

Matrix definition:

- rows: the top 20 Hallmark terms selected across clusters
- columns: `cluster_final` labels
- cell value:
  - `annotation_score` if annotation scores are available,
  - otherwise `-log10(FDR)`
- when multiple rows exist for the same cluster and Hallmark label, the maximum score is used

Legend:

- the heatmap color bar represents the numeric cell value
- larger values indicate stronger support for that Hallmark in that cluster

### `umap_cluster_final.pdf`, `.png`, `.tiff`

UMAP embedding colored by `cluster_final`.

Data source:

- the Seurat object after cluster merging

Legend:

- point color encodes `cluster_final`
- text labels identify clusters on the embedding

### `umap_cluster_final_annotation_primary.pdf`, `.png`, `.tiff`

UMAP embedding colored by the primary Hallmark annotation assigned to each cluster.

Data source:

- the Seurat object after writing `cluster_final_annotation_primary` into metadata

Legend:

- point color encodes `cluster_final_annotation_primary`
- text labels identify the primary annotation categories on the embedding

## 04_objects

### `integrated_sct_cca_seurat_cluster_final_annotation.rds`

Annotated Seurat object written at the end of the workflow.

New metadata columns added by this step:

- `cluster_final`
- `cluster_final_annotation_primary`
- `cluster_final_annotation_multi`

This object is the direct bridge between the cluster annotation results and any downstream single-cell analysis that requires the final cluster labels or Hallmark-based labels.
