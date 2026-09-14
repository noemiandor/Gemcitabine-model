#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
)

resolve_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(hit) != 1L) stop("Cannot resolve script directory.", call. = FALSE)
  dirname(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = TRUE))
}

script_dir <- resolve_script_dir()
source(file.path(script_dir, "Utils.R"))
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")
required_packages <- c("yaml", "Seurat", "Matrix", "future")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
input_rds <- file.path(results_root, "05j_myeloid_sct_cca_recluster", "mouse_myeloid_sct_cca_reclustered.rds")
out_dir <- file.path(results_root, "05k_myeloid_cluster_degs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds)) stop("Missing 05j input: ", input_rds, call. = FALSE)

icfg <- cfg$integration
configure_future_for_seurat(
  max_size_gb = as.numeric(icfg$future_globals_maxsize_gb),
  strategy = as.character(icfg$future_strategy),
  workers = as.integer(icfg$future_workers)
)
options(future.rng.onMisuse = "error")

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("05j input is not a Seurat object.", call. = FALSE)
if (!all(c("sample", "seurat_clusters") %in% colnames(obj@meta.data))) stop("05j input lacks sample or seurat_clusters.", call. = FALSE)
if (sum(startsWith(rownames(obj[["RNA"]]), as.character(cfg$species$human_prefix))) != 0L) stop("Human features detected before myeloid DEG analysis.", call. = FALSE)

Seurat::DefaultAssay(obj) <- "RNA"
obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Seurat::Idents(obj) <- obj$seurat_clusters
markers_file <- file.path(out_dir, "myeloid_cluster_positive_degs_full.csv")
if (file.exists(markers_file) && file.info(markers_file)$size > 0L) {
  message("[05k] Reusing completed DEG table: ", markers_file)
  markers <- utils::read.csv(markers_file, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  markers <- Seurat::FindAllMarkers(
    obj, assay = "RNA", slot = "data", only.pos = TRUE,
    min.pct = as.numeric(cfg$markers$min_pct), logfc.threshold = as.numeric(cfg$markers$logfc_threshold),
    verbose = FALSE
  )
}
if (nrow(markers) == 0L) stop("No positive myeloid cluster DEGs were identified.", call. = FALSE)
fc_column <- intersect(c("avg_log2FC", "avg_logFC"), colnames(markers))
if (length(fc_column) != 1L) stop("Myeloid DEG output lacks exactly one fold-change column.", call. = FALSE)
markers$cluster <- as.character(markers$cluster)
markers$gene_upper <- toupper(trimws(as.character(markers$gene)))
markers$rank_fc <- as.numeric(markers[[fc_column]])
markers <- markers[order(markers$cluster, -markers$rank_fc, markers$gene_upper), , drop = FALSE]
utils::write.csv(markers, markers_file, row.names = FALSE, na = "NA")

ranked_markers <- do.call(rbind, lapply(split(markers, markers$cluster), function(z) {
  z <- z[!duplicated(z$gene_upper), , drop = FALSE]
  z$marker_rank <- seq_len(nrow(z))
  z
}))
rownames(ranked_markers) <- NULL
for (top_n in c(20L, 50L, 100L)) {
  local <- ranked_markers[ranked_markers$marker_rank <= top_n, , drop = FALSE]
  utils::write.csv(local, file.path(out_dir, paste0("myeloid_cluster_top", top_n, "_positive_degs.csv")), row.names = FALSE, na = "NA")
}

cluster_ids <- sort(unique(as.character(obj$seurat_clusters)))
pair_grid <- if (length(cluster_ids) >= 2L) t(utils::combn(cluster_ids, 2L)) else matrix(character(), ncol = 2L)
jaccard_rows <- list()
for (top_n in c(20L, 50L, 100L, 200L, 500L)) {
  for (i in seq_len(nrow(pair_grid))) {
    a <- pair_grid[i, 1L]
    b <- pair_grid[i, 2L]
    genes_a <- unique(ranked_markers$gene_upper[ranked_markers$cluster == a & ranked_markers$marker_rank <= top_n])
    genes_b <- unique(ranked_markers$gene_upper[ranked_markers$cluster == b & ranked_markers$marker_rank <= top_n])
    overlap <- intersect(genes_a, genes_b)
    union_genes <- union(genes_a, genes_b)
    jaccard_rows[[length(jaccard_rows) + 1L]] <- data.frame(
      top_n = top_n, cluster_a = a, cluster_b = b, n_a = length(genes_a), n_b = length(genes_b),
      overlap_n = length(overlap), jaccard = if (length(union_genes) == 0L) NA_real_ else length(overlap) / length(union_genes),
      overlap_genes = paste(overlap, collapse = ";"), stringsAsFactors = FALSE
    )
  }
}
pairwise_jaccard <- if (length(jaccard_rows) > 0L) do.call(rbind, jaccard_rows) else data.frame()
utils::write.csv(pairwise_jaccard, file.path(out_dir, "myeloid_cluster_deg_pairwise_jaccard.csv"), row.names = FALSE, na = "NA")

# Sample-aware pseudobulk audit. Counts are summed within each sample x cluster,
# converted to log1p(CPM), and correlations are calculated on the union of the
# top 500 cluster DEGs. This is a validation product, not the annotation score.
counts <- Seurat::GetAssayData(obj, assay = "RNA", slot = "counts")
group_id <- paste(as.character(obj$sample), as.character(obj$seurat_clusters), sep = "|||")
group_factor <- factor(group_id, levels = sort(unique(group_id)))
indicator <- Matrix::sparseMatrix(
  i = seq_along(group_factor), j = as.integer(group_factor), x = 1,
  dims = c(length(group_factor), nlevels(group_factor)),
  dimnames = list(colnames(obj), levels(group_factor))
)
pb_counts <- counts %*% indicator
pb_library <- as.numeric(Matrix::colSums(pb_counts))
if (any(pb_library <= 0)) stop("At least one sample x cluster pseudobulk has zero counts.", call. = FALSE)
pb_cpm <- pb_counts %*% Matrix::Diagonal(x = 1e6 / pb_library)
dimnames(pb_cpm) <- dimnames(pb_counts)
pb_logcpm <- log1p(pb_cpm)
group_parts <- strsplit(colnames(pb_logcpm), "|||", fixed = TRUE)
group_meta <- data.frame(
  pseudobulk_id = colnames(pb_logcpm),
  sample = vapply(group_parts, `[[`, character(1), 1L),
  cluster = vapply(group_parts, `[[`, character(1), 2L),
  library_size = pb_library, stringsAsFactors = FALSE
)
marker_union <- unique(ranked_markers$gene[ranked_markers$marker_rank <= 500L])
marker_union <- intersect(marker_union, rownames(pb_logcpm))
cluster_pb_mean <- vapply(cluster_ids, function(cluster_id) {
  idx <- which(group_meta$cluster == cluster_id)
  as.numeric(Matrix::rowMeans(pb_logcpm[marker_union, idx, drop = FALSE]))
}, numeric(length(marker_union)))
rownames(cluster_pb_mean) <- marker_union
colnames(cluster_pb_mean) <- cluster_ids
pb_cor <- stats::cor(cluster_pb_mean, method = "spearman", use = "pairwise.complete.obs")
utils::write.csv(cbind(gene = rownames(cluster_pb_mean), as.data.frame(cluster_pb_mean, check.names = FALSE)), file.path(out_dir, "myeloid_cluster_pseudobulk_logcpm_marker_union.csv"), row.names = FALSE)
utils::write.csv(cbind(cluster = rownames(pb_cor), as.data.frame(pb_cor, check.names = FALSE)), file.path(out_dir, "myeloid_cluster_pseudobulk_spearman_correlation.csv"), row.names = FALSE)
utils::write.csv(group_meta, file.path(out_dir, "myeloid_pseudobulk_sample_cluster_metadata.csv"), row.names = FALSE)
saveRDS(list(logcpm = pb_logcpm, metadata = group_meta, marker_union = marker_union), file.path(out_dir, "myeloid_sample_cluster_pseudobulk.rds"), compress = FALSE)

# Sample, Dose, and Ploidy dominance audits.
sample_metadata <- do.call(rbind, lapply(names(cfg$immune_composition$sample_metadata), function(sample_id) {
  z <- cfg$immune_composition$sample_metadata[[sample_id]]
  data.frame(sample = sample_id, Ploidy = as.character(z$Ploidy), Dose = as.character(z$Dose), stringsAsFactors = FALSE)
}))
cell_meta <- data.frame(cell = colnames(obj), sample = as.character(obj$sample), cluster = as.character(obj$seurat_clusters), stringsAsFactors = FALSE)
idx <- match(cell_meta$sample, sample_metadata$sample)
if (anyNA(idx)) stop("A myeloid sample is absent from configured Ploidy/Dose metadata.", call. = FALSE)
cell_meta$Ploidy <- sample_metadata$Ploidy[idx]
cell_meta$Dose <- sample_metadata$Dose[idx]

cluster_sample_counts <- as.data.frame(table(cluster = cell_meta$cluster, sample = cell_meta$sample), stringsAsFactors = FALSE)
colnames(cluster_sample_counts)[3] <- "n_cells"
cluster_totals <- tapply(cluster_sample_counts$n_cells, cluster_sample_counts$cluster, sum)
cluster_sample_counts$fraction_within_cluster <- cluster_sample_counts$n_cells / unname(cluster_totals[cluster_sample_counts$cluster])
utils::write.csv(cluster_sample_counts, file.path(out_dir, "myeloid_cluster_sample_dominance_counts.csv"), row.names = FALSE)

dominance <- do.call(rbind, lapply(split(cluster_sample_counts, cluster_sample_counts$cluster), function(z) {
  p <- z$fraction_within_cluster[z$fraction_within_cluster > 0]
  data.frame(
    cluster = z$cluster[[1]], n_cells = sum(z$n_cells), samples_detected = sum(z$n_cells > 0),
    max_single_sample_fraction = max(z$fraction_within_cluster),
    normalized_shannon_entropy = if (length(p) <= 1L) 0 else -sum(p * log(p)) / log(nrow(z)),
    dominant_sample = z$sample[[which.max(z$fraction_within_cluster)]], stringsAsFactors = FALSE
  )
}))
utils::write.csv(dominance, file.path(out_dir, "myeloid_cluster_sample_dominance_summary.csv"), row.names = FALSE)

cluster_strata <- as.data.frame(table(cluster = cell_meta$cluster, Ploidy = cell_meta$Ploidy, Dose = cell_meta$Dose), stringsAsFactors = FALSE)
colnames(cluster_strata)[4] <- "n_cells"
utils::write.csv(cluster_strata, file.path(out_dir, "myeloid_cluster_by_ploidy_dose_counts.csv"), row.names = FALSE)

obj@misc$mouse05_myeloid_deg_contract <- list(
  annotation_level = "reclustered myeloid cluster", only_positive = TRUE,
  min_pct = as.numeric(cfg$markers$min_pct), logfc_threshold = as.numeric(cfg$markers$logfc_threshold),
  top_deg_exports = c(20L, 50L, 100L), overlap_audit_top_n = c(20L, 50L, 100L, 200L, 500L),
  pseudobulk_unit = "sample x myeloid cluster", pseudobulk_correlation = "Spearman on cluster-mean log1p(CPM)",
  sample_dominance_audited = TRUE
)
saveRDS(obj, file.path(out_dir, "mouse_myeloid_reclustered_with_degs.rds"), compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05k.txt"))
writeLines(c(
  "05k myeloid cluster DEG analysis completed.", paste0("Cells: ", ncol(obj)),
  paste0("Clusters: ", length(cluster_ids)), paste0("Positive DEG rows: ", nrow(markers)),
  "Cluster similarity, pseudobulk correlation, and sample/Ploidy/Dose dominance audits were exported."
), file.path(out_dir, "completion.txt"))
message("[05k] Completed. Output: ", out_dir)
