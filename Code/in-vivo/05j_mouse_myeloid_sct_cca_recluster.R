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

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  a <- sum(choose2(tab))
  row_sum <- sum(choose2(rowSums(tab)))
  col_sum <- sum(choose2(colSums(tab)))
  total <- choose2(sum(tab))
  expected <- if (total == 0) 0 else row_sum * col_sum / total
  maximum <- (row_sum + col_sum) / 2
  if (maximum == expected) return(1)
  (a - expected) / (maximum - expected)
}

script_dir <- resolve_script_dir()
source(file.path(script_dir, "Utils.R"))
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")
required_packages <- c("yaml", "Seurat", "future", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
input_rds <- file.path(results_root, "05i_myeloid_extraction", "mouse_myeloid_rna_only_preintegration.rds")
out_dir <- file.path(results_root, "05j_myeloid_sct_cca_recluster")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds)) stop("Missing 05i input: ", input_rds, call. = FALSE)

icfg <- cfg$integration
seed <- as.integer(icfg$seed)
nfeatures <- as.integer(icfg$nfeatures)
reference_samples <- as.character(unlist(icfg$reference_samples, use.names = FALSE))
npcs_compute <- as.integer(icfg$npcs_compute)
dims_requested <- as.integer(icfg$dims_use)
resolutions <- sort(unique(as.numeric(unlist(cfg$myeloid$cluster_resolutions, use.names = FALSE))))
primary_resolution <- as.numeric(cfg$myeloid$primary_resolution)
if (!primary_resolution %in% resolutions) stop("myeloid.primary_resolution must be included in cluster_resolutions.", call. = FALSE)

myeloid <- readRDS(input_rds)
if (!inherits(myeloid, "Seurat")) stop("05i input is not a Seurat object.", call. = FALSE)
if (!identical(names(myeloid@assays), "RNA")) stop("05i input must contain only the RNA assay.", call. = FALSE)
if (sum(startsWith(rownames(myeloid), as.character(cfg$species$human_prefix))) != 0L) stop("Human features detected before myeloid integration.", call. = FALSE)
seurat_list <- Seurat::SplitObject(myeloid, split.by = "sample")
if (length(seurat_list) != length(cfg$sample_sources)) stop("Expected all configured tissue samples before myeloid CCA.", call. = FALSE)
reference_indices <- match(reference_samples, names(seurat_list))
if (anyNA(reference_indices)) stop("Configured CCA references are absent from the myeloid compartment.", call. = FALSE)

future_plan <- configure_future_for_seurat(
  max_size_gb = as.numeric(icfg$future_globals_maxsize_gb),
  strategy = as.character(icfg$future_strategy),
  workers = as.integer(icfg$future_workers)
)
options(future.rng.onMisuse = "error")
prepared_checkpoint <- file.path(out_dir, "checkpoint_myeloid_sct_prepared.rds")
if (file.exists(prepared_checkpoint)) {
  checkpoint <- readRDS(prepared_checkpoint)
  expected_counts <- vapply(seurat_list, ncol, integer(1))
  if (!identical(checkpoint$sample_names, names(seurat_list)) || !identical(checkpoint$cell_counts, expected_counts) || as.integer(checkpoint$nfeatures) != nfeatures) {
    stop("05j SCT checkpoint does not match the current 05i input/configuration.", call. = FALSE)
  }
  seurat_list <- checkpoint$seurat_list
  integration_features <- checkpoint$integration_features
  rm(checkpoint)
  invisible(gc())
} else {
  set.seed(seed)
  seurat_list <- lapply(seurat_list, function(x) Seurat::SCTransform(
    x, assay = "RNA", new.assay.name = "SCT", vars.to.regress = "percent.mt",
    return.only.var.genes = FALSE, verbose = FALSE
  ))
  integration_features <- Seurat::SelectIntegrationFeatures(seurat_list, nfeatures = nfeatures)
  seurat_list <- Seurat::PrepSCTIntegration(seurat_list, anchor.features = integration_features, verbose = FALSE)
  saveRDS(list(
    sample_names = names(seurat_list), cell_counts = vapply(seurat_list, ncol, integer(1)),
    nfeatures = nfeatures, integration_features = integration_features, seurat_list = seurat_list
  ), prepared_checkpoint, compress = FALSE)
}

minimum_sample_cells <- min(vapply(seurat_list, ncol, integer(1)))
if (minimum_sample_cells < 6L) stop("Every sample requires at least six myeloid cells for CCA.", call. = FALSE)
k_weight <- min(100L, minimum_sample_cells - 1L)
anchor_checkpoint <- file.path(out_dir, "checkpoint_myeloid_reference_cca_anchors.rds")
if (file.exists(anchor_checkpoint)) {
  checkpoint <- readRDS(anchor_checkpoint)
  if (!identical(checkpoint$sample_names, names(seurat_list)) || !identical(checkpoint$reference_samples, reference_samples) || as.integer(checkpoint$nfeatures) != nfeatures) {
    stop("05j anchor checkpoint does not match the current input/configuration.", call. = FALSE)
  }
  anchors <- checkpoint$anchors
  rm(checkpoint)
  invisible(gc())
} else {
  set.seed(seed)
  anchors <- Seurat::FindIntegrationAnchors(
    object.list = seurat_list, normalization.method = "SCT", anchor.features = integration_features,
    reduction = "cca", reference = reference_indices, verbose = FALSE
  )
  saveRDS(list(
    sample_names = names(seurat_list), reference_samples = reference_samples,
    nfeatures = nfeatures, anchors = anchors
  ), anchor_checkpoint, compress = FALSE)
}

integrated_checkpoint <- file.path(out_dir, "checkpoint_myeloid_integrated_precluster.rds")
if (file.exists(integrated_checkpoint)) {
  integrated <- readRDS(integrated_checkpoint)
  if (ncol(integrated) != ncol(myeloid) || !setequal(colnames(integrated), colnames(myeloid))) stop("05j integrated checkpoint cell universe mismatch.", call. = FALSE)
} else {
  # Seurat 4.4 calls future_lapply() inside IntegrateData without exposing a
  # future.seed argument. Keep this stochastic/internal step sequential under
  # the explicit stage seed so future.rng.onMisuse="error" remains an active
  # reproducibility guard instead of being disabled.
  integration_plan <- future::plan()
  future::plan(future::sequential)
  set.seed(seed)
  integrated <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = k_weight, verbose = FALSE)
  future::plan(integration_plan)
  saveRDS(integrated, integrated_checkpoint, compress = FALSE)
}

Seurat::DefaultAssay(integrated) <- "integrated"
parallel_plan <- future::plan()
future::plan(future::sequential)
set.seed(seed)
integrated <- Seurat::RunPCA(integrated, npcs = npcs_compute, seed.use = seed, verbose = FALSE)
dims_use <- seq_len(min(dims_requested, ncol(Seurat::Embeddings(integrated, "pca"))))
integrated <- Seurat::FindNeighbors(integrated, dims = dims_use, verbose = FALSE)
set.seed(seed)
integrated <- Seurat::RunUMAP(integrated, dims = dims_use, seed.use = seed, verbose = FALSE)

resolution_columns <- character(length(resolutions))
for (i in seq_along(resolutions)) {
  resolution <- resolutions[[i]]
  set.seed(seed)
  integrated <- Seurat::FindClusters(integrated, resolution = resolution, random.seed = seed, verbose = FALSE)
  resolution_columns[[i]] <- paste0("integrated_snn_res.", format(resolution, trim = TRUE, scientific = FALSE))
  if (!resolution_columns[[i]] %in% colnames(integrated@meta.data)) stop("FindClusters did not create expected column: ", resolution_columns[[i]], call. = FALSE)
}
primary_column <- paste0("integrated_snn_res.", format(primary_resolution, trim = TRUE, scientific = FALSE))
integrated$seurat_clusters <- factor(as.character(integrated@meta.data[[primary_column]]), levels = sort(unique(as.character(integrated@meta.data[[primary_column]]))))
Seurat::Idents(integrated) <- integrated$seurat_clusters
future::plan(parallel_plan)

stability <- data.frame(
  resolution = resolutions,
  metadata_column = resolution_columns,
  n_clusters = vapply(resolution_columns, function(x) length(unique(integrated@meta.data[[x]])), integer(1)),
  ari_to_primary = vapply(resolution_columns, function(x) adjusted_rand_index(integrated@meta.data[[x]], integrated@meta.data[[primary_column]]), numeric(1)),
  ari_to_previous = NA_real_, stringsAsFactors = FALSE
)
if (nrow(stability) > 1L) {
  for (i in 2:nrow(stability)) stability$ari_to_previous[[i]] <- adjusted_rand_index(integrated@meta.data[[resolution_columns[[i - 1L]]]], integrated@meta.data[[resolution_columns[[i]]]])
}
utils::write.csv(stability, file.path(out_dir, "myeloid_resolution_stability.csv"), row.names = FALSE)

cluster_counts <- as.data.frame(table(cluster = integrated$seurat_clusters), stringsAsFactors = FALSE)
colnames(cluster_counts)[2] <- "n_cells"
utils::write.csv(cluster_counts, file.path(out_dir, "myeloid_primary_cluster_counts.csv"), row.names = FALSE)
cluster_sample <- as.data.frame.matrix(table(cluster = integrated$seurat_clusters, sample = integrated$sample))
cluster_sample <- cbind(cluster = rownames(cluster_sample), cluster_sample)
rownames(cluster_sample) <- NULL
utils::write.csv(cluster_sample, file.path(out_dir, "myeloid_cluster_by_sample_counts.csv"), row.names = FALSE)

p_cluster <- Seurat::DimPlot(integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, raster = FALSE) +
  ggplot2::labs(title = "Mouse myeloid compartment: SCT-CCA reclustering", subtitle = paste0("Primary resolution = ", primary_resolution))
p_sample <- Seurat::DimPlot(integrated, reduction = "umap", group.by = "sample", raster = FALSE) +
  ggplot2::labs(title = "Mouse myeloid compartment by sample")
ggplot2::ggsave(file.path(out_dir, "umap_mouse_myeloid_clusters.pdf"), p_cluster, width = 10, height = 8)
ggplot2::ggsave(file.path(out_dir, "umap_mouse_myeloid_clusters.png"), p_cluster, width = 10, height = 8, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "umap_mouse_myeloid_samples.pdf"), p_sample, width = 12, height = 9)
ggplot2::ggsave(file.path(out_dir, "umap_mouse_myeloid_samples.png"), p_sample, width = 12, height = 9, dpi = 300)

integrated@misc$mouse05_myeloid_integration_contract <- list(
  normalization_method = "SCT", reduction = "cca", integration_mode = "reference_based",
  reference_samples = reference_samples, integration_features = nfeatures, k_weight = k_weight,
  pca_npcs_computed = npcs_compute, dims_used = dims_use, resolutions = resolutions,
  primary_resolution = primary_resolution, primary_cluster_column = primary_column,
  stochastic_future_strategy = "sequential for IntegrateData/PCA/UMAP/clustering", seed = seed, configured_future = future_plan,
  human_features_retained_in_RNA = sum(startsWith(rownames(integrated[["RNA"]]), as.character(cfg$species$human_prefix)))
)
if (integrated@misc$mouse05_myeloid_integration_contract$human_features_retained_in_RNA != 0L) stop("Human features detected after myeloid CCA.", call. = FALSE)
saveRDS(integrated, file.path(out_dir, "mouse_myeloid_sct_cca_reclustered.rds"), compress = FALSE)
utils::write.csv(data.frame(
  parameter = c("samples", "cells", "nfeatures", "reference_samples", "k_weight", "npcs_compute", "dims_use", "resolution_grid", "primary_resolution", "seed"),
  value = c(length(seurat_list), ncol(integrated), nfeatures, paste(reference_samples, collapse = ";"), k_weight, npcs_compute, paste(dims_use, collapse = ","), paste(resolutions, collapse = ","), primary_resolution, seed)
), file.path(out_dir, "myeloid_integration_clustering_parameters.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05j.txt"))
writeLines(c(
  "05j mouse myeloid SCT-CCA integration and reclustering completed.",
  paste0("Cells: ", ncol(integrated)), paste0("Samples: ", length(seurat_list)),
  paste0("Primary clusters: ", length(unique(integrated$seurat_clusters))), paste0("Primary resolution: ", primary_resolution)
), file.path(out_dir, "completion.txt"))
message("[05j] Completed. Output: ", out_dir)
