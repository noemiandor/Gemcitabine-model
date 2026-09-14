#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220)
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
required_packages <- c("yaml", "Seurat", "future", "ggplot2", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
input_rds <- file.path(results_root, "05c_cca_integration", "mouse_integrated_sct_cca.rds")
out_dir <- file.path(results_root, "05d_clustering")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds)) stop("Missing 05c input: ", input_rds, call. = FALSE)

icfg <- cfg$integration
seed <- as.integer(icfg$seed)
npcs_compute <- as.integer(icfg$npcs_compute)
dims_requested <- as.integer(icfg$dims_use)
resolution <- as.numeric(icfg$cluster_resolution)
configure_future_for_seurat(
  max_size_gb = as.numeric(icfg$future_globals_maxsize_gb),
  strategy = as.character(icfg$future_strategy),
  workers = as.integer(icfg$future_workers)
)
parallel_plan <- future::plan()
options(future.rng.onMisuse = "error")

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("05c input is not a Seurat object.", call. = FALSE)
if (sum(startsWith(rownames(obj[["RNA"]]), cfg$species$human_prefix)) != 0L) stop("Human features detected before clustering.", call. = FALSE)

Seurat::DefaultAssay(obj) <- "integrated"
# Seurat 4.4 parallelizes FindClusters() with future_lapply() without declaring
# future.seed. Run all stochastic embedding/graph/clustering steps sequentially
# with explicit seeds, then restore the configured parallel plan for the
# deterministic normalization and marker calculations below.
future::plan(future::sequential)
set.seed(seed)
obj <- Seurat::RunPCA(obj, npcs = npcs_compute, seed.use = seed, verbose = FALSE)
dims_use <- seq_len(min(dims_requested, ncol(Seurat::Embeddings(obj, "pca"))))
obj <- Seurat::FindNeighbors(obj, dims = dims_use, verbose = FALSE)
set.seed(seed)
obj <- Seurat::FindClusters(obj, resolution = resolution, random.seed = seed, verbose = FALSE)
set.seed(seed)
obj <- Seurat::RunUMAP(obj, dims = dims_use, seed.use = seed, verbose = FALSE)
future::plan(parallel_plan)

# Normalize mouse-only RNA explicitly for markers and reference annotation.
Seurat::DefaultAssay(obj) <- "RNA"
obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Seurat::Idents(obj) <- obj$seurat_clusters
markers <- Seurat::FindAllMarkers(
  obj,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  min.pct = as.numeric(cfg$markers$min_pct),
  logfc.threshold = as.numeric(cfg$markers$logfc_threshold),
  verbose = FALSE
)
utils::write.csv(markers, file.path(out_dir, "cluster_positive_markers.csv"), row.names = FALSE, na = "NA")

cluster_sample <- as.data.frame.matrix(table(cluster = obj$seurat_clusters, sample = obj$sample))
cluster_sample <- cbind(cluster = rownames(cluster_sample), cluster_sample)
rownames(cluster_sample) <- NULL
utils::write.csv(cluster_sample, file.path(out_dir, "cluster_by_sample_counts.csv"), row.names = FALSE)

cluster_counts <- data.frame(
  cluster = names(table(obj$seurat_clusters)),
  n_cells = as.integer(table(obj$seurat_clusters)),
  stringsAsFactors = FALSE
)
utils::write.csv(cluster_counts, file.path(out_dir, "cluster_counts.csv"), row.names = FALSE)

p_cluster <- Seurat::DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, raster = FALSE) +
  ggplot2::labs(title = "Mouse cells: SCT-CCA clusters")
p_sample <- Seurat::DimPlot(obj, reduction = "umap", group.by = "sample", raster = FALSE) +
  ggplot2::labs(title = "Mouse cells by sample")
ggplot2::ggsave(file.path(out_dir, "umap_mouse_clusters.pdf"), p_cluster, width = 9, height = 7)
ggplot2::ggsave(file.path(out_dir, "umap_mouse_clusters.png"), p_cluster, width = 9, height = 7, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "umap_mouse_samples.pdf"), p_sample, width = 11, height = 8)
ggplot2::ggsave(file.path(out_dir, "umap_mouse_samples.png"), p_sample, width = 11, height = 8, dpi = 300)

obj@misc$mouse05_clustering_contract <- list(
  pca_npcs_computed = npcs_compute,
  dims_used = dims_use,
  resolution = resolution,
  seed = seed,
  stochastic_future_strategy = "sequential",
  explicit_seed_arguments = c("RunPCA.seed.use", "FindClusters.random.seed", "RunUMAP.seed.use"),
  future_rng_misuse = "error",
  human_features_retained_in_RNA = 0L
)
saveRDS(obj, file.path(out_dir, "mouse_integrated_sct_cca_clustered.rds"), compress = FALSE)
utils::write.csv(
  data.frame(
    parameter = c("pca_npcs_computed", "dims_used", "cluster_resolution", "seed", "stochastic_future_strategy", "future_rng_misuse", "cells", "clusters"),
    value = c(npcs_compute, paste(dims_use, collapse = ","), resolution, seed, "sequential", "error", ncol(obj), length(unique(obj$seurat_clusters)))
  ),
  file.path(out_dir, "clustering_parameters.csv"), row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05d.txt"))
writeLines(
  c("05d clustering completed.", paste0("Cells: ", ncol(obj)), paste0("Clusters: ", length(unique(obj$seurat_clusters)))),
  file.path(out_dir, "completion.txt")
)
message("[05d] Completed. Output: ", out_dir)
