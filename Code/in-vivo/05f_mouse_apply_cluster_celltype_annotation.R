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

cluster_order <- function(x) {
  x <- unique(as.character(x))
  numeric_x <- suppressWarnings(as.numeric(x))
  if (all(is.finite(numeric_x))) x[order(numeric_x)] else sort(x)
}

script_dir <- resolve_script_dir()
script_path <- file.path(script_dir, "05f_mouse_apply_cluster_celltype_annotation.R")
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")
required_packages <- c("yaml", "Seurat", "ggplot2", "ggrepel", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
input_rds <- file.path(results_root, "05d_clustering", "mouse_integrated_sct_cca_clustered.rds")
review_completion <- file.path(results_root, "05e_broad_annotation_marker_review", "completion.txt")
review_manifest <- file.path(results_root, "05e_broad_annotation_marker_review", "output_manifest_sha256.csv")
annotation_file <- file.path(script_dir, "annotation_resources", "mouse_broad_cluster_annotations.csv")
out_dir <- file.path(results_root, "05f_broad_celltype_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(input_rds, review_completion, review_manifest, annotation_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) stop("Missing required input(s): ", paste(missing_inputs, collapse = ", "), call. = FALSE)

completion_lines <- readLines(review_completion, warn = FALSE)
required_review_contract <- c(
  "Annotation level: cluster only.",
  "No automated cell-type label was calculated.",
  "Manual review is required before any cluster label is applied."
)
if (!all(required_review_contract %in% completion_lines)) {
  stop("05e completion file does not satisfy the manual-review gate.", call. = FALSE)
}

annotation <- utils::read.csv(annotation_file, stringsAsFactors = FALSE, check.names = FALSE, na.strings = character())
required_columns <- c("cluster", "broad_cell_type", "confidence", "evidence_tier", "annotation_basis", "review_status")
missing_columns <- setdiff(required_columns, colnames(annotation))
if (length(missing_columns) > 0L) stop("Annotation mapping lacks columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
annotation$cluster <- as.character(annotation$cluster)
for (column in required_columns) {
  if (any(!nzchar(trimws(as.character(annotation[[column]]))))) stop("Annotation mapping has blank ", column, ".", call. = FALSE)
}
if (anyDuplicated(annotation$cluster)) stop("Annotation mapping contains duplicate clusters.", call. = FALSE)
if (any(annotation$review_status != "manually_approved")) stop("Every mapping row must be manually_approved.", call. = FALSE)
if (!all(annotation$confidence %in% c("high", "medium", "low"))) stop("Invalid confidence value.", call. = FALSE)

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("05d input is not a Seurat object.", call. = FALSE)
if (!"seurat_clusters" %in% colnames(obj@meta.data)) stop("05d object lacks seurat_clusters.", call. = FALSE)
if (!"RNA" %in% names(obj@assays)) stop("05d object lacks the RNA assay.", call. = FALSE)
human_feature_n <- sum(startsWith(rownames(obj[["RNA"]]), as.character(cfg$species$human_prefix)))
if (human_feature_n != 0L) stop("Human features detected before broad cell-type annotation: ", human_feature_n, call. = FALSE)
if (!"umap" %in% names(obj@reductions)) stop("05d object lacks UMAP.", call. = FALSE)

object_clusters <- cluster_order(obj$seurat_clusters)
mapping_clusters <- cluster_order(annotation$cluster)
if (!identical(object_clusters, mapping_clusters)) {
  stop(
    "Annotation mapping and object clusters differ. Missing mapping: ", paste(setdiff(object_clusters, mapping_clusters), collapse = ","),
    "; extra mapping: ", paste(setdiff(mapping_clusters, object_clusters), collapse = ","), call. = FALSE
  )
}

annotation <- annotation[match(object_clusters, annotation$cluster), , drop = FALSE]
cell_cluster <- as.character(obj$seurat_clusters)
lookup <- function(column) unname(stats::setNames(annotation[[column]], annotation$cluster)[cell_cluster])

cell_type_order <- as.character(cfg$annotation$manual_broad_celltype_order)
observed_cell_types <- unique(annotation$broad_cell_type)
if (length(cell_type_order) == 0L || !setequal(cell_type_order, observed_cell_types)) {
  stop("annotation.manual_broad_celltype_order must match the manually approved cell types.", call. = FALSE)
}

obj$mouse_broad_cell_type <- factor(lookup("broad_cell_type"), levels = cell_type_order)
obj$mouse_cell_type <- obj$mouse_broad_cell_type
obj$mouse_broad_confidence <- lookup("confidence")
obj$mouse_broad_evidence_tier <- lookup("evidence_tier")
obj$mouse_broad_annotation_basis <- lookup("annotation_basis")
obj$mouse_broad_review_status <- lookup("review_status")
obj$mouse_broad_annotation_level <- "cluster"

label_n_per_cluster <- tapply(as.character(obj$mouse_broad_cell_type), cell_cluster, function(x) length(unique(x)))
if (any(label_n_per_cluster != 1L)) stop("Broad cell-type annotation is not uniform within every cluster.", call. = FALSE)
if (any(is.na(obj$mouse_broad_cell_type))) stop("At least one cell lacks a broad cell-type label.", call. = FALSE)
if (any(obj$mouse_broad_annotation_level != "cluster")) stop("Annotation-level metadata is inconsistent.", call. = FALSE)

cluster_counts <- as.integer(table(factor(cell_cluster, levels = object_clusters)))
cluster_annotation <- annotation
cluster_annotation$n_cells <- cluster_counts
cluster_annotation$percent_mouse_object <- 100 * cluster_annotation$n_cells / ncol(obj)

celltype_counts <- as.data.frame(table(cell_type = obj$mouse_broad_cell_type), stringsAsFactors = FALSE)
colnames(celltype_counts)[2] <- "n_cells"
celltype_counts$percent_mouse_object <- 100 * celltype_counts$n_cells / ncol(obj)

sample_column <- if ("sample" %in% colnames(obj@meta.data)) "sample" else "orig.ident"
sample_celltype <- as.data.frame.matrix(table(
  sample = as.character(obj@meta.data[[sample_column]]),
  cell_type = as.character(obj$mouse_broad_cell_type)
))
sample_celltype <- cbind(sample = rownames(sample_celltype), sample_celltype)
rownames(sample_celltype) <- NULL

colors <- unlist(cfg$annotation$manual_broad_celltype_colors, use.names = TRUE)
if (is.null(names(colors)) || !setequal(names(colors), cell_type_order)) {
  stop("annotation.manual_broad_celltype_colors must contain one named color for every cell type.", call. = FALSE)
}
colors <- colors[cell_type_order]

umap_coordinates <- as.data.frame(Seurat::Embeddings(obj, reduction = "umap"))
umap_coordinates$mouse_broad_cell_type <- as.character(obj$mouse_broad_cell_type)
label_coordinates <- stats::aggregate(
  umap_coordinates[, c("UMAP_1", "UMAP_2"), drop = FALSE],
  by = list(mouse_broad_cell_type = umap_coordinates$mouse_broad_cell_type),
  FUN = stats::median
)

p_celltype <- Seurat::DimPlot(
  obj, reduction = "umap", group.by = "mouse_broad_cell_type",
  cols = colors, label = FALSE, raster = FALSE
) +
  ggrepel::geom_label_repel(
    data = label_coordinates,
    ggplot2::aes(x = UMAP_1, y = UMAP_2, label = mouse_broad_cell_type),
    inherit.aes = FALSE,
    seed = 5826,
    size = 3.4,
    box.padding = 0.6,
    point.padding = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf,
    force = 3,
    max.time = 5,
    fill = "white",
    alpha = 0.85,
    label.size = 0.2,
    segment.color = "grey45"
  ) +
  ggplot2::labs(
    title = "Mouse broad cell types: manually approved cluster-level annotation",
    subtitle = "One label per 05d cluster; no cell-level classifier was used",
    color = "Broad cell type"
  ) +
  ggplot2::theme_bw(base_size = 11)

p_cluster <- Seurat::DimPlot(
  obj, reduction = "umap", group.by = "seurat_clusters",
  label = TRUE, repel = TRUE, raster = FALSE
) +
  ggplot2::labs(
    title = "Mouse 05d clusters retained for broad cell-type annotation",
    subtitle = "Cluster IDs provide the audit trail for the manually approved mapping"
  ) +
  ggplot2::theme_bw(base_size = 11)

ggplot2::ggsave(file.path(out_dir, "umap_mouse_broad_celltype.pdf"), p_celltype, width = 11, height = 8, limitsize = FALSE)
ggplot2::ggsave(file.path(out_dir, "umap_mouse_broad_celltype.png"), p_celltype, width = 11, height = 8, dpi = 300, limitsize = FALSE)
ggplot2::ggsave(file.path(out_dir, "umap_mouse_clusters_annotation_audit.pdf"), p_cluster, width = 9, height = 7, limitsize = FALSE)
ggplot2::ggsave(file.path(out_dir, "umap_mouse_clusters_annotation_audit.png"), p_cluster, width = 9, height = 7, dpi = 300, limitsize = FALSE)

utils::write.csv(cluster_annotation, file.path(out_dir, "cluster_broad_celltype_annotation.csv"), row.names = FALSE)
utils::write.csv(celltype_counts, file.path(out_dir, "broad_celltype_counts.csv"), row.names = FALSE)
utils::write.csv(sample_celltype, file.path(out_dir, "sample_by_broad_celltype_counts.csv"), row.names = FALSE)
utils::write.csv(annotation, file.path(out_dir, basename(annotation_file)), row.names = FALSE)

obj@misc$mouse05_broad_celltype_annotation_contract <- list(
  annotation_level = "cluster",
  annotation_method = "manually approved mapping after literature-marker and cluster-DEG review",
  mapping_source = normalizePath(annotation_file, mustWork = TRUE),
  review_completion = normalizePath(review_completion, mustWork = TRUE),
  cluster_count = length(object_clusters),
  cell_type_count = length(cell_type_order),
  human_features_retained_in_RNA = human_feature_n,
  cell_level_classifier_used = FALSE,
  cluster_uniformity_check_passed = TRUE
)

output_rds <- file.path(out_dir, "mouse_clustered_broad_celltype_annotated.rds")
saveRDS(obj, output_rds, compress = FALSE)

input_files <- c(
  input_rds, review_completion, review_manifest, normalizePath(config_path, mustWork = TRUE),
  script_path, annotation_file
)
input_manifest <- data.frame(
  role = c("05d_clustered_object", "05e_review_completion", "05e_output_manifest", "analysis_config", "analysis_script", "approved_cluster_mapping"),
  path = normalizePath(input_files, mustWork = TRUE),
  size_bytes = file.info(input_files)$size,
  sha256 = vapply(input_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(input_manifest, file.path(out_dir, "input_manifest_sha256.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05f_broad_celltype.txt"))
writeLines(
  c(
    "05f manually approved broad cell-type annotation completed.",
    paste0("Cells: ", ncol(obj)),
    paste0("Clusters: ", length(object_clusters)),
    paste0("Broad cell types: ", length(cell_type_order)),
    paste0("Human features in RNA assay: ", human_feature_n),
    "Annotation level: cluster only.",
    "Exactly one broad cell-type label was applied per 05d cluster.",
    "No cell-level classifier was used.",
    "The ambiguous mixed-myeloid cluster was retained and not forced into a lineage.",
    "No immune state or subtype labels were applied in this step."
  ),
  file.path(out_dir, "completion.txt")
)

output_files <- setdiff(list.files(out_dir, full.names = TRUE), file.path(out_dir, "output_manifest_sha256.csv"))
output_manifest <- data.frame(
  file = basename(output_files),
  size_bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(output_manifest, file.path(out_dir, "output_manifest_sha256.csv"), row.names = FALSE)

message("[05f] Manually approved broad cell-type annotation completed. Output: ", out_dir)
