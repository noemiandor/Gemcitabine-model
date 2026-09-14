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

zscore_safe <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(rep(NA_real_, length(x)))
  finite_x <- x[is.finite(x)]
  if (length(finite_x) < 2L || stats::sd(finite_x) == 0) return(ifelse(is.finite(x), 0, NA_real_))
  z <- (x - mean(finite_x)) / stats::sd(finite_x)
  pmax(-2.5, pmin(2.5, z))
}

validate_marker_resource <- function(x, resource_name) {
  required <- c(
    "annotation_level", "parent_cell_type", "candidate_cell_type", "plot_panel",
    "marker_gene", "marker_rank", "marker_role", "expected_in_nsg", "reference",
    "pmid", "doi", "marker_rationale", "caveat"
  )
  missing_columns <- setdiff(required, colnames(x))
  if (length(missing_columns) > 0L) {
    stop(resource_name, " lacks required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  for (column in c("candidate_cell_type", "plot_panel", "marker_gene", "marker_role", "reference")) {
    if (any(!nzchar(trimws(as.character(x[[column]]))))) stop(resource_name, " has blank ", column, ".", call. = FALSE)
  }
  if (any(tolower(x$marker_role) != "positive")) stop(resource_name, " contains a non-positive marker row.", call. = FALSE)
  marker_counts <- table(x$candidate_cell_type)
  invalid_counts <- marker_counts[marker_counts < 3L | marker_counts > 5L]
  if (length(invalid_counts) > 0L) {
    stop(resource_name, " must contain 3-5 markers per candidate; invalid: ",
         paste(names(invalid_counts), invalid_counts, sep = "=", collapse = ", "), call. = FALSE)
  }
  key <- paste(toupper(x$candidate_cell_type), toupper(x$marker_gene), sep = "||")
  if (anyDuplicated(key)) stop(resource_name, " contains duplicate candidate/marker pairs.", call. = FALSE)
  if (any(!is.finite(as.numeric(x$marker_rank)))) stop(resource_name, " has invalid marker_rank values.", call. = FALSE)
  x$marker_rank <- as.integer(x$marker_rank)
  x
}

match_resource_features <- function(resource, assay_features) {
  feature_upper <- toupper(assay_features)
  if (anyDuplicated(feature_upper)) stop("RNA assay has case-insensitive duplicate feature names.", call. = FALSE)
  resource$gene_upper <- toupper(trimws(resource$marker_gene))
  resource$matched_feature <- assay_features[match(resource$gene_upper, feature_upper)]
  resource$is_detected_in_rna_assay <- !is.na(resource$matched_feature)
  resource
}

summarize_marker_expression <- function(resource, expression_matrix, clusters, positive_degs) {
  cluster_ids <- cluster_order(clusters)
  found <- unique(resource[resource$is_detected_in_rna_assay, c("marker_gene", "gene_upper", "matched_feature"), drop = FALSE])

  expression_rows <- vector("list", length(cluster_ids))
  for (i in seq_along(cluster_ids)) {
    cluster_id <- cluster_ids[[i]]
    cell_index <- which(as.character(clusters) == cluster_id)
    if (nrow(found) > 0L) {
      z <- expression_matrix[found$matched_feature, cell_index, drop = FALSE]
      expression_rows[[i]] <- data.frame(
        cluster = cluster_id,
        marker_gene = found$marker_gene,
        gene_upper = found$gene_upper,
        matched_feature = found$matched_feature,
        average_log_normalized_expression = as.numeric(Matrix::rowMeans(z)),
        percent_cells_expressing = 100 * as.numeric(Matrix::rowMeans(z > 0)),
        stringsAsFactors = FALSE
      )
    } else {
      expression_rows[[i]] <- data.frame()
    }
  }
  observed <- do.call(rbind, expression_rows)
  if (is.null(observed)) observed <- data.frame()
  if (nrow(observed) > 0L) {
    observed$scaled_average_expression <- ave(
      observed$average_log_normalized_expression, observed$gene_upper,
      FUN = zscore_safe
    )
  }

  resource_key <- unique(resource[, c(
    "annotation_level", "parent_cell_type", "candidate_cell_type", "plot_panel",
    "marker_gene", "marker_rank", "marker_role", "expected_in_nsg", "reference",
    "pmid", "doi", "marker_rationale", "caveat", "gene_upper", "matched_feature",
    "is_detected_in_rna_assay"
  ), drop = FALSE])
  grid <- merge(
    expand.grid(cluster = cluster_ids, marker_gene = unique(resource_key$marker_gene), stringsAsFactors = FALSE),
    resource_key, by = "marker_gene", all.x = TRUE, sort = FALSE
  )
  if (nrow(observed) > 0L) {
    grid <- merge(
      grid,
      observed[, c("cluster", "marker_gene", "average_log_normalized_expression", "percent_cells_expressing", "scaled_average_expression"), drop = FALSE],
      by = c("cluster", "marker_gene"), all.x = TRUE, sort = FALSE
    )
  } else {
    grid$average_log_normalized_expression <- NA_real_
    grid$percent_cells_expressing <- NA_real_
    grid$scaled_average_expression <- NA_real_
  }
  grid$percent_cells_expressing[is.na(grid$percent_cells_expressing)] <- 0

  deg_keep <- positive_degs[, c(
    "cluster", "gene_upper", "deg_avg_log2FC", "deg_pct_in", "deg_pct_out", "deg_p_adj"
  ), drop = FALSE]
  grid <- merge(grid, deg_keep, by = c("cluster", "gene_upper"), all.x = TRUE, sort = FALSE)
  grid$is_positive_cluster_deg <- !is.na(grid$deg_avg_log2FC)

  candidate_levels <- unique(resource$candidate_cell_type)
  marker_keys <- unique(paste(resource$candidate_cell_type, resource$marker_gene, sep = "|||"))
  grid$candidate_cell_type <- factor(grid$candidate_cell_type, levels = candidate_levels)
  grid$cluster <- factor(grid$cluster, levels = cluster_ids)
  grid$marker_key <- factor(
    paste(as.character(grid$candidate_cell_type), grid$marker_gene, sep = "|||"),
    levels = rev(marker_keys)
  )
  grid
}

make_dotplot <- function(plot_data, title, subtitle) {
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = cluster, y = marker_key, size = percent_cells_expressing, color = scaled_average_expression)
  ) +
    ggplot2::geom_point(alpha = 0.92) +
    ggplot2::facet_grid(candidate_cell_type ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_y_discrete(labels = function(x) sub("^.*\\|\\|\\|", "", x)) +
    ggplot2::scale_size_continuous(name = "% expressing", range = c(0, 7), limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    ggplot2::scale_color_gradient2(
      name = "Mean expression\n(row z-score)", low = "#3B6FB6", mid = "#F2F2F2", high = "#D95F02",
      midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish, na.value = "#F7F7F7"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste(strwrap(subtitle, width = 118), collapse = "\n"),
      x = "05d cluster", y = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, face = "bold"),
      panel.grid.major = ggplot2::element_line(color = "#E6E6E6", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title.position = "plot",
      legend.position = "right"
    )
}

write_plot_pair <- function(plot, stem, out_dir, width, height) {
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = width, height = height, limitsize = FALSE)
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300, limitsize = FALSE)
}

script_dir <- resolve_script_dir()
script_path <- file.path(script_dir, "05e_mouse_broad_annotation.R")
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")
required_packages <- c("yaml", "Seurat", "Matrix", "ggplot2", "scales", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
input_rds <- file.path(results_root, "05d_clustering", "mouse_integrated_sct_cca_clustered.rds")
deg_file <- file.path(results_root, "05d_clustering", "cluster_positive_markers.csv")
resource_dir <- file.path(script_dir, "annotation_resources")
broad_resource_file <- file.path(resource_dir, "mouse_broad_celltype_markers.csv")
immune_resource_file <- file.path(resource_dir, "mouse_immune_subtype_markers.csv")
out_dir <- file.path(results_root, "05e_broad_annotation_marker_review")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(input_rds, deg_file, broad_resource_file, immune_resource_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) stop("Missing required input(s): ", paste(missing_inputs, collapse = ", "), call. = FALSE)

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("05d input is not a Seurat object.", call. = FALSE)
if (!"seurat_clusters" %in% colnames(obj@meta.data)) stop("05d object lacks seurat_clusters.", call. = FALSE)
if (!"RNA" %in% names(obj@assays)) stop("05d object lacks an RNA assay.", call. = FALSE)
human_feature_n <- sum(startsWith(rownames(obj[["RNA"]]), as.character(cfg$species$human_prefix)))
if (human_feature_n != 0L) stop("Human features detected before marker review: ", human_feature_n, call. = FALSE)
if (!"umap" %in% names(obj@reductions)) stop("05d object lacks the UMAP reduction.", call. = FALSE)

broad_resource <- validate_marker_resource(
  utils::read.csv(broad_resource_file, stringsAsFactors = FALSE, check.names = FALSE, na.strings = character()),
  basename(broad_resource_file)
)
immune_resource <- validate_marker_resource(
  utils::read.csv(immune_resource_file, stringsAsFactors = FALSE, check.names = FALSE, na.strings = character()),
  basename(immune_resource_file)
)
broad_resource <- match_resource_features(broad_resource, rownames(obj[["RNA"]]))
immune_resource <- match_resource_features(immune_resource, rownames(obj[["RNA"]]))

degs <- utils::read.csv(deg_file, stringsAsFactors = FALSE, check.names = FALSE)
fc_column <- intersect(c("avg_log2FC", "avg_logFC"), colnames(degs))
if (!all(c("cluster", "gene") %in% colnames(degs)) || length(fc_column) != 1L) {
  stop("05d marker table must contain cluster, gene, and exactly one avg_log2FC/avg_logFC column.", call. = FALSE)
}
if (nrow(degs) == 0L) stop("05d marker table is empty.", call. = FALSE)
degs$cluster <- as.character(degs$cluster)
degs$gene_upper <- toupper(trimws(as.character(degs$gene)))
degs$deg_avg_log2FC <- as.numeric(degs[[fc_column]])
degs$deg_pct_in <- if ("pct.1" %in% colnames(degs)) as.numeric(degs$pct.1) else NA_real_
degs$deg_pct_out <- if ("pct.2" %in% colnames(degs)) as.numeric(degs$pct.2) else NA_real_
degs$deg_p_adj <- if ("p_val_adj" %in% colnames(degs)) as.numeric(degs$p_val_adj) else NA_real_
degs <- degs[nzchar(degs$gene_upper) & is.finite(degs$deg_avg_log2FC), , drop = FALSE]
degs <- degs[order(degs$cluster, -degs$deg_avg_log2FC, degs$gene_upper), , drop = FALSE]
degs <- degs[!duplicated(paste(degs$cluster, degs$gene_upper, sep = "||")), , drop = FALSE]

clusters <- as.character(obj$seurat_clusters)
cluster_ids <- cluster_order(clusters)
missing_deg_clusters <- setdiff(cluster_ids, unique(degs$cluster))
if (length(missing_deg_clusters) > 0L) stop("No positive DEGs for clusters: ", paste(missing_deg_clusters, collapse = ", "), call. = FALSE)

Seurat::DefaultAssay(obj) <- "RNA"
expression_matrix <- Seurat::GetAssayData(obj, assay = "RNA", slot = "data")
if (nrow(expression_matrix) == 0L || ncol(expression_matrix) != ncol(obj)) stop("RNA data slot is empty or misaligned.", call. = FALSE)

broad_expression <- summarize_marker_expression(broad_resource, expression_matrix, clusters, degs)
immune_expression <- summarize_marker_expression(immune_resource, expression_matrix, clusters, degs)

broad_audit <- broad_resource[, c(
  "candidate_cell_type", "marker_gene", "marker_rank", "matched_feature", "is_detected_in_rna_assay",
  "reference", "pmid", "doi", "marker_rationale", "caveat"
), drop = FALSE]
immune_audit <- immune_resource[, c(
  "annotation_level", "parent_cell_type", "candidate_cell_type", "plot_panel", "marker_gene", "marker_rank",
  "expected_in_nsg", "matched_feature", "is_detected_in_rna_assay", "reference", "pmid", "doi", "marker_rationale", "caveat"
), drop = FALSE]

top_n <- as.integer(cfg$annotation$top_cluster_markers)
if (!is.finite(top_n) || top_n < 10L) stop("annotation.top_cluster_markers must be at least 10.", call. = FALSE)
top_deg_list <- lapply(cluster_ids, function(cluster_id) {
  z <- degs[degs$cluster == cluster_id, , drop = FALSE]
  z <- utils::head(z, top_n)
  data.frame(
    cluster = cluster_id,
    marker_rank = seq_len(nrow(z)),
    gene = z$gene,
    avg_log2FC = z$deg_avg_log2FC,
    pct_in = z$deg_pct_in,
    pct_out = z$deg_pct_out,
    p_val_adj = z$deg_p_adj,
    stringsAsFactors = FALSE
  )
})
top_degs <- do.call(rbind, top_deg_list)
rownames(top_degs) <- NULL

cluster_counts <- table(factor(clusters, levels = cluster_ids))
review_template <- data.frame(
  cluster = cluster_ids,
  n_cells = as.integer(cluster_counts),
  top_10_positive_degs = vapply(cluster_ids, function(cluster_id) {
    paste(utils::head(top_degs$gene[top_degs$cluster == cluster_id], 10L), collapse = ";")
  }, character(1)),
  reviewer_broad_cell_type = rep("", length(cluster_ids)),
  reviewer_confidence = rep("", length(cluster_ids)),
  reviewer_notes = rep("", length(cluster_ids)),
  review_status = rep("pending_manual_review", length(cluster_ids)),
  stringsAsFactors = FALSE
)

utils::write.csv(broad_resource[, setdiff(colnames(broad_resource), c("gene_upper", "matched_feature", "is_detected_in_rna_assay")), drop = FALSE],
                 file.path(out_dir, basename(broad_resource_file)), row.names = FALSE, na = "")
utils::write.csv(immune_resource[, setdiff(colnames(immune_resource), c("gene_upper", "matched_feature", "is_detected_in_rna_assay")), drop = FALSE],
                 file.path(out_dir, basename(immune_resource_file)), row.names = FALSE, na = "")
utils::write.csv(broad_audit, file.path(out_dir, "mouse_broad_marker_detection_audit.csv"), row.names = FALSE, na = "")
utils::write.csv(immune_audit, file.path(out_dir, "mouse_immune_marker_detection_audit.csv"), row.names = FALSE, na = "")
utils::write.csv(broad_expression, file.path(out_dir, "mouse_broad_marker_cluster_expression.csv"), row.names = FALSE, na = "")
utils::write.csv(immune_expression, file.path(out_dir, "mouse_immune_marker_cluster_expression.csv"), row.names = FALSE, na = "")
utils::write.csv(top_degs, file.path(out_dir, "cluster_top_positive_degs.csv"), row.names = FALSE, na = "")
utils::write.csv(review_template, file.path(out_dir, "cluster_manual_review_template.csv"), row.names = FALSE, na = "")

cluster_plot <- Seurat::DimPlot(
  obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, raster = FALSE
) +
  ggplot2::labs(
    title = "Mouse cells: 05d cluster map for manual annotation review",
    subtitle = "Cluster IDs only; no cell-type labels have been assigned in this stage"
  ) +
  ggplot2::theme_bw(base_size = 11)
write_plot_pair(cluster_plot, "umap_mouse_clusters_manual_review", out_dir, width = 9, height = 7)

common_subtitle <- paste0(
  length(cluster_ids), " clusters; dot size = percent expressing; color = cluster mean expression scaled within each gene. ",
  "Markers support manual cluster-level review and do not assign labels."
)
broad_plot <- make_dotplot(
  broad_expression,
  "Literature-curated mouse broad cell-type markers across 05d clusters",
  paste0(common_subtitle, " Mammary-tissue markers were prioritized.")
)
write_plot_pair(broad_plot, "dotplot_mouse_broad_celltype_markers", out_dir, width = 13, height = 16)

panel_order <- unique(immune_resource$plot_panel)
for (panel_name in panel_order) {
  panel_data <- immune_expression[immune_expression$plot_panel == panel_name, , drop = FALSE]
  panel_data$candidate_cell_type <- droplevels(panel_data$candidate_cell_type)
  type_n <- length(unique(as.character(panel_data$candidate_cell_type)))
  plot_height <- max(7, 1.55 * type_n + 3.5)
  panel_plot <- make_dotplot(
    panel_data,
    paste0("Literature-curated mouse immune markers: ", gsub("_", " ", panel_name)),
    paste0(
      common_subtitle,
      if (panel_name == "lymphoid_qc") " NSG mice should lack mature mouse T/B/NK cells; this panel is a QC diagnostic." else ""
    )
  )
  write_plot_pair(panel_plot, paste0("dotplot_mouse_immune_markers_", panel_name), out_dir, width = 13, height = plot_height)
}

input_files <- c(input_rds, deg_file, normalizePath(config_path, mustWork = TRUE), script_path, broad_resource_file, immune_resource_file)
input_manifest <- data.frame(
  role = c("05d_clustered_object", "05d_positive_degs", "analysis_config", "analysis_script", "broad_marker_resource", "immune_marker_resource"),
  path = normalizePath(input_files, mustWork = TRUE),
  size_bytes = file.info(input_files)$size,
  sha256 = vapply(input_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(input_manifest, file.path(out_dir, "input_manifest_sha256.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05e_marker_review.txt"))
writeLines(
  c(
    "05e literature-marker review artifacts completed.",
    paste0("Cells: ", ncol(obj)),
    paste0("Clusters: ", length(cluster_ids)),
    paste0("Human features in RNA assay: ", human_feature_n),
    paste0("Broad marker rows: ", nrow(broad_resource)),
    paste0("Broad markers detected: ", sum(broad_resource$is_detected_in_rna_assay), "/", nrow(broad_resource)),
    paste0("Immune marker rows: ", nrow(immune_resource)),
    paste0("Immune markers detected: ", sum(immune_resource$is_detected_in_rna_assay), "/", nrow(immune_resource)),
    "Annotation level: cluster only.",
    "No automated cell-type label was calculated.",
    "No annotation label was added to Seurat metadata.",
    "No annotated Seurat object was saved.",
    "Manual review is required before any cluster label is applied.",
    "Downstream steps 05f and later were not invoked."
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

message("[05e] Manual marker-review artifacts completed. Output: ", out_dir)
