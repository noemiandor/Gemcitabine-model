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
  for (column in c("annotation_level", "parent_cell_type", "candidate_cell_type", "plot_panel", "marker_gene", "marker_role", "reference")) {
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

summarize_marker_expression <- function(resource, expression_matrix, clusters, positive_degs, cluster_ids) {
  found <- unique(resource[resource$is_detected_in_rna_assay, c("marker_gene", "gene_upper", "matched_feature"), drop = FALSE])
  expression_rows <- vector("list", length(cluster_ids))
  for (i in seq_along(cluster_ids)) {
    cluster_id <- cluster_ids[[i]]
    cell_index <- which(as.character(clusters) == cluster_id)
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
  }
  observed <- do.call(rbind, expression_rows)
  observed$scaled_average_expression <- ave(
    observed$average_log_normalized_expression, observed$gene_upper,
    FUN = zscore_safe
  )

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
  grid <- merge(
    grid,
    observed[, c("cluster", "marker_gene", "average_log_normalized_expression", "percent_cells_expressing", "scaled_average_expression"), drop = FALSE],
    by = c("cluster", "marker_gene"), all.x = TRUE, sort = FALSE
  )
  grid$percent_cells_expressing[is.na(grid$percent_cells_expressing)] <- 0

  deg_keep <- positive_degs[, c("cluster", "gene_upper", "deg_avg_log2FC", "deg_pct_in", "deg_pct_out", "deg_p_adj"), drop = FALSE]
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

summarize_candidate_panels <- function(expression_table) {
  keys <- unique(expression_table[, c("cluster", "annotation_level", "parent_cell_type", "candidate_cell_type", "plot_panel"), drop = FALSE])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    z <- expression_table[
      as.character(expression_table$cluster) == as.character(keys$cluster[[i]]) &
        expression_table$annotation_level == keys$annotation_level[[i]] &
        expression_table$parent_cell_type == keys$parent_cell_type[[i]] &
        as.character(expression_table$candidate_cell_type) == as.character(keys$candidate_cell_type[[i]]) &
        expression_table$plot_panel == keys$plot_panel[[i]],
      , drop = FALSE
    ]
    z <- z[order(z$marker_rank), , drop = FALSE]
    data.frame(
      cluster = as.character(keys$cluster[[i]]),
      annotation_level = keys$annotation_level[[i]],
      parent_cell_type = keys$parent_cell_type[[i]],
      candidate_cell_type = as.character(keys$candidate_cell_type[[i]]),
      plot_panel = keys$plot_panel[[i]],
      marker_n = nrow(z),
      markers_expressed_ge_10pct = sum(z$percent_cells_expressing >= 10, na.rm = TRUE),
      positive_cluster_deg_marker_n = sum(z$is_positive_cluster_deg, na.rm = TRUE),
      mean_percent_cells_expressing = mean(z$percent_cells_expressing, na.rm = TRUE),
      mean_log_normalized_expression = mean(z$average_log_normalized_expression, na.rm = TRUE),
      marker_expression = paste0(z$marker_gene, "=", sprintf("%.1f%%", z$percent_cells_expressing), collapse = ";"),
      positive_deg_markers = paste(z$marker_gene[z$is_positive_cluster_deg], collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_dotplot <- function(plot_data, title, subtitle, x_title = "Existing 05d cluster") {
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
      x = x_title, y = NULL
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

prepare_combined_panel <- function(expression_table, lineage_name, state_parent) {
  lineage_rows <- expression_table[
    expression_table$annotation_level == "lineage" &
      as.character(expression_table$candidate_cell_type) == lineage_name,
    , drop = FALSE
  ]
  state_rows <- expression_table[
    expression_table$annotation_level == "state" & expression_table$parent_cell_type == state_parent,
    , drop = FALSE
  ]
  lineage_rows$candidate_cell_type <- paste0("LINEAGE | ", lineage_name)
  state_rows$candidate_cell_type <- paste0("STATE | ", as.character(state_rows$candidate_cell_type))
  out <- rbind(lineage_rows, state_rows)
  candidate_levels <- unique(c(paste0("LINEAGE | ", lineage_name), unique(as.character(state_rows$candidate_cell_type))))
  marker_keys <- unique(paste(as.character(out$candidate_cell_type), out$marker_gene, sep = "|||"))
  out$candidate_cell_type <- factor(out$candidate_cell_type, levels = candidate_levels)
  out$marker_key <- factor(paste(as.character(out$candidate_cell_type), out$marker_gene, sep = "|||"), levels = rev(marker_keys))
  out
}

script_dir <- resolve_script_dir()
script_path <- file.path(script_dir, "05g_review_existing_mouse_immune_clusters.R")
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")
required_packages <- c("yaml", "Seurat", "Matrix", "ggplot2", "scales", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
input_rds <- file.path(results_root, "05f_broad_celltype_annotation", "mouse_clustered_broad_celltype_annotated.rds")
deg_file <- file.path(results_root, "05d_clustering", "cluster_positive_markers.csv")
resource_file <- file.path(script_dir, as.character(cfg$annotation$immune_marker_resource))
out_dir <- file.path(results_root, as.character(cfg$immune_fine_annotation$existing_review_output_directory))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(input_rds, deg_file, resource_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) stop("Missing required input(s): ", paste(missing_inputs, collapse = ", "), call. = FALSE)

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("05f input is not a Seurat object.", call. = FALSE)
required_meta <- c(
  "sample", "barcode_raw", "seurat_clusters", "mouse_fraction", "scDblFinder.class",
  "mouse_broad_cell_type", "mouse_broad_annotation_level", "mouse_broad_review_status"
)
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0L) stop("05f object lacks metadata: ", paste(missing_meta, collapse = ", "), call. = FALSE)
if (!"RNA" %in% names(obj@assays)) stop("05f object lacks the RNA assay.", call. = FALSE)
if (!"umap" %in% names(obj@reductions)) stop("05f object lacks UMAP.", call. = FALSE)
human_feature_n <- sum(startsWith(rownames(obj[["RNA"]]), as.character(cfg$species$human_prefix)))
if (human_feature_n != 0L) stop("Human features detected before immune marker review: ", human_feature_n, call. = FALSE)
if (any(as.character(obj$mouse_broad_annotation_level) != "cluster")) stop("05f broad annotation is not cluster-level.", call. = FALSE)
if (any(as.character(obj$mouse_broad_review_status) != "manually_approved")) stop("05f broad annotations are not all manually approved.", call. = FALSE)

cluster_vector <- as.character(obj$seurat_clusters)
broad_vector <- as.character(obj$mouse_broad_cell_type)
cluster_broad_n <- tapply(broad_vector, cluster_vector, function(x) length(unique(x)))
if (any(cluster_broad_n != 1L)) stop("05f broad annotation is not uniform within every cluster.", call. = FALSE)

confirmed_types <- as.character(unlist(cfg$immune_fine_annotation$confirmed_broad_cell_types, use.names = FALSE))
candidate_types <- as.character(unlist(cfg$immune_fine_annotation$candidate_broad_cell_types, use.names = FALSE))
if (length(confirmed_types) == 0L) stop("No confirmed broad immune types are configured.", call. = FALSE)
review_mask <- broad_vector %in% c(confirmed_types, candidate_types)
review_cells <- colnames(obj)[review_mask]
review_clusters <- cluster_order(cluster_vector[review_mask])
if (length(review_cells) == 0L || length(review_clusters) == 0L) stop("No cells qualified for existing immune-cluster review.", call. = FALSE)
if (any(as.character(obj$scDblFinder.class[review_mask]) != "singlet")) stop("The immune review universe contains non-singlets.", call. = FALSE)
if (any(as.numeric(obj$mouse_fraction[review_mask]) < as.numeric(cfg$species$mouse_fraction_min), na.rm = TRUE)) {
  stop("The immune review universe violates the mouse-fraction threshold.", call. = FALSE)
}

resource <- validate_marker_resource(
  utils::read.csv(resource_file, stringsAsFactors = FALSE, check.names = FALSE, na.strings = character()),
  basename(resource_file)
)
resource <- match_resource_features(resource, rownames(obj[["RNA"]]))
if (any(!resource$is_detected_in_rna_assay)) {
  stop("Immune marker resource has undetected RNA features: ", paste(resource$marker_gene[!resource$is_detected_in_rna_assay], collapse = ", "), call. = FALSE)
}

degs <- utils::read.csv(deg_file, stringsAsFactors = FALSE, check.names = FALSE)
fc_column <- intersect(c("avg_log2FC", "avg_logFC"), colnames(degs))
if (!all(c("cluster", "gene") %in% colnames(degs)) || length(fc_column) != 1L) {
  stop("05d marker table must contain cluster, gene, and exactly one avg_log2FC/avg_logFC column.", call. = FALSE)
}
degs$cluster <- as.character(degs$cluster)
degs$gene_upper <- toupper(trimws(as.character(degs$gene)))
degs$deg_avg_log2FC <- as.numeric(degs[[fc_column]])
degs$deg_pct_in <- if ("pct.1" %in% colnames(degs)) as.numeric(degs$pct.1) else NA_real_
degs$deg_pct_out <- if ("pct.2" %in% colnames(degs)) as.numeric(degs$pct.2) else NA_real_
degs$deg_p_adj <- if ("p_val_adj" %in% colnames(degs)) as.numeric(degs$p_val_adj) else NA_real_
degs <- degs[nzchar(degs$gene_upper) & is.finite(degs$deg_avg_log2FC), , drop = FALSE]
degs <- degs[order(degs$cluster, -degs$deg_avg_log2FC, degs$gene_upper), , drop = FALSE]
degs <- degs[!duplicated(paste(degs$cluster, degs$gene_upper, sep = "||")), , drop = FALSE]
missing_deg_clusters <- setdiff(review_clusters, unique(degs$cluster))
if (length(missing_deg_clusters) > 0L) stop("No positive DEGs for review clusters: ", paste(missing_deg_clusters, collapse = ", "), call. = FALSE)

Seurat::DefaultAssay(obj) <- "RNA"
expression_matrix <- Seurat::GetAssayData(obj, assay = "RNA", slot = "data")
if (nrow(expression_matrix) == 0L || ncol(expression_matrix) != ncol(obj)) stop("RNA data slot is empty or misaligned.", call. = FALSE)
marker_expression <- summarize_marker_expression(resource, expression_matrix, cluster_vector, degs, review_clusters)
cluster_to_broad <- tapply(broad_vector[review_mask], cluster_vector[review_mask], unique)
marker_expression$source_broad_cell_type <- unname(cluster_to_broad[as.character(marker_expression$cluster)])
panel_summary <- summarize_candidate_panels(marker_expression)
panel_summary$source_broad_cell_type <- unname(cluster_to_broad[panel_summary$cluster])

top_exports <- list()
for (top_n in c(20L, 50L, 100L)) {
  local <- do.call(rbind, lapply(review_clusters, function(cluster_id) {
    z <- utils::head(degs[degs$cluster == cluster_id, , drop = FALSE], top_n)
    data.frame(
      cluster = cluster_id, marker_rank = seq_len(nrow(z)), gene = z$gene,
      avg_log2FC = z$deg_avg_log2FC, pct_in = z$deg_pct_in, pct_out = z$deg_pct_out,
      p_val_adj = z$deg_p_adj, stringsAsFactors = FALSE
    )
  }))
  rownames(local) <- NULL
  top_exports[[as.character(top_n)]] <- local
  utils::write.csv(local, file.path(out_dir, paste0("existing_immune_cluster_top", top_n, "_positive_degs.csv")), row.names = FALSE, na = "")
}

cluster_counts <- table(factor(cluster_vector[review_mask], levels = review_clusters))
cluster_sample_counts <- as.data.frame(table(
  cluster = factor(cluster_vector[review_mask], levels = review_clusters),
  sample = as.character(obj$sample[review_mask])
), stringsAsFactors = FALSE)
colnames(cluster_sample_counts)[3L] <- "n_cells"
cluster_sample_counts <- cluster_sample_counts[cluster_sample_counts$n_cells > 0L, , drop = FALSE]
cluster_totals <- tapply(cluster_sample_counts$n_cells, cluster_sample_counts$cluster, sum)
cluster_sample_counts$fraction_within_cluster <- cluster_sample_counts$n_cells / unname(cluster_totals[as.character(cluster_sample_counts$cluster)])
cluster_summary <- do.call(rbind, lapply(review_clusters, function(cluster_id) {
  sample_z <- cluster_sample_counts[as.character(cluster_sample_counts$cluster) == cluster_id, , drop = FALSE]
  broad_type <- unname(cluster_to_broad[[cluster_id]])
  data.frame(
    cluster = cluster_id,
    n_cells = as.integer(cluster_counts[[cluster_id]]),
    source_broad_cell_type = broad_type,
    review_scope = if (broad_type %in% confirmed_types) "confirmed_immune_broad_type" else "candidate_ambiguous_mixed_myeloid",
    samples_detected = nrow(sample_z),
    dominant_sample = as.character(sample_z$sample[[which.max(sample_z$n_cells)]]),
    max_single_sample_fraction = max(sample_z$fraction_within_cluster),
    top_10_positive_degs = paste(utils::head(top_exports[["100"]]$gene[top_exports[["100"]]$cluster == cluster_id], 10L), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
rownames(cluster_summary) <- NULL

review_template <- cluster_summary
review_template$reviewer_immune_lineage <- ""
review_template$reviewer_lineage_confidence <- ""
review_template$reviewer_immune_state <- ""
review_template$reviewer_state_confidence <- ""
review_template$lineage_support_markers <- ""
review_template$state_support_markers <- ""
review_template$reviewer_notes <- ""
review_template$review_status <- "pending_manual_review"

resource_export <- resource[, setdiff(colnames(resource), c("gene_upper", "matched_feature", "is_detected_in_rna_assay")), drop = FALSE]
utils::write.csv(resource_export, file.path(out_dir, basename(resource_file)), row.names = FALSE, na = "")
utils::write.csv(resource[, c(
  "annotation_level", "parent_cell_type", "candidate_cell_type", "plot_panel", "marker_gene", "marker_rank",
  "expected_in_nsg", "matched_feature", "is_detected_in_rna_assay", "reference", "pmid", "doi", "marker_rationale", "caveat"
), drop = FALSE], file.path(out_dir, "mouse_immune_marker_detection_audit.csv"), row.names = FALSE, na = "")
utils::write.csv(marker_expression, file.path(out_dir, "existing_immune_cluster_marker_expression.csv"), row.names = FALSE, na = "")
utils::write.csv(panel_summary, file.path(out_dir, "existing_immune_cluster_marker_panel_summary.csv"), row.names = FALSE, na = "")
utils::write.csv(cluster_summary, file.path(out_dir, "existing_immune_cluster_review_summary.csv"), row.names = FALSE, na = "")
utils::write.csv(cluster_sample_counts, file.path(out_dir, "existing_immune_cluster_by_sample_counts.csv"), row.names = FALSE, na = "")
utils::write.csv(review_template, file.path(out_dir, "existing_immune_cluster_manual_review_template.csv"), row.names = FALSE, na = "")

review_obj <- subset(obj, cells = review_cells)
review_obj$existing_immune_review_cluster <- factor(as.character(review_obj$seurat_clusters), levels = review_clusters)
review_obj$existing_immune_review_broad_type <- factor(
  as.character(review_obj$mouse_broad_cell_type),
  levels = c(confirmed_types, candidate_types)
)
broad_colors <- unlist(cfg$annotation$manual_broad_celltype_colors, use.names = TRUE)
if (anyNA(broad_colors[levels(review_obj$existing_immune_review_broad_type)])) stop("Broad-type palette is incomplete for immune review.", call. = FALSE)
p_clusters <- Seurat::DimPlot(
  review_obj, reduction = "umap", group.by = "existing_immune_review_cluster",
  label = TRUE, repel = TRUE, raster = FALSE
) +
  ggplot2::labs(
    title = "Existing 05d immune-related clusters selected for fine-annotation review",
    subtitle = "Cluster IDs only; no fine immune label has been assigned"
  ) +
  ggplot2::theme_bw(base_size = 11)
p_broad <- Seurat::DimPlot(
  review_obj, reduction = "umap", group.by = "existing_immune_review_broad_type",
  cols = unname(broad_colors[levels(review_obj$existing_immune_review_broad_type)]),
  label = TRUE, repel = TRUE, raster = FALSE
) +
  ggplot2::labs(
    title = "05f broad labels defining the existing immune-review universe",
    subtitle = "Ambiguous mixed myeloid is retained as a review candidate, not a confirmed immune lineage",
    color = "05f broad type"
  ) +
  ggplot2::theme_bw(base_size = 11)
write_plot_pair(p_clusters, "umap_existing_immune_clusters_review", out_dir, width = 10, height = 8)
write_plot_pair(p_broad, "umap_existing_immune_review_source_broad_types", out_dir, width = 11, height = 8)

common_subtitle <- paste0(
  length(review_clusters), " existing clusters and ", length(review_cells),
  " cells; dot size = percent expressing; color = cluster mean expression scaled within each gene. ",
  "The plots support manual cluster-level review and do not assign labels."
)
for (panel_name in c("myeloid_lineage", "granulocyte_lineage", "lymphoid_qc")) {
  panel_data <- marker_expression[marker_expression$plot_panel == panel_name, , drop = FALSE]
  panel_data$candidate_cell_type <- droplevels(panel_data$candidate_cell_type)
  panel_plot <- make_dotplot(
    panel_data,
    paste0("Existing immune clusters: ", gsub("_", " ", panel_name), " markers"),
    paste0(
      common_subtitle,
      if (panel_name == "lymphoid_qc") " NSG mice should lack mature mouse T/B/NK cells; these markers are QC signals only." else ""
    )
  )
  type_n <- length(unique(as.character(panel_data$candidate_cell_type)))
  write_plot_pair(panel_plot, paste0("dotplot_existing_immune_", panel_name), out_dir, width = 13, height = max(8, 1.55 * type_n + 3.5))
}

neutrophil_panel <- prepare_combined_panel(marker_expression, "Neutrophil", "Neutrophil")
neutrophil_plot <- make_dotplot(
  neutrophil_panel,
  "Existing immune clusters: Neutrophil lineage core and candidate states",
  paste0(common_subtitle, " A Neutrophil state is interpretable only after the Neutrophil lineage core is supported.")
)
neutrophil_panel_n <- length(unique(as.character(neutrophil_panel$candidate_cell_type)))
write_plot_pair(
  neutrophil_plot, "dotplot_existing_immune_neutrophil_lineage_and_states", out_dir,
  width = 13, height = max(15, 1.6 * neutrophil_panel_n + 4)
)

macrophage_panel <- prepare_combined_panel(marker_expression, "Macrophage", "Macrophage")
macrophage_plot <- make_dotplot(
  macrophage_panel,
  "Existing immune clusters: Macrophage lineage core and candidate states",
  paste0(common_subtitle, " A Macrophage state is interpretable only after the Macrophage lineage core is supported.")
)
write_plot_pair(macrophage_plot, "dotplot_existing_immune_macrophage_lineage_and_states", out_dir, width = 13, height = 12)

input_files <- c(input_rds, deg_file, normalizePath(config_path, mustWork = TRUE), script_path, resource_file)
input_manifest <- data.frame(
  role = c("05f_manually_approved_broad_object", "05d_positive_cluster_degs", "analysis_config", "analysis_script", "immune_marker_resource"),
  path = normalizePath(input_files, mustWork = TRUE),
  size_bytes = file.info(input_files)$size,
  sha256 = vapply(input_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(input_manifest, file.path(out_dir, "input_manifest_sha256.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05g_existing_immune_review.txt"))
writeLines(
  c(
    "05g existing immune-cluster marker review completed.",
    paste0("Full 05f mouse cells: ", ncol(obj)),
    paste0("Review cells: ", length(review_cells)),
    paste0("Review clusters: ", length(review_clusters), " [", paste(review_clusters, collapse = ","), "]"),
    paste0("Confirmed broad immune cells: ", sum(broad_vector %in% confirmed_types)),
    paste0("Candidate ambiguous mixed-myeloid cells: ", sum(broad_vector %in% candidate_types)),
    paste0("Human features in RNA assay: ", human_feature_n),
    paste0("Immune markers detected: ", sum(resource$is_detected_in_rna_assay), "/", nrow(resource)),
    "Annotation level under review: existing 05d cluster only.",
    "Lineage markers must be established before any lineage-specific state is considered.",
    "NSG T/B/NK markers are QC diagnostics and were not treated as expected lineages.",
    "No automated immune lineage or state label was calculated.",
    "No fine-annotation label was added to Seurat metadata.",
    "No annotated Seurat object was saved.",
    "Manual review is required before 05h can apply any cluster label.",
    "05h and later steps were not invoked."
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

message("[05g] Existing immune-cluster review artifacts completed. Output: ", out_dir)
