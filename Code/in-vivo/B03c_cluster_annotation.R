#!/usr/bin/env Rscript

resolve_in_vivo_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_match <- grep("--file=", cmd_args, value = TRUE)
  candidate_files <- character(0)
  if (length(file_match) > 0) {
    candidate_files <- c(candidate_files, sub("--file=", "", file_match[1]))
  }

  frame_files <- vapply(
    sys.frames(),
    function(x) {
      if (!is.null(x$ofile)) x$ofile else NA_character_
    },
    character(1)
  )
  candidate_files <- c(candidate_files, frame_files[!is.na(frame_files)])

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(active_path)) candidate_files <- c(candidate_files, active_path)
  }

  candidate_files <- candidate_files[!is.na(candidate_files) & nzchar(candidate_files)]
  candidate_dirs <- unique(dirname(normalizePath(candidate_files, mustWork = FALSE)))
  cwd <- normalizePath(getwd(), mustWork = FALSE)
  cwd_parts <- strsplit(cwd, .Platform$file.sep, fixed = TRUE)[[1]]
  parent_dirs <- vapply(
    seq_along(cwd_parts),
    function(i) {
      paste(c(cwd_parts[seq_len(length(cwd_parts) - i + 1)]), collapse = .Platform$file.sep)
    },
    character(1)
  )
  parent_dirs <- parent_dirs[nzchar(parent_dirs)]
  parent_dirs <- if (grepl("^/", cwd)) paste0("/", sub("^/+", "", parent_dirs)) else parent_dirs
  candidate_dirs <- unique(c(
    candidate_dirs,
    cwd,
    file.path(cwd, "Code", "in-vivo"),
    parent_dirs,
    file.path(parent_dirs, "Code", "in-vivo")
  ))

  utils_paths <- file.path(candidate_dirs, "Utils.R")
  hit <- candidate_dirs[file.exists(utils_paths)]
  if (length(hit) > 0) return(normalizePath(hit[1], mustWork = TRUE))
  stop("Cannot locate Code/in-vivo/Utils.R from script path or working directory: ", cwd, call. = FALSE)
}

script_dir <- resolve_in_vivo_script_dir()
source(file.path(script_dir, "Utils.R"))
config <- load_in_vivo_config(file.path(script_dir, "in_vivo_config.yaml"))
results_root <- get_results_root(config)

required_packages <- c(
  "Seurat",
  "dplyr",
  "readr",
  "ggplot2",
  "msigdbr",
  "pheatmap"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(msigdbr)
  library(pheatmap)
})

set.seed(12345)

get_env_scalar <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(trimws(value))) return(default)
  value
}

parse_optional_expected_integer <- function(x, source_name) {
  if (is.null(x) || length(x) == 0) return(NA_integer_)
  x <- trimws(as.character(x)[1])
  if (!nzchar(x)) return(NA_integer_)
  value <- suppressWarnings(as.integer(x))
  if (is.na(value) || value < 1L) {
    stop(source_name, " must be empty or a positive integer.", call. = FALSE)
  }
  value
}

annotation_script_name <- get_env_scalar("CLUSTER_ANNOTATION_SCRIPT_NAME", "03c_cluster_annotation")
analysis_cluster_col <- get_env_scalar("CLUSTER_ANNOTATION_CLUSTER_COL", "cluster_final")
analysis_cluster_label <- get_env_scalar("CLUSTER_ANNOTATION_CLUSTER_LABEL", analysis_cluster_col)
analysis_cluster_stub <- sanitize_path_component(analysis_cluster_label)[1]

input_rds <- file.path(
  results_root,
  "03b_manual_cluster_merge",
  "objects",
  "integrated_sct_cca_seurat_final_manual_merge.rds"
)
input_rds <- get_env_scalar("CLUSTER_ANNOTATION_INPUT_RDS", input_rds)
output_root <- get_env_scalar(
  "CLUSTER_ANNOTATION_OUTPUT_ROOT",
  file.path(results_root, annotation_script_name)
)
expected_cluster_count <- parse_optional_expected_integer(
  Sys.getenv("CLUSTER_ANNOTATION_EXPECTED_CLUSTER_COUNT", unset = "9"),
  "Environment variable CLUSTER_ANNOTATION_EXPECTED_CLUSTER_COUNT"
)
annotation_primary_col <- paste0(analysis_cluster_col, "_annotation_primary")
annotation_multi_col <- paste0(analysis_cluster_col, "_annotation_multi")

deg_min_pct <- 0.10
deg_logfc_threshold <- 0
deg_padj_cutoff <- 0.05
deg_lfc_cutoff_for_ora <- 0.25
deg_abs_delta_pct_cutoff_for_ora <- 0.05
deg_top_n_for_ora <- 100L
ora_annotation_fdr_cutoff <- 0.05
annotation_weight_expression <- 1 / 3
annotation_weight_detection <- 1 / 3
annotation_weight_overlap <- 1 / 3
min_hallmark_size <- 15
max_hallmark_size <- 500
ora_min_overlap <- 3
top_terms_per_cluster <- 3L
top_plot_terms <- 15L
future_globals_maxsize_gb <- 30
future_plan_strategy <- Sys.getenv("SEURAT_FUTURE_PLAN", unset = "auto")
future_workers_env <- trimws(Sys.getenv("SEURAT_FUTURE_WORKERS", unset = ""))

parse_optional_positive_integer <- function(x, source_name) {
  if (is.null(x) || length(x) == 0) return(NULL)
  x <- trimws(as.character(x)[1])
  if (!nzchar(x)) return(NULL)
  value <- suppressWarnings(as.integer(x))
  if (is.na(value) || value < 1L) {
    stop(source_name, " must be a positive integer when set.", call. = FALSE)
  }
  value
}

future_workers_config <- parse_optional_positive_integer(config$Seurat_future_workers, "Config field 'Seurat_future_workers'")
future_workers_override <- parse_optional_positive_integer(future_workers_env, "Environment variable SEURAT_FUTURE_WORKERS")
future_workers <- if (!is.null(future_workers_override)) future_workers_override else future_workers_config

message("Preparing output directories.")
.ensure_dir(output_root)
out_summary <- .ensure_dir(file.path(output_root, "00_summary"))
out_deg <- .ensure_dir(file.path(output_root, paste0("01_", analysis_cluster_stub, "_vs_rest_DEG")))
out_ora <- .ensure_dir(file.path(output_root, "02_hallmark_ORA"))
out_plots <- .ensure_dir(file.path(output_root, "03_plots"))
out_objects <- .ensure_dir(file.path(output_root, "04_objects"))

if (!file.exists(input_rds)) {
  stop("Input Seurat object does not exist: ", input_rds, call. = FALSE)
}

message("Reading Seurat object: ", input_rds)
obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("Input file is not a Seurat object.", call. = FALSE)
}
assert_required_column(obj@meta.data, analysis_cluster_col)
if (!("RNA" %in% names(obj@assays))) {
  stop("RNA assay is missing in Seurat object.", call. = FALSE)
}

future_plan_info <- configure_future_for_seurat(
  max_size_gb = future_globals_maxsize_gb,
  strategy = future_plan_strategy,
  workers = future_workers
)
message(
  "Configured Seurat future plan: ",
  future_plan_info$strategy,
  " (workers=", future_plan_info$workers, ")."
)

obj <- prepare_obj_for_markers(obj, assay = "RNA")

cluster_values <- as.character(obj@meta.data[[analysis_cluster_col]])
if (any(is.na(cluster_values) | cluster_values == "")) {
  stop("Missing values were found in ", analysis_cluster_col, ".", call. = FALSE)
}
cluster_levels <- if (is.factor(obj@meta.data[[analysis_cluster_col]])) {
  as.character(levels(obj@meta.data[[analysis_cluster_col]]))
} else {
  sort_maybe_numeric(cluster_values)
}
cluster_levels <- cluster_levels[cluster_levels %in% unique(cluster_values)]
obj@meta.data[[analysis_cluster_col]] <- factor(cluster_values, levels = cluster_levels)
Idents(obj) <- obj@meta.data[[analysis_cluster_col]]

if (!is.na(expected_cluster_count) && length(cluster_levels) != expected_cluster_count) {
  stop(
    "Expected ",
    expected_cluster_count,
    " ",
    analysis_cluster_label,
    " groups, but found ",
    length(cluster_levels),
    ": ",
    paste(cluster_levels, collapse = ", "),
    call. = FALSE
  )
}

deg_dirs <- setNames(
  lapply(cluster_levels, function(cluster_id) {
    .ensure_dir(file.path(out_deg, paste0("cluster_", cluster_id)))
  }),
  cluster_levels
)
ora_dirs <- setNames(
  lapply(cluster_levels, function(cluster_id) {
    .ensure_dir(file.path(out_ora, paste0("cluster_", cluster_id)))
  }),
  cluster_levels
)

cluster_count_df <- data.frame(
  cluster = names(table(obj@meta.data[[analysis_cluster_col]])),
  n_cells = as.integer(table(obj@meta.data[[analysis_cluster_col]])),
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(match(cluster, cluster_levels))
write_table_csv(cluster_count_df, file.path(out_summary, paste0(analysis_cluster_stub, "_cell_counts.csv")))

message("Preparing Hallmark gene sets.")
hallmark_sets <- get_hallmark_sets(species = "Homo sapiens")
universe_genes <- unique(stats::na.omit(clean_gene_symbols(rownames(obj[["RNA"]]))))
write_table_csv(
  data.frame(gene_symbol = universe_genes, stringsAsFactors = FALSE),
  file.path(out_summary, "hallmark_ora_universe_genes.csv")
)

deg_summary_list <- list()
ora_all_list <- list()
annotation_summary_list <- list()

message("Running ", analysis_cluster_label, " DEG and Hallmark ORA.")
for (cluster_id in cluster_levels) {
  message("[DEG] ", analysis_cluster_label, " ", cluster_id, " vs rest")

  deg_dir_cluster <- deg_dirs[[cluster_id]]
  ora_dir_cluster <- ora_dirs[[cluster_id]]
  deg_markers_file <- file.path(deg_dir_cluster, paste0("cluster_", cluster_id, "_vs_rest_markers.csv"))

  if (file.exists(deg_markers_file)) {
    message("  Reusing existing DEG file: ", deg_markers_file)
    de <- tryCatch(
      readr::read_csv(deg_markers_file, show_col_types = FALSE),
      error = function(e) {
        stop(
          "Failed to read existing DEG file for ",
          analysis_cluster_label,
          " ",
          cluster_id,
          ": ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    de <- as.data.frame(de, stringsAsFactors = FALSE)
  } else {
    de <- tryCatch(
      FindMarkers(
        object = obj,
        ident.1 = cluster_id,
        assay = "RNA",
        slot = "data",
        min.pct = deg_min_pct,
        logfc.threshold = deg_logfc_threshold,
        verbose = FALSE
      ),
      error = function(e) {
        msg <- paste0("FindMarkers failed for ", analysis_cluster_label, " ", cluster_id, ": ", conditionMessage(e))
        warning(msg)
        writeLines(msg, con = file.path(deg_dir_cluster, "FindMarkers_error.txt"))
        NULL
      }
    )
  }

  n_cells_cluster <- cluster_count_df$n_cells[match(cluster_id, cluster_count_df$cluster)]
  if (is.null(de) || nrow(de) == 0) {
    deg_summary_list[[cluster_id]] <- data.frame(
      cluster = as.character(cluster_id),
      n_cells = n_cells_cluster,
      tested_genes = 0L,
      sig_genes = 0L,
      deg_for_ora = 0L,
      up_genes_for_ora = 0L,
      down_genes_for_ora = 0L,
      stringsAsFactors = FALSE
    )
    annotation_summary_list[[cluster_id]] <- summarize_cluster_annotation(
      cluster_id = cluster_id,
      n_cells = n_cells_cluster,
      deg_df = data.frame(),
      ora_df = data.frame(),
      top_n = top_terms_per_cluster,
      fdr_cutoff = ora_annotation_fdr_cutoff
    )
    next
  }

  if (!("gene" %in% colnames(de))) {
    de$gene <- rownames(de)
  }
  if (all(is.na(de$gene) | de$gene == "")) {
    stop(
      "DEG result for ",
      analysis_cluster_label,
      " ",
      cluster_id,
      " does not contain a usable 'gene' column.",
      call. = FALSE
    )
  }

  if (!("gene_symbol" %in% colnames(de))) {
    de$gene_symbol <- clean_gene_symbols(de$gene)
  } else {
    de$gene_symbol <- clean_gene_symbols(de$gene_symbol)
    missing_gene_symbol <- is.na(de$gene_symbol) | de$gene_symbol == ""
    de$gene_symbol[missing_gene_symbol] <- clean_gene_symbols(de$gene[missing_gene_symbol])
  }

  lfc_col <- resolve_lfc_col(de)
  if (!("p_val_adj" %in% colnames(de))) {
    stop("DEG result for ", analysis_cluster_label, " ", cluster_id, " is missing required column 'p_val_adj'.", call. = FALSE)
  }
  required_ora_cols <- c("pct.1", "pct.2")
  missing_ora_cols <- setdiff(required_ora_cols, colnames(de))
  if (length(missing_ora_cols) > 0) {
    stop(
      "DEG result for ",
      analysis_cluster_label,
      " ",
      cluster_id,
      " is missing required columns for top100-up ORA selection: ",
      paste(missing_ora_cols, collapse = ", "),
      call. = FALSE
    )
  }

  de <- de[order(de$p_val_adj, -abs(de[[lfc_col]])), , drop = FALSE]
  rownames(de) <- NULL
  if (!file.exists(deg_markers_file)) {
    write_table_csv(de, deg_markers_file)
  }

  ora_deg <- select_deg_top_up_for_ora(
    df = de,
    lfc_col = lfc_col,
    padj_max = deg_padj_cutoff,
    abs_logfc_min = deg_lfc_cutoff_for_ora,
    abs_delta_pct_min = deg_abs_delta_pct_cutoff_for_ora,
    top_n = deg_top_n_for_ora
  )

  write_table_csv(
    ora_deg,
    file.path(deg_dir_cluster, paste0("cluster_", cluster_id, "_vs_rest_top100_up_for_Hallmark_ORA.csv"))
  )

  ora_res <- run_ora_hypergeom(
    query_genes = ora_deg$gene_symbol,
    universe_genes = universe_genes,
    pathways = hallmark_sets,
    min_size = min_hallmark_size,
    max_size = max_hallmark_size,
    min_overlap = ora_min_overlap
  )
  ora_res <- score_ora_for_annotation(
    ora_df = ora_res,
    deg_df = ora_deg,
    weight_expression = annotation_weight_expression,
    weight_detection = annotation_weight_detection,
    weight_overlap = annotation_weight_overlap,
    fdr_cutoff = ora_annotation_fdr_cutoff
  )

  if (nrow(ora_res) > 0) {
    ora_res$cluster <- as.character(cluster_id)
    ora_res <- ora_res %>%
      dplyr::select(
        cluster,
        annotation_rank,
        pathway,
        hallmark_label,
        set_size,
        query_size,
        overlap,
        overlap_query_fraction,
        overlap_set_fraction,
        overlap_score_raw,
        marker_mean_abs_logfc,
        marker_mean_pct1,
        marker_mean_abs_delta_pct,
        n_overlap_up,
        n_overlap_down,
        expression_score_raw,
        detection_score_raw,
        expression_score_scaled,
        detection_score_scaled,
        overlap_score_scaled,
        annotation_score,
        gene_ratio,
        bg_ratio,
        odds_ratio,
        p_value,
        p_adj,
        overlap_genes
      )
    ora_all_list[[cluster_id]] <- ora_res
  }

  write_table_csv(ora_res, file.path(ora_dir_cluster, paste0("cluster_", cluster_id, "_Hallmark_ORA.csv")))
  plot_ora_top(
    ora_df = ora_res,
    out_pdf = file.path(ora_dir_cluster, paste0("cluster_", cluster_id, "_Hallmark_ORA_top.pdf")),
    out_png = file.path(ora_dir_cluster, paste0("cluster_", cluster_id, "_Hallmark_ORA_top.png")),
    title = paste0(analysis_cluster_label, " ", cluster_id, " vs rest | Hallmark ORA"),
    top_n = top_plot_terms
  )

  deg_summary_list[[cluster_id]] <- data.frame(
    cluster = as.character(cluster_id),
    n_cells = n_cells_cluster,
    tested_genes = nrow(de),
    sig_genes = sum(!is.na(de$p_val_adj) & de$p_val_adj < deg_padj_cutoff),
    deg_for_ora = nrow(ora_deg),
    up_genes_for_ora = sum(ora_deg$direction == "up", na.rm = TRUE),
    down_genes_for_ora = 0L,
    stringsAsFactors = FALSE
  )

  annotation_summary_list[[cluster_id]] <- summarize_cluster_annotation(
    cluster_id = cluster_id,
    n_cells = n_cells_cluster,
    deg_df = ora_deg,
    ora_df = ora_res,
    top_n = top_terms_per_cluster,
    fdr_cutoff = ora_annotation_fdr_cutoff
  )
}

deg_summary_df <- dplyr::bind_rows(deg_summary_list)
write_table_csv(deg_summary_df, file.path(out_summary, paste0(analysis_cluster_stub, "_deg_summary.csv")))

ora_all_df <- if (length(ora_all_list) > 0) dplyr::bind_rows(ora_all_list) else data.frame()
write_table_csv(ora_all_df, file.path(out_summary, paste0(analysis_cluster_stub, "_Hallmark_ORA_all.csv")))

annotation_summary_df <- dplyr::bind_rows(annotation_summary_list) %>%
  dplyr::arrange(match(cluster, cluster_levels))
write_table_csv(annotation_summary_df, file.path(out_summary, paste0(analysis_cluster_stub, "_annotation_summary.csv")))

message("Writing cluster-level annotation metadata back to object.")
primary_map <- stats::setNames(annotation_summary_df$annotation_primary, annotation_summary_df$cluster)
multi_map <- stats::setNames(annotation_summary_df$annotation_multi, annotation_summary_df$cluster)
obj@meta.data[[annotation_primary_col]] <- unname(primary_map[as.character(obj@meta.data[[analysis_cluster_col]])])
obj@meta.data[[annotation_multi_col]] <- unname(multi_map[as.character(obj@meta.data[[analysis_cluster_col]])])

message("Writing summary plots.")
stackfig_outputs <- write_stackfig_outputs(
  obj = obj,
  output_dir = out_plots,
  group_col_candidates = c("sample_type", "sample.type", "SampleType", "sampleType"),
  fallback_sample_col_candidates = c("sample", "Sample"),
  cluster_col_candidates = c(analysis_cluster_col),
  group_label = "sample_type",
  cluster_label = analysis_cluster_col,
  group_title = paste0(analysis_cluster_label, " proportion within each sample_type"),
  cluster_title = paste0("sample_type proportion within each ", analysis_cluster_label),
  group_by_cluster_stub = paste0("stack_sample_type_by_", analysis_cluster_stub),
  cluster_by_group_stub = paste0("stack_", analysis_cluster_stub, "_by_sample_type"),
  group_fill_colors = c(
    "2N-cellline" = "#BFD7EA",
    "4N-cellline" = "#2F6690",
    "2N-tumor" = "#F7C59F",
    "4N-tumor" = "#C66B3D"
  ),
  cluster_order_by_group_sum = c("2N-tumor", "4N-tumor"),
  cluster_order_tiebreak_groups = "4N-tumor",
  width = 10,
  height = 6
)

write_stackfig_outputs(
  obj = obj,
  output_dir = out_plots,
  group_col_candidates = c("sample", "Sample", "orig.ident", "sample_id", "SampleID"),
  fallback_sample_col_candidates = character(0),
  cluster_col_candidates = c(analysis_cluster_col),
  group_label = "sample",
  cluster_label = analysis_cluster_col,
  group_title = paste0(analysis_cluster_label, " proportion within each sample"),
  cluster_title = paste0("sample proportion within each ", analysis_cluster_label),
  group_by_cluster_stub = paste0("stack_sample_by_", analysis_cluster_stub),
  cluster_by_group_stub = paste0("stack_", analysis_cluster_stub, "_by_sample"),
  width = 12,
  height = 6
)

heatmap_cluster_levels <- levels(stackfig_outputs$cluster_group[[analysis_cluster_col]])
if (is.null(heatmap_cluster_levels) || length(heatmap_cluster_levels) == 0) {
  heatmap_cluster_levels <- cluster_levels
}

if (nrow(ora_all_df) > 0) {
  use_annotation_heatmap <- "annotation_score" %in% colnames(ora_all_df) && any(is.finite(ora_all_df$annotation_score))
  heatmap_df <- if (use_annotation_heatmap) {
    ora_all_df %>% dplyr::mutate(score = annotation_score)
  } else {
    ora_all_df %>% dplyr::mutate(score = -log10(pmax(p_adj, 1e-300)))
  }

  top_pathways <- if (use_annotation_heatmap) {
    heatmap_df %>%
      dplyr::group_by(hallmark_label) %>%
      dplyr::summarise(best_score = max(score, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(best_score), hallmark_label) %>%
      dplyr::slice_head(n = 20L) %>%
      dplyr::pull(hallmark_label)
  } else {
    heatmap_df %>%
      dplyr::group_by(hallmark_label) %>%
      dplyr::summarise(best_fdr = min(p_adj, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(best_fdr, hallmark_label) %>%
      dplyr::slice_head(n = 20L) %>%
      dplyr::pull(hallmark_label)
  }

  heatmap_use <- heatmap_df %>%
    dplyr::filter(hallmark_label %in% top_pathways) %>%
    dplyr::group_by(hallmark_label, cluster) %>%
    dplyr::summarise(score = max(score, na.rm = TRUE), .groups = "drop")

  heatmap_mat <- matrix(
    0,
    nrow = length(top_pathways),
    ncol = length(heatmap_cluster_levels),
    dimnames = list(top_pathways, heatmap_cluster_levels)
  )

  if (nrow(heatmap_use) > 0) {
    row_idx <- match(heatmap_use$hallmark_label, rownames(heatmap_mat))
    col_idx <- match(heatmap_use$cluster, colnames(heatmap_mat))
    keep_idx <- !is.na(row_idx) & !is.na(col_idx)
    heatmap_mat[cbind(row_idx[keep_idx], col_idx[keep_idx])] <- heatmap_use$score[keep_idx]
  }

  if (nrow(heatmap_mat) > 0 && ncol(heatmap_mat) > 0) {
    pdf(file.path(out_plots, "Hallmark_annotation_heatmap.pdf"), width = 10, height = 8)
    pheatmap::pheatmap(
      heatmap_mat,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      border_color = NA,
      main = if (use_annotation_heatmap) {
        paste0(analysis_cluster_label, " Hallmark annotation (integrated score)")
      } else {
        paste0(analysis_cluster_label, " Hallmark annotation (-log10 FDR)")
      }
    )
    dev.off()
  }
}

if ("umap" %in% names(obj@reductions)) {
  p1 <- DimPlot(
    obj,
    reduction = "umap",
    group.by = analysis_cluster_col,
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) +
    labs(title = analysis_cluster_label)
  ggsave(file.path(out_plots, paste0("umap_", analysis_cluster_stub, ".pdf")), p1, width = 9, height = 7)
  ggsave(file.path(out_plots, paste0("umap_", analysis_cluster_stub, ".png")), p1, width = 9, height = 7, dpi = 300)

  p2 <- DimPlot(
    obj,
    reduction = "umap",
    group.by = annotation_primary_col,
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) +
    labs(title = paste0("Primary Hallmark annotation by ", analysis_cluster_label))
  ggsave(file.path(out_plots, paste0("umap_", analysis_cluster_stub, "_annotation_primary.pdf")), p2, width = 11, height = 7)
  ggsave(file.path(out_plots, paste0("umap_", analysis_cluster_stub, "_annotation_primary.png")), p2, width = 11, height = 7, dpi = 300)
}

summary_lines <- c(
  paste0(annotation_script_name, " cluster annotation completed."),
  paste0("Input object: ", input_rds),
  paste0("Output root: ", output_root),
  paste0("Analysis cluster column: ", analysis_cluster_col),
  paste0(analysis_cluster_label, " count: ", length(cluster_levels)),
  paste0(analysis_cluster_label, " labels: ", paste(cluster_levels, collapse = ", ")),
  paste0("DEG test: Seurat default FindMarkers test."),
  paste0(
    "Hallmark ORA DEG rule: top ",
    deg_top_n_for_ora,
    " up genes after p_adj < ",
    deg_padj_cutoff,
    ", abs(logFC) >= ",
    deg_lfc_cutoff_for_ora,
    ", abs(pct.1 - pct.2) >= ",
    deg_abs_delta_pct_cutoff_for_ora,
    "."
  ),
  paste0(
    "Annotation rule: integrated score = ",
    annotation_weight_expression,
    " * marker_strength + ",
    annotation_weight_detection,
    " * marker_detection_fraction + ",
    annotation_weight_overlap,
    " * overlap_degree; only Hallmark ORA terms with FDR < ",
    ora_annotation_fdr_cutoff,
    " enter scoring after within-cluster min-max scaling."
  ),
  "Stack figure outputs:",
  paste0("  03_plots/stack_sample_type_by_", analysis_cluster_stub, ".pdf(.png)"),
  paste0("  03_plots/stack_", analysis_cluster_stub, "_by_sample_type.pdf(.png)"),
  paste0("  03_plots/stack_sample_by_", analysis_cluster_stub, ".pdf(.png)"),
  paste0("  03_plots/stack_", analysis_cluster_stub, "_by_sample.pdf(.png)"),
  "",
  "UMAP rendering:",
  "  All UMAP outputs are written with raster = FALSE (no point downsampling)."
)
writeLines(summary_lines, con = file.path(out_summary, "run_summary.txt"))

message("Saving annotated Seurat object.")
saveRDS(obj, file.path(out_objects, paste0("integrated_sct_cca_seurat_", analysis_cluster_stub, "_annotation.rds")))

message(annotation_script_name, " cluster annotation finished.")
message("Output root: ", output_root)
