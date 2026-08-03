#!/usr/bin/env Rscript

script_path <- NULL
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
file_match <- grep(file_arg, cmd_args, value = TRUE)
if (length(file_match) > 0) {
  script_path <- normalizePath(sub(file_arg, "", file_match[1]), mustWork = FALSE)
}
if (is.null(script_path) || !nzchar(script_path)) {
  frame_files <- vapply(
    sys.frames(),
    function(x) {
      if (!is.null(x$ofile)) normalizePath(x$ofile, mustWork = FALSE) else NA_character_
    },
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0) {
    script_path <- frame_files[length(frame_files)]
  }
}
script_dir <- if (!is.null(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
source(file.path(script_dir, "Utils.R"))
config <- load_in_vivo_config(file.path(script_dir, "in_vivo_config.yaml"))
results_root <- get_results_root(config)

required_packages <- c(
  "dplyr",
  "readr",
  "ggplot2",
  "msigdbr",
  "pheatmap",
  "fgsea"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(msigdbr)
  library(pheatmap)
  library(fgsea)
})

set.seed(12345)

input_rds <- file.path(results_root, "02b_cluster_refine", "objects", "integrated_sct_cca_seurat_cluster_refine.rds")
deg_root <- file.path(results_root, "02e_cluster_annotation", "01_cluster_final_vs_rest_DEG")
cluster_count_file <- file.path(results_root, "02e_cluster_annotation", "00_summary", "cluster_final_cell_counts.csv")
output_root <- file.path(results_root, "02f_prerank_GSEA")

base_cluster_col <- "seurat_cluster_refine"
analysis_cluster_col <- "cluster_final"
merge_rules <- c(
  "1" = "0",
  "7" = "0",
  "11" = "10",
  "12" = "10"
)

deg_test_use <- "wilcox"
deg_min_pct <- 0.10
deg_logfc_threshold <- 0
future_globals_maxsize_gb <- 30
future_plan_strategy <- Sys.getenv("SEURAT_FUTURE_PLAN", unset = "auto")
future_workers_env <- trimws(Sys.getenv("SEURAT_FUTURE_WORKERS", unset = ""))

gsea_species <- "Homo sapiens"
gsea_min_size <- 15L
gsea_max_size <- 500L
gsea_padj_cutoff <- 0.05
top_plot_terms_each_direction <- 10L
heatmap_top_terms <- 20L
kegg_subcollection_candidates <- c("CP:KEGG_MEDICUS", "CP:KEGG_LEGACY", "CP:KEGG")

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

clean_gene_symbols <- function(genes) {
  g <- as.character(genes)
  g <- trimws(g)
  g <- sub("^GRCh[0-9]+[-_]", "", g, ignore.case = TRUE)
  g <- sub("^GRCm39[-_]", "", g, ignore.case = TRUE)
  g <- sub("^hg38[-_]", "", g, ignore.case = TRUE)
  g <- sub("\\.[0-9]+$", "", g)
  g[g == ""] <- NA_character_
  g
}

normalize_gene_key <- function(x) {
  x <- clean_gene_symbols(x)
  x <- toupper(x)
  x[x == ""] <- NA_character_
  x
}

resolve_deg_gene_label <- function(df) {
  gene_symbol <- if ("gene_symbol" %in% colnames(df)) as.character(df$gene_symbol) else rep(NA_character_, nrow(df))
  gene <- if ("gene" %in% colnames(df)) as.character(df$gene) else rownames(df)
  out <- ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, gene)
  out[is.na(out) | out == ""] <- gene[is.na(out) | out == ""]
  out
}

beautify_hallmark_name <- function(x) {
  x <- as.character(x)
  x <- sub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

beautify_kegg_name <- function(x) {
  x <- as.character(x)
  x <- sub("^KEGG_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

call_msigdbr_collection <- function(species, collection, subcollection = NULL) {
  attempts <- if (is.null(subcollection)) {
    list(
      list(species = species, collection = collection),
      list(species = species, category = collection)
    )
  } else {
    list(
      list(species = species, collection = collection, subcollection = subcollection),
      list(species = species, category = collection, subcategory = subcollection)
    )
  }

  for (args in attempts) {
    res <- tryCatch(
      do.call(msigdbr::msigdbr, args),
      error = function(e) NULL
    )
    if (!is.null(res) && nrow(res) > 0) {
      return(as.data.frame(res, stringsAsFactors = FALSE))
    }
  }

  NULL
}

build_pathway_bundle <- function(df, collection_name, collection_source, label_fun = NULL) {
  if (is.null(df) || nrow(df) == 0) {
    stop("Gene-set table is empty for ", collection_name, ".", call. = FALSE)
  }

  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df$gene_key <- normalize_gene_key(df$gene_symbol)
  df <- df %>%
    dplyr::filter(!is.na(gene_key), !is.na(gs_name), gs_name != "") %>%
    dplyr::distinct(gs_name, gene_key, .keep_all = TRUE)

  if (nrow(df) == 0) {
    stop("No valid genes remain after normalization for ", collection_name, ".", call. = FALSE)
  }

  if (is.null(label_fun)) {
    label_vec <- df$gs_name
  } else {
    label_vec <- label_fun(df$gs_name)
  }

  term_meta <- data.frame(
    pathway = df$gs_name,
    pathway_label = label_vec,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::distinct(pathway, .keep_all = TRUE)

  pathways <- split(df$gene_key, df$gs_name)
  pathways <- lapply(pathways, unique)

  list(
    collection_name = collection_name,
    collection_source = collection_source,
    pathways = pathways,
    term_meta = term_meta
  )
}

get_hallmark_bundle <- function(species) {
  hallmark_df <- call_msigdbr_collection(species = species, collection = "H")
  build_pathway_bundle(
    df = hallmark_df,
    collection_name = "Hallmark",
    collection_source = "MSigDB H",
    label_fun = beautify_hallmark_name
  )
}

get_kegg_bundle <- function(species, subcollection_candidates) {
  for (subcollection in subcollection_candidates) {
    kegg_df <- call_msigdbr_collection(
      species = species,
      collection = "C2",
      subcollection = subcollection
    )
    if (!is.null(kegg_df) && nrow(kegg_df) > 0) {
      return(
        build_pathway_bundle(
          df = kegg_df,
          collection_name = "KEGG",
          collection_source = paste0("MSigDB C2:", subcollection),
          label_fun = beautify_kegg_name
        )
      )
    }
  }

  stop(
    "Unable to retrieve KEGG gene sets from msigdbr. Tried: ",
    paste(subcollection_candidates, collapse = ", "),
    call. = FALSE
  )
}

prepare_prerank_stats <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!("gene" %in% colnames(df))) {
    df$gene <- rownames(df)
  }

  if (nrow(df) == 0) {
    return(data.frame())
  }

  lfc_col <- resolve_lfc_col(df)
  p_col <- if ("p_val" %in% colnames(df)) "p_val" else if ("p_val_adj" %in% colnames(df)) "p_val_adj" else NULL
  if (is.null(p_col)) {
    stop("DEG table is missing both 'p_val' and 'p_val_adj'.", call. = FALSE)
  }

  df$gene_label <- resolve_deg_gene_label(df)
  df$gene_key <- normalize_gene_key(df$gene_label)
  df$gene_symbol <- clean_gene_symbols(df$gene_label)
  missing_gene_symbol <- is.na(df$gene_symbol) | df$gene_symbol == ""
  df$gene_symbol[missing_gene_symbol] <- df$gene_key[missing_gene_symbol]
  df$lfc_value <- suppressWarnings(as.numeric(df[[lfc_col]]))
  df$p_value_num <- suppressWarnings(as.numeric(df[[p_col]]))
  df$p_value_num[!is.finite(df$p_value_num) | is.na(df$p_value_num) | df$p_value_num <= 0] <- .Machine$double.xmin
  df$rank_stat <- sign(df$lfc_value) * (-log10(df$p_value_num))
  tie_break <- suppressWarnings(as.numeric(df$lfc_value))
  tie_break[!is.finite(tie_break)] <- 0
  df$rank_stat <- df$rank_stat + tie_break * 1e-9
  df$abs_rank_stat <- abs(df$rank_stat)
  df$abs_lfc <- abs(df$lfc_value)

  out <- df %>%
    dplyr::filter(
      !is.na(gene_key),
      !is.na(gene_symbol),
      is.finite(lfc_value),
      is.finite(rank_stat)
    ) %>%
    dplyr::arrange(dplyr::desc(abs_rank_stat), dplyr::desc(abs_lfc), p_value_num) %>%
    dplyr::group_by(gene_key) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(rank_stat), dplyr::desc(abs_lfc), p_value_num) %>%
    dplyr::select(dplyr::any_of(c("gene_key", "gene_symbol", "rank_stat", "lfc_value", "p_value_num", "p_val_adj", "pct.1", "pct.2")))

  rownames(out) <- NULL
  out
}

run_prerank_fgsea <- function(rank_tbl, bundle, cluster_id, min_size, max_size) {
  if (is.null(rank_tbl) || nrow(rank_tbl) == 0) {
    return(data.frame())
  }

  stats <- rank_tbl$rank_stat
  names(stats) <- rank_tbl$gene_key
  stats <- stats[is.finite(stats) & !is.na(names(stats)) & names(stats) != ""]
  stats <- sort(stats, decreasing = TRUE)
  stats <- stats[!duplicated(names(stats))]

  if (length(stats) < min_size) {
    return(data.frame())
  }

  fg <- fgsea::fgseaMultilevel(
    pathways = bundle$pathways,
    stats = stats,
    minSize = min_size,
    maxSize = max_size,
    eps = 0
  )

  fg <- as.data.frame(fg, stringsAsFactors = FALSE)
  if (nrow(fg) == 0) {
    return(data.frame())
  }

  fg$leadingEdge_genes <- vapply(
    fg$leadingEdge,
    function(x) paste(x, collapse = ";"),
    character(1)
  )

  fg <- fg %>%
    dplyr::left_join(bundle$term_meta, by = c("pathway" = "pathway")) %>%
    dplyr::mutate(
      cluster = as.character(cluster_id),
      collection = bundle$collection_name,
      collection_source = bundle$collection_source,
      pathway_label = ifelse(is.na(pathway_label) | pathway_label == "", pathway, pathway_label),
      direction = ifelse(NES >= 0, "up", "down")
    ) %>%
    dplyr::select(dplyr::any_of(c(
      "cluster",
      "collection",
      "collection_source",
      "pathway",
      "pathway_label",
      "size",
      "ES",
      "NES",
      "pval",
      "padj",
      "log2err",
      "direction",
      "leadingEdge_genes"
    ))) %>%
    dplyr::arrange(padj, dplyr::desc(abs(NES)), pathway)

  rownames(fg) <- NULL
  fg
}

plot_top_fgsea <- function(gsea_df, file_stub, title, top_n_each_direction = 10L) {
  if (is.null(gsea_df) || nrow(gsea_df) == 0) {
    grDevices::pdf(paste0(file_stub, ".pdf"), width = 10, height = 5)
    plot.new()
    text(0.5, 0.5, paste0(title, "\nNo GSEA result"))
    dev.off()
    grDevices::png(paste0(file_stub, ".png"), width = 3000, height = 1500, res = 300)
    plot.new()
    text(0.5, 0.5, paste0(title, "\nNo GSEA result"))
    dev.off()
    grDevices::tiff(paste0(file_stub, ".tiff"), width = 3000, height = 1500, res = 300, compression = "lzw")
    plot.new()
    text(0.5, 0.5, paste0(title, "\nNo GSEA result"))
    dev.off()
    return(invisible(NULL))
  }

  sig_df <- gsea_df %>% dplyr::filter(!is.na(padj), padj < gsea_padj_cutoff)
  plot_df <- if (nrow(sig_df) > 0) sig_df else gsea_df

  pos_df <- plot_df %>%
    dplyr::filter(NES > 0) %>%
    dplyr::arrange(padj, dplyr::desc(NES), pathway) %>%
    dplyr::slice_head(n = top_n_each_direction)

  neg_df <- plot_df %>%
    dplyr::filter(NES < 0) %>%
    dplyr::arrange(padj, NES, pathway) %>%
    dplyr::slice_head(n = top_n_each_direction)

  plot_df <- dplyr::bind_rows(neg_df, pos_df)
  if (nrow(plot_df) == 0) {
    plot_df <- gsea_df %>%
      dplyr::arrange(padj, dplyr::desc(abs(NES)), pathway) %>%
      dplyr::slice_head(n = min(2L * top_n_each_direction, nrow(gsea_df)))
  }

  plot_df$pathway_label <- factor(plot_df$pathway_label, levels = plot_df$pathway_label[order(plot_df$NES)])

  p <- ggplot(plot_df, aes(x = pathway_label, y = NES, fill = direction)) +
    geom_col(width = 0.8) +
    coord_flip() +
    scale_fill_manual(values = c(down = "#2166ac", up = "#b2182b")) +
    labs(title = title, x = NULL, y = "NES", fill = "Direction") +
    theme_classic(base_size = 11)

  save_plot_pdf_png(
    p,
    file_stub,
    width = 10,
    height = max(5, 0.28 * nrow(plot_df)),
    dpi = 300,
    include_tiff = TRUE
  )

  invisible(plot_df)
}

plot_gsea_heatmap <- function(gsea_all_df, collection_name, file_stub, top_n_terms = 20L) {
  if (is.null(gsea_all_df) || nrow(gsea_all_df) == 0) {
    grDevices::pdf(paste0(file_stub, ".pdf"), width = 10, height = 5)
    plot.new()
    text(0.5, 0.5, paste0(collection_name, " GSEA\nNo results"))
    dev.off()
    grDevices::png(paste0(file_stub, ".png"), width = 3000, height = 1500, res = 300)
    plot.new()
    text(0.5, 0.5, paste0(collection_name, " GSEA\nNo results"))
    dev.off()
    grDevices::tiff(paste0(file_stub, ".tiff"), width = 3000, height = 1500, res = 300, compression = "lzw")
    plot.new()
    text(0.5, 0.5, paste0(collection_name, " GSEA\nNo results"))
    dev.off()
    return(invisible(NULL))
  }

  sig_df <- gsea_all_df %>%
    dplyr::filter(!is.na(padj), padj < gsea_padj_cutoff)

  if (nrow(sig_df) == 0) {
    grDevices::pdf(paste0(file_stub, ".pdf"), width = 10, height = 5)
    plot.new()
    text(0.5, 0.5, paste0(collection_name, " GSEA\nNo pathway passes FDR < ", gsea_padj_cutoff))
    dev.off()
    grDevices::png(paste0(file_stub, ".png"), width = 3000, height = 1500, res = 300)
    plot.new()
    text(0.5, 0.5, paste0(collection_name, " GSEA\nNo pathway passes FDR < ", gsea_padj_cutoff))
    dev.off()
    grDevices::tiff(paste0(file_stub, ".tiff"), width = 3000, height = 1500, res = 300, compression = "lzw")
    plot.new()
    text(0.5, 0.5, paste0(collection_name, " GSEA\nNo pathway passes FDR < ", gsea_padj_cutoff))
    dev.off()
    return(invisible(NULL))
  }

  top_terms <- sig_df %>%
    dplyr::group_by(pathway_label) %>%
    dplyr::summarise(best_padj = min(padj, na.rm = TRUE), best_abs_nes = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(best_padj, dplyr::desc(best_abs_nes), pathway_label) %>%
    dplyr::slice_head(n = top_n_terms) %>%
    dplyr::pull(pathway_label)

  clusters <- sort_maybe_numeric(sig_df$cluster)
  heatmap_mat <- matrix(
    0,
    nrow = length(top_terms),
    ncol = length(clusters),
    dimnames = list(top_terms, clusters)
  )

  sig_use <- sig_df %>%
    dplyr::filter(pathway_label %in% top_terms) %>%
    dplyr::group_by(pathway_label, cluster) %>%
    dplyr::slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  row_idx <- match(sig_use$pathway_label, rownames(heatmap_mat))
  col_idx <- match(sig_use$cluster, colnames(heatmap_mat))
  keep_idx <- !is.na(row_idx) & !is.na(col_idx)
  heatmap_mat[cbind(row_idx[keep_idx], col_idx[keep_idx])] <- sig_use$NES[keep_idx]

  grDevices::pdf(paste0(file_stub, ".pdf"), width = 10, height = max(6, 0.28 * nrow(heatmap_mat)))
  pheatmap::pheatmap(
    heatmap_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    main = paste0(collection_name, " pre-ranked GSEA NES")
  )
  dev.off()

  grDevices::png(paste0(file_stub, ".png"), width = 3000, height = max(1800, as.integer(84 * nrow(heatmap_mat))), res = 300)
  pheatmap::pheatmap(
    heatmap_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    main = paste0(collection_name, " pre-ranked GSEA NES")
  )
  dev.off()

  grDevices::tiff(
    paste0(file_stub, ".tiff"),
    width = 10,
    height = max(6, 0.28 * nrow(heatmap_mat)),
    units = "in",
    res = 300,
    compression = "lzw"
  )
  pheatmap::pheatmap(
    heatmap_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    main = paste0(collection_name, " pre-ranked GSEA NES")
  )
  dev.off()

  invisible(heatmap_mat)
}

load_cluster_final_object <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required to recompute missing DEG files.", call. = FALSE)
  }

  suppressPackageStartupMessages(library(Seurat))

  if (!file.exists(input_rds)) {
    stop("Input Seurat object does not exist: ", input_rds, call. = FALSE)
  }

  future_plan_info <- configure_future_for_seurat(
    max_size_gb = future_globals_maxsize_gb,
    strategy = future_plan_strategy,
    workers = future_workers
  )
  message(
    "Configured Seurat future plan for fallback FindMarkers: ",
    future_plan_info$strategy,
    " (workers=", future_plan_info$workers, ")."
  )

  obj <- readRDS(input_rds)
  if (!inherits(obj, "Seurat")) {
    stop("Input file is not a Seurat object.", call. = FALSE)
  }
  if (!(base_cluster_col %in% colnames(obj@meta.data))) {
    stop("Metadata column '", base_cluster_col, "' is missing.", call. = FALSE)
  }
  if (!("RNA" %in% names(obj@assays))) {
    stop("RNA assay is missing in Seurat object.", call. = FALSE)
  }

  Seurat::DefaultAssay(obj) <- "RNA"
  obj <- maybe_join_layers(obj, assay = "RNA")

  rna_data <- get_assay_data_slot(obj, assay = "RNA", slot_name = "data")
  if (is.null(rna_data) || nrow(rna_data) == 0 || ncol(rna_data) == 0) {
    message("RNA data slot is empty. Running NormalizeData on RNA assay.")
    obj <- Seurat::NormalizeData(obj, assay = "RNA", verbose = FALSE)
  }

  original_clusters <- as.character(obj@meta.data[[base_cluster_col]])
  if (any(is.na(original_clusters) | original_clusters == "")) {
    stop("Missing values were found in ", base_cluster_col, ".", call. = FALSE)
  }

  merged_clusters <- as.character(original_clusters)
  idx <- merged_clusters %in% names(merge_rules)
  merged_clusters[idx] <- unname(merge_rules[merged_clusters[idx]])
  merged_levels <- sort_maybe_numeric(merged_clusters)
  obj@meta.data[[analysis_cluster_col]] <- factor(merged_clusters, levels = merged_levels)
  Seurat::Idents(obj) <- obj@meta.data[[analysis_cluster_col]]

  list(obj = obj, merged_levels = merged_levels)
}

ensure_marker_file <- function(cluster_id, obj) {
  cluster_dir <- .ensure_dir(file.path(deg_root, paste0("cluster_", cluster_id)))
  marker_file <- file.path(cluster_dir, paste0("cluster_", cluster_id, "_vs_rest_markers.csv"))
  if (file.exists(marker_file)) {
    return(marker_file)
  }

  message("Recomputing missing DEG file for cluster_final ", cluster_id)
  de <- tryCatch(
    Seurat::FindMarkers(
      object = obj,
      ident.1 = cluster_id,
      assay = "RNA",
      slot = "data",
      test.use = deg_test_use,
      min.pct = deg_min_pct,
      logfc.threshold = deg_logfc_threshold,
      verbose = FALSE
    ),
    error = function(e) {
      stop("FindMarkers failed for cluster_final ", cluster_id, ": ", conditionMessage(e), call. = FALSE)
    }
  )

  de <- as.data.frame(de, stringsAsFactors = FALSE)
  if (!("gene" %in% colnames(de))) {
    de$gene <- rownames(de)
  }
  if (!("gene_symbol" %in% colnames(de))) {
    de$gene_symbol <- clean_gene_symbols(de$gene)
  } else {
    de$gene_symbol <- clean_gene_symbols(de$gene_symbol)
    missing_gene_symbol <- is.na(de$gene_symbol) | de$gene_symbol == ""
    de$gene_symbol[missing_gene_symbol] <- clean_gene_symbols(de$gene[missing_gene_symbol])
  }

  if (nrow(de) > 0) {
    lfc_col <- resolve_lfc_col(de)
    de <- de[order(de$p_val_adj, -abs(de[[lfc_col]])), , drop = FALSE]
  }
  rownames(de) <- NULL
  readr::write_csv(de, marker_file)
  marker_file
}

.ensure_dir(output_root)
out_summary <- .ensure_dir(file.path(output_root, "00_summary"))
out_rankings <- .ensure_dir(file.path(out_summary, "rankings"))
out_hallmark <- .ensure_dir(file.path(output_root, "01_hallmark_prerank_GSEA"))
out_kegg <- .ensure_dir(file.path(output_root, "02_kegg_prerank_GSEA"))
out_plots <- .ensure_dir(file.path(output_root, "03_plots"))

message("Loading Hallmark and KEGG gene sets for pre-ranked GSEA.")
hallmark_bundle <- get_hallmark_bundle(species = gsea_species)
kegg_bundle <- get_kegg_bundle(species = gsea_species, subcollection_candidates = kegg_subcollection_candidates)

if (!dir.exists(deg_root)) {
  .ensure_dir(deg_root)
}

existing_cluster_dirs <- list.dirs(deg_root, recursive = FALSE, full.names = FALSE)
existing_cluster_ids <- sub("^cluster_", "", existing_cluster_dirs)
existing_cluster_ids <- existing_cluster_ids[nzchar(existing_cluster_ids)]
existing_cluster_ids <- sort_maybe_numeric(existing_cluster_ids)
expected_cluster_ids <- existing_cluster_ids

if (file.exists(cluster_count_file)) {
  cluster_count_df <- tryCatch(
    readr::read_csv(cluster_count_file, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (!is.null(cluster_count_df) && "cluster" %in% colnames(cluster_count_df)) {
    expected_cluster_ids <- sort_maybe_numeric(as.character(cluster_count_df$cluster))
  }
}

obj_info <- NULL
cluster_ids <- expected_cluster_ids

if (length(cluster_ids) == 0) {
  message("No existing DEG directories found. Falling back to FindMarkers.")
  obj_info <- load_cluster_final_object()
  cluster_ids <- obj_info$merged_levels
}

marker_files <- setNames(
  file.path(deg_root, paste0("cluster_", cluster_ids), paste0("cluster_", cluster_ids, "_vs_rest_markers.csv")),
  cluster_ids
)
missing_marker_clusters <- names(marker_files)[!file.exists(marker_files)]

if (length(missing_marker_clusters) > 0) {
  if (is.null(obj_info)) {
    obj_info <- load_cluster_final_object()
    cluster_ids <- obj_info$merged_levels
    marker_files <- setNames(
      file.path(deg_root, paste0("cluster_", cluster_ids), paste0("cluster_", cluster_ids, "_vs_rest_markers.csv")),
      cluster_ids
    )
    missing_marker_clusters <- names(marker_files)[!file.exists(marker_files)]
  }

  for (cluster_id in missing_marker_clusters) {
    marker_files[[cluster_id]] <- ensure_marker_file(cluster_id, obj_info$obj)
  }
}

cluster_ids <- sort_maybe_numeric(names(marker_files))
marker_files <- marker_files[cluster_ids]

hallmark_all_list <- list()
kegg_all_list <- list()
cluster_summary_list <- list()

message("Running Hallmark and KEGG pre-ranked GSEA.")
for (cluster_id in cluster_ids) {
  message("[GSEA] cluster_final ", cluster_id)

  marker_file <- marker_files[[cluster_id]]
  if (!file.exists(marker_file)) {
    stop("Missing marker file for cluster_final ", cluster_id, ": ", marker_file, call. = FALSE)
  }

  de <- readr::read_csv(marker_file, show_col_types = FALSE)
  de <- as.data.frame(de, stringsAsFactors = FALSE)
  rank_tbl <- prepare_prerank_stats(de)

  if (nrow(rank_tbl) == 0) {
    warning("No valid ranked genes for cluster_final ", cluster_id)
    next
  }

  readr::write_csv(rank_tbl, file.path(out_rankings, paste0("cluster_", cluster_id, "_prerank_stats.csv")))

  hallmark_dir <- .ensure_dir(file.path(out_hallmark, paste0("cluster_", cluster_id)))
  kegg_dir <- .ensure_dir(file.path(out_kegg, paste0("cluster_", cluster_id)))

  hallmark_res <- run_prerank_fgsea(
    rank_tbl = rank_tbl,
    bundle = hallmark_bundle,
    cluster_id = cluster_id,
    min_size = gsea_min_size,
    max_size = gsea_max_size
  )
  kegg_res <- run_prerank_fgsea(
    rank_tbl = rank_tbl,
    bundle = kegg_bundle,
    cluster_id = cluster_id,
    min_size = gsea_min_size,
    max_size = gsea_max_size
  )

  readr::write_csv(
    hallmark_res,
    file.path(hallmark_dir, paste0("cluster_", cluster_id, "_Hallmark_prerank_GSEA.csv"))
  )
  readr::write_csv(
    kegg_res,
    file.path(kegg_dir, paste0("cluster_", cluster_id, "_KEGG_prerank_GSEA.csv"))
  )

  plot_top_fgsea(
    hallmark_res,
    file_stub = file.path(hallmark_dir, paste0("cluster_", cluster_id, "_Hallmark_prerank_GSEA_top")),
    title = paste0("cluster_final ", cluster_id, " | Hallmark pre-ranked GSEA"),
    top_n_each_direction = top_plot_terms_each_direction
  )
  plot_top_fgsea(
    kegg_res,
    file_stub = file.path(kegg_dir, paste0("cluster_", cluster_id, "_KEGG_prerank_GSEA_top")),
    title = paste0("cluster_final ", cluster_id, " | KEGG pre-ranked GSEA"),
    top_n_each_direction = top_plot_terms_each_direction
  )

  hallmark_all_list[[cluster_id]] <- hallmark_res
  kegg_all_list[[cluster_id]] <- kegg_res

  top_hallmark_up <- hallmark_res %>%
    dplyr::filter(!is.na(padj), padj < gsea_padj_cutoff, NES > 0) %>%
    dplyr::arrange(padj, dplyr::desc(NES), pathway) %>%
    dplyr::slice_head(n = 1)
  top_hallmark_down <- hallmark_res %>%
    dplyr::filter(!is.na(padj), padj < gsea_padj_cutoff, NES < 0) %>%
    dplyr::arrange(padj, NES, pathway) %>%
    dplyr::slice_head(n = 1)
  top_kegg_up <- kegg_res %>%
    dplyr::filter(!is.na(padj), padj < gsea_padj_cutoff, NES > 0) %>%
    dplyr::arrange(padj, dplyr::desc(NES), pathway) %>%
    dplyr::slice_head(n = 1)
  top_kegg_down <- kegg_res %>%
    dplyr::filter(!is.na(padj), padj < gsea_padj_cutoff, NES < 0) %>%
    dplyr::arrange(padj, NES, pathway) %>%
    dplyr::slice_head(n = 1)

  cluster_summary_list[[cluster_id]] <- data.frame(
    cluster = as.character(cluster_id),
    n_ranked_genes = nrow(rank_tbl),
    n_hallmark_sig = sum(!is.na(hallmark_res$padj) & hallmark_res$padj < gsea_padj_cutoff),
    top_hallmark_up = if (nrow(top_hallmark_up) > 0) top_hallmark_up$pathway_label[1] else NA_character_,
    top_hallmark_up_nes = if (nrow(top_hallmark_up) > 0) top_hallmark_up$NES[1] else NA_real_,
    top_hallmark_down = if (nrow(top_hallmark_down) > 0) top_hallmark_down$pathway_label[1] else NA_character_,
    top_hallmark_down_nes = if (nrow(top_hallmark_down) > 0) top_hallmark_down$NES[1] else NA_real_,
    n_kegg_sig = sum(!is.na(kegg_res$padj) & kegg_res$padj < gsea_padj_cutoff),
    top_kegg_up = if (nrow(top_kegg_up) > 0) top_kegg_up$pathway_label[1] else NA_character_,
    top_kegg_up_nes = if (nrow(top_kegg_up) > 0) top_kegg_up$NES[1] else NA_real_,
    top_kegg_down = if (nrow(top_kegg_down) > 0) top_kegg_down$pathway_label[1] else NA_character_,
    top_kegg_down_nes = if (nrow(top_kegg_down) > 0) top_kegg_down$NES[1] else NA_real_,
    stringsAsFactors = FALSE
  )
}

hallmark_all_df <- if (length(hallmark_all_list) > 0) dplyr::bind_rows(hallmark_all_list) else data.frame()
kegg_all_df <- if (length(kegg_all_list) > 0) dplyr::bind_rows(kegg_all_list) else data.frame()
cluster_summary_df <- if (length(cluster_summary_list) > 0) dplyr::bind_rows(cluster_summary_list) else data.frame()

readr::write_csv(hallmark_all_df, file.path(out_summary, "Hallmark_prerank_GSEA_all.csv"))
readr::write_csv(kegg_all_df, file.path(out_summary, "KEGG_prerank_GSEA_all.csv"))
readr::write_csv(cluster_summary_df, file.path(out_summary, "cluster_prerank_GSEA_summary.csv"))

plot_gsea_heatmap(
  gsea_all_df = hallmark_all_df,
  collection_name = "Hallmark",
  file_stub = file.path(out_plots, "Hallmark_prerank_GSEA_NES_heatmap"),
  top_n_terms = heatmap_top_terms
)
plot_gsea_heatmap(
  gsea_all_df = kegg_all_df,
  collection_name = "KEGG",
  file_stub = file.path(out_plots, "KEGG_prerank_GSEA_NES_heatmap"),
  top_n_terms = heatmap_top_terms
)

summary_lines <- c(
  "02f pre-ranked GSEA completed.",
  paste0("Input Seurat object (fallback for missing DEG only): ", input_rds),
  paste0("Input DEG root: ", deg_root),
  paste0("Output root: ", output_root),
  paste0("cluster_final labels: ", paste(cluster_ids, collapse = ", ")),
  paste0("Hallmark source: ", hallmark_bundle$collection_source),
  paste0("KEGG source: ", kegg_bundle$collection_source),
  paste0("Rank metric: sign(logFC) * -log10(p_val), with duplicate genes collapsed by max absolute rank statistic."),
  paste0("GSEA FDR cutoff for summary/heatmap: ", gsea_padj_cutoff),
  paste0("GSEA min/max gene-set size: ", gsea_min_size, "/", gsea_max_size),
  if (length(missing_marker_clusters) > 0) {
    paste0("Fallback FindMarkers rerun for missing clusters: ", paste(missing_marker_clusters, collapse = ", "))
  } else {
    "Fallback FindMarkers rerun for missing clusters: none"
  },
  "Main outputs:",
  "  00_summary/Hallmark_prerank_GSEA_all.csv",
  "  00_summary/KEGG_prerank_GSEA_all.csv",
  "  00_summary/cluster_prerank_GSEA_summary.csv",
  "  00_summary/rankings/cluster_*_prerank_stats.csv",
  "  01_hallmark_prerank_GSEA/cluster_*/cluster_*_Hallmark_prerank_GSEA.csv",
  "  02_kegg_prerank_GSEA/cluster_*/cluster_*_KEGG_prerank_GSEA.csv",
  "  03_plots/Hallmark_prerank_GSEA_NES_heatmap.pdf(.png/.tiff)",
  "  03_plots/KEGG_prerank_GSEA_NES_heatmap.pdf(.png/.tiff)"
)
writeLines(summary_lines, con = file.path(out_summary, "run_summary.txt"))

message("02f pre-ranked GSEA finished.")
message("Output root: ", output_root)
