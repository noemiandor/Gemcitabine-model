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
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
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

input_rds <- file.path(results_root, "02b_cluster_refine", "objects", "integrated_sct_cca_seurat_cluster_refine.rds")
output_root <- file.path(results_root, "02e_cluster_annotation")
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
deg_padj_cutoff <- 0.05
deg_lfc_cutoff_for_ora <- 0.25
min_hallmark_size <- 15
max_hallmark_size <- 500
ora_min_overlap <- 3
top_terms_per_cluster <- 3L
top_plot_terms <- 15L
future_globals_maxsize_gb <- 30

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

beautify_hallmark_name <- function(x) {
  x <- as.character(x)
  x <- sub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

get_hallmark_sets <- function(species = "Homo sapiens") {
  df_h <- tryCatch(
    msigdbr::msigdbr(species = species, collection = "H"),
    error = function(e) msigdbr::msigdbr(species = species, category = "H")
  )
  sets <- split(df_h$gene_symbol, df_h$gs_name)
  lapply(sets, unique)
}

run_ora_hypergeom <- function(query_genes, universe_genes, pathways, min_size = 15, max_size = 500, min_overlap = 3) {
  query <- unique(na.omit(as.character(query_genes)))
  universe <- unique(na.omit(as.character(universe_genes)))
  if (length(query) == 0 || length(universe) == 0) return(data.frame())
  query <- intersect(query, universe)
  if (length(query) == 0) return(data.frame())

  out <- lapply(names(pathways), function(pw_name) {
    pw_genes <- unique(intersect(as.character(pathways[[pw_name]]), universe))
    M <- length(pw_genes)
    if (M < min_size || M > max_size) return(NULL)
    overlap <- intersect(query, pw_genes)
    k <- length(overlap)
    if (k < min_overlap) return(NULL)

    U <- length(universe)
    N <- length(query)
    pval <- stats::phyper(q = k - 1, m = M, n = U - M, k = N, lower.tail = FALSE)
    odds <- suppressWarnings((k / max(1, N - k)) / (M / max(1, U - M)))

    data.frame(
      pathway = pw_name,
      hallmark_label = beautify_hallmark_name(pw_name),
      set_size = M,
      query_size = N,
      overlap = k,
      gene_ratio = sprintf("%d/%d", k, N),
      bg_ratio = sprintf("%d/%d", M, U),
      odds_ratio = odds,
      p_value = pval,
      p_adj = NA_real_,
      overlap_genes = paste(sort(overlap), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  if (is.null(out) || nrow(out) == 0) return(data.frame())
  out$p_adj <- p.adjust(out$p_value, method = "BH")
  out <- out[order(out$p_adj, out$p_value, -out$overlap), , drop = FALSE]
  rownames(out) <- NULL
  out
}

plot_ora_top <- function(ora_df, out_pdf, out_png, title, top_n = 15) {
  if (is.null(ora_df) || nrow(ora_df) == 0) {
    pdf(out_pdf, width = 10, height = 5)
    plot.new()
    text(0.5, 0.5, paste0(title, "\nNo Hallmark ORA results"))
    dev.off()
    png(out_png, width = 3000, height = 1500, res = 300)
    plot.new()
    text(0.5, 0.5, paste0(title, "\nNo Hallmark ORA results"))
    dev.off()
    return(invisible(NULL))
  }

  df <- ora_df
  df <- df[is.finite(df$p_adj) & !is.na(df$p_adj), , drop = FALSE]
  if (nrow(df) == 0) {
    pdf(out_pdf, width = 10, height = 5)
    plot.new()
    text(0.5, 0.5, paste0(title, "\nNo valid Hallmark ORA rows"))
    dev.off()
    png(out_png, width = 3000, height = 1500, res = 300)
    plot.new()
    text(0.5, 0.5, paste0(title, "\nNo valid Hallmark ORA rows"))
    dev.off()
    return(invisible(NULL))
  }

  df_sig <- df[df$p_adj < 0.05, , drop = FALSE]
  if (nrow(df_sig) > 0) {
    df <- df_sig
  }
  df <- head(df, top_n)
  df$score <- -log10(pmax(df$p_adj, 1e-300))
  df$hallmark_label <- factor(df$hallmark_label, levels = rev(df$hallmark_label))

  p <- ggplot(df, aes(x = hallmark_label, y = score)) +
    geom_col(fill = "#2c7fb8", width = 0.8) +
    coord_flip() +
    labs(title = title, x = NULL, y = "-log10(FDR)") +
    theme_classic(base_size = 11)

  ggsave(out_pdf, p, width = 10, height = max(5, 0.28 * nrow(df)))
  ggsave(out_png, p, width = 10, height = max(5, 0.28 * nrow(df)), dpi = 300)
  invisible(df)
}

build_merged_clusters <- function(cluster_vec, rules) {
  out <- as.character(cluster_vec)
  idx <- out %in% names(rules)
  out[idx] <- unname(rules[out[idx]])
  out
}

make_merge_mapping_table <- function(original_clusters, rules) {
  original_clusters <- sort_maybe_numeric(original_clusters)
  data.frame(
    original_cluster = original_clusters,
    analysis_cluster = ifelse(original_clusters %in% names(rules), unname(rules[original_clusters]), original_clusters),
    retained_cluster_label = ifelse(original_clusters %in% names(rules), "", original_clusters),
    merge_status = ifelse(original_clusters %in% names(rules), "absorbed", "retained"),
    stringsAsFactors = FALSE
  )
}

summarize_cluster_annotation <- function(cluster_id, n_cells, deg_df, ora_df, top_n = 3L) {
  up_genes_n <- if (is.null(deg_df) || nrow(deg_df) == 0) 0L else nrow(deg_df)

  if (is.null(ora_df) || nrow(ora_df) == 0) {
    return(data.frame(
      cluster = as.character(cluster_id),
      n_cells = n_cells,
      n_up_deg_for_ora = up_genes_n,
      n_significant_hallmarks = 0L,
      annotation_primary = NA_character_,
      annotation_secondary = NA_character_,
      annotation_tertiary = NA_character_,
      annotation_multi = NA_character_,
      top_hallmark_1_fdr = NA_real_,
      top_hallmark_2_fdr = NA_real_,
      top_hallmark_3_fdr = NA_real_,
      note = "no Hallmark ORA result",
      stringsAsFactors = FALSE
    ))
  }

  ora_sig <- ora_df %>%
    dplyr::filter(!is.na(p_adj), p_adj < 0.05) %>%
    dplyr::arrange(p_adj, dplyr::desc(overlap), pathway)

  note <- if (nrow(ora_sig) > 0) {
    "annotation based on significant Hallmark ORA"
  } else {
    "no Hallmark passes FDR < 0.05; top nominal terms are listed"
  }

  ora_use <- if (nrow(ora_sig) > 0) ora_sig else ora_df %>% dplyr::arrange(p_adj, dplyr::desc(overlap), pathway)
  ora_use <- head(ora_use, top_n)

  labels <- ora_use$hallmark_label
  padj_values <- ora_use$p_adj
  labels <- c(labels, rep(NA_character_, max(0, top_n - length(labels))))
  padj_values <- c(padj_values, rep(NA_real_, max(0, top_n - length(padj_values))))

  data.frame(
    cluster = as.character(cluster_id),
    n_cells = n_cells,
    n_up_deg_for_ora = up_genes_n,
    n_significant_hallmarks = nrow(ora_sig),
    annotation_primary = labels[1],
    annotation_secondary = labels[2],
    annotation_tertiary = labels[3],
    annotation_multi = paste(labels[!is.na(labels)], collapse = "; "),
    top_hallmark_1_fdr = padj_values[1],
    top_hallmark_2_fdr = padj_values[2],
    top_hallmark_3_fdr = padj_values[3],
    note = note,
    stringsAsFactors = FALSE
  )
}

.ensure_dir(output_root)
out_summary <- .ensure_dir(file.path(output_root, "00_summary"))
out_deg <- .ensure_dir(file.path(output_root, "01_cluster_final_vs_rest_DEG"))
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

  if (!(base_cluster_col %in% colnames(obj@meta.data))) {
    stop("Metadata column '", base_cluster_col, "' is missing.", call. = FALSE)
  }

  if (!("RNA" %in% names(obj@assays))) {
    stop("RNA assay is missing in Seurat object.", call. = FALSE)
  }

  configure_future_for_seurat(max_size_gb = future_globals_maxsize_gb)
  DefaultAssay(obj) <- "RNA"
  obj <- maybe_join_layers(obj, assay = "RNA")

  rna_data <- get_assay_data_slot(obj, assay = "RNA", slot_name = "data")
  if (is.null(rna_data) || nrow(rna_data) == 0 || ncol(rna_data) == 0) {
    message("RNA data slot is empty. Running NormalizeData on RNA assay.")
    obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
  }

  original_clusters <- as.character(obj@meta.data[[base_cluster_col]])
  if (any(is.na(original_clusters) | original_clusters == "")) {
    stop("Missing values were found in ", base_cluster_col, ".", call. = FALSE)
  }

  merged_clusters <- build_merged_clusters(original_clusters, merge_rules)
  merged_levels <- sort_maybe_numeric(merged_clusters)
  obj@meta.data[[analysis_cluster_col]] <- factor(merged_clusters, levels = merged_levels)
  Idents(obj) <- obj@meta.data[[analysis_cluster_col]]

  deg_dirs <- setNames(
    lapply(merged_levels, function(cluster_id) {
      .ensure_dir(file.path(out_deg, paste0("cluster_", cluster_id)))
    }),
    merged_levels
  )
  ora_dirs <- setNames(
    lapply(merged_levels, function(cluster_id) {
      .ensure_dir(file.path(out_ora, paste0("cluster_", cluster_id)))
    }),
    merged_levels
  )

  merge_mapping_df <- make_merge_mapping_table(original_clusters, merge_rules)
  write_csv(merge_mapping_df, file.path(out_summary, "cluster_merge_mapping_for_annotation.csv"))

  merged_count_df <- data.frame(
    cluster = names(table(obj@meta.data[[analysis_cluster_col]])),
    n_cells = as.integer(table(obj@meta.data[[analysis_cluster_col]])),
    stringsAsFactors = FALSE
  )
  write_csv(merged_count_df, file.path(out_summary, "cluster_final_cell_counts.csv"))

  message("Preparing Hallmark gene sets.")
  hallmark_sets <- get_hallmark_sets(species = "Homo sapiens")
  universe_genes <- unique(na.omit(clean_gene_symbols(rownames(obj[["RNA"]]))))
  write_csv(
    data.frame(gene_symbol = universe_genes, stringsAsFactors = FALSE),
    file.path(out_summary, "hallmark_ora_universe_genes.csv")
  )

  deg_summary_list <- list()
  ora_all_list <- list()
  annotation_summary_list <- list()

  message("Running cluster_final DEG and Hallmark ORA.")
  for (cluster_id in merged_levels) {
    message("[DEG] cluster_final ", cluster_id, " vs rest")

    deg_dir_cluster <- deg_dirs[[cluster_id]]
    ora_dir_cluster <- ora_dirs[[cluster_id]]
    deg_markers_file <- file.path(deg_dir_cluster, paste0("cluster_", cluster_id, "_vs_rest_markers.csv"))

    if (file.exists(deg_markers_file)) {
      message("  Reusing existing DEG file: ", deg_markers_file)
      de <- tryCatch(
        readr::read_csv(deg_markers_file, show_col_types = FALSE),
        error = function(e) {
          stop(
            "Failed to read existing DEG file for cluster_final ",
            cluster_id, ": ", conditionMessage(e),
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
          test.use = deg_test_use,
          min.pct = deg_min_pct,
          logfc.threshold = deg_logfc_threshold,
          verbose = FALSE
        ),
        error = function(e) {
          msg <- paste0("FindMarkers failed for cluster_final ", cluster_id, ": ", conditionMessage(e))
          warning(msg)
          writeLines(msg, con = file.path(deg_dir_cluster, "FindMarkers_error.txt"))
          NULL
        }
      )
    }

    if (is.null(de) || nrow(de) == 0) {
      deg_summary_list[[cluster_id]] <- data.frame(
        cluster = as.character(cluster_id),
        n_cells = merged_count_df$n_cells[match(cluster_id, merged_count_df$cluster)],
        tested_genes = 0L,
        sig_genes = 0L,
        up_genes_for_ora = 0L,
        stringsAsFactors = FALSE
      )
      annotation_summary_list[[cluster_id]] <- summarize_cluster_annotation(
        cluster_id = cluster_id,
        n_cells = merged_count_df$n_cells[match(cluster_id, merged_count_df$cluster)],
        deg_df = data.frame(),
        ora_df = data.frame(),
        top_n = top_terms_per_cluster
      )
      next
    }

    if (!("gene" %in% colnames(de))) {
      de$gene <- rownames(de)
    }
    if (all(is.na(de$gene) | de$gene == "")) {
      stop(
        "DEG result for cluster_final ", cluster_id,
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
      stop(
        "DEG result for cluster_final ", cluster_id,
        " is missing required column 'p_val_adj'.",
        call. = FALSE
      )
    }

    de <- de[order(de$p_val_adj, -abs(de[[lfc_col]])), , drop = FALSE]
    rownames(de) <- NULL
    if (!file.exists(deg_markers_file)) {
      write_csv(de, deg_markers_file)
    }

    sig_up <- de %>%
      dplyr::filter(
        !is.na(p_val_adj),
        p_val_adj < deg_padj_cutoff,
        !is.na(gene_symbol),
        .data[[lfc_col]] >= deg_lfc_cutoff_for_ora
      ) %>%
      dplyr::arrange(p_val_adj, dplyr::desc(.data[[lfc_col]]), gene_symbol) %>%
      dplyr::distinct(gene_symbol, .keep_all = TRUE)

    write_csv(sig_up, file.path(deg_dir_cluster, paste0("cluster_", cluster_id, "_vs_rest_up_for_Hallmark_ORA.csv")))

    ora_res <- run_ora_hypergeom(
      query_genes = sig_up$gene_symbol,
      universe_genes = universe_genes,
      pathways = hallmark_sets,
      min_size = min_hallmark_size,
      max_size = max_hallmark_size,
      min_overlap = ora_min_overlap
    )

    if (nrow(ora_res) > 0) {
      ora_res$cluster <- as.character(cluster_id)
      ora_res <- ora_res %>%
        dplyr::select(cluster, pathway, hallmark_label, set_size, query_size, overlap, gene_ratio, bg_ratio, odds_ratio, p_value, p_adj, overlap_genes)
      ora_all_list[[cluster_id]] <- ora_res
    }

    write_csv(ora_res, file.path(ora_dir_cluster, paste0("cluster_", cluster_id, "_Hallmark_ORA.csv")))
    plot_ora_top(
      ora_df = ora_res,
      out_pdf = file.path(ora_dir_cluster, paste0("cluster_", cluster_id, "_Hallmark_ORA_top.pdf")),
      out_png = file.path(ora_dir_cluster, paste0("cluster_", cluster_id, "_Hallmark_ORA_top.png")),
      title = paste0("cluster_final ", cluster_id, " vs rest | Hallmark ORA"),
      top_n = top_plot_terms
    )

    deg_summary_list[[cluster_id]] <- data.frame(
      cluster = as.character(cluster_id),
      n_cells = merged_count_df$n_cells[match(cluster_id, merged_count_df$cluster)],
      tested_genes = nrow(de),
      sig_genes = sum(!is.na(de$p_val_adj) & de$p_val_adj < deg_padj_cutoff),
      up_genes_for_ora = nrow(sig_up),
      stringsAsFactors = FALSE
    )

    annotation_summary_list[[cluster_id]] <- summarize_cluster_annotation(
      cluster_id = cluster_id,
      n_cells = merged_count_df$n_cells[match(cluster_id, merged_count_df$cluster)],
      deg_df = sig_up,
      ora_df = ora_res,
      top_n = top_terms_per_cluster
    )
  }

  deg_summary_df <- bind_rows(deg_summary_list)
  write_csv(deg_summary_df, file.path(out_summary, "cluster_final_deg_summary.csv"))

  ora_all_df <- if (length(ora_all_list) > 0) bind_rows(ora_all_list) else data.frame()
  write_csv(ora_all_df, file.path(out_summary, "cluster_final_Hallmark_ORA_all.csv"))

  annotation_summary_df <- bind_rows(annotation_summary_list) %>%
    dplyr::arrange(match(cluster, merged_levels))
  write_csv(annotation_summary_df, file.path(out_summary, "cluster_final_annotation_summary.csv"))

  message("Writing cluster-level annotation metadata back to object.")
  primary_map <- setNames(annotation_summary_df$annotation_primary, annotation_summary_df$cluster)
  multi_map <- setNames(annotation_summary_df$annotation_multi, annotation_summary_df$cluster)

  obj@meta.data[["cluster_final_annotation_primary"]] <- unname(primary_map[as.character(obj@meta.data[[analysis_cluster_col]])])
  obj@meta.data[["cluster_final_annotation_multi"]] <- unname(multi_map[as.character(obj@meta.data[[analysis_cluster_col]])])

  saveRDS(obj, file.path(out_objects, "integrated_sct_cca_seurat_cluster_final_annotation.rds"))

  message("Writing summary plots.")
  if (nrow(ora_all_df) > 0) {
    heatmap_df <- ora_all_df %>%
      dplyr::mutate(score = -log10(pmax(p_adj, 1e-300)))

    top_pathways <- heatmap_df %>%
      dplyr::group_by(hallmark_label) %>%
      dplyr::summarise(best_fdr = min(p_adj, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(best_fdr, hallmark_label) %>%
      dplyr::slice_head(n = 20L) %>%
      dplyr::pull(hallmark_label)

    heatmap_use <- heatmap_df %>%
      dplyr::filter(hallmark_label %in% top_pathways) %>%
      dplyr::group_by(hallmark_label, cluster) %>%
      dplyr::summarise(score = max(score, na.rm = TRUE), .groups = "drop")

    heatmap_mat <- matrix(
      0,
      nrow = length(top_pathways),
      ncol = length(merged_levels),
      dimnames = list(top_pathways, merged_levels)
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
        main = "cluster_final Hallmark annotation (-log10 FDR)"
      )
      dev.off()
    }
  }

  if ("umap" %in% names(obj@reductions)) {
    p1 <- DimPlot(obj, reduction = "umap", group.by = analysis_cluster_col, label = TRUE, repel = TRUE, raster = TRUE) +
      labs(title = "cluster_final for annotation")
    ggsave(file.path(out_plots, "umap_cluster_final.pdf"), p1, width = 9, height = 7)
    ggsave(file.path(out_plots, "umap_cluster_final.png"), p1, width = 9, height = 7, dpi = 300)

    p2 <- DimPlot(obj, reduction = "umap", group.by = "cluster_final_annotation_primary", label = TRUE, repel = TRUE, raster = TRUE) +
      labs(title = "Primary Hallmark annotation by cluster_final")
    ggsave(file.path(out_plots, "umap_cluster_final_annotation_primary.pdf"), p2, width = 11, height = 7)
    ggsave(file.path(out_plots, "umap_cluster_final_annotation_primary.png"), p2, width = 11, height = 7, dpi = 300)
  }

  summary_lines <- c(
    "02d cluster annotation completed.",
    paste0("Input object: ", input_rds),
    paste0("Output root: ", output_root),
    paste0("Base cluster column: ", base_cluster_col),
    paste0("Analysis cluster column: ", analysis_cluster_col),
    paste0("cluster_final count: ", length(merged_levels)),
    paste0("cluster_final labels: ", paste(merged_levels, collapse = ", ")),
    paste0("Hard-coded merge rules: ", paste(paste(names(merge_rules), "->", unname(merge_rules)), collapse = "; ")),
    "",
    "Interpretation reminder:",
    "  Clusters 1 and 7 are absorbed into 0 for downstream annotation.",
    "  Clusters 11 and 12 are absorbed into 10 for downstream annotation."
  )
  writeLines(summary_lines, con = file.path(out_summary, "run_summary.txt"))

message("02d cluster annotation finished.")
message("Output root: ", output_root)
