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

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(Matrix)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  KMP_DUPLICATE_LIB_OK = "TRUE",
  KMP_INIT_AT_FORK = "FALSE"
)

set.seed(12345)

input_rds <- file.path(results_root, "02b_cluster_refine", "objects", "integrated_sct_cca_seurat_cluster_refine.rds")
output_root <- file.path(results_root, "02c_cluster_quality")
base_cluster_col <- "seurat_cluster_refine"
assay_use <- "RNA"
top_n_genes_fraction <- 50L
focus_clusters <- c("9", "9c")
flag_threshold_abs_robust_z <- 1.5

safe_mad <- function(x) {
  out <- stats::mad(x, center = stats::median(x, na.rm = TRUE), constant = 1, na.rm = TRUE)
  if (!is.finite(out) || out == 0) return(NA_real_)
  out
}

plot_simple_heatmap <- function(mat, file_pdf, title_text) {
  if (is.null(dim(mat)) || nrow(mat) == 0 || ncol(mat) == 0) {
    return(invisible(NULL))
  }

  cols <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)
  breaks <- seq(-3, 3, length.out = 102)
  mat_plot <- t(mat[nrow(mat):1, , drop = FALSE])
  z <- mat_plot
  z[!is.finite(z)] <- NA_real_

  pdf(file_pdf, width = 8, height = max(5, 0.35 * nrow(mat) + 2))
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  layout(matrix(c(1, 2), nrow = 1), widths = c(6, 1))
  par(mar = c(7, 9, 3, 2))
  graphics::image(
    x = seq_len(nrow(z)),
    y = seq_len(ncol(z)),
    z = z,
    col = cols,
    breaks = breaks,
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    main = title_text,
    useRaster = TRUE
  )
  axis(1, at = seq_len(nrow(z)), labels = rownames(z), las = 2, cex.axis = 0.8)
  axis(2, at = seq_len(ncol(z)), labels = colnames(z), las = 2, cex.axis = 0.85)

  na_idx <- which(is.na(z), arr.ind = TRUE)
  if (nrow(na_idx) > 0) {
    rect(
      xleft = na_idx[, 1] - 0.5,
      ybottom = na_idx[, 2] - 0.5,
      xright = na_idx[, 1] + 0.5,
      ytop = na_idx[, 2] + 0.5,
      col = "grey90",
      border = NA
    )
  }
  abline(h = seq_len(ncol(z)) + 0.5, v = seq_len(nrow(z)) + 0.5, col = "white", lwd = 0.4)

  par(mar = c(7, 1, 3, 4))
  legend_vals <- seq(-3, 3, length.out = 101)
  legend_mat <- rbind(legend_vals, legend_vals)
  graphics::image(
    x = c(0, 1),
    y = legend_vals,
    z = legend_mat,
    col = cols,
    breaks = breaks,
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = ""
  )
  axis(4, at = c(-3, -1.5, 0, 1.5, 3), las = 2, cex.axis = 0.8)
  mtext("robust z", side = 4, line = 2.2, cex = 0.85)

  invisible(NULL)
}

compute_top_expression_metrics <- function(counts, top_n = 50L) {
  total_counts <- Matrix::colSums(counts)
  dominant_gene_fraction <- rep(NA_real_, ncol(counts))
  percent_topn <- rep(NA_real_, ncol(counts))

  if (inherits(counts, "dgCMatrix")) {
    ptr <- counts@p
    vals_all <- counts@x
    for (j in seq_len(ncol(counts))) {
      total_j <- total_counts[j]
      if (!is.finite(total_j) || total_j <= 0) {
        dominant_gene_fraction[j] <- 0
        percent_topn[j] <- 0
        next
      }
      start <- ptr[j] + 1L
      end <- ptr[j + 1L]
      vals <- if (start <= end) vals_all[start:end] else numeric(0)
      if (length(vals) == 0) {
        dominant_gene_fraction[j] <- 0
        percent_topn[j] <- 0
        next
      }
      vals_sorted <- sort(vals, decreasing = TRUE)
      dominant_gene_fraction[j] <- vals_sorted[1] / total_j
      percent_topn[j] <- 100 * sum(head(vals_sorted, n = min(top_n, length(vals_sorted)))) / total_j
    }
  } else {
    for (j in seq_len(ncol(counts))) {
      vals <- as.numeric(counts[, j, drop = TRUE])
      total_j <- sum(vals, na.rm = TRUE)
      if (!is.finite(total_j) || total_j <= 0) {
        dominant_gene_fraction[j] <- 0
        percent_topn[j] <- 0
        next
      }
      vals <- vals[vals > 0]
      if (length(vals) == 0) {
        dominant_gene_fraction[j] <- 0
        percent_topn[j] <- 0
        next
      }
      vals_sorted <- sort(vals, decreasing = TRUE)
      dominant_gene_fraction[j] <- vals_sorted[1] / total_j
      percent_topn[j] <- 100 * sum(head(vals_sorted, n = min(top_n, length(vals_sorted)))) / total_j
    }
  }

  data.frame(
    cell = colnames(counts),
    dominant_gene_fraction = dominant_gene_fraction,
    percent.top50 = percent_topn,
    stringsAsFactors = FALSE
  )
}

build_cluster_summary <- function(df, cluster_col, metrics) {
  parts <- lapply(metrics, function(metric_i) {
    df %>%
      group_by(cluster = .data[[cluster_col]]) %>%
      summarise(
        metric = metric_i,
        n_cells = n(),
        mean_value = mean(.data[[metric_i]], na.rm = TRUE),
        median_value = median(.data[[metric_i]], na.rm = TRUE),
        sd_value = stats::sd(.data[[metric_i]], na.rm = TRUE),
        iqr_value = IQR(.data[[metric_i]], na.rm = TRUE),
        q05 = stats::quantile(.data[[metric_i]], probs = 0.05, na.rm = TRUE, names = FALSE),
        q25 = stats::quantile(.data[[metric_i]], probs = 0.25, na.rm = TRUE, names = FALSE),
        q75 = stats::quantile(.data[[metric_i]], probs = 0.75, na.rm = TRUE, names = FALSE),
        q95 = stats::quantile(.data[[metric_i]], probs = 0.95, na.rm = TRUE, names = FALSE),
        .groups = "drop"
      )
  })
  bind_rows(parts)
}

compute_global_metric_tests <- function(df, cluster_col, metrics) {
  bind_rows(lapply(metrics, function(metric_i) {
    sub_df <- df %>%
      select(cluster = all_of(cluster_col), value = all_of(metric_i)) %>%
      filter(!is.na(value), !is.na(cluster))
    if (nrow(sub_df) == 0 || dplyr::n_distinct(sub_df$cluster) < 2) {
      return(data.frame(
        metric = metric_i,
        n_cells = nrow(sub_df),
        n_clusters = dplyr::n_distinct(sub_df$cluster),
        kruskal_statistic = NA_real_,
        kruskal_df = NA_real_,
        kruskal_p = NA_real_,
        epsilon_squared = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    kt <- stats::kruskal.test(value ~ cluster, data = sub_df)
    h <- as.numeric(kt$statistic)
    n <- nrow(sub_df)
    k <- dplyr::n_distinct(sub_df$cluster)
    epsilon_squared <- if (n > k) max(0, (h - k + 1) / (n - k)) else NA_real_
    data.frame(
      metric = metric_i,
      n_cells = n,
      n_clusters = k,
      kruskal_statistic = h,
      kruskal_df = as.numeric(kt$parameter),
      kruskal_p = as.numeric(kt$p.value),
      epsilon_squared = epsilon_squared,
      stringsAsFactors = FALSE
    )
  }))
}

compute_cluster_metric_association <- function(df, cluster_col, metrics) {
  cluster_levels <- unique(as.character(df[[cluster_col]]))
  out <- list()
  idx <- 1L

  for (cluster_i in cluster_levels) {
    in_cluster <- as.character(df[[cluster_col]]) == cluster_i
    for (metric_i in metrics) {
      x <- df[[metric_i]][in_cluster]
      y <- df[[metric_i]][!in_cluster]
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]

      if (length(x) == 0 || length(y) == 0) {
        out[[idx]] <- data.frame(
          cluster = cluster_i,
          metric = metric_i,
          n_cluster = length(x),
          n_rest = length(y),
          median_cluster = NA_real_,
          median_rest = NA_real_,
          mean_cluster = NA_real_,
          mean_rest = NA_real_,
          median_diff = NA_real_,
          wilcox_p = NA_real_,
          cles = NA_real_,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
        next
      }

      wt <- suppressWarnings(stats::wilcox.test(x, y, exact = FALSE))
      w_stat <- as.numeric(wt$statistic)
      u_stat <- w_stat - length(x) * (length(x) + 1) / 2
      cles <- u_stat / (length(x) * length(y))

      out[[idx]] <- data.frame(
        cluster = cluster_i,
        metric = metric_i,
        n_cluster = length(x),
        n_rest = length(y),
        median_cluster = median(x, na.rm = TRUE),
        median_rest = median(y, na.rm = TRUE),
        mean_cluster = mean(x, na.rm = TRUE),
        mean_rest = mean(y, na.rm = TRUE),
        median_diff = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
        wilcox_p = as.numeric(wt$p.value),
        cles = cles,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  bind_rows(out)
}

compute_cluster_logistic_qc <- function(df, cluster_col, metrics) {
  cluster_levels <- unique(as.character(df[[cluster_col]]))
  out <- list()
  idx <- 1L

  model_df <- df %>%
    select(cluster = all_of(cluster_col), all_of(metrics)) %>%
    mutate(cluster = as.character(cluster))

  for (cluster_i in cluster_levels) {
    sub_df <- model_df %>%
      mutate(target = as.integer(cluster == cluster_i)) %>%
      select(-cluster)

    keep <- stats::complete.cases(sub_df)
    sub_df <- sub_df[keep, , drop = FALSE]
    if (nrow(sub_df) < 50 || sum(sub_df$target == 1) < 10 || sum(sub_df$target == 0) < 10) {
      out[[idx]] <- data.frame(
        cluster = cluster_i,
        n_cells_used = nrow(sub_df),
        n_target = sum(sub_df$target == 1),
        n_rest = sum(sub_df$target == 0),
        mcfadden_r2 = NA_real_,
        model_aic = NA_real_,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
      next
    }

    for (metric_i in metrics) {
      sub_df[[metric_i]] <- as.numeric(scale(sub_df[[metric_i]]))
    }

    full_formula <- stats::as.formula(
      paste("target ~", paste(metrics, collapse = " + "))
    )

    fit_full <- tryCatch(
      stats::glm(full_formula, data = sub_df, family = stats::binomial()),
      error = function(e) NULL,
      warning = function(w) suppressWarnings(stats::glm(full_formula, data = sub_df, family = stats::binomial()))
    )
    fit_null <- tryCatch(
      stats::glm(target ~ 1, data = sub_df, family = stats::binomial()),
      error = function(e) NULL
    )

    if (is.null(fit_full) || is.null(fit_null)) {
      out[[idx]] <- data.frame(
        cluster = cluster_i,
        n_cells_used = nrow(sub_df),
        n_target = sum(sub_df$target == 1),
        n_rest = sum(sub_df$target == 0),
        mcfadden_r2 = NA_real_,
        model_aic = NA_real_,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
      next
    }

    ll_full <- as.numeric(stats::logLik(fit_full))
    ll_null <- as.numeric(stats::logLik(fit_null))
    mcfadden_r2 <- if (is.finite(ll_full) && is.finite(ll_null) && ll_null != 0) 1 - (ll_full / ll_null) else NA_real_

    out[[idx]] <- data.frame(
      cluster = cluster_i,
      n_cells_used = nrow(sub_df),
      n_target = sum(sub_df$target == 1),
      n_rest = sum(sub_df$target == 0),
      mcfadden_r2 = mcfadden_r2,
      model_aic = stats::AIC(fit_full),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }

  bind_rows(out)
}

make_cluster_flag_table <- function(cluster_summary_long, logistic_df, threshold_abs_z = 1.5) {
  suspicious_rules <- c(
    nCount_RNA = "low",
    nFeature_RNA = "low",
    percent.mt = "high",
    log10_genes_per_umi = "low",
    umi_per_gene = "high",
    dominant_gene_fraction = "high",
    percent.top50 = "high",
    percent.ribo = "high"
  )

  flag_long <- cluster_summary_long %>%
    mutate(
      suspicious_direction = suspicious_rules[metric]
    ) %>%
    group_by(metric) %>%
    mutate(
      across_cluster_median = median(median_value, na.rm = TRUE),
      across_cluster_mad = safe_mad(median_value),
      robust_z = ifelse(
        is.na(across_cluster_mad),
        NA_real_,
        (median_value - across_cluster_median) / across_cluster_mad
      ),
      suspicious_flag = case_when(
        suspicious_direction == "high" & !is.na(robust_z) & robust_z >= threshold_abs_z ~ TRUE,
        suspicious_direction == "low" & !is.na(robust_z) & robust_z <= -threshold_abs_z ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    ungroup()

  flag_summary <- flag_long %>%
    group_by(cluster) %>%
    summarise(
      n_flagged_metrics = sum(suspicious_flag, na.rm = TRUE),
      flagged_metrics = paste(metric[suspicious_flag], collapse = ";"),
      .groups = "drop"
    ) %>%
    left_join(logistic_df, by = "cluster") %>%
    mutate(
      flagged_metrics = ifelse(flagged_metrics == "", NA_character_, flagged_metrics),
      concern_level = case_when(
        n_flagged_metrics >= 3 ~ "High",
        n_flagged_metrics >= 2 ~ "Moderate",
        !is.na(mcfadden_r2) & mcfadden_r2 >= 0.10 ~ "Moderate",
        TRUE ~ "Low"
      ),
      interpretation = case_when(
        concern_level == "High" ~ "Multiple QC metrics are shifted in suspicious directions; cluster separation may be influenced by sequencing quality.",
        concern_level == "Moderate" ~ "Some QC metrics differ enough to warrant caution when interpreting this cluster.",
        TRUE ~ "No strong QC-driven separation signal was detected from the selected metrics."
      )
    )

  list(flag_long = flag_long, flag_summary = flag_summary)
}

plot_umap_metric <- function(plot_df, metric_name, file_stub) {
  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[metric_name]])) +
    geom_point(size = 0.30, alpha = 0.85) +
    scale_color_gradient(low = "grey90", high = "#b2182b", na.value = "grey80") +
    theme_bw(base_size = 11) +
    labs(title = paste("UMAP colored by", metric_name), color = metric_name)
  save_plot_pdf_png(p, file_stub, width = 8, height = 6.5)
}

.ensure_dir(output_root)
out_plots <- .ensure_dir(file.path(output_root, "plots"))
out_summary <- .ensure_dir(file.path(output_root, "summaries"))
out_objects <- .ensure_dir(file.path(output_root, "objects"))

  if (!file.exists(input_rds)) {
    stop("Input RDS does not exist: ", input_rds, call. = FALSE)
  }

  message("Reading Seurat object: ", input_rds)
  obj <- readRDS(input_rds)
  if (!inherits(obj, "Seurat")) {
    stop("Input file is not a Seurat object.", call. = FALSE)
  }
  assert_required_column(obj@meta.data, base_cluster_col)
  if (!("umap" %in% names(obj@reductions))) {
    stop("UMAP reduction not found in Seurat object.", call. = FALSE)
  }
  if (!(assay_use %in% names(obj@assays))) {
    stop("Assay not found in Seurat object: ", assay_use, call. = FALSE)
  }

  obj <- maybe_join_layers(obj, assay = assay_use)
  counts <- get_assay_matrix(obj, assay = assay_use, slot_name = "counts")
  if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) {
    stop("Counts matrix is unavailable in assay ", assay_use, ".", call. = FALSE)
  }

  message("Computing per-cell QC metrics.")
  cluster_vec <- as.character(obj@meta.data[[base_cluster_col]])
  cluster_levels <- if (is.factor(obj@meta.data[[base_cluster_col]])) as.character(levels(obj@meta.data[[base_cluster_col]])) else sort_maybe_numeric(cluster_vec)
  obj@meta.data[[base_cluster_col]] <- factor(cluster_vec, levels = cluster_levels)

  if (!("percent.mt" %in% colnames(obj@meta.data))) {
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, assay = assay_use, pattern = "^MT-")
  }

  percent_ribo <- if (any(grepl("^RPL|^RPS|^MRPL|^MRPS", rownames(counts), ignore.case = TRUE))) {
    Seurat::PercentageFeatureSet(obj, assay = assay_use, pattern = "^RPL|^RPS|^MRPL|^MRPS")
  } else {
    rep(NA_real_, ncol(obj))
  }

  nCount_RNA <- if ("nCount_RNA" %in% colnames(obj@meta.data)) as.numeric(obj@meta.data[["nCount_RNA"]]) else as.numeric(Matrix::colSums(counts))
  nFeature_RNA <- if ("nFeature_RNA" %in% colnames(obj@meta.data)) as.numeric(obj@meta.data[["nFeature_RNA"]]) else as.numeric(Matrix::colSums(counts > 0))

  top_metrics_df <- compute_top_expression_metrics(counts, top_n = top_n_genes_fraction)
  umap_df <- as.data.frame(Embeddings(obj, reduction = "umap")[, 1:2, drop = FALSE])
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")

  qc_df <- data.frame(
    cell = colnames(obj),
    cluster = as.character(obj@meta.data[[base_cluster_col]]),
    nCount_RNA = nCount_RNA,
    nFeature_RNA = nFeature_RNA,
    percent.mt = as.numeric(obj@meta.data[["percent.mt"]]),
    percent.ribo = as.numeric(percent_ribo),
    log10_genes_per_umi = ifelse(nCount_RNA > 1 & nFeature_RNA > 1, log10(nFeature_RNA) / log10(nCount_RNA), NA_real_),
    umi_per_gene = ifelse(nFeature_RNA > 0, nCount_RNA / nFeature_RNA, NA_real_),
    stringsAsFactors = FALSE
  ) %>%
    left_join(top_metrics_df, by = "cell") %>%
    bind_cols(umap_df)

  if ("S.Score" %in% colnames(obj@meta.data)) {
    qc_df$S.Score <- as.numeric(obj@meta.data[["S.Score"]])
  }
  if ("G2M.Score" %in% colnames(obj@meta.data)) {
    qc_df$G2M.Score <- as.numeric(obj@meta.data[["G2M.Score"]])
  }
  if ("Phase" %in% colnames(obj@meta.data)) {
    qc_df$Phase <- as.character(obj@meta.data[["Phase"]])
  }

  obj@meta.data[["percent.ribo"]] <- qc_df$percent.ribo
  obj@meta.data[["log10_genes_per_umi"]] <- qc_df$log10_genes_per_umi
  obj@meta.data[["umi_per_gene"]] <- qc_df$umi_per_gene
  obj@meta.data[["dominant_gene_fraction"]] <- qc_df$dominant_gene_fraction
  obj@meta.data[["percent.top50"]] <- qc_df$percent.top50

  readr::write_csv(qc_df, file.path(out_summary, "cell_qc_metrics.csv"))

  qc_metrics <- c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt",
    "log10_genes_per_umi",
    "umi_per_gene",
    "dominant_gene_fraction",
    "percent.top50"
  )
  if (all(is.na(qc_df$percent.ribo)) == FALSE) {
    qc_metrics <- c(qc_metrics, "percent.ribo")
  }

  message("Building cluster-level summaries and tests.")
  cluster_summary_long <- build_cluster_summary(qc_df, "cluster", qc_metrics)
  readr::write_csv(cluster_summary_long, file.path(out_summary, "cluster_qc_summary_long.csv"))

  cluster_summary_wide <- cluster_summary_long %>%
    select(cluster, metric, n_cells, mean_value, median_value, sd_value, iqr_value) %>%
    pivot_wider(
      names_from = metric,
      values_from = c(mean_value, median_value, sd_value, iqr_value),
      names_sep = "__"
    )
  readr::write_csv(cluster_summary_wide, file.path(out_summary, "cluster_qc_summary.csv"))

  global_metric_tests <- compute_global_metric_tests(qc_df, "cluster", qc_metrics)
  readr::write_csv(global_metric_tests, file.path(out_summary, "cluster_qc_global_metric_tests.csv"))

  cluster_metric_assoc <- compute_cluster_metric_association(qc_df, "cluster", qc_metrics)
  readr::write_csv(cluster_metric_assoc, file.path(out_summary, "cluster_vs_qc_association_long.csv"))

  logistic_qc_df <- compute_cluster_logistic_qc(qc_df, "cluster", qc_metrics)
  readr::write_csv(logistic_qc_df, file.path(out_summary, "cluster_qc_logistic_multivariate.csv"))

  flag_tables <- make_cluster_flag_table(cluster_summary_long, logistic_qc_df, threshold_abs_z = flag_threshold_abs_robust_z)
  readr::write_csv(flag_tables$flag_long, file.path(out_summary, "cluster_qc_metric_flags_long.csv"))
  readr::write_csv(flag_tables$flag_summary, file.path(out_summary, "cluster_qc_outlier_flags.csv"))

  message("Writing QC plots.")
  p_cluster <- DimPlot(
    obj,
    reduction = "umap",
    group.by = base_cluster_col,
    label = TRUE,
    repel = TRUE,
    raster = TRUE,
    pt.size = 0.30
  ) + labs(title = paste("UMAP by", base_cluster_col))
  save_plot_pdf_png(p_cluster, file.path(out_plots, "umap_by_refined_cluster"), width = 9, height = 7)

  for (metric_i in qc_metrics) {
    plot_umap_metric(qc_df, metric_i, file.path(out_plots, paste0("umap_", metric_i)))
  }

  violin_df <- qc_df %>%
    select(cluster, all_of(qc_metrics)) %>%
    pivot_longer(cols = all_of(qc_metrics), names_to = "metric", values_to = "value") %>%
    filter(is.finite(value))

  p_violin <- ggplot(violin_df, aes(x = factor(cluster, levels = cluster_levels), y = value, fill = factor(cluster, levels = cluster_levels))) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.7) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(title = "QC metric distributions by refined cluster", x = "Cluster", y = "Value")
  save_plot_pdf_png(p_violin, file.path(out_plots, "vln_qc_metrics_by_cluster"), width = max(12, 0.55 * length(cluster_levels) + 7), height = 12)

  heatmap_df <- flag_tables$flag_long %>%
    select(cluster, metric, robust_z) %>%
    mutate(cluster = factor(cluster, levels = cluster_levels)) %>%
    arrange(cluster, metric)
  heatmap_mat <- heatmap_df %>%
    pivot_wider(names_from = metric, values_from = robust_z) %>%
    tibble::column_to_rownames("cluster") %>%
    as.matrix()
  plot_simple_heatmap(
    mat = heatmap_mat,
    file_pdf = file.path(out_plots, "heatmap_cluster_qc_robust_zscores.pdf"),
    title_text = "Cluster QC median shifts (robust z-scores)"
  )

  focus_clusters_present <- intersect(focus_clusters, unique(as.character(qc_df$cluster)))
  if (length(focus_clusters_present) > 0) {
    focus_diag_df <- bind_rows(lapply(focus_clusters_present, function(cluster_i) {
      qc_df %>%
        mutate(group = ifelse(cluster == cluster_i, cluster_i, "rest")) %>%
        group_by(group) %>%
        summarise(
          n_cells = n(),
          across(all_of(qc_metrics), list(mean = ~ mean(.x, na.rm = TRUE), median = ~ median(.x, na.rm = TRUE))),
          .groups = "drop"
        ) %>%
        mutate(focus_cluster = cluster_i)
    }))
    readr::write_csv(focus_diag_df, file.path(out_summary, "cluster9_qc_diagnostic.csv"))

    focus_plot_df <- qc_df %>%
      filter(cluster %in% focus_clusters_present | !(cluster %in% focus_clusters_present)) %>%
      mutate(group = ifelse(cluster %in% focus_clusters_present, as.character(cluster), "rest")) %>%
      filter(group %in% c("rest", focus_clusters_present)) %>%
      pivot_longer(cols = all_of(qc_metrics), names_to = "metric", values_to = "value") %>%
      filter(is.finite(value))

    p_focus <- ggplot(focus_plot_df, aes(x = group, y = value, fill = group)) +
      geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
      geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.7) +
      facet_wrap(~ metric, scales = "free_y", ncol = 2) +
      theme_bw(base_size = 10) +
      theme(legend.position = "none") +
      labs(title = "Cluster 9 / 9c QC diagnostics versus rest", x = NULL, y = "Value")
    save_plot_pdf_png(p_focus, file.path(out_plots, "cluster9_vs_rest_qc_comparison"), width = 11, height = 12)
  }

  saveRDS(obj, file.path(out_objects, "integrated_sct_cca_seurat_cluster_refine_qc.rds"))

  summary_lines <- c(
    "Cluster quality QC workflow completed.",
    paste0("Input object: ", input_rds),
    paste0("Output root: ", output_root),
    paste0("Base cluster column: ", base_cluster_col),
    paste0("Core QC metrics: ", paste(qc_metrics, collapse = ", ")),
    paste0("Top-gene concentration metric: percent.top", top_n_genes_fraction),
    "",
    "Interpretation guide:",
    "Clusters with lower nCount_RNA, lower nFeature_RNA, lower log10_genes_per_umi, and higher percent.mt / umi_per_gene / dominant_gene_fraction / percent.top50 are more suspicious for quality-driven separation.",
    "The file cluster_qc_outlier_flags.csv provides a first-pass concern level for each cluster.",
    "The file cluster_qc_logistic_multivariate.csv reports how well QC metrics alone predict cluster membership in one-vs-rest models."
  )
  writeLines(summary_lines, con = file.path(out_summary, "run_summary.txt"))

message("Cluster quality workflow finished.")
message("Outputs written to: ", normalizePath(output_root))
