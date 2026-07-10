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
  candidate_files <- candidate_files[!is.na(candidate_files) & nzchar(candidate_files)]
  candidate_dirs <- unique(dirname(normalizePath(candidate_files, mustWork = FALSE)))
  cwd <- normalizePath(getwd(), mustWork = FALSE)
  candidate_dirs <- unique(c(candidate_dirs, cwd, file.path(cwd, "Code", "in-vivo")))
  hit <- candidate_dirs[file.exists(file.path(candidate_dirs, "Utils.R"))]
  if (length(hit) > 0) return(normalizePath(hit[1], mustWork = TRUE))
  stop("Cannot locate Code/in-vivo/Utils.R from script path or working directory: ", cwd, call. = FALSE)
}

script_dir <- resolve_in_vivo_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
source(file.path(script_dir, "Utils.R"))

required_packages <- c("dplyr", "ggplot2", "readr", "tidyr", "tibble", "patchwork", "sandwich", "lmtest", "jsonlite", "hdf5r")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(sandwich)
  library(lmtest)
  library(jsonlite)
  library(hdf5r)
})

set.seed(1234)

get_env_scalar <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(trimws(value))) return(default)
  value
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

as_clean_chr_04f <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "None")] <- NA_character_
  x
}

first_present_value <- function(x) {
  x <- as_clean_chr_04f(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) NA_character_ else x[1]
}

first_present_num <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else x[1]
}

ensure_dir_04f <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = FALSE)
}

write_csv_04f <- function(df, file) {
  ensure_dir_04f(dirname(file))
  if (tolower(tools::file_ext(file)) %in% c("tsv", "txt")) {
    readr::write_tsv(as.data.frame(df, stringsAsFactors = FALSE), file)
  } else {
    readr::write_csv(as.data.frame(df, stringsAsFactors = FALSE), file)
  }
  invisible(file)
}

save_plot_04f <- function(plot_obj, file, width = 8, height = 6) {
  ensure_dir_04f(dirname(file))
  ggplot2::ggsave(file, plot_obj, width = width, height = height, limitsize = FALSE)
  invisible(file)
}

save_plot_both_04f <- function(plot_obj, file_stub, width = 8, height = 6) {
  save_plot_04f(plot_obj, paste0(file_stub, ".pdf"), width = width, height = height)
  save_plot_04f(plot_obj, paste0(file_stub, ".png"), width = width, height = height)
  invisible(file_stub)
}

theme_04f <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(color = "grey35"),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = "grey94", color = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold")
    )
}

make_cluster_palette_04g <- function(clusters) {
  clusters <- as.character(clusters)
  if (length(clusters) == 0) return(stats::setNames(character(0), character(0)))
  colors <- grDevices::hcl.colors(max(length(clusters), 3), palette = "Dark 3")
  stats::setNames(colors[seq_along(clusters)], clusters)
}

cluster_facet_nrow_04g <- function(n_clusters) {
  if (n_clusters <= 6) 1 else 2
}

standardize_ploidy_04f <- function(x) {
  x <- toupper(as_clean_chr_04f(x))
  x <- gsub("\\s+", "", x)
  x[x %in% c("2", "2.0", "2N")] <- "2N"
  x[x %in% c("4", "4.0", "4N")] <- "4N"
  x[!(x %in% c("2N", "4N"))] <- NA_character_
  x
}

standardize_dose_04f <- function(x) {
  if (exists("standardize_in_vivo_dose", mode = "function")) {
    out <- standardize_in_vivo_dose(x)
  } else {
    out <- as_clean_chr_04f(x)
  }
  out <- as_clean_chr_04f(out)
  out[out %in% c("0", "0mg", "0 mg/kg", "0mg/kg")] <- "0mg/kg"
  out[out %in% c("30", "30mg", "30 mg/kg", "30mg/kg")] <- "30mg/kg"
  out[out %in% c("120", "120mg", "120 mg/kg", "120mg/kg")] <- "120mg/kg"
  out
}

dose_panel_label <- function(x) {
  out <- standardize_dose_04f(x)
  out <- sub("mg/kg$", "", out)
  out[!(out %in% c("0", "30", "120"))] <- NA_character_
  out
}

choose_time_bin_width_04f <- function(x, n_bins = 200L) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (length(x) < 2) return(0.01)
  rng <- range(x, na.rm = TRUE)
  width <- diff(rng) / n_bins
  if (!is.finite(width) || width <= 0) width <- 0.01
  width
}

make_sample_histogram_data <- function(df, sample_info, sample_order, n_bins = 200L) {
  df <- df %>%
    dplyr::mutate(.time = safe_num(.data$pseudotime)) %>%
    dplyr::filter(is.finite(.data$.time), !is.na(.data$sampleID), nzchar(.data$sampleID))
  if (nrow(df) == 0) return(data.frame())

  x_range <- range(df$.time, na.rm = TRUE)
  if (!all(is.finite(x_range)) || diff(x_range) <= 0) {
    x_range <- c(min(df$.time, na.rm = TRUE) - 0.005, max(df$.time, na.rm = TRUE) + 0.005)
  }
  bin_width <- choose_time_bin_width_04f(df$.time, n_bins = n_bins)
  breaks <- seq(x_range[1], x_range[2] + bin_width, by = bin_width)
  if (length(breaks) < 3) breaks <- seq(x_range[1], x_range[2] + 0.01, length.out = 3)
  n_bin <- length(breaks) - 1L
  df$.bin <- findInterval(df$.time, breaks, rightmost.closed = TRUE, all.inside = TRUE)
  df$.bin <- pmax(1L, pmin(n_bin, df$.bin))

  bin_df <- data.frame(
    .bin = seq_len(n_bin),
    time_bin_left = breaks[seq_len(n_bin)],
    time_bin_right = breaks[seq_len(n_bin) + 1L],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      time_bin_mid = (.data$time_bin_left + .data$time_bin_right) / 2,
      bin_width = .data$time_bin_right - .data$time_bin_left
    )

  sample_meta <- sample_info %>%
    dplyr::filter(.data$sampleID %in% unique(df$sampleID)) %>%
    dplyr::mutate(sampleID = factor(as.character(.data$sampleID), levels = sample_order)) %>%
    dplyr::arrange(.data$sampleID) %>%
    dplyr::mutate(sampleID = as.character(.data$sampleID))

  counts <- df %>%
    dplyr::count(.data$sampleID, .data$Ploidy, .data$Dose, .data$DosePanel, .data$.bin, name = "n_cells")
  totals <- df %>%
    dplyr::count(.data$sampleID, name = "total_cells_in_sample")

  tidyr::crossing(
    sample_meta %>% dplyr::select("sampleID", "Ploidy", "Dose", "DosePanel", "sampleFacetLabel", "mean_cell_ploidy"),
    bin_df
  ) %>%
    dplyr::left_join(counts, by = c("sampleID", "Ploidy", "Dose", "DosePanel", ".bin")) %>%
    dplyr::left_join(totals, by = "sampleID") %>%
    dplyr::mutate(
      n_cells = dplyr::coalesce(.data$n_cells, 0L),
      cell_percent = 100 * .data$n_cells / .data$total_cells_in_sample,
      sampleID = factor(as.character(.data$sampleID), levels = sample_order),
      sampleFacetLabel = factor(as.character(.data$sampleFacetLabel), levels = sample_info$sampleFacetLabel[match(sample_order, sample_info$sampleID)]),
      Ploidy = factor(as.character(.data$Ploidy), levels = c("2N", "4N")),
      DosePanel = factor(as.character(.data$DosePanel), levels = c("0", "30", "120"))
    )
}

plot_sample_histogram_16x1 <- function(hist_df, title, file_stub, height = NULL) {
  if (nrow(hist_df) == 0) return(invisible(NULL))
  y_limit <- max(hist_df$cell_percent, na.rm = TRUE)
  if (!is.finite(y_limit) || y_limit <= 0) y_limit <- 1
  x_limits <- range(c(hist_df$time_bin_left, hist_df$time_bin_right), na.rm = TRUE)
  n_samples <- length(unique(hist_df$sampleID))
  if (is.null(height)) height <- max(18, n_samples * 1.15)

  p <- ggplot(hist_df, aes(x = .data$time_bin_mid, y = .data$cell_percent, fill = .data$Ploidy, color = .data$Ploidy)) +
    geom_col(aes(width = .data$bin_width), position = "identity", alpha = 0.38, linewidth = 0.08) +
    geom_step(linewidth = 0.45) +
    facet_grid(rows = vars(sampleFacetLabel), scales = "fixed", drop = TRUE) +
    scale_x_continuous(limits = x_limits, expand = ggplot2::expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0, y_limit * 1.05), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = c("2N" = "#3B73B9", "4N" = "#D95F02"), drop = FALSE) +
    scale_color_manual(values = c("2N" = "#1D4E89", "4N" = "#A84600"), drop = FALSE) +
    labs(
      title = title,
      subtitle = paste0(analysis_scope_label, ". Y axis is percent of selected cells within each sample."),
      x = "scVelo pseudotime",
      y = "Cell percent",
      fill = "Initial Ploidy",
      color = "Initial Ploidy"
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE), color = "none") +
    theme_04f(6) +
    theme(
      legend.position = "bottom",
      strip.placement = "outside",
      strip.text.y.right = element_text(angle = 0, hjust = 0.5, size = 5.3),
      panel.spacing.y = grid::unit(0.28, "lines")
    )
  save_plot_both_04f(p, file_stub, width = 7.2, height = height)
  invisible(p)
}

plot_dose_histogram_panels <- function(hist_df, title, file_stub, sample_info, sample_order) {
  if (nrow(hist_df) == 0) return(invisible(NULL))
  x_limits <- range(c(hist_df$time_bin_left, hist_df$time_bin_right), na.rm = TRUE)
  y_limit <- max(hist_df$cell_percent, na.rm = TRUE)
  if (!is.finite(y_limit) || y_limit <= 0) y_limit <- 1

  dose_levels <- c("0", "30", "120")
  dose_plots <- lapply(dose_levels, function(dose_value) {
    dose_samples <- sample_info %>%
      dplyr::filter(.data$DosePanel == dose_value) %>%
      dplyr::filter(.data$sampleID %in% sample_order) %>%
      dplyr::mutate(sampleID = factor(.data$sampleID, levels = sample_order)) %>%
      dplyr::arrange(.data$sampleID) %>%
      dplyr::pull(.data$sampleID) %>%
      as.character()
    plot_df <- hist_df %>%
      dplyr::filter(as.character(.data$DosePanel) == dose_value, as.character(.data$sampleID) %in% dose_samples) %>%
      dplyr::mutate(
        sampleID = factor(as.character(.data$sampleID), levels = dose_samples),
        sampleFacetLabel = factor(as.character(.data$sampleFacetLabel), levels = sample_info$sampleFacetLabel[match(dose_samples, sample_info$sampleID)])
      )
    if (nrow(plot_df) == 0) {
      return(ggplot() + theme_void() + labs(title = paste0(dose_value, " mg/kg")))
    }
    ggplot(plot_df, aes(x = .data$time_bin_mid, y = .data$cell_percent, fill = .data$Ploidy, color = .data$Ploidy)) +
      geom_col(aes(width = .data$bin_width), position = "identity", alpha = 0.38, linewidth = 0.08) +
      geom_step(linewidth = 0.45) +
      facet_grid(rows = vars(sampleFacetLabel), scales = "fixed", drop = TRUE) +
      scale_x_continuous(limits = x_limits, expand = ggplot2::expansion(mult = c(0, 0))) +
      scale_y_continuous(limits = c(0, y_limit * 1.05), expand = ggplot2::expansion(mult = c(0, 0.02))) +
      scale_fill_manual(values = c("2N" = "#3B73B9", "4N" = "#D95F02"), drop = FALSE) +
      scale_color_manual(values = c("2N" = "#1D4E89", "4N" = "#A84600"), drop = FALSE) +
      labs(title = paste0(dose_value, " mg/kg"), x = "scVelo pseudotime", y = "Cell percent", fill = "Initial Ploidy", color = "Initial Ploidy") +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE), color = "none") +
      theme_04f(6) +
      theme(
        legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.right = element_text(angle = 0, hjust = 0.5, size = 5.2),
        panel.spacing.y = grid::unit(0.28, "lines")
      )
  })

  p <- patchwork::wrap_plots(dose_plots, nrow = 1, guides = "collect") +
    patchwork::plot_annotation(
      title = title,
      subtitle = analysis_scope_label
    )
  save_plot_both_04f(p, file_stub, width = 14.5, height = 11)
  invisible(p)
}

make_density_data <- function(df, sample_info, n = 512L) {
  df <- df %>%
    dplyr::mutate(.time = safe_num(.data$pseudotime)) %>%
    dplyr::filter(is.finite(.data$.time), !is.na(.data$sampleID), nzchar(.data$sampleID))
  if (nrow(df) == 0) return(data.frame())
  x_range <- range(df$.time, na.rm = TRUE)
  if (!all(is.finite(x_range)) || diff(x_range) <= 0) return(data.frame())
  rows <- lapply(sort(unique(as.character(df$sampleID))), function(sid) {
    x <- df$.time[as.character(df$sampleID) == sid]
    x <- x[is.finite(x)]
    if (length(x) < 2 || diff(range(x)) <= 0) return(NULL)
    den <- stats::density(x, from = x_range[1], to = x_range[2], n = n, na.rm = TRUE)
    data.frame(sampleID = sid, pseudotime = den$x, density = den$y, stringsAsFactors = FALSE)
  })
  dplyr::bind_rows(rows) %>%
    dplyr::left_join(
      sample_info %>%
        dplyr::select("sampleID", "Ploidy", "Dose", "DosePanel", "mean_cell_ploidy", "sampleFacetLabel"),
      by = "sampleID"
    ) %>%
    dplyr::mutate(
      Ploidy = factor(as.character(.data$Ploidy), levels = c("2N", "4N")),
      DosePanel = factor(as.character(.data$DosePanel), levels = c("0", "30", "120"))
    )
}

plot_density_curves <- function(density_df, title, file_stub, facet_by_dose = FALSE) {
  if (nrow(density_df) == 0) return(invisible(NULL))
  p <- ggplot(density_df, aes(x = .data$pseudotime, y = .data$density, group = .data$sampleID, color = .data$mean_cell_ploidy)) +
    geom_line(linewidth = 0.65, alpha = 0.88) +
    scale_color_gradientn(
      colors = c("#313695", "#74ADD1", "#FEE08B", "#F46D43", "#A50026"),
      name = "Sample mean\nCell Ploidy"
    ) +
    labs(
      title = title,
      subtitle = paste0("One density curve per sample; ", analysis_scope_label, "."),
      x = "scVelo pseudotime",
      y = "Density"
    ) +
    theme_04f(10) +
    theme(legend.position = "right")
  if (isTRUE(facet_by_dose)) {
    p <- p + facet_wrap(~DosePanel, nrow = 1, labeller = as_labeller(c("0" = "0 mg/kg", "30" = "30 mg/kg", "120" = "120 mg/kg")))
  }
  save_plot_both_04f(p, file_stub, width = ifelse(isTRUE(facet_by_dose), 13, 8.8), height = 5.8)
  invisible(p)
}

safe_kruskal <- function(df, value_col, group_col, min_groups = 2L) {
  df <- df %>%
    dplyr::filter(is.finite(.data[[value_col]]), !is.na(.data[[group_col]]), nzchar(as.character(.data[[group_col]])))
  groups <- unique(as.character(df[[group_col]]))
  if (nrow(df) == 0 || length(groups) < min_groups) {
    return(data.frame(
      statistic = NA_real_,
      p_value = NA_real_,
      n = nrow(df),
      n_groups = length(groups),
      status = "skipped_insufficient_groups",
      stringsAsFactors = FALSE
    ))
  }
  kt <- tryCatch(stats::kruskal.test(df[[value_col]] ~ as.factor(df[[group_col]])), error = function(e) NULL)
  if (is.null(kt)) {
    return(data.frame(statistic = NA_real_, p_value = NA_real_, n = nrow(df), n_groups = length(groups), status = "error", stringsAsFactors = FALSE))
  }
  data.frame(
    statistic = unname(kt$statistic),
    p_value = kt$p.value,
    n = nrow(df),
    n_groups = length(groups),
    status = "ok",
    stringsAsFactors = FALSE
  )
}

dose_position_04f <- function(x) {
  x <- as.character(x)
  x <- dplyr::case_when(
    x %in% c("0", "0mg/kg") ~ "0mg/kg",
    x %in% c("30", "30mg/kg") ~ "30mg/kg",
    x %in% c("120", "120mg/kg") ~ "120mg/kg",
    TRUE ~ x
  )
  unname(match(x, c("0mg/kg", "30mg/kg", "120mg/kg")))
}

add_dose_x_04f <- function(df, dose_col = "dose") {
  df %>% dplyr::mutate(dose_x = dose_position_04f(.data[[dose_col]]))
}

dose_scale_04f <- function(labels = c("0mg/kg", "30mg/kg", "120mg/kg")) {
  ggplot2::scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = labels,
    limits = c(0.55, 3.45),
    expand = ggplot2::expansion(mult = c(0.02, 0.03))
  )
}

significance_label_04f <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.0001 ~ "****",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

safe_wilcox_p_04f <- function(x, y, min_n = 2L) {
  x <- safe_num(x)
  y <- safe_num(y)
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < min_n || length(y) < min_n) return(NA_real_)
  if (length(unique(c(x, y))) < 2) return(1)
  out <- tryCatch(
    stats::wilcox.test(x, y, exact = FALSE)$p.value,
    warning = function(w) suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )
  as.numeric(out)
}

group_offset_04f <- function(group_value, group_levels = c("2N", "4N"), dodge_width = 0.56) {
  idx <- match(as.character(group_value), group_levels)
  n <- length(group_levels)
  ifelse(is.na(idx), 0, (idx - (n + 1) / 2) * dodge_width / n)
}

make_group_significance_04f <- function(
  df,
  value_col,
  dose_col = "dose",
  group_col = "initial_ploidy_group",
  facet_cols = character(0),
  group_levels = c("2N", "4N"),
  group_label = "Initial Ploidy",
  hide_ns = FALSE,
  dodge_width = 0.56,
  y_pad_frac = 0.065,
  y_step_frac = 0.07
) {
  if (nrow(df) == 0) return(data.frame())
  stopifnot(value_col %in% names(df), dose_col %in% names(df), group_col %in% names(df))
  facet_cols <- intersect(facet_cols, names(df))
  plot_df <- as.data.frame(df, stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      .value = safe_num(.data[[value_col]]),
      .dose_chr = as.character(.data[[dose_col]]),
      .group_chr = as.character(.data[[group_col]]),
      .dose_pos = dose_position_04f(.data$.dose_chr)
    ) %>%
    dplyr::filter(is.finite(.data$.value), is.finite(.data$.dose_pos), .data$.group_chr %in% group_levels)
  if (nrow(plot_df) == 0) return(data.frame())
  plot_df$.facet_key <- if (length(facet_cols) == 0) {
    "all"
  } else {
    interaction(plot_df[, facet_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
  }

  rows <- list()
  row_i <- 1L
  for (facet_key in unique(as.character(plot_df$.facet_key))) {
    facet_df <- plot_df[as.character(plot_df$.facet_key) == facet_key, , drop = FALSE]
    if (nrow(facet_df) == 0) next
    y_range <- range(facet_df$.value, na.rm = TRUE)
    y_span <- diff(y_range)
    if (!is.finite(y_span) || y_span <= 0) y_span <- max(abs(y_range), na.rm = TRUE) * 0.1
    if (!is.finite(y_span) || y_span <= 0) y_span <- 1
    facet_values <- if (length(facet_cols) == 0) data.frame(stringsAsFactors = FALSE) else facet_df[1, facet_cols, drop = FALSE]

    comparisons <- list()
    comp_i <- 1L
    for (dose_value in c("0mg/kg", "30mg/kg", "120mg/kg")) {
      dose_matches <- facet_df$.dose_pos == dose_position_04f(dose_value)
      x <- facet_df$.value[dose_matches & facet_df$.group_chr == group_levels[1]]
      y <- facet_df$.value[dose_matches & facet_df$.group_chr == group_levels[2]]
      comparisons[[comp_i]] <- data.frame(
        comparison_type = paste0(group_label, "_within_dose"),
        comparison = paste(group_levels[1], "vs", group_levels[2], "at", dose_value),
        group_1 = group_levels[1],
        group_2 = group_levels[2],
        dose_1 = dose_value,
        dose_2 = dose_value,
        x = dose_position_04f(dose_value) + group_offset_04f(group_levels[1], group_levels, dodge_width),
        xend = dose_position_04f(dose_value) + group_offset_04f(group_levels[2], group_levels, dodge_width),
        xmid = dose_position_04f(dose_value),
        p_value = safe_wilcox_p_04f(x, y),
        n_1 = length(x[is.finite(x)]),
        n_2 = length(y[is.finite(y)]),
        stringsAsFactors = FALSE
      )
      comp_i <- comp_i + 1L
    }
    for (group_value in group_levels) {
      for (dose_value in c("30mg/kg", "120mg/kg")) {
        x <- facet_df$.value[facet_df$.dose_pos == dose_position_04f("0mg/kg") & facet_df$.group_chr == group_value]
        y <- facet_df$.value[facet_df$.dose_pos == dose_position_04f(dose_value) & facet_df$.group_chr == group_value]
        x0 <- dose_position_04f("0mg/kg") + group_offset_04f(group_value, group_levels, dodge_width)
        x1 <- dose_position_04f(dose_value) + group_offset_04f(group_value, group_levels, dodge_width)
        comparisons[[comp_i]] <- data.frame(
          comparison_type = "dose_vs_0mgkg_within_group",
          comparison = paste(group_value, dose_value, "vs 0mg/kg"),
          group_1 = group_value,
          group_2 = group_value,
          dose_1 = "0mg/kg",
          dose_2 = dose_value,
          x = x0,
          xend = x1,
          xmid = mean(c(x0, x1)),
          p_value = safe_wilcox_p_04f(x, y),
          n_1 = length(x[is.finite(x)]),
          n_2 = length(y[is.finite(y)]),
          stringsAsFactors = FALSE
        )
        comp_i <- comp_i + 1L
      }
    }
    comp_df <- dplyr::bind_rows(comparisons)
    if (nrow(comp_df) == 0) next
    comp_df$p_adjust_BH <- stats::p.adjust(comp_df$p_value, method = "BH")
    comp_df$label <- significance_label_04f(comp_df$p_adjust_BH)
    if (isTRUE(hide_ns)) comp_df <- comp_df %>% dplyr::filter(!is.na(.data$label), .data$label != "ns")
    if (nrow(comp_df) == 0) next
    comp_df <- comp_df %>%
      dplyr::mutate(
        annotation_index = dplyr::row_number(),
        y = y_range[2] + y_span * (y_pad_frac + (.data$annotation_index - 1) * y_step_frac),
        ytick = .data$y - y_span * 0.018,
        ylabel = .data$y + y_span * 0.018,
        value_col = value_col
      )
    if (length(facet_cols) > 0) {
      comp_df <- dplyr::bind_cols(facet_values[rep(1, nrow(comp_df)), , drop = FALSE], comp_df)
    }
    rows[[row_i]] <- comp_df
    row_i <- row_i + 1L
  }
  dplyr::bind_rows(rows)
}

make_dose_significance_04f <- function(
  df,
  value_col,
  dose_col = "dose",
  facet_cols = character(0),
  hide_ns = FALSE,
  y_pad_frac = 0.065,
  y_step_frac = 0.07
) {
  if (nrow(df) == 0) return(data.frame())
  stopifnot(value_col %in% names(df), dose_col %in% names(df))
  facet_cols <- intersect(facet_cols, names(df))
  plot_df <- as.data.frame(df, stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      .value = safe_num(.data[[value_col]]),
      .dose_chr = as.character(.data[[dose_col]]),
      .dose_pos = dose_position_04f(.data$.dose_chr)
    ) %>%
    dplyr::filter(is.finite(.data$.value), is.finite(.data$.dose_pos))
  if (nrow(plot_df) == 0) return(data.frame())
  plot_df$.facet_key <- if (length(facet_cols) == 0) {
    "all"
  } else {
    interaction(plot_df[, facet_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
  }
  dose_pairs <- list(c("0mg/kg", "30mg/kg"), c("0mg/kg", "120mg/kg"), c("30mg/kg", "120mg/kg"))
  rows <- list()
  row_i <- 1L
  for (facet_key in unique(as.character(plot_df$.facet_key))) {
    facet_df <- plot_df[as.character(plot_df$.facet_key) == facet_key, , drop = FALSE]
    y_range <- range(facet_df$.value, na.rm = TRUE)
    y_span <- diff(y_range)
    if (!is.finite(y_span) || y_span <= 0) y_span <- max(abs(y_range), na.rm = TRUE) * 0.1
    if (!is.finite(y_span) || y_span <= 0) y_span <- 1
    facet_values <- if (length(facet_cols) == 0) data.frame(stringsAsFactors = FALSE) else facet_df[1, facet_cols, drop = FALSE]
    comp_df <- dplyr::bind_rows(lapply(dose_pairs, function(pair) {
      x <- facet_df$.value[facet_df$.dose_pos == dose_position_04f(pair[1])]
      y <- facet_df$.value[facet_df$.dose_pos == dose_position_04f(pair[2])]
      data.frame(
        comparison_type = "dose_pairwise",
        comparison = paste(pair[2], "vs", pair[1]),
        dose_1 = pair[1],
        dose_2 = pair[2],
        x = dose_position_04f(pair[1]),
        xend = dose_position_04f(pair[2]),
        xmid = mean(c(dose_position_04f(pair[1]), dose_position_04f(pair[2]))),
        p_value = safe_wilcox_p_04f(x, y),
        n_1 = length(x[is.finite(x)]),
        n_2 = length(y[is.finite(y)]),
        stringsAsFactors = FALSE
      )
    }))
    if (nrow(comp_df) == 0) next
    comp_df$p_adjust_BH <- stats::p.adjust(comp_df$p_value, method = "BH")
    comp_df$label <- significance_label_04f(comp_df$p_adjust_BH)
    if (isTRUE(hide_ns)) comp_df <- comp_df %>% dplyr::filter(!is.na(.data$label), .data$label != "ns")
    if (nrow(comp_df) == 0) next
    comp_df <- comp_df %>%
      dplyr::mutate(
        annotation_index = dplyr::row_number(),
        y = y_range[2] + y_span * (y_pad_frac + (.data$annotation_index - 1) * y_step_frac),
        ytick = .data$y - y_span * 0.018,
        ylabel = .data$y + y_span * 0.018,
        value_col = value_col
      )
    if (length(facet_cols) > 0) {
      comp_df <- dplyr::bind_cols(facet_values[rep(1, nrow(comp_df)), , drop = FALSE], comp_df)
    }
    rows[[row_i]] <- comp_df
    row_i <- row_i + 1L
  }
  dplyr::bind_rows(rows)
}

add_significance_layers_04f <- function(plot_obj, sig_df, label_size = 2.2, linewidth = 0.22) {
  if (is.null(sig_df) || nrow(sig_df) == 0) return(plot_obj)
  sig_df <- sig_df %>%
    dplyr::filter(!is.na(.data$label), is.finite(.data$x), is.finite(.data$xend), is.finite(.data$y), is.finite(.data$ytick), is.finite(.data$ylabel))
  if (nrow(sig_df) == 0) return(plot_obj)
  plot_obj +
    geom_segment(data = sig_df, aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y), inherit.aes = FALSE, linewidth = linewidth) +
    geom_segment(data = sig_df, aes(x = .data$x, xend = .data$x, y = .data$ytick, yend = .data$y), inherit.aes = FALSE, linewidth = linewidth) +
    geom_segment(data = sig_df, aes(x = .data$xend, xend = .data$xend, y = .data$ytick, yend = .data$y), inherit.aes = FALSE, linewidth = linewidth) +
    geom_text(data = sig_df, aes(x = .data$xmid, y = .data$ylabel, label = .data$label), inherit.aes = FALSE, size = label_size, vjust = 0)
}

tidy_lm_coefficients_04f <- function(model, model_name) {
  if (is.null(model)) return(data.frame())
  co <- as.data.frame(summary(model)$coefficients)
  co$term <- rownames(co)
  rownames(co) <- NULL
  names(co) <- c("estimate", "std_error", "statistic", "p_value", "term")
  co %>%
    dplyr::mutate(model = model_name) %>%
    dplyr::select("model", "term", "estimate", "std_error", "statistic", "p_value")
}

tidy_anova_04f <- function(model, model_name) {
  if (is.null(model)) return(data.frame())
  av <- as.data.frame(stats::anova(model))
  av$term <- rownames(av)
  rownames(av) <- NULL
  names(av) <- gsub(" ", "_", tolower(names(av)), fixed = TRUE)
  names(av) <- gsub("\\(>f\\)", "p_value", names(av))
  av %>%
    dplyr::mutate(model = model_name) %>%
    dplyr::relocate("model", "term")
}

fit_lm_safe <- function(formula, data) {
  tryCatch(stats::lm(formula, data = data), error = function(e) NULL)
}

build_sample_info <- function(target_df, all_tumor_df) {
  selected_sample <- target_df %>%
    dplyr::group_by(.data$sampleID) %>%
    dplyr::summarise(
      Ploidy = first_present_value(.data$Ploidy),
      Dose = first_present_value(.data$Dose),
      DosePanel = first_present_value(.data$DosePanel),
      n_target_cells = dplyr::n(),
      mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
      median_cell_ploidy = stats::median(.data$cell_ploidy, na.rm = TRUE),
      mean_pseudotime = mean(.data$pseudotime, na.rm = TRUE),
      median_pseudotime = stats::median(.data$pseudotime, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      mean_cell_ploidy = ifelse(is.nan(.data$mean_cell_ploidy), NA_real_, .data$mean_cell_ploidy),
      median_cell_ploidy = ifelse(is.nan(.data$median_cell_ploidy), NA_real_, .data$median_cell_ploidy),
      mean_pseudotime = ifelse(is.nan(.data$mean_pseudotime), NA_real_, .data$mean_pseudotime),
      median_pseudotime = ifelse(is.nan(.data$median_pseudotime), NA_real_, .data$median_pseudotime)
    )

  tumor_sample <- all_tumor_df %>%
    dplyr::group_by(.data$sampleID) %>%
    dplyr::summarise(n_tumor_cells = dplyr::n(), .groups = "drop")

  selected_sample %>%
    dplyr::left_join(tumor_sample, by = "sampleID") %>%
    dplyr::mutate(
      target_cluster_fraction = .data$n_target_cells / .data$n_tumor_cells,
      target_cluster_percent = 100 * .data$target_cluster_fraction,
      Ploidy = factor(as.character(.data$Ploidy), levels = c("2N", "4N")),
      DosePanel = factor(as.character(.data$DosePanel), levels = c("0", "30", "120")),
      Dose = factor(as.character(.data$Dose), levels = c("0mg/kg", "30mg/kg", "120mg/kg"))
    )
}

add_sample_labels <- function(sample_info, sample_order) {
  sample_info %>%
    dplyr::mutate(
      sampleID = as.character(.data$sampleID),
      mean_cell_ploidy_label = ifelse(is.finite(.data$mean_cell_ploidy), sprintf("%.2f", .data$mean_cell_ploidy), "NA"),
      sampleFacetLabel = paste0(
        as.character(.data$Ploidy), " | ", as.character(.data$Dose), "\n",
        .data$sampleID, " | mean CP=", .data$mean_cell_ploidy_label
      ),
      sampleID = factor(.data$sampleID, levels = sample_order),
      sampleFacetLabel = factor(.data$sampleFacetLabel, levels = .data$sampleFacetLabel[order(.data$sampleID)])
    ) %>%
    dplyr::arrange(.data$sampleID) %>%
    dplyr::mutate(
      sampleID = as.character(.data$sampleID),
      sampleFacetLabel = as.character(.data$sampleFacetLabel)
    )
}

plot_initial_ploidy_cell_response <- function(df, file_stub) {
  plot_df <- add_dose_x_04f(df, "DosePanel")
  sig_df <- make_group_significance_04f(plot_df, "pseudotime", dose_col = "DosePanel", group_col = "Ploidy", group_levels = c("2N", "4N"))
  p <- ggplot(plot_df, aes(x = .data$dose_x, y = .data$pseudotime, fill = .data$Ploidy, group = interaction(.data$DosePanel, .data$Ploidy))) +
    geom_violin(position = position_dodge(width = 0.56), alpha = 0.25, trim = TRUE, linewidth = 0.2) +
    geom_boxplot(position = position_dodge(width = 0.56), width = 0.18, outlier.size = 0.25, alpha = 0.78) +
    dose_scale_04f(labels = c("0", "30", "120")) +
    scale_fill_manual(values = c("2N" = "#3B73B9", "4N" = "#D95F02"), drop = FALSE) +
    labs(
      title = "scVelo all_cells-derived non-cell-cycle Tumor clusters pseudotime by Initial Ploidy and Dose",
      x = "Dose (mg/kg)",
      y = "scVelo pseudotime",
      fill = "Initial Ploidy"
    ) +
    theme_04f(10) +
    theme(legend.position = "bottom")
  p <- add_significance_layers_04f(p, sig_df, label_size = 2.3)
  save_plot_both_04f(p, file_stub, width = 8.6, height = 5.5)
  invisible(p)
}

plot_sample_level_response <- function(sample_df, file_stub) {
  plot_df <- sample_df %>%
    dplyr::select("sampleID", "Ploidy", "DosePanel", "n_target_cells", "target_cluster_percent", "median_pseudotime", "mean_cell_ploidy") %>%
    tidyr::pivot_longer(
      cols = c("target_cluster_percent", "median_pseudotime", "mean_cell_ploidy"),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = factor(
        .data$metric,
        levels = c("target_cluster_percent", "median_pseudotime", "mean_cell_ploidy"),
        labels = c("Target cluster percent", "Sample median pseudotime", "Sample mean Cell Ploidy")
      )
    )
  plot_df <- add_dose_x_04f(plot_df, "DosePanel")
  sig_df <- make_group_significance_04f(plot_df, "value", dose_col = "DosePanel", group_col = "Ploidy", facet_cols = "metric", group_levels = c("2N", "4N"))
  p <- ggplot(plot_df, aes(x = .data$dose_x, y = .data$value, color = .data$Ploidy, group = interaction(.data$DosePanel, .data$Ploidy))) +
    geom_boxplot(position = position_dodge(width = 0.56), width = 0.5, outlier.shape = NA, alpha = 0.15) +
    geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.58), size = 2.2, alpha = 0.9) +
    facet_wrap(~metric, scales = "free_y", nrow = 1) +
    dose_scale_04f(labels = c("0", "30", "120")) +
    scale_color_manual(values = c("2N" = "#1D4E89", "4N" = "#A84600"), drop = FALSE) +
    labs(
      title = "Sample-level Dose response in non-cell-cycle Tumor clusters",
      subtitle = "Each point is one sample; target cluster percent uses all Tumor cells from the scVelo all_cells object as denominator.",
      x = "Dose (mg/kg)",
      y = NULL,
      color = "Initial Ploidy"
    ) +
    theme_04f(9) +
    theme(legend.position = "bottom")
  p <- add_significance_layers_04f(p, sig_df, label_size = 2.2)
  save_plot_both_04f(p, file_stub, width = 12.5, height = 4.8)
  invisible(p)
}

plot_high_low_response <- function(df, sample_high_df, file_stub_prefix) {
  cell_plot_df <- add_dose_x_04f(df, "DosePanel")
  cell_sig_df <- make_group_significance_04f(
    cell_plot_df,
    "pseudotime",
    dose_col = "DosePanel",
    group_col = "cell_ploidy_group",
    facet_cols = "Ploidy",
    group_levels = c("Low Cell Ploidy", "High Cell Ploidy"),
    group_label = "Cell Ploidy group"
  )
  p_cell <- ggplot(cell_plot_df, aes(x = .data$dose_x, y = .data$pseudotime, fill = .data$cell_ploidy_group, group = interaction(.data$DosePanel, .data$cell_ploidy_group))) +
    geom_violin(position = position_dodge(width = 0.56), alpha = 0.24, trim = TRUE, linewidth = 0.2) +
    geom_boxplot(position = position_dodge(width = 0.56), width = 0.18, outlier.size = 0.25, alpha = 0.75) +
    facet_wrap(~Ploidy, nrow = 1, drop = FALSE) +
    dose_scale_04f(labels = c("0", "30", "120")) +
    scale_fill_manual(values = c("Low Cell Ploidy" = "#4C78A8", "High Cell Ploidy" = "#E45756"), drop = FALSE) +
    labs(
      title = "Cell Ploidy high/low pseudotime response by Dose",
      subtitle = "High/low split uses the median Cell Ploidy across selected Tumor cells.",
      x = "Dose (mg/kg)",
      y = "scVelo pseudotime",
      fill = "Cell Ploidy group"
    ) +
    theme_04f(9) +
    theme(legend.position = "bottom")
  p_cell <- add_significance_layers_04f(p_cell, cell_sig_df, label_size = 2.1)
  save_plot_both_04f(p_cell, paste0(file_stub_prefix, "_cell_pseudotime"), width = 10, height = 5.2)

  sample_plot_df <- add_dose_x_04f(sample_high_df, "DosePanel")
  sample_sig_df <- make_group_significance_04f(sample_plot_df, "high_cell_ploidy_percent", dose_col = "DosePanel", group_col = "Ploidy", group_levels = c("2N", "4N"))
  p_sample <- ggplot(sample_plot_df, aes(x = .data$dose_x, y = .data$high_cell_ploidy_percent, color = .data$Ploidy, group = interaction(.data$DosePanel, .data$Ploidy))) +
    geom_boxplot(position = position_dodge(width = 0.56), width = 0.48, outlier.shape = NA, alpha = 0.1) +
    geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.55), size = 2.2, alpha = 0.9) +
    dose_scale_04f(labels = c("0", "30", "120")) +
    scale_color_manual(values = c("2N" = "#1D4E89", "4N" = "#A84600"), drop = FALSE) +
    labs(
      title = "Sample-level high Cell Ploidy fraction by Dose",
      x = "Dose (mg/kg)",
      y = "High Cell Ploidy cells (%)",
      color = "Initial Ploidy"
    ) +
    theme_04f(9) +
    theme(legend.position = "bottom")
  p_sample <- add_significance_layers_04f(p_sample, sample_sig_df, label_size = 2.3)
  save_plot_both_04f(p_sample, paste0(file_stub_prefix, "_sample_high_fraction"), width = 7.8, height = 5.2)
  invisible(list(cell = p_cell, sample = p_sample))
}

save_prompt_figure <- function(plot_obj, pdf_dir, png_dir, name, width = 8, height = 6) {
  save_plot_04f(plot_obj, file.path(pdf_dir, paste0(name, ".pdf")), width = width, height = height)
  save_plot_04f(plot_obj, file.path(png_dir, paste0(name, ".png")), width = width, height = height)
  invisible(name)
}

write_text_04f <- function(lines, file) {
  ensure_dir_04f(dirname(file))
  writeLines(as.character(lines), con = file)
  invisible(file)
}

logit_fraction_04f <- function(count, total) {
  stats::qlogis((count + 0.5) / (total + 1))
}

entropy_04f <- function(p) {
  p <- p[is.finite(p) & p > 0]
  if (length(p) == 0) return(NA_real_)
  -sum(p * log(p))
}

empirical_wasserstein_04f <- function(x, y, n_grid = 1000L) {
  x <- safe_num(x)
  y <- safe_num(y)
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  probs <- seq(0, 1, length.out = n_grid)
  mean(abs(stats::quantile(x, probs, names = FALSE, type = 8) - stats::quantile(y, probs, names = FALSE, type = 8)))
}

standard_model_table <- function(model, data, response_variable, analysis_name, model_formula, vcov_type = c("HC3", "CL"), cluster = NULL, notes = "") {
  vcov_type <- match.arg(vcov_type)
  if (is.null(model)) return(data.frame())
  vc <- tryCatch({
    if (identical(vcov_type, "CL") && !is.null(cluster)) {
      sandwich::vcovCL(model, cluster = cluster, type = "HC1")
    } else {
      sandwich::vcovHC(model, type = "HC3")
    }
  }, error = function(e) stats::vcov(model))
  ct <- tryCatch(lmtest::coeftest(model, vcov. = vc), error = function(e) NULL)
  if (is.null(ct)) return(data.frame())
  out <- data.frame(
    analysis_name = analysis_name,
    response_variable = response_variable,
    model_formula = model_formula,
    n_samples = if ("sampleID" %in% names(data)) {
      dplyr::n_distinct(data$sampleID)
    } else if ("sample_id" %in% names(data)) {
      dplyr::n_distinct(data$sample_id)
    } else {
      NA_integer_
    },
    n_cells = if ("n_cells" %in% names(data)) sum(data$n_cells, na.rm = TRUE) else nrow(data),
    coefficient = rownames(ct),
    estimate = as.numeric(ct[, 1]),
    std_error = as.numeric(ct[, 2]),
    statistic = as.numeric(ct[, 3]),
    p_value = as.numeric(ct[, 4]),
    notes = notes,
    stringsAsFactors = FALSE
  )
  out$conf_low <- out$estimate - 1.96 * out$std_error
  out$conf_high <- out$estimate + 1.96 * out$std_error
  out$effect_direction <- ifelse(out$estimate > 0, "positive", ifelse(out$estimate < 0, "negative", "zero"))
  out$q_value_BH <- stats::p.adjust(out$p_value, method = "BH")
  out %>%
    dplyr::select(
      "analysis_name", "response_variable", "model_formula", "n_samples", "n_cells",
      "coefficient", "estimate", "std_error", "conf_low", "conf_high", "p_value",
      "q_value_BH", "effect_direction", "notes"
    )
}

fit_robust_lm_table <- function(formula, data, response_variable, analysis_name, weights = NULL, cluster_col = NULL, notes = "") {
  model_data <- data
  if (!is.null(weights)) {
    model_data$.__weights_04f <- as.numeric(weights)
  }
  model <- tryCatch({
    if (is.null(weights)) stats::lm(formula, data = model_data) else stats::lm(formula, data = model_data, weights = .__weights_04f)
  }, error = function(e) NULL)
  cluster <- if (!is.null(cluster_col) && cluster_col %in% names(model_data)) model_data[[cluster_col]] else NULL
  vcov_type <- if (is.null(cluster)) "HC3" else "CL"
  list(
    model = model,
    table = standard_model_table(
      model,
      data = model_data,
      response_variable = response_variable,
      analysis_name = analysis_name,
      model_formula = paste(deparse(formula), collapse = " "),
      vcov_type = vcov_type,
      cluster = cluster,
      notes = notes
    )
  )
}

write_model_summary <- function(model, file) {
  if (is.null(model)) {
    write_text_04f("Model could not be fit.", file)
  } else {
    write_text_04f(capture.output(summary(model)), file)
  }
}

make_input_inventory <- function(input_root, selected_input_file, out_data_dir) {
  patterns <- "\\.(h5ad|loom|csv|tsv|parquet|rds|rda|pkl)$"
  files <- list.files(input_root, pattern = patterns, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  inventory <- data.frame(
    path = files,
    extension = tolower(tools::file_ext(files)),
    size_bytes = file.info(files)$size,
    modified_time = as.character(file.info(files)$mtime),
    selected = normalizePath(files, mustWork = FALSE) == normalizePath(selected_input_file, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
  write_csv_04f(inventory, file.path(out_data_dir, "input_file_inventory.tsv"))
  selected <- list(
    selected_input_file = selected_input_file,
    selected_reason = "Project-confirmed merged scVelo all_cells metadata with Cell Ploidy from 04e output; raw scVelo CSV lacks cell_ploidy.",
    raw_scvelo_input_root = input_root
  )
  writeLines(jsonlite::toJSON(selected, pretty = TRUE, auto_unbox = TRUE), file.path(out_data_dir, "selected_input_file.json"))
  inventory
}

clean_gene_symbol_04f <- function(x) {
  x <- as.character(x)
  x <- sub("^GRCh[0-9]+_", "", x, ignore.case = TRUE)
  x <- sub("^hg[0-9]+_", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(x)
}

signature_definitions_04f <- function() {
  list(
    gemcitabine_metabolism_transport = c("SLC29A1", "SLC29A2", "SLC28A1", "SLC28A3", "DCK", "CDA", "DCTD", "CMPK1", "NME1", "NME2", "NT5C2", "NT5C3A", "RRM1", "RRM2"),
    s_phase_dna_replication = c("MKI67", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "TYMS", "TK1", "DHFR", "CDC6", "ORC1", "RPA1", "RPA2", "CLSPN"),
    replication_stress_atr_chk1 = c("ATR", "ATRIP", "CHEK1", "CHEK2", "WEE1", "CLSPN", "RPA1", "RPA2", "TIMELESS", "TIPIN"),
    dna_damage_response = c("H2AFX", "GADD45A", "GADD45B", "GADD45G", "RAD51", "BRCA1", "BRCA2", "PARP1", "XRCC5", "XRCC6", "TP53BP1", "MDC1"),
    g2m_checkpoint = c("CDK1", "CCNB1", "CCNB2", "CDC20", "CDC25C", "PLK1", "AURKA", "AURKB", "TOP2A", "UBE2C", "BIRC5"),
    apoptosis = c("BAX", "BAK1", "BCL2", "BCL2L1", "BBC3", "PMAIP1", "CASP3", "CASP7", "CASP8", "CASP9"),
    senescence_p53_response = c("CDKN1A", "CDKN2A", "MDM2", "SERPINE1", "LMNB1", "TP53I3", "BTG2", "DDB2")
  )
}

read_h5ad_signature_expression <- function(h5ad_file, cell_ids, signatures) {
  if (!file.exists(h5ad_file)) {
    return(list(status = "missing_h5ad", expr = NULL, matching = data.frame(), source_info = data.frame()))
  }
  h <- hdf5r::H5File$new(h5ad_file, mode = "r")
  on.exit(h$close_all(), add = TRUE)
  if (!(all(c("X", "obs", "var") %in% names(h)))) {
    return(list(status = "missing_required_h5ad_groups", expr = NULL, matching = data.frame(), source_info = data.frame()))
  }
  obs_cells <- h[["obs/cell"]]$read()
  var_genes <- if ("Gene" %in% names(h[["var"]])) h[["var/Gene"]]$read() else h[["var/_index"]]$read()
  var_clean <- clean_gene_symbol_04f(var_genes)
  query <- unique(unlist(signatures, use.names = FALSE))
  query_clean <- clean_gene_symbol_04f(query)
  match_idx <- match(query_clean, var_clean)
  matching <- data.frame(
    query_gene = query,
    query_gene_clean = query_clean,
    matched = !is.na(match_idx),
    h5ad_gene = ifelse(is.na(match_idx), NA_character_, var_genes[match_idx]),
    h5ad_gene_index_1based = match_idx,
    stringsAsFactors = FALSE
  )
  matched_idx <- unique(match_idx[!is.na(match_idx)])
  if (length(matched_idx) == 0) {
    return(list(status = "no_signature_genes_matched", expr = NULL, matching = matching, source_info = data.frame()))
  }
  row_idx <- match(cell_ids, obs_cells)
  keep <- !is.na(row_idx)
  if (!all(keep)) {
    cell_ids <- cell_ids[keep]
    row_idx <- row_idx[keep]
  }
  indptr <- h[["X/indptr"]]$read()
  shape <- h[["X"]]$attr_open("shape")$read()
  wanted0 <- matched_idx - 1L
  expr <- matrix(0, nrow = length(row_idx), ncol = length(matched_idx))
  colnames(expr) <- var_genes[matched_idx]
  rownames(expr) <- cell_ids
  names_wanted <- as.character(wanted0)
  col_map <- seq_along(matched_idx)
  names(col_map) <- names_wanted
  for (i in seq_along(row_idx)) {
    r <- row_idx[i]
    start <- indptr[r] + 1L
    end <- indptr[r + 1L]
    if (!is.finite(start) || !is.finite(end) || end < start) next
    idx0 <- h[["X/indices"]]$read(args = list(seq.int(start, end)))
    hit <- idx0 %in% wanted0
    if (!any(hit)) next
    vals <- h[["X/data"]]$read(args = list(seq.int(start, end)))
    expr[i, col_map[as.character(idx0[hit])]] <- vals[hit]
  }
  sample_vals <- as.vector(expr[seq_len(min(nrow(expr), 50)), seq_len(min(ncol(expr), 20)), drop = FALSE])
  sample_vals <- sample_vals[is.finite(sample_vals) & sample_vals != 0]
  integer_like <- length(sample_vals) > 0 && mean(abs(sample_vals - round(sample_vals)) < 1e-6) > 0.98
  source_info <- data.frame(
    h5ad_file = h5ad_file,
    expression_source = "X",
    n_obs = shape[1],
    n_vars = shape[2],
    n_target_cells_matched = nrow(expr),
    n_signature_genes_matched = ncol(expr),
    expression_integer_like = integer_like,
    notes = ifelse(integer_like, "X appears integer-like.", "X does not appear integer-like; DESeq2 is skipped and expression analyses are exploratory."),
    stringsAsFactors = FALSE
  )
  list(status = "ok", expr = expr, matching = matching, source_info = source_info)
}

score_signatures_04f <- function(expr, signatures, matching) {
  if (is.null(expr) || nrow(expr) == 0 || ncol(expr) == 0) return(data.frame())
  expr_clean <- clean_gene_symbol_04f(colnames(expr))
  z <- scale(expr)
  z[, !is.finite(colSums(z, na.rm = TRUE))] <- NA_real_
  out <- data.frame(cell_id = rownames(expr), stringsAsFactors = FALSE)
  for (sig in names(signatures)) {
    genes <- clean_gene_symbol_04f(signatures[[sig]])
    use <- which(expr_clean %in% genes)
    if (length(use) >= 3) {
      out[[sig]] <- rowMeans(z[, use, drop = FALSE], na.rm = TRUE)
    } else {
      out[[sig]] <- NA_real_
    }
  }
  out
}

read_h5ad_group_expression_summaries <- function(h5ad_file, meta_df, group_cols) {
  if (!file.exists(h5ad_file)) {
    return(list(status = "missing_h5ad", group_sums = NULL, group_meta = data.frame(), gene_names = character(0), source_info = data.frame()))
  }
  h <- hdf5r::H5File$new(h5ad_file, mode = "r")
  on.exit(h$close_all(), add = TRUE)
  if (!(all(c("X", "obs", "var") %in% names(h)))) {
    return(list(status = "missing_required_h5ad_groups", group_sums = NULL, group_meta = data.frame(), gene_names = character(0), source_info = data.frame()))
  }
  if (!all(group_cols %in% names(meta_df)) || !"cell_id" %in% names(meta_df)) {
    return(list(status = "missing_required_metadata_columns", group_sums = NULL, group_meta = data.frame(), gene_names = character(0), source_info = data.frame()))
  }
  obs_cells <- h[["obs/cell"]]$read()
  var_genes <- if ("Gene" %in% names(h[["var"]])) h[["var/Gene"]]$read() else h[["var/_index"]]$read()
  shape <- h[["X"]]$attr_open("shape")$read()

  meta_use <- meta_df %>%
    dplyr::mutate(cell_id = as.character(.data$cell_id)) %>%
    dplyr::filter(!is.na(.data$cell_id), nzchar(.data$cell_id))
  row_idx <- match(meta_use$cell_id, obs_cells)
  keep <- !is.na(row_idx)
  meta_use <- meta_use[keep, , drop = FALSE]
  row_idx <- row_idx[keep]
  if (nrow(meta_use) == 0) {
    return(list(status = "no_target_cells_matched_h5ad", group_sums = NULL, group_meta = data.frame(), gene_names = var_genes, source_info = data.frame()))
  }

  group_df <- meta_use %>%
    dplyr::select(dplyr::all_of(group_cols)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  group_id <- do.call(paste, c(group_df, sep = "|"))
  group_levels <- unique(group_id)
  group_index <- match(group_id, group_levels)
  group_meta <- group_df[match(group_levels, group_id), , drop = FALSE]
  group_meta$group_id <- group_levels
  group_meta$n_cells <- as.integer(tabulate(group_index, nbins = length(group_levels)))
  group_meta <- group_meta %>% dplyr::relocate("group_id", "n_cells")

  n_groups <- length(group_levels)
  n_vars <- as.integer(shape[2])
  group_sums <- matrix(0, nrow = n_groups, ncol = n_vars)
  colnames(group_sums) <- var_genes
  rownames(group_sums) <- group_levels
  total_sums <- numeric(n_vars)
  indptr <- h[["X/indptr"]]$read()
  for (i in seq_along(row_idx)) {
    r <- row_idx[i]
    start <- indptr[r] + 1L
    end <- indptr[r + 1L]
    if (!is.finite(start) || !is.finite(end) || end < start) next
    idx0 <- h[["X/indices"]]$read(args = list(seq.int(start, end)))
    vals <- h[["X/data"]]$read(args = list(seq.int(start, end)))
    idx <- idx0 + 1L
    g <- group_index[i]
    group_sums[g, idx] <- group_sums[g, idx] + vals
    total_sums[idx] <- total_sums[idx] + vals
  }
  sample_vals <- as.vector(group_sums[seq_len(min(nrow(group_sums), 10)), seq_len(min(ncol(group_sums), 200)), drop = FALSE])
  sample_vals <- sample_vals[is.finite(sample_vals) & sample_vals != 0]
  integer_like <- length(sample_vals) > 0 && mean(abs(sample_vals - round(sample_vals)) < 1e-6) > 0.98
  source_info <- data.frame(
    h5ad_file = h5ad_file,
    expression_source = "X",
    n_obs = shape[1],
    n_vars = shape[2],
    n_target_cells_matched = length(row_idx),
    n_groups = n_groups,
    expression_integer_like = integer_like,
    notes = "All h5ad X genes were aggregated by target Tumor Dose x Initial Ploidy x pseudotime-bin groups for preranked GSEA.",
    stringsAsFactors = FALSE
  )
  list(status = "ok", group_sums = group_sums, total_sums = total_sums, group_meta = group_meta, gene_names = var_genes, source_info = source_info)
}

collapse_rank_stats_by_gene_04f <- function(stats, genes) {
  genes <- clean_gene_symbol_04f(genes)
  stats <- safe_num(stats)
  keep <- is.finite(stats) & !is.na(genes) & nzchar(genes)
  stats <- stats[keep]
  genes <- genes[keep]
  if (length(stats) == 0) return(numeric(0))
  collapsed <- stats::aggregate(stats, by = list(gene = genes), FUN = mean, na.rm = TRUE)
  out <- collapsed$x
  names(out) <- collapsed$gene
  out[order(out, decreasing = TRUE)]
}

gsea_es_from_positions_04f <- function(sorted_stats, positions, exponent = 1) {
  n <- length(sorted_stats)
  positions <- sort(unique(as.integer(positions[is.finite(positions)])))
  positions <- positions[positions >= 1L & positions <= n]
  k <- length(positions)
  if (n == 0 || k == 0 || k >= n) return(NA_real_)
  weights <- abs(sorted_stats[positions])^exponent
  if (!all(is.finite(weights)) || sum(weights) <= 0) weights <- rep(1, k)
  hit_cum <- cumsum(weights / sum(weights))
  miss_before <- (positions - seq_len(k)) / (n - k)
  running_before <- c(0, head(hit_cum, -1)) - miss_before
  running_after <- hit_cum - miss_before
  max_es <- max(running_after, na.rm = TRUE)
  min_es <- min(running_before, na.rm = TRUE)
  if (abs(max_es) >= abs(min_es)) max_es else min_es
}

compute_preranked_gsea_04f <- function(gsea_input, signatures, nperm = 1000L, min_size = 3L, exponent = 1, seed = 1234L) {
  if (!identical(gsea_input$status, "ok") || is.null(gsea_input$group_sums) || nrow(gsea_input$group_sums) == 0) return(data.frame())
  set.seed(seed)
  group_sums <- gsea_input$group_sums
  total_sums <- gsea_input$total_sums
  group_meta <- gsea_input$group_meta
  gene_names <- colnames(group_sums)
  total_n <- sum(group_meta$n_cells)
  nperm <- as.integer(nperm)
  if (!is.finite(nperm) || nperm < 1L) nperm <- 1000L
  min_size <- as.integer(min_size)
  if (!is.finite(min_size) || min_size < 1L) min_size <- 3L

  rows <- list()
  row_i <- 1L
  for (g in seq_len(nrow(group_sums))) {
    n_group <- group_meta$n_cells[g]
    n_other <- total_n - n_group
    if (!is.finite(n_group) || !is.finite(n_other) || n_group < 1 || n_other < 1) next
    mean_group <- group_sums[g, ] / n_group
    mean_other <- (total_sums - group_sums[g, ]) / n_other
    rank_stats <- collapse_rank_stats_by_gene_04f(mean_group - mean_other, gene_names)
    rank_stats <- rank_stats[is.finite(rank_stats)]
    rank_stats <- rank_stats[order(rank_stats, decreasing = TRUE)]
    if (length(rank_stats) < 10) next
    stat_names <- names(rank_stats)
    for (sig in names(signatures)) {
      genes <- unique(clean_gene_symbol_04f(signatures[[sig]]))
      positions <- match(genes, stat_names)
      positions <- positions[!is.na(positions)]
      set_size <- length(unique(positions))
      if (set_size < min_size || set_size >= length(rank_stats)) {
        es <- NA_real_
        nes <- NA_real_
        p_value <- NA_real_
      } else {
        es <- gsea_es_from_positions_04f(rank_stats, positions, exponent = exponent)
        null_es <- replicate(nperm, gsea_es_from_positions_04f(rank_stats, sample.int(length(rank_stats), set_size), exponent = exponent))
        null_es <- null_es[is.finite(null_es)]
        if (!is.finite(es) || length(null_es) == 0) {
          nes <- NA_real_
          p_value <- NA_real_
        } else if (es >= 0) {
          denom <- mean(null_es[null_es >= 0], na.rm = TRUE)
          if (!is.finite(denom) || denom == 0) denom <- mean(abs(null_es), na.rm = TRUE)
          nes <- es / denom
          p_value <- mean(null_es >= es)
        } else {
          denom <- abs(mean(null_es[null_es < 0], na.rm = TRUE))
          if (!is.finite(denom) || denom == 0) denom <- mean(abs(null_es), na.rm = TRUE)
          nes <- es / denom
          p_value <- mean(null_es <= es)
        }
      }
      rows[[row_i]] <- data.frame(
        group_meta[g, , drop = FALSE],
        signature = sig,
        set_size = set_size,
        ES = es,
        NES = nes,
        p_value = p_value,
        nperm = nperm,
        rank_statistic = "mean_X_in_group_minus_mean_X_in_other_target_cells",
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }
  out <- dplyr::bind_rows(rows)
  if (nrow(out) > 0) {
    out <- out %>%
      dplyr::group_by(.data$signature) %>%
      dplyr::mutate(p_adjust_BH = stats::p.adjust(.data$p_value, method = "BH")) %>%
      dplyr::ungroup()
  }
  out
}

scenario_id <- get_env_scalar("SCENARIO", "ROOT_6_END_NULL")
results_root <- get_env_scalar("EXTRA_RESULTS_ROOT", "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results")
non_cell_cycle_results_root <- get_env_scalar(
  "NON_CELL_CYCLE_RESULTS_SHARED_ROOT",
  file.path(results_root, "04g_nonCellCycle_results")
)
growth_curve_file <- get_env_scalar(
  "IN_VIVO_GROWTH_CURVE_FILE",
  file.path(project_root, "Data", "in-vivo", "dt_Gem_VT_20241223_v4.xlsx")
)
input_file <- get_env_scalar(
  "SCVELO_ALL_CELLS_METRICS",
  file.path(results_root, "04e_velocity_results_all_cells", scenario_id, "scVelo", "All_cells", "cell_metrics_with_cell_ploidy.csv")
)
output_base <- get_env_scalar(
  "NON_CELL_CYCLE_RESULTS_ROOT",
  file.path(results_root, "04g_nonCellCycle_results", scenario_id)
)
excluded_clusters <- strsplit(get_env_scalar("EXCLUDE_CLUSTERS", "4c,6,10"), ",", fixed = TRUE)[[1]]
excluded_clusters <- trimws(excluded_clusters)

if (!file.exists(input_file)) {
  stop("Missing scVelo all_cells input file: ", input_file, call. = FALSE)
}

invisible(ensure_dir_04f(output_base))
out_inputs <- ensure_dir_04f(file.path(output_base, "00_inputs"))
out_hist <- ensure_dir_04f(file.path(output_base, "01_sample_pseudotime_histograms"))
out_density <- ensure_dir_04f(file.path(output_base, "02_sample_pseudotime_density"))
out_dose_hist <- ensure_dir_04f(file.path(output_base, "03_dose_pseudotime_histograms"))
out_dose_density <- ensure_dir_04f(file.path(output_base, "04_dose_pseudotime_density"))
out_response <- ensure_dir_04f(file.path(output_base, "05_dose_response"))

raw_scvelo_input_root <- get_env_scalar(
  "RAW_SCVELO_ALL_CELLS_ROOT",
  file.path(results_root, "04_trajectory", scenario_id, "scvelo", "01_velocity_groups", "All_cells")
)
h5ad_file <- get_env_scalar(
  "SCVELO_ALL_CELLS_H5AD",
  file.path(raw_scvelo_input_root, "scvelo_result.h5ad")
)
prompt_data_dir <- ensure_dir_04f(file.path(output_base, "data"))
prompt_fig_dir <- ensure_dir_04f(file.path(output_base, "figures"))
prompt_pdf_dir <- ensure_dir_04f(file.path(prompt_fig_dir, "pdf"))
prompt_png_dir <- ensure_dir_04f(file.path(prompt_fig_dir, "png"))
prompt_stats_dir <- ensure_dir_04f(file.path(output_base, "stats"))
prompt_models_dir <- ensure_dir_04f(file.path(output_base, "models"))
prompt_logs_dir <- ensure_dir_04f(file.path(output_base, "logs"))
prompt_report_dir <- ensure_dir_04f(file.path(output_base, "report"))
prompt_scripts_dir <- ensure_dir_04f(file.path(output_base, "scripts"))
non_cell_cycle_results_root <- ensure_dir_04f(non_cell_cycle_results_root)
fig4_tgi_export_file <- file.path(
  non_cell_cycle_results_root,
  "fig4_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
)
fig4_tgi_legacy_export_file <- file.path(
  non_cell_cycle_results_root,
  "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
)
tgi_process_file <- file.path(non_cell_cycle_results_root, "TGI_calculation_process.md")

message("04g non-cell-cycle results")
message("  scenario: ", scenario_id)
message("  input all_cells scVelo metrics: ", input_file)
message("  output_base: ", output_base)
message("  growth curve file for TGI: ", growth_curve_file)
message("  fig4 TGI export: ", fig4_tgi_export_file)
message("  fig4 TGI legacy export: ", fig4_tgi_legacy_export_file)
message("  excluded clusters: ", paste(excluded_clusters, collapse = ", "))

raw_df <- readr::read_csv(
  input_file,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE,
  progress = FALSE
) %>%
  as.data.frame(stringsAsFactors = FALSE)

required_cols <- c("cell", "TN", "clusters", "sample", "Ploidy", "Dose", "pseudotime", "cell_ploidy")
missing_cols <- setdiff(required_cols, names(raw_df))
if (length(missing_cols) > 0) {
  stop("Input all_cells metrics file is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

all_cells <- raw_df %>%
  dplyr::mutate(
    cell = as_clean_chr_04f(.data$cell),
    TN = as_clean_chr_04f(.data$TN),
    clusters = as_clean_chr_04f(.data$clusters),
    sampleID = as_clean_chr_04f(.data$sample),
    Ploidy = standardize_ploidy_04f(.data$Ploidy),
    Dose = standardize_dose_04f(.data$Dose),
    DosePanel = dose_panel_label(.data$Dose),
    pseudotime = safe_num(.data$pseudotime),
    cell_ploidy = safe_num(.data$cell_ploidy),
    dose_numeric = safe_num(.data$DosePanel)
  ) %>%
  dplyr::filter(!is.na(.data$cell), nzchar(.data$cell), is.finite(.data$pseudotime))

tumor_df <- all_cells %>%
  dplyr::filter(.data$TN == "Tumor", .data$Ploidy %in% c("2N", "4N"), .data$DosePanel %in% c("0", "30", "120"))

target_clusters_env <- get_env_scalar("TARGET_CLUSTERS", "")
if (nzchar(target_clusters_env)) {
  target_clusters <- strsplit(target_clusters_env, ",", fixed = TRUE)[[1]]
  target_clusters <- trimws(target_clusters)
} else {
  tumor_clusters <- unique(as.character(tumor_df$clusters))
  tumor_clusters <- tumor_clusters[!is.na(tumor_clusters) & nzchar(tumor_clusters)]
  target_clusters <- setdiff(tumor_clusters, excluded_clusters)
}
target_clusters <- target_clusters[!is.na(target_clusters) & nzchar(target_clusters)]
target_clusters <- target_clusters[target_clusters %in% unique(as.character(tumor_df$clusters))]
if (length(target_clusters) == 0) {
  stop("No target Tumor clusters remain after excluding clusters: ", paste(excluded_clusters, collapse = ", "), call. = FALSE)
}
target_clusters <- target_clusters[order(match(target_clusters, unique(as.character(tumor_df$clusters))))]
target_cluster_label <- paste(target_clusters, collapse = ", ")
excluded_cluster_label <- paste(excluded_clusters, collapse = ", ")
analysis_scope_label <- paste0("Tumor cells from scVelo all_cells; excluding clusters ", excluded_cluster_label)
cluster_palette <- make_cluster_palette_04g(target_clusters)
cluster_facet_nrow <- cluster_facet_nrow_04g(length(target_clusters))

message("  target clusters: ", target_cluster_label)

target_df <- tumor_df %>%
  dplyr::filter(.data$clusters %in% target_clusters) %>%
  dplyr::filter(!is.na(.data$sampleID), nzchar(.data$sampleID), is.finite(.data$pseudotime)) %>%
  dplyr::mutate(
    Ploidy = factor(as.character(.data$Ploidy), levels = c("2N", "4N")),
    Dose = factor(as.character(.data$Dose), levels = c("0mg/kg", "30mg/kg", "120mg/kg")),
    DosePanel = factor(as.character(.data$DosePanel), levels = c("0", "30", "120")),
    clusters = factor(as.character(.data$clusters), levels = target_clusters)
  )

if (nrow(target_df) == 0) {
  stop("No Tumor cells in non-cell-cycle target clusters were found after filtering all_cells scVelo metrics.", call. = FALSE)
}

sample_info_base <- build_sample_info(target_df, tumor_df)
sample_tgi <- calculate_in_vivo_tgi(
  growth_curve_file,
  control_dose = "0mg/kg",
  match_control_by = "initial_ploidy"
)
tgi_extra_cols <- grep(
  "^(tumor_volume_Day_|tumor_volume_delta_Day_|matched_control_mean_delta_Day_|matched_control_n_Day_|TGI_percent_Day_|tumor_volume_auc|tumor_volume_auc_delta|matched_control_mean_auc_delta|matched_control_n_auc|TGI_percent_auc|auc_)",
  names(sample_tgi),
  value = TRUE
)
tgi_endpoint_days <- sub("^TGI_percent_", "", grep("^TGI_percent_Day_", names(sample_tgi), value = TRUE))
tgi_endpoint_days <- tgi_endpoint_days[order(suppressWarnings(as.numeric(sub("^Day_", "", tgi_endpoint_days))))]
sample_tgi_join <- sample_tgi %>%
  dplyr::transmute(
    sampleID = .data$sample_id,
    growth_curve_sample_id_raw = .data$sample_id_raw,
    growth_curve_harvest = .data$harvest,
    tgi_initial_ploidy = .data$initial_ploidy,
    tgi_dose = .data$dose,
    tumor_volume_baseline_day = .data$baseline_day,
    tumor_volume_baseline = .data$baseline_volume,
    tumor_volume_final_day = .data$final_day,
    tumor_volume_final = .data$final_volume,
    tumor_volume_delta = .data$tumor_volume_delta,
    matched_control_mean_delta = .data$matched_control_mean_delta,
    matched_control_n = .data$matched_control_n,
    TGI_percent = .data$TGI_percent
  )
if (length(tgi_extra_cols) > 0) {
  sample_tgi_join <- sample_tgi_join %>%
    dplyr::left_join(
      sample_tgi %>%
        dplyr::select("sample_id", dplyr::all_of(tgi_extra_cols)) %>%
        dplyr::rename(sampleID = .data$sample_id),
      by = "sampleID"
    )
}
tgi_missing_samples <- setdiff(unique(as.character(sample_info_base$sampleID)), sample_tgi_join$sampleID)
if (length(tgi_missing_samples) > 0) {
  warning("TGI was not matched for target samples: ", paste(tgi_missing_samples, collapse = ", "), call. = FALSE)
}
reference_sample_order <- sample_info_base %>%
  dplyr::mutate(
    .ploidy_order = match(as.character(.data$Ploidy), c("2N", "4N")),
    .dose_order = match(as.character(.data$DosePanel), c("0", "30", "120"))
  ) %>%
  dplyr::arrange(.data$.ploidy_order, .data$.dose_order, .data$sampleID) %>%
  dplyr::pull(.data$sampleID) %>%
  as.character()
mean_cell_ploidy_sample_order <- sample_info_base %>%
  dplyr::arrange(.data$mean_cell_ploidy, .data$sampleID) %>%
  dplyr::pull(.data$sampleID) %>%
  as.character()

sample_info_reference <- add_sample_labels(sample_info_base, reference_sample_order)
sample_info_ploidy_order <- add_sample_labels(sample_info_base, mean_cell_ploidy_sample_order)

analysis_manifest <- data.frame(
  scenario = scenario_id,
  input_file = input_file,
  output_base = output_base,
  growth_curve_file = growth_curve_file,
  fig4_tgi_export_file = fig4_tgi_export_file,
  fig4_tgi_legacy_export_file = fig4_tgi_legacy_export_file,
  source_object = "scVelo All_cells",
  extracted_tn = "Tumor",
  excluded_clusters = excluded_cluster_label,
  target_clusters = target_cluster_label,
  n_all_cells_finite_pseudotime = nrow(all_cells),
  n_tumor_cells_finite_pseudotime = nrow(tumor_df),
  n_target_tumor_cells = nrow(target_df),
  n_target_samples = length(unique(target_df$sampleID)),
  stringsAsFactors = FALSE
)
write_csv_04f(analysis_manifest, file.path(out_inputs, "analysis_manifest.csv"))

filter_summary <- dplyr::bind_rows(
  all_cells %>% dplyr::summarise(filter = "all_cells_scVelo", n_cells = dplyr::n(), n_samples = dplyr::n_distinct(.data$sampleID)),
  tumor_df %>% dplyr::summarise(filter = "all_cells_scVelo_Tumor", n_cells = dplyr::n(), n_samples = dplyr::n_distinct(.data$sampleID)),
  target_df %>% dplyr::summarise(filter = "all_cells_scVelo_Tumor_excluding_clusters_4c_6_10", n_cells = dplyr::n(), n_samples = dplyr::n_distinct(.data$sampleID))
)
write_csv_04f(filter_summary, file.path(out_inputs, "filter_summary.csv"))

cluster_sample_summary <- target_df %>%
  dplyr::group_by(.data$clusters, .data$Ploidy, .data$DosePanel, .data$sampleID) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean_pseudotime = mean(.data$pseudotime, na.rm = TRUE),
    median_pseudotime = stats::median(.data$pseudotime, na.rm = TRUE),
    mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_04f(cluster_sample_summary, file.path(out_inputs, "target_cluster_sample_summary.csv"))
write_csv_04f(sample_info_base, file.path(out_inputs, "target_sample_summary.csv"))
write_csv_04f(sample_tgi, file.path(out_inputs, "sample_tgi_from_growth_curve.csv"))
write_csv_04f(sample_tgi, file.path(prompt_data_dir, "sample_tgi_from_growth_curve.tsv"))
write_csv_04f(
  data.frame(
    n_target_samples = length(unique(as.character(sample_info_base$sampleID))),
    n_tgi_rows = nrow(sample_tgi),
    n_target_samples_without_tgi = length(tgi_missing_samples),
    target_samples_without_tgi = paste(tgi_missing_samples, collapse = ";"),
    tgi_formula = "100 * (1 - tumor_volume_delta / matched_control_mean_delta)",
    endpoint_tgi_days = paste(tgi_endpoint_days, collapse = ";"),
    endpoint_tgi_formula = "100 * (1 - tumor_volume_delta_Day_X / matched_control_mean_delta_Day_X)",
    auc_tgi_formula = "100 * (1 - tumor_volume_auc_delta / matched_control_mean_auc_delta)",
    auc_method = "Trapezoidal AUC over all finite Day_* tumor-volume points after subtracting baseline volume.",
    control_matching = "initial_ploidy-specific controls: 2N vs 2N 0mg/kg, 4N vs 4N 0mg/kg",
    stringsAsFactors = FALSE
  ),
  file.path(prompt_logs_dir, "tgi_export_status.tsv")
)
write_csv_04f(sample_info_reference, file.path(out_inputs, "sample_order_reference.csv"))
write_csv_04f(sample_info_ploidy_order, file.path(out_inputs, "sample_order_by_mean_cell_ploidy.csv"))

hist_reference <- make_sample_histogram_data(target_df, sample_info_reference, reference_sample_order, n_bins = 200L)
hist_ploidy_order <- make_sample_histogram_data(target_df, sample_info_ploidy_order, mean_cell_ploidy_sample_order, n_bins = 200L)
write_csv_04f(hist_reference, file.path(out_hist, "all_samples_pseudotime_frequency_histograms_16x1_data.csv"))
write_csv_04f(hist_ploidy_order, file.path(out_hist, "all_samples_pseudotime_frequency_histograms_16x1_by_mean_cell_ploidy_data.csv"))

fig4_tgi_export <- target_df %>%
  dplyr::transmute(
    cell_id = .data$cell,
    sampleID = as.character(.data$sampleID),
    cluster = as.character(.data$clusters),
    initial_ploidy = as.character(.data$Ploidy),
    gemcitabine_dose = as.character(.data$Dose),
    gemcitabine_dose_mg_per_kg = safe_num(as.character(.data$DosePanel)),
    pseudotime = .data$pseudotime,
    cell_ploidy = .data$cell_ploidy
  ) %>%
  dplyr::left_join(
    sample_info_base %>%
      dplyr::transmute(
        sampleID = as.character(.data$sampleID),
        average_ploidy = .data$mean_cell_ploidy,
        median_ploidy = .data$median_cell_ploidy,
        sample_mean_pseudotime = .data$mean_pseudotime,
        sample_median_pseudotime = .data$median_pseudotime,
        n_target_cells = .data$n_target_cells,
        n_tumor_cells = .data$n_tumor_cells,
        target_cluster_fraction = .data$target_cluster_fraction,
        target_cluster_percent = .data$target_cluster_percent
      ),
    by = "sampleID"
  ) %>%
  dplyr::left_join(sample_tgi_join, by = "sampleID") %>%
  dplyr::rename(sample_id = sampleID) %>%
  dplyr::arrange(.data$initial_ploidy, .data$gemcitabine_dose_mg_per_kg, .data$sample_id, .data$pseudotime, .data$cell_id)
write_csv_04f(fig4_tgi_export, fig4_tgi_export_file)
write_csv_04f(fig4_tgi_export, fig4_tgi_legacy_export_file)
write_csv_04f(fig4_tgi_export, file.path(prompt_data_dir, "fig4_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.tsv"))

write_text_04f(
  make_in_vivo_tgi_process_markdown(
    growth_curve_file = growth_curve_file,
    output_csv_file = fig4_tgi_export_file,
    analysis_label = "04g non-cell-cycle results",
    target_cell_description = "target Tumor cell from all clusters except 4c, 6, and 10",
    cluster_description = paste0("04g non-cell-cycle Tumor clusters: ", target_cluster_label, "; excluded clusters: ", excluded_cluster_label),
    code_location = "/Users/4482173/Documents/GitHub/Gemcitabine-model/Code/in-vivo/04g_nonCellCycle_results.R",
    sample_tgi = sample_tgi
  ),
  tgi_process_file
)

plot_sample_histogram_16x1(
  hist_reference,
  "scVelo all_cells-derived non-cell-cycle Tumor cluster sample pseudotime histograms",
  file.path(out_hist, "all_samples_pseudotime_frequency_histograms_16x1")
)
plot_sample_histogram_16x1(
  hist_ploidy_order,
  "scVelo all_cells-derived non-cell-cycle Tumor cluster sample pseudotime histograms ordered by mean Cell Ploidy",
  file.path(out_hist, "all_samples_pseudotime_frequency_histograms_16x1_by_mean_cell_ploidy")
)

density_df <- make_density_data(target_df, sample_info_reference, n = 512L)
write_csv_04f(density_df, file.path(out_density, "all_samples_pseudotime_density_data.csv"))
write_csv_04f(
  sample_info_reference %>% dplyr::select("sampleID", "Ploidy", "Dose", "DosePanel", "mean_cell_ploidy", "sampleFacetLabel"),
  file.path(out_density, "sample_density_color_mapping.csv")
)
plot_density_curves(
  density_df,
  "scVelo all_cells-derived non-cell-cycle Tumor cluster sample pseudotime density",
  file.path(out_density, "all_samples_pseudotime_density_by_mean_cell_ploidy"),
  facet_by_dose = FALSE
)

plot_dose_histogram_panels(
  hist_reference,
  "Dose-stratified sample pseudotime histograms",
  file.path(out_dose_hist, "all_samples_pseudotime_frequency_histograms_by_dose_16x1"),
  sample_info_reference,
  reference_sample_order
)
plot_dose_histogram_panels(
  hist_ploidy_order,
  "Dose-stratified sample pseudotime histograms ordered by mean Cell Ploidy",
  file.path(out_dose_hist, "all_samples_pseudotime_frequency_histograms_by_dose_16x1_by_mean_cell_ploidy"),
  sample_info_ploidy_order,
  mean_cell_ploidy_sample_order
)
write_csv_04f(hist_reference, file.path(out_dose_hist, "all_samples_pseudotime_frequency_histograms_by_dose_16x1_data.csv"))

plot_density_curves(
  density_df,
  "Dose-stratified sample pseudotime density",
  file.path(out_dose_density, "all_samples_pseudotime_density_by_dose_mean_cell_ploidy"),
  facet_by_dose = TRUE
)
write_csv_04f(density_df, file.path(out_dose_density, "all_samples_pseudotime_density_by_dose_data.csv"))

cell_summary <- target_df %>%
  dplyr::group_by(.data$Ploidy, .data$DosePanel) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    n_samples = dplyr::n_distinct(.data$sampleID),
    mean_pseudotime = mean(.data$pseudotime, na.rm = TRUE),
    median_pseudotime = stats::median(.data$pseudotime, na.rm = TRUE),
    mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
    median_cell_ploidy = stats::median(.data$cell_ploidy, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_04f(cell_summary, file.path(out_response, "initial_ploidy_dose_cell_level_summary.csv"))

cell_kruskal <- target_df %>%
  dplyr::group_by(.data$Ploidy) %>%
  dplyr::group_modify(~ safe_kruskal(.x, "pseudotime", "DosePanel")) %>%
  dplyr::ungroup()
write_csv_04f(cell_kruskal, file.path(out_response, "initial_ploidy_dose_cell_level_kruskal_tests.csv"))

cell_model_df <- target_df %>%
  dplyr::filter(.data$Ploidy %in% c("2N", "4N"), .data$DosePanel %in% c("0", "30", "120"), is.finite(.data$pseudotime)) %>%
  dplyr::mutate(Ploidy = stats::relevel(factor(.data$Ploidy), ref = "2N"), DosePanel = stats::relevel(factor(.data$DosePanel), ref = "0"))
cell_model <- fit_lm_safe(pseudotime ~ Ploidy * DosePanel, cell_model_df)
write_csv_04f(tidy_lm_coefficients_04f(cell_model, "cell_level_pseudotime_initial_ploidy_dose"), file.path(out_response, "initial_ploidy_dose_cell_level_lm_coefficients.csv"))
write_csv_04f(tidy_anova_04f(cell_model, "cell_level_pseudotime_initial_ploidy_dose"), file.path(out_response, "initial_ploidy_dose_cell_level_lm_anova.csv"))
plot_initial_ploidy_cell_response(target_df, file.path(out_response, "initial_ploidy_dose_cell_level_pseudotime_violin_box"))

write_csv_04f(sample_info_base, file.path(out_response, "sample_level_dose_response_summary.csv"))
sample_metric_models <- lapply(c("target_cluster_percent", "median_pseudotime", "mean_cell_ploidy"), function(metric) {
  df <- sample_info_base %>%
    dplyr::filter(is.finite(.data[[metric]]), .data$Ploidy %in% c("2N", "4N"), .data$DosePanel %in% c("0", "30", "120")) %>%
    dplyr::mutate(Ploidy = stats::relevel(factor(.data$Ploidy), ref = "2N"), DosePanel = stats::relevel(factor(.data$DosePanel), ref = "0"))
  model <- fit_lm_safe(stats::as.formula(paste(metric, "~ Ploidy * DosePanel")), df)
  list(
    coefficients = tidy_lm_coefficients_04f(model, paste0("sample_level_", metric, "_initial_ploidy_dose")),
    anova = tidy_anova_04f(model, paste0("sample_level_", metric, "_initial_ploidy_dose"))
  )
})
write_csv_04f(dplyr::bind_rows(lapply(sample_metric_models, `[[`, "coefficients")), file.path(out_response, "sample_level_dose_response_lm_coefficients.csv"))
write_csv_04f(dplyr::bind_rows(lapply(sample_metric_models, `[[`, "anova")), file.path(out_response, "sample_level_dose_response_lm_anova.csv"))
plot_sample_level_response(sample_info_base, file.path(out_response, "sample_level_dose_response"))

ploidy_threshold <- stats::median(target_df$cell_ploidy[is.finite(target_df$cell_ploidy)], na.rm = TRUE)
target_high_low <- target_df %>%
  dplyr::filter(is.finite(.data$cell_ploidy)) %>%
  dplyr::mutate(
    cell_ploidy_group = ifelse(.data$cell_ploidy > ploidy_threshold, "High Cell Ploidy", "Low Cell Ploidy"),
    cell_ploidy_group = factor(.data$cell_ploidy_group, levels = c("Low Cell Ploidy", "High Cell Ploidy")),
    high_cell_ploidy = as.integer(.data$cell_ploidy > ploidy_threshold)
  )
write_csv_04f(
  data.frame(cell_ploidy_median_threshold = ploidy_threshold, n_cells_with_finite_cell_ploidy = nrow(target_high_low)),
  file.path(out_response, "cell_ploidy_high_low_threshold.csv")
)

high_low_cell_summary <- target_high_low %>%
  dplyr::group_by(.data$Ploidy, .data$DosePanel, .data$cell_ploidy_group) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    n_samples = dplyr::n_distinct(.data$sampleID),
    mean_pseudotime = mean(.data$pseudotime, na.rm = TRUE),
    median_pseudotime = stats::median(.data$pseudotime, na.rm = TRUE),
    mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_04f(high_low_cell_summary, file.path(out_response, "cell_ploidy_high_low_cell_level_summary.csv"))

sample_high_low <- target_high_low %>%
  dplyr::group_by(.data$sampleID, .data$Ploidy, .data$DosePanel) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    n_high_cell_ploidy = sum(.data$high_cell_ploidy, na.rm = TRUE),
    high_cell_ploidy_fraction = .data$n_high_cell_ploidy / .data$n_cells,
    high_cell_ploidy_percent = 100 * .data$high_cell_ploidy_fraction,
    mean_pseudotime = mean(.data$pseudotime, na.rm = TRUE),
    median_pseudotime = stats::median(.data$pseudotime, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_04f(sample_high_low, file.path(out_response, "cell_ploidy_high_low_sample_level_summary.csv"))

high_low_cell_model_df <- target_high_low %>%
  dplyr::filter(.data$Ploidy %in% c("2N", "4N"), .data$DosePanel %in% c("0", "30", "120"), is.finite(.data$pseudotime)) %>%
  dplyr::mutate(
    Ploidy = stats::relevel(factor(.data$Ploidy), ref = "2N"),
    DosePanel = stats::relevel(factor(.data$DosePanel), ref = "0"),
    cell_ploidy_group = stats::relevel(factor(.data$cell_ploidy_group), ref = "Low Cell Ploidy")
  )
high_low_cell_model <- fit_lm_safe(pseudotime ~ Ploidy * DosePanel * cell_ploidy_group, high_low_cell_model_df)
write_csv_04f(tidy_lm_coefficients_04f(high_low_cell_model, "cell_level_pseudotime_initial_ploidy_dose_cell_ploidy_high_low"), file.path(out_response, "cell_ploidy_high_low_cell_level_lm_coefficients.csv"))
write_csv_04f(tidy_anova_04f(high_low_cell_model, "cell_level_pseudotime_initial_ploidy_dose_cell_ploidy_high_low"), file.path(out_response, "cell_ploidy_high_low_cell_level_lm_anova.csv"))

high_fraction_model_df <- sample_high_low %>%
  dplyr::filter(.data$Ploidy %in% c("2N", "4N"), .data$DosePanel %in% c("0", "30", "120"), is.finite(.data$high_cell_ploidy_percent)) %>%
  dplyr::mutate(Ploidy = stats::relevel(factor(.data$Ploidy), ref = "2N"), DosePanel = stats::relevel(factor(.data$DosePanel), ref = "0"))
high_fraction_model <- fit_lm_safe(high_cell_ploidy_percent ~ Ploidy * DosePanel, high_fraction_model_df)
write_csv_04f(tidy_lm_coefficients_04f(high_fraction_model, "sample_level_high_cell_ploidy_fraction_initial_ploidy_dose"), file.path(out_response, "cell_ploidy_high_low_sample_level_lm_coefficients.csv"))
write_csv_04f(tidy_anova_04f(high_fraction_model, "sample_level_high_cell_ploidy_fraction_initial_ploidy_dose"), file.path(out_response, "cell_ploidy_high_low_sample_level_lm_anova.csv"))
plot_high_low_response(target_high_low, sample_high_low, file.path(out_response, "cell_ploidy_high_low_dose_response"))

message("Generating full prompt downstream outputs.")

inventory <- make_input_inventory(raw_scvelo_input_root, input_file, prompt_data_dir)
metadata_detected <- data.frame(
  standardized_role = c("cell_id", "sample_id", "cell_type", "cluster", "dose", "initial_ploidy_group", "cell_ploidy", "pseudotime", "umap1", "umap2"),
  selected_column = c("cell", "sample", "TN", "clusters", "Dose", "Ploidy", "cell_ploidy", "pseudotime/velocity_pseudotime", "UMAP_1", "UMAP_2"),
  source_file = input_file,
  notes = c(
    "Merged 04e scVelo all_cells table.",
    "Explicit sample metadata column.",
    "Tumor extracted from all_cells using TN == Tumor.",
    paste0("Target values are all Tumor clusters excluding ", excluded_cluster_label, "."),
    "Explicit Dose metadata, standardized to 0/30/120mg/kg.",
    "Explicit initial Ploidy metadata, standardized to 2N/4N.",
    "Cell Ploidy merged in 04e from all_ploidy.tsv.",
    "scVelo velocity pseudotime carried as pseudotime in 04e output.",
    "UMAP coordinate from scVelo all_cells metadata.",
    "UMAP coordinate from scVelo all_cells metadata."
  ),
  stringsAsFactors = FALSE
)
write_csv_04f(metadata_detected, file.path(prompt_data_dir, "metadata_columns_detected.tsv"))

pt_range <- range(target_df$pseudotime, na.rm = TRUE)
pt_span <- diff(pt_range)
if (!is.finite(pt_span) || pt_span <= 0) pt_span <- 1
target_meta <- target_df %>%
  dplyr::mutate(
    cell_id = as.character(.data$cell),
    sample_id = as.character(.data$sampleID),
    cell_type_original = as.character(.data$TN),
    cluster_original = as.character(.data$clusters),
    cluster_target = as.character(.data$clusters),
    initial_ploidy_group = factor(as.character(.data$Ploidy), levels = c("2N", "4N")),
    initial_ploidy_numeric = ifelse(as.character(.data$Ploidy) == "2N", 2, ifelse(as.character(.data$Ploidy) == "4N", 4, NA_real_)),
    dose = factor(as.character(.data$Dose), levels = c("0mg/kg", "30mg/kg", "120mg/kg")),
    dose_numeric = safe_num(as.character(.data$DosePanel)),
    pseudotime_raw = safe_num(.data$pseudotime),
    pseudotime = (.data$pseudotime_raw - pt_range[1]) / pt_span,
    relative_cell_ploidy = .data$cell_ploidy / .data$initial_ploidy_numeric,
    umap1 = safe_num(.data$UMAP_1),
    umap2 = safe_num(.data$UMAP_2),
    pseudotime_bin_index = pmax(1L, pmin(10L, findInterval(.data$pseudotime, seq(0, 1, by = 0.1), rightmost.closed = TRUE, all.inside = TRUE))),
    pseudotime_bin_10 = factor(sprintf("bin_%02d", .data$pseudotime_bin_index), levels = sprintf("bin_%02d", 1:10)),
    bin_start = (.data$pseudotime_bin_index - 1) / 10,
    bin_end = .data$pseudotime_bin_index / 10,
    bin_mid = (.data$bin_start + .data$bin_end) / 2,
    pseudotime_stage = dplyr::case_when(
      .data$pseudotime < 1 / 3 ~ "early",
      .data$pseudotime < 2 / 3 ~ "middle",
      TRUE ~ "late"
    ),
    pseudotime_stage = factor(.data$pseudotime_stage, levels = c("early", "middle", "late")),
    cluster_target = factor(.data$cluster_target, levels = target_clusters)
  ) %>%
  dplyr::select(
    "cell_id", "sample_id", "cell_type_original", "cluster_original", "cluster_target",
    "initial_ploidy_group", "initial_ploidy_numeric", "dose", "dose_numeric",
    "pseudotime_raw", "pseudotime", "cell_ploidy", "relative_cell_ploidy",
    "umap1", "umap2", "pseudotime_bin_10", "bin_start", "bin_end", "bin_mid", "pseudotime_stage"
  )

missingness <- data.frame(
  field = c("sample_id", "dose", "initial_ploidy_group", "pseudotime", "cell_ploidy", "cluster_target", "umap1", "umap2"),
  n_missing = c(
    sum(is.na(target_meta$sample_id) | !nzchar(target_meta$sample_id)),
    sum(is.na(target_meta$dose)),
    sum(is.na(target_meta$initial_ploidy_group)),
    sum(!is.finite(target_meta$pseudotime)),
    sum(!is.finite(target_meta$cell_ploidy)),
    sum(is.na(target_meta$cluster_target)),
    sum(!is.finite(target_meta$umap1)),
    sum(!is.finite(target_meta$umap2))
  ),
  n_total = nrow(target_meta),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(fraction_missing = .data$n_missing / .data$n_total)
write_csv_04f(target_meta, file.path(prompt_data_dir, "target_cells_metadata.tsv"))
write_csv_04f(missingness, file.path(prompt_data_dir, "missingness_summary.tsv"))

counts_by_sample <- target_meta %>%
  dplyr::count(.data$sample_id, .data$initial_ploidy_group, .data$dose, name = "n_target_cells") %>%
  dplyr::arrange(.data$initial_ploidy_group, .data$dose, .data$sample_id)
counts_by_group_cluster <- target_meta %>%
  dplyr::count(.data$initial_ploidy_group, .data$dose, .data$cluster_target, name = "n_target_cells")
counts_by_sample_cluster <- target_meta %>%
  dplyr::count(.data$sample_id, .data$initial_ploidy_group, .data$dose, .data$cluster_target, name = "n_target_cells")
counts_by_sample_bin <- target_meta %>%
  dplyr::count(.data$sample_id, .data$initial_ploidy_group, .data$dose, .data$pseudotime_bin_10, name = "n_target_cells")
write_csv_04f(counts_by_sample, file.path(prompt_data_dir, "counts_by_sample.tsv"))
write_csv_04f(counts_by_group_cluster, file.path(prompt_data_dir, "counts_by_initial_ploidy_dose_cluster.tsv"))
write_csv_04f(counts_by_sample_cluster, file.path(prompt_data_dir, "counts_by_sample_cluster.tsv"))
write_csv_04f(counts_by_sample_bin, file.path(prompt_data_dir, "counts_by_sample_pseudotime_bin.tsv"))

umap_base <- target_meta %>% dplyr::filter(is.finite(.data$umap1), is.finite(.data$umap2))
p_fig1 <- ggplot(umap_base, aes(x = .data$umap1, y = .data$umap2, color = .data$pseudotime)) +
  geom_point(size = 0.2, alpha = 0.8) +
  facet_wrap(~initial_ploidy_group, nrow = 1) +
  scale_color_viridis_c(option = "magma", name = "Pseudotime") +
  labs(title = "Tumor cells only; non-cell-cycle clusters; scVelo all-cells-derived pseudotime", x = "UMAP 1", y = "UMAP 2") +
  coord_equal() +
  theme_04f(9)
save_prompt_figure(p_fig1, prompt_pdf_dir, prompt_png_dir, "fig1_umap_pseudotime_target_cells", width = 8.5, height = 4.2)

p_fig1_dose <- ggplot(umap_base, aes(x = .data$umap1, y = .data$umap2, color = .data$pseudotime)) +
  geom_point(size = 0.2, alpha = 0.8) +
  facet_wrap(~dose, nrow = 1) +
  scale_color_viridis_c(option = "magma", name = "Pseudotime") +
  labs(title = "Target Tumor cells pseudotime by dose", x = "UMAP 1", y = "UMAP 2") +
  coord_equal() +
  theme_04f(9)
save_prompt_figure(p_fig1_dose, prompt_pdf_dir, prompt_png_dir, "fig1b_umap_pseudotime_target_cells_by_dose", width = 10.5, height = 4.1)

p_fig1_dose_ploidy <- ggplot(umap_base, aes(x = .data$umap1, y = .data$umap2, color = .data$cell_ploidy)) +
  geom_point(size = 0.2, alpha = 0.8) +
  facet_wrap(~dose, nrow = 1) +
  scale_color_gradientn(colors = c("#313695", "#74ADD1", "#FEE08B", "#F46D43", "#A50026"), name = "Cell Ploidy") +
  labs(title = "Target Tumor cells Cell Ploidy by dose", x = "UMAP 1", y = "UMAP 2") +
  coord_equal() +
  theme_04f(9)
save_prompt_figure(p_fig1_dose_ploidy, prompt_pdf_dir, prompt_png_dir, "fig1c_umap_cell_ploidy_target_cells_by_dose", width = 10.5, height = 4.1)
save_prompt_figure(p_fig1_dose_ploidy, prompt_pdf_dir, prompt_png_dir, "umap_cell_ploidy_target_cells_by_dose", width = 10.5, height = 4.1)

p_fig2 <- ggplot(umap_base, aes(x = .data$umap1, y = .data$umap2, color = .data$cluster_target)) +
  geom_point(size = 0.24, alpha = 0.82) +
  facet_wrap(~initial_ploidy_group, nrow = 1) +
  scale_color_manual(values = cluster_palette, drop = FALSE) +
  labs(title = "Tumor cells only; non-cell-cycle clusters; cluster annotation", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
  coord_equal() +
  theme_04f(9)
save_prompt_figure(p_fig2, prompt_pdf_dir, prompt_png_dir, "fig2_umap_clusters_target_cells", width = 8.5, height = 4.2)

p_fig3 <- ggplot(umap_base, aes(x = .data$umap1, y = .data$umap2, color = .data$cell_ploidy)) +
  geom_point(size = 0.22, alpha = 0.82) +
  facet_grid(initial_ploidy_group ~ dose) +
  scale_color_gradientn(colors = c("#313695", "#74ADD1", "#FEE08B", "#F46D43", "#A50026"), name = "Cell Ploidy") +
  labs(title = "Tumor cells only; non-cell-cycle clusters; Cell Ploidy by dose", x = "UMAP 1", y = "UMAP 2") +
  coord_equal() +
  theme_04f(8)
save_prompt_figure(p_fig3, prompt_pdf_dir, prompt_png_dir, "fig3_umap_cell_ploidy_by_dose", width = 9.6, height = 6.4)

p_fig3_time <- ggplot(umap_base, aes(x = .data$umap1, y = .data$umap2, color = .data$pseudotime)) +
  geom_point(size = 0.22, alpha = 0.82) +
  facet_grid(initial_ploidy_group ~ dose) +
  scale_color_viridis_c(option = "magma", name = "Pseudotime") +
  labs(title = "Tumor cells only; non-cell-cycle clusters; scVelo pseudotime by dose", x = "UMAP 1", y = "UMAP 2") +
  coord_equal() +
  theme_04f(8)
save_prompt_figure(p_fig3_time, prompt_pdf_dir, prompt_png_dir, "fig3b_umap_pseudotime_by_dose", width = 9.6, height = 6.4)
save_prompt_figure(p_fig3_time, prompt_pdf_dir, prompt_png_dir, "umap_pseudotime_by_initial_ploidy_dose", width = 9.6, height = 6.4)

p_fig4_counts <- ggplot(counts_by_sample_cluster, aes(x = .data$sample_id, y = .data$n_target_cells, fill = .data$cluster_target)) +
  geom_col(width = 0.72) +
  facet_grid(initial_ploidy_group ~ dose, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cluster_palette, drop = FALSE) +
  labs(title = "Target-cell counts per sample", x = NULL, y = "Cells", fill = "Cluster") +
  theme_04f(8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_prompt_figure(p_fig4_counts, prompt_pdf_dir, prompt_png_dir, "fig4_cell_counts_per_sample", width = 11.5, height = 5.5)

invisible(file.copy(file.path(out_hist, "all_samples_pseudotime_frequency_histograms_16x1.pdf"), file.path(prompt_pdf_dir, "fig4_pseudotime_distribution_per_sample.pdf"), overwrite = TRUE))
invisible(file.copy(file.path(out_hist, "all_samples_pseudotime_frequency_histograms_16x1.png"), file.path(prompt_png_dir, "fig4_pseudotime_distribution_per_sample.png"), overwrite = TRUE))

sample_bin_counts <- target_meta %>%
  dplyr::count(.data$sample_id, .data$pseudotime_bin_10, name = "n_bin") %>%
  dplyr::group_by(.data$sample_id) %>%
  dplyr::mutate(total = sum(.data$n_bin), fraction = .data$n_bin / .data$total) %>%
  dplyr::ungroup()
sample_pseudotime_metrics <- target_meta %>%
  dplyr::group_by(.data$sample_id, .data$initial_ploidy_group, .data$dose) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean_pseudotime = mean(.data$pseudotime, na.rm = TRUE),
    median_pseudotime = stats::median(.data$pseudotime, na.rm = TRUE),
    sd_pseudotime = stats::sd(.data$pseudotime, na.rm = TRUE),
    iqr_pseudotime = stats::IQR(.data$pseudotime, na.rm = TRUE),
    early_fraction = mean(.data$pseudotime_stage == "early", na.rm = TRUE),
    middle_fraction = mean(.data$pseudotime_stage == "middle", na.rm = TRUE),
    late_fraction = mean(.data$pseudotime_stage == "late", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    sample_bin_counts %>%
      dplyr::group_by(.data$sample_id) %>%
      dplyr::summarise(
        pseudotime_entropy = entropy_04f(.data$fraction),
        peak_bin = as.character(.data$pseudotime_bin_10[which.max(.data$fraction)]),
        .groups = "drop"
      ),
    by = "sample_id"
  ) %>%
  dplyr::left_join(
    target_meta %>%
      dplyr::count(.data$sample_id, .data$cluster_target, name = "n_cluster") %>%
      dplyr::group_by(.data$sample_id) %>%
      dplyr::mutate(frac = .data$n_cluster / sum(.data$n_cluster)) %>%
      dplyr::ungroup() %>%
      dplyr::select("sample_id", "cluster_target", "frac") %>%
      tidyr::pivot_wider(names_from = "cluster_target", values_from = "frac", names_prefix = "fraction_cluster_", values_fill = 0),
    by = "sample_id"
  )
write_csv_04f(sample_pseudotime_metrics, file.path(prompt_data_dir, "sample_pseudotime_metrics.tsv"))

wasserstein_rows <- lapply(seq_len(nrow(sample_pseudotime_metrics)), function(i) {
  sid <- sample_pseudotime_metrics$sample_id[i]
  grp <- as.character(sample_pseudotime_metrics$initial_ploidy_group[i])
  dose_i <- as.character(sample_pseudotime_metrics$dose[i])
  x <- target_meta$pseudotime[target_meta$sample_id == sid]
  control_pool <- target_meta %>% dplyr::filter(.data$initial_ploidy_group == grp, .data$dose == "0mg/kg")
  if (identical(dose_i, "0mg/kg") && dplyr::n_distinct(control_pool$sample_id) > 1) {
    control_pool <- control_pool %>% dplyr::filter(.data$sample_id != sid)
  }
  data.frame(
    sample_id = sid,
    initial_ploidy_group = grp,
    dose = dose_i,
    wasserstein_to_control = empirical_wasserstein_04f(x, control_pool$pseudotime),
    n_sample_cells = length(x),
    n_control_cells = nrow(control_pool),
    stringsAsFactors = FALSE
  )
})
sample_wasserstein <- dplyr::bind_rows(wasserstein_rows)
write_csv_04f(sample_wasserstein, file.path(prompt_data_dir, "sample_pseudotime_wasserstein.tsv"))

sample_wasserstein_plot <- add_dose_x_04f(sample_wasserstein, "dose")
sample_wasserstein_sig <- make_group_significance_04f(sample_wasserstein_plot, "wasserstein_to_control", dose_col = "dose", group_col = "initial_ploidy_group")
write_csv_04f(sample_wasserstein_sig, file.path(prompt_stats_dir, "pseudotime_wasserstein_pairwise_significance.tsv"))
p_wasserstein <- ggplot(sample_wasserstein_plot, aes(x = .data$dose_x, y = .data$wasserstein_to_control, color = .data$initial_ploidy_group, group = interaction(.data$dose, .data$initial_ploidy_group))) +
  geom_boxplot(position = position_dodge(width = 0.56), outlier.shape = NA, width = 0.48) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.55), size = 2.2) +
  dose_scale_04f() +
  scale_color_manual(values = c("2N" = "#1D4E89", "4N" = "#A84600"), drop = FALSE) +
  labs(title = "Sample pseudotime Wasserstein distance to matched 0mg/kg controls", x = "Dose", y = "Wasserstein distance", color = "Initial Ploidy") +
  theme_04f(9) +
  theme(legend.position = "bottom")
p_wasserstein <- add_significance_layers_04f(p_wasserstein, sample_wasserstein_sig, label_size = 2.3)
save_prompt_figure(p_wasserstein, prompt_pdf_dir, prompt_png_dir, "pseudotime_wasserstein_by_group", width = 7.8, height = 5.2)

cluster_fraction_cols <- intersect(paste0("fraction_cluster_", target_clusters), names(sample_pseudotime_metrics))
sample_metric_cols <- c("mean_pseudotime", "median_pseudotime", "early_fraction", "middle_fraction", "late_fraction", "pseudotime_entropy", cluster_fraction_cols)
sample_metric_long <- sample_pseudotime_metrics %>%
  dplyr::select("sample_id", "initial_ploidy_group", "dose", "n_cells", dplyr::all_of(sample_metric_cols)) %>%
  tidyr::pivot_longer(cols = -c("sample_id", "initial_ploidy_group", "dose", "n_cells"), names_to = "metric", values_to = "value")
sample_metric_long_plot <- add_dose_x_04f(sample_metric_long, "dose")
sample_metric_sig <- make_group_significance_04f(sample_metric_long_plot, "value", dose_col = "dose", group_col = "initial_ploidy_group", facet_cols = "metric")
write_csv_04f(sample_metric_sig, file.path(prompt_stats_dir, "pseudotime_sample_metrics_pairwise_significance.tsv"))
p_fig5 <- ggplot(sample_metric_long_plot, aes(x = .data$dose_x, y = .data$value, color = .data$initial_ploidy_group, group = interaction(.data$dose, .data$initial_ploidy_group))) +
  geom_boxplot(position = position_dodge(width = 0.56), outlier.shape = NA, width = 0.48, alpha = 0.1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.55), size = 2) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  dose_scale_04f() +
  scale_color_manual(values = c("2N" = "#1D4E89", "4N" = "#A84600"), drop = FALSE) +
  labs(title = "Sample-level pseudotime and composition metrics", x = "Dose", y = NULL, color = "Initial Ploidy") +
  theme_04f(8) +
  theme(legend.position = "bottom")
p_fig5 <- add_significance_layers_04f(p_fig5, sample_metric_sig, label_size = 1.8)
save_prompt_figure(p_fig5, prompt_pdf_dir, prompt_png_dir, "fig5_sample_pseudotime_metrics", width = 12, height = 8)
save_prompt_figure(p_fig5, prompt_pdf_dir, prompt_png_dir, "pseudotime_sample_metrics_boxplots", width = 12, height = 8)

target_meta_dose_plot <- add_dose_x_04f(target_meta, "dose")
fig6_sig <- make_dose_significance_04f(target_meta_dose_plot, "cell_ploidy", dose_col = "dose", facet_cols = "initial_ploidy_group")
write_csv_04f(fig6_sig, file.path(prompt_stats_dir, "fig6_ploidy_violin_pairwise_significance.tsv"))
p_fig6 <- ggplot(target_meta_dose_plot, aes(x = .data$dose_x, y = .data$cell_ploidy, fill = .data$dose, group = .data$dose)) +
  geom_violin(alpha = 0.22, trim = TRUE, linewidth = 0.2) +
  geom_boxplot(width = 0.18, outlier.size = 0.15, alpha = 0.65) +
  geom_point(
    data = sample_info_base %>% dplyr::mutate(dose = factor(as.character(.data$Dose), levels = c("0mg/kg", "30mg/kg", "120mg/kg")), initial_ploidy_group = .data$Ploidy) %>% add_dose_x_04f("dose"),
    aes(x = .data$dose_x, y = .data$mean_cell_ploidy),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.07),
    size = 2,
    color = "black"
  ) +
  facet_wrap(~initial_ploidy_group, nrow = 1) +
  dose_scale_04f() +
  scale_fill_manual(values = c("0mg/kg" = "#4C78A8", "30mg/kg" = "#F58518", "120mg/kg" = "#54A24B"), drop = FALSE) +
  labs(title = "Cell Ploidy by dose with sample-level overlays", x = "Dose", y = "Cell Ploidy", fill = "Dose") +
  theme_04f(9) +
  theme(legend.position = "bottom")
p_fig6 <- add_significance_layers_04f(p_fig6, fig6_sig, label_size = 2.2)
save_prompt_figure(p_fig6, prompt_pdf_dir, prompt_png_dir, "fig6_ploidy_violin_box_sample_overlay", width = 8.8, height = 5.2)

cluster_ploidy_sig <- make_group_significance_04f(target_meta_dose_plot, "cell_ploidy", dose_col = "dose", group_col = "initial_ploidy_group", facet_cols = "cluster_target")
write_csv_04f(cluster_ploidy_sig, file.path(prompt_stats_dir, "ploidy_by_cluster_and_dose_pairwise_significance.tsv"))
p_fig7_cluster_ploidy <- ggplot(target_meta_dose_plot, aes(x = .data$dose_x, y = .data$cell_ploidy, fill = .data$initial_ploidy_group, group = interaction(.data$dose, .data$initial_ploidy_group))) +
  geom_violin(position = position_dodge(width = 0.56), alpha = 0.24, trim = TRUE, linewidth = 0.2) +
  geom_boxplot(position = position_dodge(width = 0.56), width = 0.18, outlier.size = 0.12, alpha = 0.72) +
  facet_wrap(~cluster_target, nrow = cluster_facet_nrow) +
  dose_scale_04f() +
  scale_fill_manual(values = c("2N" = "#1D4E89", "4N" = "#A84600"), drop = FALSE) +
  labs(title = "Cell Ploidy by dose and target cluster", x = "Dose", y = "Cell Ploidy", fill = "Initial Ploidy") +
  theme_04f(8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
p_fig7_cluster_ploidy <- add_significance_layers_04f(p_fig7_cluster_ploidy, cluster_ploidy_sig, label_size = 2)
save_prompt_figure(p_fig7_cluster_ploidy, prompt_pdf_dir, prompt_png_dir, "ploidy_by_cluster_and_dose", width = 11.2, height = 5.6)

pseudo_models <- lapply(
  setdiff(unique(sample_metric_long$metric), character(0)),
  function(metric_name) {
    df <- sample_metric_long %>%
      dplyr::filter(.data$metric == metric_name, is.finite(.data$value)) %>%
      dplyr::mutate(initial_ploidy_group = stats::relevel(factor(.data$initial_ploidy_group), ref = "2N"), dose = stats::relevel(factor(.data$dose), ref = "0mg/kg"))
    fit_robust_lm_table(value ~ initial_ploidy_group * dose, df, metric_name, "pseudotime_sample_metric", notes = "Sample-level OLS with HC3 robust standard errors.")$table
  }
)
pseudotime_metric_models <- dplyr::bind_rows(pseudo_models)
write_csv_04f(pseudotime_metric_models, file.path(prompt_stats_dir, "pseudotime_metric_models.tsv"))

ploidy_thresholds <- target_meta %>%
  dplyr::filter(.data$dose == "0mg/kg") %>%
  dplyr::group_by(.data$initial_ploidy_group) %>%
  dplyr::summarise(control_q90_cell_ploidy = stats::quantile(.data$cell_ploidy, 0.9, na.rm = TRUE), .groups = "drop")
sample_ploidy_metrics <- target_meta %>%
  dplyr::left_join(ploidy_thresholds, by = "initial_ploidy_group") %>%
  dplyr::group_by(.data$sample_id, .data$initial_ploidy_group, .data$dose) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
    median_cell_ploidy = stats::median(.data$cell_ploidy, na.rm = TRUE),
    sd_cell_ploidy = stats::sd(.data$cell_ploidy, na.rm = TRUE),
    iqr_cell_ploidy = stats::IQR(.data$cell_ploidy, na.rm = TRUE),
    cv_cell_ploidy = .data$sd_cell_ploidy / .data$mean_cell_ploidy,
    q10_cell_ploidy = stats::quantile(.data$cell_ploidy, 0.10, na.rm = TRUE),
    q25_cell_ploidy = stats::quantile(.data$cell_ploidy, 0.25, na.rm = TRUE),
    q75_cell_ploidy = stats::quantile(.data$cell_ploidy, 0.75, na.rm = TRUE),
    q90_cell_ploidy = stats::quantile(.data$cell_ploidy, 0.90, na.rm = TRUE),
    upper_tail_fraction = mean(.data$cell_ploidy > first_present_num(.data$control_q90_cell_ploidy), na.rm = TRUE),
    mean_relative_cell_ploidy = mean(.data$relative_cell_ploidy, na.rm = TRUE),
    median_relative_cell_ploidy = stats::median(.data$relative_cell_ploidy, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_04f(sample_ploidy_metrics, file.path(prompt_data_dir, "sample_ploidy_metrics.tsv"))

ploidy_metric_long <- sample_ploidy_metrics %>%
  dplyr::select("sample_id", "initial_ploidy_group", "dose", "n_cells", "mean_cell_ploidy", "median_cell_ploidy", "upper_tail_fraction", "mean_relative_cell_ploidy", "median_relative_cell_ploidy") %>%
  tidyr::pivot_longer(cols = -c("sample_id", "initial_ploidy_group", "dose", "n_cells"), names_to = "metric", values_to = "value")
ploidy_metric_long_plot <- add_dose_x_04f(ploidy_metric_long, "dose")
ploidy_metric_sig <- make_group_significance_04f(ploidy_metric_long_plot, "value", dose_col = "dose", group_col = "initial_ploidy_group", facet_cols = "metric")
write_csv_04f(ploidy_metric_sig, file.path(prompt_stats_dir, "ploidy_sample_metrics_pairwise_significance.tsv"))
ploidy_metric_models <- lapply(unique(ploidy_metric_long$metric), function(metric_name) {
  df <- ploidy_metric_long %>%
    dplyr::filter(.data$metric == metric_name, is.finite(.data$value)) %>%
    dplyr::mutate(initial_ploidy_group = stats::relevel(factor(.data$initial_ploidy_group), ref = "2N"), dose = stats::relevel(factor(.data$dose), ref = "0mg/kg"))
  fit_robust_lm_table(value ~ initial_ploidy_group * dose, df, metric_name, "ploidy_sample_metric", notes = "Sample-level OLS with HC3 robust standard errors.")$table
})
write_csv_04f(dplyr::bind_rows(ploidy_metric_models), file.path(prompt_stats_dir, "ploidy_sample_metric_models.tsv"))

p_ploidy_metrics <- ggplot(ploidy_metric_long_plot, aes(x = .data$dose_x, y = .data$value, color = .data$initial_ploidy_group, group = interaction(.data$dose, .data$initial_ploidy_group))) +
  geom_boxplot(position = position_dodge(width = 0.56), outlier.shape = NA, width = 0.48, alpha = 0.1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.55), size = 2) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  dose_scale_04f() +
  scale_color_manual(values = c("2N" = "#1D4E89", "4N" = "#A84600"), drop = FALSE) +
  labs(title = "Sample-level Cell Ploidy metrics", x = "Dose", y = NULL, color = "Initial Ploidy") +
  theme_04f(8) +
  theme(legend.position = "bottom")
p_ploidy_metrics <- add_significance_layers_04f(p_ploidy_metrics, ploidy_metric_sig, label_size = 1.9)
save_prompt_figure(p_ploidy_metrics, prompt_pdf_dir, prompt_png_dir, "ploidy_sample_metrics_boxplots", width = 11.5, height = 7.2)

upper_tail_plot <- ploidy_metric_long_plot %>% dplyr::filter(.data$metric == "upper_tail_fraction")
upper_tail_sig <- make_group_significance_04f(upper_tail_plot, "value", dose_col = "dose", group_col = "initial_ploidy_group")
write_csv_04f(upper_tail_sig, file.path(prompt_stats_dir, "ploidy_upper_tail_fraction_pairwise_significance.tsv"))
p_upper_tail <- ggplot(upper_tail_plot, aes(x = .data$dose_x, y = .data$value, color = .data$initial_ploidy_group, group = interaction(.data$dose, .data$initial_ploidy_group))) +
  geom_boxplot(position = position_dodge(width = 0.56), outlier.shape = NA, width = 0.48, alpha = 0.1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.55), size = 2.3) +
  dose_scale_04f() +
  scale_color_manual(values = c("2N" = "#1D4E89", "4N" = "#A84600"), drop = FALSE) +
  labs(title = "Sample-level upper-tail Cell Ploidy fraction", x = "Dose", y = "Upper-tail fraction", color = "Initial Ploidy") +
  theme_04f(9) +
  theme(legend.position = "bottom")
p_upper_tail <- add_significance_layers_04f(p_upper_tail, upper_tail_sig, label_size = 2.3)
save_prompt_figure(p_upper_tail, prompt_pdf_dir, prompt_png_dir, "ploidy_upper_tail_fraction", width = 7.8, height = 5.2)

cell_model_prompt_df <- target_meta %>%
  dplyr::mutate(initial_ploidy_group = stats::relevel(factor(.data$initial_ploidy_group), ref = "2N"), dose = stats::relevel(factor(.data$dose), ref = "0mg/kg"), cluster_target = factor(.data$cluster_target))
cell_level_ploidy_model <- fit_robust_lm_table(
  cell_ploidy ~ initial_ploidy_group * dose + pseudotime + cluster_target,
  cell_model_prompt_df,
  "cell_ploidy",
  "exploratory_cell_level_ploidy",
  cluster_col = "sample_id",
  notes = "Exploratory cell-level OLS with sample-cluster robust standard errors; sample-aware summaries are primary."
)
write_model_summary(cell_level_ploidy_model$model, file.path(prompt_models_dir, "ploidy_cell_level_model_summary.txt"))
write_csv_04f(cell_level_ploidy_model$table, file.path(prompt_stats_dir, "ploidy_cell_level_model_coefficients.tsv"))

sample_bin_ploidy <- target_meta %>%
  dplyr::group_by(.data$sample_id, .data$initial_ploidy_group, .data$dose, .data$cluster_target, .data$pseudotime_bin_10, .data$bin_start, .data$bin_end, .data$bin_mid) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
    median_cell_ploidy = stats::median(.data$cell_ploidy, na.rm = TRUE),
    mean_relative_cell_ploidy = mean(.data$relative_cell_ploidy, na.rm = TRUE),
    median_relative_cell_ploidy = stats::median(.data$relative_cell_ploidy, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_04f(sample_bin_ploidy, file.path(prompt_data_dir, "pseudotime_bin_ploidy_by_sample.tsv"))

bin_model_df <- sample_bin_ploidy %>%
  dplyr::mutate(initial_ploidy_group = stats::relevel(factor(.data$initial_ploidy_group), ref = "2N"), dose = stats::relevel(factor(.data$dose), ref = "0mg/kg"), pseudotime_bin_10 = factor(.data$pseudotime_bin_10))
spline_fit <- fit_robust_lm_table(
  mean_cell_ploidy ~ initial_ploidy_group * dose * splines::bs(bin_mid, df = 5, degree = 3),
  bin_model_df,
  "mean_cell_ploidy",
  "ploidy_pseudotime_spline",
  weights = bin_model_df$n_cells,
  cluster_col = "sample_id",
  notes = "Sample-bin WLS weighted by n_cells with sample-cluster robust standard errors."
)
bin_fit <- fit_robust_lm_table(
  mean_cell_ploidy ~ initial_ploidy_group * dose * pseudotime_bin_10,
  bin_model_df,
  "mean_cell_ploidy",
  "ploidy_pseudotime_bin_interaction",
  weights = bin_model_df$n_cells,
  cluster_col = "sample_id",
  notes = "Sample-bin WLS weighted by n_cells with sample-cluster robust standard errors."
)
write_model_summary(spline_fit$model, file.path(prompt_models_dir, "ploidy_pseudotime_spline_model_summary.txt"))
write_csv_04f(spline_fit$table, file.path(prompt_stats_dir, "ploidy_pseudotime_spline_model_coefficients.tsv"))
write_csv_04f(bin_fit$table, file.path(prompt_stats_dir, "ploidy_pseudotime_bin_interaction_model.tsv"))

pred_grid <- tidyr::expand_grid(
  initial_ploidy_group = factor(c("2N", "4N"), levels = levels(bin_model_df$initial_ploidy_group)),
  dose = factor(c("0mg/kg", "30mg/kg", "120mg/kg"), levels = levels(bin_model_df$dose)),
  bin_mid = seq(min(bin_model_df$bin_mid, na.rm = TRUE), max(bin_model_df$bin_mid, na.rm = TRUE), length.out = 100)
)
if (!is.null(spline_fit$model)) {
  pr <- predict(spline_fit$model, newdata = pred_grid, se.fit = TRUE)
  pred_grid$predicted_mean_cell_ploidy <- as.numeric(pr$fit)
  pred_grid$conf_low <- as.numeric(pr$fit - 1.96 * pr$se.fit)
  pred_grid$conf_high <- as.numeric(pr$fit + 1.96 * pr$se.fit)
} else {
  pred_grid$predicted_mean_cell_ploidy <- NA_real_
  pred_grid$conf_low <- NA_real_
  pred_grid$conf_high <- NA_real_
}
write_csv_04f(pred_grid, file.path(prompt_data_dir, "ploidy_pseudotime_model_predictions.tsv"))

p_fig7 <- ggplot() +
  geom_point(data = sample_bin_ploidy, aes(x = .data$bin_mid, y = .data$mean_cell_ploidy, size = .data$n_cells, color = .data$dose), alpha = 0.28) +
  geom_ribbon(data = pred_grid, aes(x = .data$bin_mid, ymin = .data$conf_low, ymax = .data$conf_high, fill = .data$dose), alpha = 0.12) +
  geom_line(data = pred_grid, aes(x = .data$bin_mid, y = .data$predicted_mean_cell_ploidy, color = .data$dose), linewidth = 0.9) +
  facet_wrap(~initial_ploidy_group, nrow = 1) +
  scale_color_manual(values = c("0mg/kg" = "#4C78A8", "30mg/kg" = "#F58518", "120mg/kg" = "#54A24B"), drop = FALSE) +
  scale_fill_manual(values = c("0mg/kg" = "#4C78A8", "30mg/kg" = "#F58518", "120mg/kg" = "#54A24B"), drop = FALSE) +
  labs(title = "Cell Ploidy vs pseudotime smooth model", x = "Normalized scVelo pseudotime", y = "Mean Cell Ploidy", color = "Dose", fill = "Dose", size = "Cells") +
  theme_04f(9) +
  theme(legend.position = "bottom")
save_prompt_figure(p_fig7, prompt_pdf_dir, prompt_png_dir, "fig7_ploidy_vs_pseudotime_smooth", width = 9.5, height = 5.2)
save_prompt_figure(p_fig7, prompt_pdf_dir, prompt_png_dir, "ploidy_vs_pseudotime_smooth_curves", width = 9.5, height = 5.2)

p_by_cluster <- ggplot(sample_bin_ploidy, aes(x = .data$bin_mid, y = .data$mean_cell_ploidy, color = .data$dose, weight = .data$n_cells)) +
  geom_point(aes(size = .data$n_cells), alpha = 0.35) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8) +
  facet_grid(initial_ploidy_group ~ cluster_target) +
  scale_color_manual(values = c("0mg/kg" = "#4C78A8", "30mg/kg" = "#F58518", "120mg/kg" = "#54A24B"), drop = FALSE) +
  labs(title = "Cell Ploidy vs pseudotime by target cluster", x = "Normalized scVelo pseudotime", y = "Mean Cell Ploidy", color = "Dose", size = "Cells") +
  theme_04f(8) +
  theme(legend.position = "bottom")
save_prompt_figure(p_by_cluster, prompt_pdf_dir, prompt_png_dir, "ploidy_vs_pseudotime_by_cluster", width = 11, height = 7)

relative_fit <- fit_robust_lm_table(
  mean_relative_cell_ploidy ~ initial_ploidy_group * dose * splines::bs(bin_mid, df = 5, degree = 3),
  bin_model_df,
  "mean_relative_cell_ploidy",
  "relative_ploidy_pseudotime_spline",
  weights = bin_model_df$n_cells,
  cluster_col = "sample_id",
  notes = "Sensitivity model using relative_cell_ploidy."
)
write_csv_04f(relative_fit$table, file.path(prompt_stats_dir, "relative_ploidy_pseudotime_spline_model_coefficients.tsv"))
p_relative <- ggplot(sample_bin_ploidy, aes(x = .data$bin_mid, y = .data$mean_relative_cell_ploidy, color = .data$dose, weight = .data$n_cells)) +
  geom_point(aes(size = .data$n_cells), alpha = 0.32) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  facet_wrap(~initial_ploidy_group, nrow = 1) +
  scale_color_manual(values = c("0mg/kg" = "#4C78A8", "30mg/kg" = "#F58518", "120mg/kg" = "#54A24B"), drop = FALSE) +
  labs(title = "Relative Cell Ploidy vs pseudotime", x = "Normalized scVelo pseudotime", y = "Mean relative Cell Ploidy", color = "Dose", size = "Cells") +
  theme_04f(9) +
  theme(legend.position = "bottom")
save_prompt_figure(p_relative, prompt_pdf_dir, prompt_png_dir, "relative_ploidy_vs_pseudotime_smooth_curves", width = 9.5, height = 5.2)

baseline_bin <- target_meta %>%
  dplyr::filter(.data$dose == "0mg/kg") %>%
  dplyr::group_by(.data$initial_ploidy_group, .data$pseudotime_bin_10) %>%
  dplyr::summarise(control_mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE), .groups = "drop")
cell_delta <- target_meta %>%
  dplyr::left_join(baseline_bin, by = c("initial_ploidy_group", "pseudotime_bin_10")) %>%
  dplyr::mutate(delta_ploidy_vs_control_bin = .data$cell_ploidy - .data$control_mean_cell_ploidy)
sample_bin_delta <- cell_delta %>%
  dplyr::group_by(.data$sample_id, .data$initial_ploidy_group, .data$dose, .data$pseudotime_bin_10, .data$bin_mid) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean_delta_ploidy_vs_control_bin = mean(.data$delta_ploidy_vs_control_bin, na.rm = TRUE),
    median_delta_ploidy_vs_control_bin = stats::median(.data$delta_ploidy_vs_control_bin, na.rm = TRUE),
    .groups = "drop"
)
write_csv_04f(cell_delta, file.path(prompt_data_dir, "cell_level_delta_ploidy.tsv"))
write_csv_04f(sample_bin_delta, file.path(prompt_data_dir, "sample_bin_delta_ploidy.tsv"))
sample_bin_delta_plot <- add_dose_x_04f(sample_bin_delta, "dose")
sample_bin_delta_sig <- make_group_significance_04f(sample_bin_delta_plot, "mean_delta_ploidy_vs_control_bin", dose_col = "dose", group_col = "initial_ploidy_group", facet_cols = "pseudotime_bin_10")
write_csv_04f(sample_bin_delta_sig, file.path(prompt_stats_dir, "pseudotime_matched_delta_ploidy_pairwise_significance.tsv"))
delta_model_df <- sample_bin_delta %>%
  dplyr::filter(is.finite(.data$mean_delta_ploidy_vs_control_bin)) %>%
  dplyr::mutate(initial_ploidy_group = stats::relevel(factor(.data$initial_ploidy_group), ref = "2N"), dose = stats::relevel(factor(.data$dose), ref = "0mg/kg"), pseudotime_bin_10 = factor(.data$pseudotime_bin_10))
delta_fit <- fit_robust_lm_table(
  mean_delta_ploidy_vs_control_bin ~ initial_ploidy_group * dose * pseudotime_bin_10,
  delta_model_df,
  "mean_delta_ploidy_vs_control_bin",
  "delta_ploidy_bin_model",
  weights = delta_model_df$n_cells,
  cluster_col = "sample_id",
  notes = "Sample-bin WLS of pseudotime-matched delta Cell Ploidy."
)
write_csv_04f(delta_fit$table, file.path(prompt_stats_dir, "delta_ploidy_models.tsv"))

delta_heat <- sample_bin_delta %>%
  dplyr::group_by(.data$initial_ploidy_group, .data$dose, .data$pseudotime_bin_10) %>%
  dplyr::summarise(mean_delta = weighted.mean(.data$mean_delta_ploidy_vs_control_bin, w = .data$n_cells, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(row_label = paste(.data$initial_ploidy_group, .data$dose, sep = " | "))
p_fig8 <- ggplot(delta_heat, aes(x = .data$pseudotime_bin_10, y = .data$row_label, fill = .data$mean_delta)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#313695", mid = "white", high = "#A50026", midpoint = 0, name = "Delta Cell Ploidy") +
  labs(title = "Pseudotime-matched delta Cell Ploidy vs matched 0mg/kg control", x = "Pseudotime bin", y = NULL) +
  theme_04f(8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_prompt_figure(p_fig8, prompt_pdf_dir, prompt_png_dir, "fig8_delta_ploidy_heatmap", width = 8.8, height = 4.8)
save_prompt_figure(p_fig8, prompt_pdf_dir, prompt_png_dir, "pseudotime_matched_delta_ploidy_heatmap", width = 8.8, height = 4.8)

p_delta_box <- ggplot(sample_bin_delta_plot, aes(x = .data$dose_x, y = .data$mean_delta_ploidy_vs_control_bin, color = .data$initial_ploidy_group, group = interaction(.data$dose, .data$initial_ploidy_group))) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
  geom_boxplot(position = position_dodge(width = 0.56), outlier.shape = NA, width = 0.48) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.55), alpha = 0.75, size = 1.7) +
  facet_wrap(~pseudotime_bin_10, nrow = 2) +
  dose_scale_04f() +
  scale_color_manual(values = c("2N" = "#1D4E89", "4N" = "#A84600"), drop = FALSE) +
  labs(title = "Pseudotime-matched delta Cell Ploidy by sample-bin", x = "Dose", y = "Delta Cell Ploidy", color = "Initial Ploidy") +
  theme_04f(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
p_delta_box <- add_significance_layers_04f(p_delta_box, sample_bin_delta_sig, label_size = 1.6)
save_prompt_figure(p_delta_box, prompt_pdf_dir, prompt_png_dir, "pseudotime_matched_delta_ploidy_boxplots", width = 12, height = 7.5)

decompose_once <- function(df, group_value, dose_value) {
  bins <- levels(target_meta$pseudotime_bin_10)
  control <- df %>% dplyr::filter(.data$initial_ploidy_group == group_value, .data$dose == "0mg/kg")
  treat <- df %>% dplyr::filter(.data$initial_ploidy_group == group_value, .data$dose == dose_value)
  if (nrow(control) == 0 || nrow(treat) == 0) return(data.frame())
  summarize_bins <- function(x) {
    x %>%
      dplyr::group_by(.data$pseudotime_bin_10) %>%
      dplyr::summarise(n = dplyr::n(), mu = mean(.data$cell_ploidy, na.rm = TRUE), .groups = "drop") %>%
      tidyr::complete(pseudotime_bin_10 = factor(bins, levels = bins), fill = list(n = 0, mu = NA_real_)) %>%
      dplyr::mutate(p = .data$n / sum(.data$n), mu = ifelse(is.na(.data$mu), mean(x$cell_ploidy, na.rm = TRUE), .data$mu))
  }
  ctab <- summarize_bins(control)
  ttab <- summarize_bins(treat)
  total <- sum(ttab$p * ttab$mu, na.rm = TRUE) - sum(ctab$p * ctab$mu, na.rm = TRUE)
  composition <- sum((ttab$p - ctab$p) * ctab$mu, na.rm = TRUE)
  within <- sum(ttab$p * (ttab$mu - ctab$mu), na.rm = TRUE)
  data.frame(
    initial_ploidy_group = group_value,
    comparison = paste0(dose_value, "_vs_0mg/kg"),
    dose = dose_value,
    total_effect = total,
    composition_effect = composition,
    within_state_effect = within,
    residual = total - composition - within,
    n_control_cells = nrow(control),
    n_treated_cells = nrow(treat),
    stringsAsFactors = FALSE
  )
}
decomp <- dplyr::bind_rows(lapply(c("2N", "4N"), function(g) dplyr::bind_rows(lapply(c("30mg/kg", "120mg/kg"), function(dv) decompose_once(target_meta, g, dv)))))
write_csv_04f(decomp, file.path(prompt_data_dir, "composition_within_state_decomposition.tsv"))

bootstrap_n <- as.integer(get_env_scalar("DECOMPOSITION_BOOTSTRAP_N", "1000"))
set.seed(1234)
boot_rows <- list()
boot_i <- 1L
for (g in c("2N", "4N")) {
  for (dv in c("30mg/kg", "120mg/kg")) {
    sample_pool <- target_meta %>% dplyr::filter(.data$initial_ploidy_group == g, .data$dose %in% c("0mg/kg", dv))
    for (b in seq_len(bootstrap_n)) {
      boot_df <- dplyr::bind_rows(lapply(split(sample_pool, sample_pool$dose), function(xx) {
        sids <- unique(xx$sample_id)
        sampled <- sample(sids, length(sids), replace = TRUE)
        dplyr::bind_rows(lapply(seq_along(sampled), function(j) xx %>% dplyr::filter(.data$sample_id == sampled[j]) %>% dplyr::mutate(sample_id = paste0(.data$sample_id, "_boot", j))))
      }))
      tmp <- decompose_once(boot_df, g, dv)
      if (nrow(tmp) > 0) {
        tmp$bootstrap_id <- b
        boot_rows[[boot_i]] <- tmp
        boot_i <- boot_i + 1L
      }
    }
  }
}
decomp_boot <- dplyr::bind_rows(boot_rows)
write_csv_04f(decomp_boot, file.path(prompt_data_dir, "composition_within_state_decomposition_bootstrap.tsv"))
decomp_ci <- decomp_boot %>%
  tidyr::pivot_longer(cols = c("total_effect", "composition_effect", "within_state_effect"), names_to = "component", values_to = "effect") %>%
  dplyr::group_by(.data$initial_ploidy_group, .data$comparison, .data$component) %>%
  dplyr::summarise(conf_low = stats::quantile(.data$effect, 0.025, na.rm = TRUE), conf_high = stats::quantile(.data$effect, 0.975, na.rm = TRUE), .groups = "drop")
decomp_long <- decomp %>%
  tidyr::pivot_longer(cols = c("total_effect", "composition_effect", "within_state_effect"), names_to = "component", values_to = "effect") %>%
  dplyr::left_join(decomp_ci, by = c("initial_ploidy_group", "comparison", "component"))
p_fig9 <- ggplot(decomp_long, aes(x = .data$comparison, y = .data$effect, fill = .data$component)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_errorbar(aes(ymin = .data$conf_low, ymax = .data$conf_high), position = position_dodge(width = 0.75), width = 0.18) +
  facet_wrap(~initial_ploidy_group, nrow = 1) +
  labs(title = "Composition vs within-state Cell Ploidy decomposition", x = NULL, y = "Effect vs 0mg/kg", fill = "Component") +
  theme_04f(9) +
  theme(legend.position = "bottom")
save_prompt_figure(p_fig9, prompt_pdf_dir, prompt_png_dir, "fig9_composition_within_state_decomposition", width = 9.5, height = 4.8)
save_prompt_figure(p_fig9, prompt_pdf_dir, prompt_png_dir, "decomposition_total_composition_within_barplot", width = 9.5, height = 4.8)
save_prompt_figure(p_fig9, prompt_pdf_dir, prompt_png_dir, "decomposition_bootstrap_confidence_intervals", width = 9.5, height = 4.8)

cluster_abundance <- counts_by_sample_cluster %>%
  dplyr::group_by(.data$sample_id, .data$initial_ploidy_group, .data$dose) %>%
  dplyr::mutate(total_target_cells = sum(.data$n_target_cells), fraction = .data$n_target_cells / .data$total_target_cells, logit_fraction = logit_fraction_04f(.data$n_target_cells, .data$total_target_cells)) %>%
  dplyr::ungroup()
write_csv_04f(cluster_abundance, file.path(prompt_data_dir, "cluster_fraction_by_sample.tsv"))
cluster_abundance_plot <- add_dose_x_04f(cluster_abundance, "dose")
cluster_abundance_sig <- make_group_significance_04f(cluster_abundance_plot, "fraction", dose_col = "dose", group_col = "initial_ploidy_group", facet_cols = "cluster_target")
write_csv_04f(cluster_abundance_sig, file.path(prompt_stats_dir, "cluster_fraction_pairwise_significance.tsv"))
cluster_models <- lapply(levels(target_meta$cluster_target), function(cl) {
  df <- cluster_abundance %>%
    dplyr::filter(.data$cluster_target == cl) %>%
    dplyr::mutate(initial_ploidy_group = stats::relevel(factor(.data$initial_ploidy_group), ref = "2N"), dose = stats::relevel(factor(.data$dose), ref = "0mg/kg"))
  fit_robust_lm_table(logit_fraction ~ initial_ploidy_group * dose, df, paste0("logit_fraction_cluster_", cl), "cluster_abundance", notes = "Sample-level logit fraction model with HC3 robust standard errors.")$table
})
write_csv_04f(dplyr::bind_rows(cluster_models), file.path(prompt_stats_dir, "cluster_abundance_models.tsv"))

bin_abundance <- counts_by_sample_bin %>%
  dplyr::group_by(.data$sample_id, .data$initial_ploidy_group, .data$dose) %>%
  dplyr::mutate(total_target_cells = sum(.data$n_target_cells), fraction = .data$n_target_cells / .data$total_target_cells, logit_bin_fraction = logit_fraction_04f(.data$n_target_cells, .data$total_target_cells)) %>%
  dplyr::ungroup()
bin_abundance_model_df <- bin_abundance %>%
  dplyr::mutate(initial_ploidy_group = stats::relevel(factor(.data$initial_ploidy_group), ref = "2N"), dose = stats::relevel(factor(.data$dose), ref = "0mg/kg"), pseudotime_bin_10 = factor(.data$pseudotime_bin_10))
bin_abundance_fit <- fit_robust_lm_table(
  logit_bin_fraction ~ initial_ploidy_group * dose * pseudotime_bin_10,
  bin_abundance_model_df,
  "logit_bin_fraction",
  "pseudotime_bin_abundance",
  cluster_col = "sample_id",
  notes = "Sample-bin logit fraction model with sample-cluster robust standard errors."
)
write_csv_04f(bin_abundance_fit$table, file.path(prompt_stats_dir, "pseudotime_bin_abundance_models.tsv"))

p_cluster_fraction <- ggplot(cluster_abundance_plot, aes(x = .data$dose_x, y = .data$fraction, color = .data$initial_ploidy_group, group = interaction(.data$dose, .data$initial_ploidy_group))) +
  geom_boxplot(position = position_dodge(width = 0.56), outlier.shape = NA, width = 0.48) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.55), size = 2) +
  facet_wrap(~cluster_target, nrow = cluster_facet_nrow) +
  dose_scale_04f() +
  scale_color_manual(values = c("2N" = "#1D4E89", "4N" = "#A84600"), drop = FALSE) +
  labs(title = "Cluster fraction by sample", x = "Dose", y = "Fraction among target cells", color = "Initial Ploidy") +
  theme_04f(9) +
  theme(legend.position = "bottom")
p_cluster_fraction <- add_significance_layers_04f(p_cluster_fraction, cluster_abundance_sig, label_size = 2.1)
save_prompt_figure(p_cluster_fraction, prompt_pdf_dir, prompt_png_dir, "cluster_fraction_by_sample", width = 9.5, height = 4.8)

bin_heat <- bin_abundance %>%
  dplyr::group_by(.data$initial_ploidy_group, .data$dose, .data$pseudotime_bin_10) %>%
  dplyr::summarise(mean_fraction = mean(.data$fraction, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(row_label = paste(.data$initial_ploidy_group, .data$dose, sep = " | "))
p_bin_heat <- ggplot(bin_heat, aes(x = .data$pseudotime_bin_10, y = .data$row_label, fill = .data$mean_fraction)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(name = "Mean fraction") +
  labs(title = "Pseudotime-bin fraction by group", x = "Pseudotime bin", y = NULL) +
  theme_04f(8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_prompt_figure(p_bin_heat, prompt_pdf_dir, prompt_png_dir, "pseudotime_bin_fraction_heatmap", width = 8.8, height = 4.8)
save_prompt_figure(p_bin_heat, prompt_pdf_dir, prompt_png_dir, "pseudotime_bin_fraction_by_group_heatmap", width = 8.8, height = 4.8)
save_prompt_figure(p_bin_heat, prompt_pdf_dir, prompt_png_dir, "pseudotime_bin_abundance_effects", width = 8.8, height = 4.8)

p_fig10 <- p_cluster_fraction / p_bin_heat + patchwork::plot_annotation(title = "Cluster and pseudotime-bin abundance")
save_prompt_figure(p_fig10, prompt_pdf_dir, prompt_png_dir, "fig10_cluster_and_pseudotime_bin_abundance", width = 10.5, height = 9)

sensitivity_tables <- list()
median_fit <- fit_robust_lm_table(
  median_cell_ploidy ~ initial_ploidy_group * dose * pseudotime_bin_10,
  bin_model_df,
  "median_cell_ploidy",
  "sensitivity_median_cell_ploidy_bin_model",
  weights = bin_model_df$n_cells,
  cluster_col = "sample_id",
  notes = "Sensitivity using median Cell Ploidy."
)
sensitivity_tables[["median"]] <- median_fit$table
cluster_cov_fit <- fit_robust_lm_table(
  mean_cell_ploidy ~ initial_ploidy_group * dose * splines::bs(bin_mid, df = 5, degree = 3) + cluster_target,
  bin_model_df,
  "mean_cell_ploidy",
  "sensitivity_spline_plus_cluster",
  weights = bin_model_df$n_cells,
  cluster_col = "sample_id",
  notes = "Sensitivity with cluster_target covariate."
)
sensitivity_tables[["cluster_covariate"]] <- cluster_cov_fit$table
stage_df <- target_meta %>%
  dplyr::group_by(.data$sample_id, .data$initial_ploidy_group, .data$dose, .data$pseudotime_stage) %>%
  dplyr::summarise(n_cells = dplyr::n(), mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(initial_ploidy_group = stats::relevel(factor(.data$initial_ploidy_group), ref = "2N"), dose = stats::relevel(factor(.data$dose), ref = "0mg/kg"), pseudotime_stage = factor(.data$pseudotime_stage))
stage_fit <- fit_robust_lm_table(
  mean_cell_ploidy ~ initial_ploidy_group * dose * pseudotime_stage,
  stage_df,
  "mean_cell_ploidy",
  "sensitivity_early_middle_late",
  weights = stage_df$n_cells,
  cluster_col = "sample_id",
  notes = "Sensitivity using early/middle/late bins."
)
sensitivity_tables[["stage"]] <- stage_fit$table
numeric_dose_df <- sample_ploidy_metrics %>%
  dplyr::mutate(initial_ploidy_group = stats::relevel(factor(.data$initial_ploidy_group), ref = "2N"), dose_numeric = as.numeric(sub("mg/kg", "", as.character(.data$dose), fixed = TRUE)))
numeric_dose_fit <- fit_robust_lm_table(
  mean_cell_ploidy ~ initial_ploidy_group * dose_numeric + initial_ploidy_group * I(dose_numeric^2),
  numeric_dose_df,
  "mean_cell_ploidy",
  "sensitivity_numeric_quadratic_dose",
  notes = "Secondary numeric dose and quadratic dose model."
)
sensitivity_tables[["numeric_dose"]] <- numeric_dose_fit$table
sensitivity_summary <- dplyr::bind_rows(sensitivity_tables)
write_csv_04f(sensitivity_summary, file.path(prompt_stats_dir, "sensitivity_analysis_summary.tsv"))

signature_status <- data.frame(status = "not_run", reason = NA_character_, stringsAsFactors = FALSE)
gsea_status <- data.frame(status = "not_run", reason = NA_character_, stringsAsFactors = FALSE)
signature_scores_joined <- data.frame()
signature_gene_matching <- data.frame()
signature_source_info <- data.frame()
signatures <- signature_definitions_04f()
h5res <- tryCatch(read_h5ad_signature_expression(h5ad_file, target_meta$cell_id, signatures), error = function(e) list(status = "error", error = conditionMessage(e), expr = NULL, matching = data.frame(), source_info = data.frame()))
signature_status <- data.frame(status = h5res$status, reason = if (!is.null(h5res$error)) h5res$error else NA_character_, stringsAsFactors = FALSE)
if (nrow(h5res$matching) > 0) {
  signature_gene_matching <- h5res$matching
  signature_gene_matching$signature <- vapply(signature_gene_matching$query_gene_clean, function(g) {
    hit <- names(signatures)[vapply(signatures, function(gs) g %in% clean_gene_symbol_04f(gs), logical(1))]
    paste(hit, collapse = ";")
  }, character(1))
  write_csv_04f(signature_gene_matching, file.path(prompt_data_dir, "signature_gene_matching.tsv"))
}
if (nrow(h5res$source_info) > 0) {
  signature_source_info <- h5res$source_info
  write_csv_04f(signature_source_info, file.path(prompt_logs_dir, "expression_source_info.tsv"))
}
if (identical(h5res$status, "ok")) {
  sig_scores <- score_signatures_04f(h5res$expr, signatures, h5res$matching)
  signature_scores_joined <- target_meta %>%
    dplyr::left_join(sig_scores, by = "cell_id")
  write_csv_04f(sig_scores, file.path(prompt_data_dir, "cell_signature_scores.tsv"))
  sig_cols <- intersect(names(signatures), names(signature_scores_joined))
  sig_sample_bin <- signature_scores_joined %>%
    dplyr::group_by(.data$sample_id, .data$initial_ploidy_group, .data$dose, .data$cluster_target, .data$pseudotime_bin_10, .data$bin_mid) %>%
    dplyr::summarise(n_cells = dplyr::n(), dplyr::across(dplyr::all_of(sig_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  sig_sample_cluster <- signature_scores_joined %>%
    dplyr::group_by(.data$sample_id, .data$initial_ploidy_group, .data$dose, .data$cluster_target) %>%
    dplyr::summarise(n_cells = dplyr::n(), dplyr::across(dplyr::all_of(sig_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  write_csv_04f(sig_sample_bin, file.path(prompt_data_dir, "signature_scores_by_sample_bin.tsv"))
  write_csv_04f(sig_sample_cluster, file.path(prompt_data_dir, "signature_scores_by_sample_cluster.tsv"))
  sig_models <- lapply(sig_cols, function(sc) {
    df <- sig_sample_cluster %>%
      dplyr::filter(is.finite(.data[[sc]])) %>%
      dplyr::mutate(initial_ploidy_group = stats::relevel(factor(.data$initial_ploidy_group), ref = "2N"), dose = stats::relevel(factor(.data$dose), ref = "0mg/kg"))
    fit_robust_lm_table(stats::as.formula(paste(sc, "~ initial_ploidy_group * dose")), df, sc, "signature_score_model", weights = df$n_cells, cluster_col = "sample_id", notes = "Signature score model on sample-cluster aggregated values.")$table
  })
  sig_bin_models <- lapply(sig_cols, function(sc) {
    df <- sig_sample_bin %>%
      dplyr::filter(is.finite(.data[[sc]])) %>%
      dplyr::mutate(initial_ploidy_group = stats::relevel(factor(.data$initial_ploidy_group), ref = "2N"), dose = stats::relevel(factor(.data$dose), ref = "0mg/kg"), pseudotime_bin_10 = factor(.data$pseudotime_bin_10))
    fit_robust_lm_table(stats::as.formula(paste(sc, "~ initial_ploidy_group * dose * pseudotime_bin_10")), df, sc, "signature_pseudotime_interaction_model", weights = df$n_cells, cluster_col = "sample_id", notes = "Signature score model on sample-bin aggregated values.")$table
  })
  write_csv_04f(dplyr::bind_rows(sig_models), file.path(prompt_stats_dir, "signature_score_models.tsv"))
  write_csv_04f(dplyr::bind_rows(sig_bin_models), file.path(prompt_stats_dir, "signature_pseudotime_interaction_models.tsv"))
  sig_long <- sig_sample_cluster %>%
    tidyr::pivot_longer(cols = dplyr::all_of(sig_cols), names_to = "signature", values_to = "signature_score") %>%
    add_dose_x_04f("dose")
  sig_score_sig <- make_group_significance_04f(sig_long, "signature_score", dose_col = "dose", group_col = "initial_ploidy_group", facet_cols = "signature")
  write_csv_04f(sig_score_sig, file.path(prompt_stats_dir, "signature_scores_pairwise_significance.tsv"))
  p_sig_dose <- ggplot(sig_long, aes(x = .data$dose_x, y = .data$signature_score, color = .data$initial_ploidy_group, group = interaction(.data$dose, .data$initial_ploidy_group))) +
    geom_boxplot(position = position_dodge(width = 0.56), outlier.shape = NA, width = 0.48) +
    geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.55), size = 1.7) +
    facet_wrap(~signature, scales = "free_y", ncol = 2) +
    dose_scale_04f() +
    scale_color_manual(values = c("2N" = "#1D4E89", "4N" = "#A84600"), drop = FALSE) +
    labs(title = "Signature scores by dose and Initial Ploidy", x = "Dose", y = "Signature score", color = "Initial Ploidy") +
    theme_04f(7) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
  p_sig_dose <- add_significance_layers_04f(p_sig_dose, sig_score_sig, label_size = 1.7)
  save_prompt_figure(p_sig_dose, prompt_pdf_dir, prompt_png_dir, "signature_scores_by_dose_initial_ploidy", width = 10, height = 10)
  sig_heat <- sig_sample_bin %>%
    tidyr::pivot_longer(cols = dplyr::all_of(sig_cols), names_to = "signature", values_to = "signature_score") %>%
    dplyr::group_by(.data$initial_ploidy_group, .data$dose, .data$pseudotime_bin_10, .data$signature) %>%
    dplyr::summarise(mean_score = weighted.mean(.data$signature_score, w = .data$n_cells, na.rm = TRUE), .groups = "drop")
  sig_row_meta <- sig_heat %>%
    dplyr::distinct(.data$initial_ploidy_group, .data$signature, .data$dose) %>%
    dplyr::mutate(
      initial_ploidy_group = factor(as.character(.data$initial_ploidy_group), levels = c("2N", "4N")),
      signature = factor(as.character(.data$signature), levels = sig_cols),
      dose = factor(as.character(.data$dose), levels = c("0mg/kg", "30mg/kg", "120mg/kg"))
    ) %>%
    dplyr::arrange(.data$dose, .data$signature, .data$initial_ploidy_group) %>%
    dplyr::mutate(
      row_id = paste(.data$dose, .data$signature, .data$initial_ploidy_group, sep = " | "),
      dose_label = as.character(.data$dose)
    )
  sig_row_order <- sig_row_meta$row_id
  write_csv_04f(sig_row_meta, file.path(prompt_data_dir, "signature_heatmap_row_annotations.tsv"))
  sig_heat <- sig_heat %>%
    dplyr::left_join(sig_row_meta, by = c("initial_ploidy_group", "signature", "dose")) %>%
    dplyr::mutate(row_id = factor(.data$row_id, levels = rev(sig_row_order)))
  sig_annotation_ploidy <- sig_row_meta %>%
    dplyr::mutate(row_id = factor(.data$row_id, levels = rev(sig_row_order)), annotation = "Initial\nPloidy")
  sig_annotation_pathway <- sig_row_meta %>%
    dplyr::mutate(row_id = factor(.data$row_id, levels = rev(sig_row_order)), annotation = "Pathway")
  sig_annotation_dose <- sig_row_meta %>%
    dplyr::mutate(row_id = factor(.data$row_id, levels = rev(sig_row_order)), annotation = "Dose")
  pathway_colors <- setNames(grDevices::hcl.colors(length(sig_cols), palette = "Dark 3"), sig_cols)
  dose_colors <- c("0mg/kg" = "#4C78A8", "30mg/kg" = "#F58518", "120mg/kg" = "#54A24B")
  ploidy_annotation_colors <- c("2N" = "#A6CEE3", "4N" = "#FDBF6F")
  p_sig_annot_ploidy <- ggplot(sig_annotation_ploidy, aes(x = .data$annotation, y = .data$row_id, fill = .data$initial_ploidy_group)) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_manual(values = ploidy_annotation_colors, drop = FALSE, name = "Initial Ploidy") +
    scale_y_discrete(limits = rev(sig_row_order), expand = c(0, 0)) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_04f(6) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", plot.margin = margin(5.5, 1, 5.5, 5.5))
  p_sig_annot_pathway <- ggplot(sig_annotation_pathway, aes(x = .data$annotation, y = .data$row_id, fill = .data$signature)) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_manual(values = pathway_colors, drop = FALSE, name = "Pathway", guide = guide_legend(ncol = 1, byrow = TRUE)) +
    scale_y_discrete(limits = rev(sig_row_order), expand = c(0, 0)) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_04f(6) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", plot.margin = margin(5.5, 1, 5.5, 1))
  p_sig_annot_dose <- ggplot(sig_annotation_dose, aes(x = .data$annotation, y = .data$row_id, fill = .data$dose)) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_manual(values = dose_colors, drop = FALSE, name = "Dose") +
    scale_y_discrete(limits = rev(sig_row_order), expand = c(0, 0)) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_04f(6) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", plot.margin = margin(5.5, 1, 5.5, 1))
  p_sig_heat_main <- ggplot(sig_heat, aes(x = .data$pseudotime_bin_10, y = .data$row_id, fill = .data$mean_score)) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_gradient2(low = "#313695", mid = "white", high = "#A50026", midpoint = 0, name = "Score") +
    scale_y_discrete(limits = rev(sig_row_order), labels = setNames(sig_row_meta$dose_label, as.character(sig_row_meta$row_id)), expand = c(0, 0)) +
    labs(title = "Signature scores over pseudotime", x = "Pseudotime bin", y = NULL) +
    theme_04f(6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom", plot.margin = margin(5.5, 5.5, 5.5, 1))
  p_sig_heat <- p_sig_annot_dose + p_sig_annot_pathway + p_sig_annot_ploidy + p_sig_heat_main +
    patchwork::plot_layout(widths = c(0.42, 0.5, 0.42, 8.2), guides = "collect")
  save_prompt_figure(p_sig_heat, prompt_pdf_dir, prompt_png_dir, "fig11_signature_scores_heatmap", width = 12, height = 12)
  save_prompt_figure(p_sig_heat, prompt_pdf_dir, prompt_png_dir, "signature_scores_pseudotime_heatmap", width = 12, height = 12)
  save_prompt_figure(p_sig_heat, prompt_pdf_dir, prompt_png_dir, "signature_delta_vs_control_heatmap", width = 12, height = 12)

  gsea_nperm <- as.integer(get_env_scalar("GSEA_NPERM", "1000"))
  gsea_min_size <- as.integer(get_env_scalar("GSEA_MIN_SIZE", "3"))
  gsea_exponent <- safe_num(get_env_scalar("GSEA_EXPONENT", "1"))
  gsea_input <- tryCatch(
    read_h5ad_group_expression_summaries(h5ad_file, target_meta, c("initial_ploidy_group", "dose", "pseudotime_bin_10")),
    error = function(e) list(status = "error", reason = conditionMessage(e), group_sums = NULL, group_meta = data.frame(), gene_names = character(0), source_info = data.frame())
  )
  gsea_status <- data.frame(
    status = gsea_input$status,
    reason = if (!is.null(gsea_input$reason)) gsea_input$reason else NA_character_,
    nperm = gsea_nperm,
    min_size = gsea_min_size,
    exponent = gsea_exponent,
    ranking = "All h5ad X genes ranked by mean expression in group minus mean expression in other target Tumor cells.",
    stringsAsFactors = FALSE
  )
  if (nrow(gsea_input$source_info) > 0) {
    write_csv_04f(gsea_input$source_info, file.path(prompt_logs_dir, "gsea_expression_source_info.tsv"))
  }
  if (identical(gsea_input$status, "ok")) {
    gsea_nes <- compute_preranked_gsea_04f(
      gsea_input,
      signatures,
      nperm = gsea_nperm,
      min_size = gsea_min_size,
      exponent = gsea_exponent,
      seed = as.integer(get_env_scalar("GSEA_SEED", "1234"))
    )
    write_csv_04f(gsea_nes, file.path(prompt_data_dir, "gsea_nes_by_group_pseudotime_bin.tsv"))
    gsea_status$n_results <- nrow(gsea_nes)
    if (nrow(gsea_nes) > 0) {
      gsea_heat <- gsea_nes %>%
        dplyr::mutate(
          initial_ploidy_group = factor(as.character(.data$initial_ploidy_group), levels = c("2N", "4N")),
          signature = factor(as.character(.data$signature), levels = sig_cols),
          dose = factor(as.character(.data$dose), levels = c("0mg/kg", "30mg/kg", "120mg/kg"))
        )
      gsea_row_meta <- gsea_heat %>%
        dplyr::distinct(.data$initial_ploidy_group, .data$signature, .data$dose) %>%
        dplyr::arrange(.data$dose, .data$signature, .data$initial_ploidy_group) %>%
        dplyr::mutate(row_id = paste(.data$dose, .data$signature, .data$initial_ploidy_group, sep = " | "))
      gsea_row_order <- gsea_row_meta$row_id
      write_csv_04f(gsea_row_meta, file.path(prompt_data_dir, "gsea_nes_heatmap_row_annotations.tsv"))
      gsea_heat <- gsea_heat %>%
        dplyr::left_join(gsea_row_meta, by = c("initial_ploidy_group", "signature", "dose")) %>%
        dplyr::mutate(row_id = factor(.data$row_id, levels = rev(gsea_row_order)))
      gsea_annot_dose <- gsea_row_meta %>%
        dplyr::mutate(row_id = factor(.data$row_id, levels = rev(gsea_row_order)), annotation = "Dose")
      gsea_annot_pathway <- gsea_row_meta %>%
        dplyr::mutate(row_id = factor(.data$row_id, levels = rev(gsea_row_order)), annotation = "Pathway")
      gsea_annot_ploidy <- gsea_row_meta %>%
        dplyr::mutate(row_id = factor(.data$row_id, levels = rev(gsea_row_order)), annotation = "Initial\nPloidy")
      p_gsea_annot_dose <- ggplot(gsea_annot_dose, aes(x = .data$annotation, y = .data$row_id, fill = .data$dose)) +
        geom_tile(color = "white", linewidth = 0.15) +
        scale_fill_manual(values = dose_colors, drop = FALSE, name = "Dose") +
        scale_y_discrete(limits = rev(gsea_row_order), expand = c(0, 0)) +
        scale_x_discrete(position = "top") +
        labs(x = NULL, y = NULL) +
        theme_04f(6) +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", plot.margin = margin(5.5, 1, 5.5, 5.5))
      p_gsea_annot_pathway <- ggplot(gsea_annot_pathway, aes(x = .data$annotation, y = .data$row_id, fill = .data$signature)) +
        geom_tile(color = "white", linewidth = 0.15) +
        scale_fill_manual(values = pathway_colors, drop = FALSE, name = "Pathway", guide = guide_legend(ncol = 1, byrow = TRUE)) +
        scale_y_discrete(limits = rev(gsea_row_order), expand = c(0, 0)) +
        scale_x_discrete(position = "top") +
        labs(x = NULL, y = NULL) +
        theme_04f(6) +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", plot.margin = margin(5.5, 1, 5.5, 1))
      p_gsea_annot_ploidy <- ggplot(gsea_annot_ploidy, aes(x = .data$annotation, y = .data$row_id, fill = .data$initial_ploidy_group)) +
        geom_tile(color = "white", linewidth = 0.15) +
        scale_fill_manual(values = ploidy_annotation_colors, drop = FALSE, name = "Initial Ploidy") +
        scale_y_discrete(limits = rev(gsea_row_order), expand = c(0, 0)) +
        scale_x_discrete(position = "top") +
        labs(x = NULL, y = NULL) +
        theme_04f(6) +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", plot.margin = margin(5.5, 1, 5.5, 1))
      nes_limit <- max(abs(gsea_heat$NES), na.rm = TRUE)
      if (!is.finite(nes_limit) || nes_limit <= 0) nes_limit <- 1
      p_gsea_heat_main <- ggplot(gsea_heat, aes(x = .data$pseudotime_bin_10, y = .data$row_id, fill = .data$NES)) +
        geom_tile(color = "white", linewidth = 0.15) +
        scale_fill_gradient2(low = "#313695", mid = "white", high = "#A50026", midpoint = 0, limits = c(-nes_limit, nes_limit), name = "NES") +
        scale_y_discrete(limits = rev(gsea_row_order), expand = c(0, 0)) +
        labs(title = "Preranked GSEA NES over pseudotime", subtitle = "Ranks use all h5ad X genes: group mean expression minus other target Tumor cells.", x = "Pseudotime bin", y = NULL) +
        theme_04f(6) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom", plot.margin = margin(5.5, 5.5, 5.5, 1))
      p_gsea_heat <- p_gsea_annot_dose + p_gsea_annot_pathway + p_gsea_annot_ploidy + p_gsea_heat_main +
        patchwork::plot_layout(widths = c(0.42, 0.5, 0.42, 8.2), guides = "collect")
      save_prompt_figure(p_gsea_heat, prompt_pdf_dir, prompt_png_dir, "fig12_gsea_nes_heatmap", width = 12, height = 12)
      save_prompt_figure(p_gsea_heat, prompt_pdf_dir, prompt_png_dir, "gsea_nes_heatmap", width = 12, height = 12)
    }
  } else {
    write_csv_04f(data.frame(), file.path(prompt_data_dir, "gsea_nes_by_group_pseudotime_bin.tsv"))
  }
  write_csv_04f(gsea_status, file.path(prompt_logs_dir, "gsea_analysis_status.tsv"))

  expr <- h5res$expr
  group_cluster <- signature_scores_joined$sample_id[match(rownames(expr), signature_scores_joined$cell_id)]
  group_cluster <- paste(
    group_cluster,
    signature_scores_joined$initial_ploidy_group[match(rownames(expr), signature_scores_joined$cell_id)],
    signature_scores_joined$dose[match(rownames(expr), signature_scores_joined$cell_id)],
    signature_scores_joined$cluster_target[match(rownames(expr), signature_scores_joined$cell_id)],
    sep = "|"
  )
  pseudo_cluster <- rowsum(expr, group = group_cluster, reorder = FALSE)
  write_csv_04f(data.frame(group_id = rownames(pseudo_cluster), pseudo_cluster, check.names = FALSE), file.path(prompt_data_dir, "pseudobulk_counts_by_cluster.tsv"))
  meta_cluster <- data.frame(group_id = rownames(pseudo_cluster), stringsAsFactors = FALSE) %>%
    tidyr::separate(.data$group_id, into = c("sample_id", "initial_ploidy_group", "dose", "cluster_target"), sep = "\\|", remove = FALSE)
  write_csv_04f(meta_cluster, file.path(prompt_data_dir, "pseudobulk_metadata_by_cluster.tsv"))
  group_bin <- paste(
    signature_scores_joined$sample_id[match(rownames(expr), signature_scores_joined$cell_id)],
    signature_scores_joined$initial_ploidy_group[match(rownames(expr), signature_scores_joined$cell_id)],
    signature_scores_joined$dose[match(rownames(expr), signature_scores_joined$cell_id)],
    signature_scores_joined$pseudotime_bin_10[match(rownames(expr), signature_scores_joined$cell_id)],
    sep = "|"
  )
  pseudo_bin <- rowsum(expr, group = group_bin, reorder = FALSE)
  write_csv_04f(data.frame(group_id = rownames(pseudo_bin), pseudo_bin, check.names = FALSE), file.path(prompt_data_dir, "pseudobulk_counts_by_pseudotime_bin.tsv"))
  meta_bin <- data.frame(group_id = rownames(pseudo_bin), stringsAsFactors = FALSE) %>%
    tidyr::separate(.data$group_id, into = c("sample_id", "initial_ploidy_group", "dose", "pseudotime_bin_10"), sep = "\\|", remove = FALSE)
  write_csv_04f(meta_bin, file.path(prompt_data_dir, "pseudobulk_metadata_by_pseudotime_bin.tsv"))
  de_summary <- data.frame(
    analysis = "pseudobulk_DE",
    status = "skipped_deseq2",
    reason = "h5ad X does not appear integer-like raw counts; signature-gene pseudobulk matrices were saved and signature-score models were run.",
    stringsAsFactors = FALSE
  )
  ensure_dir_04f(file.path(prompt_stats_dir, "pseudobulk_DE_results"))
  write_csv_04f(de_summary, file.path(prompt_stats_dir, "pseudobulk_DE_summary.tsv"))
} else {
  write_csv_04f(signature_status, file.path(prompt_logs_dir, "signature_analysis_status.tsv"))
  write_csv_04f(data.frame(), file.path(prompt_data_dir, "signature_gene_matching.tsv"))
  write_csv_04f(data.frame(), file.path(prompt_data_dir, "cell_signature_scores.tsv"))
  write_csv_04f(gsea_status, file.path(prompt_logs_dir, "gsea_analysis_status.tsv"))
  write_csv_04f(data.frame(), file.path(prompt_data_dir, "gsea_nes_by_group_pseudotime_bin.tsv"))
  write_csv_04f(data.frame(analysis = "pseudobulk_DE", status = "skipped", reason = h5res$status), file.path(prompt_stats_dir, "pseudobulk_DE_summary.tsv"))
}
write_csv_04f(signature_status, file.path(prompt_logs_dir, "signature_analysis_status.tsv"))

cellrank_status <- data.frame(
  analysis = "CellRank",
  status = "skipped",
  reason = "CellRank is not available in the current R/Python environment; pseudotime-bin and abundance analyses are the primary trajectory-level analyses.",
  stringsAsFactors = FALSE
)
write_csv_04f(cellrank_status, file.path(prompt_logs_dir, "cellrank_status.tsv"))

required_figures <- c(
  "fig1_umap_pseudotime_target_cells.pdf",
  "fig1c_umap_cell_ploidy_target_cells_by_dose.pdf",
  "fig2_umap_clusters_target_cells.pdf",
  "fig3_umap_cell_ploidy_by_dose.pdf",
  "fig3b_umap_pseudotime_by_dose.pdf",
  "fig4_pseudotime_distribution_per_sample.pdf",
  "fig5_sample_pseudotime_metrics.pdf",
  "fig6_ploidy_violin_box_sample_overlay.pdf",
  "fig7_ploidy_vs_pseudotime_smooth.pdf",
  "fig8_delta_ploidy_heatmap.pdf",
  "fig9_composition_within_state_decomposition.pdf",
  "fig10_cluster_and_pseudotime_bin_abundance.pdf",
  "fig11_signature_scores_heatmap.pdf",
  "fig12_gsea_nes_heatmap.pdf"
)

model_tables_for_summary <- dplyr::bind_rows(
  pseudotime_metric_models,
  dplyr::bind_rows(ploidy_metric_models),
  cell_level_ploidy_model$table,
  spline_fit$table,
  bin_fit$table,
  delta_fit$table,
  dplyr::bind_rows(cluster_models),
  bin_abundance_fit$table,
  sensitivity_summary
)
write_csv_04f(model_tables_for_summary, file.path(prompt_stats_dir, "all_model_coefficients.tsv"))

key_results <- data.frame(
  question = c(
    "Do 2N and 4N Tumor cells occupy different pseudotime positions?",
    "Does Gemcitabine shift the pseudotime distribution?",
    "Does Cell Ploidy change with Gemcitabine dose?",
    "Are ploidy changes composition-driven or within-state?",
    "Do 2N and 4N cells show different dose-response patterns?",
    "At which pseudotime stages are strongest ploidy changes observed?",
    "Which signatures are associated with response?"
  ),
  main_result = c(
    "See sample pseudotime metric models and per-sample distributions.",
    "Dose-specific shifts are summarized by sample pseudotime metrics, Wasserstein distances, and bin abundance models.",
    "Cell Ploidy response is summarized by sample-level ploidy metrics and sample-bin ploidy models.",
    "Composition and within-state effects were decomposed relative to 0mg/kg controls.",
    "Initial Ploidy x Dose interactions were tested in sample-level, bin-level, and delta-ploidy models.",
    "Pseudotime-bin models and delta heatmaps localize effects across ten bins and early/middle/late sensitivity bins.",
    ifelse(identical(h5res$status, "ok"), "Signature scores and preranked GSEA NES were computed from h5ad X; pseudobulk DESeq2 was skipped because X was not integer-like.", "Signature and GSEA analyses were skipped; see logs/signature_analysis_status.tsv and logs/gsea_analysis_status.tsv.")
  ),
  supporting_figure = c(
    "figures/pdf/fig4_pseudotime_distribution_per_sample.pdf",
    "figures/pdf/fig5_sample_pseudotime_metrics.pdf",
    "figures/pdf/fig6_ploidy_violin_box_sample_overlay.pdf",
    "figures/pdf/fig9_composition_within_state_decomposition.pdf",
    "figures/pdf/fig7_ploidy_vs_pseudotime_smooth.pdf",
    "figures/pdf/fig8_delta_ploidy_heatmap.pdf",
    "figures/pdf/fig11_signature_scores_heatmap.pdf; figures/pdf/fig12_gsea_nes_heatmap.pdf"
  ),
  supporting_table = c(
    "stats/pseudotime_metric_models.tsv",
    "data/sample_pseudotime_wasserstein.tsv",
    "stats/ploidy_sample_metric_models.tsv",
    "data/composition_within_state_decomposition.tsv",
    "stats/ploidy_pseudotime_bin_interaction_model.tsv",
    "stats/delta_ploidy_models.tsv",
    "stats/signature_score_models.tsv; data/gsea_nes_by_group_pseudotime_bin.tsv"
  ),
  model = c(
    "Sample-level OLS HC3",
    "Sample-level OLS HC3 and Wasserstein distance",
    "Sample-level OLS HC3; exploratory cell-level robust model",
    "Pseudotime-bin decomposition with sample bootstrap",
    "WLS sample-bin interaction models",
    "WLS sample-bin interaction and early/middle/late sensitivity",
    ifelse(identical(h5res$status, "ok"), "WLS signature models; preranked GSEA", "Not run")
  ),
  effect_size = NA_character_,
  p_value = NA_real_,
  q_value = NA_real_,
  interpretation = c(
    "Use sample-aware estimates rather than cell-level p-values as primary evidence.",
    "Non-monotonic dose responses are allowed because dose is categorical.",
    "Initial Ploidy groups should be interpreted separately when interactions are present.",
    "Composition effect suggests trajectory-state redistribution; within-state effect suggests ploidy change after pseudotime matching.",
    "Interaction terms quantify 4N-specific responses relative to 2N.",
    "Bins with largest delta values identify stage-specific response.",
    "Signature results are exploratory unless raw count DE is available."
  ),
  stringsAsFactors = FALSE
)
write_csv_04f(key_results, file.path(prompt_report_dir, "key_results_summary.tsv"))
write_csv_04f(key_results, file.path(prompt_report_dir, "analysis_summary.tsv"))

report_lines <- c(
  "# Downstream scVelo Tumor Non-Cell-Cycle Trajectory Analysis",
  "",
  "## 1. Analysis Object and Input Files",
  "",
  paste0("Selected input file: `", input_file, "`."),
  paste0("Raw scVelo all-cells input root inspected: `", raw_scvelo_input_root, "`."),
  paste0("The selected table is the project-confirmed 04e merged scVelo all-cells output with Cell Ploidy annotations. Tumor cells were extracted from this all-cells table and then restricted to all clusters except ", excluded_cluster_label, "."),
  "",
  paste0("Target cells: ", nrow(target_meta), ". Target samples: ", dplyr::n_distinct(target_meta$sample_id), "."),
  "",
  "## 2. QC Results",
  "",
  paste0("Expected 16 samples found: ", dplyr::n_distinct(target_meta$sample_id) == 16, "."),
  paste0("Dose groups found: ", paste(levels(droplevels(target_meta$dose)), collapse = ", "), "."),
  paste0("Initial Ploidy groups found: ", paste(levels(droplevels(target_meta$initial_ploidy_group)), collapse = ", "), "."),
  "Missingness is saved in `data/missingness_summary.tsv`. Counts are saved in `data/counts_by_*`. QC figures are saved under `figures/pdf` and `figures/png`.",
  "",
  "## 3. Pseudotime Distribution Analysis",
  "",
  "Sample-level pseudotime metrics, entropy, peak bin, early/middle/late fractions, cluster fractions, and Wasserstein distances to matched 0mg/kg controls were computed. Models treat dose as categorical and include Initial Ploidy x Dose interactions.",
  "",
  "Key files: `data/sample_pseudotime_metrics.tsv`, `data/sample_pseudotime_wasserstein.tsv`, and `stats/pseudotime_metric_models.tsv`.",
  "",
  "## 4. Overall Ploidy Response",
  "",
  "Cell Ploidy metrics were summarized per sample, including mean, median, tail fractions, and relative Cell Ploidy sensitivity metrics. Cell-level models are explicitly labeled exploratory; sample-level summaries and sample-bin models are primary.",
  "",
  "Key files: `data/sample_ploidy_metrics.tsv`, `stats/ploidy_sample_metric_models.tsv`, and `stats/ploidy_cell_level_model_coefficients.tsv`.",
  "",
  "## 5. Joint Ploidy-Pseudotime Analysis",
  "",
  "A sample-bin WLS spline model and an interpretable ten-bin interaction model were fit for Cell Ploidy over pseudotime. Models include Initial Ploidy x Dose x pseudotime terms and cluster-robust standard errors by sample.",
  "",
  "Key files: `data/pseudotime_bin_ploidy_by_sample.tsv`, `stats/ploidy_pseudotime_spline_model_coefficients.tsv`, and `stats/ploidy_pseudotime_bin_interaction_model.tsv`.",
  "",
  "## 6. Composition vs Within-State Decomposition",
  "",
  "Ploidy shifts relative to 0mg/kg controls were decomposed into total, pseudotime-composition, and within-state components using ten pseudotime bins. Bootstrap confidence intervals were estimated by sample resampling.",
  "",
  "Key files: `data/composition_within_state_decomposition.tsv` and `data/composition_within_state_decomposition_bootstrap.tsv`.",
  "",
  "## 7. Cluster and Pseudotime-Bin Abundance",
  "",
  "Cluster fractions and pseudotime-bin fractions were modeled at sample level using logit-transformed fractions with pseudocount correction.",
  "",
  "Key files: `stats/cluster_abundance_models.tsv` and `stats/pseudotime_bin_abundance_models.tsv`.",
  "",
  "## 8. Signature and Pathway Analysis",
  "",
  ifelse(
    identical(h5res$status, "ok"),
    "Signature scores were computed from h5ad `X` for matched genes in Gemcitabine metabolism, S phase/DNA replication, replication stress, DNA damage response, G2/M checkpoint, apoptosis, and senescence/p53 response signatures. Preranked GSEA NES was also computed over pseudotime bins using all h5ad `X` genes ranked by group mean expression minus other target Tumor cells. h5ad `X` did not appear integer-like, so DESeq2 pseudobulk DE was skipped and signature-score/GSEA summaries are the primary expression-level analysis.",
    paste0("Signature/GSEA analysis was skipped. Status: ", h5res$status, ". See `logs/signature_analysis_status.tsv` and `logs/gsea_analysis_status.tsv`.")
  ),
  "",
  "Key GSEA files: `data/gsea_nes_by_group_pseudotime_bin.tsv`, `data/gsea_nes_heatmap_row_annotations.tsv`, and `figures/pdf/fig12_gsea_nes_heatmap.pdf`.",
  "",
  "## 9. Main Conclusions",
  "",
  "Use `report/key_results_summary.tsv` as the concise index linking each biological question to supporting figures, tables, and model outputs. Because the sample size is small, conclusions should prioritize sample-aware effect sizes and confidence intervals over cell-level p-values.",
  "",
  "## 10. Notes and Limitations",
  "",
  "- Pseudotime is an inferred trajectory coordinate, not chronological time.",
  "- Cell Ploidy may not equal absolute DNA content unless experimentally calibrated.",
  "- Cell-level p-values are exploratory because cells within a sample are not independent.",
  "- Validation is recommended by DNA-content flow cytometry/DAPI/PI, EdU/BrdU, pH3, gamma-H2AX, Ki67, and optionally scDNA/CNV/FISH."
)
write_text_04f(report_lines, file.path(prompt_report_dir, "analysis_report.md"))

readme_lines <- c(
  "# 04g non-cell-cycle downstream results",
  "",
  paste0("This directory contains the full prompt-driven downstream analysis for scVelo all-cells-derived Tumor cells in all clusters except ", excluded_cluster_label, "."),
  "",
  "The selected input is the 04e merged scVelo all_cells table with Cell Ploidy annotations. Outputs are organized under data, figures, stats, models, logs, report, and scripts."
)
write_text_04f(readme_lines, file.path(output_base, "README.md"))
write_text_04f(
  c(
    "Rscript Code/in-vivo/04g_nonCellCycle_results.R",
    "",
    "Required R packages: dplyr, ggplot2, readr, tidyr, tibble, patchwork, sandwich, lmtest, jsonlite, hdf5r."
  ),
  file.path(output_base, "requirements.txt")
)
write_text_04f(
  c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "Rscript Code/in-vivo/04g_nonCellCycle_results.R"
  ),
  file.path(prompt_scripts_dir, "run_all.sh")
)
write_text_04f(capture.output(sessionInfo()), file.path(prompt_logs_dir, "session_info.txt"))
write_text_04f(
  c(
    paste0("Completed at: ", as.character(Sys.time())),
    paste0("Output directory: ", output_base),
    paste0("Target cells: ", nrow(target_meta)),
    paste0("Target samples: ", dplyr::n_distinct(target_meta$sample_id)),
    "Key figures:",
    paste0("  figures/pdf/", required_figures),
    "Key model tables:",
    "  stats/pseudotime_metric_models.tsv",
    "  stats/ploidy_sample_metric_models.tsv",
    "  stats/ploidy_pseudotime_spline_model_coefficients.tsv",
    "  stats/delta_ploidy_models.tsv"
  ),
  file.path(prompt_logs_dir, "pipeline.log")
)

status <- data.frame(
  output_group = c("histograms", "density", "dose_histograms", "dose_density", "dose_response", "full_prompt_outputs"),
  status = "ok",
  n_target_cells = nrow(target_df),
  n_samples = length(unique(target_df$sampleID)),
  stringsAsFactors = FALSE
)
write_csv_04f(status, file.path(output_base, "04g_nonCellCycle_results_status.csv"))

message("04g non-cell-cycle results completed.")
message("Output root: ", output_base)
message("Target Tumor cells excluding clusters ", excluded_cluster_label, ": ", nrow(target_df))
message("Target samples: ", length(unique(target_df$sampleID)))
