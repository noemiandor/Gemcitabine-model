#!/usr/bin/env Rscript

# Plot-only generator for Supplementary Figures 4-7.
#
# The 11 committed tables under Data/in-vivo/SIfigures are the reviewed,
# plot-facing contract. This script intentionally cannot load raw Seurat data,
# run differential expression, or recompute enrichment. The canonical default
# accepts only the reviewed SI7 cache. A separately labelled generated
# human-only raw rebuild can be rendered only with an explicit opt-in and is
# never suitable for publication over Data/in-vivo/SIfigures until reviewed.

parse_args <- function(tokens) {
  result <- list()
  i <- 1L
  while (i <= length(tokens)) {
    token <- tokens[[i]]
    if (grepl("^--[^=]+=", token)) {
      key <- sub("^--([^=]+)=.*$", "\\1", token)
      result[[key]] <- sub("^--[^=]+=", "", token)
      i <- i + 1L
    } else if (grepl("^--", token)) {
      key <- sub("^--", "", token)
      if (i < length(tokens) && !grepl("^--", tokens[[i + 1L]])) {
        result[[key]] <- tokens[[i + 1L]]
        i <- i + 2L
      } else {
        result[[key]] <- "TRUE"
        i <- i + 1L
      }
    } else {
      stop("Unexpected positional argument: ", token, call. = FALSE)
    }
  }
  result
}

arg_value <- function(args, name, default = NULL) {
  value <- args[[name]]
  if (!is.null(value) && length(value) == 1L && nzchar(value)) value else default
}

arg_flag <- function(args, name, default = FALSE) {
  value <- arg_value(args, name, NULL)
  if (is.null(value)) return(default)
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript Code/in-vivo/SI_figures/generate_supplementary_figures.R \\",
      "    [--table-cache-dir Data/in-vivo/SIfigures] \\",
      "    [--config Code/in-vivo/figure7/figure7_config.yaml] \\",
      "    --output-dir Results/in-vivo/SI_figures/runs/RUN_si_figures",
      "",
      "Options:",
      "  --overwrite  Replace this generator's known outputs.",
      paste(
        "  --allow-generated-human-only-si7",
        "Render an unreviewed run-scoped cache rebuilt with the GRCh-only policy."
      ),
      "  --help       Show this help.",
      "",
      "This entrypoint is deliberately plot-only and emits four composite",
      "Supplementary Figures (4-7) as PDF/PNG pairs.",
      sep = "\n"
    ),
    "\n"
  )
}

script_path <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  match <- grep("^--file=", command, value = TRUE)
  if (!length(match)) return(NA_character_)
  normalizePath(sub("^--file=", "", match[[1L]]), mustWork = TRUE)
}

resolve_path <- function(path, root, must_work = FALSE) {
  candidate <- if (grepl("^/", path)) path else file.path(root, path)
  normalizePath(candidate, mustWork = must_work)
}

repo_relative <- function(path, root) {
  normalized <- normalizePath(path, mustWork = TRUE)
  prefix <- paste0(normalizePath(root, mustWork = TRUE), .Platform$file.sep)
  if (!startsWith(normalized, prefix)) {
    return(paste0("external:", basename(normalized)))
  }
  substring(normalized, nchar(prefix) + 1L)
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Missing required R package: ", package, call. = FALSE)
  }
}

read_table <- function(path, delimiter, label) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing or empty ", label, ": ", path, call. = FALSE)
  }
  reader <- if (delimiter == "\t") utils::read.delim else utils::read.csv
  data <- reader(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = if (delimiter == "\t") "" else "\"",
    comment.char = "",
    na.strings = c("", "NA", "NaN")
  )
  if (!nrow(data)) stop(label, " has no data rows", call. = FALSE)
  data
}

read_csv <- function(path, label) read_table(path, ",", label)
read_tsv <- function(path, label) read_table(path, "\t", label)

file_sha256 <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(file = path, algo = "sha256", serialize = FALSE))
  }
  executable <- Sys.which("shasum")
  if (!nzchar(executable)) stop("Need digest or shasum for SHA-256", call. = FALSE)
  result <- system2(executable, c("-a", "256", shQuote(path)), stdout = TRUE)
  strsplit(trimws(result[[1L]]), "[[:space:]]+")[[1L]][[1L]]
}

write_tsv <- function(data, path) {
  utils::write.table(
    data,
    path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = ""
  )
}

figure_theme <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", color = "#222222", size = base_size + 1
      ),
      plot.subtitle = ggplot2::element_text(
        color = "#444444", size = base_size - 1
      ),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 3),
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(
        fill = "grey94", color = "grey75", linewidth = 0.35
      ),
      strip.text = ggplot2::element_text(face = "bold", color = "#333333"),
      axis.title = ggplot2::element_text(color = "#333333")
    )
}

umap_theme <- function(base_size = 10) {
  figure_theme(base_size) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(color = "grey35", linewidth = 0.35),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      legend.text = ggplot2::element_text(size = base_size + 1),
      legend.key.height = grid::unit(0.85, "lines"),
      legend.key.width = grid::unit(0.90, "lines")
    )
}

add_tag <- function(plot, tag) plot + ggplot2::labs(tag = tag)

shuffle_cells <- function(data, seed) {
  set.seed(seed)
  data[sample.int(nrow(data)), , drop = FALSE]
}

make_umap_discrete <- function(
  data,
  field,
  colors,
  title,
  legend_title,
  tag,
  point_size,
  subtitle = NULL,
  labels = FALSE
) {
  data <- shuffle_cells(data, plot_seed + utf8ToInt(tag)[[1L]])
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(UMAP_1, UMAP_2, color = .data[[field]])
  ) +
    ggplot2::geom_point(size = point_size, alpha = 0.76, stroke = 0) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE, name = legend_title) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 3, alpha = 1, stroke = 0)
      )
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2") +
    umap_theme(10)
  if (labels) {
    centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, data = data, FUN = median)
    plot <- plot + ggplot2::geom_label(
      data = centers,
      ggplot2::aes(UMAP_1, UMAP_2, label = cluster),
      inherit.aes = FALSE,
      size = 2.5,
      linewidth = 0.2,
      fill = "white",
      color = "#222222",
      alpha = 0.86,
      label.padding = grid::unit(0.10, "lines")
    )
  }
  add_tag(plot, tag)
}

make_umap_continuous <- function(
  data,
  field,
  title,
  legend_title,
  tag,
  point_size,
  limits = NULL,
  subtitle = NULL,
  diverging = FALSE
) {
  data <- shuffle_cells(data, plot_seed + 1000L + utf8ToInt(tag)[[1L]])
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(UMAP_1, UMAP_2, color = .data[[field]])
  ) +
    ggplot2::geom_point(size = point_size, alpha = 0.80, stroke = 0) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2") +
    umap_theme(10)
  if (diverging) {
    plot <- plot + ggplot2::scale_color_gradient2(
      low = "#2C7BB6",
      mid = "white",
      high = "#D7191C",
      midpoint = 0,
      limits = limits,
      name = legend_title
    )
  } else {
    plot <- plot + ggplot2::scale_color_gradient(
      low = "#2C7BB6",
      high = "#D7191C",
      limits = limits,
      name = legend_title
    )
  }
  add_tag(plot, tag)
}

percent_labels <- function(values) paste0(round(100 * values), "%")

make_composition_plot <- function(
  data,
  fill_colors,
  title,
  x_title,
  y_title,
  tag,
  use_proportion
) {
  y_field <- if (use_proportion) "proportion" else "n_cells"
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(group_value, .data[[y_field]], fill = fill_value)
  ) +
    ggplot2::geom_col(width = 0.82, color = "white", linewidth = 0.12) +
    ggplot2::scale_fill_manual(values = fill_colors, drop = FALSE) +
    ggplot2::labs(title = title, x = x_title, y = y_title, fill = NULL) +
    figure_theme(9.5) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  if (use_proportion) {
    plot <- plot +
      ggplot2::scale_y_continuous(
        breaks = seq(0, 1, 0.25),
        labels = percent_labels,
        expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1))
  } else {
    plot <- plot + ggplot2::scale_y_continuous(
      labels = function(x) format(x, big.mark = ",", scientific = FALSE),
      expand = ggplot2::expansion(mult = c(0, 0.04))
    )
  }
  add_tag(plot, tag)
}

heatmap_plot <- function(matrix_data, title, diverging, tag) {
  arguments <- list(
    mat = matrix_data,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 8,
    fontsize_col = 7,
    angle_col = 45,
    main = paste0(tag, "  ", title),
    silent = TRUE
  )
  if (diverging) {
    max_abs <- max(abs(matrix_data))
    arguments$color <- grDevices::colorRampPalette(
      c("#2C7BB6", "white", "#D7191C")
    )(101)
    arguments$breaks <- seq(
      -max_abs,
      max_abs,
      length.out = length(arguments$color) + 1L
    )
  }
  heatmap <- do.call(pheatmap::pheatmap, arguments)
  patchwork::wrap_elements(full = heatmap$gtable)
}

save_composite <- function(
  plot,
  figure_dir,
  figure_id,
  filename_stub,
  width,
  height
) {
  pdf_path <- file.path(figure_dir, paste0(filename_stub, ".pdf"))
  png_path <- file.path(figure_dir, paste0(filename_stub, ".png"))
  pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf"
  ggplot2::ggsave(
    pdf_path,
    plot = plot,
    device = pdf_device,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
  ggplot2::ggsave(
    png_path,
    plot = plot,
    device = "png",
    dpi = 300,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
  data.frame(
    panel_id = c(figure_id, paste0(figure_id, "_png")),
    filename = basename(c(pdf_path, png_path)),
    variant = c("pdf", "png"),
    stringsAsFactors = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (arg_flag(args, "help")) {
  usage()
  quit(save = "no", status = 0L)
}

for (package in c("ggplot2", "patchwork", "pheatmap", "yaml")) {
  require_package(package)
}

script <- script_path()
if (is.na(script)) stop("Cannot resolve script path", call. = FALSE)
root <- normalizePath(
  dirname(dirname(dirname(dirname(script)))),
  mustWork = TRUE
)
cache_dir <- resolve_path(
  arg_value(args, "table-cache-dir", "Data/in-vivo/SIfigures"),
  root,
  must_work = TRUE
)
config_path <- resolve_path(
  arg_value(args, "config", "Code/in-vivo/figure7/figure7_config.yaml"),
  root,
  must_work = TRUE
)
output_dir <- resolve_path(
  arg_value(args, "output-dir", "Results/in-vivo/SI_figures"),
  root,
  must_work = FALSE
)
overwrite <- arg_flag(args, "overwrite")
allow_generated_human_only_si7 <- arg_flag(
  args,
  "allow-generated-human-only-si7"
)
upstream_manifest_arg <- arg_value(args, "upstream-input-manifest", NULL)
upstream_manifest_path <- if (is.null(upstream_manifest_arg)) {
  NA_character_
} else {
  resolve_path(upstream_manifest_arg, root, must_work = TRUE)
}
if (!is.na(upstream_manifest_path) && !allow_generated_human_only_si7) {
  stop(
    paste(
      "--upstream-input-manifest is valid only for a generated",
      "human-only raw-table render"
    ),
    call. = FALSE
  )
}
if (allow_generated_human_only_si7 && is.na(upstream_manifest_path)) {
  stop(
    paste(
      "--upstream-input-manifest is required for a generated",
      "human-only raw-table render"
    ),
    call. = FALSE
  )
}

validator <- file.path(root, "Code", "tools", "validate_si_figures_table_cache.py")
validator_args <- c(
  shQuote(validator),
  "--cache-dir",
  shQuote(cache_dir)
)
if (allow_generated_human_only_si7) {
  validator_args <- c(validator_args, "--si7-policy", "generated-human-only")
}
validation_status <- system2(
  "python3",
  validator_args
)
if (!identical(validation_status, 0L)) {
  stop("The SI Figures table cache failed validation", call. = FALSE)
}

config <- yaml::read_yaml(config_path)
si_config <- config$si_figures
if (is.null(si_config)) {
  stop("Figure 7 config is missing the si_figures contract", call. = FALSE)
}
cluster_levels <- as.character(unlist(si_config$cluster_order))
ploidy_levels <- as.character(unlist(si_config$ploidy_levels))
dose_levels <- as.character(unlist(si_config$dose_levels))
plot_seed <- as.integer(si_config$plot_shuffle_seed)
if (!identical(cluster_levels, c("0", "2", "4c", "5", "6", "8", "10", "13", "14"))) {
  stop("Unexpected SI Figures cluster order in config", call. = FALSE)
}
if (!is.finite(plot_seed)) stop("Invalid SI Figures plot seed", call. = FALSE)

if (dir.exists(output_dir)) {
  existing_figures <- list.files(
    file.path(output_dir, "figures"),
    all.files = TRUE,
    no.. = TRUE
  )
  existing_tables <- list.files(
    file.path(output_dir, "tables"),
    all.files = TRUE,
    no.. = TRUE
  )
  existing_metadata <- list.files(
    file.path(output_dir, "metadata"),
    all.files = TRUE,
    no.. = TRUE
  )
  if (
    (length(existing_figures) || length(existing_tables) || length(existing_metadata)) &&
      !overwrite
  ) {
    stop(
      "Output contains generated files; pass --overwrite: ",
      output_dir,
      call. = FALSE
    )
  }
}

figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
metadata_dir <- file.path(output_dir, "metadata")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
if (overwrite) {
  unlink(list.files(figure_dir, full.names = TRUE), recursive = TRUE)
  unlink(list.files(table_dir, full.names = TRUE), recursive = TRUE)
  unlink(list.files(metadata_dir, full.names = TRUE), recursive = TRUE)
}

cache_manifest <- read_tsv(
  file.path(cache_dir, "manifest.tsv"),
  "SI Figures cache manifest"
)
cache_files <- as.character(cache_manifest$filename)
if (length(cache_files) != 11L || anyDuplicated(cache_files)) {
  stop("SI Figures cache manifest must identify exactly 11 tables", call. = FALSE)
}
source_paths <- file.path(cache_dir, cache_files)
copied_paths <- file.path(table_dir, cache_files)
copied <- file.copy(source_paths, copied_paths, overwrite = TRUE, copy.mode = FALSE)
if (!all(copied)) stop("Failed to copy SI Figures tables", call. = FALSE)
if (!identical(
  unname(vapply(source_paths, file_sha256, character(1))),
  unname(vapply(copied_paths, file_sha256, character(1)))
)) {
  stop("SI Figures tables changed while being copied", call. = FALSE)
}
manifest_copied <- file.copy(
  file.path(cache_dir, "manifest.tsv"),
  file.path(table_dir, "manifest.tsv"),
  overwrite = TRUE,
  copy.mode = FALSE
)
if (!manifest_copied) stop("Failed to copy the SI Figures cache manifest", call. = FALSE)
copied_validation <- system2(
  "python3",
  c(
    shQuote(validator),
    "--cache-dir",
    shQuote(table_dir),
    if (allow_generated_human_only_si7) {
      c("--si7-policy", "generated-human-only")
    }
  )
)
if (!identical(copied_validation, 0L)) {
  stop("Copied SI Figures tables failed validation", call. = FALSE)
}

cells <- read_csv(
  file.path(table_dir, "si_figures_cell_metadata.csv"),
  "SI Figures cell metadata"
)
cluster_key <- read_tsv(
  file.path(table_dir, "si_figures_cluster_key.tsv"),
  "SI Figures cluster key"
)
endpoint_audit <- read_csv(
  file.path(table_dir, "si_figure6_endpoint_ploidy_join_audit.csv"),
  "SI Figure 6 endpoint-ploidy audit"
)
s4_context <- read_csv(
  file.path(table_dir, "si_figure4_cluster_context_composition.csv"),
  "SI Figure 4 context composition"
)
s4_ploidy <- read_csv(
  file.path(table_dir, "si_figure4_cluster_initial_ploidy_composition.csv"),
  "SI Figure 4 ploidy composition"
)
composition_table <- read_csv(
  file.path(table_dir, "si_figure5_cluster_composition_by_mouse.csv"),
  "SI Figure 5 per-mouse composition"
)
group_summary <- read_csv(
  file.path(
    table_dir,
    "si_figure5_mouse_weighted_composition_by_initial_ploidy_dose.csv"
  ),
  "SI Figure 5 mouse-weighted composition"
)
s5_dose <- read_csv(
  file.path(table_dir, "si_figure5_cluster_dose_composition.csv"),
  "SI Figure 5 dose composition"
)
s5_ploidy <- read_csv(
  file.path(table_dir, "si_figure5_cluster_initial_ploidy_composition.csv"),
  "SI Figure 5 ploidy composition"
)

cells$UMAP_1 <- as.numeric(cells$UMAP_1)
cells$UMAP_2 <- as.numeric(cells$UMAP_2)
cells$s_phase_score <- as.numeric(cells$s_phase_score)
cells$endpoint_ploidy <- suppressWarnings(as.numeric(cells$endpoint_ploidy))
cells$included_in_si_figures <- tolower(
  as.character(cells$included_in_si_figures)
) %in% c("true", "t", "1")

cluster_colors <- stats::setNames(
  as.character(cluster_key$color),
  as.character(cluster_key$cluster_id)
)
ploidy_colors <- c("2N" = "#4C78A8", "4N" = "#E45756")
dose_colors <- c(
  "0mg/kg" = "#666666",
  "30mg/kg" = "#d95f02",
  "120mg/kg" = "#1b9e77"
)
context_colors <- c("Tumor" = "#4C78A8", "CellLine" = "#F2CF5B")

seurat <- data.frame(
  cell = cells$cell_id,
  UMAP_1 = cells$UMAP_1,
  UMAP_2 = cells$UMAP_2,
  mouse = cells$sample_id,
  cluster = factor(cells$cluster_id, levels = cluster_levels),
  initial_ploidy = factor(cells$initial_ploidy, levels = ploidy_levels),
  dose = factor(cells$dose, levels = dose_levels),
  s_phase_score = cells$s_phase_score,
  endpoint_ploidy = cells$endpoint_ploidy,
  context = factor(cells$context, levels = c("Tumor", "CellLine")),
  included = cells$included_in_si_figures,
  stringsAsFactors = FALSE
)
tumor <- seurat[seurat$included, , drop = FALSE]
n_tumor <- nrow(tumor)
point_size <- if (n_tumor > 50000L) {
  0.08
} else if (n_tumor > 20000L) {
  0.14
} else {
  0.24
}

sample_metadata <- unique(
  tumor[, c("mouse", "initial_ploidy", "dose"), drop = FALSE]
)
sample_metadata <- sample_metadata[
  order(sample_metadata$initial_ploidy, sample_metadata$dose, sample_metadata$mouse),
  ,
  drop = FALSE
]
mouse_levels <- as.character(sample_metadata$mouse)
tumor$mouse <- factor(tumor$mouse, levels = mouse_levels)
n_mice <- length(mouse_levels)

factor_composition <- function(data, group_levels, fill_levels) {
  data$group_value <- factor(data$group_value, levels = group_levels)
  data$fill_value <- factor(data$fill_value, levels = fill_levels)
  data$n_cells <- as.numeric(data$n_cells)
  data$total_cells <- as.numeric(data$total_cells)
  data$proportion <- as.numeric(data$proportion)
  data
}
s4_context <- factor_composition(
  s4_context,
  cluster_levels,
  c("Tumor", "CellLine")
)
s4_ploidy <- factor_composition(s4_ploidy, cluster_levels, ploidy_levels)
s5_dose <- factor_composition(s5_dose, cluster_levels, dose_levels)
s5_ploidy <- factor_composition(s5_ploidy, cluster_levels, ploidy_levels)

message("Generating Supplementary Figure 4 composite.")
s4a <- make_umap_discrete(
  seurat,
  "cluster",
  cluster_colors,
  "UMAP by cluster",
  "Cluster",
  "A",
  point_size,
  sprintf("Tumor and CellLine cells; n = %s", format(nrow(seurat), big.mark = ",")),
  labels = TRUE
)
s4b <- make_umap_discrete(
  seurat,
  "initial_ploidy",
  ploidy_colors,
  "UMAP by initial ploidy",
  "Initial ploidy",
  "B",
  point_size,
  "Tumor and CellLine cells"
)
s4c <- make_umap_discrete(
  seurat,
  "context",
  context_colors,
  "UMAP by Tumor/CellLine context",
  "Context",
  "C",
  point_size,
  "Tumor and CellLine cells"
)
s4d <- make_umap_continuous(
  seurat,
  "s_phase_score",
  "UMAP by S phase score",
  "S phase score",
  "D",
  point_size,
  limits = range(seurat$s_phase_score),
  diverging = TRUE
)
s4e <- make_composition_plot(
  s4_context,
  context_colors,
  "Tumor and CellLine proportions by cluster",
  "Cluster",
  "Cell proportion",
  "E",
  TRUE
)
s4f <- make_composition_plot(
  s4_context,
  context_colors,
  "Tumor and CellLine cell counts by cluster",
  "Cluster",
  "Number of cells",
  "F",
  FALSE
)
s4g <- make_composition_plot(
  s4_ploidy,
  ploidy_colors,
  "Initial ploidy proportions by cluster",
  "Cluster",
  "Cell proportion",
  "G",
  TRUE
)
s4h <- make_composition_plot(
  s4_ploidy,
  ploidy_colors,
  "Initial ploidy cell counts by cluster",
  "Cluster",
  "Number of cells",
  "H",
  FALSE
)
s4 <- patchwork::wrap_plots(
  patchwork::wrap_plots(s4a, s4b, s4c, s4d, ncol = 4),
  patchwork::wrap_plots(s4e, s4f, ncol = 2),
  patchwork::wrap_plots(s4g, s4h, ncol = 2),
  ncol = 1,
  heights = c(1.0, 0.95, 0.95)
) + patchwork::plot_annotation(
  title = "Supplementary Figure 4 | Tumor and CellLine cellular landscape"
)
panel_rows <- list(
  save_composite(
    s4,
    figure_dir,
    "SuppFig4",
    "panel_SuppFig4_composite",
    20,
    15.5
  )
)

message("Generating Supplementary Figure 5 composite.")
s5a <- make_umap_discrete(
  tumor,
  "cluster",
  cluster_colors,
  "UMAP by cluster",
  "Cluster",
  "A",
  point_size,
  sprintf("Tumor cells; n = %s", format(n_tumor, big.mark = ",")),
  labels = TRUE
)
s5b <- make_umap_discrete(
  tumor,
  "initial_ploidy",
  ploidy_colors,
  "UMAP by initial tumor ploidy",
  "Initial ploidy",
  "B",
  point_size
)
s5c <- make_umap_discrete(
  tumor,
  "dose",
  dose_colors,
  "UMAP by Gemcitabine dose",
  "Gemcitabine dose",
  "C",
  point_size
)
s5d <- make_umap_continuous(
  tumor,
  "s_phase_score",
  "UMAP by S phase score",
  "S phase score",
  "D",
  point_size,
  limits = range(tumor$s_phase_score),
  diverging = TRUE
)

panel_metadata <- sample_metadata
panel_metadata$mouse_panel <- paste(
  panel_metadata$mouse,
  panel_metadata$initial_ploidy,
  panel_metadata$dose,
  sep = " | "
)
panel_lookup <- stats::setNames(
  panel_metadata$mouse_panel,
  panel_metadata$mouse
)
mouse_plot_data <- tumor
mouse_plot_data$mouse_panel <- factor(
  unname(panel_lookup[as.character(mouse_plot_data$mouse)]),
  levels = panel_metadata$mouse_panel
)
s5e <- ggplot2::ggplot(
  mouse_plot_data,
  ggplot2::aes(UMAP_1, UMAP_2, color = initial_ploidy)
) +
  ggplot2::geom_point(
    size = max(point_size * 1.35, 0.14),
    alpha = 0.88,
    stroke = 0
  ) +
  ggplot2::facet_wrap(~mouse_panel, ncol = 4, drop = FALSE) +
  ggplot2::scale_color_manual(
    values = ploidy_colors,
    drop = FALSE,
    name = "Initial ploidy"
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "UMAP by mouse of origin",
    subtitle = sprintf("%s mice; shared 2N/4N color code", n_mice),
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  umap_theme(8.5) +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_text(size = 7.2),
    panel.spacing = grid::unit(0.10, "lines"),
    plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
  )
s5e <- add_tag(s5e, "E")

composition_table$mouse <- factor(composition_table$mouse, levels = mouse_levels)
composition_table$initial_ploidy <- factor(
  composition_table$initial_ploidy,
  levels = ploidy_levels
)
composition_table$dose <- factor(composition_table$dose, levels = dose_levels)
composition_table$cluster_final <- factor(
  composition_table$cluster_final,
  levels = cluster_levels
)
composition_table$proportion <- as.numeric(composition_table$proportion)
s5f <- ggplot2::ggplot(
  composition_table,
  ggplot2::aes(mouse, proportion, fill = cluster_final)
) +
  ggplot2::geom_col(width = 0.84, color = "white", linewidth = 0.10) +
  ggplot2::facet_grid(
    . ~ initial_ploidy + dose,
    scales = "free_x",
    space = "free_x",
    drop = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = cluster_colors,
    drop = FALSE,
    name = "Cluster"
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = percent_labels,
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Cluster composition of each mouse",
    subtitle = "Each mouse sums to 100%",
    x = "Mouse",
    y = "Cluster composition"
  ) +
  figure_theme(9.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 50, hjust = 1, size = 7.5)
  )
s5f <- add_tag(s5f, "F")

group_summary$initial_ploidy <- factor(
  group_summary$initial_ploidy,
  levels = ploidy_levels
)
group_summary$dose <- factor(group_summary$dose, levels = dose_levels)
group_summary$cluster_final <- factor(
  group_summary$cluster_final,
  levels = cluster_levels
)
group_summary$n_mice <- as.integer(group_summary$n_mice)
group_summary$mean_proportion <- as.numeric(group_summary$mean_proportion)
group_n <- unique(
  group_summary[, c("initial_ploidy", "dose", "n_mice"), drop = FALSE]
)
s5g <- ggplot2::ggplot(
  group_summary,
  ggplot2::aes(dose, mean_proportion, fill = cluster_final)
) +
  ggplot2::geom_col(width = 0.72, color = "white", linewidth = 0.12) +
  ggplot2::geom_text(
    data = group_n,
    ggplot2::aes(dose, 0.985, label = paste0("n=", n_mice)),
    inherit.aes = FALSE,
    size = 2.8,
    vjust = 1
  ) +
  ggplot2::facet_wrap(~initial_ploidy, nrow = 1, drop = FALSE) +
  ggplot2::scale_fill_manual(
    values = cluster_colors,
    drop = FALSE,
    name = "Cluster"
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = percent_labels,
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Mouse-weighted composition by initial ploidy and dose",
    subtitle = paste(
      "Cluster proportions are calculated per mouse,",
      "then averaged with equal mouse weights"
    ),
    x = "Gemcitabine dose",
    y = "Mean cluster composition"
  ) +
  figure_theme(9.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
  )
s5g <- add_tag(s5g, "G")

s5h <- make_composition_plot(
  s5_dose,
  dose_colors,
  "Gemcitabine dose composition by cluster",
  "Cluster",
  "Cell proportion",
  "H",
  TRUE
)
s5i <- make_composition_plot(
  s5_ploidy,
  ploidy_colors,
  "Initial ploidy composition by cluster",
  "Cluster",
  "Cell proportion",
  "I",
  TRUE
)

s5_left <- patchwork::wrap_plots(s5a, s5b, s5c, s5d, ncol = 1)
s5_upper <- patchwork::wrap_plots(
  s5_left,
  s5e,
  ncol = 2,
  widths = c(1, 4)
)
s5 <- patchwork::wrap_plots(
  s5_upper,
  patchwork::wrap_plots(s5f, s5g, ncol = 2, widths = c(1.2, 1)),
  patchwork::wrap_plots(s5h, s5i, ncol = 2),
  ncol = 1,
  heights = c(2.0, 0.85, 0.78)
) + patchwork::plot_annotation(
  title = "Supplementary Figure 5 | Initial ploidy, dose, and sample composition"
)
panel_rows[[length(panel_rows) + 1L]] <- save_composite(
  s5,
  figure_dir,
  "SuppFig5",
  "panel_SuppFig5_composite",
  18,
  23
)

message("Generating Supplementary Figure 6 composite.")
endpoint_limits <- range(tumor$endpoint_ploidy)
tumor_2n <- tumor[tumor$initial_ploidy == "2N", , drop = FALSE]
tumor_4n <- tumor[tumor$initial_ploidy == "4N", , drop = FALSE]
s6a <- make_umap_continuous(
  tumor,
  "endpoint_ploidy",
  "Endpoint ploidy in all tumors",
  "Endpoint ploidy",
  "A",
  point_size,
  limits = endpoint_limits,
  subtitle = sprintf(
    "Range %.3f-%.3f",
    endpoint_limits[[1L]],
    endpoint_limits[[2L]]
  )
)
s6b <- make_umap_continuous(
  tumor_2n,
  "endpoint_ploidy",
  "Endpoint ploidy in initial 2N tumors",
  "Endpoint ploidy",
  "B",
  point_size,
  limits = endpoint_limits
)
s6c <- make_umap_continuous(
  tumor_4n,
  "endpoint_ploidy",
  "Endpoint ploidy in initial 4N tumors",
  "Endpoint ploidy",
  "C",
  point_size,
  limits = endpoint_limits
)
s6_mouse <- tumor
s6_mouse$mouse_panel <- factor(
  unname(panel_lookup[as.character(s6_mouse$mouse)]),
  levels = panel_metadata$mouse_panel
)
s6d <- ggplot2::ggplot(
  shuffle_cells(s6_mouse, plot_seed + 2000L),
  ggplot2::aes(UMAP_1, UMAP_2, color = endpoint_ploidy)
) +
  ggplot2::geom_point(
    size = max(point_size * 1.35, 0.14),
    alpha = 0.88,
    stroke = 0
  ) +
  ggplot2::facet_wrap(~mouse_panel, ncol = 4, nrow = 4, drop = FALSE) +
  ggplot2::scale_color_gradient(
    low = "#2C7BB6",
    high = "#D7191C",
    limits = endpoint_limits,
    name = "Endpoint ploidy"
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Endpoint ploidy by mouse of origin",
    subtitle = "All panels use the same endpoint-ploidy scale",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  umap_theme(8.5) +
  ggplot2::theme(
    legend.position = "top",
    legend.justification = "right",
    legend.box.just = "right",
    strip.text = ggplot2::element_text(size = 7.2),
    panel.spacing = grid::unit(0.10, "lines")
  )
s6d <- add_tag(s6d, "D")
s6 <- patchwork::wrap_plots(
  patchwork::wrap_plots(s6a, s6b, s6c, ncol = 1),
  s6d,
  ncol = 2,
  widths = c(1, 3)
) + patchwork::plot_annotation(
  title = "Supplementary Figure 6 | Endpoint tumor ploidy"
)
panel_rows[[length(panel_rows) + 1L]] <- save_composite(
  s6,
  figure_dir,
  "SuppFig6",
  "panel_SuppFig6_composite",
  20,
  16.5
)

message("Generating Supplementary Figure 7 composite.")
read_matrix <- function(filename, label) {
  data <- read_tsv(file.path(table_dir, filename), label)
  if (!identical(names(data), c("pathway", cluster_levels))) {
    stop(label, " does not match the canonical cluster order", call. = FALSE)
  }
  matrix_data <- as.matrix(data[, cluster_levels, drop = FALSE])
  storage.mode(matrix_data) <- "double"
  rownames(matrix_data) <- data$pathway
  if (!identical(dim(matrix_data), c(20L, length(cluster_levels))) ||
      any(!is.finite(matrix_data))) {
    stop(label, " must be a finite 20 x 9 matrix", call. = FALSE)
  }
  matrix_data
}
ora_matrix <- read_matrix(
  "si_figure7_cluster_Hallmark_ORA_annotation_score_heatmap_top20_matrix.tsv",
  "SI Figure 7 ORA matrix"
)
gsea_matrix <- read_matrix(
  "si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv",
  "SI Figure 7 GSEA matrix"
)
s7a <- heatmap_plot(
  ora_matrix,
  "Cluster Hallmark ORA annotation score",
  FALSE,
  "A"
)
s7b <- heatmap_plot(
  gsea_matrix,
  "Cluster Hallmark GSEA NES",
  TRUE,
  "B"
)
s7 <- patchwork::wrap_plots(s7a, s7b, ncol = 2) +
  patchwork::plot_annotation(
    title = "Supplementary Figure 7 | Cluster Hallmark pathway analysis"
  )
panel_rows[[length(panel_rows) + 1L]] <- save_composite(
  s7,
  figure_dir,
  "SuppFig7",
  "panel_SuppFig7_composite",
  20,
  8.8
)

panel_contract <- do.call(rbind, panel_rows)
rownames(panel_contract) <- NULL
expected_figures <- sort(panel_contract$filename)
observed_figures <- sort(list.files(figure_dir, pattern = "[.](pdf|png)$"))
if (!identical(expected_figures, observed_figures)) {
  stop(
    "Supplementary composite inventory mismatch; missing=",
    paste(setdiff(expected_figures, observed_figures), collapse = ","),
    "; unexpected=",
    paste(setdiff(observed_figures, expected_figures), collapse = ","),
    call. = FALSE
  )
}
if (
  any(
    !file.exists(file.path(figure_dir, expected_figures)) |
      file.info(file.path(figure_dir, expected_figures))$size <= 0
  )
) {
  stop("One or more Supplementary composite files are empty", call. = FALSE)
}
write_tsv(panel_contract, file.path(metadata_dir, "panel_contract.tsv"))

input_paths <- c(
  config_path,
  file.path(cache_dir, "manifest.tsv"),
  source_paths
)
input_manifest <- data.frame(
  role = c(
    "figure7_config",
    if (allow_generated_human_only_si7) {
      "si_figures_generated_cache_manifest"
    } else {
      "si_figures_cache_manifest"
    },
    rep(
      if (allow_generated_human_only_si7) {
        "si_figures_generated_table"
      } else {
        "si_figures_frozen_table"
      },
      length(source_paths)
    )
  ),
  repo_relative_path = vapply(
    input_paths,
    repo_relative,
    character(1),
    root = root
  ),
  sha256 = vapply(input_paths, file_sha256, character(1)),
  bytes = as.numeric(file.info(input_paths)$size),
  stringsAsFactors = FALSE
)
if (!is.na(upstream_manifest_path)) {
  upstream <- read_tsv(
    upstream_manifest_path,
    "raw SI build input manifest"
  )
  required_upstream_columns <- c("role", "locator", "sha256", "bytes")
  if (!identical(names(upstream), required_upstream_columns) ||
      anyNA(upstream$role) ||
      any(!nzchar(upstream$role)) ||
      anyDuplicated(upstream$role) ||
      any(!grepl("^[0-9a-f]{64}$", upstream$sha256)) ||
      any(!is.finite(as.numeric(upstream$bytes))) ||
      any(grepl("^/", upstream$locator))) {
    stop("Raw SI build input manifest is malformed or nonportable", call. = FALSE)
  }
  upstream_rows <- data.frame(
    role = paste0("raw_build_", upstream$role),
    repo_relative_path = as.character(upstream$locator),
    sha256 = as.character(upstream$sha256),
    bytes = as.numeric(upstream$bytes),
    stringsAsFactors = FALSE
  )
  upstream_manifest_row <- data.frame(
    role = "raw_build_input_manifest",
    repo_relative_path = repo_relative(upstream_manifest_path, root),
    sha256 = file_sha256(upstream_manifest_path),
    bytes = as.numeric(file.info(upstream_manifest_path)$size),
    stringsAsFactors = FALSE
  )
  input_manifest <- rbind(
    input_manifest,
    upstream_manifest_row,
    upstream_rows
  )
}
write_tsv(input_manifest, file.path(metadata_dir, "input_manifest.tsv"))

source_revision <- paste(unique(cache_manifest$source_revision), collapse = ";")
run_config <- data.frame(
  key = c(
    "schema_version",
    "module",
    "figures",
    "figure_file_count",
    "table_mode",
    "cache_file_count",
    "plot_shuffle_seed",
    "cluster_order",
    "si7_feature_species_policy",
    "si7_gene_set_database",
    "si7_canonical_publication_allowed",
    if (allow_generated_human_only_si7) {
      "generated_table_source_revision"
    } else {
      "frozen_table_source_revision"
    }
  ),
  value = c(
    "1",
    "si_figures",
    "4,5,6,7",
    as.character(nrow(panel_contract)),
    if (allow_generated_human_only_si7) {
      "run_scoped_generated_human_only_plot_tables"
    } else {
      "frozen_plot_tables_only"
    },
    as.character(length(cache_files)),
    as.character(plot_seed),
    paste(cluster_levels, collapse = ","),
    as.character(
      if (allow_generated_human_only_si7) {
        si_config$raw_rebuild_species_policy
      } else {
        si_config$si7_feature_species_policy
      }
    ),
    as.character(
      if (allow_generated_human_only_si7) {
        si_config$raw_rebuild_gene_set_database
      } else {
        si_config$si7_gene_set_database
      }
    ),
    tolower(as.character(!allow_generated_human_only_si7)),
    source_revision
  ),
  stringsAsFactors = FALSE
)
write_tsv(run_config, file.path(metadata_dir, "run_config.tsv"))

input_qc <- data.frame(
  key = c(
    "all_cells",
    "included_tumor_cells",
    "cellline_cells",
    "tumor_mice",
    "clusters",
    "endpoint_audit_rows",
    "tumor_endpoint_ploidy_min",
    "tumor_endpoint_ploidy_max",
    "composite_files"
  ),
  value = c(
    nrow(seurat),
    n_tumor,
    sum(seurat$context == "CellLine"),
    n_mice,
    length(cluster_levels),
    nrow(endpoint_audit),
    format(endpoint_limits[[1L]], digits = 16),
    format(endpoint_limits[[2L]], digits = 16),
    nrow(panel_contract)
  ),
  stringsAsFactors = FALSE
)
write_tsv(input_qc, file.path(metadata_dir, "input_qc.tsv"))

git_revision <- tryCatch(
  system2(
    "git",
    c("-C", shQuote(root), "rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
  )[[1L]],
  error = function(error) "unavailable"
)
provenance <- data.frame(
  key = c(
    "schema_version",
    "artifact",
    "entrypoint",
    "entrypoint_sha256",
    "source_code_revision",
    "table_cache",
    "table_cache_manifest_sha256",
    "analysis_mode",
    "si7_feature_species_policy",
    "si7_gene_set_database",
    "si7_canonical_publication_allowed",
    if (allow_generated_human_only_si7) {
      "si7_generated_matrix_note"
    } else {
      "si7_frozen_matrix_note"
    }
  ),
  value = c(
    "1",
    "supplementary_figures_4_7",
    repo_relative(script, root),
    file_sha256(script),
    git_revision,
    repo_relative(cache_dir, root),
    file_sha256(file.path(cache_dir, "manifest.tsv")),
    if (allow_generated_human_only_si7) {
      paste(
        "plot-only rendering of a run-scoped raw rebuild;",
        "canonical publication is prohibited"
      )
    } else {
      "plot-only; raw Seurat and enrichment dependencies are intentionally absent"
    },
    as.character(
      if (allow_generated_human_only_si7) {
        si_config$raw_rebuild_species_policy
      } else {
        si_config$si7_feature_species_policy
      }
    ),
    as.character(
      if (allow_generated_human_only_si7) {
        si_config$raw_rebuild_gene_set_database
      } else {
        si_config$si7_gene_set_database
      }
    ),
    tolower(as.character(!allow_generated_human_only_si7)),
    if (allow_generated_human_only_si7) {
      paste(
        "Generated SI7 statistics use a newly normalized object containing",
        "only exact GRCh38-prefixed RNA counts; Hallmark matching preserves",
        "symbol case. These rebuilt tables are not canonical until reviewed."
      )
    } else {
      paste(
        "Frozen matrices were recalculated from Tao's cluster DEG cache after",
        "explicitly retaining GRCh features and excluding GRCm39 features."
      )
    }
  ),
  stringsAsFactors = FALSE
)
write_tsv(
  provenance,
  file.path(metadata_dir, "si_figures_provenance.tsv")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(metadata_dir, "sessionInfo.txt")
)
writeLines(
  c(
    "Supplementary Figures 4-7 generation completed.",
    if (allow_generated_human_only_si7) {
      "Mode: generated human-only plot-facing tables (noncanonical)"
    } else {
      "Mode: frozen plot-facing tables only"
    },
    paste0("Composite files: ", nrow(panel_contract)),
    paste0("All cells: ", nrow(seurat)),
    paste0("Included tumor cells: ", n_tumor),
    paste0("Tumor mice: ", n_mice),
    paste0(
      "SI7 feature policy: ",
      as.character(
        if (allow_generated_human_only_si7) {
          si_config$raw_rebuild_species_policy
        } else {
          si_config$si7_feature_species_policy
        }
      )
    )
  ),
  file.path(metadata_dir, "run_summary.txt")
)

message("Supplementary Figures completed: ", output_dir)
