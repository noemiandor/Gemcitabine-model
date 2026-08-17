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
      "    [--all-ploidy Data/in-vivo/all_ploidy.tsv] \\",
      "    [--cbs-dir Data/in-vivo/scRNAseq_Numbat] \\",
      paste(
        "    [--cellcycle-pseudotime",
        "Data/in-vivo/figure7/processed/CellCycleCells_",
        "pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv] \\",
        sep = ""
      ),
      paste(
        "    [--injected-reference-dir",
        "Data/in-vivo/scRNAseq_Numbat/injected_reference] \\",
      ),
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

percent_labels <- function(values) paste0(round(100 * values), "%")

make_composition_plot <- function(
  data,
  fill_colors,
  title,
  x_title,
  y_title,
  tag,
  use_proportion,
  subtitle = NULL
) {
  y_field <- if (use_proportion) "proportion" else "n_cells"
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(group_value, .data[[y_field]], fill = fill_value)
  ) +
    ggplot2::geom_col(width = 0.82, color = "white", linewidth = 0.12) +
    ggplot2::scale_fill_manual(values = fill_colors, drop = FALSE) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_title,
      y = y_title,
      fill = NULL
    ) +
    shared_context_figure_theme(9.5) +
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
  shared_context_add_tag(plot, tag)
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
composition_helper_path <- file.path(
  root,
  "Code",
  "in-vivo",
  "SI_figures",
  "normalized_composition.R"
)
if (!file.exists(composition_helper_path)) {
  stop("Missing normalized-composition helper: ", composition_helper_path, call. = FALSE)
}
sys.source(composition_helper_path, envir = environment())
shared_context_helper_path <- file.path(
  root,
  "Code",
  "in-vivo",
  "SI_figures",
  "shared_context_panels.R"
)
if (!file.exists(shared_context_helper_path)) {
  stop("Missing shared context-panel helper: ", shared_context_helper_path, call. = FALSE)
}
sys.source(shared_context_helper_path, envir = environment())
copy_number_helper_path <- file.path(
  root,
  "Code",
  "in-vivo",
  "SI_figures",
  "copy_number_heatmap.R"
)
if (!file.exists(copy_number_helper_path)) {
  stop("Missing copy-number heatmap helper: ", copy_number_helper_path,
       call. = FALSE)
}
sys.source(copy_number_helper_path, envir = environment())
figure7_shared_paths <- file.path(
  root,
  "Code",
  "in-vivo",
  "figure7",
  "src",
  c("common_io.R", "tgi_statistics.R", "tgi_panels.R")
)
if (any(!file.exists(figure7_shared_paths))) {
  stop(
    "Missing Figure 7 density-localization helper(s): ",
    paste(figure7_shared_paths[!file.exists(figure7_shared_paths)], collapse = ", "),
    call. = FALSE
  )
}
for (path in figure7_shared_paths) {
  sys.source(path, envir = environment())
}
weighted_ploidy_path <- file.path(
  root, "Data", "in-vivo", "weighted_ploidy.py"
)
if (!file.exists(weighted_ploidy_path)) {
  stop("Missing endpoint-ploidy derivation helper: ", weighted_ploidy_path,
       call. = FALSE)
}
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
cellcycle_pseudotime_path <- resolve_path(
  arg_value(
    args,
    "cellcycle-pseudotime",
    paste0(
      "Data/in-vivo/figure7/processed/",
      "CellCycleCells_pseudotime_distribution_per_sample_",
      "cell_level_with_ploidy_dose_tgi.csv"
    )
  ),
  root,
  must_work = TRUE
)
output_dir <- resolve_path(
  arg_value(args, "output-dir", "Results/in-vivo/SI_figures"),
  root,
  must_work = FALSE
)
all_ploidy_path <- resolve_path(
  arg_value(args, "all-ploidy", "Data/in-vivo/all_ploidy.tsv"),
  root,
  must_work = TRUE
)
cbs_dir <- resolve_path(
  arg_value(args, "cbs-dir", "Data/in-vivo/scRNAseq_Numbat"),
  root,
  must_work = TRUE
)
injected_reference_dir <- resolve_path(
  arg_value(
    args,
    "injected-reference-dir",
    file.path(cbs_dir, "injected_reference")
  ),
  root,
  must_work = TRUE
)
cbs_manifest_path <- resolve_path(
  file.path(cbs_dir, "cbs_manifest.tsv"),
  root,
  must_work = TRUE
)
injected_reference_manifest_path <- resolve_path(
  file.path(injected_reference_dir, "reference_manifest.tsv"),
  root,
  must_work = TRUE
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
figure7_config <- figure7_read_config(config_path)
figure7_config <- figure7_attach_density_localization_config(
  figure7_config,
  file.path(root, "Code", "in-vivo", "figure7", "density_localization_config.yaml")
)
density_localization_config_path <- attr(
  figure7_config, "density_localization_config_path", exact = TRUE
)
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
joint_umap_point_size <- shared_context_publication_umap_point_size(nrow(seurat))
tumor_umap_point_size <- shared_context_publication_umap_point_size(n_tumor)
faceted_umap_point_size <- shared_context_publication_umap_point_size(
  max(table(tumor$mouse)),
  faceted = TRUE
)

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

message("Generating Supplementary Figure 4 composite.")
shared_si4 <- shared_context_build_si4_panels(
  data = seurat,
  cluster_levels = cluster_levels,
  cluster_colors = cluster_colors,
  ploidy_colors = ploidy_colors,
  context_colors = context_colors,
  point_size = joint_umap_point_size,
  plot_seed = plot_seed
)
s4a <- shared_si4$plots$cluster
s4b <- shared_si4$plots$initial_ploidy
s4c <- shared_si4$plots$context
s4e_result <- shared_si4$composition_result
s4e <- shared_si4$plots$composition
s4d <- shared_context_make_umap_continuous(
  seurat,
  "s_phase_score",
  "UMAP by S phase score",
  "S phase score",
  "D",
  max(joint_umap_point_size, 0.70),
  plot_seed,
  limits = range(seurat$s_phase_score),
  diverging = TRUE,
  mid_color = "#A6A6A6"
)
selected_cellcycle_clusters <- c("4c", "6", "10")
selected_tumor_cells <- sum(
  s4_context$n_cells[
    as.character(s4_context$group_value) %in% selected_cellcycle_clusters &
      as.character(s4_context$fill_value) == "Tumor"
  ]
)
metadata_selected_tumor_cells <- sum(
  as.character(tumor$cluster) %in% selected_cellcycle_clusters
)
if (!identical(as.numeric(selected_tumor_cells), as.numeric(metadata_selected_tumor_cells))) {
  stop("SI4F selected tumor-cell count does not reconcile to cell metadata", call. = FALSE)
}
s4f <- make_composition_plot(
  s4_context,
  context_colors,
  "Tumor and CellLine cell counts by cluster",
  "Cluster",
  "Number of cells",
  "F",
  FALSE,
  subtitle = sprintf(
    "Selected CellCycle clusters 4c + 6 + 10 contain n = %s tumor cells",
    format(selected_tumor_cells, big.mark = ",", scientific = FALSE)
  )
)
s4g_result <- make_normalized_composition_plot(
  data = seurat,
  unit_col = "mouse",
  cluster_col = "cluster",
  group_col = "initial_ploidy",
  cluster_levels = cluster_levels,
  group_levels = ploidy_levels,
  fill_colors = ploidy_colors,
  bar_axis = "cluster",
  strata_cols = "context",
  title = "Equal-sample cluster composition by initial ploidy",
  subtitle = paste(
    "Each sample contributes equally within initial ploidy;",
    "stars mark ploidy enrichment within context at BH FDR <= 0.05"
  ),
  x_title = "Cluster",
  y_title = "Mean within-sample cluster proportion",
  legend_title = "Initial ploidy",
  tag = "G",
  theme_function = shared_context_figure_theme,
  x_text_angle = 0
)
s4g <- s4g_result$plot
s4h <- make_composition_plot(
  s4_ploidy,
  ploidy_colors,
  "Initial ploidy cell counts by cluster",
  "Cluster",
  "Number of cells",
  "H",
  FALSE
)

# Supplementary Figure 4 already defines the S-phase landscape and the exact
# 2,881-cell clusters 4c/6/10 subset. The mouse-balanced localization formerly
# shown beneath main Figure 7H therefore belongs here as a full-width panel.
if (!identical(
      file_sha256(cellcycle_pseudotime_path),
      as.character(figure7_config$inputs$cellcycle_sha256)
    )) {
  stop(
    "SI4I CellCycle pseudotime input differs from the reviewed Figure 7 input",
    call. = FALSE
  )
}
cellcycle_pseudotime <- utils::read.csv(
  cellcycle_pseudotime_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_localization_columns <- c(
  "cell_id", "sample_id", "initial_ploidy",
  "gemcitabine_dose_mg_per_kg", "pseudotime"
)
if (!all(required_localization_columns %in% names(cellcycle_pseudotime)) ||
    nrow(cellcycle_pseudotime) != 2881L ||
    anyDuplicated(cellcycle_pseudotime$cell_id) ||
    any(!is.finite(cellcycle_pseudotime$pseudotime))) {
  stop(
    "SI4I requires the exact reviewed 2,881-cell pseudotime table",
    call. = FALSE
  )
}
localization_samples <- unique(data.frame(
  sample_id = as.character(cellcycle_pseudotime$sample_id),
  initial_ploidy = as.character(cellcycle_pseudotime$initial_ploidy),
  dose_mg = suppressWarnings(as.numeric(
    cellcycle_pseudotime$gemcitabine_dose_mg_per_kg
  )),
  stringsAsFactors = FALSE
))
localization_samples <- localization_samples[order(
  localization_samples$sample_id
), , drop = FALSE]
if (nrow(localization_samples) != 16L ||
    anyDuplicated(localization_samples$sample_id) ||
    any(!is.finite(localization_samples$dose_mg))) {
  stop("SI4I requires one reviewed treatment row for each of 16 mice",
       call. = FALSE)
}
localization <- figure7_density_localization(
  cellcycle_pseudotime,
  localization_samples,
  figure7_config
)
write_tsv(
  localization$grid,
  file.path(metadata_dir, "si_figure4I_density_localization_grid.tsv")
)
write_tsv(
  localization$intervals,
  file.path(metadata_dir, "si_figure4I_density_localization_intervals.tsv")
)
write_tsv(
  localization$test,
  file.path(metadata_dir, "si_figure4I_density_localization_test.tsv")
)
s4i <- figure7_panel_b_localization_plot(
  localization$grid,
  localization$intervals,
  localization$test
) +
  ggplot2::labs(
    title = "Mouse-balanced localization of gemcitabine-associated cell excess",
    x = "Cell-cycle pseudotime",
    y = "Gemcitabine - vehicle density difference"
  ) +
  shared_context_figure_theme(9)
s4i <- shared_context_add_tag(s4i, "I")
s4 <- patchwork::wrap_plots(
  patchwork::wrap_plots(s4d, s4f, ncol = 2, widths = c(1, 1.55)),
  s4i,
  ncol = 1,
  heights = c(1.12, 0.88)
) + patchwork::plot_annotation(
  title = paste(
    "Supplementary Figure 4 | Cell-cycle state definition and",
    "gemcitabine-associated pseudotime localization"
  )
)
panel_rows <- list(
  save_composite(
    s4,
    figure_dir,
    "SuppFig4",
    "panel_SuppFig4_composite",
    14,
    10.5
  )
)

message("Generating Supplementary Figure 5 composite.")
s5a <- shared_context_make_umap_discrete(
  tumor,
  "cluster",
  cluster_colors,
  "UMAP by cluster",
  "Cluster",
  "A",
  tumor_umap_point_size,
  plot_seed,
  sprintf("Tumor cells; n = %s", format(n_tumor, big.mark = ",")),
  labels = TRUE
)
s5b <- shared_context_make_umap_discrete(
  tumor,
  "initial_ploidy",
  ploidy_colors,
  "UMAP by initial tumor ploidy",
  "Initial ploidy",
  "B",
  tumor_umap_point_size,
  plot_seed
)
s5c <- shared_context_make_umap_discrete(
  tumor,
  "dose",
  dose_colors,
  "UMAP by Gemcitabine dose",
  "Gemcitabine dose",
  "C",
  tumor_umap_point_size,
  plot_seed
)
s5d <- shared_context_make_umap_continuous(
  tumor,
  "s_phase_score",
  "UMAP by S phase score",
  "S phase score",
  "D",
  tumor_umap_point_size,
  plot_seed,
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
    size = faceted_umap_point_size,
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
  shared_context_umap_theme(8.5) +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_text(size = 7.2),
    panel.spacing = grid::unit(0.10, "lines"),
    plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
  )
s5e <- shared_context_add_tag(s5e, "E")

tumor$ploidy_dose_group <- paste(
  as.character(tumor$initial_ploidy),
  as.character(tumor$dose),
  sep = " | "
)
ploidy_dose_levels <- unlist(lapply(
  ploidy_levels,
  function(ploidy) paste(ploidy, dose_levels, sep = " | ")
))

s5f_result <- make_normalized_composition_plot(
  data = tumor,
  unit_col = "mouse",
  cluster_col = "cluster",
  group_col = "mouse",
  cluster_levels = cluster_levels,
  group_levels = mouse_levels,
  fill_colors = cluster_colors,
  bar_axis = "group",
  x_col = "mouse",
  x_levels = mouse_levels,
  facet_cols = c("initial_ploidy", "dose"),
  facet_formula = stats::as.formula(". ~ initial_ploidy + dose"),
  facet_type = "grid",
  facet_scales = "free_x",
  facet_space = "free_x",
  title = "Within-mouse cluster composition",
  subtitle = paste(
    "Every mouse sums to 100%; one biological sample per bar,",
    "so this panel is descriptive and has no inferential stars"
  ),
  x_title = "Mouse",
  y_title = "Within-mouse cluster proportion",
  legend_title = "Cluster",
  tag = "F",
  theme_function = shared_context_figure_theme,
  x_text_angle = 50,
  test_mode = "descriptive",
  descriptive_reason = "one biological sample per displayed mouse"
)
s5f <- s5f_result$plot + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 50, hjust = 1, size = 7.5)
)

s5g_result <- make_normalized_composition_plot(
  data = tumor,
  unit_col = "mouse",
  cluster_col = "cluster",
  group_col = "ploidy_dose_group",
  cluster_levels = cluster_levels,
  group_levels = ploidy_dose_levels,
  fill_colors = cluster_colors,
  bar_axis = "group",
  x_col = "dose",
  x_levels = dose_levels,
  facet_cols = "initial_ploidy",
  strata_cols = "initial_ploidy",
  facet_formula = stats::as.formula("~ initial_ploidy"),
  facet_type = "wrap",
  facet_nrow = 1,
  title = "Equal-mouse composition by initial ploidy and dose",
  subtitle = paste(
    "Within-mouse proportions are averaged with equal mouse weights;",
    "stars mark dose enrichment within ploidy at BH FDR <= 0.05"
  ),
  x_title = "Gemcitabine dose",
  y_title = "Mean within-mouse cluster proportion",
  legend_title = "Cluster",
  tag = "G",
  theme_function = shared_context_figure_theme,
  x_text_angle = 25,
  show_n = TRUE
)
s5g <- s5g_result$plot

s5h_result <- make_normalized_composition_plot(
  data = tumor,
  unit_col = "mouse",
  cluster_col = "cluster",
  group_col = "dose",
  cluster_levels = cluster_levels,
  group_levels = dose_levels,
  fill_colors = dose_colors,
  bar_axis = "cluster",
  strata_cols = "initial_ploidy",
  title = "Equal-mouse cluster composition by gemcitabine dose",
  subtitle = paste(
    "Each mouse contributes equally within dose;",
    "stars mark dose enrichment within ploidy at BH FDR <= 0.05"
  ),
  x_title = "Cluster",
  y_title = "Mean within-mouse cluster proportion",
  legend_title = "Gemcitabine dose",
  tag = "H",
  theme_function = shared_context_figure_theme,
  x_text_angle = 0
)
s5h <- s5h_result$plot

s5i_result <- make_normalized_composition_plot(
  data = tumor,
  unit_col = "mouse",
  cluster_col = "cluster",
  group_col = "initial_ploidy",
  cluster_levels = cluster_levels,
  group_levels = ploidy_levels,
  fill_colors = ploidy_colors,
  bar_axis = "cluster",
  strata_cols = "dose",
  title = "Equal-mouse cluster composition by initial ploidy",
  subtitle = paste(
    "Each mouse contributes equally within initial ploidy;",
    "stars mark ploidy enrichment within dose at BH FDR <= 0.05"
  ),
  x_title = "Cluster",
  y_title = "Mean within-mouse cluster proportion",
  legend_title = "Initial ploidy",
  tag = "I",
  theme_function = shared_context_figure_theme,
  x_text_angle = 0
)
s5i <- s5i_result$plot

s5 <- patchwork::wrap_plots(
  s5e,
  s5g,
  ncol = 1,
  heights = c(1.65, 0.85)
) + patchwork::plot_annotation(
  title = paste(
    "Supplementary Figure 5 | Mouse-level tumor landscape and",
    "composition by initial ploidy and dose"
  )
)
panel_rows[[length(panel_rows) + 1L]] <- save_composite(
  s5,
  figure_dir,
  "SuppFig5",
  "panel_SuppFig5_composite",
  16,
  11.5
)

message("Generating Supplementary Figure 6 composite.")
endpoint_limits <- range(tumor$endpoint_ploidy)
tumor_2n <- tumor[tumor$initial_ploidy == "2N", , drop = FALSE]
tumor_4n <- tumor[tumor$initial_ploidy == "4N", , drop = FALSE]
s6a <- shared_context_make_umap_continuous(
  tumor,
  "endpoint_ploidy",
  "NUMBAT-derived ploidy in all tumors",
  "NUMBAT-derived ploidy",
  "A",
  tumor_umap_point_size,
  plot_seed,
  limits = endpoint_limits,
  subtitle = sprintf(
    "Range %.3f-%.3f",
    endpoint_limits[[1L]],
    endpoint_limits[[2L]]
  )
)
s6b <- shared_context_make_umap_continuous(
  tumor_2n,
  "endpoint_ploidy",
  "NUMBAT-derived ploidy in 2N-origin tumors",
  "NUMBAT-derived ploidy",
  "B",
  tumor_umap_point_size,
  plot_seed,
  limits = endpoint_limits
)
s6c <- shared_context_make_umap_continuous(
  tumor_4n,
  "endpoint_ploidy",
  "NUMBAT-derived ploidy in 4N-origin tumors",
  "NUMBAT-derived ploidy",
  "C",
  tumor_umap_point_size,
  plot_seed,
  limits = endpoint_limits
)
s6_mouse <- tumor
s6_mouse$mouse_panel <- factor(
  unname(panel_lookup[as.character(s6_mouse$mouse)]),
  levels = panel_metadata$mouse_panel
)
s6d <- ggplot2::ggplot(
  shared_context_shuffle_cells(s6_mouse, plot_seed + 2000L),
  ggplot2::aes(UMAP_1, UMAP_2, color = endpoint_ploidy)
) +
  ggplot2::geom_point(
    size = faceted_umap_point_size,
    alpha = 0.88,
    stroke = 0
  ) +
  ggplot2::facet_wrap(~mouse_panel, ncol = 4, nrow = 4, drop = FALSE) +
  ggplot2::scale_color_gradient(
    low = "#2C7BB6",
    high = "#D7191C",
    limits = endpoint_limits,
    name = "NUMBAT-derived ploidy"
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "NUMBAT-derived ploidy by mouse of origin",
    subtitle = "All panels use the same chromosome-length-weighted scale",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  shared_context_umap_theme(8.5) +
  ggplot2::theme(
    legend.position = "top",
    legend.justification = "right",
    legend.box.just = "right",
    strip.text = ggplot2::element_text(size = 7.2),
    panel.spacing = grid::unit(0.10, "lines")
  )
s6d <- shared_context_add_tag(s6d, "D")
message("Validating the Figure 7J downstream CBS copy-number heatmap inputs.")
cbs_collection <- si_copy_number_read_collection(
  cbs_dir = cbs_dir,
  all_ploidy_path = all_ploidy_path,
  endpoint_audit = endpoint_audit,
  sample_metadata = sample_metadata
)
cbs_harmonized <- si_copy_number_harmonize(cbs_collection)
cbs_sample_summary <- si_copy_number_sample_summary(
  cbs_harmonized$cell_annotations
)
cbs_injected_references <- si_copy_number_read_injected_references(
  injected_reference_dir
)
cbs_reduction_summary <- si_copy_number_reduction_summary(
  cbs_injected_references$cells,
  cbs_sample_summary,
  cbs_harmonized$cell_annotations
)
cbs_separation_summary <- si_copy_number_separation_summary(
  cbs_reduction_summary
)
cbs_4n_reduction <- cbs_reduction_summary[
  cbs_reduction_summary$injected_origin == "4N", , drop = FALSE
]
cbs_2n_reduction <- cbs_reduction_summary[
  cbs_reduction_summary$injected_origin == "2N", , drop = FALSE
]
if (nrow(cbs_2n_reduction) != 1L || nrow(cbs_4n_reduction) != 1L) {
  stop("Expected one descriptive ploidy-change row per injected origin",
       call. = FALSE)
}
cbs_endpoint_summaries <- si_copy_number_endpoint_summaries(
  cbs_sample_summary
)
write_tsv(
  cbs_harmonized$cell_annotations,
  file.path(metadata_dir, "figure7J_copy_number_cell_annotations.tsv")
)
write_tsv(
  cbs_harmonized$chromosomes,
  file.path(metadata_dir, "figure7J_copy_number_chromosomes.tsv")
)
write_tsv(
  cbs_harmonized$chromosome_schema_audit,
  file.path(metadata_dir, "figure7J_cbs_schema_chromosome_coverage.tsv")
)
write_tsv(
  cbs_harmonized$chromosome_file_audit,
  file.path(metadata_dir, "figure7J_cbs_file_chromosome_coverage.tsv")
)
write_tsv(
  cbs_endpoint_summaries$samples,
  file.path(metadata_dir, "si_figure6_postprocessed_copy_number_score_by_mouse.tsv")
)
write_tsv(
  cbs_endpoint_summaries$origin_summary,
  file.path(metadata_dir, "si_figure6_postprocessed_copy_number_score_by_injected_origin.tsv")
)
write_tsv(
  cbs_endpoint_summaries$origin_dose_summary,
  file.path(metadata_dir, "si_figure6_postprocessed_copy_number_score_by_origin_and_dose.tsv")
)
write_tsv(
  cbs_injected_references$cells,
  file.path(metadata_dir, "si_figure6_injected_cell_reference_ploidy.tsv")
)
write_tsv(
  cbs_reduction_summary,
  file.path(metadata_dir, "si_figure6_injected_to_endpoint_ploidy_change.tsv")
)
write_tsv(
  cbs_separation_summary,
  file.path(
    metadata_dir,
    "si_figure6_injected_to_endpoint_ploidy_separation.tsv"
  )
)
saveRDS(
  cbs_harmonized$matrix,
  file.path(metadata_dir, "figure7J_copy_number_heatmap_matrix.rds"),
  compress = "gzip"
)
cbs_reference_plot_data <- cbs_injected_references$cells
cbs_reference_plot_data$initial_ploidy <- cbs_reference_plot_data$injected_origin
cbs_stage_levels <- c("Injected-cell\nreference", "Endpoint\ntumor")
cbs_reference_plot_data$stage <- factor(
  "Injected-cell\nreference", levels = cbs_stage_levels
)
cbs_reference_plot_data$x_position <- 1
cbs_reference_plot_data$value <- cbs_reference_plot_data$ploidy
cbs_mouse_plot_data <- cbs_endpoint_summaries$samples
cbs_mouse_plot_data$stage <- factor(
  "Endpoint\ntumor", levels = cbs_stage_levels
)
cbs_mouse_plot_data$x_position <- 2
cbs_mouse_plot_data$value <-
  cbs_mouse_plot_data$mean_postprocessed_copy_number_score
cbs_mouse_dose_colors <- c(
  "0" = "#666666", "30" = "#D95F02", "120" = "#1B9E77"
)
cbs_reference_means <- stats::aggregate(
  value ~ initial_ploidy + stage,
  data = cbs_reference_plot_data,
  FUN = mean
)
cbs_endpoint_means <- stats::aggregate(
  value ~ initial_ploidy + stage,
  data = cbs_mouse_plot_data,
  FUN = mean
)
cbs_stage_means <- rbind(cbs_reference_means, cbs_endpoint_means)
cbs_stage_means$x_position <- ifelse(
  cbs_stage_means$stage == "Injected-cell\nreference", 1, 2
)
cbs_reduction_labels <- cbs_reduction_summary
cbs_reduction_labels$initial_ploidy <- cbs_reduction_labels$injected_origin
cbs_reduction_labels$label <- sprintf(
  "Change: %+.2f (%+.1f%%)",
  cbs_reduction_labels$absolute_change,
  cbs_reduction_labels$relative_change_percent
)
cbs_reduction_labels$x <- 1.5
cbs_reduction_labels$y <- vapply(
  cbs_reduction_labels$injected_origin,
  function(origin) {
    max(c(
      cbs_reference_plot_data$value[
        cbs_reference_plot_data$initial_ploidy == origin
      ],
      cbs_mouse_plot_data$value[cbs_mouse_plot_data$initial_ploidy == origin]
    )) + 0.13
  },
  numeric(1L)
)
cbs_reduction_arrows <- data.frame(
  initial_ploidy = cbs_reduction_summary$injected_origin,
  x = 1,
  xend = 2,
  y = cbs_reduction_summary$reference_mean_ploidy,
  yend = cbs_reduction_summary$endpoint_mouse_balanced_mean_ploidy,
  stringsAsFactors = FALSE
)
s6e_mouse <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = cbs_reduction_arrows,
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.55,
    color = "#555555",
    arrow = grid::arrow(length = grid::unit(0.07, "inches"), type = "closed"),
    inherit.aes = FALSE
  ) +
  ggplot2::geom_point(
    data = cbs_reference_plot_data,
    ggplot2::aes(x_position, value),
    shape = 4,
    size = 1.6,
    stroke = 0.55,
    color = "#333333",
    position = ggplot2::position_jitter(
      width = 0.10, height = 0, seed = plot_seed + 61L
    )
  ) +
  ggplot2::geom_point(
    data = cbs_mouse_plot_data,
    ggplot2::aes(x_position, value, fill = factor(dose_mg_per_kg)),
    shape = 21,
    size = 3.0,
    stroke = 0.45,
    color = "black",
    position = ggplot2::position_jitter(
      width = 0.10, height = 0, seed = plot_seed + 62L
    )
  ) +
  ggplot2::geom_crossbar(
    data = cbs_stage_means,
    ggplot2::aes(x_position, value, ymin = value, ymax = value),
    width = 0.50,
    linewidth = 0.6,
    color = "black",
    inherit.aes = FALSE
  ) +
  ggplot2::geom_text(
    data = cbs_reduction_labels,
    ggplot2::aes(x, y, label = label),
    size = 2.5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  ggplot2::facet_wrap(~initial_ploidy, nrow = 1) +
  ggplot2::scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Injected-cell\nreference", "Endpoint\ntumor"),
    limits = c(0.72, 2.28)
  ) +
  ggplot2::scale_fill_manual(
    values = cbs_mouse_dose_colors,
    name = "Dose (mg/kg)"
  ) +
  ggplot2::labs(
    title = "Injected-reference to endpoint ploidy",
    subtitle = paste(
      "Reference crosses: project-designated proxy cells",
      "(chr999 added in haploid-genome-equivalent units);",
      "endpoint circles: one mean per mouse; descriptive only"
    ),
    x = NULL,
    y = "Reference ploidy / endpoint copy-number score"
  ) +
  shared_context_figure_theme(9) +
  ggplot2::theme(
    legend.position = "bottom",
    plot.subtitle = ggplot2::element_text(size = 7.2, lineheight = 1.05),
    axis.text.x = ggplot2::element_text(size = 7.2),
    strip.text = ggplot2::element_text(face = "bold")
  )
s6e_mouse <- shared_context_add_tag(s6e_mouse, "E")
s6_top <- patchwork::wrap_plots(
  patchwork::wrap_plots(s6a, s6b, s6c, ncol = 1),
  s6d,
  ncol = 2,
  widths = c(1, 3)
)
s6 <- patchwork::wrap_plots(
  s6_top,
  s6e_mouse,
  ncol = 1,
  heights = c(1.15, 0.85)
) + patchwork::plot_annotation(
  title = paste(
    "Supplementary Figure 6 | Endpoint tumor ploidy and",
    "NUMBAT-derived copy-number states"
  )
)
panel_rows[[length(panel_rows) + 1L]] <- save_composite(
  s6,
  figure_dir,
  "SuppFig6",
  "panel_SuppFig6_composite",
  20,
  18
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
shared_si7 <- shared_context_build_si7_heatmap_panels(
  ora_matrix = ora_matrix,
  gsea_matrix = gsea_matrix
)
s7a <- shared_si7$plots$ora
s7b <- shared_si7$plots$gsea
s7 <- patchwork::wrap_plots(s7a, ncol = 1) +
  patchwork::plot_annotation(
    title = "Supplementary Figure 7 | Cluster Hallmark over-representation analysis"
  )
panel_rows[[length(panel_rows) + 1L]] <- save_composite(
  s7,
  figure_dir,
  "SuppFig7",
  "panel_SuppFig7_composite",
  10,
  8.8
)

composition_results <- list(
  SI4E = s4e_result,
  SI4G = s4g_result,
  SI5F = s5f_result,
  SI5G = s5g_result,
  SI5H = s5h_result,
  SI5I = s5i_result
)
composition_test_audit <- do.call(rbind, lapply(
  names(composition_results),
  function(panel_id) {
    data <- composition_results[[panel_id]]$tests
    data$panel_id <- panel_id
    data[, c(
      "panel_id", "group_value", "cluster", "n_group_samples",
      "n_rest_samples", "mean_group_proportion", "mean_rest_proportion",
      "difference", "exact_permutations", "p_value", "q_value", "enriched",
      "significance", "testable", "test", "strata", "adjustment"
    )]
  }
))
rownames(composition_test_audit) <- NULL
write_tsv(
  composition_test_audit,
  file.path(metadata_dir, "normalized_composition_enrichment_tests.tsv")
)

composition_plot_audit <- do.call(rbind, lapply(
  names(composition_results),
  function(panel_id) {
    data <- composition_results[[panel_id]]$plot_data
    data$panel_id <- panel_id
    data.frame(
      panel_id = data$panel_id,
      group_value = as.character(data$group_value),
      cluster = as.character(data$cluster),
      x_value = as.character(data$x_value),
      n_samples = data$n_samples,
      sum_n_cells = data$sum_n_cells,
      mean_proportion = data$mean_proportion,
      min_sample_proportion = data$min_sample_proportion,
      max_sample_proportion = data$max_sample_proportion,
      difference = data$difference,
      q_value = data$q_value,
      significance = data$significance,
      testable = data$testable,
      stringsAsFactors = FALSE
    )
  }
))
rownames(composition_plot_audit) <- NULL
write_tsv(
  composition_plot_audit,
  file.path(metadata_dir, "normalized_composition_plot_data.tsv")
)

panel_contract <- do.call(rbind, panel_rows)
rownames(panel_contract) <- NULL
displayed_panel_contract <- data.frame(
  figure_id = c("SuppFig4", "SuppFig5", "SuppFig6", "SuppFig7"),
  panel_ids = c("D,F,I", "E,G", "A,B,C,D,E", "A"),
  selection_policy = rep("panels_cited_in_manuscript_results", 4L),
  stringsAsFactors = FALSE
)
write_tsv(
  displayed_panel_contract,
  file.path(metadata_dir, "displayed_panel_contract.tsv")
)
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

cbs_input_paths <- file.path(
  cbs_dir,
  cbs_collection$reviewed_manifest$filename
)
injected_reference_input_paths <- file.path(
  injected_reference_dir,
  cbs_injected_references$manifest$filename
)
input_paths <- c(
  config_path,
  composition_helper_path,
  shared_context_helper_path,
  copy_number_helper_path,
  figure7_shared_paths,
  density_localization_config_path,
  cellcycle_pseudotime_path,
  weighted_ploidy_path,
  all_ploidy_path,
  cbs_manifest_path,
  cbs_input_paths,
  injected_reference_manifest_path,
  injected_reference_input_paths,
  file.path(cache_dir, "manifest.tsv"),
  source_paths
)
input_manifest <- data.frame(
  role = c(
    "figure7_config",
    "normalized_composition_helper",
    "shared_context_panels_helper",
    "copy_number_heatmap_helper",
    "figure7_common_io_helper",
    "figure7_tgi_statistics_helper",
    "figure7_tgi_panels_helper",
    "density_localization_config",
    "si4i_cellcycle_pseudotime",
    "endpoint_ploidy_derivation_helper",
    "endpoint_ploidy_source",
    "numbat_cbs_manifest",
    paste0(
      "numbat_cbs_matrix_",
      gsub("[^A-Za-z0-9]+", "_", basename(cbs_input_paths))
    ),
    "injected_cell_reference_manifest",
    paste0(
      "injected_cell_reference_matrix_",
      gsub(
        "[^A-Za-z0-9]+", "_", basename(injected_reference_input_paths)
      )
    ),
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
    "displayed_panel_sets",
    "figure_file_count",
    "table_mode",
    "cache_file_count",
    "plot_shuffle_seed",
    "cluster_order",
    "composition_normalization",
    "composition_enrichment_test",
    "composition_permutation_strata",
    "composition_multiple_testing",
    "composition_fdr_threshold",
    "si4i_cell_universe",
    "si4i_n_cells",
    "si4i_n_mice",
    "si4i_pointwise_positive_interval",
    "si4i_simultaneous_positive_interval",
    "si4i_global_max_abs_t_p_two_sided",
    "figure7j_copy_number_source",
    "figure7j_column_statistic",
    "figure7j_cbs_matrix_count",
    "figure7j_source_cell_count",
    "figure7j_qc_passed_tumor_cell_count",
    "figure7j_treated_tumor_cell_count",
    "figure7j_qc_selection_policy",
    "figure7j_chromosome_count",
    "figure7j_row_order",
    "figure7j_column_order",
    "si6_postprocessed_copy_number_score_unit",
    "si6_injected_reference_ploidy_policy",
    "si6_injected_reference_cell_counts",
    "si6_injected_reference_source_repository",
    "si6_injected_reference_source_commit",
    "si6_injected_reference_chr999_unit",
    "si6_injected_reference_chr999_interpretation_basis",
    "si6_injected_reference_designation_basis",
    "si6_endpoint_summary_analysis_type",
    "si6_2n_reference_assigned_autosomal_mean_ploidy",
    "si6_2n_reference_mean_chr999_extra_dna_haploid_genome_equivalents",
    "si6_2n_reference_mean_ploidy",
    "si6_2n_endpoint_mouse_balanced_mean_ploidy",
    "si6_2n_absolute_change",
    "si6_2n_relative_change_percent",
    "si6_4n_reference_assigned_autosomal_mean_ploidy",
    "si6_4n_reference_mean_chr999_extra_dna_haploid_genome_equivalents",
    "si6_4n_reference_mean_ploidy",
    "si6_4n_endpoint_mouse_balanced_mean_ploidy",
    "si6_4n_absolute_change",
    "si6_4n_relative_change_percent",
    "si6_reference_4n_minus_2n_mean_ploidy",
    "si6_endpoint_4n_minus_2n_mouse_balanced_mean_ploidy",
    "si6_absolute_separation_change",
    "si6_separation_contraction_percent",
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
    paste(
      paste0(
        displayed_panel_contract$figure_id,
        "=",
        displayed_panel_contract$panel_ids
      ),
      collapse = ";"
    ),
    as.character(nrow(panel_contract)),
    if (allow_generated_human_only_si7) {
      "run_scoped_generated_human_only_plot_tables"
    } else {
      "frozen_plot_tables_only"
    },
    as.character(length(cache_files)),
    as.character(plot_seed),
    paste(cluster_levels, collapse = ","),
    "within-sample cluster proportions averaged with equal sample weights within group",
    "exact independent-sample label permutation; one group versus exchangeable remaining samples",
    paste(
      "SI4E=initial_ploidy;SI4G=context;SI5F=descriptive;",
      "SI5G=initial_ploidy;SI5H=initial_ploidy;SI5I=dose",
      sep = ""
    ),
    "Benjamini-Hochberg across all group-by-cluster contrasts within each panel",
    "0.05",
    as.character(localization$test$cell_universe),
    as.character(localization$test$n_cells),
    as.character(localization$test$n_samples),
    sprintf(
      "%.3f-%.3f",
      localization$test$pointwise_start,
      localization$test$pointwise_end
    ),
    sprintf(
      "%.3f-%.3f",
      localization$test$simultaneous_start,
      localization$test$simultaneous_end
    ),
    format(localization$test$global_max_t_p_two_sided, digits = 16),
    "postprocessed NUMBAT-derived cell-by-segment CBS matrices",
    paste(
      "per-cell length-weighted mean across available CBS segments within",
      "each chromosome and file-specific schema"
    ),
    as.character(length(cbs_collection$matrices)),
    as.character(nrow(cbs_collection$ploidy)),
    as.character(nrow(cbs_harmonized$matrix)),
    as.character(sum(cbs_harmonized$cell_annotations$dose_mg_per_kg > 0)),
    paste(
      "exact cells matched by the frozen endpoint-ploidy audit from the",
      "final QC-curated Seurat object; clusters 3, 4, 9, and 9c excluded"
    ),
    as.character(ncol(cbs_harmonized$matrix)),
    paste(
      "injected origin, dose, mouse, post-processed copy-number score, cell ID;",
      "no row clustering"
    ),
    "chromosomes 1-22 in genomic order; no column clustering",
    "sequenced mouse/CBS file",
    unique(cbs_injected_references$cells$ploidy_policy),
    paste(
      paste0(
        cbs_reduction_summary$injected_origin,
        "=",
        cbs_reduction_summary$n_reference_cells
      ),
      collapse = ";"
    ),
    paste(
      unique(cbs_injected_references$cells$source_repository),
      collapse = ";"
    ),
    paste(unique(cbs_injected_references$cells$source_commit), collapse = ";"),
    paste(
      unique(cbs_injected_references$cells$chr999_unit),
      collapse = ";"
    ),
    paste(
      unique(cbs_injected_references$cells$chr999_interpretation_basis),
      collapse = ";"
    ),
    paste(
      paste0(
        cbs_reduction_summary$injected_origin,
        "=",
        cbs_reduction_summary$designation_basis
      ),
      collapse = ";"
    ),
    "descriptive_only",
    format(
      cbs_2n_reduction$reference_mean_assigned_autosomal_ploidy,
      digits = 16
    ),
    format(
      cbs_2n_reduction[[
        "reference_mean_chr999_extra_dna_haploid_genome_equivalents"
      ]],
      digits = 16
    ),
    format(cbs_2n_reduction$reference_mean_ploidy, digits = 16),
    format(
      cbs_2n_reduction$endpoint_mouse_balanced_mean_ploidy,
      digits = 16
    ),
    format(cbs_2n_reduction$absolute_change, digits = 16),
    format(cbs_2n_reduction$relative_change_percent, digits = 16),
    format(
      cbs_4n_reduction$reference_mean_assigned_autosomal_ploidy,
      digits = 16
    ),
    format(
      cbs_4n_reduction[[
        "reference_mean_chr999_extra_dna_haploid_genome_equivalents"
      ]],
      digits = 16
    ),
    format(cbs_4n_reduction$reference_mean_ploidy, digits = 16),
    format(
      cbs_4n_reduction$endpoint_mouse_balanced_mean_ploidy,
      digits = 16
    ),
    format(cbs_4n_reduction$absolute_change, digits = 16),
    format(cbs_4n_reduction$relative_change_percent, digits = 16),
    format(
      cbs_separation_summary$reference_4n_minus_2n_mean_ploidy,
      digits = 16
    ),
    format(
      cbs_separation_summary$endpoint_4n_minus_2n_mouse_balanced_mean_ploidy,
      digits = 16
    ),
    format(cbs_separation_summary$absolute_separation_change, digits = 16),
    format(
      cbs_separation_summary$separation_contraction_percent,
      digits = 16
    ),
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
    "selected_cellcycle_tumor_cells",
    "si4i_localization_cells",
    "si4i_localization_mice",
    "si4i_localization_permutations",
    "normalized_composition_panels",
    "normalized_composition_contrasts",
    "normalized_composition_tested_contrasts",
    "normalized_composition_descriptive_contrasts",
    "normalized_composition_significant_enrichments",
    "figure7j_cbs_matrices",
    "figure7j_source_cells",
    "figure7j_qc_passed_tumor_cells",
    "figure7j_qc_passed_treated_tumor_cells",
    "figure7j_chromosomes",
    "figure7j_missing_chromosome_mean_fraction",
    "figure7j_cbs_schemas",
    "si6_postprocessed_copy_number_score_mice",
    "si6_injected_reference_cells",
    "si6_endpoint_summary_analysis_type",
    "si6_2n_reference_mean_ploidy",
    "si6_2n_endpoint_mouse_balanced_mean_ploidy",
    "si6_2n_absolute_change",
    "si6_2n_relative_change_percent",
    "si6_2n_all_endpoint_cells_below_reference_min",
    "si6_4n_reference_mean_ploidy",
    "si6_4n_endpoint_mouse_balanced_mean_ploidy",
    "si6_4n_absolute_change",
    "si6_4n_relative_change_percent",
    "si6_4n_all_endpoint_cells_below_reference_min",
    "si6_separation_contraction_percent",
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
    selected_tumor_cells,
    localization$test$n_cells,
    localization$test$n_samples,
    localization$test$n_permutations,
    length(composition_results),
    nrow(composition_test_audit),
    sum(composition_test_audit$testable),
    sum(!composition_test_audit$testable),
    sum(composition_test_audit$enriched),
    length(cbs_collection$matrices),
    nrow(cbs_collection$ploidy),
    nrow(cbs_harmonized$matrix),
    sum(cbs_harmonized$cell_annotations$dose_mg_per_kg > 0),
    ncol(cbs_harmonized$matrix),
    format(mean(is.na(cbs_harmonized$matrix)), digits = 16),
    length(unique(cbs_harmonized$chromosome_schema_audit$schema_id)),
    nrow(cbs_endpoint_summaries$samples),
    nrow(cbs_injected_references$cells),
    "descriptive_only",
    format(cbs_2n_reduction$reference_mean_ploidy, digits = 16),
    format(
      cbs_2n_reduction$endpoint_mouse_balanced_mean_ploidy,
      digits = 16
    ),
    format(cbs_2n_reduction$absolute_change, digits = 16),
    format(cbs_2n_reduction$relative_change_percent, digits = 16),
    tolower(as.character(
      cbs_2n_reduction$all_endpoint_cells_below_reference_min
    )),
    format(cbs_4n_reduction$reference_mean_ploidy, digits = 16),
    format(
      cbs_4n_reduction$endpoint_mouse_balanced_mean_ploidy,
      digits = 16
    ),
    format(cbs_4n_reduction$absolute_change, digits = 16),
    format(cbs_4n_reduction$relative_change_percent, digits = 16),
    tolower(as.character(
      cbs_4n_reduction$all_endpoint_cells_below_reference_min
    )),
    format(
      cbs_separation_summary$separation_contraction_percent,
      digits = 16
    ),
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
    "normalized_composition_helper",
    "normalized_composition_helper_sha256",
    "shared_context_panels_helper",
    "shared_context_panels_helper_sha256",
    "si4i_cellcycle_pseudotime",
    "si4i_cellcycle_pseudotime_sha256",
    "si4i_density_localization_config",
    "si4i_density_localization_config_sha256",
    "si4i_density_localization_method",
    "copy_number_heatmap_helper",
    "copy_number_heatmap_helper_sha256",
    "endpoint_ploidy_derivation_helper",
    "endpoint_ploidy_derivation_helper_sha256",
    "endpoint_ploidy_source",
    "endpoint_ploidy_source_sha256",
    "numbat_cbs_manifest",
    "numbat_cbs_manifest_sha256",
    "numbat_cbs_matrix_count",
    "numbat_cbs_matrix_hashes",
    "figure7j_qc_selection",
    "injected_cell_reference_manifest",
    "injected_cell_reference_manifest_sha256",
    "injected_cell_reference_matrix_hashes",
    "figure7j_harmonization",
    "si6e_reference_source",
    "si6e_reference_chr999_interpretation",
    "si6e_reference_designation_basis",
    "si6e_reference_ploidy_policy",
    "si6e_summary_analysis_type",
    "si6e_ploidy_reduction_comparison",
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
    repo_relative(composition_helper_path, root),
    file_sha256(composition_helper_path),
    repo_relative(shared_context_helper_path, root),
    file_sha256(shared_context_helper_path),
    repo_relative(cellcycle_pseudotime_path, root),
    file_sha256(cellcycle_pseudotime_path),
    repo_relative(density_localization_config_path, root),
    file_sha256(density_localization_config_path),
    paste(
      "equal-mouse Gaussian-kernel treated-minus-vehicle density contrast;",
      "exact injected-origin-stratified pointwise and studentized max-|T|",
      "permutation inference on the reviewed 2,881-cell subset"
    ),
    repo_relative(copy_number_helper_path, root),
    file_sha256(copy_number_helper_path),
    repo_relative(weighted_ploidy_path, root),
    file_sha256(weighted_ploidy_path),
    repo_relative(all_ploidy_path, root),
    file_sha256(all_ploidy_path),
    repo_relative(cbs_manifest_path, root),
    file_sha256(cbs_manifest_path),
    as.character(length(cbs_input_paths)),
    paste(
      paste0(
        basename(cbs_input_paths),
        "=",
        vapply(cbs_input_paths, file_sha256, character(1L))
      ),
      collapse = ";"
    ),
    paste(
      "the complete 14,125-cell CBS source is checksum/value validated,",
      "then restricted by the frozen endpoint-ploidy audit to the exact",
      "9,832 QC-passed tumor cells (5,335 treated); clusters 3, 4, 9,",
      "and 9c remain excluded"
    ),
    repo_relative(injected_reference_manifest_path, root),
    file_sha256(injected_reference_manifest_path),
    paste(
      paste0(
        basename(injected_reference_input_paths),
        "=",
        vapply(
          injected_reference_input_paths,
          file_sha256,
          character(1L)
        )
      ),
      collapse = ";"
    ),
    paste(
      "each separately postprocessed file-specific CBS schema reduced",
      "independently to chr1-22 finite-segment length-weighted means;",
      "exported and finite coverage audited separately; no coordinate",
      "alignment; rows and columns unclustered"
    ),
    paste(
      unique(cbs_injected_references$cells$source_repository),
      paste(unique(cbs_injected_references$cells$source_commit), collapse = ";"),
      sep = "@"
    ),
    paste(
      paste(unique(cbs_injected_references$cells$chr999_unit), collapse = ";"),
      paste(
        unique(
          cbs_injected_references$cells$chr999_interpretation_basis
        ),
        collapse = ";"
      ),
      sep = "; "
    ),
    paste(
      paste0(
        cbs_reduction_summary$injected_origin,
        "=",
        cbs_reduction_summary$designation_basis
      ),
      collapse = ";"
    ),
    unique(cbs_injected_references$cells$ploidy_policy),
    "descriptive_only; no endpoint cross-origin or reference-to-endpoint test",
    paste(
      "project-designated lineage-matched 2N-A7M/4N-A5M karyotype",
      "reference distributions, including the chr999 unassigned-extra-DNA",
      "haploid-genome-equivalent term added to autosomal ploidy, compared",
      "descriptively with one postprocessed endpoint",
      "mean per mouse; no formal P value because each reference is one",
      "culture-level biological unit and origin-specific endpoint runs use",
      "different schemas/calibration"
    ),
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
        "Reviewed SI7 matrices are the exact outputs from raw-refit run",
        "grch_human_only_v2_20260729_raw_refit_retry3_si_figures:",
        "exact GRCh38-prefixed RNA counts were retained before fresh RNA",
        "normalization, differential expression, symbol cleanup,",
        "deduplication, ORA, and GSEA."
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
    paste0(
      "Figure 7J CBS matrices/source cells/QC-passed cells/treated cells/",
      "chromosomes: ",
      length(cbs_collection$matrices),
      "/",
      nrow(cbs_collection$ploidy),
      "/",
      nrow(cbs_harmonized$matrix),
      "/",
      sum(cbs_harmonized$cell_annotations$dose_mg_per_kg > 0),
      "/",
      ncol(cbs_harmonized$matrix)
    ),
    paste0(
      "SI6 2N injected-reference to endpoint descriptive change: ",
      format(cbs_2n_reduction$absolute_change, digits = 9),
      " (",
      format(cbs_2n_reduction$relative_change_percent, digits = 9),
      "%)"
    ),
    paste0(
      "SI6 4N injected-reference to endpoint descriptive change: ",
      format(cbs_4n_reduction$absolute_change, digits = 9),
      " (",
      format(cbs_4n_reduction$relative_change_percent, digits = 9),
      "%)"
    ),
    paste0(
      "SI6 descriptive 4N-minus-2N separation contraction: ",
      format(cbs_separation_summary$separation_contraction_percent,
             digits = 9),
      "%"
    ),
    paste0("Tumor mice: ", n_mice),
    paste0("Selected CellCycle tumor cells: ", selected_tumor_cells),
    paste0(
      "SI4I pointwise/simultaneous positive intervals and global P: ",
      sprintf(
        "%.3f-%.3f / %.3f-%.3f / %.6g",
        localization$test$pointwise_start,
        localization$test$pointwise_end,
        localization$test$simultaneous_start,
        localization$test$simultaneous_end,
        localization$test$global_max_t_p_two_sided
      )
    ),
    paste0(
      "Normalized composition panels: ",
      paste(names(composition_results), collapse = ",")
    ),
    paste0(
      "Significant composition enrichments (BH FDR <= 0.05): ",
      sum(composition_test_audit$enriched)
    ),
    paste0(
      "Tested/descriptive composition contrasts: ",
      sum(composition_test_audit$testable),
      "/",
      sum(!composition_test_audit$testable)
    ),
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
