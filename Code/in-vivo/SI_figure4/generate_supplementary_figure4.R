#!/usr/bin/env Rscript

# Supplementary Figure 4
#
# Descriptive UMAP and cluster-composition views for the in-vivo tumor cells.
# This script is intentionally downstream of the Figure 7 input workflow:
#   * seurat_metadata.csv supplies the complete Seurat cell universe and UMAP.
#   * scvelo_cell_metrics.csv supplies/validates sample-level ploidy and context.
#
# The main figure uses all tumor cells in seurat_metadata.csv. It does not
# restrict composition estimates to cells retained by scVelo.

parse_cli_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (grepl("^--[^=]+=", token)) {
      key <- sub("^--([^=]+)=.*$", "\\1", token)
      out[[key]] <- sub("^--[^=]+=", "", token)
      i <- i + 1L
    } else if (grepl("^--", token)) {
      key <- sub("^--", "", token)
      if (i < length(args) && !grepl("^--", args[[i + 1L]])) {
        out[[key]] <- args[[i + 1L]]
        i <- i + 2L
      } else {
        out[[key]] <- "TRUE"
        i <- i + 1L
      }
    } else {
      stop("Unexpected positional argument: ", token, call. = FALSE)
    }
  }
  out
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
      "  Rscript Code/in-vivo/SI_figure4/generate_supplementary_figure4.R \\",
      "    --seurat-metadata Data/in-vivo/seurat_metadata.csv \\",
      "    --scvelo-metrics Data/in-vivo/scvelo_cell_metrics.csv \\",
      "    --output-dir Results/in-vivo/SI_figure4",
      "",
      "Options:",
      "  --overwrite       Replace this script's known outputs in a nonempty output directory.",
      "  --help            Show this help message.",
      sep = "\n"
    ),
    "\n"
  )
}

resolve_script_path <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command, value = TRUE)
  if (!length(file_arg)) return(NA_character_)
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
}

is_absolute_path <- function(path) grepl("^/", path)

resolve_path <- function(path, repo_root, must_work = FALSE) {
  resolved <- if (is_absolute_path(path)) path else file.path(repo_root, path)
  normalizePath(resolved, mustWork = must_work)
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit)) hit[[1L]] else paths[[1L]]
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Missing required R package: ", package, call. = FALSE)
  }
}

read_csv_checked <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
  if (file.info(path)$size <= 0) stop(label, " is empty: ", path, call. = FALSE)
  out <- utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN")
  )
  if (!nrow(out)) stop(label, " has no data rows: ", path, call. = FALSE)
  out
}

assert_single_column <- function(df, column, label) {
  n <- sum(names(df) == column)
  if (n != 1L) {
    stop(label, " must contain exactly one `", column, "` column; found ", n, call. = FALSE)
  }
  invisible(column)
}

resolve_column <- function(df, candidates, label, required = TRUE) {
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  exact <- candidates[candidates %in% names(df)]
  if (length(exact)) {
    assert_single_column(df, exact[[1L]], label)
    return(exact[[1L]])
  }
  lower_names <- tolower(names(df))
  for (candidate in candidates) {
    index <- which(lower_names == tolower(candidate))
    if (length(index) == 1L) return(names(df)[index])
    if (length(index) > 1L) {
      stop(label, " has ambiguous case-insensitive matches for `", candidate, "`", call. = FALSE)
    }
  }
  if (required) {
    stop(label, " is missing a required column; tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  NA_character_
}

clean_character <- function(values) {
  out <- trimws(enc2utf8(as.character(values)))
  out[is.na(values) | out == "" | out %in% c("NA", "NaN", "None")] <- NA_character_
  out
}

numeric_strict <- function(values, label) {
  raw <- clean_character(values)
  out <- suppressWarnings(as.numeric(raw))
  bad <- !is.na(raw) & !is.finite(out)
  if (any(bad)) {
    stop(label, " contains nonnumeric/nonfinite values, including: ",
         paste(head(unique(raw[bad]), 10L), collapse = ", "), call. = FALSE)
  }
  out
}

standardize_dose <- function(values, label = "Dose") {
  raw <- clean_character(values)
  normalized <- tolower(gsub("[[:space:]]+", "", raw))
  out <- rep(NA_character_, length(raw))
  out[normalized %in% c("0", "0mg", "0mg/kg")] <- "0mg/kg"
  out[normalized %in% c("30", "30mg", "30mg/kg")] <- "30mg/kg"
  out[normalized %in% c("120", "120mg", "120mg/kg")] <- "120mg/kg"
  unknown <- !is.na(raw) & is.na(out)
  if (any(unknown)) {
    stop(label, " contains unexpected values: ",
         paste(sort(unique(raw[unknown])), collapse = ", "), call. = FALSE)
  }
  out
}

standardize_ploidy <- function(values) {
  raw <- clean_character(values)
  upper <- toupper(raw)
  out <- rep(NA_character_, length(raw))
  out[!is.na(upper) & grepl("2N", upper, fixed = TRUE)] <- "2N"
  out[is.na(out) & !is.na(upper) & grepl("4N", upper, fixed = TRUE)] <- "4N"
  out
}

standardize_context <- function(values) {
  raw <- tolower(clean_character(values))
  out <- rep(NA_character_, length(raw))
  out[!is.na(raw) & grepl("cell[-_ ]?(line|culture)|cellline", raw)] <- "CellLine"
  out[is.na(out) & !is.na(raw) & grepl("tumou?r", raw)] <- "Tumor"
  out
}

sort_cluster_levels <- function(values) {
  values <- sort(unique(clean_character(values)))
  numeric_prefix <- suppressWarnings(as.numeric(sub("^([0-9]+).*$", "\\1", values)))
  has_prefix <- grepl("^[0-9]+", values)
  order_key <- ifelse(has_prefix, numeric_prefix, Inf)
  values[order(order_key, values, na.last = TRUE, method = "radix")]
}

mapping_by_group <- function(group, value, value_label) {
  group <- clean_character(group)
  value <- clean_character(value)
  keep <- !is.na(group) & !is.na(value)
  pieces <- split(value[keep], group[keep])
  conflicts <- names(pieces)[vapply(pieces, function(x) length(unique(x)) > 1L, logical(1))]
  if (length(conflicts)) {
    stop(
      "Each mouse must map to one ", value_label, "; conflicts include: ",
      paste(head(conflicts, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  if (!length(pieces)) return(setNames(character(0), character(0)))
  setNames(vapply(pieces, function(x) unique(x)[[1L]], character(1)), names(pieces))
}

assert_consistent <- function(left, right, label) {
  conflict <- !is.na(left) & !is.na(right) & left != right
  if (any(conflict)) {
    stop(label, " disagrees for ", sum(conflict), " matched cells", call. = FALSE)
  }
  invisible(TRUE)
}

file_sha256 <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(file = path, algo = "sha256", serialize = FALSE))
  }
  sha256sum <- Sys.which("sha256sum")
  if (nzchar(sha256sum)) {
    result <- system2(sha256sum, shQuote(path), stdout = TRUE, stderr = TRUE)
  } else {
    shasum <- Sys.which("shasum")
    if (!nzchar(shasum)) return(NA_character_)
    result <- system2(shasum, c("-a", "256", shQuote(path)), stdout = TRUE, stderr = TRUE)
  }
  if (!length(result)) return(NA_character_)
  strsplit(trimws(result[[1L]]), "[[:space:]]+")[[1L]][[1L]]
}

write_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, na = "", quote = TRUE)
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Failed to write output table: ", path, call. = FALSE)
  }
  invisible(path)
}

write_tsv <- function(df, path) {
  utils::write.table(
    df,
    path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = ""
  )
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Failed to write output table: ", path, call. = FALSE)
  }
  invisible(path)
}

figure_theme <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", color = "#222222", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(color = "#444444", size = base_size - 1),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 3),
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey94", color = "grey75", linewidth = 0.35),
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

umap_color_guide <- function() {
  ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1, stroke = 0))
}

add_panel_tag <- function(plot, tag) {
  plot + ggplot2::labs(tag = tag)
}

save_pdf_png <- function(plot, stub, width, height, dpi = 300) {
  pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf"
  ggplot2::ggsave(
    paste0(stub, ".pdf"), plot = plot, device = pdf_device,
    width = width, height = height, units = "in", bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    paste0(stub, ".png"), plot = plot, device = "png", dpi = dpi,
    width = width, height = height, units = "in", bg = "white", limitsize = FALSE
  )
  invisible(c(paste0(stub, ".pdf"), paste0(stub, ".png")))
}

percent_labels <- function(values) paste0(round(100 * values), "%")

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
if (arg_flag(args, "help", FALSE)) {
  usage()
  quit(save = "no", status = 0L)
}

require_package("ggplot2")
require_package("patchwork")

script_path <- resolve_script_path()
if (is.na(script_path)) stop("Cannot resolve the script path", call. = FALSE)
script_dir <- dirname(script_path)
repo_root <- normalizePath(dirname(dirname(dirname(script_dir))), mustWork = TRUE)

seurat_metadata_path <- resolve_path(
  arg_value(args, "seurat-metadata", file.path("Data", "in-vivo", "seurat_metadata.csv")),
  repo_root,
  must_work = FALSE
)

default_scvelo <- first_existing(c(
  file.path(repo_root, "Data", "in-vivo", "scvelo_cell_metrics.csv"),
  file.path(repo_root, "Results", "in-vivo", "figure7", "intermediates", "scvelo_cell_metrics.csv")
))
scvelo_metrics_path <- resolve_path(
  arg_value(args, "scvelo-metrics", default_scvelo),
  repo_root,
  must_work = FALSE
)

output_dir <- resolve_path(
  arg_value(args, "output-dir", file.path("Results", "in-vivo", "SI_figure4")),
  repo_root,
  must_work = FALSE
)
overwrite <- arg_flag(args, "overwrite", FALSE)

if (dir.exists(output_dir)) {
  existing <- list.files(output_dir, all.files = TRUE, no.. = TRUE)
  manager_bootstrap <- c("logs", "metadata")
  blocking <- setdiff(existing, manager_bootstrap)
  metadata_existing <- list.files(file.path(output_dir, "metadata"), all.files = TRUE, no.. = TRUE)
  if ((length(blocking) || length(metadata_existing)) && !overwrite) {
    stop("Output directory is not empty; pass --overwrite to replace known outputs: ", output_dir, call. = FALSE)
  }
}

figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
metadata_dir <- file.path(output_dir, "metadata")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading Seurat metadata: ", seurat_metadata_path)
seurat <- read_csv_checked(seurat_metadata_path, "Seurat metadata")
message("Reading scVelo cell metrics: ", scvelo_metrics_path)
metrics <- read_csv_checked(scvelo_metrics_path, "scVelo cell metrics")

for (column in c("cell", "UMAP_1", "UMAP_2", "Dose")) {
  assert_single_column(seurat, column, "Seurat metadata")
}
assert_single_column(metrics, "cell", "scVelo cell metrics")

canonical_cluster_col <- resolve_column(
  seurat,
  c("cluster_final", "clusters"),
  "Seurat canonical cluster"
)

seurat$cell <- clean_character(seurat$cell)
metrics$cell <- clean_character(metrics$cell)
if (anyNA(seurat$cell) || anyNA(metrics$cell)) stop("Input cell identifiers cannot be missing", call. = FALSE)
if (anyDuplicated(seurat$cell)) stop("Seurat metadata contains duplicated cell identifiers", call. = FALSE)
if (anyDuplicated(metrics$cell)) stop("scVelo metrics contains duplicated cell identifiers", call. = FALSE)

missing_metrics_cells <- setdiff(metrics$cell, seurat$cell)
if (length(missing_metrics_cells)) {
  stop(
    "scVelo metrics contains cells absent from Seurat metadata, including: ",
    paste(head(missing_metrics_cells, 10L), collapse = ", "),
    call. = FALSE
  )
}

sample_col <- resolve_column(
  seurat,
  c("sample", "IDs", "ID", "orig.ident", "Sequencing.IDs"),
  "Seurat metadata sample identity"
)
metric_ploidy_col <- resolve_column(metrics, c("Ploidy", "ploidy"), "scVelo ploidy")
metric_context_col <- resolve_column(
  metrics,
  c("TN", "trajectory_context", "trajectory_group", "sample_type"),
  "scVelo tumor/cell-line context"
)
metric_dose_col <- resolve_column(metrics, c("Dose", "Dose_DEG"), "scVelo dose")

seurat$UMAP_1 <- numeric_strict(seurat$UMAP_1, "Seurat UMAP_1")
seurat$UMAP_2 <- numeric_strict(seurat$UMAP_2, "Seurat UMAP_2")
if (any(!is.finite(seurat$UMAP_1)) || any(!is.finite(seurat$UMAP_2))) {
  stop("Seurat metadata contains missing or nonfinite UMAP coordinates", call. = FALSE)
}

seurat$mouse <- clean_character(seurat[[sample_col]])
seurat$cluster <- clean_character(seurat[[canonical_cluster_col]])
seurat$dose <- standardize_dose(seurat$Dose, "Seurat Dose")
if (anyNA(seurat$mouse)) stop("Seurat metadata contains missing mouse/sample identities", call. = FALSE)
if (anyNA(seurat$cluster)) stop("Seurat metadata contains missing canonical cluster labels", call. = FALSE)

metric_index <- match(seurat$cell, metrics$cell)
seurat$in_scvelo <- !is.na(metric_index)
metrics_ploidy <- standardize_ploidy(metrics[[metric_ploidy_col]])
metrics_context <- standardize_context(metrics[[metric_context_col]])
metrics_dose <- standardize_dose(metrics[[metric_dose_col]], "scVelo Dose")
if (anyNA(metrics_ploidy)) stop("scVelo metrics contains unresolved Ploidy values", call. = FALSE)
if (anyNA(metrics_context)) stop("scVelo metrics contains unresolved TN/context values", call. = FALSE)
if (any(is.na(metrics_dose) & metrics_context == "Tumor")) {
  stop("scVelo tumor cells contain missing Dose values", call. = FALSE)
}

matched_ploidy <- rep(NA_character_, nrow(seurat))
matched_context <- rep(NA_character_, nrow(seurat))
matched_dose <- rep(NA_character_, nrow(seurat))
matched_ploidy[seurat$in_scvelo] <- metrics_ploidy[metric_index[seurat$in_scvelo]]
matched_context[seurat$in_scvelo] <- metrics_context[metric_index[seurat$in_scvelo]]
matched_dose[seurat$in_scvelo] <- metrics_dose[metric_index[seurat$in_scvelo]]
assert_consistent(seurat$dose, matched_dose, "Seurat/scVelo Dose")

metric_cluster_col <- resolve_column(
  metrics,
  canonical_cluster_col,
  "scVelo canonical cluster"
)
matched_cluster <- rep(NA_character_, nrow(seurat))
matched_cluster[seurat$in_scvelo] <- clean_character(metrics[[metric_cluster_col]])[metric_index[seurat$in_scvelo]]
assert_consistent(seurat$cluster, matched_cluster, "Seurat/scVelo canonical cluster")

sample_ploidy <- mapping_by_group(seurat$mouse, matched_ploidy, "initial ploidy")
sample_context <- mapping_by_group(seurat$mouse, matched_context, "tumor/cell-line context")

raw_ploidy_col <- resolve_column(
  seurat,
  c("Ploidy", "initial_ploidy", "sample_type"),
  "Seurat raw ploidy",
  required = FALSE
)
raw_ploidy <- if (is.na(raw_ploidy_col)) rep(NA_character_, nrow(seurat)) else standardize_ploidy(seurat[[raw_ploidy_col]])
raw_ploidy[is.na(raw_ploidy)] <- standardize_ploidy(seurat$mouse)[is.na(raw_ploidy)]

raw_context_col <- resolve_column(
  seurat,
  c("TN", "trajectory_context", "sample_type"),
  "Seurat raw context",
  required = FALSE
)
raw_context <- if (is.na(raw_context_col)) rep(NA_character_, nrow(seurat)) else standardize_context(seurat[[raw_context_col]])
raw_context[is.na(raw_context) & grepl("cell[-_ ]?culture", seurat$mouse, ignore.case = TRUE)] <- "CellLine"

seurat$initial_ploidy <- unname(sample_ploidy[seurat$mouse])
seurat$context <- unname(sample_context[seurat$mouse])
assert_consistent(seurat$initial_ploidy, raw_ploidy, "Metric-derived/raw initial ploidy")
assert_consistent(seurat$context, raw_context, "Metric-derived/raw context")
seurat$initial_ploidy[is.na(seurat$initial_ploidy)] <- raw_ploidy[is.na(seurat$initial_ploidy)]
seurat$context[is.na(seurat$context)] <- raw_context[is.na(seurat$context)]
seurat$context[is.na(seurat$context) & !is.na(seurat$initial_ploidy)] <- "Tumor"

if (anyNA(seurat$context)) {
  stop("Could not resolve tumor/cell-line context for every Seurat cell", call. = FALSE)
}

tumor <- seurat[seurat$context == "Tumor", , drop = FALSE]
if (!nrow(tumor)) stop("No tumor cells remain after context filtering", call. = FALSE)
if (anyNA(tumor$initial_ploidy)) stop("Tumor cells contain unresolved initial ploidy", call. = FALSE)
if (any(!tumor$initial_ploidy %in% c("2N", "4N"))) stop("Unexpected tumor ploidy value", call. = FALSE)
if (anyNA(tumor$dose)) stop("Tumor cells contain missing Dose values", call. = FALSE)

mouse_dose <- mapping_by_group(tumor$mouse, tumor$dose, "dose")
mouse_ploidy <- mapping_by_group(tumor$mouse, tumor$initial_ploidy, "initial ploidy")
mouse_context <- mapping_by_group(tumor$mouse, tumor$context, "context")
if (any(mouse_context != "Tumor")) stop("Tumor subset contains a non-tumor mouse", call. = FALSE)

ploidy_levels <- c("2N", "4N")
dose_levels <- c("0mg/kg", "30mg/kg", "120mg/kg")
cluster_levels <- sort_cluster_levels(tumor$cluster)
if (!length(cluster_levels)) stop("No canonical clusters remain in tumor cells", call. = FALSE)

sample_metadata <- unique(tumor[, c("mouse", "initial_ploidy", "dose"), drop = FALSE])
if (nrow(sample_metadata) != length(unique(tumor$mouse))) {
  stop("Mouse metadata are not unique after ploidy/dose resolution", call. = FALSE)
}
sample_metadata$initial_ploidy <- factor(sample_metadata$initial_ploidy, levels = ploidy_levels)
sample_metadata$dose <- factor(sample_metadata$dose, levels = dose_levels)
sample_metadata <- sample_metadata[
  order(sample_metadata$initial_ploidy, sample_metadata$dose, sample_metadata$mouse, method = "radix"),
  ,
  drop = FALSE
]
mouse_levels <- sample_metadata$mouse

tumor$mouse <- factor(tumor$mouse, levels = mouse_levels)
tumor$initial_ploidy <- factor(tumor$initial_ploidy, levels = ploidy_levels)
tumor$dose <- factor(tumor$dose, levels = dose_levels)
tumor$cluster <- factor(tumor$cluster, levels = cluster_levels)

cluster_colors <- setNames(
  grDevices::hcl.colors(max(length(cluster_levels), 3L), palette = "Dark 3")[seq_along(cluster_levels)],
  cluster_levels
)
ploidy_colors <- c("2N" = "#4C78A8", "4N" = "#E45756")
dose_colors <- c("0mg/kg" = "#666666", "30mg/kg" = "#d95f02", "120mg/kg" = "#1b9e77")

annotation_col <- resolve_column(
  tumor,
  c("cluster_final_annotation_primary"),
  "Canonical cluster annotation",
  required = FALSE
)
annotation_map <- setNames(rep(NA_character_, length(cluster_levels)), cluster_levels)
if (!is.na(annotation_col)) {
  annotation_map <- mapping_by_group(tumor$cluster, tumor[[annotation_col]], "primary annotation")
  annotation_map <- annotation_map[cluster_levels]
}

cluster_counts <- as.integer(table(factor(tumor$cluster, levels = cluster_levels)))
color_key <- data.frame(
  cluster_final = cluster_levels,
  source_field = canonical_cluster_col,
  color = unname(cluster_colors[cluster_levels]),
  annotation_primary = unname(annotation_map[cluster_levels]),
  n_tumor_cells = cluster_counts,
  stringsAsFactors = FALSE
)

composition_table <- as.data.frame(
  table(
    mouse = factor(tumor$mouse, levels = mouse_levels),
    cluster_final = factor(tumor$cluster, levels = cluster_levels)
  ),
  stringsAsFactors = FALSE
)
names(composition_table)[names(composition_table) == "Freq"] <- "n_cells"
composition_table$mouse <- as.character(composition_table$mouse)
composition_table$cluster_final <- as.character(composition_table$cluster_final)
composition_table$total_cells <- ave(composition_table$n_cells, composition_table$mouse, FUN = sum)
composition_table$proportion <- composition_table$n_cells / composition_table$total_cells
composition_table$initial_ploidy <- as.character(mouse_ploidy[composition_table$mouse])
composition_table$dose <- as.character(mouse_dose[composition_table$mouse])
composition_table <- composition_table[, c(
  "mouse", "initial_ploidy", "dose", "cluster_final", "n_cells", "total_cells", "proportion"
)]

mouse_sums <- tapply(composition_table$proportion, composition_table$mouse, sum)
if (any(abs(mouse_sums - 1) > 1e-12)) {
  stop("Cluster proportions do not sum to one within every mouse", call. = FALSE)
}

group_key <- interaction(
  composition_table$initial_ploidy,
  composition_table$dose,
  composition_table$cluster_final,
  drop = TRUE,
  lex.order = TRUE
)
group_pieces <- split(composition_table, group_key)
group_summary <- do.call(rbind, lapply(group_pieces, function(piece) {
  data.frame(
    initial_ploidy = piece$initial_ploidy[[1L]],
    dose = piece$dose[[1L]],
    cluster_final = piece$cluster_final[[1L]],
    n_mice = length(unique(piece$mouse)),
    mean_proportion = mean(piece$proportion),
    sd_proportion = if (nrow(piece) > 1L) stats::sd(piece$proportion) else NA_real_,
    min_proportion = min(piece$proportion),
    max_proportion = max(piece$proportion),
    stringsAsFactors = FALSE
  )
}))
rownames(group_summary) <- NULL
group_summary$initial_ploidy <- factor(group_summary$initial_ploidy, levels = ploidy_levels)
group_summary$dose <- factor(group_summary$dose, levels = dose_levels)
group_summary$cluster_final <- factor(group_summary$cluster_final, levels = cluster_levels)
group_summary <- group_summary[order(
  group_summary$initial_ploidy,
  group_summary$dose,
  group_summary$cluster_final
), , drop = FALSE]

group_sums <- aggregate(mean_proportion ~ initial_ploidy + dose, data = group_summary, sum)
if (any(abs(group_sums$mean_proportion - 1) > 1e-12)) {
  stop("Equal-mouse group mean proportions do not sum to one", call. = FALSE)
}

n_tumor <- nrow(tumor)
n_mice <- length(mouse_levels)
n_clusters <- length(cluster_levels)
point_size <- if (n_tumor > 50000L) 0.08 else if (n_tumor > 20000L) 0.14 else 0.24
set.seed(5826)
plot_data <- tumor[sample.int(n_tumor), , drop = FALSE]

cluster_centers <- aggregate(
  cbind(UMAP_1, UMAP_2) ~ cluster,
  data = tumor,
  FUN = stats::median
)

p_a <- ggplot2::ggplot(plot_data, ggplot2::aes(UMAP_1, UMAP_2, color = cluster)) +
  ggplot2::geom_point(size = point_size, alpha = 0.78, stroke = 0) +
  ggplot2::geom_label(
    data = cluster_centers,
    ggplot2::aes(x = UMAP_1, y = UMAP_2, label = cluster),
    inherit.aes = FALSE,
    size = 2.6,
    linewidth = 0.2,
    fill = "white",
    color = "#222222",
    alpha = 0.86,
    label.padding = grid::unit(0.12, "lines")
  ) +
  ggplot2::scale_color_manual(values = cluster_colors, drop = FALSE, name = "Cluster") +
  ggplot2::guides(color = umap_color_guide()) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "UMAP by cluster",
    subtitle = sprintf("Tumor cells; %s cells across %s clusters", format(n_tumor, big.mark = ","), n_clusters),
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  umap_theme(10)
p_a <- add_panel_tag(p_a, "A")

panel_metadata <- sample_metadata
panel_metadata$mouse_panel <- paste(
  panel_metadata$mouse,
  as.character(panel_metadata$initial_ploidy),
  as.character(panel_metadata$dose),
  sep = " | "
)
panel_levels <- panel_metadata$mouse_panel
panel_lookup <- setNames(panel_metadata$mouse_panel, panel_metadata$mouse)
mouse_plot_data <- tumor
mouse_plot_data$mouse_panel <- factor(unname(panel_lookup[as.character(mouse_plot_data$mouse)]), levels = panel_levels)
background_data <- tumor[, c("UMAP_1", "UMAP_2"), drop = FALSE]

p_b <- ggplot2::ggplot(mouse_plot_data, ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(
    data = background_data,
    color = "grey88",
    size = max(point_size * 0.70, 0.06),
    alpha = 0.30,
    stroke = 0
  ) +
  ggplot2::geom_point(color = "#1565C0", size = max(point_size * 1.45, 0.14), alpha = 0.90, stroke = 0) +
  ggplot2::facet_wrap(~mouse_panel, ncol = 4, drop = FALSE) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "UMAP by mouse of origin",
    subtitle = sprintf("Blue: focal mouse; grey: all tumor cells in the shared embedding; n = %s mice", n_mice),
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  umap_theme(8.5) +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_text(size = 7.2),
    panel.spacing = grid::unit(0.10, "lines")
  )
p_b <- add_panel_tag(p_b, "D")

p_c <- ggplot2::ggplot(plot_data, ggplot2::aes(UMAP_1, UMAP_2, color = initial_ploidy)) +
  ggplot2::geom_point(size = point_size, alpha = 0.76, stroke = 0) +
  ggplot2::scale_color_manual(values = ploidy_colors, drop = FALSE, name = "Initial ploidy") +
  ggplot2::guides(color = umap_color_guide()) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "UMAP by initial tumor ploidy",
    subtitle = "Tumor-level classification: 2N versus 4N",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  umap_theme(10)
p_c <- add_panel_tag(p_c, "B")

observed_dose_levels <- dose_levels[dose_levels %in% unique(as.character(tumor$dose))]
dose_panel_informative <- identical(observed_dose_levels, dose_levels)
if (dose_panel_informative) {
  p_d <- ggplot2::ggplot(plot_data, ggplot2::aes(UMAP_1, UMAP_2, color = dose)) +
    ggplot2::geom_point(size = point_size, alpha = 0.76, stroke = 0) +
    ggplot2::scale_color_manual(values = dose_colors, drop = FALSE, name = "Gemcitabine dose") +
    ggplot2::guides(color = umap_color_guide()) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "UMAP by treatment dose",
      subtitle = "Control, 30 mg/kg, and 120 mg/kg tumors",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    umap_theme(10)
} else {
  missing_doses <- setdiff(dose_levels, observed_dose_levels)
  p_d <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0, y = 0,
      label = paste0("Dose UMAP not included\nMissing tumor group(s): ", paste(missing_doses, collapse = ", ")),
      size = 4,
      color = "grey30"
    ) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = "UMAP by treatment dose") +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}
p_d <- add_panel_tag(p_d, "C")

composition_plot <- composition_table
composition_plot$mouse <- factor(composition_plot$mouse, levels = mouse_levels)
composition_plot$initial_ploidy <- factor(composition_plot$initial_ploidy, levels = ploidy_levels)
composition_plot$dose <- factor(composition_plot$dose, levels = dose_levels)
composition_plot$cluster_final <- factor(composition_plot$cluster_final, levels = cluster_levels)

p_e <- ggplot2::ggplot(
  composition_plot,
  ggplot2::aes(mouse, proportion, fill = cluster_final)
) +
  ggplot2::geom_col(width = 0.84, color = "white", linewidth = 0.10) +
  ggplot2::facet_grid(
    . ~ initial_ploidy + dose,
    scales = "free_x",
    space = "free_x",
    drop = FALSE
  ) +
  ggplot2::scale_fill_manual(values = cluster_colors, drop = FALSE, name = "Cluster") +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = percent_labels,
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    title = "Cluster composition of each mouse",
    subtitle = "Each bar is one mouse and sums to 100% of its tumor cells",
    x = "Mouse",
    y = "Cluster composition"
  ) +
  figure_theme(9.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 50, hjust = 1, vjust = 1, size = 7.5),
    panel.spacing.x = grid::unit(0.15, "lines"),
    legend.position = "right"
  )
p_e <- add_panel_tag(p_e, "E")

group_plot <- group_summary
group_plot$initial_ploidy <- factor(group_plot$initial_ploidy, levels = ploidy_levels)
group_plot$dose <- factor(group_plot$dose, levels = dose_levels)
group_plot$cluster_final <- factor(group_plot$cluster_final, levels = cluster_levels)
group_n <- unique(group_plot[, c("initial_ploidy", "dose", "n_mice"), drop = FALSE])

p_f <- ggplot2::ggplot(
  group_plot,
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
  ggplot2::scale_fill_manual(values = cluster_colors, drop = FALSE, name = "Cluster") +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = percent_labels,
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    title = "Mouse-weighted composition by ploidy and dose",
    subtitle = "Cluster proportions are calculated per mouse, then averaged with equal mouse weights",
    x = "Gemcitabine dose",
    y = "Mean cluster composition"
  ) +
  figure_theme(9.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    legend.position = "right"
  )
p_f <- add_panel_tag(p_f, "F")

message("Writing Supplementary Figure 4 panels.")
save_pdf_png(p_a, file.path(figure_dir, "panel_SuppFig4A_umap_cluster"), 8.2, 6.5)
save_pdf_png(p_c, file.path(figure_dir, "panel_SuppFig4B_umap_initial_ploidy"), 7.2, 6.2)
save_pdf_png(p_d, file.path(figure_dir, "panel_SuppFig4C_umap_treatment_dose"), 7.2, 6.2)
save_pdf_png(p_b, file.path(figure_dir, "panel_SuppFig4D_umap_mouse_facets"), 13.5, 10.5)
save_pdf_png(p_e, file.path(figure_dir, "panel_SuppFig4E_cluster_composition_by_mouse"), 12.8, 6.5)
save_pdf_png(p_f, file.path(figure_dir, "panel_SuppFig4F_cluster_composition_by_ploidy_dose"), 9.5, 6.5)

left_column <- patchwork::wrap_plots(p_a, p_c, p_d, ncol = 1, heights = c(1, 1, 1))
upper_block <- patchwork::wrap_plots(left_column, p_b, ncol = 2, widths = c(1, 2.15))
bottom_row <- patchwork::wrap_plots(p_e, p_f, ncol = 2, widths = c(1.35, 1))
composite <- patchwork::wrap_plots(
  upper_block,
  bottom_row,
  ncol = 1,
  heights = c(1.75, 1)
) + patchwork::plot_annotation(
  title = "Supplementary Figure 4 | In-vivo tumor-cell landscape and cluster composition"
)
save_pdf_png(composite, file.path(figure_dir, "panel_SuppFig4_composite"), 17, 16, dpi = 300)

panel_contract <- data.frame(
  panel_id = c(
    "SuppFig4A", "SuppFig4A_png", "SuppFig4B", "SuppFig4B_png",
    "SuppFig4C", "SuppFig4C_png", "SuppFig4D", "SuppFig4D_png",
    "SuppFig4E", "SuppFig4E_png", "SuppFig4F", "SuppFig4F_png",
    "SuppFig4_composite", "SuppFig4_composite_png"
  ),
  filename = c(
    "panel_SuppFig4A_umap_cluster.pdf", "panel_SuppFig4A_umap_cluster.png",
    "panel_SuppFig4B_umap_initial_ploidy.pdf", "panel_SuppFig4B_umap_initial_ploidy.png",
    "panel_SuppFig4C_umap_treatment_dose.pdf", "panel_SuppFig4C_umap_treatment_dose.png",
    "panel_SuppFig4D_umap_mouse_facets.pdf", "panel_SuppFig4D_umap_mouse_facets.png",
    "panel_SuppFig4E_cluster_composition_by_mouse.pdf", "panel_SuppFig4E_cluster_composition_by_mouse.png",
    "panel_SuppFig4F_cluster_composition_by_ploidy_dose.pdf", "panel_SuppFig4F_cluster_composition_by_ploidy_dose.png",
    "panel_SuppFig4_composite.pdf", "panel_SuppFig4_composite.png"
  ),
  variant = rep(c("pdf", "png"), 7L),
  stringsAsFactors = FALSE
)
observed_figure_files <- sort(list.files(figure_dir, pattern = "[.](pdf|png)$"))
expected_figure_files <- sort(panel_contract$filename)
if (!identical(observed_figure_files, expected_figure_files)) {
  stop(
    "Supplementary Figure 4 exact figure inventory failed; missing: ",
    paste(setdiff(expected_figure_files, observed_figure_files), collapse = ", "),
    "; unexpected: ",
    paste(setdiff(observed_figure_files, expected_figure_files), collapse = ", "),
    call. = FALSE
  )
}
write_tsv(panel_contract, file.path(metadata_dir, "panel_contract.tsv"))

write_csv(color_key, file.path(table_dir, "cluster_color_key.csv"))
write_csv(composition_table, file.path(table_dir, "cluster_composition_by_mouse.csv"))
group_summary_output <- group_summary
group_summary_output$initial_ploidy <- as.character(group_summary_output$initial_ploidy)
group_summary_output$dose <- as.character(group_summary_output$dose)
group_summary_output$cluster_final <- as.character(group_summary_output$cluster_final)
write_csv(group_summary_output, file.path(table_dir, "cluster_composition_by_ploidy_dose.csv"))

input_manifest <- data.frame(
  role = c("seurat_metadata", "scvelo_cell_metrics"),
  path = c(seurat_metadata_path, scvelo_metrics_path),
  sha256 = c(file_sha256(seurat_metadata_path), file_sha256(scvelo_metrics_path)),
  bytes = c(file.info(seurat_metadata_path)$size, file.info(scvelo_metrics_path)$size),
  rows = c(nrow(seurat), nrow(metrics)),
  stringsAsFactors = FALSE
)
write_csv(input_manifest, file.path(metadata_dir, "input_manifest.csv"))

run_config <- data.frame(
  key = c(
    "module", "figure", "panel_count", "figure_file_count", "plot_shuffle_seed",
    "seurat_metadata", "seurat_metadata_sha256", "scvelo_cell_metrics",
    "scvelo_cell_metrics_sha256", "main_cell_universe", "cluster_field", "sample_field"
  ),
  value = c(
    "si_figure4", "Supplementary", "7", as.character(nrow(panel_contract)), "5826",
    seurat_metadata_path, input_manifest$sha256[input_manifest$role == "seurat_metadata"],
    scvelo_metrics_path, input_manifest$sha256[input_manifest$role == "scvelo_cell_metrics"],
    "all_tumor_cells_in_seurat_metadata", canonical_cluster_col, sample_col
  ),
  stringsAsFactors = FALSE
)
write_tsv(run_config, file.path(metadata_dir, "run_config.tsv"))

qc <- data.frame(
  key = c(
    "generated_at", "seurat_cells", "scvelo_cells", "seurat_cells_in_scvelo",
    "scvelo_cells_missing_from_seurat", "tumor_cells", "tumor_cells_in_scvelo",
    "tumor_mice", "canonical_clusters", "initial_ploidy_levels", "dose_levels",
    "dose_panel_informative", "main_cell_universe", "cluster_field", "sample_field"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    nrow(seurat),
    nrow(metrics),
    sum(seurat$in_scvelo),
    length(missing_metrics_cells),
    nrow(tumor),
    sum(tumor$in_scvelo),
    n_mice,
    n_clusters,
    paste(ploidy_levels[ploidy_levels %in% unique(as.character(tumor$initial_ploidy))], collapse = ","),
    paste(observed_dose_levels, collapse = ","),
    tolower(as.character(dose_panel_informative)),
    "all_tumor_cells_in_seurat_metadata",
    canonical_cluster_col,
    sample_col
  ),
  stringsAsFactors = FALSE
)
write_csv(qc, file.path(metadata_dir, "input_qc.csv"))

writeLines(capture.output(sessionInfo()), file.path(metadata_dir, "sessionInfo.txt"))
writeLines(
  c(
    "Supplementary Figure 4 generation completed.",
    paste0("Seurat metadata: ", seurat_metadata_path),
    paste0("scVelo metrics: ", scvelo_metrics_path),
    paste0("Main cell universe: all ", format(n_tumor, big.mark = ","), " tumor cells in Seurat metadata"),
    paste0("Tumor mice: ", n_mice),
    paste0("Clusters: ", paste(cluster_levels, collapse = ", ")),
    paste0("Cluster source field: ", canonical_cluster_col),
    paste0("Dose panel informative: ", dose_panel_informative),
    "Ploidy/dose group composition is the equal-weight mean of mouse-level proportions.",
    paste0("Output directory: ", output_dir)
  ),
  file.path(metadata_dir, "run_summary.txt")
)

message("Supplementary Figure 4 completed: ", output_dir)
