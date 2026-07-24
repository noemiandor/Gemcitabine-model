#!/usr/bin/env Rscript

# Supplementary Figures 4-7
#
# Descriptive UMAP/composition views plus standalone Hallmark cluster analysis.
# This script is intentionally downstream of the Figure 7 input workflow:
#   * seurat_metadata.csv supplies the complete Seurat cell universe and UMAP.
#   * scvelo_cell_metrics.csv supplies/validates sample-level ploidy and context.
#   * all_ploidy.tsv supplies endpoint tumor ploidy.
#   * the raw Seurat RDS supplies cluster-vs-rest DEGs for SI Figure 7.
#   * seurat_metadata_provenance.tsv is validated when available.
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
      "  Rscript Code/in-vivo/SI_figures/generate_supplementary_figures.R \\",
      "    --seurat-metadata Data/in-vivo/seurat_metadata.csv \\",
      "    --scvelo-metrics Data/in-vivo/scvelo_cell_metrics.csv \\",
      "    --all-ploidy Data/in-vivo/all_ploidy.tsv \\",
      "    --seurat-rds /absolute/path/to/integrated_sct_cca_seurat_final_reclustered.rds \\",
      "    [--seurat-metadata-provenance Data/in-vivo/seurat_metadata_provenance.tsv] \\",
      "    --config Code/in-vivo/figure7/figure7_config.yaml \\",
      "    --output-dir Results/in-vivo/SI_figures",
      "",
      "Options:",
      "  --overwrite       Replace this script's known outputs in a nonempty output directory.",
      "  --deg-cache-dir    Reuse a complete, validated set of per-cluster SI7 DEG CSV files.",
      "  --table-cache-dir  Canonical 32-table cache (default: Data/in-vivo/SIfigures).",
      "  --force-reanalysis Ignore a valid table cache and rerun all analyses.",
      "  --skip-si7        Generate SI Figures 4-6 only (validation/debug only).",
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

read_tsv_checked <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
  if (file.info(path)$size <= 0) stop(label, " is empty: ", path, call. = FALSE)
  out <- utils::read.delim(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    quote = "", comment.char = "", na.strings = c("", "NA", "NaN")
  )
  if (!nrow(out)) stop(label, " has no data rows: ", path, call. = FALSE)
  out
}

key_value_map <- function(data, label) {
  if (!identical(names(data), c("key", "value"))) {
    stop(label, " must have exactly key and value columns", call. = FALSE)
  }
  if (anyNA(data$key) || any(!nzchar(data$key)) || anyDuplicated(data$key)) {
    stop(label, " contains missing or duplicated keys", call. = FALSE)
  }
  stats::setNames(as.character(data$value), data$key)
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
require_package("yaml")

script_path <- resolve_script_path()
if (is.na(script_path)) stop("Cannot resolve the script path", call. = FALSE)
script_dir <- dirname(script_path)
repo_root <- normalizePath(dirname(dirname(dirname(script_dir))), mustWork = TRUE)
table_cache_dir <- resolve_path(
  arg_value(args, "table-cache-dir", file.path("Data", "in-vivo", "SIfigures")),
  repo_root,
  must_work = FALSE
)
force_reanalysis <- arg_flag(args, "force-reanalysis", FALSE)
cache_validator <- file.path(repo_root, "Code", "tools", "validate_si_figures_table_cache.py")
if (!file.exists(cache_validator)) {
  stop("Missing SI Figures table-cache validator: ", cache_validator, call. = FALSE)
}
use_table_cache <- !force_reanalysis && dir.exists(table_cache_dir)
if (use_table_cache) {
  validation_status <- system2(
    "python3",
    c(shQuote(cache_validator), "--cache-dir", shQuote(table_cache_dir))
  )
  if (!identical(validation_status, 0L)) {
    stop(
      "Data/in-vivo/SIfigures exists but is not a complete valid cache; ",
      "fix it or pass --force-reanalysis",
      call. = FALSE
    )
  }
}
config_path <- resolve_path(
  arg_value(args, "config", file.path("Code", "in-vivo", "figure7", "figure7_config.yaml")),
  repo_root,
  must_work = TRUE
)
config <- yaml::read_yaml(config_path)
si_config <- config$si_figures
if (is.null(si_config)) stop("Figure 7 config is missing si_figures", call. = FALSE)

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
seurat_metadata_provenance_path <- resolve_path(
  arg_value(
    args,
    "seurat-metadata-provenance",
    file.path("Data", "in-vivo", "seurat_metadata_provenance.tsv")
  ),
  repo_root,
  must_work = FALSE
)
all_ploidy_path <- resolve_path(
  arg_value(args, "all-ploidy", file.path("Data", "in-vivo", "all_ploidy.tsv")),
  repo_root,
  must_work = FALSE
)
skip_si7 <- arg_flag(args, "skip-si7", FALSE)
deg_cache_arg <- arg_value(args, "deg-cache-dir", NULL)
deg_cache_dir <- if (is.null(deg_cache_arg)) {
  NA_character_
} else {
  resolve_path(deg_cache_arg, repo_root, must_work = FALSE)
}
if (!is.na(deg_cache_dir) && !dir.exists(deg_cache_dir)) {
  stop("--deg-cache-dir does not exist: ", deg_cache_dir, call. = FALSE)
}
if (force_reanalysis && !is.na(deg_cache_dir)) {
  stop("--force-reanalysis and --deg-cache-dir are mutually exclusive", call. = FALSE)
}

output_dir <- resolve_path(
  arg_value(args, "output-dir", file.path("Results", "in-vivo", "SI_figures")),
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

if (use_table_cache) {
  message("Using canonical SI Figures table cache; analysis is disabled: ", table_cache_dir)
  cache_files <- sort(list.files(
    table_cache_dir,
    pattern = "[.](csv|tsv)$",
    full.names = TRUE
  ))
  cache_files <- cache_files[basename(cache_files) != "manifest.tsv"]
  if (length(cache_files) != 32L) {
    stop("Validated SI Figures table cache must contain exactly 32 data files", call. = FALSE)
  }
  copied <- file.copy(
    cache_files,
    file.path(table_dir, basename(cache_files)),
    overwrite = TRUE,
    copy.mode = TRUE
  )
  if (!all(copied)) stop("Failed to copy the canonical SI Figures table cache", call. = FALSE)
  copied_paths <- file.path(table_dir, basename(cache_files))
  copied_hashes <- vapply(copied_paths, file_sha256, character(1))
  source_hashes <- vapply(cache_files, file_sha256, character(1))
  if (!identical(unname(copied_hashes), unname(source_hashes))) {
    stop("SI Figures table cache changed while being copied into the run", call. = FALSE)
  }
}

if (!use_table_cache) {
message("Reading Seurat metadata: ", seurat_metadata_path)
seurat <- read_csv_checked(seurat_metadata_path, "Seurat metadata")
message("Reading scVelo cell metrics: ", scvelo_metrics_path)
metrics <- read_csv_checked(scvelo_metrics_path, "scVelo cell metrics")
message("Reading endpoint ploidy: ", all_ploidy_path)
endpoint_ploidy <- read_tsv_checked(all_ploidy_path, "endpoint ploidy")
input_provenance <- character(0)
if (file.exists(seurat_metadata_provenance_path)) {
  message("Reading Seurat metadata provenance: ", seurat_metadata_provenance_path)
  input_provenance <- key_value_map(
    read_tsv_checked(seurat_metadata_provenance_path, "Seurat metadata provenance"),
    "Seurat metadata provenance"
  )
}
required_provenance <- c(
  "source_seurat_rds", "source_seurat_rds_sha256", "source_object_cells",
  "seurat_metadata_sha256", "scvelo_metrics_sha256", "umap_reduction",
  "cluster_id_field", "base_cluster_field", "clustering_resolution",
  "cluster_annotation_field", "sample_field", "dose_field", "ploidy_field",
  "context_field", "cellcycle_mapping", "inclusion_context",
  "inclusion_ploidy_levels", "inclusion_dose_levels",
  "inclusion_required_nonmissing", "source_qc_fields", "source_qc_policy",
  "figure7_config_sha256", "export_script_sha256", "source_code_revision"
)
missing_provenance <- if (length(input_provenance)) {
  setdiff(required_provenance, names(input_provenance))
} else {
  character(0)
}
if (length(missing_provenance)) {
  stop(
    "Seurat metadata provenance is missing required key(s): ",
    paste(missing_provenance, collapse = ", "),
    call. = FALSE
  )
}
if (length(input_provenance) &&
    (!identical(input_provenance[["seurat_metadata_sha256"]], file_sha256(seurat_metadata_path)) ||
     !identical(input_provenance[["scvelo_metrics_sha256"]], file_sha256(scvelo_metrics_path)) ||
     !identical(input_provenance[["figure7_config_sha256"]], file_sha256(config_path)))) {
  stop("Seurat metadata provenance checksum does not match the selected inputs/config", call. = FALSE)
}
if (length(input_provenance) &&
    !grepl("^[0-9a-f]{64}$", input_provenance[["source_seurat_rds_sha256"]])) {
  stop("Seurat metadata provenance has an invalid source Seurat RDS SHA-256", call. = FALSE)
}

seurat_rds_default <- if (length(input_provenance)) input_provenance[["source_seurat_rds"]] else NULL
seurat_rds_arg <- arg_value(args, "seurat-rds", seurat_rds_default)
seurat_rds_path <- if (is.null(seurat_rds_arg)) NA_character_ else resolve_path(
  seurat_rds_arg, repo_root, must_work = FALSE
)
if (!skip_si7 && (is.na(seurat_rds_path) || !file.exists(seurat_rds_path))) {
  stop("SI Figure 7 requires --seurat-rds pointing to the raw Seurat RDS", call. = FALSE)
}

for (column in c("cell", "UMAP_1", "UMAP_2", "Dose", "S.Score", "harvest", "barcode_raw")) {
  assert_single_column(seurat, column, "Seurat metadata")
}
assert_single_column(metrics, "cell", "scVelo cell metrics")
for (column in c("file", "cell_id", "ploidy")) {
  assert_single_column(endpoint_ploidy, column, "endpoint ploidy")
}

canonical_cluster_col <- resolve_column(
  seurat,
  as.character(si_config$cluster_id_field),
  "Seurat canonical cluster"
)
annotation_col <- resolve_column(
  seurat,
  as.character(si_config$cluster_annotation_field),
  "Seurat canonical cluster annotation"
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
  as.character(si_config$sample_field),
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
seurat$s_phase_score <- numeric_strict(seurat$S.Score, "Seurat S.Score")
if (any(!is.finite(seurat$UMAP_1)) || any(!is.finite(seurat$UMAP_2))) {
  stop("Seurat metadata contains missing or nonfinite UMAP coordinates", call. = FALSE)
}

seurat$mouse <- clean_character(seurat[[sample_col]])
seurat$cluster <- clean_character(seurat[[canonical_cluster_col]])
seurat$cluster_annotation <- clean_character(seurat[[annotation_col]])
seurat$dose <- standardize_dose(seurat$Dose, "Seurat Dose")
if (anyNA(seurat$mouse)) stop("Seurat metadata contains missing mouse/sample identities", call. = FALSE)
if (anyNA(seurat$cluster)) stop("Seurat metadata contains missing canonical cluster labels", call. = FALSE)
if (anyNA(seurat$cluster_annotation)) {
  stop("Seurat metadata contains missing canonical cluster annotations", call. = FALSE)
}
cellcycle_mapping <- unlist(si_config$cellcycle_mapping, use.names = TRUE)
unknown_annotations <- setdiff(unique(seurat$cluster_annotation), names(cellcycle_mapping))
if (length(unknown_annotations)) {
  stop(
    "Cluster annotations are absent from the reviewed CellCycle mapping: ",
    paste(unknown_annotations, collapse = ", "),
    call. = FALSE
  )
}
seurat$cellcycle_classification <- unname(cellcycle_mapping[seurat$cluster_annotation])

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

endpoint_ploidy$file <- clean_character(endpoint_ploidy$file)
endpoint_ploidy$cell_id <- clean_character(endpoint_ploidy$cell_id)
endpoint_ploidy$ploidy <- numeric_strict(endpoint_ploidy$ploidy, "endpoint ploidy value")
endpoint_ploidy$join_key <- paste(endpoint_ploidy$file, endpoint_ploidy$cell_id, sep = "\r")
if (anyNA(endpoint_ploidy$file) || anyNA(endpoint_ploidy$cell_id) ||
    any(!is.finite(endpoint_ploidy$ploidy)) || anyDuplicated(endpoint_ploidy$join_key)) {
  stop("Endpoint ploidy requires unique, nonmissing file + cell_id keys and finite ploidy", call. = FALSE)
}
seurat$endpoint_file <- paste0(clean_character(seurat$harvest), ".sps.cbs")
seurat$endpoint_cell_id <- clean_character(seurat$barcode_raw)
seurat$endpoint_join_key <- paste(seurat$endpoint_file, seurat$endpoint_cell_id, sep = "\r")
endpoint_index <- match(seurat$endpoint_join_key, endpoint_ploidy$join_key)
seurat$endpoint_ploidy <- endpoint_ploidy$ploidy[endpoint_index]
tumor_endpoint_missing <- seurat$context == "Tumor" & !is.finite(seurat$endpoint_ploidy)
cellline_endpoint_present <- seurat$context == "CellLine" & is.finite(seurat$endpoint_ploidy)
if (any(tumor_endpoint_missing)) {
  stop(
    "Endpoint ploidy is missing for ", sum(tumor_endpoint_missing),
    " tumor cells; first cells: ",
    paste(head(seurat$cell[tumor_endpoint_missing], 10L), collapse = ", "),
    call. = FALSE
  )
}
if (any(cellline_endpoint_present)) {
  stop("Endpoint ploidy unexpectedly matched CellLine cells", call. = FALSE)
}

ploidy_levels <- as.character(unlist(si_config$inclusion$ploidy_levels))
dose_levels <- as.character(unlist(si_config$inclusion$dose_levels))
cluster_levels <- sort_cluster_levels(seurat$cluster)
if (!length(cluster_levels)) stop("No canonical clusters remain", call. = FALSE)
cluster_order <- stats::setNames(seq_along(cluster_levels), cluster_levels)
cluster_colors <- setNames(
  grDevices::hcl.colors(max(length(cluster_levels), 3L), palette = "Dark 3")[seq_along(cluster_levels)],
  cluster_levels
)
cluster_annotation_map <- mapping_by_group(
  seurat$cluster, seurat$cluster_annotation, "cluster annotation"
)[cluster_levels]
cluster_cellcycle_map <- mapping_by_group(
  seurat$cluster, seurat$cellcycle_classification, "CellCycle classification"
)[cluster_levels]
if (anyNA(cluster_annotation_map) || anyNA(cluster_cellcycle_map)) {
  stop("Every cluster must have one nonmissing annotation and CellCycle classification", call. = FALSE)
}

seurat$cluster_order <- unname(cluster_order[seurat$cluster])
seurat$cluster_color <- unname(cluster_colors[seurat$cluster])
seurat$dose_mg_per_kg <- as.numeric(sub("mg/kg$", "", seurat$dose))
seurat$treatment <- ifelse(
  is.na(seurat$dose), NA_character_,
  ifelse(seurat$dose == "0mg/kg", "Control", "Gemcitabine")
)
required_nonmissing <- as.character(unlist(si_config$inclusion$required_nonmissing))
canonical_required <- list(
  cell_id = seurat$cell,
  UMAP_1 = seurat$UMAP_1,
  UMAP_2 = seurat$UMAP_2,
  sample_id = seurat$mouse,
  cluster_id = seurat$cluster,
  cluster_annotation = seurat$cluster_annotation,
  initial_ploidy = seurat$initial_ploidy,
  dose = seurat$dose,
  cellcycle_classification = seurat$cellcycle_classification
)
unknown_required <- setdiff(required_nonmissing, names(canonical_required))
if (length(unknown_required)) {
  stop("Unknown configured SI Figure 5 required field(s): ", paste(unknown_required, collapse = ", "), call. = FALSE)
}
nonmissing_required <- rep(TRUE, nrow(seurat))
for (field in required_nonmissing) {
  values <- canonical_required[[field]]
  nonmissing_required <- nonmissing_required & !is.na(values)
  if (is.character(values)) nonmissing_required <- nonmissing_required & nzchar(values)
}
included_context <- as.character(si_config$inclusion$context)
seurat$included_in_si_figures <- (
  seurat$context == included_context &
    seurat$initial_ploidy %in% ploidy_levels &
    seurat$dose %in% dose_levels &
    nonmissing_required
)
seurat$exclusion_reason <- ifelse(
  seurat$included_in_si_figures,
  "",
  ifelse(
    is.na(seurat$context) | seurat$context != included_context,
    "context_not_Tumor",
    ifelse(
      is.na(seurat$initial_ploidy) | !seurat$initial_ploidy %in% ploidy_levels,
      "invalid_initial_ploidy",
      ifelse(
        is.na(seurat$dose) | !seurat$dose %in% dose_levels,
        "invalid_dose",
        "missing_required_metadata"
      )
    )
  )
)

canonical_cells <- data.frame(
  cell_id = seurat$cell,
  UMAP_1 = seurat$UMAP_1,
  UMAP_2 = seurat$UMAP_2,
  sample_id = seurat$mouse,
  cluster_id = seurat$cluster,
  cluster_annotation = seurat$cluster_annotation,
  cluster_order = seurat$cluster_order,
  cluster_color = seurat$cluster_color,
  initial_ploidy = seurat$initial_ploidy,
  treatment = seurat$treatment,
  dose = seurat$dose,
  dose_mg_per_kg = seurat$dose_mg_per_kg,
  s_phase_score = seurat$s_phase_score,
  endpoint_ploidy = seurat$endpoint_ploidy,
  endpoint_file = seurat$endpoint_file,
  endpoint_cell_id = seurat$endpoint_cell_id,
  cellcycle_classification = seurat$cellcycle_classification,
  context = seurat$context,
  included_in_si_figures = seurat$included_in_si_figures,
  exclusion_reason = seurat$exclusion_reason,
  stringsAsFactors = FALSE
)
canonical_cells <- canonical_cells[
  order(canonical_cells$cell_id, method = "radix"),
  ,
  drop = FALSE
]
if (anyDuplicated(canonical_cells$cell_id)) stop("Canonical cell table contains duplicated cell IDs", call. = FALSE)

tumor <- seurat[seurat$included_in_si_figures, , drop = FALSE]
if (!nrow(tumor)) stop("No tumor cells remain after context filtering", call. = FALSE)
if (anyNA(tumor$initial_ploidy)) stop("Tumor cells contain unresolved initial ploidy", call. = FALSE)
if (any(!tumor$initial_ploidy %in% c("2N", "4N"))) stop("Unexpected tumor ploidy value", call. = FALSE)
if (anyNA(tumor$dose)) stop("Tumor cells contain missing Dose values", call. = FALSE)

mouse_dose <- mapping_by_group(tumor$mouse, tumor$dose, "dose")
mouse_ploidy <- mapping_by_group(tumor$mouse, tumor$initial_ploidy, "initial ploidy")
mouse_context <- mapping_by_group(tumor$mouse, tumor$context, "context")
if (any(mouse_context != "Tumor")) stop("Tumor subset contains a non-tumor mouse", call. = FALSE)

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

ploidy_colors <- c("2N" = "#4C78A8", "4N" = "#E45756")
dose_colors <- c("0mg/kg" = "#666666", "30mg/kg" = "#d95f02", "120mg/kg" = "#1b9e77")

cluster_key <- data.frame(
  cluster_id = cluster_levels,
  cluster_annotation = unname(cluster_annotation_map[cluster_levels]),
  cellcycle_classification = unname(cluster_cellcycle_map[cluster_levels]),
  cluster_order = seq_along(cluster_levels),
  color = unname(cluster_colors[cluster_levels]),
  n_all_cells = as.integer(table(factor(seurat$cluster, levels = cluster_levels))),
  n_included_tumor_cells = as.integer(table(factor(tumor$cluster, levels = cluster_levels))),
  stringsAsFactors = FALSE
)
if (anyNA(cluster_key) || any(!nzchar(cluster_key$cluster_annotation)) ||
    anyDuplicated(cluster_key$cluster_order) || anyDuplicated(cluster_key$color) ||
    any(!grepl("^#[0-9A-Fa-f]{6}$", cluster_key$color))) {
  stop("Cluster key violates annotation, order, or color requirements", call. = FALSE)
}

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
composition_table$cluster_annotation <- unname(cluster_annotation_map[composition_table$cluster_final])
composition_table$cluster_order <- unname(cluster_order[composition_table$cluster_final])
composition_table$denominator_definition <- "all canonical cells with included_in_si_figures=true within mouse"
composition_table <- composition_table[, c(
  "mouse", "initial_ploidy", "dose", "cluster_final", "cluster_annotation",
  "cluster_order", "n_cells", "total_cells", "proportion", "denominator_definition"
)]

mouse_sums <- tapply(composition_table$proportion, composition_table$mouse, sum)
if (any(abs(mouse_sums - 1) > 1e-12)) {
  stop("Cluster proportions do not sum to one within every mouse", call. = FALSE)
}
canonical_included_counts <- table(canonical_cells$sample_id[canonical_cells$included_in_si_figures])
composition_denominators <- stats::setNames(
  composition_table$total_cells[!duplicated(composition_table$mouse)],
  composition_table$mouse[!duplicated(composition_table$mouse)]
)
if (!identical(
  as.integer(canonical_included_counts[names(composition_denominators)]),
  as.integer(composition_denominators)
)) {
  stop("Composition denominators do not reconcile to the canonical included-cell table", call. = FALSE)
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
    cluster_annotation = piece$cluster_annotation[[1L]],
    cluster_order = piece$cluster_order[[1L]],
    n_mice = length(unique(piece$mouse)),
    sum_n_cells = sum(piece$n_cells),
    sum_total_cells = sum(piece$total_cells),
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
  group_summary$cluster_order
), , drop = FALSE]

group_sums <- aggregate(mean_proportion ~ initial_ploidy + dose, data = group_summary, sum)
if (any(abs(group_sums$mean_proportion - 1) > 1e-12)) {
  stop("Equal-mouse group mean proportions do not sum to one", call. = FALSE)
}
} else {
  canonical_cells <- read_csv_checked(
    file.path(table_dir, "si_figures_cell_metadata.csv"),
    "cached SI Figures cell metadata"
  )
  cluster_key <- read_tsv_checked(
    file.path(table_dir, "si_figures_cluster_key.tsv"),
    "cached SI Figures cluster key"
  )
  required_cached_cell_columns <- c(
    "cell_id", "UMAP_1", "UMAP_2", "sample_id", "cluster_id",
    "cluster_annotation", "cluster_order", "cluster_color", "initial_ploidy",
    "treatment", "dose", "dose_mg_per_kg", "s_phase_score",
    "endpoint_ploidy", "endpoint_file", "endpoint_cell_id",
    "cellcycle_classification", "context", "included_in_si_figures",
    "exclusion_reason"
  )
  missing_cached_columns <- setdiff(required_cached_cell_columns, names(canonical_cells))
  if (length(missing_cached_columns)) {
    stop(
      "Cached SI Figures cell metadata is missing: ",
      paste(missing_cached_columns, collapse = ", "),
      call. = FALSE
    )
  }
  canonical_cells$UMAP_1 <- numeric_strict(canonical_cells$UMAP_1, "cached UMAP_1")
  canonical_cells$UMAP_2 <- numeric_strict(canonical_cells$UMAP_2, "cached UMAP_2")
  canonical_cells$s_phase_score <- numeric_strict(
    canonical_cells$s_phase_score,
    "cached S phase score"
  )
  canonical_cells$endpoint_ploidy <- suppressWarnings(as.numeric(canonical_cells$endpoint_ploidy))
  canonical_cells$included_in_si_figures <- tolower(
    as.character(canonical_cells$included_in_si_figures)
  ) %in% c("true", "t", "1")
  if (anyDuplicated(canonical_cells$cell_id)) {
    stop("Cached SI Figures cell metadata contains duplicated cell IDs", call. = FALSE)
  }

  cluster_key$cluster_order <- as.integer(cluster_key$cluster_order)
  cluster_key <- cluster_key[order(cluster_key$cluster_order), , drop = FALSE]
  cluster_levels <- as.character(cluster_key$cluster_id)
  cluster_colors <- stats::setNames(as.character(cluster_key$color), cluster_levels)
  cluster_order <- stats::setNames(cluster_key$cluster_order, cluster_levels)
  cluster_annotation_map <- stats::setNames(
    as.character(cluster_key$cluster_annotation),
    cluster_levels
  )
  cluster_cellcycle_map <- stats::setNames(
    as.character(cluster_key$cellcycle_classification),
    cluster_levels
  )
  ploidy_levels <- as.character(unlist(si_config$inclusion$ploidy_levels))
  dose_levels <- as.character(unlist(si_config$inclusion$dose_levels))
  canonical_cluster_col <- as.character(si_config$cluster_id_field)
  sample_col <- as.character(si_config$sample_field)
  input_provenance <- character(0)

  seurat <- data.frame(
    cell = clean_character(canonical_cells$cell_id),
    UMAP_1 = canonical_cells$UMAP_1,
    UMAP_2 = canonical_cells$UMAP_2,
    mouse = clean_character(canonical_cells$sample_id),
    cluster = clean_character(canonical_cells$cluster_id),
    cluster_annotation = clean_character(canonical_cells$cluster_annotation),
    cluster_order = as.integer(canonical_cells$cluster_order),
    cluster_color = clean_character(canonical_cells$cluster_color),
    initial_ploidy = clean_character(canonical_cells$initial_ploidy),
    treatment = clean_character(canonical_cells$treatment),
    dose = clean_character(canonical_cells$dose),
    dose_mg_per_kg = suppressWarnings(as.numeric(canonical_cells$dose_mg_per_kg)),
    s_phase_score = canonical_cells$s_phase_score,
    endpoint_ploidy = canonical_cells$endpoint_ploidy,
    endpoint_file = clean_character(canonical_cells$endpoint_file),
    endpoint_cell_id = clean_character(canonical_cells$endpoint_cell_id),
    cellcycle_classification = clean_character(canonical_cells$cellcycle_classification),
    context = clean_character(canonical_cells$context),
    included_in_si_figures = canonical_cells$included_in_si_figures,
    exclusion_reason = as.character(canonical_cells$exclusion_reason),
    stringsAsFactors = FALSE
  )
  seurat$in_scvelo <- TRUE
  missing_metrics_cells <- character(0)
  tumor_endpoint_missing <- seurat$context == "Tumor" & !is.finite(seurat$endpoint_ploidy)
  cellline_endpoint_present <- seurat$context == "CellLine" & is.finite(seurat$endpoint_ploidy)
  if (any(tumor_endpoint_missing) || any(cellline_endpoint_present)) {
    stop("Cached endpoint-ploidy matching violates the SI Figure 6 contract", call. = FALSE)
  }
  metrics <- seurat
  endpoint_ploidy <- read_csv_checked(
    file.path(table_dir, "si_figure6_endpoint_ploidy_join_audit.csv"),
    "cached endpoint-ploidy audit"
  )

  tumor <- seurat[seurat$included_in_si_figures, , drop = FALSE]
  sample_metadata <- unique(tumor[, c("mouse", "initial_ploidy", "dose"), drop = FALSE])
  sample_metadata$initial_ploidy <- factor(
    sample_metadata$initial_ploidy,
    levels = ploidy_levels
  )
  sample_metadata$dose <- factor(sample_metadata$dose, levels = dose_levels)
  sample_metadata <- sample_metadata[
    order(sample_metadata$initial_ploidy, sample_metadata$dose, sample_metadata$mouse),
    ,
    drop = FALSE
  ]
  mouse_levels <- as.character(sample_metadata$mouse)
  mouse_ploidy <- stats::setNames(
    as.character(sample_metadata$initial_ploidy),
    mouse_levels
  )
  mouse_dose <- stats::setNames(as.character(sample_metadata$dose), mouse_levels)
  tumor$mouse <- factor(tumor$mouse, levels = mouse_levels)
  tumor$initial_ploidy <- factor(tumor$initial_ploidy, levels = ploidy_levels)
  tumor$dose <- factor(tumor$dose, levels = dose_levels)
  tumor$cluster <- factor(tumor$cluster, levels = cluster_levels)

  ploidy_colors <- c("2N" = "#4C78A8", "4N" = "#E45756")
  dose_colors <- c("0mg/kg" = "#666666", "30mg/kg" = "#d95f02", "120mg/kg" = "#1b9e77")

  s4_context <- read_csv_checked(
    file.path(table_dir, "si_figure4_cluster_context_composition.csv"),
    "cached SI Figure 4 context composition"
  )
  s4_ploidy <- read_csv_checked(
    file.path(table_dir, "si_figure4_cluster_initial_ploidy_composition.csv"),
    "cached SI Figure 4 ploidy composition"
  )
  composition_table <- read_csv_checked(
    file.path(table_dir, "si_figure5_cluster_composition_by_mouse.csv"),
    "cached SI Figure 5 mouse composition"
  )
  group_summary <- read_csv_checked(
    file.path(table_dir, "si_figure5_mouse_weighted_composition_by_initial_ploidy_dose.csv"),
    "cached SI Figure 5 mouse-weighted composition"
  )
  s5_dose_composition <- read_csv_checked(
    file.path(table_dir, "si_figure5_cluster_dose_composition.csv"),
    "cached SI Figure 5 dose composition"
  )
  s5_ploidy_composition <- read_csv_checked(
    file.path(table_dir, "si_figure5_cluster_initial_ploidy_composition.csv"),
    "cached SI Figure 5 ploidy composition"
  )
  endpoint_audit <- endpoint_ploidy

  factor_composition <- function(df, group_levels, fill_levels) {
    df$group_value <- factor(as.character(df$group_value), levels = group_levels)
    df$fill_value <- factor(as.character(df$fill_value), levels = fill_levels)
    df$n_cells <- as.numeric(df$n_cells)
    df$total_cells <- as.numeric(df$total_cells)
    df$proportion <- as.numeric(df$proportion)
    df
  }
  s4_context <- factor_composition(s4_context, cluster_levels, c("Tumor", "CellLine"))
  s4_ploidy <- factor_composition(s4_ploidy, cluster_levels, ploidy_levels)
  s5_dose_composition <- factor_composition(
    s5_dose_composition,
    cluster_levels,
    dose_levels
  )
  s5_ploidy_composition <- factor_composition(
    s5_ploidy_composition,
    cluster_levels,
    ploidy_levels
  )
  composition_table$n_cells <- as.numeric(composition_table$n_cells)
  composition_table$total_cells <- as.numeric(composition_table$total_cells)
  composition_table$proportion <- as.numeric(composition_table$proportion)
  group_summary$n_mice <- as.integer(group_summary$n_mice)
  group_summary$mean_proportion <- as.numeric(group_summary$mean_proportion)
}

n_tumor <- nrow(tumor)
n_mice <- length(mouse_levels)
n_clusters <- length(cluster_levels)
point_size <- if (n_tumor > 50000L) 0.08 else if (n_tumor > 20000L) 0.14 else 0.24
set.seed(5826)
plot_data <- tumor[sample.int(n_tumor), , drop = FALSE]

if (overwrite) {
  stale_panels <- list.files(
    figure_dir,
    pattern = "^panel_SuppFig[4567].*[.](pdf|png)$",
    full.names = TRUE
  )
  if (length(stale_panels)) unlink(stale_panels)
}

panel_rows <- list()
record_panel <- function(plot, panel_id, filename_stub, width, height, caption_role) {
  save_pdf_png(plot, file.path(figure_dir, filename_stub), width, height)
  panel_rows[[length(panel_rows) + 1L]] <<- data.frame(
    panel_id = c(panel_id, paste0(panel_id, "_png")),
    filename = c(paste0(filename_stub, ".pdf"), paste0(filename_stub, ".png")),
    variant = c("pdf", "png"),
    caption_role = caption_role,
    stringsAsFactors = FALSE
  )
  invisible(plot)
}

shuffle_cells <- function(df, seed) {
  set.seed(seed)
  df[sample.int(nrow(df)), , drop = FALSE]
}

make_umap_discrete <- function(df, field, colors, title, legend_title, tag, subtitle = NULL, labels = FALSE) {
  df <- shuffle_cells(df, 5826L + utf8ToInt(tag)[[1L]])
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(UMAP_1, UMAP_2, color = .data[[field]])
  ) +
    ggplot2::geom_point(size = point_size, alpha = 0.76, stroke = 0) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE, name = legend_title) +
    ggplot2::guides(color = umap_color_guide()) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2") +
    umap_theme(10)
  if (labels) {
    centers <- aggregate(
      cbind(UMAP_1, UMAP_2) ~ cluster,
      data = df,
      FUN = stats::median
    )
    p <- p + ggplot2::geom_label(
      data = centers,
      ggplot2::aes(x = UMAP_1, y = UMAP_2, label = cluster),
      inherit.aes = FALSE,
      size = 2.5,
      linewidth = 0.2,
      fill = "white",
      color = "#222222",
      alpha = 0.86,
      label.padding = grid::unit(0.10, "lines")
    )
  }
  add_panel_tag(p, tag)
}

make_umap_continuous <- function(
  df,
  field,
  title,
  legend_title,
  tag,
  limits = NULL,
  subtitle = NULL,
  diverging = FALSE
) {
  df <- shuffle_cells(df, 6826L + utf8ToInt(tag)[[1L]])
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(UMAP_1, UMAP_2, color = .data[[field]])
  ) +
    ggplot2::geom_point(size = point_size, alpha = 0.80, stroke = 0) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2") +
    umap_theme(10)
  if (diverging) {
    p <- p + ggplot2::scale_color_gradient2(
      low = "#2C7BB6", mid = "white", high = "#D7191C",
      midpoint = 0, limits = limits, name = legend_title
    )
  } else {
    p <- p + ggplot2::scale_color_gradient(
      low = "#2C7BB6", high = "#D7191C",
      limits = limits, name = legend_title
    )
  }
  add_panel_tag(p, tag)
}

composition_by <- function(df, group_field, fill_field, group_levels, fill_levels) {
  out <- as.data.frame(
    table(
      group_value = factor(df[[group_field]], levels = group_levels),
      fill_value = factor(df[[fill_field]], levels = fill_levels)
    ),
    stringsAsFactors = FALSE
  )
  names(out)[3L] <- "n_cells"
  out$group_value <- factor(as.character(out$group_value), levels = group_levels)
  out$fill_value <- factor(as.character(out$fill_value), levels = fill_levels)
  out$total_cells <- ave(out$n_cells, out$group_value, FUN = sum)
  out$proportion <- ifelse(out$total_cells > 0, out$n_cells / out$total_cells, 0)
  out$group_field <- group_field
  out$fill_field <- fill_field
  out$denominator_definition <- paste0("all plotted cells within ", group_field)
  out
}

make_composition_plot <- function(
  df,
  fill_colors,
  title,
  x_title,
  y_title,
  tag,
  use_proportion
) {
  y_field <- if (use_proportion) "proportion" else "n_cells"
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(group_value, .data[[y_field]], fill = fill_value)
  ) +
    ggplot2::geom_col(width = 0.82, color = "white", linewidth = 0.12) +
    ggplot2::scale_fill_manual(values = fill_colors, drop = FALSE) +
    ggplot2::labs(title = title, x = x_title, y = y_title, fill = NULL) +
    figure_theme(9.5) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  if (use_proportion) {
    p <- p + ggplot2::scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      labels = percent_labels,
      expand = ggplot2::expansion(mult = c(0, 0))
    ) + ggplot2::coord_cartesian(ylim = c(0, 1))
  } else {
    p <- p + ggplot2::scale_y_continuous(
      labels = function(x) format(x, big.mark = ",", scientific = FALSE),
      expand = ggplot2::expansion(mult = c(0, 0.04))
    )
  }
  add_panel_tag(p, tag)
}

all_cells <- seurat
all_cells$cluster <- factor(all_cells$cluster, levels = cluster_levels)
all_cells$context <- factor(all_cells$context, levels = c("Tumor", "CellLine"))
all_cells$initial_ploidy <- factor(all_cells$initial_ploidy, levels = ploidy_levels)
context_colors <- c("Tumor" = "#4C78A8", "CellLine" = "#F2CF5B")

message("Generating Supplementary Figure 4.")
s4a <- make_umap_discrete(
  all_cells, "cluster", cluster_colors, "UMAP by cluster", "Cluster", "A",
  sprintf("Tumor and CellLine cells; n = %s", format(nrow(all_cells), big.mark = ",")),
  labels = TRUE
)
s4b <- make_umap_discrete(
  all_cells, "initial_ploidy", ploidy_colors,
  "UMAP by initial ploidy", "Initial ploidy", "B", "Tumor and CellLine cells"
)
s4c <- make_umap_discrete(
  all_cells, "context", context_colors,
  "UMAP by Tumor/CellLine context", "Context", "C", "Tumor and CellLine cells"
)
s4d <- make_umap_continuous(
  all_cells, "s_phase_score", "UMAP by S phase score", "S phase score", "D",
  limits = range(all_cells$s_phase_score, na.rm = TRUE), diverging = TRUE
)
if (!use_table_cache) {
  s4_context <- composition_by(
    all_cells, "cluster", "context", cluster_levels, c("Tumor", "CellLine")
  )
  s4_ploidy <- composition_by(
    all_cells, "cluster", "initial_ploidy", cluster_levels, ploidy_levels
  )
}
if (any(abs(tapply(s4_context$proportion, s4_context$group_value, sum) - 1) > 1e-12) ||
    any(abs(tapply(s4_ploidy$proportion, s4_ploidy$group_value, sum) - 1) > 1e-12)) {
  stop("SI Figure 4 within-cluster proportions do not sum to one", call. = FALSE)
}
s4e <- make_composition_plot(
  s4_context, context_colors, "Tumor and CellLine proportions by cluster",
  "Cluster", "Cell proportion", "E", TRUE
)
s4f <- make_composition_plot(
  s4_context, context_colors, "Tumor and CellLine cell counts by cluster",
  "Cluster", "Number of cells", "F", FALSE
)
s4g <- make_composition_plot(
  s4_ploidy, ploidy_colors, "Initial ploidy proportions by cluster",
  "Cluster", "Cell proportion", "G", TRUE
)
s4h <- make_composition_plot(
  s4_ploidy, ploidy_colors, "Initial ploidy cell counts by cluster",
  "Cluster", "Number of cells", "H", FALSE
)
record_panel(s4a, "SuppFig4A", "panel_SuppFig4A_umap_cluster", 7.0, 5.8, "All-cell UMAP by cluster")
record_panel(s4b, "SuppFig4B", "panel_SuppFig4B_umap_initial_ploidy", 7.0, 5.8, "All-cell UMAP by initial ploidy")
record_panel(s4c, "SuppFig4C", "panel_SuppFig4C_umap_context", 7.0, 5.8, "All-cell UMAP by Tumor/CellLine context")
record_panel(s4d, "SuppFig4D", "panel_SuppFig4D_umap_s_phase_score", 7.0, 5.8, "All-cell UMAP by S phase score")
record_panel(s4e, "SuppFig4E", "panel_SuppFig4E_cluster_context_proportion", 8.5, 5.5, "Tumor and CellLine proportions by cluster")
record_panel(s4f, "SuppFig4F", "panel_SuppFig4F_cluster_context_count", 8.5, 5.5, "Tumor and CellLine counts by cluster")
record_panel(s4g, "SuppFig4G", "panel_SuppFig4G_cluster_initial_ploidy_proportion", 8.5, 5.5, "Initial ploidy proportions by cluster")
record_panel(s4h, "SuppFig4H", "panel_SuppFig4H_cluster_initial_ploidy_count", 8.5, 5.5, "Initial ploidy counts by cluster")
s4_composite <- patchwork::wrap_plots(
  patchwork::wrap_plots(s4a, s4b, s4c, s4d, ncol = 4),
  patchwork::wrap_plots(s4e, s4f, ncol = 2),
  patchwork::wrap_plots(s4g, s4h, ncol = 2),
  ncol = 1,
  heights = c(1.0, 0.95, 0.95)
) + patchwork::plot_annotation(
  title = "Supplementary Figure 4 | Tumor and CellLine cellular landscape"
)
record_panel(
  s4_composite, "SuppFig4_composite", "panel_SuppFig4_composite",
  20, 15.5, "Supplementary Figure 4 composite"
)

message("Generating Supplementary Figure 5.")
observed_dose_levels <- dose_levels[dose_levels %in% unique(as.character(tumor$dose))]
dose_panel_informative <- identical(observed_dose_levels, dose_levels)
s5a <- make_umap_discrete(
  tumor, "cluster", cluster_colors, "UMAP by cluster", "Cluster", "A",
  sprintf("Tumor cells; n = %s", format(n_tumor, big.mark = ",")), labels = TRUE
)
s5b <- make_umap_discrete(
  tumor, "initial_ploidy", ploidy_colors,
  "UMAP by initial tumor ploidy", "Initial ploidy", "B"
)
s5c <- make_umap_discrete(
  tumor, "dose", dose_colors, "UMAP by Gemcitabine dose", "Gemcitabine dose", "C"
)
s5d <- make_umap_continuous(
  tumor, "s_phase_score", "UMAP by S phase score", "S phase score", "D",
  limits = range(tumor$s_phase_score, na.rm = TRUE), diverging = TRUE
)
panel_metadata <- sample_metadata
panel_metadata$mouse_panel <- paste(
  panel_metadata$mouse,
  as.character(panel_metadata$initial_ploidy),
  as.character(panel_metadata$dose),
  sep = " | "
)
panel_lookup <- setNames(panel_metadata$mouse_panel, panel_metadata$mouse)
mouse_plot_data <- tumor
mouse_plot_data$mouse_panel <- factor(
  unname(panel_lookup[as.character(mouse_plot_data$mouse)]),
  levels = panel_metadata$mouse_panel
)
s5e <- ggplot2::ggplot(
  mouse_plot_data,
  ggplot2::aes(UMAP_1, UMAP_2, color = initial_ploidy)
) +
  ggplot2::geom_point(size = max(point_size * 1.35, 0.14), alpha = 0.88, stroke = 0) +
  ggplot2::facet_wrap(~mouse_panel, ncol = 4, drop = FALSE) +
  ggplot2::scale_color_manual(values = ploidy_colors, drop = FALSE, name = "Initial ploidy") +
  ggplot2::guides(color = umap_color_guide()) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "UMAP by mouse of origin",
    subtitle = sprintf("%s mice; shared 2N/4N color code", n_mice),
    x = "UMAP 1", y = "UMAP 2"
  ) +
  umap_theme(8.5) +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_text(size = 7.2),
    panel.spacing = grid::unit(0.10, "lines"),
    plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
  )
s5e <- add_panel_tag(s5e, "E")

composition_plot <- composition_table
composition_plot$mouse <- factor(composition_plot$mouse, levels = mouse_levels)
composition_plot$initial_ploidy <- factor(composition_plot$initial_ploidy, levels = ploidy_levels)
composition_plot$dose <- factor(composition_plot$dose, levels = dose_levels)
composition_plot$cluster_final <- factor(composition_plot$cluster_final, levels = cluster_levels)
s5f <- ggplot2::ggplot(
  composition_plot,
  ggplot2::aes(mouse, proportion, fill = cluster_final)
) +
  ggplot2::geom_col(width = 0.84, color = "white", linewidth = 0.10) +
  ggplot2::facet_grid(
    . ~ initial_ploidy + dose,
    scales = "free_x", space = "free_x", drop = FALSE
  ) +
  ggplot2::scale_fill_manual(values = cluster_colors, drop = FALSE, name = "Cluster") +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, 0.25), labels = percent_labels,
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Cluster composition of each mouse",
    subtitle = "Each mouse sums to 100%",
    x = "Mouse", y = "Cluster composition"
  ) +
  figure_theme(9.5) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 50, hjust = 1, size = 7.5))
s5f <- add_panel_tag(s5f, "F")

group_plot <- group_summary
group_plot$initial_ploidy <- factor(group_plot$initial_ploidy, levels = ploidy_levels)
group_plot$dose <- factor(group_plot$dose, levels = dose_levels)
group_plot$cluster_final <- factor(group_plot$cluster_final, levels = cluster_levels)
group_n <- unique(group_plot[, c("initial_ploidy", "dose", "n_mice"), drop = FALSE])
s5g <- ggplot2::ggplot(
  group_plot,
  ggplot2::aes(dose, mean_proportion, fill = cluster_final)
) +
  ggplot2::geom_col(width = 0.72, color = "white", linewidth = 0.12) +
  ggplot2::geom_text(
    data = group_n,
    ggplot2::aes(dose, 0.985, label = paste0("n=", n_mice)),
    inherit.aes = FALSE, size = 2.8, vjust = 1
  ) +
  ggplot2::facet_wrap(~initial_ploidy, nrow = 1, drop = FALSE) +
  ggplot2::scale_fill_manual(values = cluster_colors, drop = FALSE, name = "Cluster") +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, 0.25), labels = percent_labels,
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Mouse-weighted composition by initial ploidy and dose",
    subtitle = "Cluster proportions are calculated per mouse, then averaged with equal mouse weights",
    x = "Gemcitabine dose", y = "Mean cluster composition"
  ) +
  figure_theme(9.5) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
s5g <- add_panel_tag(s5g, "G")

if (!use_table_cache) {
  s5_dose_composition <- composition_by(
    tumor, "cluster", "dose", cluster_levels, dose_levels
  )
  s5_ploidy_composition <- composition_by(
    tumor, "cluster", "initial_ploidy", cluster_levels, ploidy_levels
  )
}
s5h <- make_composition_plot(
  s5_dose_composition, dose_colors, "Gemcitabine dose composition by cluster",
  "Cluster", "Cell proportion", "H", TRUE
)
s5i <- make_composition_plot(
  s5_ploidy_composition, ploidy_colors, "Initial ploidy composition by cluster",
  "Cluster", "Cell proportion", "I", TRUE
)
record_panel(s5a, "SuppFig5A", "panel_SuppFig5A_umap_cluster", 7.0, 5.8, "Tumor UMAP by cluster")
record_panel(s5b, "SuppFig5B", "panel_SuppFig5B_umap_initial_ploidy", 7.0, 5.8, "Tumor UMAP by initial ploidy")
record_panel(s5c, "SuppFig5C", "panel_SuppFig5C_umap_treatment_dose", 7.0, 5.8, "Tumor UMAP by Gemcitabine dose")
record_panel(s5d, "SuppFig5D", "panel_SuppFig5D_umap_s_phase_score", 7.0, 5.8, "Tumor UMAP by S phase score")
record_panel(s5e, "SuppFig5E", "panel_SuppFig5E_umap_mouse_facets", 14.0, 14.0, "Mouse-faceted UMAP colored by initial ploidy")
record_panel(s5f, "SuppFig5F", "panel_SuppFig5F_cluster_composition_by_mouse", 12.0, 6.2, "Per-mouse cluster composition")
record_panel(s5g, "SuppFig5G", "panel_SuppFig5G_mouse_weighted_cluster_composition", 10.0, 6.2, "Mouse-weighted composition by initial ploidy and dose")
record_panel(s5h, "SuppFig5H", "panel_SuppFig5H_cluster_composition_by_dose", 9.5, 5.8, "Dose composition by cluster")
record_panel(s5i, "SuppFig5I", "panel_SuppFig5I_cluster_composition_by_initial_ploidy", 9.5, 5.8, "Initial ploidy composition by cluster")
s5_left <- patchwork::wrap_plots(
  s5a, s5b, s5c, s5d,
  ncol = 1,
  heights = c(1, 1, 1, 1)
)
s5_upper <- patchwork::wrap_plots(s5_left, s5e, ncol = 2, widths = c(1, 4))
s5_composite <- patchwork::wrap_plots(
  s5_upper,
  patchwork::wrap_plots(s5f, s5g, ncol = 2, widths = c(1.2, 1)),
  patchwork::wrap_plots(s5h, s5i, ncol = 2),
  ncol = 1,
  heights = c(2.0, 0.85, 0.78)
) + patchwork::plot_annotation(
  title = "Supplementary Figure 5 | Initial ploidy, dose, and sample composition"
)
record_panel(
  s5_composite, "SuppFig5_composite", "panel_SuppFig5_composite",
  18, 23, "Supplementary Figure 5 composite"
)

message("Generating Supplementary Figure 6.")
endpoint_limits <- range(tumor$endpoint_ploidy, na.rm = TRUE)
tumor_2n <- tumor[as.character(tumor$initial_ploidy) == "2N", , drop = FALSE]
tumor_4n <- tumor[as.character(tumor$initial_ploidy) == "4N", , drop = FALSE]
s6a <- make_umap_continuous(
  tumor, "endpoint_ploidy", "Endpoint ploidy in all tumors", "Endpoint ploidy", "A",
  limits = endpoint_limits, subtitle = sprintf("Range %.3f-%.3f", endpoint_limits[[1L]], endpoint_limits[[2L]])
)
s6b <- make_umap_continuous(
  tumor_2n, "endpoint_ploidy", "Endpoint ploidy in initial 2N tumors", "Endpoint ploidy", "B",
  limits = endpoint_limits
)
s6c <- make_umap_continuous(
  tumor_4n, "endpoint_ploidy", "Endpoint ploidy in initial 4N tumors", "Endpoint ploidy", "C",
  limits = endpoint_limits
)
s6_mouse <- tumor
s6_mouse$mouse_panel <- factor(
  unname(panel_lookup[as.character(s6_mouse$mouse)]),
  levels = panel_metadata$mouse_panel
)
s6d <- ggplot2::ggplot(
  shuffle_cells(s6_mouse, 7826L),
  ggplot2::aes(UMAP_1, UMAP_2, color = endpoint_ploidy)
) +
  ggplot2::geom_point(size = max(point_size * 1.35, 0.14), alpha = 0.88, stroke = 0) +
  ggplot2::facet_wrap(~mouse_panel, ncol = 4, nrow = 4, drop = FALSE) +
  ggplot2::scale_color_gradient(
    low = "#2C7BB6", high = "#D7191C",
    limits = endpoint_limits, name = "Endpoint ploidy"
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Endpoint ploidy by mouse of origin",
    subtitle = "All panels use the same endpoint-ploidy scale",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  umap_theme(8.5) +
  ggplot2::theme(
    legend.position = "top",
    legend.justification = "right",
    legend.box.just = "right",
    strip.text = ggplot2::element_text(size = 7.2),
    panel.spacing = grid::unit(0.10, "lines")
  )
s6d <- add_panel_tag(s6d, "D")
record_panel(s6a, "SuppFig6A", "panel_SuppFig6A_umap_endpoint_ploidy_all_tumors", 7.0, 5.8, "All-tumor UMAP by endpoint ploidy")
record_panel(s6b, "SuppFig6B", "panel_SuppFig6B_umap_endpoint_ploidy_initial_2N", 7.0, 5.8, "Initial-2N tumor UMAP by endpoint ploidy")
record_panel(s6c, "SuppFig6C", "panel_SuppFig6C_umap_endpoint_ploidy_initial_4N", 7.0, 5.8, "Initial-4N tumor UMAP by endpoint ploidy")
record_panel(s6d, "SuppFig6D", "panel_SuppFig6D_umap_endpoint_ploidy_mouse_facets", 14.0, 14.0, "Mouse-faceted UMAP by endpoint ploidy")
s6_left <- patchwork::wrap_plots(
  s6a, s6b, s6c,
  ncol = 1,
  heights = c(1, 1, 1)
)
s6_composite <- patchwork::wrap_plots(
  s6_left, s6d, ncol = 2, widths = c(1, 3)
) + patchwork::plot_annotation(
  title = "Supplementary Figure 6 | Endpoint tumor ploidy"
)
record_panel(
  s6_composite, "SuppFig6_composite", "panel_SuppFig6_composite",
  20, 16.5, "Supplementary Figure 6 composite"
)

safe_cluster_stub <- function(x) {
  out <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  out <- gsub("_+", "_", out)
  sub("^_|_$", "", out)
}

clean_gene_symbol <- function(x) {
  out <- trimws(as.character(x))
  out <- sub("^GRCh[0-9]+[-_]", "", out, ignore.case = TRUE)
  out <- sub("^GRCm39[-_]", "", out, ignore.case = TRUE)
  out <- sub("^hg38[-_]", "", out, ignore.case = TRUE)
  out <- sub("\\.[0-9]+$", "", out)
  out[out == ""] <- NA_character_
  out
}

lfc_column <- function(df) {
  hit <- c("avg_log2FC", "avg_logFC", "log2FC", "logFC")
  hit <- hit[hit %in% names(df)]
  if (!length(hit)) stop("Cannot find a log-fold-change column in FindMarkers output", call. = FALSE)
  hit[[1L]]
}

hallmark_sets <- function() {
  hallmark <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  )
  split(hallmark$gene_symbol, hallmark$gs_name)
}

hallmark_label <- function(x) {
  tools::toTitleCase(tolower(gsub("_", " ", sub("^HALLMARK_", "", x))))
}

prepare_ora_markers <- function(de, lfc_col) {
  out <- de
  out$gene_symbol <- clean_gene_symbol(out$gene)
  out$gene_key <- toupper(out$gene_symbol)
  out$lfc_value <- as.numeric(out[[lfc_col]])
  out$p_val_adj_num <- as.numeric(out$p_val_adj)
  out$delta_pct <- as.numeric(out$pct.1) - as.numeric(out$pct.2)
  out$abs_logfc <- abs(out$lfc_value)
  out$abs_delta_pct <- abs(out$delta_pct)
  out <- out[
    !is.na(out$gene_key) & is.finite(out$lfc_value) &
      is.finite(out$p_val_adj_num) & is.finite(out$delta_pct) &
      out$p_val_adj_num < 0.05 & out$abs_logfc >= 0.25 &
      out$abs_delta_pct >= 0.05 & out$lfc_value > 0,
    ,
    drop = FALSE
  ]
  out <- out[order(-out$lfc_value, out$p_val_adj_num), , drop = FALSE]
  out <- out[!duplicated(out$gene_key), , drop = FALSE]
  head(out, 100L)
}

run_ora <- function(query_df, universe, pathways) {
  query <- unique(query_df$gene_symbol)
  query <- query[!is.na(query)]
  results <- lapply(names(pathways), function(pathway) {
    genes <- unique(intersect(pathways[[pathway]], universe))
    if (length(genes) < 15L || length(genes) > 500L) return(NULL)
    overlap <- intersect(query, genes)
    if (length(overlap) < 3L) return(NULL)
    p_value <- stats::phyper(
      length(overlap) - 1L,
      length(genes),
      length(universe) - length(genes),
      length(query),
      lower.tail = FALSE
    )
    data.frame(
      pathway = pathway,
      hallmark_label = hallmark_label(pathway),
      set_size = length(genes),
      query_size = length(query),
      overlap = length(overlap),
      p_value = p_value,
      overlap_genes = paste(sort(overlap), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, results)
  if (is.null(out) || !nrow(out)) return(data.frame())
  out$p_adj <- stats::p.adjust(out$p_value, method = "BH")
  out <- out[out$p_adj < 0.05, , drop = FALSE]
  if (!nrow(out)) return(out)
  overlap_stats <- lapply(strsplit(out$overlap_genes, ";", fixed = TRUE), function(genes) {
    match_rows <- query_df[toupper(query_df$gene_symbol) %in% toupper(genes), , drop = FALSE]
    data.frame(
      expression_score_raw = mean(match_rows$abs_logfc, na.rm = TRUE),
      detection_score_raw = mean(as.numeric(match_rows$pct.1), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  score <- do.call(rbind, overlap_stats)
  scale01 <- function(x) {
    if (!length(x) || !any(is.finite(x))) return(rep(0, length(x)))
    limits <- range(x[is.finite(x)])
    if (diff(limits) == 0) return(rep(1, length(x)))
    (x - limits[[1L]]) / diff(limits)
  }
  out <- cbind(out, score)
  out$overlap_score_raw <- rowMeans(cbind(
    out$overlap / pmax(1, out$query_size),
    out$overlap / pmax(1, out$set_size)
  ))
  out$expression_score_scaled <- scale01(out$expression_score_raw)
  out$detection_score_scaled <- scale01(out$detection_score_raw)
  out$overlap_score_scaled <- scale01(out$overlap_score_raw)
  out$annotation_score <- rowMeans(out[, c(
    "expression_score_scaled", "detection_score_scaled", "overlap_score_scaled"
  )])
  out <- out[order(-out$annotation_score, out$p_adj, out$pathway), , drop = FALSE]
  out$annotation_rank <- seq_len(nrow(out))
  out
}

prepare_rank_stats <- function(de, lfc_col) {
  gene <- clean_gene_symbol(de$gene)
  lfc <- as.numeric(de[[lfc_col]])
  p_value <- as.numeric(de$p_val_adj)
  finite_positive <- p_value[is.finite(p_value) & p_value > 0]
  p_floor <- if (length(finite_positive)) max(min(finite_positive) * 0.1, 1e-300) else 1e-300
  p_value[!is.finite(p_value) | p_value <= 0] <- p_floor
  rank_metric <- sign(lfc) * (-log10(p_value) + abs(lfc) * 1e-6)
  ranked <- data.frame(gene = gene, rank_metric = rank_metric, lfc = lfc)
  ranked <- ranked[
    !is.na(ranked$gene) & is.finite(ranked$rank_metric) & ranked$rank_metric != 0,
    ,
    drop = FALSE
  ]
  ranked <- ranked[order(-abs(ranked$rank_metric), -abs(ranked$lfc), ranked$gene), , drop = FALSE]
  ranked <- ranked[!duplicated(ranked$gene), , drop = FALSE]
  ranked <- ranked[order(-ranked$rank_metric, ranked$gene), , drop = FALSE]
  stats <- ranked$rank_metric
  names(stats) <- ranked$gene
  stats
}

heatmap_plot <- function(matrix_data, title, legend_title, diverging, tag) {
  pheatmap_args <- list(
    mat = matrix_data,
    cluster_rows = nrow(matrix_data) > 2L,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 8,
    fontsize_col = 7,
    angle_col = 45,
    main = paste0(tag, "  ", title),
    silent = TRUE
  )
  if (diverging) {
    max_abs <- max(abs(matrix_data), na.rm = TRUE)
    if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
    pheatmap_args$color <- grDevices::colorRampPalette(
      c("#2C7BB6", "white", "#D7191C")
    )(101)
    pheatmap_args$breaks <- seq(
      -max_abs, max_abs,
      length.out = length(pheatmap_args$color) + 1L
    )
  }
  heatmap <- do.call(pheatmap::pheatmap, pheatmap_args)
  patchwork::wrap_elements(full = heatmap$gtable)
}

load_si7_cache <- function() {
  require_package("pheatmap")
  read_matrix <- function(filename, label) {
    data <- read_tsv_checked(file.path(table_dir, filename), label)
    if (!identical(names(data), c("pathway", cluster_levels))) {
      stop(label, " does not match the canonical cluster order", call. = FALSE)
    }
    pathways <- as.character(data$pathway)
    matrix_data <- as.matrix(data[, cluster_levels, drop = FALSE])
    storage.mode(matrix_data) <- "double"
    rownames(matrix_data) <- pathways
    if (nrow(matrix_data) != 20L || any(!is.finite(matrix_data))) {
      stop(label, " must be a finite 20 x ", length(cluster_levels), " matrix", call. = FALSE)
    }
    matrix_data
  }
  ora_matrix <- read_matrix(
    "si_figure7_cluster_Hallmark_ORA_annotation_score_heatmap_top20_matrix.tsv",
    "cached SI Figure 7 ORA matrix"
  )
  gsea_matrix <- read_matrix(
    "si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv",
    "cached SI Figure 7 GSEA matrix"
  )
  list(
    ora = heatmap_plot(
      ora_matrix, "Cluster Hallmark ORA annotation score", "Annotation score", FALSE, "A"
    ),
    gsea = heatmap_plot(
      gsea_matrix, "Cluster Hallmark GSEA NES", "NES", TRUE, "B"
    ),
    nperm_simple = "cached",
    detected_cpus = "not_run",
    reserved_cpus = "not_run",
    usable_cpus = "not_run",
    cluster_workers = "not_run",
    inner_workers = "not_run",
    memory_limit_bytes = "not_run",
    object_size_bytes = "not_run",
    memory_per_worker_bytes = "not_run",
    memory_worker_budget = "not_run",
    execution_worker_budget = "not_run",
    future_global_limit_bytes = "not_run",
    deg_cache_files = character(0),
    deg_cache_reused_clusters = length(cluster_levels),
    n_object_cells = nrow(canonical_cells),
    n_genes = "not_loaded"
  )
}

parse_memory_bytes <- function(value, default_unit = 1) {
  if (is.null(value) || !length(value)) return(NA_real_)
  value <- trimws(as.character(value[[1L]]))
  if (!nzchar(value) || tolower(value) == "max") return(NA_real_)
  match <- regexec("^([0-9]+[.]?[0-9]*)([KMGT]?)$", toupper(value))
  pieces <- regmatches(toupper(value), match)[[1L]]
  if (!length(pieces)) return(NA_real_)
  multiplier <- switch(
    pieces[[3L]],
    K = 1024,
    M = 1024^2,
    G = 1024^3,
    T = 1024^4,
    default_unit
  )
  as.numeric(pieces[[2L]]) * multiplier
}

detect_memory_limit_bytes <- function(detected_cpus) {
  candidates <- numeric(0)
  slurm_node <- parse_memory_bytes(Sys.getenv("SLURM_MEM_PER_NODE"), 1024^2)
  if (is.finite(slurm_node) && slurm_node > 0) {
    candidates <- c(candidates, slurm_node)
  }
  slurm_cpu <- parse_memory_bytes(Sys.getenv("SLURM_MEM_PER_CPU"), 1024^2)
  if (is.finite(slurm_cpu) && slurm_cpu > 0) {
    candidates <- c(candidates, slurm_cpu * detected_cpus)
  }
  for (path in c("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes")) {
    if (!file.exists(path)) next
    value <- tryCatch(readLines(path, n = 1L, warn = FALSE), error = function(e) NA_character_)
    parsed <- parse_memory_bytes(value, 1)
    if (is.finite(parsed) && parsed > 0 && parsed < 2^60) {
      candidates <- c(candidates, parsed)
    }
  }
  meminfo <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) character(0))
  memtotal <- grep("^MemTotal:", meminfo, value = TRUE)
  if (length(memtotal)) {
    parsed <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", memtotal[[1L]]))) * 1024
    if (is.finite(parsed) && parsed > 0) candidates <- c(candidates, parsed)
  }
  if (!length(candidates)) return(NA_real_)
  min(candidates)
}

find_cluster_markers_parallel <- function(
  cluster_index,
  cluster_levels,
  inner_workers,
  multicore_enabled,
  seurat_object,
  table_dir
) {
  cluster_id <- cluster_levels[[cluster_index]]
  assigned_inner_workers <- inner_workers[[cluster_index]]
  worker_previous_plan <- future::plan()
  on.exit(future::plan(worker_previous_plan), add = TRUE)
  if (multicore_enabled && assigned_inner_workers > 1L) {
    future::plan(future::multicore, workers = assigned_inner_workers)
  } else {
    future::plan(future::sequential)
  }
  message("FindMarkers: cluster ", cluster_id, " versus rest")
  de <- Seurat::FindMarkers(
    object = seurat_object,
    ident.1 = cluster_id,
    ident.2 = NULL,
    assay = "RNA",
    slot = "data",
    verbose = FALSE
  )
  de$gene <- rownames(de)
  de$gene_symbol <- clean_gene_symbol(de$gene)
  de$cluster <- cluster_id
  lfc_col <- lfc_column(de)
  de <- de[order(de$p_val_adj, -abs(de[[lfc_col]])), , drop = FALSE]
  rownames(de) <- NULL
  write_csv(
    de,
    file.path(
      table_dir,
      paste0("si_figure7_cluster_", safe_cluster_stub(cluster_id), "_vs_rest_DEG.csv")
    )
  )
  list(cluster = cluster_id, de = de)
}

read_cached_cluster_markers <- function(cluster_id, cache_dir, table_dir) {
  cache_path <- file.path(
    cache_dir,
    paste0("si_figure7_cluster_", safe_cluster_stub(cluster_id), "_vs_rest_DEG.csv")
  )
  de <- read_csv_checked(cache_path, paste0("SI Figure 7 DEG cache for cluster ", cluster_id))
  required <- c("gene", "gene_symbol", "cluster", "p_val_adj", "pct.1", "pct.2")
  missing <- setdiff(required, names(de))
  if (length(missing)) {
    stop(
      "SI Figure 7 DEG cache is missing columns for cluster ", cluster_id, ": ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  lfc_column(de)
  if (anyNA(de$gene) || any(!nzchar(de$gene)) || anyDuplicated(de$gene)) {
    stop("SI Figure 7 DEG cache has invalid gene identifiers for cluster ", cluster_id, call. = FALSE)
  }
  cached_clusters <- unique(clean_character(de$cluster))
  if (!identical(cached_clusters, cluster_id)) {
    stop(
      "SI Figure 7 DEG cache cluster mismatch for ", cluster_id, ": ",
      paste(cached_clusters, collapse = ","),
      call. = FALSE
    )
  }
  target_path <- file.path(
    table_dir,
    paste0("si_figure7_cluster_", safe_cluster_stub(cluster_id), "_vs_rest_DEG.csv")
  )
  copied <- file.copy(cache_path, target_path, overwrite = TRUE, copy.mode = TRUE)
  if (!copied || !file.exists(target_path) || file.info(target_path)$size <= 0) {
    stop("Failed to materialize SI Figure 7 DEG cache for cluster ", cluster_id, call. = FALSE)
  }
  if (!identical(file_sha256(cache_path), file_sha256(target_path))) {
    stop("SI Figure 7 DEG cache hash mismatch after copy for cluster ", cluster_id, call. = FALSE)
  }
  list(cluster = cluster_id, de = de, cache_path = normalizePath(cache_path, mustWork = TRUE))
}

run_si7 <- function() {
  require_package("Seurat")
  require_package("future")
  require_package("future.apply")
  require_package("fgsea")
  require_package("msigdbr")
  require_package("pheatmap")
  previous_future_plan <- future::plan()
  previous_future_max_size <- getOption("future.globals.maxSize")
  on.exit(
    {
      future::plan(previous_future_plan)
      options(future.globals.maxSize = previous_future_max_size)
    },
    add = TRUE
  )
  options(future.globals.maxSize = 8 * 1024^3)
  message("Generating Supplementary Figure 7 from raw Seurat RDS: ", seurat_rds_path)
  object <- readRDS(seurat_rds_path)
  if (!inherits(object, "Seurat")) stop("--seurat-rds is not a Seurat object", call. = FALSE)
  future_global_limit_bytes <- max(
    8 * 1024^3,
    as.numeric(utils::object.size(object)) * 3 + 2 * 1024^3
  )
  options(future.globals.maxSize = future_global_limit_bytes)
  if (!(canonical_cluster_col %in% names(object@meta.data))) {
    stop("Raw Seurat object is missing cluster field: ", canonical_cluster_col, call. = FALSE)
  }
  if (ncol(object) != nrow(seurat)) {
    stop(
      "Raw Seurat RDS cell count differs from seurat_metadata.csv: ",
      ncol(object), " versus ", nrow(seurat),
      call. = FALSE
    )
  }
  if (!("RNA" %in% names(object@assays))) stop("Raw Seurat object is missing RNA assay", call. = FALSE)
  Seurat::DefaultAssay(object) <- "RNA"
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    object <- tryCatch(
      SeuratObject::JoinLayers(object, assay = "RNA"),
      error = function(e) object
    )
  }
  data_matrix <- tryCatch(
    SeuratObject::LayerData(object, assay = "RNA", layer = "data"),
    error = function(e) tryCatch(
      Seurat::GetAssayData(object, assay = "RNA", slot = "data"),
      error = function(e2) NULL
    )
  )
  if (is.null(data_matrix) || !nrow(data_matrix) || !ncol(data_matrix)) {
    object <- Seurat::NormalizeData(object, assay = "RNA", verbose = FALSE)
  }
  rds_clusters <- sort_cluster_levels(object@meta.data[[canonical_cluster_col]])
  if (!identical(rds_clusters, cluster_levels)) {
    stop(
      "Raw RDS and Seurat metadata cluster levels differ: RDS=",
      paste(rds_clusters, collapse = ","), "; metadata=", paste(cluster_levels, collapse = ","),
      call. = FALSE
    )
  }
  Seurat::Idents(object) <- factor(
    as.character(object@meta.data[[canonical_cluster_col]]),
    levels = cluster_levels
  )
  detected_cpus <- max(1L, as.integer(future::availableCores()))
  reserved_cpus <- min(2L, max(0L, detected_cpus - 1L))
  usable_cpus <- max(1L, detected_cpus - reserved_cpus)
  object_size_bytes <- as.numeric(utils::object.size(object))
  memory_limit_bytes <- detect_memory_limit_bytes(detected_cpus)
  memory_per_worker_bytes <- object_size_bytes * 1.5
  memory_worker_budget <- if (is.finite(memory_limit_bytes)) {
    max(
      1L,
      as.integer(floor(
        max(0, memory_limit_bytes * 0.90 - object_size_bytes) /
          max(memory_per_worker_bytes, 1)
      ))
    )
  } else {
    usable_cpus
  }
  execution_worker_budget <- max(
    1L,
    min(usable_cpus, memory_worker_budget)
  )
  multicore_enabled <- future::supportsMulticore() && execution_worker_budget > 1L
  cluster_workers <- if (multicore_enabled) {
    min(length(cluster_levels), execution_worker_budget)
  } else {
    1L
  }
  active_inner_workers <- rep(
    execution_worker_budget %/% cluster_workers,
    cluster_workers
  )
  remainder <- execution_worker_budget %% cluster_workers
  if (remainder > 0L) {
    active_inner_workers[seq_len(remainder)] <-
      active_inner_workers[seq_len(remainder)] + 1L
  }
  inner_workers <- rep(active_inner_workers, length.out = length(cluster_levels))
  names(inner_workers) <- cluster_levels
  message(
    "SI Figure 7 differential expression layout: detected_cpus=", detected_cpus,
    "; reserved_cpus=", reserved_cpus,
    "; usable_cpus=", usable_cpus,
    "; memory_limit_bytes=", format(memory_limit_bytes, scientific = FALSE),
    "; object_size_bytes=", format(object_size_bytes, scientific = FALSE),
    "; memory_per_worker_bytes=", format(memory_per_worker_bytes, scientific = FALSE),
    "; memory_worker_budget=", memory_worker_budget,
    "; execution_worker_budget=", execution_worker_budget,
    "; cluster_workers=", cluster_workers,
    "; inner_workers=", paste(inner_workers, collapse = ",")
  )

  if (!is.na(deg_cache_dir)) {
    message("Reusing complete SI Figure 7 DEG cache: ", deg_cache_dir)
    deg_results <- lapply(
      cluster_levels,
      read_cached_cluster_markers,
      cache_dir = deg_cache_dir,
      table_dir = table_dir
    )
  } else if (cluster_workers > 1L) {
    future::plan(future::multicore, workers = cluster_workers)
    deg_results <- future.apply::future_lapply(
      seq_along(cluster_levels),
      find_cluster_markers_parallel,
      cluster_levels = cluster_levels,
      inner_workers = inner_workers,
      multicore_enabled = multicore_enabled,
      seurat_object = object,
      table_dir = table_dir,
      future.seed = TRUE,
      future.scheduling = 1
    )
  } else {
    future::plan(future::sequential)
    deg_results <- lapply(
      seq_along(cluster_levels),
      find_cluster_markers_parallel,
      cluster_levels = cluster_levels,
      inner_workers = inner_workers,
      multicore_enabled = multicore_enabled,
      seurat_object = object,
      table_dir = table_dir
    )
  }
  deg_list <- stats::setNames(
    lapply(deg_results, function(result) result$de),
    vapply(deg_results, function(result) result$cluster, character(1))
  )
  deg_cache_files <- if (is.na(deg_cache_dir)) {
    character(0)
  } else {
    vapply(deg_results, function(result) result$cache_path, character(1))
  }
  universe <- sort(unique(unlist(lapply(deg_list, function(de) de$gene_symbol))))
  universe <- universe[!is.na(universe) & nzchar(universe)]
  write_csv(
    data.frame(gene_symbol = universe, stringsAsFactors = FALSE),
    file.path(table_dir, "si_figure7_hallmark_ORA_universe.csv")
  )
  pathways <- hallmark_sets()
  ora_list <- list()
  gsea_list <- list()
  nperm_simple <- as.integer(arg_value(args, "fgsea-nperm-simple", "10000"))
  if (!is.finite(nperm_simple) || nperm_simple < 1L) {
    stop("--fgsea-nperm-simple must be a positive integer", call. = FALSE)
  }
  for (cluster_id in cluster_levels) {
    de <- deg_list[[cluster_id]]
    lfc_col <- lfc_column(de)
    ora_markers <- prepare_ora_markers(de, lfc_col)
    write_csv(
      ora_markers,
      file.path(
        table_dir,
        paste0("si_figure7_cluster_", safe_cluster_stub(cluster_id), "_top100_up_ORA_input.csv")
      )
    )
    ora <- run_ora(ora_markers, universe, pathways)
    if (nrow(ora)) ora$cluster <- cluster_id
    ora_list[[cluster_id]] <- ora
    stats <- prepare_rank_stats(de, lfc_col)
    gsea <- fgsea::fgseaMultilevel(
      pathways = pathways,
      stats = stats,
      minSize = 15L,
      maxSize = 500L,
      nPermSimple = nperm_simple,
      eps = 0
    )
    gsea <- as.data.frame(gsea, stringsAsFactors = FALSE)
    if ("leadingEdge" %in% names(gsea)) {
      gsea$leadingEdge <- vapply(
        gsea$leadingEdge,
        function(x) paste(as.character(x), collapse = ";"),
        character(1)
      )
    }
    gsea$hallmark_label <- hallmark_label(gsea$pathway)
    gsea$cluster <- cluster_id
    gsea <- gsea[order(gsea$padj, -abs(gsea$NES), gsea$pathway), , drop = FALSE]
    gsea_list[[cluster_id]] <- gsea
  }
  ora_nonempty <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, ora_list)
  ora_all <- if (length(ora_nonempty)) do.call(rbind, ora_nonempty) else data.frame()
  gsea_all <- do.call(rbind, gsea_list)
  if (is.null(ora_all) || !nrow(ora_all)) stop("No significant Hallmark ORA results", call. = FALSE)
  if (is.null(gsea_all) || !nrow(gsea_all)) stop("No Hallmark GSEA results", call. = FALSE)
  rownames(ora_all) <- NULL
  rownames(gsea_all) <- NULL
  write_csv(ora_all, file.path(table_dir, "si_figure7_cluster_Hallmark_ORA_all.csv"))
  write_csv(gsea_all, file.path(table_dir, "si_figure7_cluster_Hallmark_GSEA_all.csv"))

  ora_top <- aggregate(annotation_score ~ hallmark_label, data = ora_all, FUN = max)
  ora_top <- ora_top[order(-ora_top$annotation_score, ora_top$hallmark_label), , drop = FALSE]
  ora_pathways <- head(ora_top$hallmark_label, 20L)
  ora_matrix <- matrix(
    0,
    nrow = length(ora_pathways),
    ncol = length(cluster_levels),
    dimnames = list(ora_pathways, cluster_levels)
  )
  for (i in seq_len(nrow(ora_all))) {
    row_id <- match(ora_all$hallmark_label[[i]], rownames(ora_matrix))
    col_id <- match(ora_all$cluster[[i]], colnames(ora_matrix))
    if (!is.na(row_id) && !is.na(col_id)) {
      ora_matrix[row_id, col_id] <- max(
        ora_matrix[row_id, col_id], ora_all$annotation_score[[i]], na.rm = TRUE
      )
    }
  }
  gsea_heatmap <- gsea_all[is.finite(gsea_all$NES), , drop = FALSE]
  if (!nrow(gsea_heatmap)) {
    stop("No finite Hallmark GSEA NES values are available for the heatmap", call. = FALSE)
  }
  gsea_heatmap$abs_NES <- abs(gsea_heatmap$NES)
  gsea_top <- aggregate(abs_NES ~ hallmark_label, data = gsea_heatmap, FUN = max)
  names(gsea_top)[2L] <- "best_abs_nes"
  gsea_top <- gsea_top[order(-gsea_top$best_abs_nes, gsea_top$hallmark_label), , drop = FALSE]
  gsea_pathways <- head(gsea_top$hallmark_label, 20L)
  gsea_matrix <- matrix(
    0,
    nrow = length(gsea_pathways),
    ncol = length(cluster_levels),
    dimnames = list(gsea_pathways, cluster_levels)
  )
  for (i in seq_len(nrow(gsea_heatmap))) {
    row_id <- match(gsea_heatmap$hallmark_label[[i]], rownames(gsea_matrix))
    col_id <- match(gsea_heatmap$cluster[[i]], colnames(gsea_matrix))
    if (!is.na(row_id) && !is.na(col_id)) {
      current <- gsea_matrix[row_id, col_id]
      candidate <- gsea_heatmap$NES[[i]]
      if (abs(candidate) >= abs(current)) gsea_matrix[row_id, col_id] <- candidate
    }
  }
  write_tsv(
    data.frame(pathway = rownames(ora_matrix), ora_matrix, check.names = FALSE),
    file.path(table_dir, "si_figure7_cluster_Hallmark_ORA_annotation_score_heatmap_top20_matrix.tsv")
  )
  write_tsv(
    data.frame(pathway = rownames(gsea_matrix), gsea_matrix, check.names = FALSE),
    file.path(table_dir, "si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv")
  )
  list(
    ora = heatmap_plot(
      ora_matrix, "Cluster Hallmark ORA annotation score", "Annotation score", FALSE, "A"
    ),
    gsea = heatmap_plot(
      gsea_matrix, "Cluster Hallmark GSEA NES", "NES", TRUE, "B"
    ),
    nperm_simple = nperm_simple,
    detected_cpus = detected_cpus,
    reserved_cpus = reserved_cpus,
    usable_cpus = usable_cpus,
    cluster_workers = cluster_workers,
    inner_workers = paste(inner_workers, collapse = ","),
    memory_limit_bytes = memory_limit_bytes,
    object_size_bytes = object_size_bytes,
    memory_per_worker_bytes = memory_per_worker_bytes,
    memory_worker_budget = memory_worker_budget,
    execution_worker_budget = execution_worker_budget,
    future_global_limit_bytes = future_global_limit_bytes,
    deg_cache_files = deg_cache_files,
    deg_cache_reused_clusters = length(deg_cache_files),
    n_object_cells = ncol(object),
    n_genes = nrow(object)
  )
}

si7_result <- NULL
if (!skip_si7) {
  si7_result <- if (use_table_cache) load_si7_cache() else run_si7()
  record_panel(
    si7_result$ora, "SuppFig7A",
    "panel_SuppFig7A_cluster_Hallmark_ORA_annotation_score_heatmap_top20",
    10, 8, "Cluster Hallmark ORA annotation-score heatmap"
  )
  record_panel(
    si7_result$gsea, "SuppFig7B",
    "panel_SuppFig7B_cluster_Hallmark_GSEA_NES_heatmap_top20",
    10, 8, "Cluster Hallmark GSEA NES heatmap"
  )
  s7_composite <- patchwork::wrap_plots(
    si7_result$ora, si7_result$gsea, ncol = 2
  ) + patchwork::plot_annotation(
    title = "Supplementary Figure 7 | Cluster Hallmark pathway analysis"
  )
  record_panel(
    s7_composite, "SuppFig7_composite", "panel_SuppFig7_composite",
    20, 8.8, "Supplementary Figure 7 composite"
  )
}

panel_contract <- do.call(rbind, panel_rows)
rownames(panel_contract) <- NULL
observed_figure_files <- sort(list.files(figure_dir, pattern = "[.](pdf|png)$"))
expected_figure_files <- sort(panel_contract$filename)
if (!identical(observed_figure_files, expected_figure_files)) {
  stop(
    "Supplementary Figures exact panel inventory failed; missing=",
    paste(setdiff(expected_figure_files, observed_figure_files), collapse = ","),
    "; unexpected=",
    paste(setdiff(observed_figure_files, expected_figure_files), collapse = ","),
    call. = FALSE
  )
}
write_tsv(panel_contract, file.path(metadata_dir, "panel_contract.tsv"))

canonical_table_path <- file.path(table_dir, "si_figures_cell_metadata.csv")
cluster_key_path <- file.path(table_dir, "si_figures_cluster_key.tsv")
if (!use_table_cache) {
  write_csv(canonical_cells, canonical_table_path)
  write_tsv(cluster_key, cluster_key_path)
  write_csv(s4_context, file.path(table_dir, "si_figure4_cluster_context_composition.csv"))
  write_csv(s4_ploidy, file.path(table_dir, "si_figure4_cluster_initial_ploidy_composition.csv"))
  write_csv(composition_table, file.path(table_dir, "si_figure5_cluster_composition_by_mouse.csv"))
  group_summary_output <- group_summary
  group_summary_output$initial_ploidy <- as.character(group_summary_output$initial_ploidy)
  group_summary_output$dose <- as.character(group_summary_output$dose)
  group_summary_output$cluster_final <- as.character(group_summary_output$cluster_final)
  write_csv(
    group_summary_output,
    file.path(table_dir, "si_figure5_mouse_weighted_composition_by_initial_ploidy_dose.csv")
  )
  write_csv(s5_dose_composition, file.path(table_dir, "si_figure5_cluster_dose_composition.csv"))
  write_csv(s5_ploidy_composition, file.path(table_dir, "si_figure5_cluster_initial_ploidy_composition.csv"))
  endpoint_audit <- data.frame(
    cell = seurat$cell,
    context = seurat$context,
    initial_ploidy = seurat$initial_ploidy,
    endpoint_file = seurat$endpoint_file,
    endpoint_cell_id = seurat$endpoint_cell_id,
    endpoint_ploidy = seurat$endpoint_ploidy,
    matched = is.finite(seurat$endpoint_ploidy),
    stringsAsFactors = FALSE
  )
  write_csv(endpoint_audit, file.path(table_dir, "si_figure6_endpoint_ploidy_join_audit.csv"))
}
table_validation_status <- system2(
  "python3",
  c(shQuote(cache_validator), "--cache-dir", shQuote(table_dir))
)
if (!identical(table_validation_status, 0L)) {
  stop("Generated SI Figures tables failed the canonical cache contract", call. = FALSE)
}

if (use_table_cache) {
  input_roles <- c(
    "figure7_config",
    paste0("si_figures_table_cache_", basename(cache_files))
  )
  input_paths <- c(config_path, cache_files)
  input_rows <- c(NA_integer_, rep(NA_integer_, length(cache_files)))
} else {
  input_roles <- c("seurat_metadata", "scvelo_cell_metrics", "all_ploidy", "figure7_config")
  input_paths <- c(seurat_metadata_path, scvelo_metrics_path, all_ploidy_path, config_path)
  input_rows <- c(nrow(seurat), nrow(metrics), nrow(endpoint_ploidy), NA_integer_)
}
if (!skip_si7 && !use_table_cache) {
  input_roles <- c(input_roles, "seurat_rds")
  input_paths <- c(input_paths, seurat_rds_path)
  input_rows <- c(input_rows, si7_result$n_object_cells)
  if (length(si7_result$deg_cache_files)) {
    input_roles <- c(
      input_roles,
      paste0("si_figure7_deg_cache_cluster_", cluster_levels)
    )
    input_paths <- c(input_paths, si7_result$deg_cache_files)
    input_rows <- c(
      input_rows,
      vapply(si7_result$deg_cache_files, function(path) {
        nrow(read_csv_checked(path, "SI Figure 7 DEG cache manifest input"))
      }, integer(1))
    )
  }
}
if (!use_table_cache && file.exists(seurat_metadata_provenance_path)) {
  input_roles <- c(input_roles, "seurat_metadata_provenance")
  input_paths <- c(input_paths, seurat_metadata_provenance_path)
  input_rows <- c(input_rows, length(input_provenance))
}
input_manifest <- data.frame(
  role = input_roles,
  path = normalizePath(input_paths, mustWork = TRUE),
  sha256 = vapply(input_paths, file_sha256, character(1)),
  bytes = as.numeric(file.info(input_paths)$size),
  rows = input_rows,
  stringsAsFactors = FALSE
)
write_csv(input_manifest, file.path(metadata_dir, "input_manifest.csv"))

run_config <- data.frame(
  key = c(
    "module", "figures", "logical_panels", "figure_file_count",
    "plot_shuffle_seed", "cluster_field", "sample_field", "umap_reduction",
    "endpoint_ploidy_join_key", "endpoint_ploidy_color_scale",
    "si_figure7_from_raw_rds", "fgsea_nperm_simple",
    "si_figure7_detected_cpus", "si_figure7_reserved_cpus",
    "si_figure7_usable_cpus", "si_figure7_cluster_workers",
    "si_figure7_inner_workers", "si_figure7_future_global_limit_bytes",
    "si_figure7_memory_limit_bytes", "si_figure7_object_size_bytes",
    "si_figure7_memory_per_worker_bytes",
    "si_figure7_memory_worker_budget", "si_figure7_execution_worker_budget",
    "si_figure7_deg_cache_reused_clusters",
    "table_mode", "canonical_table_root", "cache_file_count",
    "force_reanalysis", "raw_seurat_loaded", "deg_analysis_executed"
  ),
  value = c(
    "si_figures", if (skip_si7) "4,5,6" else "4,5,6,7",
    as.character(nrow(panel_contract) / 2L), as.character(nrow(panel_contract)),
    "5826", canonical_cluster_col, sample_col, as.character(si_config$umap_reduction),
    "paste0(harvest,'.sps.cbs')+barcode_raw -> file+cell_id",
    sprintf("blue-red; shared range %.15g to %.15g", endpoint_limits[[1L]], endpoint_limits[[2L]]),
    tolower(as.character(!skip_si7 && !use_table_cache)),
    if (skip_si7) "not_run" else as.character(si7_result$nperm_simple),
    if (skip_si7) "not_run" else as.character(si7_result$detected_cpus),
    if (skip_si7) "not_run" else as.character(si7_result$reserved_cpus),
    if (skip_si7) "not_run" else as.character(si7_result$usable_cpus),
    if (skip_si7) "not_run" else as.character(si7_result$cluster_workers),
    if (skip_si7) "not_run" else si7_result$inner_workers,
    if (skip_si7) "not_run" else format(si7_result$future_global_limit_bytes, scientific = FALSE),
    if (skip_si7) "not_run" else format(si7_result$memory_limit_bytes, scientific = FALSE),
    if (skip_si7) "not_run" else format(si7_result$object_size_bytes, scientific = FALSE),
    if (skip_si7) "not_run" else format(si7_result$memory_per_worker_bytes, scientific = FALSE),
    if (skip_si7) "not_run" else as.character(si7_result$memory_worker_budget),
    if (skip_si7) "not_run" else as.character(si7_result$execution_worker_budget),
    if (skip_si7) "not_run" else as.character(si7_result$deg_cache_reused_clusters),
    if (use_table_cache) "canonical_cache" else "full_reanalysis",
    table_cache_dir,
    if (use_table_cache) as.character(length(cache_files)) else "0",
    tolower(as.character(force_reanalysis)),
    tolower(as.character(!skip_si7 && !use_table_cache)),
    tolower(as.character(!skip_si7 && !use_table_cache && is.na(deg_cache_dir)))
  ),
  stringsAsFactors = FALSE
)
write_tsv(run_config, file.path(metadata_dir, "run_config.tsv"))

qc <- data.frame(
  key = c(
    "generated_at", "seurat_cells", "scvelo_cells", "seurat_cells_in_scvelo",
    "scvelo_cells_missing_from_seurat", "tumor_cells", "cellline_cells",
    "tumor_cells_with_endpoint_ploidy", "tumor_cells_missing_endpoint_ploidy",
    "cellline_cells_with_endpoint_ploidy", "endpoint_ploidy_min",
    "endpoint_ploidy_max", "endpoint_ploidy_input_rows", "tumor_mice",
    "clusters", "initial_ploidy_levels", "dose_levels", "dose_panel_informative",
    "panel_files", "si_figure7_rds_cells"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    nrow(seurat), nrow(metrics), sum(seurat$in_scvelo), length(missing_metrics_cells),
    sum(seurat$context == "Tumor"), sum(seurat$context == "CellLine"),
    sum(seurat$context == "Tumor" & is.finite(seurat$endpoint_ploidy)),
    sum(tumor_endpoint_missing), sum(cellline_endpoint_present),
    format(endpoint_limits[[1L]], digits = 16),
    format(endpoint_limits[[2L]], digits = 16),
    nrow(endpoint_ploidy), n_mice, paste(cluster_levels, collapse = ","),
    paste(ploidy_levels, collapse = ","), paste(observed_dose_levels, collapse = ","),
    tolower(as.character(dose_panel_informative)), nrow(panel_contract),
    if (skip_si7) "not_run" else as.character(si7_result$n_object_cells)
  ),
  stringsAsFactors = FALSE
)
write_tsv(qc, file.path(metadata_dir, "input_qc.tsv"))

provenance <- data.frame(
  key = c(
    "schema_version", "artifact", "script", "script_sha256",
    "seurat_metadata", "seurat_metadata_sha256", "scvelo_metrics",
    "scvelo_metrics_sha256", "all_ploidy", "all_ploidy_sha256",
    "seurat_rds", "seurat_rds_sha256", "source_provenance_present",
    "figure7_config", "figure7_config_sha256",
    "si_figure7_traced_source", "si_figure7_deg_method",
    "si_figure7_ora_method", "si_figure7_gsea_method",
    "source_code_revision"
  ),
  value = c(
    "2", "supplementary_figures_4_7", script_path, file_sha256(script_path),
    if (use_table_cache) "not_read_in_plot_only_mode" else seurat_metadata_path,
    if (use_table_cache) "not_read_in_plot_only_mode" else file_sha256(seurat_metadata_path),
    if (use_table_cache) "not_read_in_plot_only_mode" else scvelo_metrics_path,
    if (use_table_cache) "not_read_in_plot_only_mode" else file_sha256(scvelo_metrics_path),
    if (use_table_cache) "not_read_in_plot_only_mode" else all_ploidy_path,
    if (use_table_cache) "not_read_in_plot_only_mode" else file_sha256(all_ploidy_path),
    if (skip_si7 || use_table_cache) "not_run" else seurat_rds_path,
    if (skip_si7 || use_table_cache) "not_run" else file_sha256(seurat_rds_path),
    tolower(as.character(!use_table_cache && file.exists(seurat_metadata_provenance_path))),
    config_path, file_sha256(config_path),
    "Code/in-vivo/03a_DEGs.R + Code/in-vivo/03b_cluster_annotation_and_GSEA.R",
    "Seurat::FindMarkers cluster versus rest; RNA assay; default thresholds and test",
    "top-100 up markers; Hallmark hypergeometric ORA; BH FDR < 0.05; equal expression/detection/overlap annotation weights",
    if (skip_si7) "not_run" else paste0(
      "fgseaMultilevel Hallmark preranked; minSize=15; maxSize=500; nPermSimple=",
      si7_result$nperm_simple, "; eps=0"
    ),
    tryCatch(
      system2("git", c("-C", shQuote(repo_root), "rev-parse", "HEAD"), stdout = TRUE)[[1L]],
      error = function(e) "unavailable"
    )
  ),
  stringsAsFactors = FALSE
)
write_tsv(provenance, file.path(metadata_dir, "si_figures_provenance.tsv"))
writeLines(capture.output(sessionInfo()), file.path(metadata_dir, "sessionInfo.txt"))
writeLines(
  c(
    "Supplementary Figures 4-7 generation completed.",
    paste0("Output directory: ", output_dir),
    paste0("Panel files: ", nrow(panel_contract)),
    paste0("Seurat cells: ", nrow(seurat)),
    paste0("Tumor cells: ", nrow(tumor)),
    paste0("CellLine cells: ", sum(seurat$context == "CellLine")),
    paste0("Tumor endpoint-ploidy matches: ", sum(is.finite(tumor$endpoint_ploidy)), "/", nrow(tumor)),
    paste0("Endpoint ploidy range: ", paste(format(endpoint_limits, digits = 16), collapse = " to ")),
    paste0("SI Figure 7 run: ", !skip_si7)
  ),
  file.path(metadata_dir, "run_summary.txt")
)

message("Supplementary Figures completed: ", output_dir)
