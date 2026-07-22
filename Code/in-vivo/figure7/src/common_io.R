# Shared I/O and output-contract helpers for the narrow Figure 7 module.

figure7_stop <- function(...) stop(..., call. = FALSE)

figure7_parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) figure7_stop("Unexpected positional argument: ", arg)
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L || !nzchar(pieces[[1L]])) {
      figure7_stop("Arguments must use --name=value syntax: ", arg)
    }
    out[[gsub("-", "_", pieces[[1L]])]] <- pieces[[2L]]
  }
  out
}

figure7_arg <- function(args, name, default = NULL, required = FALSE) {
  value <- args[[gsub("-", "_", name)]]
  if (is.null(value) || !nzchar(value)) value <- default
  if (isTRUE(required) && (is.null(value) || !nzchar(value))) {
    figure7_stop("Missing required argument --", gsub("_", "-", name))
  }
  value
}

figure7_sha256 <- function(path) {
  if (!file.exists(path)) figure7_stop("Cannot checksum missing file: ", path)
  if (!requireNamespace("digest", quietly = TRUE)) figure7_stop("R package 'digest' is required")
  unname(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

figure7_read_config <- function(path, tgi_day = NULL) {
  if (!file.exists(path)) figure7_stop("Missing Figure 7 config: ", path)
  if (!requireNamespace("yaml", quietly = TRUE)) figure7_stop("R package 'yaml' is required")
  config <- yaml::read_yaml(path)
  required <- c("schema_version", "module", "raw_data", "inputs", "tgi", "statistics", "etp", "intervals", "state_pathways", "panels")
  missing <- setdiff(required, names(config))
  if (length(missing)) figure7_stop("Config is missing section(s): ", paste(missing, collapse = ", "))
  if (!identical(as.character(config$module), "in_vivo_figure7")) figure7_stop("Unexpected config module")
  if (!identical(as.character(config$raw_data$doi), "10.5281/zenodo.21463392") ||
      !identical(as.integer(config$raw_data$required_loom_files), 18L) ||
      !identical(as.character(config$raw_data$seurat_rds_sha256),
                 "727b8a5e5868da911c3b0873838fb1b0023377ed21ea5498dc6493acbbef6d98")) {
    figure7_stop("Figure 7 raw-data Zenodo contract is not the reviewed record")
  }
  configured_day <- suppressWarnings(as.integer(config$tgi$day))
  if (length(configured_day) != 1L || !is.finite(configured_day) || configured_day < 0L ||
      !identical(as.character(config$tgi$outcome), "day") ||
      !identical(as.character(config$tgi$matched_control_summary), "mean") ||
      !identical(as.character(config$tgi$matched_control_group), "initial_ploidy")) {
    figure7_stop("Figure 7 requires endpoint-day TGI and the mean initial-ploidy-matched control reference")
  }
  if (!is.null(tgi_day)) {
    selected_day <- suppressWarnings(as.integer(tgi_day))
    if (length(selected_day) != 1L || !is.finite(selected_day) || selected_day < 0L) {
      figure7_stop("--tgi-day must be one non-negative integer")
    }
    config$tgi$day <- selected_day
  }
  if (!isTRUE(all.equal(as.numeric(config$etp$threshold), 2.24, tolerance = 0))) {
    figure7_stop("Figure 7 requires the reference-balanced ETP threshold 2.24")
  }
  if (!identical(as.integer(unlist(config$panels$selected_direct_comparison_ids)), c(1L, 8L, 9L))) {
    figure7_stop("Figure 7 panel 7B requires comparison IDs 1, 8, and 9")
  }
  config
}

figure7_tgi_day <- function(config) as.integer(config$tgi$day)

figure7_tgi_measure <- function(config) {
  paste0("TGI_percent_Day_", figure7_tgi_day(config))
}

figure7_tgi_delta_measure <- function(config) {
  paste0("tumor_volume_delta_Day_", figure7_tgi_day(config))
}

figure7_tgi_label <- function(config) {
  paste("Day", figure7_tgi_day(config))
}

figure7_add_tgi_metadata <- function(data, config) {
  data$tgi_outcome <- "day"
  data$tgi_day <- figure7_tgi_day(config)
  data$tgi_measure <- figure7_tgi_measure(config)
  data$matched_control_summary <- "mean"
  data$matched_control_group <- "initial_ploidy"
  data
}

figure7_verify_checksum <- function(path, expected, label = basename(path)) {
  if (!file.exists(path)) figure7_stop("Missing ", label, ": ", path)
  expected <- as.character(expected)
  if (!grepl("^[0-9a-f]{64}$", expected)) {
    figure7_stop("Canonical SHA-256 has not been frozen for ", label, ": ", expected)
  }
  observed <- figure7_sha256(path)
  if (!identical(observed, expected)) {
    figure7_stop("SHA-256 mismatch for ", label, ": expected ", expected, "; observed ", observed)
  }
  invisible(observed)
}

figure7_verify_named_checksums <- function(paths, expected) {
  if (is.null(names(paths)) || is.null(names(expected)) ||
      !identical(sort(names(paths)), sort(names(expected)))) {
    figure7_stop("Named checksum paths and expectations must contain identical keys")
  }
  observed <- vapply(names(expected), function(key) {
    figure7_verify_checksum(paths[[key]], expected[[key]], key)
  }, character(1L))
  invisible(observed)
}

figure7_assert_empty_output <- function(path) {
  if (dir.exists(path)) {
    existing <- list.files(path, recursive = TRUE, all.files = TRUE, no.. = TRUE,
                           include.dirs = FALSE)
    allowed_manager_skeleton <- c("logs/stdout.log", "logs/stderr.log")
    unexpected <- setdiff(existing, allowed_manager_skeleton)
    if (length(unexpected)) {
      figure7_stop("Output directory contains pre-existing analysis files: ", paste(unexpected, collapse = ", "))
    }
  }
  invisible(TRUE)
}

figure7_prepare_output <- function(path) {
  figure7_assert_empty_output(path)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  for (subdir in c("figures", "tables", "metadata", "logs")) {
    dir.create(file.path(path, subdir), recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

figure7_read_tsv <- function(path, required = character(), label = basename(path)) {
  if (!file.exists(path)) figure7_stop("Missing ", label, ": ", path)
  data <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  character_columns <- vapply(data, is.character, logical(1L))
  data[character_columns] <- lapply(data[character_columns], function(x) {
    x <- gsub("\\n", "\n", x, fixed = TRUE)
    x <- gsub("\\r", "\r", x, fixed = TRUE)
    gsub("\\t", "\t", x, fixed = TRUE)
  })
  missing <- setdiff(required, names(data))
  if (length(missing)) figure7_stop(label, " is missing column(s): ", paste(missing, collapse = ", "))
  data
}

figure7_write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  encoded <- as.data.frame(data, stringsAsFactors = FALSE)
  character_columns <- vapply(encoded, is.character, logical(1L))
  encoded[character_columns] <- lapply(encoded[character_columns], function(x) {
    x <- gsub("\t", "\\t", x, fixed = TRUE)
    x <- gsub("\r", "\\r", x, fixed = TRUE)
    gsub("\n", "\\n", x, fixed = TRUE)
  })
  utils::write.table(encoded, path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "NA")
  invisible(path)
}

figure7_save_panel <- function(plot, pdf_path, width, height, png_dpi = 300) {
  png_path <- sub("[.]pdf$", ".png", pdf_path, ignore.case = TRUE)
  if (identical(png_path, pdf_path)) figure7_stop("Figure 7 panel path must end in .pdf: ", pdf_path)
  ggplot2::ggsave(pdf_path, plot, device = grDevices::cairo_pdf,
                  width = width, height = height, units = "in")
  ggplot2::ggsave(png_path, plot, device = "png", dpi = png_dpi,
                  width = width, height = height, units = "in", bg = "white")
  paths <- c(pdf_path, png_path)
  if (any(!file.exists(paths)) || any(file.info(paths)$size <= 0)) {
    figure7_stop("Failed to write PDF/PNG panel pair: ", pdf_path)
  }
  invisible(paths)
}

figure7_panel_ids <- function(include_panel_f = TRUE) {
  if (isTRUE(include_panel_f)) c("7A", "7B", "7C", "7D", "7E", "7F") else c("7A", "7B", "7C", "7D", "7E")
}

figure7_panel_filenames <- function(config, panel_ids = figure7_panel_ids(TRUE)) {
  filenames <- unname(unlist(config$panels$filenames[panel_ids], use.names = FALSE))
  sub("day17", paste0("day", figure7_tgi_day(config)), filenames, ignore.case = TRUE)
}

figure7_panel_asset_filenames <- function(config, panel_ids = figure7_panel_ids(TRUE)) {
  pdfs <- figure7_panel_filenames(config, panel_ids)
  c(pdfs, sub("[.]pdf$", ".png", pdfs, ignore.case = TRUE))
}

figure7_validate_figure_inventory <- function(output_dir, config, panel_ids = figure7_panel_ids(TRUE)) {
  figures_dir <- file.path(output_dir, "figures")
  all_figures <- list.files(
    output_dir, pattern = "[.](pdf|png|svg|tiff?|jpg|jpeg)$",
    recursive = TRUE, full.names = TRUE, ignore.case = TRUE
  )
  observed <- sort(basename(all_figures))
  expected <- sort(figure7_panel_asset_filenames(config, panel_ids))
  if (!identical(observed, expected) ||
      any(dirname(normalizePath(all_figures)) != normalizePath(figures_dir))) {
    figure7_stop("Figure inventory mismatch. Expected: ", paste(expected, collapse = ", "),
                 "; observed: ", paste(observed, collapse = ", "))
  }
  if (any(file.info(all_figures)$size <= 0)) figure7_stop("Every Figure 7 asset must be nonempty")
  invisible(TRUE)
}

figure7_copy_file <- function(from, to) {
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(from, to, overwrite = FALSE, copy.mode = TRUE)) figure7_stop("Could not copy ", from, " to ", to)
  invisible(to)
}
