setup_local_lib <- function(base_dir) {
  local_lib <- file.path(base_dir, ".Rlibs")
  if (dir.exists(local_lib)) {
    .libPaths(c(local_lib, .libPaths()))
  }
  invisible(.libPaths())
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required R packages: ",
      paste(missing, collapse = ", "),
      ". Install them before running the analysis.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

write_metadata_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

write_run_metadata <- function(parameters, metadata_dir) {
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  if (!all(c("key", "value") %in% colnames(parameters))) {
    stop("Run metadata must contain 'key' and 'value' columns.", call. = FALSE)
  }
  parameters <- parameters[, c("key", "value"), drop = FALSE]
  parameters$value <- as.character(parameters$value)
  write_metadata_tsv(parameters, file.path(metadata_dir, "run_config.tsv"))

  run_parameters <- parameters
  colnames(run_parameters) <- c("parameter", "value")
  write_metadata_tsv(run_parameters, file.path(metadata_dir, "run_parameters.tsv"))
}

write_session_metadata <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  capture.output(sessionInfo(), file = path)
}

write_input_manifest <- function(files, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  files <- unique(normalizePath(files[file.exists(files)], mustWork = TRUE))
  manifest <- data.frame(
    path = files,
    md5 = unname(tools::md5sum(files)),
    stringsAsFactors = FALSE
  )
  write.table(manifest, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}
