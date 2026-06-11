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

write_session_metadata <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  capture.output(sessionInfo(), file = path)
}

write_input_manifest <- function(files, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  existing <- files[file.exists(files)]
  manifest <- data.frame(
    path = existing,
    md5 = unname(tools::md5sum(existing)),
    stringsAsFactors = FALSE
  )
  write.table(manifest, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}
