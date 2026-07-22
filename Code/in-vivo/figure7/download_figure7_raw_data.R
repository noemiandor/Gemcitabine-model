#!/usr/bin/env Rscript

# Materialize the pinned Figure 7 loom and Seurat inputs from Zenodo. Downloads
# use resumable .part files and are promoted only after size/checksum validation.

parse_cli_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
      stop("Arguments must use --name=value syntax: ", arg, call. = FALSE)
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[gsub("-", "_", pieces[[1L]])]] <- paste(pieces[-1L], collapse = "=")
  }
  out
}

arg_value <- function(args, name, default = NULL, required = FALSE) {
  value <- args[[gsub("-", "_", name)]]
  if (is.null(value) || !nzchar(value)) value <- default
  if (isTRUE(required) && (is.null(value) || !nzchar(value))) {
    stop("Missing required argument --", gsub("_", "-", name), call. = FALSE)
  }
  value
}

parse_boolean <- function(value, argument) {
  normalized <- tolower(trimws(as.character(value)))
  if (!normalized %in% c("true", "false", "1", "0", "yes", "no", "y", "n")) {
    stop("--", argument, " must be true or false", call. = FALSE)
  }
  normalized %in% c("true", "1", "yes", "y")
}

script_location <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else "Code/in-vivo/figure7/download_figure7_raw_data.R"
}

read_download_manifest <- function(path, enforce_published_contract = FALSE) {
  if (!file.exists(path)) stop("Missing Zenodo download manifest: ", path, call. = FALSE)
  manifest <- utils::read.delim(
    path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "",
    colClasses = "character"
  )
  required <- c("role", "filename", "size_bytes", "md5", "sha256", "url")
  missing <- setdiff(required, names(manifest))
  if (length(missing)) stop("Zenodo manifest is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  manifest$size_bytes <- suppressWarnings(as.numeric(manifest$size_bytes))
  if (!nrow(manifest) || any(!manifest$role %in% c("loom", "seurat_rds")) ||
      any(!is.finite(manifest$size_bytes) | manifest$size_bytes <= 0) ||
      any(basename(manifest$filename) != manifest$filename) || anyDuplicated(manifest$filename) ||
      any(!grepl("^[0-9a-f]{32}$", manifest$md5)) ||
      any(!grepl("^(https://|file://)", manifest$url))) {
    stop("Zenodo manifest contains invalid or unsafe rows", call. = FALSE)
  }
  has_sha <- nzchar(manifest$sha256)
  if (any(has_sha & !grepl("^[0-9a-f]{64}$", manifest$sha256))) {
    stop("Zenodo manifest contains an invalid SHA-256", call. = FALSE)
  }
  if (isTRUE(enforce_published_contract) &&
      (sum(manifest$role == "loom") != 18L || sum(manifest$role == "seurat_rds") != 1L ||
       sum(manifest$size_bytes) != 11395115098)) {
    stop("Pinned Zenodo manifest must contain 18 loom files and one Seurat RDS", call. = FALSE)
  }
  manifest
}

manifest_target_path <- function(raw_data_dir, role, filename) {
  if (identical(role, "loom")) file.path(raw_data_dir, "velocyto_loom", filename) else file.path(raw_data_dir, filename)
}

file_checksum_status <- function(path, expected_size, expected_md5, expected_sha256 = "") {
  if (!file.exists(path)) return(list(valid = FALSE, reason = "missing", md5 = "", sha256 = ""))
  observed_size <- unname(file.info(path)$size)
  if (!is.finite(observed_size) || observed_size != expected_size) {
    return(list(valid = FALSE, reason = "size_mismatch", md5 = "", sha256 = ""))
  }
  observed_md5 <- unname(tools::md5sum(path))
  if (!identical(observed_md5, expected_md5)) {
    return(list(valid = FALSE, reason = "md5_mismatch", md5 = observed_md5, sha256 = ""))
  }
  observed_sha256 <- ""
  if (nzchar(expected_sha256)) {
    if (!requireNamespace("digest", quietly = TRUE)) stop("R package 'digest' is required", call. = FALSE)
    observed_sha256 <- unname(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
    if (!identical(observed_sha256, expected_sha256)) {
      return(list(valid = FALSE, reason = "sha256_mismatch", md5 = observed_md5, sha256 = observed_sha256))
    }
  }
  list(valid = TRUE, reason = "verified", md5 = observed_md5, sha256 = observed_sha256)
}

nearest_existing_directory <- function(path) {
  current <- normalizePath(path, mustWork = FALSE)
  while (!dir.exists(current)) {
    parent <- dirname(current)
    if (identical(parent, current)) return("")
    current <- parent
  }
  current
}

available_bytes <- function(path) {
  existing <- nearest_existing_directory(path)
  if (!nzchar(existing) || !nzchar(Sys.which("df"))) return(NA_real_)
  lines <- suppressWarnings(system2("df", c("-Pk", shQuote(existing)), stdout = TRUE, stderr = TRUE))
  if (!length(lines)) return(NA_real_)
  fields <- strsplit(trimws(tail(lines, 1L)), "[[:space:]]+")[[1L]]
  if (length(fields) < 4L) return(NA_real_)
  blocks <- suppressWarnings(as.numeric(fields[[4L]]))
  if (is.finite(blocks)) blocks * 1024 else NA_real_
}

download_resumable <- function(url, part_path, wget_bin, curl_bin, expected_size) {
  dir.create(dirname(part_path), recursive = TRUE, showWarnings = FALSE)
  if (startsWith(url, "https://") && nzchar(wget_bin) && file.exists(wget_bin)) {
    status <- system2(
      wget_bin,
      c(
        "--continue", "--tries=5", "--timeout=60", "--progress=dot:giga",
        "--output-document", shQuote(part_path), shQuote(url)
      )
    )
    if (!identical(status, 0L) || !file.exists(part_path) || file.info(part_path)$size > expected_size) {
      stop("wget failed or returned an invalid byte range for ", url, call. = FALSE)
    }
    return(invisible(part_path))
  }
  if (capabilities("libcurl")) {
    status <- 1L
    last_error <- NULL
    for (attempt in seq_len(5L)) {
      part_size <- if (file.exists(part_path)) unname(file.info(part_path)$size) else 0
      headers <- if (part_size > 0) c(Range = paste0("bytes=", format(part_size, scientific = FALSE), "-")) else NULL
      status <- tryCatch(
        utils::download.file(
          url,
          part_path,
          method = "libcurl",
          quiet = FALSE,
          mode = if (part_size > 0) "ab" else "wb",
          headers = headers
        ),
        error = function(error) {
          last_error <<- conditionMessage(error)
          1L
        }
      )
      if (identical(status, 0L)) break
      if (attempt < 5L) Sys.sleep(attempt)
    }
    if (!identical(status, 0L) || !file.exists(part_path) || file.info(part_path)$size > expected_size) {
      stop(
        "Download failed or server returned an invalid byte range for ", url,
        if (!is.null(last_error)) paste0(": ", last_error) else "",
        call. = FALSE
      )
    }
    return(invisible(part_path))
  }
  status <- system2(
    curl_bin,
    c(
      "--fail", "--location", "--retry", "5", "--retry-delay", "5",
      "--continue-at", "-", "--output", shQuote(part_path), shQuote(url)
    )
  )
  if (!identical(status, 0L)) stop("Download failed for ", url, "; curl exit status ", status, call. = FALSE)
  invisible(part_path)
}

write_download_provenance <- function(raw_data_dir, source_manifest, audit) {
  provenance_dir <- file.path(raw_data_dir, "provenance")
  dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(source_manifest, file.path(provenance_dir, "required_files_manifest.tsv"), overwrite = TRUE)
  audit_path <- file.path(provenance_dir, "downloaded_files_checksums.tsv")
  if (file.exists(audit_path)) {
    previous <- utils::read.delim(audit_path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
    if (all(names(audit) %in% names(previous))) {
      previous <- previous[!previous$filename %in% audit$filename, names(audit), drop = FALSE]
      audit <- rbind(previous, audit)
    }
  }
  utils::write.table(
    audit,
    audit_path,
    sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )
  record_json <- paste0(
    "{\n",
    "  \"doi\": \"10.5281/zenodo.21463392\",\n",
    "  \"record_id\": 21463392,\n",
    "  \"record_url\": \"https://zenodo.org/records/21463392\",\n",
    "  \"accessed_utc\": \"", format(Sys.time(), tz = "UTC", usetz = TRUE), "\"\n",
    "}\n"
  )
  writeLines(record_json, file.path(provenance_dir, "zenodo_record.json"), useBytes = TRUE)
  invisible(provenance_dir)
}

download_figure7_raw_data <- function(raw_data_dir, manifest_path, roles = c("loom", "seurat_rds"), wget_bin = Sys.which("wget"), curl_bin = Sys.which("curl"), enforce_published_contract = FALSE, allow_download = TRUE) {
  roles <- unique(roles)
  if (!length(roles) || any(!roles %in% c("loom", "seurat_rds"))) stop("--roles must select loom, seurat_rds, or all", call. = FALSE)
  has_wget <- nzchar(wget_bin) && file.exists(wget_bin)
  has_curl <- nzchar(curl_bin) && file.exists(curl_bin)
  if (!has_wget && !capabilities("libcurl") && !has_curl) {
    stop("wget, R libcurl support, or a system curl executable is required for resumable Zenodo downloads", call. = FALSE)
  }
  manifest <- read_download_manifest(manifest_path, enforce_published_contract)
  selected <- manifest[manifest$role %in% roles, , drop = FALSE]
  if (!nrow(selected)) stop("No required files selected from Zenodo manifest", call. = FALSE)
  selected$path <- mapply(manifest_target_path, raw_data_dir, selected$role, selected$filename, USE.NAMES = FALSE)

  initial <- lapply(seq_len(nrow(selected)), function(i) {
    file_checksum_status(selected$path[[i]], selected$size_bytes[[i]], selected$md5[[i]], selected$sha256[[i]])
  })
  need <- !vapply(initial, `[[`, logical(1L), "valid")
  if (any(need) && !isTRUE(allow_download)) {
    stop(
      "Required Figure 7 raw files are missing or invalid and downloading is disabled: ",
      paste(selected$filename[need], collapse = ", "),
      call. = FALSE
    )
  }
  required_bytes <- 0
  if (any(need)) {
    for (i in which(need)) {
      part <- paste0(selected$path[[i]], ".part")
      part_size <- if (file.exists(part)) unname(file.info(part)$size) else 0
      if (!is.finite(part_size) || part_size > selected$size_bytes[[i]]) {
        unlink(part, force = TRUE)
        part_size <- 0
      }
      required_bytes <- required_bytes + selected$size_bytes[[i]] - part_size
    }
    free_bytes <- available_bytes(raw_data_dir)
    if (is.finite(free_bytes) && free_bytes < required_bytes * 1.05) {
      stop(
        "Insufficient disk space for Figure 7 raw data: need approximately ",
        format(ceiling(required_bytes * 1.05), scientific = FALSE), " bytes including safety margin; available ",
        format(free_bytes, scientific = FALSE), " bytes",
        call. = FALSE
      )
    }
  }

  audit_rows <- vector("list", nrow(selected))
  for (i in seq_len(nrow(selected))) {
    row <- selected[i, , drop = FALSE]
    target <- row$path[[1L]]
    status <- initial[[i]]
    action <- "reused"
    if (!isTRUE(status$valid)) {
      part <- paste0(target, ".part")
      part_status <- file_checksum_status(part, row$size_bytes[[1L]], row$md5[[1L]], row$sha256[[1L]])
      if (!isTRUE(part_status$valid) && file.exists(part) && file.info(part)$size >= row$size_bytes[[1L]]) {
        unlink(part, force = TRUE)
      }
      if (!isTRUE(part_status$valid)) download_resumable(row$url[[1L]], part, wget_bin, curl_bin, row$size_bytes[[1L]])
      part_status <- file_checksum_status(part, row$size_bytes[[1L]], row$md5[[1L]], row$sha256[[1L]])
      if (!isTRUE(part_status$valid)) {
        stop("Downloaded file failed integrity validation: ", row$filename[[1L]], " (", part_status$reason, ")", call. = FALSE)
      }
      if (file.exists(target) && !unlink(target, force = TRUE)) stop("Could not replace invalid cached file: ", target, call. = FALSE)
      if (!file.rename(part, target)) stop("Could not atomically promote verified download: ", target, call. = FALSE)
      status <- file_checksum_status(target, row$size_bytes[[1L]], row$md5[[1L]], row$sha256[[1L]])
      action <- "downloaded"
    }
    if (!isTRUE(status$valid)) stop("Cached file failed final validation: ", target, call. = FALSE)
    observed_sha256 <- status$sha256
    if (!nzchar(observed_sha256)) {
      if (!requireNamespace("digest", quietly = TRUE)) stop("R package 'digest' is required", call. = FALSE)
      observed_sha256 <- unname(digest::digest(target, algo = "sha256", file = TRUE, serialize = FALSE))
    }
    audit_rows[[i]] <- data.frame(
      role = row$role, filename = row$filename, path = target, action = action,
      size_bytes = row$size_bytes, md5 = status$md5, sha256 = observed_sha256,
      verified_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), stringsAsFactors = FALSE
    )
    message("[", action, "] ", row$filename[[1L]])
  }
  audit <- do.call(rbind, audit_rows)
  write_download_provenance(raw_data_dir, manifest_path, audit)
  list(
    raw_data_dir = normalizePath(raw_data_dir, mustWork = TRUE),
    loom_root = normalizePath(file.path(raw_data_dir, "velocyto_loom"), mustWork = FALSE),
    seurat_rds = if (any(manifest$role == "seurat_rds")) {
      normalizePath(
        manifest_target_path(raw_data_dir, "seurat_rds", manifest$filename[manifest$role == "seurat_rds"][[1L]]),
        mustWork = FALSE
      )
    } else "",
    downloaded = sum(audit$action == "downloaded"), reused = sum(audit$action == "reused"), audit = audit
  )
}

main <- function() {
  args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  script_dir <- dirname(normalizePath(script_location(), mustWork = FALSE))
  repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE)
  raw_data_dir <- arg_value(args, "raw-data-dir", file.path(repo_root, "Data", "in-vivo", "figure7", "raw", "zenodo_21463392"))
  if (!grepl("^/", raw_data_dir)) raw_data_dir <- file.path(repo_root, raw_data_dir)
  manifest_arg <- arg_value(args, "manifest", file.path(script_dir, "zenodo_required_files.tsv"))
  if (!grepl("^/", manifest_arg)) manifest_arg <- file.path(repo_root, manifest_arg)
  roles_arg <- strsplit(arg_value(args, "roles", "all"), ",", fixed = TRUE)[[1L]]
  roles <- if ("all" %in% roles_arg) c("loom", "seurat_rds") else trimws(roles_arg)
  default_manifest <- normalizePath(file.path(script_dir, "zenodo_required_files.tsv"), mustWork = TRUE)
  result <- download_figure7_raw_data(
    raw_data_dir = normalizePath(raw_data_dir, mustWork = FALSE),
    manifest_path = normalizePath(manifest_arg, mustWork = TRUE),
    roles = roles,
    wget_bin = arg_value(args, "wget", Sys.which("wget")),
    curl_bin = arg_value(args, "curl", Sys.which("curl")),
    enforce_published_contract = identical(normalizePath(manifest_arg, mustWork = TRUE), default_manifest),
    allow_download = parse_boolean(arg_value(args, "allow-download", "true"), "allow-download")
  )
  message("Figure 7 raw data ready: ", result$raw_data_dir)
  message("Downloaded: ", result$downloaded, "; reused: ", result$reused)
}

if (identical(environment(), globalenv())) main()
