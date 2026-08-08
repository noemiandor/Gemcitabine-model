#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L, 3L)) {
  stop("usage: restore-in-vivo-cluster.R PACKAGES_TSV [START END]", call. = FALSE)
}

versions <- read.delim(args[[1L]], stringsAsFactors = FALSE, check.names = FALSE)
required <- c("ecosystem", "package", "version")
if (!identical(names(versions), required)) {
  stop("Version manifest must contain exactly: ecosystem, package, version", call. = FALSE)
}

architecture <- trimws(system2("uname", "-m", stdout = TRUE))
if (!identical(architecture, "x86_64")) {
  stop("This image supports linux/amd64 only; found architecture: ", architecture, call. = FALSE)
}

options(timeout = 3600, Ncpus = 4L)
cran_versions <- versions[versions$ecosystem == "r", c("package", "version")]
if (length(args) == 3L) {
  start <- as.integer(args[[2L]])
  end <- min(as.integer(args[[3L]]), nrow(cran_versions))
  if (is.na(start) || is.na(end) || start < 1L || end < start) {
    stop("START and END must define a valid positive row range", call. = FALSE)
  }
  cran_versions <- cran_versions[seq.int(start, end), , drop = FALSE]
}

installed_version <- function(package) {
  tryCatch(
    as.character(packageVersion(package)),
    error = function(error) NA_character_
  )
}

for (i in seq_len(nrow(cran_versions))) {
  package <- cran_versions$package[[i]]
  expected <- cran_versions$version[[i]]
  actual <- installed_version(package)
  if (!is.na(actual) && package_version(actual) == package_version(expected)) next

  message("Restoring ", package, " ", expected)
  archive <- tempfile(fileext = ".tar.gz")
  urls <- c(
    sprintf("https://cran.r-project.org/src/contrib/%s_%s.tar.gz", package, expected),
    sprintf("https://cran.r-project.org/src/contrib/Archive/%s/%s_%s.tar.gz", package, package, expected)
  )
  downloaded <- FALSE
  for (url in urls) {
    status <- suppressWarnings(
      tryCatch(
        download.file(url, archive, mode = "wb", quiet = TRUE),
        error = function(error) 1L
      )
    )
    if (identical(status, 0L) && file.exists(archive) && file.size(archive) > 0L) {
      downloaded <- TRUE
      break
    }
  }
  if (!downloaded) {
    stop("Unable to download CRAN source for ", package, " ", expected, call. = FALSE)
  }

  install.packages(archive, repos = NULL, dependencies = FALSE, type = "source")
  unlink(archive)
  actual <- installed_version(package)
  if (is.na(actual) || package_version(actual) != package_version(expected)) {
    stop(
      "Failed to restore ", package, ": expected ", expected, ", found ", actual,
      call. = FALSE
    )
  }
}
