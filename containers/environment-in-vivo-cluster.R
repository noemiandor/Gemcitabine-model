#!/usr/bin/env Rscript

read_versions <- function(path) {
  versions <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("ecosystem", "package", "version")
  if (!identical(names(versions), required)) {
    stop("Version manifest must contain exactly: ecosystem, package, version", call. = FALSE)
  }
  if (any(!nzchar(versions$ecosystem)) || any(!nzchar(versions$package)) ||
      any(!nzchar(versions$version)) || anyDuplicated(versions[c("ecosystem", "package")])) {
    stop("Version manifest contains blank or duplicate records", call. = FALSE)
  }
  versions
}

expected_value <- function(versions, ecosystem, package) {
  value <- versions$version[versions$ecosystem == ecosystem & versions$package == package]
  if (length(value) != 1L) {
    stop("Missing unique version record for ", ecosystem, "/", package, call. = FALSE)
  }
  value
}

assert_platform <- function(versions) {
  architecture <- trimws(system2("uname", "-m", stdout = TRUE))
  if (!identical(architecture, "x86_64")) {
    stop("This image supports linux/amd64 only; found architecture: ", architecture, call. = FALSE)
  }
  expected_r <- expected_value(versions, "system", "R")
  if (getRversion() != package_version(expected_r)) {
    stop("R version mismatch: expected ", expected_r, ", found ", getRversion(), call. = FALSE)
  }
  architecture
}

validate_direct_packages <- function(versions) {
  locked_versions <- versions[
    versions$ecosystem %in% c("r", "recommended", "bioc"),
    c("package", "version")
  ]
  direct_packages <- c(
    "BiocManager", "pak", "Seurat", "scDblFinder", "SingleCellExperiment",
    "sctransform", "hdf5r", "readxl", "xml2", "dplyr", "tidyr", "readr",
    "tibble", "ggplot2", "pheatmap", "patchwork", "EnhancedVolcano", "RANN",
    "Matrix", "msigdbr", "fgsea", "yaml", "future", "scales", "rlang"
  )
  if (!all(direct_packages %in% locked_versions$package)) {
    stop("R package manifest is missing an audited in-vivo-cluster direct dependency", call. = FALSE)
  }
  list(locked_versions = locked_versions, direct_packages = direct_packages)
}

prepare_install <- function(versions_path) {
  versions <- read_versions(versions_path)
  assert_platform(versions)
  options(timeout = 3600, Ncpus = 4L)
  validation <- validate_direct_packages(versions)
  c(list(versions = versions), validation)
}

install_tools <- function(versions_path) {
  environment <- prepare_install(versions_path)
  versions <- environment$versions
  locked_versions <- environment$locked_versions
  bioc_manager_version <- locked_versions$version[locked_versions$package == "BiocManager"]
  pak_version <- locked_versions$version[locked_versions$package == "pak"]
  install.packages(
    sprintf("https://cran.r-project.org/src/contrib/Archive/BiocManager/BiocManager_%s.tar.gz", bioc_manager_version),
    repos = NULL,
    type = "source"
  )
  BiocManager::install(
    version = expected_value(versions, "system", "Bioconductor"),
    ask = FALSE,
    update = FALSE
  )
  install.packages(
    sprintf("https://cran.r-project.org/src/contrib/Archive/pak/pak_%s.tar.gz", pak_version),
    repos = NULL,
    type = "source"
  )
}

install_cran_bootstrap <- function(versions_path) {
  environment <- prepare_install(versions_path)
  versions <- environment$versions
  direct_packages <- environment$direct_packages
  bioc_versions <- versions[versions$ecosystem == "bioc", c("package", "version")]
  bootstrap_packages <- setdiff(
    direct_packages,
    c("BiocManager", "pak", "Matrix", bioc_versions$package)
  )
  pak::pkg_install(
    bootstrap_packages,
    dependencies = c("Depends", "Imports", "LinkingTo"),
    upgrade = FALSE,
    ask = FALSE
  )
}

install_bioc_packages <- function(versions_path) {
  environment <- prepare_install(versions_path)
  versions <- environment$versions
  bioc_versions <- versions[versions$ecosystem == "bioc", c("package", "version")]
  BiocManager::install(
    bioc_versions$package,
    version = expected_value(versions, "system", "Bioconductor"),
    ask = FALSE,
    update = FALSE,
    type = "source"
  )
}

restore_cran_versions <- function(versions_path) {
  environment <- prepare_install(versions_path)
  versions <- environment$versions
  # The local environment accumulated packages over time, so a single current
  # CRAN solve cannot recreate it. Restore the audited dependency closure in
  # dependency-first manifest order without asking the solver to change it.
  cran_versions <- versions[versions$ecosystem == "r", c("package", "version")]
  for (i in seq_len(nrow(cran_versions))) {
    package <- cran_versions$package[[i]]
    expected <- cran_versions$version[[i]]
    actual <- if (requireNamespace(package, quietly = TRUE)) {
      as.character(packageVersion(package))
    } else {
      NA_character_
    }
    if (!is.na(actual) && package_version(actual) == package_version(expected)) next
    message("Restoring ", package, " ", expected)
    pak::pkg_install(
      sprintf("%s@%s", package, expected),
      dependencies = FALSE,
      upgrade = FALSE,
      ask = FALSE
    )
  }
}

install_environment <- function(versions_path) {
  install_tools(versions_path)
  install_cran_bootstrap(versions_path)
  install_bioc_packages(versions_path)
  restore_cran_versions(versions_path)
}

verify_environment <- function(versions_path) {
  versions <- read_versions(versions_path)
  architecture <- assert_platform(versions)
  r_versions <- versions[
    versions$ecosystem %in% c("r", "recommended", "bioc"),
    c("package", "version")
  ]

  installed <- vapply(
    r_versions$package,
    function(package) {
      if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
      as.character(packageVersion(package))
    },
    character(1)
  )
  mismatch <- is.na(installed) | !mapply(
    function(actual, expected) {
      if (is.na(actual)) return(FALSE)
      package_version(actual) == package_version(expected)
    },
    installed,
    r_versions$version,
    USE.NAMES = FALSE
  )
  if (any(mismatch)) {
    print(
      data.frame(
        package = r_versions$package[mismatch],
        expected = r_versions$version[mismatch],
        actual = unname(installed[mismatch]),
        row.names = NULL
      ),
      row.names = FALSE
    )
    stop(sum(mismatch), " required R package versions do not match", call. = FALSE)
  }

  expected_bioc <- expected_value(versions, "system", "Bioconductor")
  actual_bioc <- as.character(BiocManager::version())
  if (package_version(actual_bioc) != package_version(expected_bioc)) {
    stop("Bioconductor version mismatch: expected ", expected_bioc, ", found ", actual_bioc, call. = FALSE)
  }

  python_code <- paste(
    "import platform, sys, venv",
    "expected = sys.argv[1]",
    "actual = '.'.join(map(str, sys.version_info[:2]))",
    "assert actual == expected, f'Python version mismatch: expected {expected}, found {platform.python_version()}'",
    "print(f'Python={platform.python_version()}')",
    sep = "\n"
  )
  python_output <- system2(
    "python3",
    c("-c", shQuote(python_code), shQuote(expected_value(versions, "system", "Python"))),
    stdout = TRUE,
    stderr = TRUE
  )
  python_status <- attr(python_output, "status")
  cat(python_output, sep = "\n")
  cat("\n")
  if (!is.null(python_status) && python_status != 0L) {
    stop("Python verification failed", call. = FALSE)
  }

  counts <- Matrix::sparseMatrix(
    i = c(1L, 2L), j = c(1L, 2L), x = c(1, 2),
    dims = c(2L, 2L),
    dimnames = list(c("gene1", "gene2"), c("cell1", "cell2"))
  )
  object <- Seurat::CreateSeuratObject(counts = counts)
  if (!inherits(object, "Seurat")) stop("Seurat smoke test failed", call. = FALSE)
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("SingleCellExperiment smoke test failed", call. = FALSE)
  }
  if (!capabilities("tiff")) stop("R TIFF output capability is unavailable", call. = FALSE)

  cat("R=", as.character(getRversion()), "\n", sep = "")
  cat("Bioconductor=", actual_bioc, "\n", sep = "")
  cat("locked_R_packages=", nrow(r_versions), "\n", sep = "")
  cat("architecture=", architecture, "\n", sep = "")
  cat("seurat_smoke=PASS\n")
  cat("single_cell_experiment_smoke=PASS\n")
  cat("tiff_capability=PASS\n")
  cat("environment_verification=in-vivo-cluster_PASS\n")
}

actions <- c(
  "install", "install-tools", "install-cran-bootstrap", "install-bioc",
  "restore-cran", "verify"
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !args[[1L]] %in% actions) {
  stop(
    paste0("usage: environment-in-vivo-cluster.R [", paste(actions, collapse = "|"), "] PACKAGES_TSV"),
    call. = FALSE
  )
}

switch(
  args[[1L]],
  install = install_environment(args[[2L]]),
  `install-tools` = install_tools(args[[2L]]),
  `install-cran-bootstrap` = install_cran_bootstrap(args[[2L]]),
  `install-bioc` = install_bioc_packages(args[[2L]]),
  `restore-cran` = restore_cran_versions(args[[2L]]),
  verify = verify_environment(args[[2L]])
)
