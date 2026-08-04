#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: bootstrap_dependencies.R PRIVATE_R_LIBRARY VENDOR_DIR", call. = FALSE)
}

private_library <- normalizePath(args[[1]], mustWork = FALSE)
vendor_dir <- normalizePath(args[[2]], mustWork = TRUE)
if (!dir.exists(private_library) &&
    !dir.create(private_library, recursive = TRUE, showWarnings = FALSE)) {
  stop("Cannot create private R library: ", private_library, call. = FALSE)
}
private_library <- normalizePath(private_library, mustWork = TRUE)
.libPaths(c(private_library, .libPaths()))

sources <- data.frame(
  package = c(
    "RcppAnnoy", "irlba", "uwot", "sctransform",
    "assorthead", "xgboost", "BiocNeighbors"
  ),
  version = c(
    "0.0.22", "2.3.5.1", "0.2.3", "0.4.2",
    "1.2.0", "1.7.11.1", "2.2.0"
  ),
  archive = file.path(
    vendor_dir,
    c(
      "RcppAnnoy_0.0.22.tar.gz",
      "irlba_2.3.5.1.tar.gz",
      "uwot_0.2.3.tar.gz",
      "sctransform_0.4.2.tar.gz",
      "assorthead_1.2.0.tar.gz",
      "xgboost_1.7.11.1.tar.gz",
      "BiocNeighbors_2.2.0.tar.gz"
    )
  ),
  md5 = c(
    "8ae334c4634f6fdbc066f6a0e93dcb95",
    "f738200d5272c7258ee59f7074dd4a6a",
    "a66c56f91d39b984e0bf02399441004b",
    "02e80f35d5be2993cbc1740c628cc05a",
    "07ba54164ad32cde38585327e94cafc0",
    "fcb32ea43b53faf9625515a4f5197750",
    "a23b27aadce180b89bd71792608006f4"
  ),
  stringsAsFactors = FALSE
)

missing_archives <- sources$archive[!file.exists(sources$archive)]
if (length(missing_archives) > 0L) {
  stop("Vendored source archive(s) missing: ", paste(missing_archives, collapse = ", "), call. = FALSE)
}
observed_md5 <- unname(tools::md5sum(sources$archive))
if (!identical(observed_md5, sources$md5)) {
  bad <- sources$archive[observed_md5 != sources$md5]
  stop("Vendored source checksum mismatch: ", paste(bad, collapse = ", "), call. = FALSE)
}

private_version <- function(package) {
  description_file <- file.path(private_library, package, "DESCRIPTION")
  if (!file.exists(description_file)) return(NA_character_)
  as.character(read.dcf(description_file, fields = "Version")[[1]])
}

install_source <- function(package, version, archive) {
  message("Installing pinned ", package, " ", version, " into ", private_library)
  utils::install.packages(
    archive,
    repos = NULL,
    type = "source",
    lib = private_library,
    INSTALL_opts = c("--preclean", "--no-multiarch")
  )
  observed <- private_version(package)
  if (!identical(observed, version)) {
    stop(
      "Pinned package installation failed for ", package,
      ": expected ", version, ", observed ", observed,
      call. = FALSE
    )
  }
}

# xgboost must be compiled with the container defaults. The custom floating-
# point flags below are only part of the patched BiocNeighbors/Annoy contract.
Sys.unsetenv("R_MAKEVARS_USER")

for (package in c("RcppAnnoy", "irlba")) {
  expected_version <- sources$version[sources$package == package]
  if (!identical(private_version(package), expected_version)) {
    install_source(
      package,
      expected_version,
      sources$archive[sources$package == package]
    )
  }
}

uwot_source_md5 <- sources$md5[sources$package == "uwot"]
uwot_contract <- c(
  "uwot=0.2.3",
  "gradient_powf=macos_arm64_libsystem_m_powf_compat",
  "powf_reference_samples=200000",
  "powf_reference_bitwise_matches=200000",
  "powf_reference_max_ulp=0",
  "spectral_backend=irlba_with_reference_blas_lapack"
)
uwot_installed_contract <- file.path(
  private_library, "uwot", "cluster_standalone_contract.txt"
)
uwot_header <- file.path(private_library, "uwot", "include", "uwot", "apple_powf.h")
uwot_header_md5 <- "69f469487fb4d13ada4caa9c5f9473f1"
uwot_marker_contract <- c(
  uwot_contract,
  paste0("source_md5=", uwot_source_md5),
  paste0("apple_powf_header_md5=", uwot_header_md5)
)
uwot_marker <- file.path(private_library, "uwot.cluster_standalone_contract.txt")
uwot_ok <- identical(private_version("uwot"), "0.2.3") &&
  file.exists(uwot_installed_contract) &&
  identical(readLines(uwot_installed_contract, warn = FALSE), uwot_contract) &&
  file.exists(uwot_header) &&
  identical(unname(tools::md5sum(uwot_header)), uwot_header_md5) &&
  file.exists(uwot_marker) &&
  identical(readLines(uwot_marker, warn = FALSE), uwot_marker_contract)
if (!uwot_ok) {
  install_source("uwot", "0.2.3", sources$archive[sources$package == "uwot"])
  writeLines(uwot_marker_contract, uwot_marker)
}
if (!file.exists(uwot_installed_contract) ||
    !identical(readLines(uwot_installed_contract, warn = FALSE), uwot_contract) ||
    !file.exists(uwot_header) ||
    !identical(unname(tools::md5sum(uwot_header)), uwot_header_md5)) {
  stop("Installed uwot compatibility contract mismatch.", call. = FALSE)
}

sctransform_source_md5 <- sources$md5[sources$package == "sctransform"]
sctransform_contract <- paste(
  "sctransform=0.4.2",
  paste0("source_md5=", sctransform_source_md5),
  "row_gmean_math=macos_arm64_libsystem_m_log_exp_compat",
  "log10=macos_arm64_libsystem_m_compat",
  "density_sampling=R_4.5.1_apple_silicon_fma_compat",
  sep = "\n"
)
sctransform_marker <- file.path(
  private_library, "sctransform.cluster_standalone_contract.txt"
)
sctransform_marker_ok <- file.exists(sctransform_marker) &&
  identical(
    readLines(sctransform_marker, warn = FALSE),
    strsplit(sctransform_contract, "\n", fixed = TRUE)[[1]]
  )
if (!identical(private_version("sctransform"), "0.4.2") ||
    !sctransform_marker_ok) {
  install_source(
    "sctransform", "0.4.2",
    sources$archive[sources$package == "sctransform"]
  )
  writeLines(sctransform_contract, sctransform_marker)
}

assort_header <- file.path(
  private_library, "assorthead", "include", "annoy", "annoylib.h"
)
expected_assort_header_md5 <- "3706ec09e55d95fac8306d48f2bd0ba6"
assort_header_ok <- file.exists(assort_header) &&
  identical(unname(tools::md5sum(assort_header)), expected_assort_header_md5)
if (!identical(private_version("assorthead"), "1.2.0") || !assort_header_ok) {
  install_source("assorthead", "1.2.0", sources$archive[sources$package == "assorthead"])
}
if (!identical(unname(tools::md5sum(assort_header)), expected_assort_header_md5)) {
  stop("Installed assorthead compatibility header checksum mismatch.", call. = FALSE)
}

if (!identical(private_version("xgboost"), "1.7.11.1")) {
  install_source("xgboost", "1.7.11.1", sources$archive[sources$package == "xgboost"])
}

machine <- tolower(paste(R.version$arch, Sys.info()[["machine"]], sep = ":"))
makevars_file <- file.path(private_library, "Makevars.cluster_standalone")
if (grepl("x86_64|amd64", machine)) {
  makevars_lines <- "CXX17FLAGS = -g -O2 -fno-tree-vectorize -mfma -ffp-contract=fast"
} else {
  makevars_lines <- "# Native compiler defaults; the compatibility arithmetic is explicit in annoylib.h"
}
writeLines(makevars_lines, makevars_file)
Sys.setenv(R_MAKEVARS_USER = makevars_file)

build_contract <- paste(
  "BiocNeighbors=2.2.0",
  paste0("source_md5=", sources$md5[sources$package == "BiocNeighbors"]),
  paste0("assorthead_header_md5=", expected_assort_header_md5),
  paste0("machine=", machine),
  paste0("makevars=", paste(makevars_lines, collapse = " ")),
  sep = "\n"
)
build_marker <- file.path(private_library, "BiocNeighbors.cluster_standalone_contract.txt")
marker_ok <- file.exists(build_marker) &&
  identical(readLines(build_marker, warn = FALSE), strsplit(build_contract, "\n", fixed = TRUE)[[1]])
if (!identical(private_version("BiocNeighbors"), "2.2.0") || !marker_ok) {
  install_source(
    "BiocNeighbors", "2.2.0",
    sources$archive[sources$package == "BiocNeighbors"]
  )
  writeLines(build_contract, build_marker)
}

expected_versions <- c(
  RcppAnnoy = "0.0.22",
  irlba = "2.3.5.1",
  sctransform = "0.4.2",
  assorthead = "1.2.0",
  xgboost = "1.7.11.1",
  BiocNeighbors = "2.2.0",
  uwot = "0.2.3"
)
observed_versions <- vapply(names(expected_versions), private_version, character(1))
if (!identical(observed_versions, expected_versions)) {
  stop(
    "Private dependency version mismatch: ",
    paste(names(observed_versions), observed_versions, sep = "=", collapse = ", "),
    call. = FALSE
  )
}

for (package in names(expected_versions)) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Pinned package cannot be loaded: ", package, call. = FALSE)
  }
}

provenance_dir <- file.path(dirname(private_library), "00_provenance")
dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- data.frame(
  package = sources$package,
  version = sources$version,
  source_archive = basename(sources$archive),
  source_md5 = sources$md5,
  private_library = private_library,
  stringsAsFactors = FALSE
)
utils::write.table(
  manifest,
  file.path(provenance_dir, "dependency_bootstrap.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(build_contract, file.path(provenance_dir, "biocneighbors_build_contract.txt"))
writeLines(
  sctransform_contract,
  file.path(provenance_dir, "sctransform_build_contract.txt")
)
writeLines(uwot_marker_contract, file.path(provenance_dir, "uwot_build_contract.txt"))

message(
  "Pinned dependency bootstrap complete: ",
  paste(names(observed_versions), observed_versions, sep = "=", collapse = ", ")
)
