#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: run_compatible_umap.R PCA_RDS OUTPUT_RDS SEED", call. = FALSE)
}

private_library <- Sys.getenv("CLUSTER_STANDALONE_R_LIBRARY", unset = "")
if (nzchar(private_library) && dir.exists(private_library)) {
  .libPaths(c(private_library, .libPaths()))
}
if (!requireNamespace("uwot", quietly = TRUE)) {
  stop("The vendored uwot package is unavailable.", call. = FALSE)
}
if (!identical(as.character(utils::packageVersion("uwot")), "0.2.3")) {
  stop("Compatible UMAP requires uwot 0.2.3.", call. = FALSE)
}

contract_path <- system.file("cluster_standalone_contract.txt", package = "uwot")
required_contract <- c(
  "uwot=0.2.3",
  "gradient_powf=macos_arm64_libsystem_m_powf_compat",
  "powf_reference_samples=200000",
  "powf_reference_bitwise_matches=200000",
  "powf_reference_max_ulp=0",
  "spectral_backend=irlba_with_reference_blas_lapack"
)
if (!nzchar(contract_path) ||
    !identical(readLines(contract_path, warn = FALSE), required_contract)) {
  stop("Installed uwot does not satisfy the cluster compatibility contract.", call. = FALSE)
}

pca <- readRDS(normalizePath(args[[1]], mustWork = TRUE))
if (!is.matrix(pca) || !is.numeric(pca) || ncol(pca) != 30L ||
    is.null(rownames(pca))) {
  stop("PCA input must be a named numeric matrix with exactly 30 columns.", call. = FALSE)
}
seed <- as.integer(args[[3]])
if (length(seed) != 1L || is.na(seed)) stop("SEED must be one integer.", call. = FALSE)

# The reference macOS session did not have RSpectra available, so uwot used
# irlba. Force that same branch even though RSpectra is present in the SIF.
uwot_namespace <- asNamespace("uwot")
original_rspectra_check <- get("rspectra_is_installed", envir = uwot_namespace)
unlockBinding("rspectra_is_installed", uwot_namespace)
assign("rspectra_is_installed", function() FALSE, envir = uwot_namespace)
lockBinding("rspectra_is_installed", uwot_namespace)
on.exit({
  unlockBinding("rspectra_is_installed", uwot_namespace)
  assign("rspectra_is_installed", original_rspectra_check, envir = uwot_namespace)
  lockBinding("rspectra_is_installed", uwot_namespace)
}, add = TRUE)

set.seed(seed)
result <- uwot::umap(
  X = pca,
  n_neighbors = 30L,
  n_components = 2L,
  metric = "cosine",
  n_epochs = NULL,
  learning_rate = 1,
  min_dist = 0.3,
  spread = 1,
  set_op_mix_ratio = 1,
  local_connectivity = 1L,
  repulsion_strength = 1,
  negative_sample_rate = 5L,
  fast_sgd = FALSE,
  n_threads = 1L,
  ret_extra = "fgraph",
  verbose = TRUE
)
embedding <- result$embedding
rownames(embedding) <- rownames(pca)
colnames(embedding) <- c("UMAP_1", "UMAP_2")

session <- utils::sessionInfo()
saveRDS(
  list(
    embedding = embedding,
    seed = seed,
    uwot_version = as.character(utils::packageVersion("uwot")),
    uwot_built = unname(read.dcf(
      file.path(find.package("uwot"), "DESCRIPTION"), fields = "Built"
    )[[1]]),
    blas = session$BLAS,
    lapack = session$LAPACK,
    ld_preload = Sys.getenv("LD_PRELOAD", unset = ""),
    force_irlba = TRUE
  ),
  normalizePath(args[[2]], mustWork = FALSE),
  compress = "xz"
)
