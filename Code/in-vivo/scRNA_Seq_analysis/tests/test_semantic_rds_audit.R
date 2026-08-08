#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_arg <- sub(
  "^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep(
    "^--file=", commandArgs(trailingOnly = FALSE)
  )]
)
script_dir <- normalizePath(file.path(dirname(script_arg), ".."), mustWork = TRUE)
source(file.path(script_dir, "audit_generated_seurat_rds.R"), local = FALSE)

for (package in c("Seurat", "Matrix")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Audit regression test requires package '", package, "'.", call. = FALSE)
  }
}
suppressPackageStartupMessages(library(Seurat))

make_fixture <- function(n_cells = 300L) {
  set.seed(5826)
  cells <- sprintf("cell_%03d", seq_len(n_cells))
  genes <- sprintf("gene_%02d", seq_len(12L))
  counts <- Matrix::Matrix(
    matrix(
      stats::rpois(length(cells) * length(genes), lambda = 3),
      nrow = length(genes),
      dimnames = list(genes, cells)
    ),
    sparse = TRUE
  )
  object <- Seurat::CreateSeuratObject(counts, project = "semantic-audit-fixture")
  object$sample <- factor(rep(c("sample_a", "sample_b"), length.out = n_cells))
  object$clusters <- factor(rep(c("0", "1", "2"), length.out = n_cells))
  object$seurat_clusters <- object$clusters
  Seurat::Idents(object) <- object$clusters
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- Seurat::FindVariableFeatures(
    object,
    selection.method = "vst",
    nfeatures = 8L,
    verbose = FALSE
  )
  object <- Seurat::ScaleData(object, features = genes, verbose = FALSE)

  embeddings <- matrix(
    stats::rnorm(n_cells * 3L),
    nrow = n_cells,
    dimnames = list(cells, paste0("PC_", 1:3))
  )
  loadings <- matrix(
    stats::rnorm(length(genes) * 3L),
    nrow = length(genes),
    dimnames = list(genes, paste0("PC_", 1:3))
  )
  object[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = embeddings,
    loadings = loadings,
    stdev = c(3, 2, 1),
    key = "PC_",
    assay = "RNA"
  )
  umap <- cbind(
    UMAP_1 = seq(-3, 3, length.out = n_cells),
    UMAP_2 = sin(seq(0, 8 * pi, length.out = n_cells))
  )
  rownames(umap) <- cells
  object[["umap"]] <- Seurat::CreateDimReducObject(
    embeddings = umap,
    key = "UMAP_",
    assay = "RNA"
  )
  graph_matrix <- Matrix::bandSparse(
    n_cells,
    n_cells,
    k = c(-1L, 0L, 1L),
    diagonals = list(rep(0.25, n_cells - 1L), rep(1, n_cells), rep(0.25, n_cells - 1L))
  )
  dimnames(graph_matrix) <- list(cells, cells)
  object[["RNA_snn"]] <- SeuratObject::as.Graph(graph_matrix)
  object
}

clone_object <- function(object) unserialize(serialize(object, NULL))

run_case <- function(name, generated, reference, expected_status, failed_check = NULL) {
  root <- tempfile(paste0("semantic_audit_", name, "_"))
  dir.create(root)
  generated_path <- file.path(root, "generated.rds")
  reference_path <- file.path(root, "reference.rds")
  saveRDS(generated, generated_path)
  saveRDS(reference, reference_path)
  result <- audit_generated_seurat_rds(
    generated_path,
    reference_path,
    file.path(root, "audit"),
    stop_on_failure = FALSE
  )
  if (!identical(result$status, expected_status)) {
    stop(
      name, ": expected ", expected_status, " but observed ", result$status,
      call. = FALSE
    )
  }
  if (!is.null(failed_check) && !any(
    grepl(failed_check, result$summary$check) & result$summary$status == "FAIL"
  )) {
    stop(name, ": expected failing check matching ", failed_check, call. = FALSE)
  }
  expected_files <- c(
    "audit_summary.tsv", "metadata_comparison.tsv", "cluster_comparison.tsv",
    "graph_comparison.tsv", "assay_numeric_comparison.tsv", "pca_comparison.tsv",
    "umap_displacement_summary.tsv", "umap_largest_displacements.tsv",
    "command_comparison.tsv", "rds_file_identity.tsv", "AUDIT_COMPLETE.txt"
  )
  if (!all(file.exists(file.path(root, "audit", expected_files)))) {
    stop(name, ": audit output set is incomplete", call. = FALSE)
  }
  invisible(result)
}

reference <- make_fixture()
generated <- clone_object(reference)
generated@commands[[1L]]@time.stamp <-
  generated@commands[[1L]]@time.stamp + 60
generated@assays$RNA@data[1L, 1L] <-
  generated@assays$RNA@data[1L, 1L] + 5e-8
generated@reductions$pca@cell.embeddings[1L, 1L] <-
  generated@reductions$pca@cell.embeddings[1L, 1L] + 5e-8
generated@reductions$umap@cell.embeddings <-
  generated@reductions$umap@cell.embeddings + 0.01
pass_result <- run_case("tolerated_numeric_and_timestamp", generated, reference, "PASS")
if (identical(
      pass_result$identities$md5[[1L]],
      pass_result$identities$md5[[2L]]
    )) {
  stop("PASS fixture must demonstrate that byte-identical RDS MD5 is not required")
}

mutated <- clone_object(reference)
mutated$sample[[1L]] <- "sample_b"
run_case("metadata", mutated, reference, "FAIL", "column:sample")

mutated <- clone_object(reference)
levels(mutated$clusters) <- c(levels(mutated$clusters), "unused")
run_case("cluster_levels", mutated, reference, "FAIL", "metadata_cluster:clusters")

mutated <- clone_object(reference)
mutated@graphs$RNA_snn@x[[1L]] <- mutated@graphs$RNA_snn@x[[1L]] + 0.01
run_case("graph", mutated, reference, "FAIL", "graph:RNA_snn")

mutated <- clone_object(reference)
mutated@assays$RNA@counts@x[[1L]] <- mutated@assays$RNA@counts@x[[1L]] + 1
run_case("counts", mutated, reference, "FAIL", "assay:RNA:counts")

mutated <- clone_object(reference)
mutated@assays$RNA@data[1L, 1L] <- mutated@assays$RNA@data[1L, 1L] + 2e-6
run_case("assay_tolerance", mutated, reference, "FAIL", "assay:RNA:data")

mutated <- clone_object(reference)
mutated@reductions$pca@cell.embeddings[1L, 1L] <-
  mutated@reductions$pca@cell.embeddings[1L, 1L] + 2e-6
run_case("pca_tolerance", mutated, reference, "FAIL", "pca:cell.embeddings")

mutated <- clone_object(reference)
mutated@reductions$umap@cell.embeddings[, 1L] <-
  mutated@reductions$umap@cell.embeddings[, 1L] + 0.2
run_case("umap_median", mutated, reference, "FAIL", "umap:median")

mutated <- clone_object(reference)
mutated@reductions$umap@cell.embeddings[1:2, 1L] <-
  mutated@reductions$umap@cell.embeddings[1:2, 1L] + 1.1
run_case("umap_outlier_fraction", mutated, reference, "FAIL", "fraction_cell")

mutated <- clone_object(reference)
command_name <- names(mutated@commands)[[1L]]
parameter_name <- names(mutated@commands[[command_name]]@params)[[1L]]
mutated@commands[[command_name]]@params[[parameter_name]] <- "changed"
run_case("command_parameter", mutated, reference, "FAIL", "command:")

message("Semantic RDS audit regression tests passed")
