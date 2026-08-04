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
  stop(
    "Usage: validate_reproduction.R GENERATED_OUTPUT_DIR REFERENCE_REFINED_RDS REFERENCE_FINAL_RDS",
    call. = FALSE
  )
}

generated_root <- normalizePath(args[[1]], mustWork = TRUE)
reference_refined <- normalizePath(args[[2]], mustWork = TRUE)
reference_final <- normalizePath(args[[3]], mustWork = TRUE)
generated_refined <- file.path(
  generated_root,
  "02b_cluster_refine", "objects", "integrated_sct_cca_seurat_cluster_refine.rds"
)
generated_final <- file.path(
  generated_root,
  "03_final_cluster", "03_objects", "integrated_sct_cca_seurat_final_reclustered.rds"
)

for (path in c(generated_refined, generated_final)) {
  if (!file.exists(path)) stop("Generated object is missing: ", path, call. = FALSE)
}
for (package in c("Seurat", "Matrix", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Validation requires package '", package, "'.", call. = FALSE)
  }
}

suppressPackageStartupMessages(library(Seurat))

ensure_dir <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create validation directory: ", path, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

hash_atomic_chunks <- function(x, chunk_size = 5000000L) {
  n <- length(x)
  if (n == 0L) return(digest::digest(list(typeof(x), 0L), algo = "xxhash64"))
  starts <- seq.int(1L, n, by = chunk_size)
  chunk_hashes <- vapply(starts, function(start) {
    end <- min(n, start + chunk_size - 1L)
    digest::digest(x[start:end], algo = "xxhash64", serialize = TRUE)
  }, character(1))
  digest::digest(list(typeof(x), n, chunk_hashes), algo = "xxhash64")
}

hash_character <- function(x) {
  hash_atomic_chunks(enc2utf8(as.character(x)), chunk_size = 1000000L)
}

hash_matrix_exact <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(digest::digest(list(dim(x), "empty"), algo = "xxhash64"))
  }
  dim_hash <- digest::digest(dim(x), algo = "xxhash64")
  dimname_hash <- digest::digest(
    list(hash_character(rownames(x)), hash_character(colnames(x))),
    algo = "xxhash64"
  )
  if (inherits(x, "sparseMatrix")) {
    x <- methods::as(x, "dgCMatrix")
    payload_hash <- digest::digest(
      c(hash_atomic_chunks(x@p), hash_atomic_chunks(x@i), hash_atomic_chunks(x@x)),
      algo = "xxhash64"
    )
  } else {
    row_starts <- seq.int(1L, nrow(x), by = 100L)
    row_hashes <- vapply(row_starts, function(start) {
      end <- min(nrow(x), start + 99L)
      block <- x[start:end, , drop = FALSE]
      digest::digest(block, algo = "xxhash64", serialize = TRUE)
    }, character(1))
    payload_hash <- digest::digest(row_hashes, algo = "xxhash64")
  }
  digest::digest(c(dim_hash, dimname_hash, payload_hash), algo = "xxhash64")
}

assay_signature <- function(obj, assay_name) {
  assay <- obj[[assay_name]]
  slots <- intersect(c("counts", "data", "scale.data"), methods::slotNames(assay))
  slot_hashes <- stats::setNames(vapply(slots, function(slot_name) {
    matrix_value <- methods::slot(assay, slot_name)
    hash_matrix_exact(matrix_value)
  }, character(1)), slots)
  list(
    class = class(assay),
    dimensions = dim(assay),
    feature_names = rownames(assay),
    variable_features = Seurat::VariableFeatures(assay),
    slot_hashes = slot_hashes
  )
}

object_signature <- function(path) {
  message("Reading for validation: ", path)
  obj <- readRDS(path)
  if (!inherits(obj, "Seurat")) stop("Not a Seurat object: ", path, call. = FALSE)
  assay_names <- names(obj@assays)
  reduction_names <- names(obj@reductions)
  reduction_embeddings <- stats::setNames(lapply(reduction_names, function(name) {
    Seurat::Embeddings(obj, name)
  }), reduction_names)
  reduction_loadings_hash <- stats::setNames(vapply(reduction_names, function(name) {
    hash_matrix_exact(obj[[name]]@feature.loadings)
  }, character(1)), reduction_names)
  metadata <- obj@meta.data
  factor_levels <- lapply(metadata, function(x) if (is.factor(x)) levels(x) else NULL)
  signature <- list(
    file = path,
    object_class = class(obj),
    cells = colnames(obj),
    assays = stats::setNames(lapply(assay_names, function(name) assay_signature(obj, name)), assay_names),
    default_assay = Seurat::DefaultAssay(obj),
    metadata = metadata,
    metadata_classes = vapply(metadata, function(x) class(x)[1], character(1)),
    metadata_factor_levels = factor_levels,
    active_ident = Seurat::Idents(obj),
    reductions = reduction_embeddings,
    reduction_loadings_hash = reduction_loadings_hash,
    graph_names = names(obj@graphs),
    neighbor_names = names(obj@neighbors),
    command_names = names(obj@commands)
  )
  rm(obj)
  invisible(gc())
  signature
}

new_check_collector <- function() {
  rows <- list()
  index <- 1L
  add <- function(object_stage, component, passed, observed, expected, tolerance = NA_character_) {
    rows[[index]] <<- data.frame(
      object_stage = object_stage,
      component = component,
      passed = isTRUE(passed),
      observed = paste(observed, collapse = ";"),
      expected = paste(expected, collapse = ";"),
      tolerance = tolerance,
      stringsAsFactors = FALSE
    )
    index <<- index + 1L
  }
  collect <- function() {
    if (length(rows) == 0L) data.frame() else do.call(rbind, rows)
  }
  list(add = add, collect = collect)
}

numeric_difference <- function(x, y) {
  if (!identical(dim(x), dim(y))) return(Inf)
  if (length(x) == 0L) return(0)
  max(abs(as.numeric(x) - as.numeric(y)), na.rm = TRUE)
}

compare_metadata <- function(stage, generated, reference, collector) {
  collector$add(
    stage, "metadata_column_names",
    identical(colnames(generated$metadata), colnames(reference$metadata)),
    colnames(generated$metadata), colnames(reference$metadata)
  )
  collector$add(
    stage, "metadata_column_classes",
    identical(generated$metadata_classes, reference$metadata_classes),
    paste(names(generated$metadata_classes), generated$metadata_classes, sep = "="),
    paste(names(reference$metadata_classes), reference$metadata_classes, sep = "=")
  )
  common <- intersect(colnames(generated$metadata), colnames(reference$metadata))
  for (column_name in common) {
    x <- generated$metadata[[column_name]]
    y <- reference$metadata[[column_name]]
    if (is.numeric(x) && is.numeric(y)) {
      difference <- if (length(x) == length(y)) max(abs(x - y), na.rm = TRUE) else Inf
      if (!is.finite(difference) && all(is.na(x)) && all(is.na(y))) difference <- 0
      passed <- length(x) == length(y) && identical(is.na(x), is.na(y)) && difference <= 1e-10
      collector$add(
        stage, paste0("metadata_values__", column_name), passed,
        format(difference, digits = 16), "max_abs_difference<=1e-10", "1e-10"
      )
    } else {
      passed <- identical(as.character(x), as.character(y))
      collector$add(
        stage, paste0("metadata_values__", column_name), passed,
        hash_character(as.character(x)), hash_character(as.character(y))
      )
    }
    if (is.factor(x) || is.factor(y)) {
      collector$add(
        stage, paste0("metadata_factor_levels__", column_name),
        identical(levels(x), levels(y)), levels(x), levels(y)
      )
    }
  }
}

compare_signatures <- function(stage, generated, reference) {
  collector <- new_check_collector()
  collector$add(stage, "object_class", identical(generated$object_class, reference$object_class), generated$object_class, reference$object_class)
  collector$add(stage, "cell_count", length(generated$cells) == length(reference$cells), length(generated$cells), length(reference$cells))
  collector$add(stage, "cell_names_and_order", identical(generated$cells, reference$cells), hash_character(generated$cells), hash_character(reference$cells))
  collector$add(stage, "assay_names", identical(names(generated$assays), names(reference$assays)), names(generated$assays), names(reference$assays))
  collector$add(stage, "default_assay", identical(generated$default_assay, reference$default_assay), generated$default_assay, reference$default_assay)
  collector$add(stage, "graph_names", identical(generated$graph_names, reference$graph_names), generated$graph_names, reference$graph_names)
  collector$add(stage, "neighbor_names", identical(generated$neighbor_names, reference$neighbor_names), generated$neighbor_names, reference$neighbor_names)
  collector$add(stage, "command_names", identical(generated$command_names, reference$command_names), generated$command_names, reference$command_names)

  for (assay_name in intersect(names(generated$assays), names(reference$assays))) {
    x <- generated$assays[[assay_name]]
    y <- reference$assays[[assay_name]]
    collector$add(stage, paste0("assay_class__", assay_name), identical(x$class, y$class), x$class, y$class)
    collector$add(stage, paste0("assay_dimensions__", assay_name), identical(x$dimensions, y$dimensions), x$dimensions, y$dimensions)
    collector$add(stage, paste0("assay_features__", assay_name), identical(x$feature_names, y$feature_names), hash_character(x$feature_names), hash_character(y$feature_names))
    collector$add(stage, paste0("variable_features__", assay_name), identical(x$variable_features, y$variable_features), hash_character(x$variable_features), hash_character(y$variable_features))
    all_slots <- union(names(x$slot_hashes), names(y$slot_hashes))
    for (slot_name in all_slots) {
      x_hash <- unname(x$slot_hashes[slot_name])
      y_hash <- unname(y$slot_hashes[slot_name])
      collector$add(
        stage, paste0("assay_exact_hash__", assay_name, "__", slot_name),
        identical(x_hash, y_hash), x_hash, y_hash
      )
    }
  }

  compare_metadata(stage, generated, reference, collector)
  collector$add(
    stage, "active_ident_values",
    identical(as.character(generated$active_ident), as.character(reference$active_ident)),
    hash_character(as.character(generated$active_ident)),
    hash_character(as.character(reference$active_ident))
  )
  collector$add(
    stage, "active_ident_levels",
    identical(levels(generated$active_ident), levels(reference$active_ident)),
    levels(generated$active_ident), levels(reference$active_ident)
  )
  collector$add(stage, "reduction_names", identical(names(generated$reductions), names(reference$reductions)), names(generated$reductions), names(reference$reductions))
  for (reduction_name in intersect(names(generated$reductions), names(reference$reductions))) {
    x <- generated$reductions[[reduction_name]]
    y <- reference$reductions[[reduction_name]]
    difference <- numeric_difference(x, y)
    collector$add(
      stage, paste0("reduction_dimensions__", reduction_name),
      identical(dim(x), dim(y)), dim(x), dim(y)
    )
    collector$add(
      stage, paste0("reduction_cell_order__", reduction_name),
      identical(rownames(x), rownames(y)), hash_character(rownames(x)), hash_character(rownames(y))
    )
    collector$add(
      stage, paste0("reduction_values__", reduction_name),
      is.finite(difference) && difference <= 1e-8,
      format(difference, digits = 16), "max_abs_difference<=1e-8", "1e-8"
    )
    collector$add(
      stage, paste0("reduction_loadings_exact_hash__", reduction_name),
      identical(generated$reduction_loadings_hash[[reduction_name]], reference$reduction_loadings_hash[[reduction_name]]),
      generated$reduction_loadings_hash[[reduction_name]], reference$reduction_loadings_hash[[reduction_name]]
    )
  }
  collector$collect()
}

validation_root <- ensure_dir(file.path(generated_root, "00_validation"))
all_reports <- list()
object_pairs <- list(
  refined = c(generated_refined, reference_refined),
  final = c(generated_final, reference_final)
)

for (stage in names(object_pairs)) {
  generated_signature <- object_signature(object_pairs[[stage]][[1]])
  reference_signature <- object_signature(object_pairs[[stage]][[2]])
  all_reports[[stage]] <- compare_signatures(stage, generated_signature, reference_signature)
  rm(generated_signature, reference_signature)
  invisible(gc())
}

report <- do.call(rbind, all_reports)
utils::write.csv(report, file.path(validation_root, "reproduction_validation_report.csv"), row.names = FALSE)
failed <- report[!report$passed, , drop = FALSE]
utils::write.csv(failed, file.path(validation_root, "reproduction_validation_failures.csv"), row.names = FALSE)
summary_lines <- c(
  paste0("status=", if (nrow(failed) == 0L) "PASS" else "FAIL"),
  paste0("checks_total=", nrow(report)),
  paste0("checks_passed=", sum(report$passed)),
  paste0("checks_failed=", nrow(failed)),
  paste0("generated_refined=", generated_refined),
  paste0("reference_refined=", reference_refined),
  paste0("generated_final=", generated_final),
  paste0("reference_final=", reference_final),
  paste0("validated_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
)
writeLines(summary_lines, file.path(validation_root, "validation_summary.txt"))

if (nrow(failed) > 0L) {
  message("Reproduction validation FAILED with ", nrow(failed), " failed checks.")
  quit(save = "no", status = 1L)
}
writeLines(summary_lines, file.path(validation_root, "VALIDATION_PASS.txt"))
message("Reproduction validation PASSED.")
