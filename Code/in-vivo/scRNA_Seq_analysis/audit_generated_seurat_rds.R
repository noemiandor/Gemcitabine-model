#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

audit_schema_version <- "semantic_rds_audit_v1"
numeric_max_abs_tolerance <- 1e-6
numeric_rmse_tolerance <- 1e-8
numeric_correlation_tolerance <- 0.999999
umap_correlation_tolerance <- 0.99
umap_median_displacement_tolerance <- 0.10
umap_q99_displacement_tolerance <- 0.50
umap_fraction_gt_one_tolerance <- 0.005
umap_max_displacement_tolerance <- 15

audit_stop <- function(...) stop(..., call. = FALSE)

audit_require_packages <- function() {
  for (package in c("Seurat", "Matrix", "digest")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      audit_stop("Semantic RDS audit requires package '", package, "'.")
    }
  }
}

audit_ensure_dir <- function(path) {
  if (!dir.exists(path) &&
      !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    audit_stop("Cannot create semantic audit directory: ", path)
  }
  normalizePath(path, mustWork = TRUE)
}

audit_sha256 <- function(path) {
  tolower(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

audit_file_identity <- function(path, artifact, role) {
  path <- normalizePath(path, mustWork = TRUE)
  data.frame(
    artifact = artifact,
    path = path,
    size_bytes = as.character(file.info(path)$size),
    md5 = tolower(unname(tools::md5sum(path))),
    sha256 = audit_sha256(path),
    role = role,
    stringsAsFactors = FALSE
  )
}

audit_text <- function(x) {
  if (length(x) == 0L || is.null(x)) return("")
  value <- paste(as.character(x), collapse = ";")
  gsub("[\t\r\n]+", " ", value)
}

audit_row <- function(
  check,
  passed,
  observed = "",
  reference = "",
  threshold = "",
  details = ""
) {
  data.frame(
    check = as.character(check),
    status = if (isTRUE(passed)) "PASS" else "FAIL",
    observed = audit_text(observed),
    reference = audit_text(reference),
    threshold = audit_text(threshold),
    details = audit_text(details),
    stringsAsFactors = FALSE
  )
}

audit_bind_rows <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    return(data.frame(
      check = character(), status = character(), observed = character(),
      reference = character(), threshold = character(), details = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

audit_write_tsv <- function(x, path) {
  utils::write.table(
    x,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

audit_exact_summary <- function(x, y) {
  identical(x, y)
}

audit_numeric_vectors <- function(x, y) {
  if (!identical(length(x), length(y)) ||
      !identical(is.na(x), is.na(y)) ||
      !identical(is.nan(x), is.nan(y)) ||
      !identical(is.infinite(x), is.infinite(y))) {
    return(list(
      compatible = FALSE,
      exact = FALSE,
      n = max(length(x), length(y)),
      max_abs = Inf,
      rmse = Inf,
      correlation = NA_real_
    ))
  }
  keep <- is.finite(x) & is.finite(y)
  x <- as.numeric(x[keep])
  y <- as.numeric(y[keep])
  if (!length(x)) {
    return(list(
      compatible = TRUE,
      exact = TRUE,
      n = 0,
      max_abs = 0,
      rmse = 0,
      correlation = 1
    ))
  }
  delta <- x - y
  max_abs <- max(abs(delta))
  rmse <- sqrt(mean(delta * delta))
  exact <- identical(x, y)
  x_centered <- x - mean(x)
  y_centered <- y - mean(y)
  denominator <- sqrt(sum(x_centered * x_centered) * sum(y_centered * y_centered))
  correlation <- if (is.finite(denominator) && denominator > 0) {
    sum(x_centered * y_centered) / denominator
  } else if (max_abs <= numeric_max_abs_tolerance) {
    1
  } else {
    NA_real_
  }
  list(
    compatible = TRUE,
    exact = exact,
    n = length(x),
    max_abs = max_abs,
    rmse = rmse,
    correlation = correlation
  )
}

audit_numeric_matrix <- function(x, y, block_columns = 256L) {
  if (!identical(dim(x), dim(y)) ||
      !identical(dimnames(x), dimnames(y))) {
    return(list(
      compatible = FALSE,
      exact = FALSE,
      n = max(length(x), length(y)),
      max_abs = Inf,
      rmse = Inf,
      correlation = NA_real_
    ))
  }
  if (!length(x)) return(audit_numeric_vectors(numeric(), numeric()))
  sparse_x <- inherits(x, "sparseMatrix")
  sparse_y <- inherits(y, "sparseMatrix")
  if (sparse_x && sparse_y) {
    x_c <- methods::as(x, "dgCMatrix")
    y_c <- methods::as(y, "dgCMatrix")
    if (!identical(x_c@p, y_c@p) || !identical(x_c@i, y_c@i)) {
      return(list(
        compatible = FALSE,
        exact = FALSE,
        n = length(x),
        max_abs = Inf,
        rmse = Inf,
        correlation = NA_real_
      ))
    }
    stats <- audit_numeric_vectors(x_c@x, y_c@x)
    stats$n <- length(x)
    return(stats)
  }
  if (xor(sparse_x, sparse_y)) {
    return(list(
      compatible = FALSE,
      exact = FALSE,
      n = length(x),
      max_abs = Inf,
      rmse = Inf,
      correlation = NA_real_
    ))
  }

  n_columns <- ncol(x)
  starts <- seq.int(1L, n_columns, by = block_columns)
  total_n <- 0
  sum_x <- sum_y <- sum_x2 <- sum_y2 <- sum_xy <- sum_delta2 <- 0
  max_abs <- 0
  exact <- TRUE
  for (start in starts) {
    end <- min(n_columns, start + block_columns - 1L)
    xv <- as.numeric(x[, start:end, drop = FALSE])
    yv <- as.numeric(y[, start:end, drop = FALSE])
    if (!identical(is.na(xv), is.na(yv)) ||
        !identical(is.nan(xv), is.nan(yv)) ||
        !identical(is.infinite(xv), is.infinite(yv))) {
      return(list(
        compatible = FALSE,
        exact = FALSE,
        n = length(x),
        max_abs = Inf,
        rmse = Inf,
        correlation = NA_real_
      ))
    }
    keep <- is.finite(xv) & is.finite(yv)
    xv <- xv[keep]
    yv <- yv[keep]
    if (!length(xv)) next
    delta <- xv - yv
    exact <- exact && identical(xv, yv)
    total_n <- total_n + length(xv)
    max_abs <- max(max_abs, max(abs(delta)))
    sum_delta2 <- sum_delta2 + sum(delta * delta)
    sum_x <- sum_x + sum(xv)
    sum_y <- sum_y + sum(yv)
    sum_x2 <- sum_x2 + sum(xv * xv)
    sum_y2 <- sum_y2 + sum(yv * yv)
    sum_xy <- sum_xy + sum(xv * yv)
  }
  if (!total_n) return(audit_numeric_vectors(numeric(), numeric()))
  variance_x <- sum_x2 - (sum_x * sum_x / total_n)
  variance_y <- sum_y2 - (sum_y * sum_y / total_n)
  denominator <- sqrt(max(0, variance_x) * max(0, variance_y))
  correlation <- if (is.finite(denominator) && denominator > 0) {
    (sum_xy - (sum_x * sum_y / total_n)) / denominator
  } else if (max_abs <= numeric_max_abs_tolerance) {
    1
  } else {
    NA_real_
  }
  list(
    compatible = TRUE,
    exact = exact,
    n = total_n,
    max_abs = max_abs,
    rmse = sqrt(sum_delta2 / total_n),
    correlation = correlation
  )
}

audit_numeric_pass <- function(stats) {
  isTRUE(stats$compatible) &&
    is.finite(stats$max_abs) &&
    stats$max_abs <= numeric_max_abs_tolerance &&
    is.finite(stats$rmse) &&
    stats$rmse <= numeric_rmse_tolerance &&
    is.finite(stats$correlation) &&
    stats$correlation >= numeric_correlation_tolerance
}

audit_numeric_row <- function(check, stats, details = "") {
  audit_row(
    check,
    audit_numeric_pass(stats),
    sprintf(
      "compatible=%s;n=%s;exact=%s;max_abs=%.17g;rmse=%.17g;correlation=%.17g",
      tolower(as.character(stats$compatible)),
      stats$n,
      tolower(as.character(stats$exact)),
      stats$max_abs,
      stats$rmse,
      stats$correlation
    ),
    "numerically equivalent",
    paste0(
      "max_abs<=", numeric_max_abs_tolerance,
      ";rmse<=", numeric_rmse_tolerance,
      ";correlation>=", numeric_correlation_tolerance
    ),
    details
  )
}

audit_metadata <- function(generated, reference) {
  x <- generated@meta.data
  y <- reference@meta.data
  rows <- list(
    audit_row("row_names_and_order", identical(rownames(x), rownames(y))),
    audit_row(
      "column_names_and_order",
      identical(colnames(x), colnames(y)),
      colnames(x),
      colnames(y)
    )
  )
  for (name in union(colnames(x), colnames(y))) {
    present <- name %in% colnames(x) && name %in% colnames(y)
    if (!present) {
      rows[[length(rows) + 1L]] <- audit_row(
        paste0("column:", name), FALSE,
        name %in% colnames(x), name %in% colnames(y),
        details = "column presence"
      )
      next
    }
    gx <- x[[name]]
    gy <- y[[name]]
    rows[[length(rows) + 1L]] <- audit_row(
      paste0("column:", name),
      identical(gx, gy),
      paste0(
        "class=", paste(class(gx), collapse = "/"),
        ";factor_levels=", paste(if (is.factor(gx)) levels(gx) else character(), collapse = ","),
        ";na=", sum(is.na(gx))
      ),
      paste0(
        "class=", paste(class(gy), collapse = "/"),
        ";factor_levels=", paste(if (is.factor(gy)) levels(gy) else character(), collapse = ","),
        ";na=", sum(is.na(gy))
      ),
      "exact type, attributes, NA mask, and values"
    )
  }
  audit_bind_rows(rows)
}

audit_clusters <- function(generated, reference) {
  x_meta <- generated@meta.data
  y_meta <- reference@meta.data
  cluster_columns <- union(
    grep("cluster", colnames(x_meta), ignore.case = TRUE, value = TRUE),
    grep("cluster", colnames(y_meta), ignore.case = TRUE, value = TRUE)
  )
  rows <- list(
    audit_row(
      "active_ident",
      identical(generated@active.ident, reference@active.ident),
      paste0("levels=", paste(levels(generated@active.ident), collapse = ",")),
      paste0("levels=", paste(levels(reference@active.ident), collapse = ",")),
      "exact names, order, type, levels, NA mask, and values"
    ),
    audit_row(
      "cluster_column_inventory",
      identical(cluster_columns, intersect(colnames(x_meta), cluster_columns)) &&
        identical(cluster_columns, intersect(colnames(y_meta), cluster_columns)),
      intersect(colnames(x_meta), cluster_columns),
      intersect(colnames(y_meta), cluster_columns)
    )
  )
  for (name in cluster_columns) {
    present <- name %in% colnames(x_meta) && name %in% colnames(y_meta)
    rows[[length(rows) + 1L]] <- audit_row(
      paste0("metadata_cluster:", name),
      present && identical(x_meta[[name]], y_meta[[name]]),
      if (name %in% colnames(x_meta)) paste0("levels=", paste(levels(x_meta[[name]]), collapse = ",")) else "missing",
      if (name %in% colnames(y_meta)) paste0("levels=", paste(levels(y_meta[[name]]), collapse = ",")) else "missing",
      "exact type, levels, NA mask, and values"
    )
  }
  audit_bind_rows(rows)
}

audit_graphs <- function(generated, reference) {
  x_names <- names(generated@graphs)
  y_names <- names(reference@graphs)
  rows <- list(audit_row(
    "graph_inventory", identical(x_names, y_names), x_names, y_names
  ))
  for (name in union(x_names, y_names)) {
    present <- name %in% x_names && name %in% y_names
    rows[[length(rows) + 1L]] <- audit_row(
      paste0("graph:", name),
      present && identical(generated@graphs[[name]], reference@graphs[[name]]),
      if (present) paste0("class=", class(generated@graphs[[name]])[[1L]], ";dim=", paste(dim(generated@graphs[[name]]), collapse = "x")) else "missing",
      if (present) paste0("class=", class(reference@graphs[[name]])[[1L]], ";dim=", paste(dim(reference@graphs[[name]]), collapse = "x")) else "missing",
      "exact graph structure, names, and values"
    )
  }
  audit_bind_rows(rows)
}

audit_assays <- function(generated, reference) {
  x_names <- names(generated@assays)
  y_names <- names(reference@assays)
  rows <- list(
    audit_row("assay_inventory", identical(x_names, y_names), x_names, y_names),
    audit_row(
      "active_assay",
      identical(Seurat::DefaultAssay(generated), Seurat::DefaultAssay(reference)),
      Seurat::DefaultAssay(generated),
      Seurat::DefaultAssay(reference)
    )
  )
  for (name in union(x_names, y_names)) {
    present <- name %in% x_names && name %in% y_names
    if (!present) {
      rows[[length(rows) + 1L]] <- audit_row(
        paste0("assay:", name, ":presence"), FALSE,
        name %in% x_names, name %in% y_names
      )
      next
    }
    x <- generated@assays[[name]]
    y <- reference@assays[[name]]
    rows[[length(rows) + 1L]] <- audit_row(
      paste0("assay:", name, ":structure"),
      identical(class(x), class(y)) &&
        identical(methods::slotNames(x), methods::slotNames(y)) &&
        identical(dim(x), dim(y)) &&
        identical(dimnames(x), dimnames(y)) &&
        identical(x@key, y@key) &&
        identical(x@var.features, y@var.features) &&
        identical(x@meta.features, y@meta.features),
      paste0("class=", paste(class(x), collapse = "/"), ";dim=", paste(dim(x), collapse = "x")),
      paste0("class=", paste(class(y), collapse = "/"), ";dim=", paste(dim(y), collapse = "x")),
      "exact class, slots, dimensions, names, key, variable features, and feature metadata"
    )
    if ("SCTModel.list" %in% methods::slotNames(x)) {
      rows[[length(rows) + 1L]] <- audit_row(
        paste0("assay:", name, ":sct_model_structure"),
        identical(names(x@SCTModel.list), names(y@SCTModel.list)) &&
          identical(
            vapply(x@SCTModel.list, class, character(1L)),
            vapply(y@SCTModel.list, class, character(1L))
          ),
        names(x@SCTModel.list), names(y@SCTModel.list),
        "exact SCT model inventory and classes"
      )
    }
    numeric_slots <- intersect(
      c("counts", "data", "scale.data"),
      methods::slotNames(x)
    )
    for (slot_name in numeric_slots) {
      x_value <- methods::slot(x, slot_name)
      y_value <- methods::slot(y, slot_name)
      check <- paste0("assay:", name, ":", slot_name)
      if (identical(slot_name, "counts")) {
        rows[[length(rows) + 1L]] <- audit_row(
          check,
          identical(x_value, y_value),
          paste0("class=", class(x_value)[[1L]], ";dim=", paste(dim(x_value), collapse = "x")),
          paste0("class=", class(y_value)[[1L]], ";dim=", paste(dim(y_value), collapse = "x")),
          "exact count matrix representation, names, and values"
        )
      } else {
        rows[[length(rows) + 1L]] <- audit_numeric_row(
          check,
          audit_numeric_matrix(x_value, y_value),
          "row and column names/order must also be exact"
        )
      }
    }
  }
  audit_bind_rows(rows)
}

audit_pca <- function(generated, reference) {
  x_names <- names(generated@reductions)
  y_names <- names(reference@reductions)
  rows <- list(audit_row(
    "reduction_inventory", identical(x_names, y_names), x_names, y_names
  ))
  if (!"pca" %in% x_names || !"pca" %in% y_names) {
    rows[[length(rows) + 1L]] <- audit_row(
      "pca:presence", FALSE, "pca" %in% x_names, "pca" %in% y_names
    )
    return(audit_bind_rows(rows))
  }
  x <- generated@reductions[["pca"]]
  y <- reference@reductions[["pca"]]
  rows[[length(rows) + 1L]] <- audit_row(
    "pca:structure",
    identical(class(x), class(y)) &&
      identical(x@assay.used, y@assay.used) &&
      identical(x@key, y@key),
    paste0("class=", class(x)[[1L]], ";assay=", x@assay.used, ";key=", x@key),
    paste0("class=", class(y)[[1L]], ";assay=", y@assay.used, ";key=", y@key),
    "exact reduction class, assay, and key"
  )
  for (slot_name in intersect(
    c("cell.embeddings", "feature.loadings", "feature.loadings.projected", "stdev"),
    methods::slotNames(x)
  )) {
    x_value <- methods::slot(x, slot_name)
    y_value <- methods::slot(y, slot_name)
    stats <- if (is.null(dim(x_value))) {
      if (!identical(names(x_value), names(y_value))) {
        list(compatible = FALSE, exact = FALSE, n = max(length(x_value), length(y_value)), max_abs = Inf, rmse = Inf, correlation = NA_real_)
      } else {
        audit_numeric_vectors(x_value, y_value)
      }
    } else {
      audit_numeric_matrix(x_value, y_value)
    }
    rows[[length(rows) + 1L]] <- audit_numeric_row(
      paste0("pca:", slot_name), stats,
      "dimensions and names/order must also be exact"
    )
  }
  audit_bind_rows(rows)
}

audit_umap <- function(generated, reference) {
  empty_top <- data.frame(
    cell = character(),
    displacement = numeric(),
    generated_UMAP_1 = numeric(),
    generated_UMAP_2 = numeric(),
    reference_UMAP_1 = numeric(),
    reference_UMAP_2 = numeric(),
    delta_UMAP_1 = numeric(),
    delta_UMAP_2 = numeric(),
    stringsAsFactors = FALSE
  )
  if (!"umap" %in% names(generated@reductions) ||
      !"umap" %in% names(reference@reductions)) {
    return(list(
      summary = audit_row(
        "umap:presence", FALSE,
        "umap" %in% names(generated@reductions),
        "umap" %in% names(reference@reductions)
      ),
      largest = empty_top
    ))
  }
  x <- generated@reductions[["umap"]]@cell.embeddings
  y <- reference@reductions[["umap"]]@cell.embeddings
  if (!identical(dim(x), dim(y)) ||
      !identical(dimnames(x), dimnames(y)) || ncol(x) != 2L) {
    return(list(
      summary = audit_row(
        "umap:structure", FALSE,
        paste0("dim=", paste(dim(x), collapse = "x")),
        paste0("dim=", paste(dim(y), collapse = "x")),
        "exact 2D dimensions and cell/axis order"
      ),
      largest = empty_top
    ))
  }
  delta <- x - y
  displacement <- sqrt(rowSums(delta * delta))
  correlations <- vapply(seq_len(ncol(x)), function(index) {
    stats::cor(x[, index], y[, index])
  }, numeric(1L))
  median_displacement <- stats::median(displacement)
  q99_displacement <- unname(stats::quantile(displacement, 0.99, type = 7))
  fraction_gt_one <- mean(displacement > 1)
  max_displacement <- max(displacement)
  rows <- list(
    audit_row(
      "umap:structure", TRUE,
      paste0("dim=", paste(dim(x), collapse = "x")),
      paste0("dim=", paste(dim(y), collapse = "x")),
      "exact dimensions and cell/axis order"
    ),
    audit_row(
      "umap:coordinate_correlation",
      all(is.finite(correlations)) && min(correlations) >= umap_correlation_tolerance,
      paste(sprintf("%s=%.17g", colnames(x), correlations), collapse = ";"),
      "same-axis positive correlation",
      paste0("each_axis>=", umap_correlation_tolerance)
    ),
    audit_row(
      "umap:median_cell_displacement",
      median_displacement <= umap_median_displacement_tolerance,
      sprintf("%.17g", median_displacement), "0",
      paste0("<=", umap_median_displacement_tolerance)
    ),
    audit_row(
      "umap:q99_cell_displacement",
      q99_displacement <= umap_q99_displacement_tolerance,
      sprintf("%.17g", q99_displacement), "0",
      paste0("<=", umap_q99_displacement_tolerance)
    ),
    audit_row(
      "umap:fraction_cell_displacement_gt_1",
      fraction_gt_one <= umap_fraction_gt_one_tolerance,
      sprintf("%.17g", fraction_gt_one), "0",
      paste0("<=", umap_fraction_gt_one_tolerance)
    ),
    audit_row(
      "umap:max_cell_displacement",
      max_displacement <= umap_max_displacement_tolerance,
      sprintf("%.17g", max_displacement), "0",
      paste0("<=", umap_max_displacement_tolerance)
    )
  )
  order_index <- head(order(displacement, decreasing = TRUE), 100L)
  largest <- data.frame(
    cell = rownames(x)[order_index],
    displacement = displacement[order_index],
    generated_UMAP_1 = x[order_index, 1L],
    generated_UMAP_2 = x[order_index, 2L],
    reference_UMAP_1 = y[order_index, 1L],
    reference_UMAP_2 = y[order_index, 2L],
    delta_UMAP_1 = delta[order_index, 1L],
    delta_UMAP_2 = delta[order_index, 2L],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  metadata <- generated@meta.data[rownames(x)[order_index], , drop = FALSE]
  colnames(metadata) <- paste0("meta.", colnames(metadata))
  largest <- cbind(largest, metadata)
  list(summary = audit_bind_rows(rows), largest = largest)
}

audit_normalize_command <- function(command) {
  slots <- setdiff(methods::slotNames(command), "time.stamp")
  stats::setNames(lapply(slots, function(name) methods::slot(command, name)), slots)
}

audit_commands <- function(generated, reference) {
  x_names <- names(generated@commands)
  y_names <- names(reference@commands)
  rows <- list(audit_row(
    "command_inventory", identical(x_names, y_names), x_names, y_names
  ))
  for (name in union(x_names, y_names)) {
    present <- name %in% x_names && name %in% y_names
    timestamp_equal <- present && identical(
      generated@commands[[name]]@time.stamp,
      reference@commands[[name]]@time.stamp
    )
    rows[[length(rows) + 1L]] <- audit_row(
      paste0("command:", name),
      present && identical(
        audit_normalize_command(generated@commands[[name]]),
        audit_normalize_command(reference@commands[[name]])
      ),
      if (present) paste0("timestamp_equal=", tolower(as.character(timestamp_equal))) else "missing",
      "parameters, call string, assay, and recorded seeds",
      "exact after excluding time.stamp",
      "time.stamp is intentionally non-gating"
    )
  }
  audit_bind_rows(rows)
}

audit_object_contract <- function(generated, reference) {
  rows <- list(
    audit_row(
      "object_class",
      identical(class(generated), class(reference)),
      class(generated), class(reference)
    ),
    audit_row(
      "object_version",
      identical(generated@version, reference@version),
      as.character(generated@version), as.character(reference@version)
    ),
    audit_row(
      "project_name",
      identical(generated@project.name, reference@project.name),
      generated@project.name, reference@project.name
    ),
    audit_row(
      "cell_names_and_order",
      identical(colnames(generated), colnames(reference)),
      length(colnames(generated)), length(colnames(reference)),
      "exact cell names and order"
    ),
    audit_row(
      "misc",
      identical(generated@misc, reference@misc),
      names(generated@misc), names(reference@misc), "exact"
    ),
    audit_row(
      "tools",
      identical(generated@tools, reference@tools),
      names(generated@tools), names(reference@tools), "exact"
    )
  )
  audit_bind_rows(rows)
}

audit_generated_seurat_rds <- function(
  generated_rds,
  zenodo_reference_rds,
  output_dir,
  stop_on_failure = TRUE
) {
  audit_require_packages()
  suppressPackageStartupMessages(library(Seurat))
  generated_rds <- normalizePath(generated_rds, mustWork = TRUE)
  zenodo_reference_rds <- normalizePath(zenodo_reference_rds, mustWork = TRUE)
  output_dir <- audit_ensure_dir(output_dir)

  identities <- rbind(
    audit_file_identity(generated_rds, "generated", "candidate_downstream_input"),
    audit_file_identity(zenodo_reference_rds, "zenodo_reference", "semantic_reference")
  )
  audit_write_tsv(identities, file.path(output_dir, "rds_file_identity.tsv"))

  message("Reading generated Seurat RDS: ", generated_rds)
  generated <- readRDS(generated_rds)
  message("Reading Zenodo reference Seurat RDS: ", zenodo_reference_rds)
  reference <- readRDS(zenodo_reference_rds)
  if (!inherits(generated, "Seurat") || !inherits(reference, "Seurat")) {
    audit_stop("Both semantic-audit inputs must be Seurat objects")
  }

  reports <- list(
    object = audit_object_contract(generated, reference),
    metadata = audit_metadata(generated, reference),
    cluster = audit_clusters(generated, reference),
    graph = audit_graphs(generated, reference),
    assay = audit_assays(generated, reference),
    pca = audit_pca(generated, reference)
  )
  umap <- audit_umap(generated, reference)
  reports$umap <- umap$summary
  reports$command <- audit_commands(generated, reference)

  audit_write_tsv(reports$metadata, file.path(output_dir, "metadata_comparison.tsv"))
  audit_write_tsv(reports$cluster, file.path(output_dir, "cluster_comparison.tsv"))
  audit_write_tsv(reports$graph, file.path(output_dir, "graph_comparison.tsv"))
  audit_write_tsv(reports$assay, file.path(output_dir, "assay_numeric_comparison.tsv"))
  audit_write_tsv(reports$pca, file.path(output_dir, "pca_comparison.tsv"))
  audit_write_tsv(reports$umap, file.path(output_dir, "umap_displacement_summary.tsv"))
  audit_write_tsv(umap$largest, file.path(output_dir, "umap_largest_displacements.tsv"))
  audit_write_tsv(reports$command, file.path(output_dir, "command_comparison.tsv"))

  summary <- do.call(rbind, lapply(names(reports), function(group) {
    data.frame(group = group, reports[[group]], stringsAsFactors = FALSE)
  }))
  rownames(summary) <- NULL
  audit_write_tsv(summary, file.path(output_dir, "audit_summary.tsv"))
  failed <- summary$status != "PASS"
  status <- if (any(failed)) "FAIL" else "PASS"
  summary_path <- normalizePath(
    file.path(output_dir, "audit_summary.tsv"), mustWork = TRUE
  )
  report_names <- c(
    "audit_summary.tsv",
    "metadata_comparison.tsv",
    "cluster_comparison.tsv",
    "graph_comparison.tsv",
    "assay_numeric_comparison.tsv",
    "pca_comparison.tsv",
    "umap_displacement_summary.tsv",
    "umap_largest_displacements.tsv",
    "command_comparison.tsv",
    "rds_file_identity.tsv"
  )
  report_hash_lines <- vapply(report_names, function(name) {
    paste0(
      "report_sha256.", name, "=",
      audit_sha256(file.path(output_dir, name))
    )
  }, character(1L))
  marker <- c(
    paste0("schema_version=", audit_schema_version),
    paste0("status=", status),
    paste0("checks_total=", nrow(summary)),
    paste0("checks_passed=", sum(!failed)),
    paste0("checks_failed=", sum(failed)),
    paste0("generated_rds=", generated_rds),
    paste0("generated_rds_size_bytes=", identities$size_bytes[[1L]]),
    paste0("generated_rds_md5=", identities$md5[[1L]]),
    paste0("generated_rds_sha256=", identities$sha256[[1L]]),
    paste0("zenodo_reference_rds=", zenodo_reference_rds),
    paste0("zenodo_reference_rds_size_bytes=", identities$size_bytes[[2L]]),
    paste0("zenodo_reference_rds_md5=", identities$md5[[2L]]),
    paste0("zenodo_reference_rds_sha256=", identities$sha256[[2L]]),
    paste0("downstream_rds=", generated_rds),
    paste0("audit_summary=", summary_path),
    paste0("audit_summary_sha256=", audit_sha256(summary_path)),
    report_hash_lines,
    paste0("validated_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  )
  marker_path <- file.path(output_dir, "AUDIT_COMPLETE.txt")
  writeLines(marker, marker_path)

  rm(generated, reference)
  invisible(gc())
  if (identical(status, "FAIL") && isTRUE(stop_on_failure)) {
    message(
      "Semantic RDS audit FAILED with ", sum(failed),
      " failed checks; see ", summary_path
    )
    stop("Generated Seurat RDS failed semantic audit", call. = FALSE)
  }
  message("Semantic RDS audit ", status, ": ", marker_path)
  invisible(list(status = status, summary = summary, identities = identities))
}

audit_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L) {
    audit_stop(
      "Usage: audit_generated_seurat_rds.R GENERATED_RDS ",
      "ZENODO_REFERENCE_RDS OUTPUT_DIR"
    )
  }
  audit_generated_seurat_rds(args[[1L]], args[[2L]], args[[3L]])
}

if (sys.nframe() == 0L) audit_cli()
