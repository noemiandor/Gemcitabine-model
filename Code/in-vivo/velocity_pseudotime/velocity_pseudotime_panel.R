velocity_stop <- function(...) {
  stop(..., call. = FALSE)
}

velocity_sha256 <- function(path) {
  if (!file.exists(path)) velocity_stop("Cannot checksum missing file: ", path)
  if (!requireNamespace("digest", quietly = TRUE)) {
    velocity_stop("R package 'digest' is required")
  }
  unname(digest::digest(
    path,
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  ))
}

velocity_read_tsv <- function(path) {
  if (!file.exists(path)) velocity_stop("Missing input: ", path)
  utils::read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    na.strings = c("", "NA", "NaN")
  )
}

velocity_read_csv <- function(path) {
  if (!file.exists(path)) velocity_stop("Missing input: ", path)
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN")
  )
}

velocity_write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    data,
    path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = ""
  )
  invisible(path)
}

velocity_canonical_table_locator <- function() {
  "Data/in-vivo/figure7/processed/cellcycle_velocity_pseudotime_umap.tsv"
}

velocity_canonical_provenance_locator <- function() {
  paste0(
    "Data/in-vivo/figure7/processed/",
    "cellcycle_velocity_pseudotime_umap.provenance.tsv"
  )
}

velocity_path_is_within <- function(path, root) {
  normalized <- normalizePath(path, mustWork = FALSE)
  normalized_root <- normalizePath(root, mustWork = TRUE)
  identical(normalized, normalized_root) || startsWith(
    normalized,
    paste0(normalized_root, .Platform$file.sep)
  )
}

velocity_portable_locator <- function(path, repo_root) {
  normalized <- normalizePath(path, mustWork = TRUE)
  normalized_root <- normalizePath(repo_root, mustWork = TRUE)
  if (velocity_path_is_within(normalized, normalized_root)) {
    if (identical(normalized, normalized_root)) return(".")
    return(substring(
      normalized,
      nchar(paste0(normalized_root, .Platform$file.sep)) + 1L
    ))
  }
  paste0("external:", basename(normalized))
}

velocity_resolve_repo_locator <- function(locator, repo_root, must_work = TRUE) {
  locator <- trimws(as.character(locator))
  if (length(locator) != 1L || is.na(locator) || !nzchar(locator) ||
      startsWith(locator, "external:") || grepl("^[/\\\\]", locator) ||
      grepl("^[A-Za-z]:[/\\\\]", locator) ||
      any(strsplit(locator, "[/\\\\]")[[1L]] %in% c("", ".", ".."))) {
    velocity_stop("Unsafe or non-repository velocity locator: ", locator)
  }
  candidate <- file.path(repo_root, locator)
  if (!velocity_path_is_within(candidate, repo_root)) {
    velocity_stop("Velocity locator escapes the repository: ", locator)
  }
  normalizePath(candidate, mustWork = must_work)
}

velocity_is_canonical_bundle <- function(table_path, provenance_path, repo_root) {
  identical(
    normalizePath(table_path, mustWork = TRUE),
    normalizePath(
      file.path(repo_root, velocity_canonical_table_locator()),
      mustWork = FALSE
    )
  ) && identical(
    normalizePath(provenance_path, mustWork = TRUE),
    normalizePath(
      file.path(repo_root, velocity_canonical_provenance_locator()),
      mustWork = FALSE
    )
  )
}

velocity_parse_logical <- function(value, field = "logical value") {
  normalized <- tolower(trimws(as.character(value)))
  out <- normalized %in% c("true", "t", "1", "yes")
  invalid <- is.na(normalized) |
    !normalized %in% c("true", "t", "1", "yes", "false", "f", "0", "no")
  if (any(invalid)) velocity_stop("Invalid ", field)
  out
}

velocity_required_columns <- function() {
  c(
    "cell_id", "sample_id", "cluster", "UMAP_1", "UMAP_2",
    "velocity_UMAP_1", "velocity_UMAP_2", "velocity_pseudotime",
    "is_root"
  )
}

velocity_expected_cluster_counts <- function() {
  c("4c" = 89L, "6" = 686L, "10" = 2106L)
}

velocity_validate_panel_table <- function(
  table,
  expected_cluster_counts = velocity_expected_cluster_counts(),
  root_cluster = "6"
) {
  required <- velocity_required_columns()
  if (!identical(names(table), required)) {
    velocity_stop(
      "Velocity panel table columns must be exactly: ",
      paste(required, collapse = ", ")
    )
  }
  cells <- trimws(as.character(table$cell_id))
  samples <- trimws(as.character(table$sample_id))
  clusters <- trimws(as.character(table$cluster))
  if (!nrow(table) || anyNA(cells) || any(!nzchar(cells)) ||
      anyDuplicated(cells) || anyNA(samples) || any(!nzchar(samples)) ||
      anyNA(clusters) || any(!nzchar(clusters))) {
    velocity_stop(
      "Velocity panel table requires unique cell IDs and complete sample/cluster fields"
    )
  }
  expected_cluster_names <- names(expected_cluster_counts)
  expected_cluster_counts <- as.integer(expected_cluster_counts)
  names(expected_cluster_counts) <- expected_cluster_names
  if (is.null(expected_cluster_names) ||
      anyNA(expected_cluster_names) ||
      any(!nzchar(expected_cluster_names)) ||
      anyDuplicated(expected_cluster_names)) {
    velocity_stop("Expected cluster counts require unique cluster names")
  }
  observed <- table(factor(clusters, levels = names(expected_cluster_counts)))
  if (!setequal(unique(clusters), names(expected_cluster_counts)) ||
      !identical(as.integer(observed), unname(expected_cluster_counts))) {
    velocity_stop(
      "Velocity panel cluster counts differ from the exact CellCycle universe"
    )
  }
  numeric_fields <- c(
    "UMAP_1", "UMAP_2", "velocity_UMAP_1",
    "velocity_UMAP_2", "velocity_pseudotime"
  )
  for (field in numeric_fields) {
    values <- suppressWarnings(as.numeric(table[[field]]))
    if (any(!is.finite(values))) {
      velocity_stop("Velocity panel field must be finite: ", field)
    }
    table[[field]] <- values
  }
  if (any(table$velocity_pseudotime < 0 | table$velocity_pseudotime > 1)) {
    velocity_stop("Velocity pseudotime must lie within [0,1]")
  }
  vector_length <- sqrt(
    table$velocity_UMAP_1^2 + table$velocity_UMAP_2^2
  )
  if (!any(vector_length > 0)) {
    velocity_stop("Velocity panel requires at least one non-zero RNA-velocity vector")
  }
  roots <- velocity_parse_logical(table$is_root, "is_root value")
  expected_roots <- clusters == root_cluster
  if (!root_cluster %in% clusters ||
      !identical(roots, expected_roots)) {
    velocity_stop(
      "is_root must identify every and only cluster ", root_cluster,
      " cell"
    )
  }
  table$cell_id <- cells
  table$sample_id <- samples
  table$cluster <- clusters
  table$is_root <- roots
  table
}

velocity_validate_provenance <- function(
  path,
  table_path,
  root_cluster = "6",
  expected_n_cells = 2881L
) {
  provenance <- velocity_read_tsv(path)
  if (!identical(names(provenance), c("key", "value")) ||
      anyNA(provenance$key) || any(!nzchar(provenance$key)) ||
      anyDuplicated(provenance$key)) {
    velocity_stop("Malformed velocity panel provenance: ", path)
  }
  values <- stats::setNames(as.character(provenance$value), provenance$key)
  required <- c(
    "schema_version", "artifact", "table_sha256", "root_cluster",
    "cellcycle_clusters", "n_cells", "cluster_4c_cells",
    "cluster_6_cells", "cluster_10_cells", "scvelo_mode",
    "velocity_basis", "velocity_source", "embedding_source_sha256",
    "embedding_lineage_sha256", "embedding_inventory_sha256",
    "cellcycle_source_sha256", "si_metadata_source_sha256",
    "builder_sha256", "panel_logic_sha256", "scvelo_generator_sha256"
  )
  if (!all(required %in% names(values))) {
    velocity_stop(
      "Velocity panel provenance is missing: ",
      paste(setdiff(required, names(values)), collapse = ", ")
    )
  }
  if (!identical(values[["schema_version"]], "1") ||
      !identical(values[["artifact"]], "figure7_cellcycle_velocity_umap") ||
      !identical(values[["table_sha256"]], velocity_sha256(table_path)) ||
      !identical(values[["root_cluster"]], root_cluster) ||
      !identical(values[["cellcycle_clusters"]], "4c,6,10") ||
      !identical(values[["n_cells"]], as.character(expected_n_cells)) ||
      !identical(values[["cluster_4c_cells"]], "89") ||
      !identical(values[["cluster_6_cells"]], "686") ||
      !identical(values[["cluster_10_cells"]], "2106") ||
      !identical(values[["scvelo_mode"]], "stochastic") ||
      !identical(values[["velocity_basis"]], "reviewed_seurat_umap") ||
      !identical(
        values[["velocity_source"]],
        "scvelo.tl.velocity_embedding_from_velocity_graph"
      )) {
    velocity_stop("Velocity panel provenance disagrees with the reviewed contract")
  }
  hash_keys <- c(
    "table_sha256", "embedding_source_sha256", "cellcycle_source_sha256",
    "embedding_lineage_sha256", "embedding_inventory_sha256",
    "si_metadata_source_sha256", "builder_sha256", "panel_logic_sha256",
    "scvelo_generator_sha256"
  )
  if (any(!grepl("^[0-9a-f]{64}$", unname(values[hash_keys])))) {
    velocity_stop("Velocity panel provenance contains an invalid SHA-256 value")
  }
  provenance
}

velocity_key_values <- function(path, label) {
  table <- velocity_read_tsv(path)
  if (!identical(names(table), c("key", "value")) || !nrow(table) ||
      anyNA(table$key) || any(!nzchar(table$key)) || anyDuplicated(table$key) ||
      anyNA(table$value) || any(!nzchar(table$value))) {
    velocity_stop("Malformed ", label, ": ", path)
  }
  stats::setNames(as.character(table$value), as.character(table$key))
}

velocity_loom_inventory <- function(loom_root) {
  normalized_root <- normalizePath(loom_root, mustWork = TRUE)
  paths <- sort(normalizePath(
    list.files(
      normalized_root,
      pattern = "[.]loom$",
      recursive = TRUE,
      full.names = TRUE
    ),
    mustWork = TRUE
  ))
  if (length(paths) != 18L || anyDuplicated(paths)) {
    velocity_stop(
      "Velocity lineage requires exactly 18 unique loom files; observed ",
      length(paths)
    )
  }
  prefix <- paste0(normalized_root, .Platform$file.sep)
  relative <- substring(paths, nchar(prefix) + 1L)
  if (any(!startsWith(paths, prefix)) || any(!nzchar(relative)) ||
      anyDuplicated(relative) || any(grepl("(^|/)[.][.]?(/|$)", relative))) {
    velocity_stop("Velocity loom inventory contains an unsafe relative path")
  }
  data.frame(
    relative_path = relative,
    bytes = as.character(file.info(paths)$size),
    sha256 = vapply(paths, velocity_sha256, character(1L)),
    stringsAsFactors = FALSE
  )
}

velocity_write_embedding_lineage <- function(
  embedding_path,
  lineage_path,
  inventory_path,
  seurat_rds,
  loom_root,
  scvelo_generator_path,
  config_path,
  environment_lock_path,
  repo_root,
  root_cluster = "6"
) {
  inventory <- velocity_loom_inventory(loom_root)
  velocity_write_tsv(inventory, inventory_path)
  embedding <- velocity_read_tsv(embedding_path)
  if (nrow(embedding) != 35513L) {
    velocity_stop(
      "Velocity embedding must contain the complete reviewed 35,513-cell universe"
    )
  }
  lineage <- data.frame(
    key = c(
      "schema_version", "artifact", "embedding_sha256", "n_cells",
      "root_cluster", "scvelo_mode", "seurat_rds_locator",
      "seurat_rds_sha256", "loom_count", "loom_inventory_sha256",
      "scvelo_generator_sha256", "figure7_config_sha256",
      "environment_lock_sha256"
    ),
    value = c(
      "1", "figure7_scvelo_velocity_umap", velocity_sha256(embedding_path),
      as.character(nrow(embedding)), root_cluster, "stochastic",
      velocity_portable_locator(seurat_rds, repo_root),
      velocity_sha256(seurat_rds), as.character(nrow(inventory)),
      velocity_sha256(inventory_path), velocity_sha256(scvelo_generator_path),
      velocity_sha256(config_path), velocity_sha256(environment_lock_path)
    ),
    stringsAsFactors = FALSE
  )
  velocity_write_tsv(lineage, lineage_path)
  invisible(list(lineage = lineage_path, inventory = inventory_path))
}

velocity_validate_embedding_lineage <- function(
  embedding_path,
  lineage_path,
  inventory_path,
  seurat_rds,
  loom_root,
  scvelo_generator_path,
  config_path,
  environment_lock_path,
  repo_root,
  root_cluster = "6"
) {
  required_paths <- c(
    embedding_path, lineage_path, inventory_path, seurat_rds,
    scvelo_generator_path, config_path, environment_lock_path
  )
  if (any(!file.exists(required_paths)) || !dir.exists(loom_root)) {
    velocity_stop("Velocity embedding lineage has missing current dependencies")
  }
  lineage <- velocity_key_values(lineage_path, "velocity embedding lineage")
  required_keys <- c(
    "schema_version", "artifact", "embedding_sha256", "n_cells",
    "root_cluster", "scvelo_mode", "seurat_rds_locator",
    "seurat_rds_sha256", "loom_count", "loom_inventory_sha256",
    "scvelo_generator_sha256", "figure7_config_sha256",
    "environment_lock_sha256"
  )
  if (!setequal(names(lineage), required_keys) ||
      !identical(lineage[["schema_version"]], "1") ||
      !identical(lineage[["artifact"]], "figure7_scvelo_velocity_umap") ||
      !identical(lineage[["embedding_sha256"]], velocity_sha256(embedding_path)) ||
      !identical(lineage[["n_cells"]], "35513") ||
      !identical(lineage[["root_cluster"]], root_cluster) ||
      !identical(lineage[["scvelo_mode"]], "stochastic") ||
      !identical(
        lineage[["seurat_rds_locator"]],
        velocity_portable_locator(seurat_rds, repo_root)
      ) ||
      !identical(lineage[["seurat_rds_sha256"]], velocity_sha256(seurat_rds)) ||
      !identical(lineage[["loom_count"]], "18") ||
      !identical(
        lineage[["loom_inventory_sha256"]],
        velocity_sha256(inventory_path)
      ) ||
      !identical(
        lineage[["scvelo_generator_sha256"]],
        velocity_sha256(scvelo_generator_path)
      ) ||
      !identical(
        lineage[["figure7_config_sha256"]],
        velocity_sha256(config_path)
      ) ||
      !identical(
        lineage[["environment_lock_sha256"]],
        velocity_sha256(environment_lock_path)
      )) {
    velocity_stop("Velocity embedding lineage disagrees with current dependencies")
  }
  inventory <- velocity_read_tsv(inventory_path)
  if (!identical(names(inventory), c("relative_path", "bytes", "sha256")) ||
      nrow(inventory) != 18L || anyNA(inventory$relative_path) ||
      any(!nzchar(inventory$relative_path)) || anyDuplicated(inventory$relative_path) ||
      any(grepl("(^|/)[.][.]?(/|$)", inventory$relative_path))) {
    velocity_stop("Malformed 18-file velocity loom inventory")
  }
  current <- velocity_loom_inventory(loom_root)
  inventory[] <- lapply(inventory, as.character)
  current[] <- lapply(current, as.character)
  current <- current[order(current$relative_path), , drop = FALSE]
  inventory <- inventory[order(inventory$relative_path), , drop = FALSE]
  rownames(current) <- NULL
  rownames(inventory) <- NULL
  if (!identical(current, inventory)) {
    velocity_stop("Velocity loom inventory differs from the current 18 loom files")
  }
  invisible(lineage)
}

velocity_generated_bundle_matches_dependencies <- function(
  table_path,
  provenance_path,
  embedding_path,
  embedding_lineage_path,
  embedding_inventory_path,
  cellcycle_path,
  si_metadata_path,
  builder_path,
  panel_logic_path,
  scvelo_generator_path,
  seurat_rds,
  loom_root,
  config_path,
  environment_lock_path,
  repo_root
) {
  table <- velocity_validate_panel_table(velocity_read_tsv(table_path))
  provenance <- velocity_validate_provenance(provenance_path, table_path)
  values <- stats::setNames(
    as.character(provenance$value),
    as.character(provenance$key)
  )
  expected_hashes <- c(
    embedding_source_sha256 = velocity_sha256(embedding_path),
    embedding_lineage_sha256 = velocity_sha256(embedding_lineage_path),
    embedding_inventory_sha256 = velocity_sha256(embedding_inventory_path),
    cellcycle_source_sha256 = velocity_sha256(cellcycle_path),
    si_metadata_source_sha256 = velocity_sha256(si_metadata_path),
    builder_sha256 = velocity_sha256(builder_path),
    panel_logic_sha256 = velocity_sha256(panel_logic_path),
    scvelo_generator_sha256 = velocity_sha256(scvelo_generator_path)
  )
  if (any(values[names(expected_hashes)] != expected_hashes)) {
    velocity_stop("Generated velocity panel cache has stale input or code hashes")
  }
  velocity_validate_embedding_lineage(
    embedding_path,
    embedding_lineage_path,
    embedding_inventory_path,
    seurat_rds,
    loom_root,
    scvelo_generator_path,
    config_path,
    environment_lock_path,
    repo_root
  )
  invisible(table)
}

velocity_build_panel_table <- function(
  embedding_path,
  cellcycle_path,
  si_metadata_path,
  root_cluster = "6",
  pseudotime_tolerance = 1e-7,
  umap_tolerance = 1e-7
) {
  embedding <- velocity_read_tsv(embedding_path)
  expected_embedding <- c(
    "cell", "UMAP_1", "UMAP_2", "velocity_UMAP_1",
    "velocity_UMAP_2", "velocity_pseudotime", "cluster",
    "context", "sample", "is_root"
  )
  if (!identical(names(embedding), expected_embedding)) {
    velocity_stop(
      "Full scVelo embedding columns must be exactly: ",
      paste(expected_embedding, collapse = ", ")
    )
  }
  cellcycle <- velocity_read_csv(cellcycle_path)
  required_cellcycle <- c(
    "cell_id", "sample_id", "cluster", "pseudotime"
  )
  if (!all(required_cellcycle %in% names(cellcycle))) {
    velocity_stop("CellCycle input lacks cell/sample/cluster/pseudotime fields")
  }
  si_metadata <- velocity_read_csv(si_metadata_path)
  required_si <- c("cell_id", "UMAP_1", "UMAP_2", "cluster_id", "context")
  if (!all(required_si %in% names(si_metadata))) {
    velocity_stop("SI metadata lacks reviewed UMAP/cluster/context fields")
  }
  if (anyDuplicated(embedding$cell) || anyDuplicated(cellcycle$cell_id) ||
      anyDuplicated(si_metadata$cell_id)) {
    velocity_stop("Velocity source inputs contain duplicated cell IDs")
  }
  if (nrow(cellcycle) != 2881L) {
    velocity_stop("CellCycle input must contain exactly 2,881 tumor cells")
  }
  embedding_index <- match(cellcycle$cell_id, embedding$cell)
  si_index <- match(cellcycle$cell_id, si_metadata$cell_id)
  if (anyNA(embedding_index) || anyNA(si_index)) {
    velocity_stop("Every CellCycle cell must match the scVelo and SI metadata tables")
  }
  selected <- embedding[embedding_index, , drop = FALSE]
  selected_si <- si_metadata[si_index, , drop = FALSE]
  if (any(as.character(selected$context) != "Tumor") ||
      any(as.character(selected_si$context) != "Tumor") ||
      !identical(as.character(selected$cluster), as.character(cellcycle$cluster)) ||
      !identical(as.character(selected_si$cluster_id), as.character(cellcycle$cluster)) ||
      !identical(as.character(selected$sample), as.character(cellcycle$sample_id))) {
    velocity_stop(
      "CellCycle membership, context, cluster, or sample differs across source tables"
    )
  }
  source_pseudotime <- suppressWarnings(as.numeric(selected$velocity_pseudotime))
  reviewed_pseudotime <- suppressWarnings(as.numeric(cellcycle$pseudotime))
  if (any(!is.finite(source_pseudotime)) ||
      any(!is.finite(reviewed_pseudotime)) ||
      max(abs(source_pseudotime - reviewed_pseudotime)) > pseudotime_tolerance) {
    velocity_stop(
      "scVelo embedding pseudotime differs from the reviewed CellCycle values"
    )
  }
  embedding_umap <- cbind(
    suppressWarnings(as.numeric(selected$UMAP_1)),
    suppressWarnings(as.numeric(selected$UMAP_2))
  )
  reviewed_umap <- cbind(
    suppressWarnings(as.numeric(selected_si$UMAP_1)),
    suppressWarnings(as.numeric(selected_si$UMAP_2))
  )
  if (any(!is.finite(embedding_umap)) || any(!is.finite(reviewed_umap)) ||
      max(abs(embedding_umap - reviewed_umap)) > umap_tolerance) {
    velocity_stop(
      "scVelo embedding coordinates differ from the reviewed Seurat UMAP"
    )
  }
  table <- data.frame(
    cell_id = as.character(cellcycle$cell_id),
    sample_id = as.character(cellcycle$sample_id),
    cluster = as.character(cellcycle$cluster),
    UMAP_1 = embedding_umap[, 1L],
    UMAP_2 = embedding_umap[, 2L],
    velocity_UMAP_1 = suppressWarnings(as.numeric(selected$velocity_UMAP_1)),
    velocity_UMAP_2 = suppressWarnings(as.numeric(selected$velocity_UMAP_2)),
    velocity_pseudotime = reviewed_pseudotime,
    is_root = velocity_parse_logical(selected$is_root, "source is_root value"),
    stringsAsFactors = FALSE
  )
  table <- table[order(table$cell_id, method = "radix"), , drop = FALSE]
  rownames(table) <- NULL
  velocity_validate_panel_table(table, root_cluster = root_cluster)
}

velocity_write_panel_bundle <- function(
  table,
  table_path,
  provenance_path,
  embedding_path,
  embedding_lineage_path,
  embedding_inventory_path,
  cellcycle_path,
  si_metadata_path,
  builder_path,
  panel_logic_path,
  scvelo_generator_path,
  root_cluster = "6"
) {
  table <- velocity_validate_panel_table(table, root_cluster = root_cluster)
  velocity_write_tsv(table, table_path)
  provenance <- data.frame(
    key = c(
      "schema_version", "artifact", "table_sha256", "root_cluster",
      "cellcycle_clusters", "n_cells", "cluster_4c_cells",
      "cluster_6_cells", "cluster_10_cells", "scvelo_mode",
      "velocity_basis", "velocity_source", "embedding_source_sha256",
      "embedding_lineage_sha256", "embedding_inventory_sha256",
      "cellcycle_source_sha256", "si_metadata_source_sha256",
      "builder_sha256", "panel_logic_sha256", "scvelo_generator_sha256"
    ),
    value = c(
      "1", "figure7_cellcycle_velocity_umap", velocity_sha256(table_path),
      root_cluster, "4c,6,10", as.character(nrow(table)),
      as.character(sum(table$cluster == "4c")),
      as.character(sum(table$cluster == "6")),
      as.character(sum(table$cluster == "10")), "stochastic",
      "reviewed_seurat_umap",
      "scvelo.tl.velocity_embedding_from_velocity_graph",
      velocity_sha256(embedding_path), velocity_sha256(embedding_lineage_path),
      velocity_sha256(embedding_inventory_path), velocity_sha256(cellcycle_path),
      velocity_sha256(si_metadata_path), velocity_sha256(builder_path),
      velocity_sha256(panel_logic_path), velocity_sha256(scvelo_generator_path)
    ),
    stringsAsFactors = FALSE
  )
  velocity_write_tsv(provenance, provenance_path)
  velocity_validate_provenance(
    provenance_path,
    table_path,
    root_cluster = root_cluster,
    expected_n_cells = nrow(table)
  )
  invisible(list(table = table_path, provenance = provenance_path))
}

velocity_build_grid <- function(
  table,
  nx = 24L,
  ny = 24L,
  min_bin_n = 4L
) {
  table <- velocity_validate_panel_table(table)
  vector_length <- sqrt(
    table$velocity_UMAP_1^2 + table$velocity_UMAP_2^2
  )
  cap <- stats::quantile(
    vector_length,
    probs = 0.995,
    na.rm = TRUE,
    names = FALSE
  )
  keep <- vector_length > 0 & vector_length <= cap
  vectors <- table[keep, , drop = FALSE]
  x_range <- range(table$UMAP_1, finite = TRUE)
  y_range <- range(table$UMAP_2, finite = TRUE)
  if (diff(x_range) <= 0 || diff(y_range) <= 0) {
    velocity_stop("Velocity UMAP must span both dimensions")
  }
  x_breaks <- seq(x_range[[1L]], x_range[[2L]], length.out = nx + 1L)
  y_breaks <- seq(y_range[[1L]], y_range[[2L]], length.out = ny + 1L)
  vectors$grid_x <- findInterval(vectors$UMAP_1, x_breaks, all.inside = TRUE)
  vectors$grid_y <- findInterval(vectors$UMAP_2, y_breaks, all.inside = TRUE)
  split_key <- interaction(vectors$grid_x, vectors$grid_y, drop = TRUE)
  groups <- split(seq_len(nrow(vectors)), split_key)
  rows <- lapply(groups, function(index) {
    if (length(index) < min_bin_n) return(NULL)
    data.frame(
      UMAP_1 = mean(vectors$UMAP_1[index]),
      UMAP_2 = mean(vectors$UMAP_2[index]),
      velocity_UMAP_1 = mean(vectors$velocity_UMAP_1[index]),
      velocity_UMAP_2 = mean(vectors$velocity_UMAP_2[index]),
      n_cells = length(index),
      stringsAsFactors = FALSE
    )
  })
  grid <- do.call(rbind, rows[!vapply(rows, is.null, logical(1L))])
  if (is.null(grid) || nrow(grid) < 12L) {
    velocity_stop("Too few populated bins to display an RNA-velocity field")
  }
  grid$vector_length <- sqrt(
    grid$velocity_UMAP_1^2 + grid$velocity_UMAP_2^2
  )
  grid <- grid[is.finite(grid$vector_length) & grid$vector_length > 0, , drop = FALSE]
  reference <- stats::quantile(
    grid$vector_length,
    probs = 0.90,
    na.rm = TRUE,
    names = FALSE
  )
  step <- min(diff(x_range) / nx, diff(y_range) / ny)
  if (!is.finite(reference) || reference <= 0 || !is.finite(step) || step <= 0) {
    velocity_stop("Cannot scale the RNA-velocity field")
  }
  scale <- 0.80 * step / reference
  grid$UMAP_1_to <- grid$UMAP_1 + grid$velocity_UMAP_1 * scale
  grid$UMAP_2_to <- grid$UMAP_2 + grid$velocity_UMAP_2 * scale
  grid[order(grid$UMAP_1, grid$UMAP_2), , drop = FALSE]
}

velocity_cluster_labels <- function(table) {
  rows <- lapply(split(table, table$cluster), function(cluster) {
    data.frame(
      cluster = cluster$cluster[[1L]],
      UMAP_1 = stats::median(cluster$UMAP_1),
      UMAP_2 = stats::median(cluster$UMAP_2),
      stringsAsFactors = FALSE
    )
  })
  labels <- do.call(rbind, rows)
  labels$label <- ifelse(
    labels$cluster == "6",
    "6  ROOT",
    labels$cluster
  )
  rownames(labels) <- NULL
  labels
}

velocity_root_hull <- function(table, fraction = 0.90) {
  root <- table[table$is_root, , drop = FALSE]
  center <- c(stats::median(root$UMAP_1), stats::median(root$UMAP_2))
  distance <- sqrt(
    (root$UMAP_1 - center[[1L]])^2 +
      (root$UMAP_2 - center[[2L]])^2
  )
  cutoff <- stats::quantile(
    distance,
    probs = fraction,
    na.rm = TRUE,
    names = FALSE
  )
  core <- root[distance <= cutoff, , drop = FALSE]
  if (nrow(core) < 3L) velocity_stop("Cluster 6 root hull has fewer than three cells")
  hull <- grDevices::chull(core$UMAP_1, core$UMAP_2)
  core[c(hull, hull[[1L]]), c("UMAP_1", "UMAP_2"), drop = FALSE]
}

velocity_render_panel <- function(
  table_path,
  provenance_path,
  output_dir,
  width = 6.8,
  height = 5.8,
  dpi = 600L
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    velocity_stop("R package 'ggplot2' is required")
  }
  table <- velocity_validate_panel_table(velocity_read_tsv(table_path))
  velocity_validate_provenance(provenance_path, table_path)
  grid <- velocity_build_grid(table)
  labels <- velocity_cluster_labels(table)
  root_hull <- velocity_root_hull(table)
  dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "metadata"), recursive = TRUE, showWarnings = FALSE)
  plot_table_path <- file.path(
    output_dir,
    "tables",
    "panel_velocity_pseudotime_plot_data.tsv"
  )
  if (!file.copy(table_path, plot_table_path, overwrite = TRUE)) {
    velocity_stop("Could not preserve the exact velocity plot-facing table")
  }
  if (!identical(velocity_sha256(plot_table_path), velocity_sha256(table_path))) {
    velocity_stop("Rendered velocity plot table is not byte-identical to its source")
  }
  velocity_write_tsv(
    grid,
    file.path(output_dir, "tables", "panel_velocity_pseudotime_vector_grid.tsv")
  )
  plot <- ggplot2::ggplot(
    table,
    ggplot2::aes(x = UMAP_1, y = UMAP_2, color = velocity_pseudotime)
  ) +
    ggplot2::geom_point(size = 0.48, alpha = 0.78) +
    ggplot2::geom_path(
      data = root_hull,
      ggplot2::aes(x = UMAP_1, y = UMAP_2),
      inherit.aes = FALSE,
      color = "#D73027",
      linewidth = 0.75,
      linetype = "22"
    ) +
    ggplot2::geom_segment(
      data = grid,
      ggplot2::aes(
        x = UMAP_1,
        y = UMAP_2,
        xend = UMAP_1_to,
        yend = UMAP_2_to
      ),
      inherit.aes = FALSE,
      linewidth = 0.30,
      color = "black",
      alpha = 0.82,
      arrow = grid::arrow(
        length = grid::unit(0.045, "inches"),
        type = "closed"
      )
    ) +
    ggplot2::geom_label(
      data = labels[labels$cluster != "6", , drop = FALSE],
      ggplot2::aes(x = UMAP_1, y = UMAP_2, label = label),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold",
      linewidth = 0.25,
      label.padding = grid::unit(0.14, "lines"),
      fill = "white",
      color = "black",
      alpha = 0.90
    ) +
    ggplot2::geom_label(
      data = labels[labels$cluster == "6", , drop = FALSE],
      ggplot2::aes(x = UMAP_1, y = UMAP_2, label = label),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold",
      linewidth = 0.25,
      label.padding = grid::unit(0.14, "lines"),
      fill = "white",
      color = "#B2182B",
      alpha = 0.90
    ) +
    ggplot2::scale_color_gradientn(
      colors = c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"),
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      name = "Velocity\npseudotime"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Stochastic RNA velocity across the tumor CellCycle trajectory",
      subtitle = "2,881 cells; dashed outline marks cluster 6, the pre-specified root set",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 9),
      legend.position = "right",
      aspect.ratio = 1
    )
  pdf_path <- file.path(
    output_dir,
    "figures",
    "panel_SuppFig10_velocity_pseudotime.pdf"
  )
  png_path <- file.path(
    output_dir,
    "figures",
    "panel_SuppFig10_velocity_pseudotime.png"
  )
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, units = "in")
  ggplot2::ggsave(
    png_path,
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  source_provenance <- velocity_read_tsv(provenance_path)
  run_provenance <- rbind(
    source_provenance,
    data.frame(
      key = c(
        "rendered_pdf_sha256", "rendered_png_sha256",
        "plot_table_sha256", "vector_grid_sha256", "vector_grid_rows",
        "root_display_policy"
      ),
      value = c(
        velocity_sha256(pdf_path),
        velocity_sha256(png_path),
        velocity_sha256(plot_table_path),
        velocity_sha256(file.path(
          output_dir, "tables", "panel_velocity_pseudotime_vector_grid.tsv"
        )),
        as.character(nrow(grid)),
        "cluster_6_label_plus_90_percent_core_hull"
      ),
      stringsAsFactors = FALSE
    )
  )
  velocity_write_tsv(
    run_provenance,
    file.path(output_dir, "metadata", "velocity_pseudotime_provenance.tsv")
  )
  invisible(list(pdf = pdf_path, png = png_path, plot = plot, grid = grid))
}
