# Reviewed SI4/SI7 panels promoted into the assembled main Figure 7.
#
# The plotting implementation lives in SI_figures/shared_context_panels.R so
# the main and supplementary copies cannot silently diverge.  This layer owns
# only the reviewed-cache contract and Figure 7-specific A-L mapping.

figure7_read_delimited <- function(path, delimiter, label) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    figure7_stop("Missing or empty ", label, ": ", path)
  }
  reader <- if (identical(delimiter, "\t")) utils::read.delim else utils::read.csv
  data <- reader(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = if (identical(delimiter, "\t")) "" else "\"",
    comment.char = "",
    na.strings = c("", "NA", "NaN")
  )
  if (!nrow(data)) figure7_stop(label, " has no data rows")
  data
}

figure7_context_key_values <- function(path, label) {
  data <- figure7_read_delimited(path, "\t", label)
  if (!identical(names(data), c("key", "value")) ||
      anyDuplicated(data$key) ||
      anyNA(data$key) ||
      any(!nzchar(data$key))) {
    figure7_stop(label, " must contain unique key/value rows")
  }
  stats::setNames(as.character(data$value), as.character(data$key))
}

figure7_validate_si_cache <- function(
  cache_dir,
  repo_root,
  config,
  policy = "reviewed",
  upstream_input_manifest = "",
  source_run_config = "",
  source_provenance = ""
) {
  if (!policy %in% c("reviewed", "generated-human-only")) {
    figure7_stop("Unknown Figure 7 SI context-cache policy: ", policy)
  }
  cache_label <- if (identical(policy, "reviewed")) {
    "Reviewed SI cache"
  } else {
    "Generated SI context cache"
  }
  cache_dir <- normalizePath(cache_dir, mustWork = TRUE)
  manifest_path <- file.path(cache_dir, "manifest.tsv")
  if (!file.exists(manifest_path)) {
    figure7_stop(cache_label, " is missing manifest.tsv: ", cache_dir)
  }
  manifest_hash <- figure7_sha256(manifest_path)
  if (identical(policy, "reviewed")) {
    figure7_verify_checksum(
      manifest_path,
      as.character(config$si_figures$reviewed_manifest_sha256),
      "reviewed SI Figures cache manifest"
    )
  }
  manifest <- figure7_read_delimited(
    manifest_path,
    "\t",
    "SI Figures context-cache manifest"
  )
  if (!identical(names(manifest), c(
        "filename", "bytes", "sha256", "source_revision", "notes"
      )) ||
      nrow(manifest) != as.integer(config$si_figures$cache_file_count) ||
      anyDuplicated(manifest$filename)) {
    figure7_stop("SI context-cache manifest contract is invalid")
  }
  paths <- file.path(cache_dir, manifest$filename)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    figure7_stop(
      cache_label, " is incomplete: ",
      paste(basename(missing), collapse = ", ")
    )
  }
  observed <- vapply(paths, figure7_sha256, character(1L))
  if (!identical(unname(observed), as.character(manifest$sha256))) {
    bad <- manifest$filename[observed != as.character(manifest$sha256)]
    figure7_stop(
      cache_label, " checksum mismatch: ",
      paste(bad, collapse = ", ")
    )
  }
  validator <- file.path(
    repo_root,
    "Code",
    "tools",
    "validate_si_figures_table_cache.py"
  )
  validator_args <- c(shQuote(validator), "--cache-dir", shQuote(cache_dir))
  if (identical(policy, "generated-human-only")) {
    validator_args <- c(
      validator_args,
      "--si7-policy",
      "generated-human-only"
    )
  }
  status <- system2(
    "python3",
    validator_args,
    stdout = FALSE,
    stderr = FALSE
  )
  if (!identical(status, 0L)) {
    figure7_stop(
      "SI context cache failed the strict ",
      policy,
      " table validator"
    )
  }

  descriptor <- list(
    root = cache_dir,
    manifest = manifest,
    manifest_path = manifest_path,
    manifest_sha256 = manifest_hash,
    policy = policy,
    kind = if (identical(policy, "reviewed")) {
      "reviewed_human_only_frozen"
    } else {
      "generated_human_only_run_scoped"
    },
    canonical_publication_allowed = identical(policy, "reviewed"),
    upstream_input_manifest = "",
    upstream_input_manifest_sha256 = "not_applicable",
    source_run_config = "",
    source_run_config_sha256 = "not_applicable",
    source_provenance = "",
    source_provenance_sha256 = "not_applicable"
  )

  if (identical(policy, "reviewed")) {
    if (any(nzchar(c(
          upstream_input_manifest,
          source_run_config,
          source_provenance
        )))) {
      figure7_stop(
        "Reviewed SI context cache cannot carry generated-run lineage arguments"
      )
    }
    return(descriptor)
  }

  lineage_paths <- c(
    upstream_input_manifest = upstream_input_manifest,
    source_run_config = source_run_config,
    source_provenance = source_provenance
  )
  if (any(!nzchar(lineage_paths))) {
    figure7_stop(
      "Generated SI context cache requires its analysis manifest, run config, ",
      "and provenance"
    )
  }
  lineage_paths <- vapply(
    lineage_paths,
    normalizePath,
    character(1L),
    mustWork = TRUE
  )
  upstream <- figure7_read_delimited(
    lineage_paths[["upstream_input_manifest"]],
    "\t",
    "generated SI analysis input manifest"
  )
  if (!identical(
        names(upstream),
        c("role", "repo_relative_path", "sha256", "bytes")
      ) ||
      anyNA(upstream$role) ||
      any(!nzchar(upstream$role)) ||
      any(!grepl("^[0-9a-f]{64}$", upstream$sha256)) ||
      any(!is.finite(suppressWarnings(as.numeric(upstream$bytes)))) ||
      any(grepl("^/", upstream$repo_relative_path))) {
    figure7_stop("Generated SI analysis input manifest is malformed or nonportable")
  }
  role_count <- table(upstream$role)
  required_single_roles <- c(
    "si_figures_generated_cache_manifest",
    "raw_build_input_manifest",
    "raw_build_si_raw_scientific_code_contract",
    "raw_build_si_raw_config_contract",
    "raw_build_si_raw_runtime_contract",
    "raw_build_versioned_endpoint_ploidy"
  )
  single_counts <- unname(role_count[required_single_roles])
  generated_table_count <- unname(role_count["si_figures_generated_table"])
  if (anyNA(single_counts) || any(single_counts != 1L) ||
      is.na(generated_table_count) ||
      generated_table_count != nrow(manifest) ||
      !any(startsWith(upstream$role, "raw_build_source_"))) {
    figure7_stop("Generated SI analysis lineage is incomplete")
  }
  manifest_row <- upstream[
    upstream$role == "si_figures_generated_cache_manifest",
    ,
    drop = FALSE
  ]
  if (!identical(manifest_row$sha256[[1L]], manifest_hash) ||
      as.numeric(manifest_row$bytes[[1L]]) != file.info(manifest_path)$size) {
    figure7_stop("Generated SI analysis manifest does not bind its cache manifest")
  }
  table_rows <- upstream[
    upstream$role == "si_figures_generated_table",
    ,
    drop = FALSE
  ]
  table_names <- sub(
    "^external:",
    "",
    basename(as.character(table_rows$repo_relative_path))
  )
  if (anyDuplicated(table_names) || !setequal(table_names, manifest$filename)) {
    figure7_stop("Generated SI analysis manifest has the wrong table inventory")
  }
  table_rows <- table_rows[match(manifest$filename, table_names), , drop = FALSE]
  if (!identical(as.character(table_rows$sha256), as.character(manifest$sha256)) ||
      !identical(
        as.numeric(table_rows$bytes),
        as.numeric(manifest$bytes)
      )) {
    figure7_stop("Generated SI analysis manifest does not bind exact cache bytes")
  }

  run_values <- figure7_context_key_values(
    lineage_paths[["source_run_config"]],
    "generated SI source run config"
  )
  source_revision <- paste(unique(manifest$source_revision), collapse = ";")
  if (!identical(
        run_values[["table_mode"]],
        "run_scoped_generated_human_only_plot_tables"
      ) ||
      !identical(run_values[["si7_canonical_publication_allowed"]], "false") ||
      !identical(run_values[["cache_file_count"]], as.character(nrow(manifest))) ||
      !identical(run_values[["generated_table_source_revision"]], source_revision)) {
    figure7_stop("Generated SI source run config contradicts its cache")
  }
  provenance_values <- figure7_context_key_values(
    lineage_paths[["source_provenance"]],
    "generated SI source provenance"
  )
  renderer <- file.path(
    repo_root,
    "Code",
    "in-vivo",
    "SI_figures",
    "generate_supplementary_figures.R"
  )
  composition_helper <- file.path(
    repo_root,
    "Code",
    "in-vivo",
    "SI_figures",
    "normalized_composition.R"
  )
  shared_helper <- file.path(
    repo_root,
    "Code",
    "in-vivo",
    "SI_figures",
    "shared_context_panels.R"
  )
  if (!identical(provenance_values[["entrypoint_sha256"]], figure7_sha256(renderer)) ||
      !identical(
        provenance_values[["normalized_composition_helper_sha256"]],
        figure7_sha256(composition_helper)
      ) ||
      !identical(
        provenance_values[["shared_context_panels_helper_sha256"]],
        figure7_sha256(shared_helper)
      ) ||
      !identical(provenance_values[["table_cache_manifest_sha256"]], manifest_hash) ||
      !identical(provenance_values[["si7_canonical_publication_allowed"]], "false")) {
    figure7_stop("Generated SI source provenance is stale or contradictory")
  }

  descriptor$upstream_input_manifest <-
    lineage_paths[["upstream_input_manifest"]]
  descriptor$upstream_input_manifest_sha256 <-
    figure7_sha256(descriptor$upstream_input_manifest)
  descriptor$source_run_config <- lineage_paths[["source_run_config"]]
  descriptor$source_run_config_sha256 <-
    figure7_sha256(descriptor$source_run_config)
  descriptor$source_provenance <- lineage_paths[["source_provenance"]]
  descriptor$source_provenance_sha256 <-
    figure7_sha256(descriptor$source_provenance)
  descriptor
}

figure7_validate_reviewed_si_cache <- function(cache_dir, repo_root, config) {
  figure7_validate_si_cache(cache_dir, repo_root, config, policy = "reviewed")
}

figure7_context_matrix <- function(cache_dir, filename, cluster_levels, label) {
  data <- figure7_read_delimited(file.path(cache_dir, filename), "\t", label)
  if (!identical(names(data), c("pathway", cluster_levels))) {
    figure7_stop(label, " does not match the reviewed cluster order")
  }
  matrix_data <- as.matrix(data[, cluster_levels, drop = FALSE])
  storage.mode(matrix_data) <- "double"
  rownames(matrix_data) <- data$pathway
  if (!identical(dim(matrix_data), c(20L, length(cluster_levels))) ||
      any(!is.finite(matrix_data))) {
    figure7_stop(label, " must be a finite 20 x 9 matrix")
  }
  matrix_data
}

figure7_build_context_panels <- function(
  cache_dir,
  repo_root,
  config,
  output_dir = NULL,
  policy = "reviewed",
  upstream_input_manifest = "",
  source_run_config = "",
  source_provenance = ""
) {
  cache <- figure7_validate_si_cache(
    cache_dir,
    repo_root,
    config,
    policy,
    upstream_input_manifest,
    source_run_config,
    source_provenance
  )
  si <- config$si_figures
  cluster_levels <- as.character(unlist(si$cluster_order))
  ploidy_levels <- as.character(unlist(si$ploidy_levels))
  dose_levels <- as.character(unlist(si$dose_levels))
  plot_seed <- as.integer(si$plot_shuffle_seed)

  cells <- figure7_read_delimited(
    file.path(cache$root, "si_figures_cell_metadata.csv"),
    ",",
    "SI Figures cell metadata"
  )
  cluster_key <- figure7_read_delimited(
    file.path(cache$root, "si_figures_cluster_key.tsv"),
    "\t",
    "SI Figures cluster key"
  )
  cells$UMAP_1 <- figure7_numeric(cells$UMAP_1)
  cells$UMAP_2 <- figure7_numeric(cells$UMAP_2)
  cells$s_phase_score <- figure7_numeric(cells$s_phase_score)
  cells$endpoint_ploidy <- suppressWarnings(figure7_numeric(cells$endpoint_ploidy))
  cells$included_in_si_figures <- tolower(
    as.character(cells$included_in_si_figures)
  ) %in% c("true", "t", "1")
  cell_data <- data.frame(
    cell = cells$cell_id,
    UMAP_1 = cells$UMAP_1,
    UMAP_2 = cells$UMAP_2,
    mouse = cells$sample_id,
    cluster = factor(cells$cluster_id, levels = cluster_levels),
    initial_ploidy = factor(cells$initial_ploidy, levels = ploidy_levels),
    dose = factor(cells$dose, levels = dose_levels),
    s_phase_score = cells$s_phase_score,
    endpoint_ploidy = cells$endpoint_ploidy,
    context = factor(cells$context, levels = c("Tumor", "CellLine")),
    included = cells$included_in_si_figures,
    stringsAsFactors = FALSE
  )
  n_tumor <- sum(cell_data$included)
  point_size <- if (n_tumor > 50000L) {
    0.08
  } else if (n_tumor > 20000L) {
    0.14
  } else {
    0.24
  }
  cluster_colors <- stats::setNames(
    as.character(cluster_key$color),
    as.character(cluster_key$cluster_id)
  )
  ploidy_colors <- c("2N" = "#4C78A8", "4N" = "#E45756")
  context_colors <- c("Tumor" = "#4C78A8", "CellLine" = "#F2CF5B")

  si4 <- shared_context_build_si4_panels(
    data = cell_data,
    cluster_levels = cluster_levels,
    cluster_colors = cluster_colors,
    ploidy_colors = ploidy_colors,
    context_colors = context_colors,
    point_size = point_size,
    plot_seed = plot_seed,
    tags = c(
      cluster = "",
      initial_ploidy = "",
      context = "",
      composition = ""
    ),
    repel_cluster_labels = TRUE
  )
  significance_layers <- which(vapply(
    si4$plots$composition$layers,
    function(layer) {
      label_mapping <- rlang::as_label(layer$mapping$label)
      identical(label_mapping, "significance")
    },
    logical(1L)
  ))
  if (length(significance_layers) != 1L) {
    figure7_stop("Main Figure 7 composition plot lacks one significance layer")
  }
  star_layer <- si4$plots$composition$layers[[significance_layers[[1L]]]]
  star_layer$aes_params$size <- 3.8
  star_layer$aes_params$fontface <- "bold"
  si4$plots$composition$layers[[significance_layers[[1L]]]] <- star_layer
  gsea_matrix <- figure7_context_matrix(
    cache$root,
    "si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv",
    cluster_levels,
    "SI Figure 7 GSEA matrix"
  )
  gsea <- shared_context_heatmap_plot(
    gsea_matrix,
    "",
    TRUE,
    "",
    fontsize_row = 7.25,
    fontsize_col = 7,
    reader_labels = TRUE,
    mark_zero_missing = FALSE,
    treeheight_row = 20,
    treeheight_col = 20
  )

  if (!is.null(output_dir)) {
    figure7_write_tsv(
      si4$composition_result$plot_data,
      file.path(output_dir, "tables", "main_composite_panel_F_plot_data.tsv")
    )
    figure7_write_tsv(
      si4$composition_result$tests,
      file.path(
        output_dir,
        "tables",
        "main_composite_panel_F_enrichment_tests.tsv"
      )
    )
  }

  list(
    plots = list(
      C = si4$plots$cluster,
      D = si4$plots$initial_ploidy,
      E = si4$plots$context,
      F = si4$plots$composition,
      G = gsea
    ),
    composition_result = si4$composition_result,
    gsea_matrix = gsea_matrix,
    cache_manifest_sha256 = cache$manifest_sha256,
    cache_manifest_path = cache$manifest_path,
    cache_policy = cache$policy,
    cache_kind = cache$kind,
    cache_canonical_publication_allowed =
      cache$canonical_publication_allowed,
    upstream_input_manifest = cache$upstream_input_manifest,
    upstream_input_manifest_sha256 =
      cache$upstream_input_manifest_sha256,
    source_run_config = cache$source_run_config,
    source_run_config_sha256 = cache$source_run_config_sha256,
    source_provenance = cache$source_provenance,
    source_provenance_sha256 = cache$source_provenance_sha256,
    composite_filename = figure7_main_composite_filename(config, cache$policy)
  )
}

figure7_main_composite_plots <- function(
  source_plots,
  context_plots,
  interval_treatment_plot,
  copy_number_plot,
  config
) {
  if (!identical(names(source_plots), LETTERS[1:5]) ||
      !identical(names(context_plots), LETTERS[3:7]) ||
      !inherits(interval_treatment_plot, "ggplot") ||
      !inherits(copy_number_plot, c("ggplot", "patchwork", "wrapped_patch"))) {
    figure7_stop("Cannot assemble Figure 7 from an incomplete plot set")
  }
  expected_mapping <- c(
    A = "7A", B = "7C", C = "SI4A", D = "SI4B", E = "SI4C",
    F = "SI4E", G = "SI7B", H = "7B", I = "7I", J = "7J",
    K = "7D", L = "7E"
  )
  configured_mapping <- unlist(
    config$panels$main_composite$panel_order,
    use.names = TRUE
  )
  if (!identical(configured_mapping, expected_mapping)) {
    figure7_stop("Main Figure 7 panel mapping is not the reviewed A-L order")
  }

  panel_b_components <- attr(
    source_plots$B, "figure7_panel_b_components", exact = TRUE
  )
  if (is.null(panel_b_components) ||
      !inherits(panel_b_components$ecdf, "ggplot") ||
      !inherits(panel_b_components$localization, "ggplot")) {
    figure7_stop(
      "Figure 7 source panel 7B must expose separate ECDF and localization plots"
    )
  }

  list(
    A = source_plots$A,
    B = source_plots$C,
    C = context_plots$C,
    D = context_plots$D,
    E = context_plots$E,
    F = context_plots$F,
    G = context_plots$G,
    H = panel_b_components$ecdf,
    I = interval_treatment_plot,
    J = copy_number_plot,
    K = source_plots$D,
    L = source_plots$E
  )
}
