testthat::test_that("main Figure 7 config and both render paths use exact first-citation A-K order", {
  config_path <- file.path(module_dir, "figure7_config.yaml")
  config <- figure7_test_inputs()$config
  expected_mapping <- c(
    A = "7A", B = "7C", C = "SI4A", D = "SI4B", E = "SI4C",
    F = "SI4E", G = "SI7B", H = "7B", I = "7F", J = "7D", K = "7E"
  )
  observed_mapping <- unlist(
    config$panels$main_composite$panel_order,
    use.names = TRUE
  )

  testthat::expect_identical(
    as.character(observed_mapping),
    unname(expected_mapping)
  )
  testthat::expect_identical(names(observed_mapping), names(expected_mapping))
  testthat::expect_identical(
    as.character(config$panels$main_composite$filename),
    "Figure7_reviewed_GRCh.png"
  )
  testthat::expect_identical(
    as.character(config$panels$main_composite$generated_candidate_filename),
    "Figure7_generated_GRCh_candidate.png"
  )
  testthat::expect_identical(
    figure7_panel_contract(config)$panel_id,
    c(figure7_panel_ids(TRUE), "7A-7K_composite")
  )
  testthat::expect_identical(
    tail(figure7_panel_contract(config)$filename, 1L),
    "Figure7_reviewed_GRCh.png"
  )

  tampered_config <- yaml::read_yaml(config_path)
  tampered_config$panels$main_composite$panel_order$F <- "SI7B"
  tampered_config$panels$main_composite$panel_order$G <- "SI4E"
  tampered_path <- tempfile(fileext = ".yaml")
  yaml::write_yaml(tampered_config, tampered_path)
  testthat::expect_error(
    figure7_read_config(tampered_path),
    "reviewed A-K first-citation panel mapping"
  )

  find_calls <- function(node, function_name) {
    found <- list()
    visit <- function(value) {
      if (is.call(value)) {
        if (identical(value[[1L]], as.name(function_name))) {
          found[[length(found) + 1L]] <<- value
        }
        parts <- as.list(value)
        for (index in seq_along(parts)) {
          if (!rlang::is_missing(parts[[index]])) visit(parts[[index]])
        }
      } else if (is.expression(value) || is.pairlist(value)) {
        parts <- as.list(value)
        for (index in seq_along(parts)) {
          if (!rlang::is_missing(parts[[index]])) visit(parts[[index]])
        }
      }
      invisible(NULL)
    }
    visit(node)
    found
  }
  composite_calls <- find_calls(
    parse(file.path(module_dir, "run_figure7.R"), keep.source = FALSE),
    "figure7_main_composite_plots"
  )
  testthat::expect_length(composite_calls, 2L)
  call_arguments <- lapply(composite_calls, function(call) {
    vapply(as.list(call)[-1L], deparse1, character(1L))
  })
  testthat::expect_true(any(vapply(
    call_arguments,
    identical,
    logical(1L),
    c(
      "legacy_plots[LETTERS[1:5]]",
      "context_cache$plots",
      "legacy_plots$F",
      "config"
    )
  )))
  testthat::expect_true(any(vapply(
    call_arguments,
    identical,
    logical(1L),
    c("ae$plots", "context_cache$plots", "plot_f", "config")
  )))

  source_plots <- stats::setNames(
    lapply(
      paste0("7", LETTERS[1:5]),
      function(source_id) {
        ggplot2::ggplot() + ggplot2::labs(caption = source_id)
      }
    ),
    LETTERS[1:5]
  )
  context_plots <- stats::setNames(
    lapply(
      c("SI4A", "SI4B", "SI4C", "SI4E", "SI7B"),
      function(source_id) {
        ggplot2::ggplot() + ggplot2::labs(caption = source_id)
      }
    ),
    LETTERS[3:7]
  )
  state_pathway_plot <- ggplot2::ggplot() + ggplot2::labs(caption = "7F")
  assembled <- figure7_main_composite_plots(
    source_plots,
    context_plots,
    state_pathway_plot,
    config
  )
  testthat::expect_identical(names(assembled), LETTERS[1:11])
  testthat::expect_identical(
    vapply(assembled, function(plot) plot$labels$caption, character(1L)),
    expected_mapping
  )
})

testthat::test_that("promoted SI4 panels reproduce reviewed cache plots and normalized SI4E audit", {
  config <- figure7_test_inputs()$config
  cache_dir <- file.path(repo_root, "Data/in-vivo/SIfigures")
  context <- figure7_build_context_panels(cache_dir, repo_root, config)
  repeated <- figure7_build_context_panels(cache_dir, repo_root, config)

  testthat::expect_named(context$plots, c("C", "D", "E", "F", "G"))
  testthat::expect_true(all(vapply(
    context$plots[c("C", "D", "E", "F")],
    inherits,
    logical(1L),
    what = "ggplot"
  )))
  testthat::expect_identical(
    vapply(
      context$plots[c("C", "D", "E", "F")],
      function(plot) plot$labels$title,
      character(1L)
    ),
    c(
      C = "UMAP by cluster",
      D = "UMAP by initial ploidy",
      E = "UMAP by Tumor/CellLine context",
      F = "Equal-sample cluster composition by context"
    )
  )
  testthat::expect_identical(
    vapply(
      context$plots[c("C", "D", "E")],
      function(plot) rlang::as_label(plot$mapping$colour),
      character(1L)
    ),
    c(C = "cluster", D = "initial_ploidy", E = "context")
  )
  for (panel in c("C", "D", "E")) {
    testthat::expect_null(context$plots[[panel]]$labels$tag)
    testthat::expect_identical(
      context$plots[[panel]]$data$cell,
      repeated$plots[[panel]]$data$cell
    )
  }
  testthat::expect_null(context$plots$F$labels$tag)

  composition <- context$composition_result
  tests <- composition$tests
  enriched <- tests[tests$enriched, , drop = FALSE]
  testthat::expect_equal(nrow(tests), 18L)
  testthat::expect_true(all(tests$exact_permutations == 81L))
  testthat::expect_true(all(tests$strata == "initial_ploidy"))
  testthat::expect_true(all(tests$significance[tests$difference <= 0] == ""))
  testthat::expect_setequal(
    paste(enriched$group_value, enriched$cluster, sep = "/"),
    c(
      "Tumor/6", "Tumor/10", "Tumor/13",
      "CellLine/0", "CellLine/8", "CellLine/14"
    )
  )
  testthat::expect_equal(
    enriched$q_value,
    rep(0.037037037037037, 6L),
    tolerance = 1e-14
  )
  group_sums <- stats::aggregate(
    mean_proportion ~ group_value,
    data = composition$plot_data,
    FUN = sum
  )
  testthat::expect_equal(
    group_sums$mean_proportion,
    rep(1, nrow(group_sums)),
    tolerance = 1e-12
  )
  comparison_columns <- setdiff(names(composition$plot_data), "panel_tag")
  testthat::expect_equal(
    composition$plot_data[, comparison_columns, drop = FALSE],
    repeated$composition_result$plot_data[, comparison_columns, drop = FALSE],
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  testthat::expect_identical(
    context$cache_manifest_sha256,
    as.character(config$si_figures$reviewed_manifest_sha256)
  )
})

testthat::test_that("promoted SI7B deterministically clusters rows and columns with visible dendrograms", {
  config <- figure7_test_inputs()$config
  context <- figure7_build_context_panels(
    file.path(repo_root, "Data/in-vivo/SIfigures"),
    repo_root,
    config
  )
  first <- shared_context_build_heatmap(
    context$gsea_matrix,
    "Cluster Hallmark GSEA NES",
    TRUE,
    ""
  )
  second <- shared_context_build_heatmap(
    context$gsea_matrix,
    "Cluster Hallmark GSEA NES",
    TRUE,
    ""
  )

  testthat::expect_s3_class(first$tree_row, "hclust")
  testthat::expect_s3_class(first$tree_col, "hclust")
  testthat::expect_identical(first$tree_row$order, second$tree_row$order)
  testthat::expect_identical(first$tree_col$order, second$tree_col$order)
  testthat::expect_setequal(first$tree_row$labels, rownames(context$gsea_matrix))
  testthat::expect_setequal(first$tree_col$labels, colnames(context$gsea_matrix))
  testthat::expect_identical(
    colnames(context$gsea_matrix)[first$tree_col$order],
    c("0", "14", "13", "2", "5", "10", "4c", "6", "8")
  )

  promoted_gtable <- attr(context$plots$G, "grobs")$full
  testthat::expect_s3_class(promoted_gtable, "gtable")
  testthat::expect_setequal(
    intersect(promoted_gtable$layout$name, c("row_tree", "col_tree")),
    c("row_tree", "col_tree")
  )
  row_tree <- promoted_gtable$grobs[[
    which(promoted_gtable$layout$name == "row_tree")
  ]]
  col_tree <- promoted_gtable$grobs[[
    which(promoted_gtable$layout$name == "col_tree")
  ]]
  testthat::expect_s3_class(row_tree, "polyline")
  testthat::expect_s3_class(col_tree, "polyline")
})

testthat::test_that("full inventory requires composite while A-E inventory excludes it", {
  config <- figure7_test_inputs()$config
  composite <- as.character(config$panels$main_composite$filename)

  full_out <- tempfile("figure7_full_inventory_")
  dir.create(file.path(full_out, "figures"), recursive = TRUE)
  full_assets <- figure7_panel_asset_filenames(config, figure7_panel_ids(TRUE))
  testthat::expect_true(composite %in% full_assets)
  for (filename in setdiff(full_assets, composite)) {
    writeLines("fixture", file.path(full_out, "figures", filename))
  }
  testthat::expect_error(
    figure7_validate_figure_inventory(
      full_out,
      config,
      figure7_panel_ids(TRUE)
    ),
    "Figure inventory mismatch"
  )
  writeLines("fixture", file.path(full_out, "figures", composite))
  testthat::expect_silent(figure7_validate_figure_inventory(
    full_out,
    config,
    figure7_panel_ids(TRUE)
  ))

  generated_out <- tempfile("figure7_generated_inventory_")
  dir.create(file.path(generated_out, "figures"), recursive = TRUE)
  generated_composite <- figure7_main_composite_filename(
    config,
    "generated-human-only"
  )
  generated_assets <- figure7_panel_asset_filenames(
    config,
    figure7_panel_ids(TRUE),
    generated_composite
  )
  for (filename in generated_assets) {
    writeLines("fixture", file.path(generated_out, "figures", filename))
  }
  testthat::expect_silent(figure7_validate_figure_inventory(
    generated_out,
    config,
    figure7_panel_ids(TRUE),
    generated_composite
  ))
  testthat::expect_false(composite %in% generated_assets)

  ae_out <- tempfile("figure7_ae_inventory_")
  dir.create(file.path(ae_out, "figures"), recursive = TRUE)
  ae_ids <- figure7_panel_ids(FALSE)
  ae_assets <- figure7_panel_asset_filenames(config, ae_ids)
  testthat::expect_false(composite %in% ae_assets)
  testthat::expect_false(
    "7A-7K_composite" %in% figure7_panel_contract(config, ae_ids)$panel_id
  )
  for (filename in ae_assets) {
    writeLines("fixture", file.path(ae_out, "figures", filename))
  }
  testthat::expect_silent(figure7_validate_figure_inventory(
    ae_out,
    config,
    ae_ids
  ))
  writeLines("fixture", file.path(ae_out, "figures", composite))
  testthat::expect_error(
    figure7_validate_figure_inventory(ae_out, config, ae_ids),
    "Figure inventory mismatch"
  )
})

testthat::test_that("tampered reviewed SI cache fails before Figure 7 creates output", {
  source_cache <- file.path(repo_root, "Data/in-vivo/SIfigures")
  tampered_cache <- tempfile("figure7_si_cache_tampered_")
  dir.create(tampered_cache)
  source_files <- list.files(source_cache, full.names = TRUE)
  testthat::expect_true(all(file.copy(source_files, tampered_cache)))
  tampered_file <- file.path(tampered_cache, "si_figures_cluster_key.tsv")
  writeLines(
    c(readLines(tampered_file, warn = FALSE), "tampered"),
    tampered_file,
    useBytes = TRUE
  )

  config <- figure7_test_inputs()$config
  testthat::expect_error(
    figure7_validate_reviewed_si_cache(tampered_cache, repo_root, config),
    "Reviewed SI cache checksum mismatch"
  )

  output_dir <- tempfile("figure7_tampered_no_output_")
  command_output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      file.path(module_dir, "run_figure7.R"),
      "--mode=standard",
      "--panel-set=a-f",
      paste0("--config=", file.path(module_dir, "figure7_config.yaml")),
      paste0("--si-table-cache-dir=", tampered_cache),
      paste0("--output-dir=", output_dir)
    ),
    stdout = TRUE,
    stderr = TRUE
  ))
  testthat::expect_true(
    !is.null(attr(command_output, "status")) &&
      attr(command_output, "status") != 0L
  )
  testthat::expect_match(
    paste(command_output, collapse = "\n"),
    "Reviewed SI cache checksum mismatch"
  )
  testthat::expect_false(dir.exists(output_dir))
})

testthat::test_that("generated SI context cache is reproducible but explicitly noncanonical", {
  config <- figure7_test_inputs()$config
  source_cache <- file.path(repo_root, "Data/in-vivo/SIfigures")
  fixture_root <- tempfile("figure7_generated_context_")
  cache_dir <- file.path(fixture_root, "tables")
  metadata_dir <- file.path(fixture_root, "metadata")
  dir.create(cache_dir, recursive = TRUE)
  dir.create(metadata_dir, recursive = TRUE)
  testthat::expect_true(all(file.copy(
    list.files(source_cache, full.names = TRUE),
    cache_dir
  )))

  manifest_path <- file.path(cache_dir, "manifest.tsv")
  manifest <- utils::read.delim(
    manifest_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  manifest$source_revision <- "raw-generated-human-only@test-fixture"
  si7_rows <- startsWith(
    manifest$filename,
    "si_figure7_cluster_Hallmark"
  )
  manifest$notes[si7_rows] <- paste(
    "generated human-only GRCh policy;",
    "not approved for canonical publication"
  )
  utils::write.table(
    manifest,
    manifest_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  cache_paths <- file.path(cache_dir, manifest$filename)
  input_manifest <- data.frame(
    role = c(
      "figure7_config",
      "normalized_composition_helper",
      "shared_context_panels_helper",
      "si_figures_generated_cache_manifest",
      rep("si_figures_generated_table", nrow(manifest)),
      "raw_build_input_manifest",
      "raw_build_si_raw_scientific_code_contract",
      "raw_build_si_raw_config_contract",
      "raw_build_si_raw_runtime_contract",
      "raw_build_versioned_endpoint_ploidy",
      "raw_build_source_seurat_rds"
    ),
    repo_relative_path = c(
      "Code/in-vivo/figure7/figure7_config.yaml",
      "Code/in-vivo/SI_figures/normalized_composition.R",
      "Code/in-vivo/SI_figures/shared_context_panels.R",
      "external:manifest.tsv",
      paste0("external:", manifest$filename),
      "external:raw_input_manifest.tsv",
      "contract:si_raw_scientific_code_contract",
      "contract:si_raw_config_contract",
      "contract:si_raw_runtime_contract",
      "Data/in-vivo/all_ploidy.tsv",
      "external:integrated.rds"
    ),
    sha256 = c(
      figure7_sha256(file.path(
        repo_root,
        "Code/in-vivo/figure7/figure7_config.yaml"
      )),
      figure7_sha256(file.path(
        repo_root,
        "Code/in-vivo/SI_figures/normalized_composition.R"
      )),
      figure7_sha256(file.path(
        repo_root,
        "Code/in-vivo/SI_figures/shared_context_panels.R"
      )),
      figure7_sha256(manifest_path),
      vapply(cache_paths, figure7_sha256, character(1L)),
      strrep("a", 64L),
      strrep("b", 64L),
      strrep("c", 64L),
      strrep("d", 64L),
      strrep("e", 64L),
      strrep("f", 64L)
    ),
    bytes = c(
      file.info(file.path(
        repo_root,
        "Code/in-vivo/figure7/figure7_config.yaml"
      ))$size,
      file.info(file.path(
        repo_root,
        "Code/in-vivo/SI_figures/normalized_composition.R"
      ))$size,
      file.info(file.path(
        repo_root,
        "Code/in-vivo/SI_figures/shared_context_panels.R"
      ))$size,
      file.info(manifest_path)$size,
      file.info(cache_paths)$size,
      rep(0, 6L)
    ),
    stringsAsFactors = FALSE
  )
  upstream_path <- file.path(metadata_dir, "analysis_input_manifest.tsv")
  utils::write.table(
    input_manifest,
    upstream_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  source_revision <- paste(unique(manifest$source_revision), collapse = ";")
  run_config_path <- file.path(metadata_dir, "run_config.tsv")
  figure7_write_tsv(
    data.frame(
      key = c(
        "table_mode",
        "si7_canonical_publication_allowed",
        "cache_file_count",
        "generated_table_source_revision"
      ),
      value = c(
        "run_scoped_generated_human_only_plot_tables",
        "false",
        as.character(nrow(manifest)),
        source_revision
      ),
      stringsAsFactors = FALSE
    ),
    run_config_path
  )
  provenance_path <- file.path(metadata_dir, "si_figures_provenance.tsv")
  figure7_write_tsv(
    data.frame(
      key = c(
        "entrypoint_sha256",
        "normalized_composition_helper_sha256",
        "shared_context_panels_helper_sha256",
        "table_cache_manifest_sha256",
        "si7_canonical_publication_allowed"
      ),
      value = c(
        figure7_sha256(file.path(
          repo_root,
          "Code/in-vivo/SI_figures/generate_supplementary_figures.R"
        )),
        figure7_sha256(file.path(
          repo_root,
          "Code/in-vivo/SI_figures/normalized_composition.R"
        )),
        figure7_sha256(file.path(
          repo_root,
          "Code/in-vivo/SI_figures/shared_context_panels.R"
        )),
        figure7_sha256(manifest_path),
        "false"
      ),
      stringsAsFactors = FALSE
    ),
    provenance_path
  )

  generated <- figure7_build_context_panels(
    cache_dir,
    repo_root,
    config,
    policy = "generated-human-only",
    upstream_input_manifest = upstream_path,
    source_run_config = run_config_path,
    source_provenance = provenance_path
  )
  testthat::expect_named(generated$plots, LETTERS[3:7])
  testthat::expect_identical(
    generated$cache_kind,
    "generated_human_only_run_scoped"
  )
  testthat::expect_false(generated$cache_canonical_publication_allowed)
  testthat::expect_identical(
    generated$composite_filename,
    "Figure7_generated_GRCh_candidate.png"
  )
  testthat::expect_error(
    figure7_validate_reviewed_si_cache(cache_dir, repo_root, config),
    "reviewed SI Figures cache manifest"
  )

  incomplete <- input_manifest[
    input_manifest$role != "raw_build_si_raw_config_contract",
    ,
    drop = FALSE
  ]
  utils::write.table(
    incomplete,
    upstream_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  testthat::expect_error(
    figure7_validate_si_cache(
      cache_dir,
      repo_root,
      config,
      "generated-human-only",
      upstream_path,
      run_config_path,
      provenance_path
    ),
    "lineage is incomplete"
  )
})
