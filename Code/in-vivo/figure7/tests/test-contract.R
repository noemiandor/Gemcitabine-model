testthat::test_that("panel-F compact reference contract validates and selector is audited", {
  fixture <- figure7_test_reviewed_state_reference()
  reference <- figure7_validate_reviewed_state_reference(
    fixture$path,
    fixture$config
  )
  testthat::expect_equal(nrow(reference$pathways), 21L)
  testthat::expect_s3_class(figure7_panel_f_plot(reference$activity, fixture$config), "ggplot")
  collision_keys <- unlist(lapply(
    c("H", "C2:CP:REACTOME"),
    function(collection) {
      reference$pathways$pathway_id[
        match(collection, reference$pathways$collection_id)
      ]
    }
  ))
  collision_activity <- reference$activity[
    reference$activity$pathway_id %in% collision_keys,
    ,
    drop = FALSE
  ]
  collision_activity$pathway_label <- "Shared display label"
  collision_plot <- figure7_panel_f_plot(
    collision_activity,
    fixture$config
  )
  testthat::expect_silent(ggplot2::ggplot_build(collision_plot))
  testthat::expect_equal(
    nlevels(collision_plot$data$pathway_plot_key),
    2L
  )
  bad <- reference$selected
  bad$selected_rank_within_direction[[1L]] <- 99L
  figure7_write_tsv(bad, file.path(fixture$path, "panel_7F_selected_pathway_gsea.tsv"))
  testthat::expect_error(
    figure7_validate_reviewed_state_reference(
      fixture$path,
      fixture$config,
      verify_checksums = FALSE
    ),
    "metadata/order|configured selector"
  )
})

testthat::test_that("tracked v2 is byte-pinned but superseded after the endpoint-score confound", {
  input <- figure7_test_inputs()
  reference_path <- file.path(
    repo_root,
    input$config$state_pathways$reviewed_reference_root,
    "state_pathway_grch_human_only_etp2_24_day17_v2"
  )
  reference <- figure7_validate_superseded_v2_state_reference(
    reference_path,
    input$config
  )
  provenance <- stats::setNames(
    as.character(reference$provenance$value),
    reference$provenance$key
  )
  collection_counts <- table(reference$pathways$collection_id)

  testthat::expect_identical(
    reference$reference_id,
    "state_pathway_grch_human_only_etp2_24_day17_v2"
  )
  testthat::expect_identical(
    reference$reference_kind,
    "superseded_human_only_endpoint_cn_score_frozen"
  )
  testthat::expect_false(reference$canonical_publication_allowed)
  testthat::expect_equal(nrow(reference$selected), 21L)
  testthat::expect_equal(nrow(reference$activity), 21L * 501L)
  testthat::expect_identical(
    as.integer(collection_counts[c("H", "C2:CP:REACTOME", "C5:GO:BP")]),
    c(5L, 8L, 8L)
  )
  testthat::expect_lte(max(figure7_numeric(reference$selected$padj)), 0.05)
  testthat::expect_identical(
    provenance[["reviewed_source_provenance_sha256"]],
    "c325359b1d0fe67c3ad1f13523e993cb80e5168d164d06e1b2fa1b1b584be888"
  )
  testthat::expect_identical(
    provenance[["canonical_publication_allowed"]],
    "true"
  )
  testthat::expect_false(any(grepl(
    "^(/|[A-Za-z]:[\\\\/])|/(private|share|Users)/",
    provenance
  )))
})

testthat::test_that("tracked pointwise v4 candidate uses injected initial ploidy and GRCh-only features", {
  input <- figure7_test_inputs()
  reference_path <- file.path(
    repo_root,
    input$config$state_pathways$reviewed_reference_root,
    input$config$state_pathways$generated_reference_id
  )
  reference <- figure7_validate_generated_state_reference(
    reference_path,
    input$config,
    expected_inputs = list(
      cellcycle = file.path(
        repo_root,
        "Data/in-vivo/figure7/processed",
        "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
      ),
      noncellcycle = file.path(
        repo_root,
        "Data/in-vivo/figure7/processed",
        "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
      ),
      seurat_rds = NA_character_
    ),
    config_path = NULL
  )
  provenance <- stats::setNames(
    as.character(reference$provenance$value),
    reference$provenance$key
  )
  design <- figure7_read_tsv(
    file.path(reference_path, "state_pathway_design_qc.tsv")
  )

  testthat::expect_identical(
    reference$reference_id,
    "state_pathway_grch_human_only_initial_ploidy_day17_pointwise_v4_candidate"
  )
  testthat::expect_identical(
    reference$reference_kind,
    "generated_human_only_initial_ploidy_computed_pointwise_interval_candidate"
  )
  testthat::expect_false(reference$canonical_publication_allowed)
  testthat::expect_equal(nrow(reference$selected), 21L)
  testthat::expect_lte(max(figure7_numeric(reference$selected$padj)), 0.05)
  testthat::expect_identical(
    provenance[["nuisance_policy"]],
    "injected_initial_ploidy_only_no_endpoint_cn_score"
  )
  testthat::expect_identical(
    provenance[["nuisance_terms"]],
    "dose_mg_factor,initial_ploidy_factor"
  )
  testthat::expect_identical(
    provenance[["endpoint_cn_score_covariate_prohibited"]],
    "true"
  )
  testthat::expect_identical(
    provenance[["feature_species_policy_id"]],
    "grch_human_tumor_only_v2"
  )
  testthat::expect_identical(design$covariate_mode, "initial_ploidy")
  testthat::expect_identical(design$initial_ploidy_levels, "2N;4N")
  testthat::expect_match(
    design$retained_design_columns,
    "(^|;)initial_ploidy_factor4N($|;)"
  )
  testthat::expect_false(any(grepl(
    "etp|endpoint|cn_score",
    design$retained_design_columns,
    ignore.case = TRUE
  )))
})

testthat::test_that("reviewed pointwise v4 is the exact promoted initial-ploidy-adjusted result", {
  input <- figure7_test_inputs()
  reviewed_path <- file.path(
    repo_root,
    input$config$state_pathways$reviewed_reference_root,
    input$config$state_pathways$reviewed_reference_id
  )
  candidate_path <- file.path(
    repo_root,
    input$config$state_pathways$reviewed_reference_root,
    input$config$state_pathways$generated_reference_id
  )
  reviewed <- figure7_validate_reviewed_state_reference(
    reviewed_path,
    input$config
  )
  candidate <- figure7_validate_generated_state_reference(
    candidate_path,
    input$config,
    expected_inputs = NULL,
    config_path = NULL
  )
  provenance <- stats::setNames(
    as.character(reviewed$provenance$value),
    reviewed$provenance$key
  )
  keys <- c(
    "collection_id", "pathway_id", "NES", "padj",
    "selected_direction", "selected_rank_within_direction"
  )

  testthat::expect_identical(
    reviewed$reference_id,
    "state_pathway_grch_human_only_initial_ploidy_day17_pointwise_v4"
  )
  testthat::expect_identical(
    reviewed$reference_kind,
    "reviewed_human_only_initial_ploidy_computed_pointwise_interval"
  )
  testthat::expect_true(reviewed$canonical_publication_allowed)
  testthat::expect_identical(reviewed$selected[, keys], candidate$selected[, keys])
  testthat::expect_lte(max(figure7_numeric(reviewed$selected$padj)), 0.05)
  testthat::expect_identical(
    provenance[["reviewed_source_run_id"]],
    "20260802_figure7_pointwise_v4_candidate_review_figure7"
  )
  testthat::expect_identical(
    provenance[["endpoint_cn_score_covariate_prohibited"]],
    "true"
  )
})

testthat::test_that("reviewed pointwise v4 rejects provenance or selected-table tampering", {
  input <- figure7_test_inputs()
  source <- file.path(
    repo_root,
    input$config$state_pathways$reviewed_reference_root,
    input$config$state_pathways$reviewed_reference_id
  )
  parent <- tempfile("figure7_reviewed_tamper_")
  tampered <- file.path(
    parent,
    input$config$state_pathways$reviewed_reference_id
  )
  dir.create(tampered, recursive = TRUE)
  testthat::expect_true(all(file.copy(
    file.path(source, figure7_state_required_files()),
    file.path(tampered, figure7_state_required_files())
  )))

  provenance_path <- file.path(tampered, "state_pathway_provenance.tsv")
  provenance <- figure7_read_tsv(provenance_path, c("key", "value"))
  provenance$value[provenance$key == "reviewed_on"] <- "2026-08-01"
  figure7_write_tsv(provenance, provenance_path)
  testthat::expect_error(
    figure7_validate_reviewed_state_reference(tampered, input$config),
    "SHA-256 mismatch"
  )
  testthat::expect_error(
    figure7_validate_reviewed_state_reference(
      tampered,
      input$config,
      verify_checksums = FALSE
    ),
    "initial-ploidy-adjusted pointwise-interval v4"
  )
})

testthat::test_that("canonical provenance is location-independent and always checksummed", {
  fixture <- figure7_test_reviewed_state_reference()
  provenance_path <- file.path(fixture$path, "state_pathway_provenance.tsv")
  provenance <- figure7_read_tsv(provenance_path, c("key", "value"))
  provenance <- rbind(
    provenance,
    data.frame(
      key = "runtime_source_root",
      value = "/runtime/generated/results",
      stringsAsFactors = FALSE
    )
  )
  figure7_write_tsv(provenance, provenance_path)

  testthat::expect_error(
    figure7_validate_reviewed_state_reference(
      fixture$path,
      fixture$config
    ),
    "SHA-256 mismatch"
  )
  testthat::expect_error(
    figure7_validate_reviewed_state_reference(
      fixture$path,
      fixture$config,
      verify_checksums = FALSE
    ),
    "runtime filesystem paths"
  )
})

testthat::test_that("named source checksum verification rejects any changed table", {
  paths <- c(first = tempfile(), second = tempfile())
  writeLines("first", paths[["first"]], useBytes = TRUE)
  writeLines("second", paths[["second"]], useBytes = TRUE)
  expected <- vapply(paths, figure7_sha256, character(1L))
  testthat::expect_silent(figure7_verify_named_checksums(paths, expected))
  writeLines("changed", paths[["second"]], useBytes = TRUE)
  testthat::expect_error(figure7_verify_named_checksums(paths, expected), "SHA-256 mismatch")
})

testthat::test_that("TSV helpers round-trip multiline annotations without malformed rows", {
  path <- tempfile(fileext = ".tsv")
  original <- data.frame(id = 1:2, annotation = c("line one\nline two", "plain"), stringsAsFactors = FALSE)
  figure7_write_tsv(original, path)
  testthat::expect_length(readLines(path, warn = FALSE), 3L)
  observed <- figure7_read_tsv(path, c("id", "annotation"))
  testthat::expect_identical(observed$annotation, original$annotation)
})

testthat::test_that("missing F, wrong checksums, and nonempty outputs fail clearly", {
  input <- figure7_test_inputs()
  missing <- file.path(
    tempdir(),
    input$config$state_pathways$reviewed_reference_id
  )
  testthat::expect_error(
    figure7_validate_reviewed_state_reference(missing, input$config),
    "Missing reviewed human-only panel-7F"
  )
  wrong <- tempfile(); writeLines("wrong", wrong)
  testthat::expect_error(figure7_verify_checksum(wrong, paste(rep("0", 64), collapse = "")), "SHA-256 mismatch")
  nonempty <- tempfile(); dir.create(nonempty); writeLines("x", file.path(nonempty, "existing.txt"))
  testthat::expect_error(figure7_prepare_output(nonempty), "pre-existing analysis files")
})

testthat::test_that("entire-run image inventory is exact", {
  config <- figure7_test_inputs()$config; out <- tempfile(); dir.create(file.path(out, "figures"), recursive = TRUE)
  for (file in figure7_panel_asset_filenames(config)) writeLines("figure fixture", file.path(out, "figures", file))
  testthat::expect_silent(figure7_validate_figure_inventory(out, config))
  dir.create(file.path(out, "tables")); file.create(file.path(out, "tables", "unexpected.svg"))
  testthat::expect_error(figure7_validate_figure_inventory(out, config), "Figure inventory mismatch")
})

testthat::test_that("publication compositor pins dimensions, mapping, and vector asset names", {
  config <- figure7_read_config(file.path(module_dir, "figure7_config.yaml"))
  spec <- figure7_publication_spec()
  expected_mapping <- c(
    A = "7A", B = "7C", C = "SI4A", D = "SI4B", E = "SI4C",
    F = "SI4E", G = "SI7B", H = "7B", I = "7F", J = "7J",
    K = "7D", L = "7E"
  )
  configured_mapping <- unlist(
    config$panels$main_composite$panel_order,
    use.names = TRUE
  )

  testthat::expect_identical(spec$width_in, 7.1)
  testthat::expect_identical(spec$height_in, 10.645)
  testthat::expect_identical(
    spec$row_heights,
    c(1.20, 1.35, 2.45, 1.45, 2.745, 1.45)
  )
  testthat::expect_equal(spec$row_heights[[2L]] / 0.95, 1.4210526, tolerance = 1e-7)
  testthat::expect_true(all(spec$row_heights > 0))
  testthat::expect_equal(sum(spec$row_heights), spec$height_in, tolerance = 0)
  testthat::expect_identical(
    strsplit(spec$design, "\n", fixed = TRUE)[[1L]],
    c(
      paste0(strrep("A", 18L), strrep("B", 10L)),
      paste0(strrep("C", 9L), strrep("D", 10L), strrep("E", 9L)),
      paste0(strrep("F", 10L), strrep("G", 18L)),
      strrep("H", 28L),
      paste0(strrep("I", 11L), strrep("J", 17L)),
      paste0(strrep("K", 14L), strrep("L", 14L))
    )
  )
  testthat::expect_identical(spec$content_left_npc, 0.020)
  testthat::expect_identical(spec$content_right_npc, 0.005)
  testthat::expect_equal(
    spec$width_in * (1 - spec$content_left_npc - spec$content_right_npc),
    6.9225,
    tolerance = 1e-12
  )
  testthat::expect_identical(configured_mapping, expected_mapping)

  make_plot <- function(tag) {
    ggplot2::ggplot(
      data.frame(x = 1, y = 1),
      ggplot2::aes(x = x, y = y)
    ) +
      ggplot2::geom_point() +
      ggplot2::labs(tag = tag)
  }
  source_plots <- stats::setNames(
    lapply(paste0("source-", LETTERS[1:5]), make_plot),
    LETTERS[1:5]
  )
  attr(source_plots$B, "figure7_panel_b_components") <- list(
    ecdf = source_plots$B,
    localization = make_plot("localization")
  )
  context_plots <- stats::setNames(
    lapply(paste0("context-", LETTERS[3:7]), make_plot),
    LETTERS[3:7]
  )
  mapped <- figure7_main_composite_plots(
    source_plots,
    context_plots,
    make_plot("state-pathway"),
    make_plot("copy-number"),
    config
  )
  testthat::expect_identical(names(mapped), LETTERS[1:12])
  testthat::expect_identical(
    vapply(mapped, function(plot) plot$labels$tag, character(1L)),
    c(
      A = "source-A", B = "source-C", C = "context-C",
      D = "context-D", E = "context-E", F = "context-F",
      G = "context-G", H = "source-B", I = "state-pathway",
      J = "copy-number", K = "source-D", L = "source-E"
    )
  )

  composite_png <- figure7_main_composite_filename(config)
  composite_pdf <- figure7_main_composite_pdf_filename(composite_png)
  assets <- figure7_panel_asset_filenames(config)
  testthat::expect_identical(composite_png, "Figure7_reviewed_GRCh.png")
  testthat::expect_identical(composite_pdf, "Figure7_reviewed_GRCh.pdf")
  testthat::expect_true(all(c(composite_png, composite_pdf) %in% assets))
  testthat::expect_equal(sum(assets == composite_png), 1L)
  testthat::expect_equal(sum(assets == composite_pdf), 1L)
  testthat::expect_error(
    figure7_main_composite_pdf_filename("Figure7_reviewed_GRCh.tiff"),
    "must end in .png"
  )
})

testthat::test_that("publication styling removes internal prose and uses reader-facing labels", {
  config <- figure7_read_config(file.path(module_dir, "figure7_config.yaml"))
  plots <- stats::setNames(lapply(LETTERS[1:12], function(panel) {
    ggplot2::ggplot(
      data.frame(x = 1:2, y = 1:2),
      ggplot2::aes(x = x, y = y)
    ) +
      ggplot2::geom_point() +
      ggplot2::labs(
        title = paste("code-facing title", panel),
        subtitle = "implementation formula",
        caption = "methodological prose",
        tag = panel
      )
  }), LETTERS[1:12])

  styled <- figure7_publication_clean_plots(plots, config)
  for (panel in names(styled)) {
    testthat::expect_null(styled[[panel]]$labels$title, info = panel)
    testthat::expect_null(styled[[panel]]$labels$subtitle, info = panel)
    testthat::expect_null(styled[[panel]]$labels$caption, info = panel)
    testthat::expect_null(styled[[panel]]$labels$tag, info = panel)
  }
  testthat::expect_identical(styled$A$labels$x, "Days since first treatment")
  testthat::expect_identical(styled$B$labels$x, "Injected origin")
  testthat::expect_identical(styled$F$labels$y, "Mean proportion per sample")
  testthat::expect_identical(styled$H$labels$y, "Mean ECDF")
  testthat::expect_identical(styled$I$labels$fill, "Mean gene z score")
  testthat::expect_identical(styled$D$theme$legend.position, "inside")
  testthat::expect_identical(styled$E$theme$legend.position, "inside")
  testthat::expect_identical(
    styled$K$labels$x,
    "Pseudotime-distribution shift\n(dose-centered ECDF RMSE)"
  )
  testthat::expect_identical(
    styled$L$labels$x,
    "Mean endpoint tumor-cell ploidy"
  )

  composite <- suppressWarnings(figure7_main_composite_object(plots, config))
  testthat::expect_s3_class(composite, "gTree")
  testthat::expect_true(all(
    paste0("figure7_tag_", setdiff(LETTERS[1:12], "J")) %in%
      names(composite$children)
  ))
  testthat::expect_false("figure7_tag_J" %in% names(composite$children))
  locally_tagged_j <- figure7_publication_local_panel_tag(
    plots$J, "J", figure7_publication_spec()
  )
  testthat::expect_identical(locally_tagged_j$labels$tag, "J")
})

testthat::test_that("main H keeps ECDFs while localization remains a separate source component", {
  input <- figure7_test_inputs()
  panel <- figure7_panel_b(input$data, input$samples, input$config)
  panel_h <- figure7_panel_b_plot(
    panel$data,
    panel$tests,
    panel$localization_grid,
    panel$localization_intervals,
    panel$localization_test
  )
  components <- attr(panel_h, "figure7_panel_b_components", exact = TRUE)

  testthat::expect_s3_class(panel_h, "patchwork")
  testthat::expect_identical(names(components), c("ecdf", "localization"))
  testthat::expect_identical(
    components$localization$labels$y,
    "Gemcitabine - vehicle\ndensity difference"
  )
  testthat::expect_equal(
    panel$localization_intervals$start,
    c(0.296, 0.414),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    panel$localization_intervals$end,
    c(0.486, 0.426),
    tolerance = 1e-12
  )

  make_plot <- function() {
    ggplot2::ggplot(
      data.frame(x = 1:2, y = 1:2),
      ggplot2::aes(x, y)
    ) + ggplot2::geom_point()
  }
  source_plots <- stats::setNames(
    lapply(LETTERS[1:5], function(...) make_plot()),
    LETTERS[1:5]
  )
  source_plots$B <- panel_h
  context_plots <- stats::setNames(
    lapply(LETTERS[3:7], function(...) make_plot()),
    LETTERS[3:7]
  )
  mapped <- figure7_main_composite_plots(
    source_plots,
    context_plots,
    make_plot(),
    make_plot(),
    input$config
  )
  testthat::expect_s3_class(mapped$H, "ggplot")
  testthat::expect_null(attr(
    mapped$H, "figure7_panel_b_components", exact = TRUE
  ))
  styled <- figure7_publication_clean_plots(mapped, input$config)
  testthat::expect_identical(
    styled$H$labels$y,
    "Mean ECDF"
  )
  testthat::expect_identical(
    styled$H$labels$x,
    "Cell-cycle pseudotime"
  )
})

testthat::test_that("full source panels plus the A-L composite satisfy the exact inventory", {
  input <- figure7_test_inputs()
  fixture <- figure7_test_reviewed_state_reference()
  reference <- figure7_validate_reviewed_state_reference(
    fixture$path,
    fixture$config
  )
  out <- tempfile("figure7_six_"); figure7_prepare_output(out)
  ae <- figure7_build_ae(
    input$cellcycle,
    input$data,
    input$samples,
    out,
    fixture$config
  )
  plot_f <- figure7_build_f(reference, out, fixture$config)
  context <- figure7_build_context_panels(
    file.path(repo_root, "Data/in-vivo/SIfigures"),
    repo_root,
    fixture$config
  )
  copy_number <- figure7_build_copy_number_panel(repo_root)
  testthat::expect_identical(copy_number$column_width_multiplier, 2)
  testthat::expect_identical(
    copy_number$heatmap$column_width_multiplier,
    2
  )
  testthat::expect_identical(
    copy_number$row_ordering_policy,
    "hierarchical_clustering_separately_within_each_mouse"
  )
  testthat::expect_identical(copy_number$row_distance_method, "euclidean")
  testthat::expect_identical(copy_number$row_linkage_method, "ward.D2")
  testthat::expect_identical(
    copy_number$row_tie_break_method,
    "canonical_heatmap_row_id_input_order"
  )
  testthat::expect_equal(nrow(copy_number$row_order_audit), 9832L)
  testthat::expect_identical(
    unique(copy_number$row_order_audit$sample_id),
    names(copy_number$heatmap$annotation_colors$Mouse)
  )
  testthat::expect_equal(
    sum(vapply(
      context$plots$C$layers,
      function(layer) inherits(layer$geom, "GeomLabelRepel"),
      logical(1L)
    )),
    1L
  )
  figure7_save_main_composite(
    figure7_main_composite_plots(
      ae$plots,
      context$plots,
      plot_f,
      copy_number$plot,
      fixture$config
    ),
    out,
    fixture$config,
    width = 12,
    height = 14,
    png_dpi = 72
  )
  testthat::expect_silent(figure7_validate_figure_inventory(out, fixture$config))
  pdfs <- list.files(file.path(out, "figures"), pattern = "[.]pdf$", full.names = TRUE)
  pngs <- list.files(file.path(out, "figures"), pattern = "[.]png$", full.names = TRUE)
  testthat::expect_length(pdfs, 7L)
  testthat::expect_length(pngs, 7L)
  testthat::expect_true(file.exists(file.path(
    out,
    "figures",
    "Figure7_reviewed_GRCh.png"
  )))
  testthat::expect_true(file.exists(file.path(
    out,
    "figures",
    "Figure7_reviewed_GRCh.pdf"
  )))
  testthat::expect_true(all(file.info(pngs)$size > 0))
  if (nzchar(Sys.which("pdfinfo"))) {
    statuses <- vapply(pdfs, function(file) system2("pdfinfo", file, stdout = FALSE, stderr = FALSE), integer(1L))
    testthat::expect_true(all(statuses == 0L))
  }
})

testthat::test_that("entrypoint fails before output when canonical F is unavailable", {
  input <- figure7_test_inputs(); out <- tempfile("figure7_no_partial_")
  args <- c(file.path(module_dir, "run_figure7.R"), "--mode=standard",
    paste0("--config=", file.path(module_dir, "figure7_config.yaml")),
    paste0("--cellcycle-input=", file.path(repo_root, "Data/in-vivo/figure7/processed/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv")),
    paste0("--non-cellcycle-input=", file.path(repo_root, "Data/in-vivo/figure7/processed/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv")),
    paste0(
      "--saved-state-pathway-dir=",
      file.path(
        tempdir(),
        input$config$state_pathways$reviewed_reference_id
      )
    ),
    paste0("--output-dir=", out))
  status <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), args, stdout = TRUE, stderr = TRUE))
  testthat::expect_true(!is.null(attr(status, "status")) && attr(status, "status") != 0L)
  testthat::expect_false(dir.exists(out))
})

testthat::test_that("render-only rejects a source without frozen run metadata before output", {
  source <- tempfile(); dir.create(source); out <- tempfile("figure7_render_no_partial_")
  args <- c(file.path(module_dir, "run_figure7.R"), "--mode=render-only",
    paste0("--config=", file.path(module_dir, "figure7_config.yaml")),
    paste0("--source-run-dir=", source), paste0("--output-dir=", out))
  status <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), args, stdout = TRUE, stderr = TRUE))
  testthat::expect_true(!is.null(attr(status, "status")) && attr(status, "status") != 0L)
  testthat::expect_false(dir.exists(out))
})

testthat::test_that("SI8 compact provenance is complete and relocatable", {
  validation_output <- file.path(
    repo_root,
    "figures",
    paste0(".si8_contract_", basename(tempfile()))
  )
  on.exit(unlink(validation_output, recursive = TRUE, force = TRUE), add = TRUE)
  command <- c(
    file.path(module_dir, "assemble_tgi_sensitivity.R"),
    paste0("--output-dir=", validation_output)
  )
  status <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    command,
    stdout = TRUE,
    stderr = TRUE
  ))
  testthat::expect_null(attr(status, "status"), info = paste(status, collapse = "\n"))
  provenance_path <- file.path(
    validation_output,
    "Figure7_Supplement_provenance.tsv"
  )
  provenance <- figure7_read_tsv(provenance_path, c("key", "value"))
  expected_roles <- c("cellcycle", "noncellcycle", unlist(lapply(c("day24", "day31"), function(day) {
    paste0(
      day,
      c(
        "_run_config", "_endpoint_ploidy", "_panel_a",
        "_panel_c_data", "_panel_c_test", "_panel_e_data",
        "_panel_e_test"
      )
    )
  }), use.names = FALSE))

  testthat::expect_identical(anyDuplicated(provenance$key), 0L)
  testthat::expect_setequal(
    sub("^source_file:", "", provenance$key[startsWith(
      provenance$key,
      "source_file:"
    )]),
    expected_roles
  )
  testthat::expect_setequal(
    sub("^source_sha256:", "", provenance$key[startsWith(
      provenance$key,
      "source_sha256:"
    )]),
    expected_roles
  )
  source_locators <- provenance$value[startsWith(
    provenance$key,
    "source_file:"
  )]
  keyed <- stats::setNames(as.character(provenance$value), provenance$key)
  bundle_locator <- paste0(
    "Data/in-vivo/figure7/saved_tgi_sensitivity/",
    "tgi_day24_day31_curated_cbs_v3_raw_pearson"
  )
  expected_bundle_files <- unlist(lapply(c("day24", "day31"), function(day) {
    c(
      file.path(day, "metadata", "run_config.tsv"),
      file.path(
        day,
        "tables",
        c(
          "panel_7A_plot_data.tsv",
          "panel_7C_plot_data.tsv",
          "panel_7C_test.tsv",
          "panel_7E_plot_data.tsv",
          "panel_7E_test.tsv"
        )
      )
    )
  }), use.names = FALSE)
  bundle_path <- file.path(repo_root, bundle_locator)

  testthat::expect_identical(
    keyed[["source_bundle_id"]],
    "tgi_day24_day31_curated_cbs_v3_raw_pearson"
  )
  testthat::expect_identical(keyed[["source_bundle"]], bundle_locator)
  testthat::expect_identical(keyed[["source_bundle_file_count"]], "12")
  testthat::expect_identical(
    keyed[["endpoint_ploidy_inventory_n_cells"]],
    "14125"
  )
  testthat::expect_identical(
    keyed[["endpoint_ploidy_score_universe_n_cells"]],
    "9832"
  )
  testthat::expect_identical(
    keyed[["treated_endpoint_ploidy_n_cells"]],
    "5335"
  )
  testthat::expect_equal(
    as.numeric(keyed[["day24_pearson_r"]]),
    -0.387348978978662,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    as.numeric(keyed[["day24_exact_unrestricted_permutation_p"]]),
    0.346924603174603,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    as.numeric(keyed[["day31_pearson_r"]]),
    -0.296948455789288,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    as.numeric(keyed[["day31_exact_unrestricted_permutation_p"]]),
    0.485714285714286,
    tolerance = 1e-12
  )
  testthat::expect_false(any(c(
    "day24_adjusted_slope_per_origin_sd", "day24_partial_r",
    "day31_adjusted_slope_per_origin_sd", "day31_partial_r"
  ) %in% provenance$key))
  testthat::expect_setequal(
    list.files(bundle_path, recursive = TRUE, all.files = FALSE, no.. = TRUE),
    expected_bundle_files
  )
  testthat::expect_equal(
    sum(grepl("/tables/", expected_bundle_files)),
    10L
  )
  testthat::expect_equal(
    sum(grepl("/metadata/run_config[.]tsv$", expected_bundle_files)),
    2L
  )
  endpoint_roles <- grepl("_endpoint_ploidy$", expected_roles)
  testthat::expect_true(all(
    source_locators[match(expected_roles[endpoint_roles], sub(
      "^source_file:", "", provenance$key[startsWith(
        provenance$key,
        "source_file:"
      )]
    ))] == "Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv"
  ))
  processed_roles <- expected_roles %in% c("cellcycle", "noncellcycle")
  testthat::expect_setequal(
    source_locators[match(expected_roles[processed_roles], sub(
      "^source_file:", "", provenance$key[startsWith(
        provenance$key,
        "source_file:"
      )]
    ))],
    c(
      "Data/in-vivo/figure7/processed/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv",
      "Data/in-vivo/figure7/processed/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
    )
  )
  bundle_source_locators <- source_locators[!grepl(
    "all_ploidy[.]csv$|Data/in-vivo/figure7/processed/",
    source_locators
  )]
  testthat::expect_true(all(startsWith(
    bundle_source_locators,
    paste0(bundle_locator, "/")
  )))

  path_keys <- c(
    "assembly_script", "source_bundle", "day24_source_run",
    "day31_source_run", "png", "pdf"
  )
  declared_path_locators <- c(keyed[path_keys], source_locators)
  testthat::expect_false(any(grepl(
    "^(/|[A-Za-z]:[\\\\/])|^external:|^Results/",
    declared_path_locators
  )))
  testthat::expect_true(all(vapply(
    declared_path_locators,
    function(locator) file.exists(file.path(repo_root, locator)),
    logical(1L)
  )))

  source_role_keys <- sub(
    "^source_file:",
    "",
    provenance$key[startsWith(provenance$key, "source_file:")]
  )
  for (role in source_role_keys) {
    locator <- keyed[[paste0("source_file:", role)]]
    testthat::expect_identical(
      keyed[[paste0("source_sha256:", role)]],
      figure7_sha256(file.path(repo_root, locator)),
      info = role
    )
  }
  testthat::expect_identical(
    keyed[["assembly_script_sha256"]],
    figure7_sha256(file.path(module_dir, "assemble_tgi_sensitivity.R"))
  )
  testthat::expect_identical(
    keyed[["png_sha256"]],
    figure7_sha256(file.path(
      validation_output,
      "Figure7_Supplement.png"
    ))
  )
  testthat::expect_identical(
    keyed[["pdf_sha256"]],
    figure7_sha256(file.path(
      validation_output,
      "Figure7_Supplement.pdf"
    ))
  )

  relocated_root <- tempfile("figure7_si8_relocated_")
  dir.create(relocated_root, recursive = TRUE)
  file_locators <- unique(c(
    keyed[["assembly_script"]], source_locators,
    keyed[["png"]], keyed[["pdf"]]
  ))
  for (locator in file_locators) {
    source <- file.path(repo_root, locator)
    destination <- file.path(relocated_root, locator)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    testthat::expect_true(file.copy(source, destination, overwrite = TRUE))
  }
  for (locator in keyed[c(
    "source_bundle", "day24_source_run", "day31_source_run"
  )]) {
    testthat::expect_true(dir.exists(file.path(relocated_root, locator)))
  }
  relocated_provenance <- file.path(
    relocated_root,
    substring(
      provenance_path,
      nchar(paste0(repo_root, .Platform$file.sep)) + 1L
    )
  )
  dir.create(dirname(relocated_provenance), recursive = TRUE, showWarnings = FALSE)
  testthat::expect_true(file.copy(
    provenance_path,
    relocated_provenance,
    overwrite = TRUE
  ))
  relocated <- figure7_read_tsv(relocated_provenance, c("key", "value"))
  relocated_keyed <- stats::setNames(
    as.character(relocated$value),
    relocated$key
  )
  testthat::expect_identical(anyDuplicated(relocated$key), 0L)
  testthat::expect_true(all(vapply(
    declared_path_locators,
    function(locator) file.exists(file.path(relocated_root, locator)),
    logical(1L)
  )))
  for (role in source_role_keys) {
    locator <- relocated_keyed[[paste0("source_file:", role)]]
    testthat::expect_identical(
      relocated_keyed[[paste0("source_sha256:", role)]],
      figure7_sha256(file.path(relocated_root, locator)),
      info = paste("relocated", role)
    )
  }

  manifest_columns <- c(
    "figure", "panel", "asset_path", "source_file", "source_kind",
    "generated_by", "command", "input_data", "result_run_dir", "run_id",
    "caption_role", "asset_status", "not_regenerated_reason",
    "local_provenance_path", "citation_or_uri", "notes"
  )
  si8_manifest <- figure7_read_tsv(file.path(
    validation_output,
    "si8_manifest.tsv"
  ))
  source_manifest <- figure7_read_tsv(file.path(
    repo_root,
    "figures/Figure7_Supplement/manifest.tsv"
  ))
  testthat::expect_identical(names(si8_manifest), manifest_columns)
  testthat::expect_equal(nrow(si8_manifest), 2L)
  testthat::expect_setequal(
    si8_manifest$panel,
    c("SuppFig8", "SuppFig8_png")
  )
  testthat::expect_true(all(
    si8_manifest$local_provenance_path ==
      substring(
        provenance_path,
        nchar(paste0(repo_root, .Platform$file.sep)) + 1L
      )
  ))
  testthat::expect_true(all(file.exists(file.path(
    repo_root,
    si8_manifest$asset_path
  ))))
  testthat::expect_equal(nrow(source_manifest), 13L)
  testthat::expect_false(any(startsWith(source_manifest$panel, "SuppFig8")))
})
