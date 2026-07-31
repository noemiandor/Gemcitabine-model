testthat::test_that("panel-F compact reference contract validates and selector is audited", {
  fixture <- figure7_test_state_reference()
  reference <- figure7_validate_historical_state_reference(
    fixture$path,
    fixture$config
  )
  testthat::expect_equal(nrow(reference$pathways), 24L)
  testthat::expect_s3_class(figure7_panel_f_plot(reference$activity, fixture$config), "ggplot")
  collision_activity <- reference$activity[
    reference$activity$pathway_id %in% c(
      "H_positive_1",
      "C2:CP:REACTOME_positive_1"
    ),
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
  bad <- reference$selected; bad$selected_rank_within_direction[[1L]] <- 4L
  figure7_write_tsv(bad, file.path(fixture$path, "panel_7F_selected_pathway_gsea.tsv"))
  testthat::expect_error(figure7_validate_historical_state_reference(fixture$path, fixture$config, verify_checksums = FALSE),
                         "metadata/order|top-four")
})

testthat::test_that("tracked historical 04i reference validates exact report lineage", {
  input <- figure7_test_inputs()
  reference_path <- file.path(
    repo_root,
    input$config$state_pathways$reference_root,
    input$config$state_pathways$reference_id
  )
  reference <- figure7_validate_historical_state_reference(
    reference_path,
    input$config
  )
  provenance <- stats::setNames(as.character(reference$provenance$value), reference$provenance$key)

  testthat::expect_equal(nrow(reference$activity), 24L * 501L)
  testthat::expect_equal(nrow(reference$selected), 24L)
  testthat::expect_identical(provenance[["canonical_reference_id"]], "taoli_04i_etp2_24_day17_v1")
  testthat::expect_identical(provenance[["workflow_id"]], "binning")
  testthat::expect_identical(provenance[["model_id"]], "ETP_reference_balanced_threshold_2_24")
  testthat::expect_identical(provenance[["accumulated_interval"]], "[0.30,0.49]")
  testthat::expect_identical(
    provenance[["report_html_sha256"]],
    "b9644b1da0399043a6aba28178a2b61375b661780c4fb08a148c7724da00bfa1"
  )
  testthat::expect_identical(
    provenance[["code_revision_04i"]],
    "dc751eab928bc40f3edb063baec447fe32a69d73"
  )
  testthat::expect_false(any(c("report_html_path", "export_source_results_root") %in% names(provenance)))
  testthat::expect_identical(
    provenance[["source_primary_coverage_sha256"]],
    "9fecc5339a86cf2e072dccfb05358580fc0db25fd3a8e9f71af301a8ea423141"
  )
})

testthat::test_that("tracked reviewed v2 is the exact approved human-only FDR reference", {
  input <- figure7_test_inputs()
  reference_path <- file.path(
    repo_root,
    input$config$state_pathways$reviewed_reference_root,
    input$config$state_pathways$reviewed_reference_id
  )
  reference <- figure7_validate_reviewed_state_reference(
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
    "reviewed_human_only_frozen"
  )
  testthat::expect_true(reference$canonical_publication_allowed)
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

testthat::test_that("reviewed v2 rejects provenance or selected-table tampering", {
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
  provenance$value[provenance$key == "reviewed_on"] <- "2026-07-31"
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
    "exact approved retry7"
  )
})

testthat::test_that("canonical provenance is location-independent and always checksummed", {
  fixture <- figure7_test_state_reference()
  provenance_path <- file.path(fixture$path, "state_pathway_provenance.tsv")
  provenance <- figure7_read_tsv(provenance_path, c("key", "value"))
  provenance <- rbind(
    provenance,
    data.frame(key = "export_source_results_root", value = "/runtime/04i/results", stringsAsFactors = FALSE)
  )
  figure7_write_tsv(provenance, provenance_path)

  testthat::expect_error(
    figure7_validate_historical_state_reference(
      fixture$path,
      fixture$config
    ),
    "SHA-256 mismatch"
  )
  testthat::expect_error(
    figure7_validate_historical_state_reference(
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
  missing <- file.path(tempdir(), input$config$state_pathways$reference_id)
  testthat::expect_error(
    figure7_validate_historical_state_reference(missing, input$config),
    "Missing historical panel-7F"
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

testthat::test_that("full source panels plus the A-K composite satisfy the exact inventory", {
  input <- figure7_test_inputs(); fixture <- figure7_test_state_reference()
  reference <- figure7_validate_historical_state_reference(
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
  figure7_save_main_composite(
    list(
      A = ae$plots$A,
      B = ae$plots$C,
      C = context$plots$C,
      D = context$plots$D,
      E = context$plots$E,
      F = context$plots$F,
      G = context$plots$G,
      H = ae$plots$B,
      I = plot_f,
      J = ae$plots$D,
      K = ae$plots$E
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
  testthat::expect_length(pdfs, 6L)
  testthat::expect_length(pngs, 7L)
  testthat::expect_true(file.exists(file.path(
    out,
    "figures",
    "Figure7_reviewed_GRCh.png"
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
    paste0("--saved-state-pathway-dir=", file.path(tempdir(), input$config$state_pathways$reference_id)),
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
