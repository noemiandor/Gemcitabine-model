testthat::test_that("panel-F compact reference contract validates and selector is audited", {
  fixture <- figure7_test_state_reference()
  reference <- figure7_validate_state_reference(fixture$path, fixture$config)
  testthat::expect_equal(nrow(reference$pathways), 24L)
  testthat::expect_s3_class(figure7_panel_f_plot(reference$activity, fixture$config), "ggplot")
  bad <- reference$selected; bad$selected_rank_within_direction[[1L]] <- 4L
  figure7_write_tsv(bad, file.path(fixture$path, "panel_7F_selected_pathway_gsea.tsv"))
  testthat::expect_error(figure7_validate_state_reference(fixture$path, fixture$config, verify_checksums = FALSE),
                         "metadata/order|top-four")
})

testthat::test_that("missing F, wrong checksums, and nonempty outputs fail clearly", {
  input <- figure7_test_inputs()
  missing <- file.path(tempdir(), input$config$state_pathways$reference_id)
  testthat::expect_error(figure7_validate_state_reference(missing, input$config), "Missing canonical panel-7F")
  wrong <- tempfile(); writeLines("wrong", wrong)
  testthat::expect_error(figure7_verify_checksum(wrong, paste(rep("0", 64), collapse = "")), "SHA-256 mismatch")
  nonempty <- tempfile(); dir.create(nonempty); writeLines("x", file.path(nonempty, "existing.txt"))
  testthat::expect_error(figure7_prepare_output(nonempty), "pre-existing analysis files")
})

testthat::test_that("entire-run image inventory is exact", {
  config <- figure7_test_inputs()$config; out <- tempfile(); dir.create(file.path(out, "figures"), recursive = TRUE)
  for (file in figure7_panel_filenames(config)) writeLines("%PDF fixture", file.path(out, "figures", file))
  testthat::expect_silent(figure7_validate_figure_inventory(out, config))
  dir.create(file.path(out, "tables")); file.create(file.path(out, "tables", "unexpected.png"))
  testthat::expect_error(figure7_validate_figure_inventory(out, config), "Unexpected non-PDF")
})

testthat::test_that("A-E plus a validated fixture F satisfy the exact six-panel inventory", {
  input <- figure7_test_inputs(); fixture <- figure7_test_state_reference()
  reference <- figure7_validate_state_reference(fixture$path, fixture$config)
  out <- tempfile("figure7_six_"); figure7_prepare_output(out)
  figure7_build_ae(input$cellcycle, input$data, input$samples, out, fixture$config)
  figure7_build_f(reference, out, fixture$config)
  testthat::expect_silent(figure7_validate_figure_inventory(out, fixture$config))
  pdfs <- list.files(file.path(out, "figures"), pattern = "[.]pdf$", full.names = TRUE)
  testthat::expect_length(pdfs, 6L)
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
