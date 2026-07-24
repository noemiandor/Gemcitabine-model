testthat::test_that("A-E frozen numerical results and treated scope are reproduced", {
  input <- figure7_test_inputs()
  testthat::expect_equal(nrow(input$samples[input$samples$dose_mg > 0, ]), 8L)
  testthat::expect_setequal(input$samples$dose_mg[input$samples$dose_mg > 0], c(30, 120))
  c_test <- figure7_panel_c_test(input$samples, input$config)
  testthat::expect_equal(c_test$dose_adjusted_difference_high_minus_low, -39.5232168, tolerance = 1e-6)
  testthat::expect_equal(c_test$permutation_p_two_sided, 0.0555556, tolerance = 1e-7)
  testthat::expect_equal(c_test$n_permutations, 36L)
  d <- figure7_panel_d(figure7_shift_metrics(input$data, input$samples, input$config), input$config)
  testthat::expect_equal(d$test$estimate, 0.8105173, tolerance = 1e-7)
  testthat::expect_equal(d$test$permutation_p_two_sided, 0.0104167, tolerance = 1e-7)
  testthat::expect_equal(d$test$n_permutations, 576L)
  e <- figure7_panel_e(input$samples, input$config)
  testthat::expect_equal(e$test$estimate, -0.6984010, tolerance = 1e-7)
  testthat::expect_equal(e$test$permutation_p_two_sided, 0.0591270, tolerance = 1e-7)
  testthat::expect_equal(e$test$n_permutations, 40320L)
})

testthat::test_that("selected ECDF panel IDs and linetypes are frozen", {
  input <- figure7_test_inputs(); panel <- figure7_panel_b(input$data, input$samples)
  testthat::expect_identical(unique(panel$tests$comparison_id), c(1L, 8L, 9L))
  testthat::expect_identical(unique(panel$tests$panel),
    c("1. 0 vs treated", "8. 4N: 0 vs treated", "9. 2N: 0 vs treated"))
  testthat::expect_identical(unique(panel$data$line_group[panel$data$comparison_id == 8L]), "4N")
  testthat::expect_identical(unique(panel$data$line_group[panel$data$comparison_id == 9L]), "2N")
})

testthat::test_that("Day-24 override recomputes TGI and emits an explicit Day-24 contract", {
  config <- figure7_read_config(file.path(module_dir, "figure7_config.yaml"), 24L)
  cellcycle <- figure7_read_cell_table(file.path(
    repo_root, "Data/in-vivo/figure7/processed",
    "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
  ), "CellCycle", config)
  noncellcycle <- figure7_read_cell_table(file.path(
    repo_root, "Data/in-vivo/figure7/processed",
    "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
  ), "NonCellCycle", config)
  samples <- figure7_sample_table(cellcycle, noncellcycle, config)
  treated <- samples$dose_mg > 0

  testthat::expect_identical(figure7_tgi_measure(config), "TGI_percent_Day_24")
  testthat::expect_identical(unique(samples$tgi_day), 24L)
  testthat::expect_equal(
    samples$TGI_percent_Day_24[treated],
    samples$embedded_TGI_percent_Day_24[treated],
    tolerance = as.numeric(config$tgi$numerical_tolerance)
  )
  testthat::expect_true(all(is.finite(figure7_panel_c_test(samples, config)$permutation_p_two_sided)))
  testthat::expect_identical(
    figure7_panel_filenames(config, c("7A", "7B", "7C", "7D", "7E")),
    c(
      "panel_7A_day24_tgi_calculation.pdf",
      "panel_7B_cellcycle_selected_ecdf_comparisons.pdf",
      "panel_7C_day24_tgi_by_initial_ploidy.pdf",
      "panel_7D_day24_tgi_vs_centered_ecdf_shift.pdf",
      "panel_7E_day24_tgi_vs_mean_etp.pdf"
    )
  )
})

testthat::test_that("all module R files parse and A-E builders emit five PDF/PNG pairs", {
  top_level <- c("run_figure7.R", "download_figure7_raw_data.R", "export_state_pathway_reference.R")
  for (file in c(top_level, list.files(file.path(module_dir, "src"), pattern = "[.]R$", full.names = FALSE))) {
    path <- if (file %in% top_level) file.path(module_dir, file) else file.path(module_dir, "src", file)
    testthat::expect_silent(parse(file = path))
  }
  input <- figure7_test_inputs(); out <- tempfile("figure7_ae_"); figure7_prepare_output(out)
  result <- figure7_build_ae(input$cellcycle, input$data, input$samples, out, input$config)
  testthat::expect_length(list.files(file.path(out, "figures"), pattern = "[.]pdf$"), 5L)
  testthat::expect_length(list.files(file.path(out, "figures"), pattern = "[.]png$"), 5L)
  testthat::expect_silent(figure7_validate_figure_inventory(out, input$config, figure7_panel_ids(FALSE)))
  testthat::expect_identical(unique(result$panel_b$tests$comparison_id), c(1L, 8L, 9L))
})

testthat::test_that("state-pathway exporter requires an explicit results root", {
  exporter <- file.path(module_dir, "export_state_pathway_reference.R")
  output_parent <- tempfile("figure7_export_requires_root_")
  output_dir <- file.path(output_parent, "taoli_state_pathway_etp2_24_day17_v1")
  status <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(exporter, paste0("--output-dir=", output_dir)),
    stdout = TRUE,
    stderr = TRUE
  ))
  testthat::expect_true(!is.null(attr(status, "status")) && attr(status, "status") != 0L)
  testthat::expect_match(paste(status, collapse = "\n"), "Missing required argument --results-root")
  testthat::expect_false(dir.exists(output_dir))
})
