testthat::test_that("A-E frozen numerical results and treated scope are reproduced", {
  input <- figure7_test_inputs()
  testthat::expect_equal(nrow(input$samples[input$samples$dose_mg > 0, ]), 8L)
  testthat::expect_setequal(input$samples$dose_mg[input$samples$dose_mg > 0], c(30, 120))
  c_test <- figure7_panel_c_test(input$samples)
  testthat::expect_equal(c_test$dose_adjusted_difference_high_minus_low, -39.5232168, tolerance = 1e-6)
  testthat::expect_equal(c_test$permutation_p_two_sided, 0.0555556, tolerance = 1e-7)
  testthat::expect_equal(c_test$n_permutations, 36L)
  d <- figure7_panel_d(figure7_shift_metrics(input$data, input$samples))
  testthat::expect_equal(d$test$estimate, 0.8105173, tolerance = 1e-7)
  testthat::expect_equal(d$test$permutation_p_two_sided, 0.0104167, tolerance = 1e-7)
  testthat::expect_equal(d$test$n_permutations, 576L)
  e <- figure7_panel_e(input$samples)
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

testthat::test_that("all module R files parse and A-E builders emit five PDF/PNG pairs", {
  top_level <- c("run_figure7.R", "export_04i_state_pathway_reference.R")
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
