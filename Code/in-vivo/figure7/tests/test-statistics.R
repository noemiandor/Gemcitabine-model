testthat::test_that("A-E frozen numerical results and treated scope are reproduced", {
  input <- figure7_test_inputs()
  testthat::expect_equal(nrow(input$samples[input$samples$dose_mg > 0, ]), 8L)
  testthat::expect_setequal(input$samples$dose_mg[input$samples$dose_mg > 0], c(30, 120))
  c_test <- figure7_panel_c_test(input$samples, input$config)
  testthat::expect_equal(c_test$dose_adjusted_difference_high_minus_low, -39.5232168, tolerance = 1e-6)
  testthat::expect_equal(c_test$permutation_p_two_sided, 0.0555556, tolerance = 1e-7)
  testthat::expect_equal(c_test$n_permutations, 36L)
  d <- figure7_panel_d(figure7_shift_metrics(input$data, input$samples, input$config), input$config)
  testthat::expect_equal(d$test$estimate, 0.7399454530471112, tolerance = 1e-12)
  testthat::expect_equal(d$test$permutation_p_two_sided, 0.01736111111111111, tolerance = 1e-12)
  testthat::expect_equal(d$test$n_permutations, 576L)
  testthat::expect_identical(d$test$ecdf_reference_group_column, "initial_ploidy")
  e <- figure7_panel_e(input$samples, input$config)
  testthat::expect_equal(
    e$test$effect_per_within_origin_sd,
    -8.8374975584089,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    e$test$partial_correlation,
    -0.344659307863932,
    tolerance = 1e-12
  )
  testthat::expect_equal(e$test$estimate, e$test$partial_correlation)
  testthat::expect_equal(e$test$permutation_p_two_sided, 0.6875, tolerance = 1e-12)
  testthat::expect_equal(e$test$n_permutations, 16L)
})

testthat::test_that("panel 7D/J uses only injected-origin-matched untreated ECDF references", {
  input <- figure7_test_inputs()
  shifts <- figure7_shift_metrics(input$data, input$samples, input$config)
  testthat::expect_true(all(shifts$reference_group_column == "initial_ploidy"))
  testthat::expect_true(all(shifts$reference_group_value == shifts$initial_ploidy))
  for (i in seq_len(nrow(shifts))) {
    reference_ids <- strsplit(shifts$reference_sample_ids[[i]], ";", fixed = TRUE)[[1L]]
    reference <- input$samples[
      match(reference_ids, input$samples$sample_id),
      ,
      drop = FALSE
    ]
    testthat::expect_true(all(reference$dose_mg == 0))
    testthat::expect_true(all(reference$initial_ploidy == shifts$initial_ploidy[[i]]))
  }
  testthat::expect_error(
    figure7_shift_metrics(
      input$data,
      input$samples,
      input$config,
      group_column = "etp_group"
    ),
    "endpoint-CN-score groups are not permitted"
  )
})

testthat::test_that("panel 7E/K is standardized within origin and adjusted without endpoint grouping", {
  input <- figure7_test_inputs()
  e <- figure7_panel_e(input$samples, input$config)
  standardized <- split(
    e$data$terminal_cn_score_within_origin_z,
    e$data$initial_ploidy
  )
  testthat::expect_true(all(vapply(
    standardized,
    function(values) abs(mean(values)) <= 1e-12 &&
      abs(stats::sd(values) - 1) <= 1e-12,
    logical(1L)
  )))
  testthat::expect_equal(
    e$data$terminal_postprocessed_cn_score,
    e$data$sample_mean_endpoint_ploidy,
    tolerance = 1e-12
  )
  testthat::expect_identical(
    unique(e$data$endpoint_ploidy_source_total_cells),
    14125L
  )
  testthat::expect_identical(
    unique(e$data$endpoint_ploidy_source_file_count),
    16L
  )
  testthat::expect_equal(sum(e$data$n_endpoint_ploidy_cells), 7623L)
  testthat::expect_equal(
    sum(input$samples$n_endpoint_ploidy_cells),
    14125L
  )
  testthat::expect_equal(
    sum(input$samples$n_plot_table_endpoint_ploidy_cells),
    9832L
  )
  testthat::expect_true(any(
    abs(
      input$samples$sample_mean_endpoint_ploidy -
        input$samples$sample_mean_plot_table_endpoint_ploidy
    ) > 1e-6
  ))
  testthat::expect_identical(
    e$data$permutation_stratum,
    paste(e$data$initial_ploidy, e$data$dose_mg, sep = "|")
  )
  testthat::expect_identical(
    e$test$permutation_mode,
    "exact_TGI_label_enumeration_within_initial_ploidy_x_dose"
  )
  testthat::expect_identical(e$test$permutation_strata, "initial_ploidy:dose_mg")
  testthat::expect_identical(
    e$test$score_variable,
    "sample_mean_all_canonical_cbs_cell_ploidy"
  )
  testthat::expect_identical(e$test$score_source_n_cells, 14125L)
  testthat::expect_identical(e$test$score_source_n_files, 16L)
  testthat::expect_identical(e$test$treated_score_n_cells, 7623L)
  testthat::expect_identical(e$test$score_standardization, "z_score_within_initial_ploidy")
  testthat::expect_identical(e$test$adjustment_terms, "initial_ploidy+dose_mg")
  testthat::expect_identical(e$test$outcome_variable, "TGI_percent_Day_17")
  testthat::expect_false(any(grepl("etp_group", names(e$data), fixed = TRUE)))

  plot <- figure7_adjusted_cn_plot(e$data, e$test, input$config)
  testthat::expect_match(plot$labels$title, "all-cell terminal postprocessed CN score", fixed = TRUE)
  testthat::expect_match(plot$labels$subtitle, "adjusted for injected origin", fixed = TRUE)
  testthat::expect_match(plot$labels$x, "within-origin z score", fixed = TRUE)
  testthat::expect_match(plot$labels$y, "origin- and dose-adjusted", fixed = TRUE)
  testthat::expect_match(deparse(plot$mapping$shape), "initial_ploidy", fixed = TRUE)
  testthat::expect_match(deparse(plot$mapping$colour), "dose", fixed = TRUE)
})

testthat::test_that("panel 7E/K exposes the frozen Day-17/24/31 sensitivity results", {
  expected <- data.frame(
    day = c(17L, 24L, 31L),
    effect = c(
      -8.8374975584089,
      -6.77351906237078,
      -2.77364130329883
    ),
    partial_correlation = c(
      -0.344659307863932,
      -0.30605912538695,
      -0.173200987198346
    ),
    permutation_p = c(0.6875, 0.5625, 0.6875)
  )
  for (i in seq_len(nrow(expected))) {
    day <- expected$day[[i]]
    config <- figure7_read_config(
      file.path(module_dir, "figure7_config.yaml"),
      day
    )
    cellcycle <- figure7_read_cell_table(file.path(
      repo_root, "Data/in-vivo/figure7/processed",
      "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
    ), "CellCycle", config)
    noncellcycle <- figure7_read_cell_table(file.path(
      repo_root, "Data/in-vivo/figure7/processed",
      "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
    ), "NonCellCycle", config)
    endpoint_ploidy <- figure7_read_endpoint_ploidy_table(
      file.path(repo_root, "Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv"),
      config
    )
    samples <- figure7_sample_table(
      cellcycle,
      noncellcycle,
      config,
      endpoint_ploidy
    )
    result <- figure7_panel_e(samples, config)$test
    testthat::expect_equal(
      result$effect_per_within_origin_sd,
      expected$effect[[i]],
      tolerance = 1e-12
    )
    testthat::expect_equal(
      result$partial_correlation,
      expected$partial_correlation[[i]],
      tolerance = 1e-12
    )
    testthat::expect_equal(
      result$permutation_p_two_sided,
      expected$permutation_p[[i]],
      tolerance = 1e-12
    )
    testthat::expect_equal(result$n_permutations, 16L)
    testthat::expect_identical(
      result$outcome_variable,
      paste0("TGI_percent_Day_", day)
    )
  }
})

testthat::test_that("selected ECDF panel IDs and linetypes are frozen", {
  input <- figure7_test_inputs(); panel <- figure7_panel_b(input$data, input$samples)
  testthat::expect_identical(unique(panel$tests$comparison_id), c(1L, 8L, 9L))
  testthat::expect_identical(unique(panel$tests$panel),
    c("1. 0 vs treated", "8. 4N: 0 vs treated", "9. 2N: 0 vs treated"))
  testthat::expect_identical(unique(panel$data$line_group[panel$data$comparison_id == 8L]), "4N")
  testthat::expect_identical(unique(panel$data$line_group[panel$data$comparison_id == 9L]), "2N")
})

testthat::test_that("panel 7A uses treatment-relative time with the Day_0 baseline", {
  input <- figure7_test_inputs()
  growth <- figure7_growth_trajectory_data(input$cellcycle, input$samples, input$config)
  panel <- figure7_panel_a_plot(growth, input$config)
  individual <- growth[growth$series_type == "individual mouse", , drop = FALSE]

  testthat::expect_identical(panel$labels$x, "Days since first treatment")
  testthat::expect_identical(unique(input$cellcycle$tumor_volume_baseline_day), "Day_0")
  testthat::expect_equal(sum(individual$day == 0), length(unique(individual$sample_id)))
  testthat::expect_true(all(individual$tumor_volume_change[individual$day == 0] == 0))
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
  endpoint_ploidy <- figure7_read_endpoint_ploidy_table(
    file.path(repo_root, "Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv"),
    config
  )
  samples <- figure7_sample_table(
    cellcycle,
    noncellcycle,
    config,
    endpoint_ploidy
  )
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

testthat::test_that("state-pathway exporter requires an explicit results root", {
  exporter <- file.path(module_dir, "export_04i_state_pathway_reference.R")
  output_parent <- tempfile("figure7_export_requires_root_")
  output_dir <- file.path(output_parent, "taoli_04i_etp2_24_day17_v1")
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
