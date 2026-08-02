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
    e$test$estimate,
    -0.698401019253014,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    e$test$asymptotic_p,
    0.0540069781511847,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    e$test$permutation_p_two_sided,
    0.0591269841269841,
    tolerance = 1e-12
  )
  testthat::expect_equal(e$test$n_permutations, factorial(8L))
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

testthat::test_that("panel 7E/K restores the raw mouse-level Pearson definition on the QC universe", {
  input <- figure7_test_inputs()
  e <- figure7_panel_e(input$samples, input$config)
  testthat::expect_identical(
    unique(e$data$endpoint_ploidy_source_total_cells),
    14125L
  )
  testthat::expect_identical(
    unique(e$data$endpoint_ploidy_score_universe_total_cells),
    9832L
  )
  testthat::expect_identical(
    unique(e$data$endpoint_ploidy_source_file_count),
    16L
  )
  testthat::expect_equal(sum(e$data$n_endpoint_ploidy_cells), 5335L)
  testthat::expect_equal(sum(input$samples$n_endpoint_ploidy_cells), 9832L)
  testthat::expect_identical(
    stats::setNames(
      input$samples$n_endpoint_ploidy_cells,
      input$samples$sample_id
    )[names(figure7_curated_endpoint_counts())],
    figure7_curated_endpoint_counts()
  )
  testthat::expect_false(any(c(
    "sample_mean_plot_table_endpoint_ploidy",
    "sample_median_plot_table_endpoint_ploidy",
    "n_plot_table_endpoint_ploidy_cells",
    "sample_mean_all_cbs_endpoint_ploidy",
    "sample_median_all_cbs_endpoint_ploidy",
    "n_all_cbs_endpoint_ploidy_cells",
    "n_endpoint_ploidy_cells_excluded_from_score"
  ) %in% names(input$samples)))
  testthat::expect_identical(
    e$test$permutation_mode,
    "exact_TGI_label_enumeration"
  )
  testthat::expect_identical(e$test$permutation_strata, "none")
  testthat::expect_identical(
    e$test$association_type,
    "unadjusted_mouse_level_pearson"
  )
  testthat::expect_identical(
    e$test$score_variable,
    "sample_mean_qc_passed_curated_cbs_cell_ploidy"
  )
  testthat::expect_identical(e$test$score_source_n_cells, 9832L)
  testthat::expect_identical(e$test$score_inventory_n_cells, 14125L)
  testthat::expect_identical(e$test$score_source_n_files, 16L)
  testthat::expect_identical(e$test$treated_score_n_cells, 5335L)
  testthat::expect_identical(e$test$score_standardization, "none")
  testthat::expect_identical(e$test$adjustment_terms, "none")
  testthat::expect_identical(e$test$outcome_variable, "TGI_percent_Day_17")
  testthat::expect_identical(e$test$plot_x, "sample_mean_endpoint_ploidy")
  testthat::expect_identical(e$test$plot_y, "TGI_percent_Day_17")
  testthat::expect_equal(e$test$n_permutations, factorial(8L))
  testthat::expect_false(any(c(
    "terminal_postprocessed_cn_score",
    "terminal_cn_score_within_origin_z",
    "terminal_cn_score_nuisance_residual",
    "tgi_origin_dose_residual",
    "permutation_stratum"
  ) %in% names(e$data)))
  testthat::expect_false(any(grepl("etp_group", names(e$data), fixed = TRUE)))

  plot <- figure7_endpoint_ploidy_plot(e$data, e$test, input$config)
  testthat::expect_match(plot$labels$title, "mean endpoint tumor-cell ploidy", fixed = TRUE)
  testthat::expect_match(plot$labels$subtitle, "unadjusted mouse-level", fixed = TRUE)
  testthat::expect_identical(plot$labels$x, "Mean endpoint tumor-cell ploidy")
  testthat::expect_identical(plot$labels$y, "Day 17 TGI (%)")
  testthat::expect_match(deparse(plot$mapping$shape), "initial_ploidy", fixed = TRUE)
  testthat::expect_match(deparse(plot$mapping$colour), "dose", fixed = TRUE)
})

testthat::test_that("panel 7E/K exposes the frozen Day-17/24/31 sensitivity results", {
  expected <- data.frame(
    day = c(17L, 24L, 31L),
    estimate = c(
      -0.698401019253014,
      -0.387348978978662,
      -0.296948455789288
    ),
    asymptotic_p = c(
      0.0540069781511847,
      0.343097626516948,
      0.475086351653848
    ),
    permutation_p = c(
      0.0591269841269841,
      0.346924603174603,
      0.485714285714286
    )
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
      result$estimate,
      expected$estimate[[i]],
      tolerance = 1e-12
    )
    testthat::expect_equal(
      result$asymptotic_p,
      expected$asymptotic_p[[i]],
      tolerance = 1e-12
    )
    testthat::expect_equal(
      result$permutation_p_two_sided,
      expected$permutation_p[[i]],
      tolerance = 1e-12
    )
    testthat::expect_equal(result$n_permutations, factorial(8L))
    testthat::expect_identical(
      result$outcome_variable,
      paste0("TGI_percent_Day_", day)
    )
  }
})

testthat::test_that("selected ECDF panel IDs and linetypes are frozen", {
  input <- figure7_test_inputs()
  panel <- figure7_panel_b(input$data, input$samples, input$config)
  testthat::expect_identical(unique(panel$tests$comparison_id), c(1L, 8L, 9L))
  testthat::expect_identical(unique(panel$tests$panel),
    c("1. 0 vs treated", "8. 4N: 0 vs treated", "9. 2N: 0 vs treated"))
  testthat::expect_identical(unique(panel$data$line_group[panel$data$comparison_id == 8L]), "4N")
  testthat::expect_identical(unique(panel$data$line_group[panel$data$comparison_id == 9L]), "2N")
})

testthat::test_that("treated-cell excess is localized with equal-mouse exact inference", {
  input <- figure7_test_inputs()
  result <- figure7_density_localization(
    input$data, input$samples, input$config
  )
  test <- result$test

  testthat::expect_identical(test$cell_universe, "reviewed_qc_retained_cellcycle_2881")
  testthat::expect_identical(test$n_cells, 2881L)
  testthat::expect_identical(test$n_samples, 16L)
  testthat::expect_identical(test$n_vehicle_samples, 8L)
  testthat::expect_identical(test$n_treated_samples, 8L)
  testthat::expect_identical(test$sample_weighting, "equal_mouse")
  testthat::expect_identical(test$permutation_strata, "initial_ploidy")
  testthat::expect_identical(test$n_permutations, 4900L)
  testthat::expect_identical(test$grid_points, 501L)
  testthat::expect_equal(
    test$bandwidth, 0.0506470455660707, tolerance = 1e-13
  )
  testthat::expect_equal(test$pointwise_start, 0.296, tolerance = 1e-12)
  testthat::expect_equal(test$pointwise_end, 0.486, tolerance = 1e-12)
  testthat::expect_equal(test$simultaneous_start, 0.414, tolerance = 1e-12)
  testthat::expect_equal(test$simultaneous_end, 0.426, tolerance = 1e-12)
  testthat::expect_equal(
    test$simultaneous_critical,
    2.6737373166169203,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    test$raw_excess_peak_pseudotime, 0.452, tolerance = 1e-12
  )
  testthat::expect_equal(
    test$raw_excess_peak_density_difference,
    0.283948126904316,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    test$max_abs_t_pseudotime, 0.420, tolerance = 1e-12
  )
  testthat::expect_equal(
    test$max_abs_t_observed_statistic,
    2.68138128387289,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    test$global_max_t_p_two_sided,
    232 / 4900,
    tolerance = 1e-12
  )
  testthat::expect_identical(
    result$intervals$support_type,
    c("positive_pointwise_two_sided", "positive_simultaneous_max_abs_t")
  )
  ordered_samples <- input$samples[order(input$samples$sample_id), , drop = FALSE]
  density_matrix <- t(vapply(ordered_samples$sample_id, function(id) {
    stats::density(
      input$data$pseudotime[input$data$sample_id == id],
      bw = test$bandwidth, from = 0, to = 1, n = 501
    )$y
  }, numeric(501L)))
  vehicle <- ordered_samples$dose_mg == 0
  treated <- !vehicle
  testthat::expect_equal(
    result$grid$vehicle_mean_density,
    colMeans(density_matrix[vehicle, , drop = FALSE]),
    tolerance = 1e-13
  )
  testthat::expect_equal(
    result$grid$treated_mean_density,
    colMeans(density_matrix[treated, , drop = FALSE]),
    tolerance = 1e-13
  )
  cell_weighted_vehicle <- stats::density(
    input$data$pseudotime[
      input$data$sample_id %in% ordered_samples$sample_id[vehicle]
    ],
    bw = test$bandwidth, from = 0, to = 1, n = 501
  )$y
  testthat::expect_gt(
    max(abs(result$grid$vehicle_mean_density - cell_weighted_vehicle)),
    1e-3
  )

  set.seed(71)
  shuffled <- figure7_density_localization(
    input$data[sample(seq_len(nrow(input$data))), , drop = FALSE],
    input$samples[sample(seq_len(nrow(input$samples))), , drop = FALSE],
    input$config
  )
  testthat::expect_equal(
    shuffled$grid,
    result$grid,
    tolerance = 1e-13,
    check.attributes = FALSE
  )
  testthat::expect_equal(
    shuffled$test,
    result$test,
    tolerance = 1e-13,
    check.attributes = FALSE
  )
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
  top_level <- "run_figure7.R"
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
  testthat::expect_true(all(file.exists(file.path(
    out,
    "tables",
    c(
      "panel_7B_density_localization_grid.tsv",
      "panel_7B_density_localization_intervals.tsv",
      "panel_7B_density_localization_test.tsv"
    )
  ))))
})
