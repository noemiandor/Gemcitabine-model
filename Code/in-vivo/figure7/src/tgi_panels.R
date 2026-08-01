# Canonical A-E panel builders, adapted from the reviewed 04h plotting layer.

figure7_theme <- function() {
  ggplot2::theme_bw(base_size = 11) + ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey92", color = "grey70"),
    legend.position = "bottom"
  )
}

figure7_dose_colors <- function() c(
  "0mg/kg" = "#666666", "30mg/kg" = "#d95f02", "120mg/kg" = "#1b9e77",
  "treated" = "#377eb8"
)

figure7_add_metadata <- function(data, config) {
  figure7_add_tgi_metadata(data, config)
}

figure7_growth_trajectory_data <- function(cellcycle, samples, config) {
  tgi_day <- figure7_tgi_day(config)
  columns <- grep("^tumor_volume_Day_[0-9]+$", names(cellcycle), value = TRUE)
  days <- figure7_numeric(sub("^tumor_volume_Day_", "", columns))
  columns <- columns[order(days)]; days <- sort(days)
  rows <- lapply(samples$sample_id, function(id) {
    local <- cellcycle[cellcycle$sample_id == id, , drop = FALSE]
    baseline_day <- figure7_numeric(sub("^Day_", "", figure7_unique_sample_value(local, "tumor_volume_baseline_day", id)))
    baseline <- figure7_unique_sample_value(local, "tumor_volume_baseline", id, TRUE)
    volumes <- vapply(columns, function(column) figure7_unique_sample_value(local, column, id, TRUE), numeric(1L))
    meta <- samples[samples$sample_id == id, , drop = FALSE]
    data.frame(
      sample_id = id, initial_ploidy = meta$initial_ploidy, dose = meta$dose, dose_mg = meta$dose_mg,
      day = c(baseline_day, days), tumor_volume = c(baseline, volumes),
      tumor_volume_change = c(0, volumes - baseline), series_type = "individual mouse",
      reference_n = NA_integer_, highlight_day = tgi_day, is_highlight_day = c(baseline_day, days) == tgi_day,
      stringsAsFactors = FALSE
    )
  })
  individual <- do.call(rbind, rows)
  controls <- individual[individual$dose_mg == 0, , drop = FALSE]
  references <- lapply(split(controls, interaction(controls$initial_ploidy, controls$day, drop = TRUE)), function(local) {
    data.frame(
      sample_id = paste0(local$initial_ploidy[[1L]], " matched-control reference"),
      initial_ploidy = local$initial_ploidy[[1L]], dose = "matched-control reference", dose_mg = 0,
      day = local$day[[1L]], tumor_volume = NA_real_, tumor_volume_change = mean(local$tumor_volume_change),
      series_type = "matched-control reference", reference_n = nrow(local), highlight_day = tgi_day,
      is_highlight_day = local$day[[1L]] == tgi_day, stringsAsFactors = FALSE
    )
  })
  figure7_add_metadata(rbind(individual, do.call(rbind, references)), config)
}

figure7_panel_a_plot <- function(data, config) {
  tgi_day <- figure7_tgi_day(config)
  individual <- data[data$series_type == "individual mouse", , drop = FALSE]
  reference <- data[data$series_type == "matched-control reference", , drop = FALSE]
  individual$dose <- factor(individual$dose, levels = c("0mg/kg", "30mg/kg", "120mg/kg"))
  ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = tgi_day, color = "#B2182B", linewidth = 0.9, alpha = 0.8) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::geom_line(data = individual,
      ggplot2::aes(day, tumor_volume_change, group = sample_id, color = dose), linewidth = 0.55, alpha = 0.72) +
    ggplot2::geom_line(data = reference,
      ggplot2::aes(day, tumor_volume_change, group = initial_ploidy), color = "black", linetype = "22", linewidth = 1.15) +
    ggplot2::geom_point(data = individual[individual$is_highlight_day, ],
      ggplot2::aes(day, tumor_volume_change, color = dose), shape = 21, fill = "white", stroke = 1, size = 2.5) +
    ggplot2::geom_point(data = reference[reference$is_highlight_day, ],
      ggplot2::aes(day, tumor_volume_change), shape = 23, fill = "#B2182B", color = "black", size = 3.1) +
    ggplot2::facet_wrap(~initial_ploidy, nrow = 1) +
    ggplot2::scale_color_manual(values = figure7_dose_colors(), breaks = c("0mg/kg", "30mg/kg", "120mg/kg"), name = "Dose") +
    ggplot2::scale_x_continuous(breaks = sort(unique(individual$day)), minor_breaks = NULL) +
    ggplot2::labs(
      title = paste("How Day", tgi_day, "TGI is calculated from tumor growth"),
      subtitle = paste0("TGI = 100 x (1 - treated Day-", tgi_day,
                        " growth delta / mean matched-control Day-", tgi_day, " growth delta)"),
      x = "Days since first treatment", y = expression(Delta * " tumor volume from Day 0 (mm"^3 * ")"),
      caption = paste0("Thin lines: individual mice. Dashed black: mean initial-ploidy-matched untreated controls. Red: Day ",
                       tgi_day, ".")
    ) + figure7_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), plot.caption = ggplot2::element_text(hjust = 0))
}

figure7_panel_b_plot <- function(data, tests) {
  levels <- c("1. 0 vs treated", "8. 4N: 0 vs treated", "9. 2N: 0 vs treated")
  data$panel <- factor(data$panel, levels = levels); tests$panel <- factor(tests$panel, levels = levels)
  data$line_group <- factor(data$line_group, levels = c("All", "2N", "4N"))
  ggplot2::ggplot(data, ggplot2::aes(pseudotime, mean_ecdf, color = color_group,
                                     linetype = line_group, group = curve_label)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_text(data = tests, ggplot2::aes(x = 0.98, y = 0.05, label = annotation),
      inherit.aes = FALSE, hjust = 1, vjust = 0, size = 2.35, lineheight = 0.92) +
    ggplot2::facet_wrap(~panel, ncol = 3, drop = FALSE) +
    ggplot2::scale_color_manual(values = figure7_dose_colors(), name = "Group") +
    ggplot2::scale_linetype_manual(values = c("All" = "solid", "2N" = "solid", "4N" = "22"),
                                   name = "Initial ploidy") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.05), clip = "off") +
    ggplot2::labs(title = "CellCycle: selected direct group mean ECDF comparisons (initial_ploidy)",
                  subtitle = "Equal-sample mean ECDFs", x = "Pseudotime", y = "Mean ECDF") + figure7_theme()
}

figure7_panel_c_plot <- function(data, test, config) {
  tgi_day <- figure7_tgi_day(config)
  tgi_measure <- figure7_tgi_measure(config)
  data$initial_ploidy <- factor(data$initial_ploidy, levels = c("2N", "4N"))
  data$dose <- factor(data$dose, levels = c("30mg/kg", "120mg/kg"))
  label <- sprintf(paste0(
    "Independent tumors; no one-to-one mouse pairing\nExact dose-stratified permutation\n",
    "Adjusted delta (4N - 2N) = %.3f percentage points; P = %.3g\n",
    "n = %d (2N) and %d (4N)"),
    test$dose_adjusted_difference_high_minus_low, test$permutation_p_two_sided,
    test$n_group_low, test$n_group_high)
  ggplot2::ggplot(data, ggplot2::aes(initial_ploidy, .data[[tgi_measure]], fill = initial_ploidy)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::geom_boxplot(width = 0.58, alpha = 0.35, outlier.shape = NA) +
    ggplot2::geom_point(ggplot2::aes(color = dose, shape = dose),
      position = ggplot2::position_jitter(width = 0.075, height = 0, seed = 1), size = 3.1) +
    ggplot2::scale_fill_manual(values = c("2N" = "#4C78A8", "4N" = "#E45756"), guide = "none") +
    ggplot2::scale_color_manual(values = figure7_dose_colors(), breaks = c("30mg/kg", "120mg/kg"), name = "Dose") +
    ggplot2::scale_shape_manual(values = c("30mg/kg" = 16, "120mg/kg" = 17),
                                breaks = c("30mg/kg", "120mg/kg"), name = "Dose") +
    ggplot2::labs(title = paste("Treated CellCycle tumors: Day", tgi_day, "TGI by Initial ploidy"),
                  subtitle = label, x = "Initial ploidy", y = paste("Day", tgi_day, "TGI (%)")) + figure7_theme() +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9, color = "grey25", lineheight = 1.05))
}

figure7_scatter_plot <- function(
  data,
  x,
  y,
  test,
  title,
  x_label,
  y_label,
  annotation_p_label = "Permutation P"
) {
  data$initial_ploidy <- factor(data$initial_ploidy, levels = c("2N", "4N"))
  data$dose <- factor(data$dose, levels = c("30mg/kg", "120mg/kg"))
  annotation <- sprintf(
    "Pearson r = %.3f\n%s = %.3g\nn = %d",
    test$estimate,
    annotation_p_label,
    test$permutation_p_two_sided,
    test$n
  )
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[[x]], y = .data[[y]], color = dose,
      shape = initial_ploidy
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::geom_smooth(data = data, ggplot2::aes(x = .data[[x]], y = .data[[y]], group = 1),
                         inherit.aes = FALSE, method = "lm", se = TRUE, color = "black", linewidth = 0.55) +
    ggplot2::geom_point(size = 2.9) +
    ggrepel::geom_text_repel(ggplot2::aes(label = sample_id), size = 2.4, seed = 1,
      min.segment.length = 0, segment.color = "grey65", max.overlaps = Inf, show.legend = FALSE) +
    ggplot2::annotate("label", x = Inf, y = Inf, label = annotation, hjust = 1.05, vjust = 1.1, size = 3) +
    ggplot2::scale_color_manual(values = figure7_dose_colors(), name = "Dose") +
    ggplot2::scale_shape_manual(
      values = c("2N" = 16, "4N" = 17),
      name = "Injected origin"
    ) +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = y_label,
      shape = "Injected origin"
    ) +
    ggplot2::coord_cartesian(clip = "off") + figure7_theme() +
    ggplot2::theme(legend.position = "right")
}

figure7_adjusted_cn_plot <- function(data, test, config) {
  required_data <- c(
    "sample_id", "initial_ploidy", "dose",
    "endpoint_ploidy_file", "n_endpoint_ploidy_cells",
    "endpoint_ploidy_source_total_cells",
    "terminal_cn_score_nuisance_residual", "tgi_origin_dose_residual"
  )
  required_test <- c(
    "n", "partial_correlation", "effect_per_within_origin_sd",
    "permutation_p_two_sided", "permutation_strata",
    "score_variable", "score_source_n_cells", "score_source_n_files",
    "treated_score_n_cells", "score_aggregation_policy",
    "sample_mapping_policy", "score_standardization", "adjustment_terms"
  )
  if (length(setdiff(required_data, names(data))) ||
      length(setdiff(required_test, names(test)))) {
    figure7_stop("Panel 7E/K adjusted CN-score plot received an incomplete analysis contract")
  }
  if (!identical(as.character(test$permutation_strata[[1L]]), "initial_ploidy:dose_mg") ||
      !identical(
        as.character(test$score_variable[[1L]]),
        "sample_mean_all_canonical_cbs_cell_ploidy"
      ) ||
      as.integer(test$score_source_n_cells[[1L]]) != 14125L ||
      as.integer(test$score_source_n_files[[1L]]) != 16L ||
      as.integer(test$treated_score_n_cells[[1L]]) != 7623L ||
      !identical(as.character(test$score_standardization[[1L]]), "z_score_within_initial_ploidy") ||
      !identical(as.character(test$adjustment_terms[[1L]]), "initial_ploidy+dose_mg")) {
    figure7_stop("Panel 7E/K adjusted CN-score plot received an incompatible analysis contract")
  }

  data$initial_ploidy <- factor(data$initial_ploidy, levels = c("2N", "4N"))
  data$dose <- factor(data$dose, levels = c("30mg/kg", "120mg/kg"))
  annotation <- sprintf(
    paste0(
      "Adjusted slope = %.2f TGI points / score SD\n",
      "Partial r = %.3f\n",
      "Exact origin x dose permutation P = %.3g\n",
      "n = %d"
    ),
    test$effect_per_within_origin_sd,
    test$partial_correlation,
    test$permutation_p_two_sided,
    test$n
  )
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = terminal_cn_score_nuisance_residual,
      y = tgi_origin_dose_residual,
      color = dose,
      shape = initial_ploidy
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::geom_smooth(
      data = data,
      ggplot2::aes(
        x = terminal_cn_score_nuisance_residual,
        y = tgi_origin_dose_residual,
        group = 1
      ),
      inherit.aes = FALSE,
      method = "lm", se = TRUE, color = "black", linewidth = 0.55
    ) +
    ggplot2::geom_point(size = 2.9) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = sample_id), size = 2.4, seed = 1,
      min.segment.length = 0, segment.color = "grey65",
      max.overlaps = Inf, show.legend = FALSE
    ) +
    ggplot2::annotate(
      "label", x = Inf, y = Inf, label = annotation,
      hjust = 1.05, vjust = 1.1, size = 2.8
    ) +
    ggplot2::scale_color_manual(values = figure7_dose_colors(), name = "Dose") +
    ggplot2::scale_shape_manual(
      values = c("2N" = 16, "4N" = 17),
      name = "Injected origin"
    ) +
    ggplot2::labs(
      title = paste(
        "Day", figure7_tgi_day(config),
        "TGI vs all-cell terminal postprocessed CN score"
      ),
      subtitle = paste0(
        "All 7,623 treated-tumor CBS cells; within-origin z score;\n",
        "association adjusted for injected origin and dose"
      ),
      x = paste0(
        "Terminal CN-score residual\n",
        "(within-origin z score; dose-adjusted)"
      ),
      y = paste0(
        "Day ", figure7_tgi_day(config), " TGI residual (%)\n",
        "(origin- and dose-adjusted)"
      ),
      shape = "Injected origin"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    figure7_theme() +
    ggplot2::theme(legend.position = "right")
}

figure7_build_ae <- function(cellcycle, data, samples, output_dir, config) {
  tgi_day <- figure7_tgi_day(config)
  tgi_measure <- figure7_tgi_measure(config)
  figures <- file.path(output_dir, "figures"); tables <- file.path(output_dir, "tables")
  trajectory <- figure7_growth_trajectory_data(cellcycle, samples, config)
  panel_b <- figure7_panel_b(data, samples)
  treated <- figure7_add_metadata(samples[samples$dose_mg > 0, , drop = FALSE], config)
  test_c <- figure7_add_metadata(figure7_panel_c_test(samples, config), config)
  shifts <- figure7_shift_metrics(data, samples, config); panel_d <- figure7_panel_d(shifts, config)
  panel_d$data <- figure7_add_metadata(panel_d$data, config); panel_d$test <- figure7_add_metadata(panel_d$test, config)
  panel_e <- figure7_panel_e(samples, config)
  panel_e$data <- figure7_add_metadata(panel_e$data, config); panel_e$test <- figure7_add_metadata(panel_e$test, config)
  panel_b$data <- figure7_add_metadata(panel_b$data, config); panel_b$tests <- figure7_add_metadata(panel_b$tests, config)

  figure7_write_tsv(trajectory, file.path(tables, "panel_7A_plot_data.tsv"))
  figure7_write_tsv(panel_b$data, file.path(tables, "panel_7B_plot_data.tsv"))
  figure7_write_tsv(panel_b$tests, file.path(tables, "panel_7B_tests.tsv"))
  figure7_write_tsv(treated, file.path(tables, "panel_7C_plot_data.tsv"))
  figure7_write_tsv(test_c, file.path(tables, "panel_7C_test.tsv"))
  figure7_write_tsv(panel_d$data, file.path(tables, "panel_7D_plot_data.tsv"))
  figure7_write_tsv(panel_d$test, file.path(tables, "panel_7D_test.tsv"))
  figure7_write_tsv(panel_e$data, file.path(tables, "panel_7E_plot_data.tsv"))
  figure7_write_tsv(panel_e$test, file.path(tables, "panel_7E_test.tsv"))

  filenames <- stats::setNames(figure7_panel_filenames(config), figure7_panel_ids(TRUE))
  plots <- list(
    A = figure7_panel_a_plot(trajectory, config),
    B = figure7_panel_b_plot(panel_b$data, panel_b$tests),
    C = figure7_panel_c_plot(treated, test_c, config),
    D = figure7_scatter_plot(
      panel_d$data,
      "shift_centered",
      "tgi_centered",
      panel_d$test,
      "CellCycle TGI association using injected-origin-matched controls",
      "Dose-centered ECDF RMSE (origin-matched untreated reference)",
      paste("Dose-centered Day", tgi_day, "TGI (%)"),
      "Exact within-dose permutation P"
    ) +
      ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35),
    E = figure7_adjusted_cn_plot(panel_e$data, panel_e$test, config)
  )
  sizes <- list(
    A = c(10, 6.5),
    B = c(15, 5.5),
    C = c(6.8, 6.4),
    D = c(6.6, 6.6),
    E = c(6.8, 6.8)
  )
  for (panel in names(plots)) {
    dimensions <- sizes[[panel]]
    figure7_save_panel(
      plots[[panel]],
      file.path(figures, filenames[[paste0("7", panel)]]),
      dimensions[[1L]],
      dimensions[[2L]]
    )
  }
  invisible(list(
    plots = plots,
    panel_b = panel_b,
    panel_c = test_c,
    panel_d = panel_d,
    panel_e = panel_e
  ))
}
