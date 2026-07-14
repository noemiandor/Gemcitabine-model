# Visual explanations of the tumor-growth-inhibition (TGI) calculation.
#
# This module is intentionally separate from 04h_pseudotime_TGI_essential_util.R.
# Source the essential utility first so that `tgi_spec` can be created with
# essential_tgi_spec(), then source this file. The public entry point is
# essential_tgi_calculation_visualizations().

essential_tgi_calc_numeric <- function(x) suppressWarnings(as.numeric(x))

essential_tgi_calc_clean <- function(x) {
  out <- trimws(as.character(x))
  out[out %in% c("", "NA", "NaN", "NULL", "None")] <- NA_character_
  out
}

essential_tgi_calc_unique <- function(data, column, sample_id, numeric = FALSE) {
  values <- data[data$sample_id == sample_id, column]
  if (numeric) {
    values <- essential_tgi_calc_numeric(values)
    values <- unique(values[is.finite(values)])
  } else {
    values <- essential_tgi_calc_clean(values)
    values <- unique(values[!is.na(values)])
  }
  if (length(values) == 0L) return(if (numeric) NA_real_ else NA_character_)
  if (length(values) > 1L) {
    stop(
      "Sample ", sample_id, " has inconsistent values in ", column, ": ",
      paste(values, collapse = ", "),
      call. = FALSE
    )
  }
  values[[1L]]
}

essential_tgi_calc_validate_spec <- function(tgi_spec) {
  required <- c(
    "outcome", "control_summary", "day", "measure", "delta_column",
    "short_label", "title_label", "axis_label", "reference_label"
  )
  missing <- setdiff(required, names(tgi_spec))
  if (length(missing) > 0L) {
    stop("tgi_spec is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!tgi_spec$outcome %in% c("auc", "day")) {
    stop("Unsupported TGI outcome: ", tgi_spec$outcome, call. = FALSE)
  }
  if (!tgi_spec$control_summary %in% c("mean", "median", "max")) {
    stop("Unsupported matched-control summary: ", tgi_spec$control_summary, call. = FALSE)
  }
  if (identical(tgi_spec$outcome, "day") && !is.finite(tgi_spec$day)) {
    stop("A finite selected day is required for endpoint TGI", call. = FALSE)
  }
  invisible(TRUE)
}

essential_tgi_calc_summary <- function(x, method) {
  x <- essential_tgi_calc_numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  switch(
    method,
    mean = mean(x),
    median = stats::median(x),
    max = max(x),
    stop("Unsupported matched-control summary: ", method, call. = FALSE)
  )
}

essential_tgi_calc_day_columns <- function(data) {
  columns <- grep("^tumor_volume_Day_[0-9]+$", names(data), value = TRUE)
  days <- essential_tgi_calc_numeric(sub("^tumor_volume_Day_", "", columns))
  columns[order(days)]
}

essential_tgi_calc_dose_label <- function(dose_mg) {
  ifelse(is.finite(dose_mg), paste0(format(dose_mg, trim = TRUE), "mg/kg"), NA_character_)
}

essential_tgi_calc_sample_table <- function(data, tgi_spec) {
  essential_tgi_calc_validate_spec(tgi_spec)
  day_columns <- essential_tgi_calc_day_columns(data)
  required <- c(
    "sample_id", "initial_ploidy", "gemcitabine_dose_mg_per_kg",
    "tumor_volume_baseline_day", "tumor_volume_baseline",
    tgi_spec$delta_column, day_columns
  )
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      "TGI calculation plots require input column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  sample_ids <- sort(unique(essential_tgi_calc_clean(data$sample_id)))
  sample_ids <- sample_ids[!is.na(sample_ids)]
  if (length(sample_ids) == 0L) stop("No samples found in TGI input", call. = FALSE)

  extract_row <- function(sample_id) {
    row <- data.frame(
      sample_id = sample_id,
      initial_ploidy = essential_tgi_calc_unique(data, "initial_ploidy", sample_id),
      dose_mg = essential_tgi_calc_unique(
        data, "gemcitabine_dose_mg_per_kg", sample_id, numeric = TRUE
      ),
      baseline_day = essential_tgi_calc_unique(
        data, "tumor_volume_baseline_day", sample_id
      ),
      baseline_volume = essential_tgi_calc_unique(
        data, "tumor_volume_baseline", sample_id, numeric = TRUE
      ),
      selected_delta = essential_tgi_calc_unique(
        data, tgi_spec$delta_column, sample_id, numeric = TRUE
      ),
      stringsAsFactors = FALSE
    )
    for (column in day_columns) {
      row[[column]] <- essential_tgi_calc_unique(data, column, sample_id, numeric = TRUE)
    }
    if (tgi_spec$measure %in% names(data)) {
      row$input_tgi_percent <- essential_tgi_calc_unique(
        data, tgi_spec$measure, sample_id, numeric = TRUE
      )
    } else {
      row$input_tgi_percent <- NA_real_
    }
    row
  }
  samples <- do.call(rbind, lapply(sample_ids, extract_row))
  samples$dose <- essential_tgi_calc_dose_label(samples$dose_mg)
  samples$treatment_status <- ifelse(samples$dose_mg > 0, "treated", "untreated control")
  samples <- samples[order(samples$initial_ploidy, samples$dose_mg, samples$sample_id), , drop = FALSE]
  rownames(samples) <- NULL
  samples
}

essential_tgi_calculation_components <- function(data, tgi_spec, control_dose_mg = 0) {
  samples <- essential_tgi_calc_sample_table(data, tgi_spec)
  controls <- samples[
    samples$dose_mg == control_dose_mg & is.finite(samples$selected_delta),
    ,
    drop = FALSE
  ]
  treated <- samples[
    samples$dose_mg > control_dose_mg & is.finite(samples$selected_delta),
    ,
    drop = FALSE
  ]
  if (nrow(treated) == 0L) stop("No treated mice have a finite selected TGI delta", call. = FALSE)

  rows <- lapply(seq_len(nrow(treated)), function(index) {
    mouse <- treated[index, , drop = FALSE]
    matched <- controls[controls$initial_ploidy == mouse$initial_ploidy, , drop = FALSE]
    if (nrow(matched) == 0L) {
      stop(
        "No finite untreated-control delta for initial ploidy ",
        mouse$initial_ploidy,
        call. = FALSE
      )
    }
    reference_delta <- essential_tgi_calc_summary(
      matched$selected_delta,
      tgi_spec$control_summary
    )
    tgi_percent <- 100 * (1 - mouse$selected_delta / reference_delta)
    endpoint_volume <- if (identical(tgi_spec$outcome, "day")) {
      column <- paste0("tumor_volume_Day_", tgi_spec$day)
      if (!column %in% names(mouse)) NA_real_ else mouse[[column]][[1L]]
    } else {
      NA_real_
    }
    data.frame(
      sample_id = mouse$sample_id,
      initial_ploidy = mouse$initial_ploidy,
      dose_mg = mouse$dose_mg,
      dose = mouse$dose,
      baseline_day = mouse$baseline_day,
      baseline_volume = mouse$baseline_volume,
      selected_day = if (identical(tgi_spec$outcome, "day")) tgi_spec$day else NA_real_,
      endpoint_volume = endpoint_volume,
      treated_delta = mouse$selected_delta,
      matched_control_reference_delta = reference_delta,
      matched_control_summary = tgi_spec$control_summary,
      matched_control_n = nrow(matched),
      matched_control_sample_ids = paste(matched$sample_id, collapse = ";"),
      matched_control_deltas = paste(
        format(matched$selected_delta, digits = 10, trim = TRUE),
        collapse = ";"
      ),
      tgi_percent = if (is.finite(tgi_percent)) tgi_percent else NA_real_,
      input_tgi_percent = mouse$input_tgi_percent,
      recalculation_difference = if (
        is.finite(tgi_percent) && is.finite(mouse$input_tgi_percent)
      ) {
        tgi_percent - mouse$input_tgi_percent
      } else {
        NA_real_
      },
      tgi_measure = tgi_spec$measure,
      tgi_outcome = tgi_spec$outcome,
      tgi_control_summary = tgi_spec$control_summary,
      delta_column = tgi_spec$delta_column,
      delta_units = if (identical(tgi_spec$outcome, "auc")) {
        "mm^3 * day"
      } else {
        "mm^3"
      },
      formula = paste0(
        "100 * (1 - ", format(mouse$selected_delta, digits = 7, trim = TRUE),
        " / ", format(reference_delta, digits = 7, trim = TRUE), ")"
      ),
      stringsAsFactors = FALSE
    )
  })
  components <- do.call(rbind, rows)
  components <- components[
    order(components$initial_ploidy, components$dose_mg, components$sample_id),
    ,
    drop = FALSE
  ]
  rownames(components) <- NULL
  components
}

essential_tgi_trajectory_plot_data <- function(
  data,
  tgi_spec,
  highlight_day = NULL,
  control_dose_mg = 0
) {
  samples <- essential_tgi_calc_sample_table(data, tgi_spec)
  day_columns <- essential_tgi_calc_day_columns(samples)
  observed_days <- essential_tgi_calc_numeric(sub("^tumor_volume_Day_", "", day_columns))
  baseline_days <- essential_tgi_calc_numeric(sub("^Day_", "", samples$baseline_day))
  if (is.null(highlight_day)) {
    highlight_day <- if (identical(tgi_spec$outcome, "day")) tgi_spec$day else 17
  }
  highlight_day <- essential_tgi_calc_numeric(highlight_day)[[1L]]
  if (!is.finite(highlight_day)) stop("highlight_day must be finite", call. = FALSE)
  if (identical(tgi_spec$outcome, "day") && highlight_day != tgi_spec$day) {
    stop(
      "For endpoint TGI, highlight_day must equal the selected TGI day (Day ",
      tgi_spec$day, ")",
      call. = FALSE
    )
  }
  all_days <- sort(unique(c(baseline_days[is.finite(baseline_days)], observed_days)))
  if (!highlight_day %in% all_days) {
    stop("Day ", highlight_day, " is not present in the tumor-volume input", call. = FALSE)
  }

  rows <- lapply(seq_len(nrow(samples)), function(index) {
    mouse <- samples[index, , drop = FALSE]
    baseline_day <- essential_tgi_calc_numeric(sub("^Day_", "", mouse$baseline_day))
    day <- c(baseline_day, observed_days)
    volume <- c(
      mouse$baseline_volume,
      essential_tgi_calc_numeric(mouse[1L, day_columns, drop = TRUE])
    )
    keep <- is.finite(day) & is.finite(volume) & is.finite(mouse$baseline_volume)
    data.frame(
      sample_id = mouse$sample_id,
      initial_ploidy = mouse$initial_ploidy,
      dose_mg = mouse$dose_mg,
      dose = mouse$dose,
      treatment_status = mouse$treatment_status,
      day = day[keep],
      tumor_volume = volume[keep],
      tumor_volume_change = volume[keep] - mouse$baseline_volume,
      series_type = "individual mouse",
      reference_n = NA_integer_,
      control_summary = tgi_spec$control_summary,
      highlight_day = highlight_day,
      is_highlight_day = day[keep] == highlight_day,
      tgi_measure = tgi_spec$measure,
      stringsAsFactors = FALSE
    )
  })
  individual <- do.call(rbind, rows)
  individual <- individual[!duplicated(individual[c("sample_id", "day")]), , drop = FALSE]

  control <- individual[
    individual$dose_mg == control_dose_mg & is.finite(individual$tumor_volume_change),
    ,
    drop = FALSE
  ]
  reference_rows <- lapply(
    split(control, interaction(control$initial_ploidy, control$day, drop = TRUE)),
    function(local) {
      data.frame(
        sample_id = paste0(local$initial_ploidy[[1L]], " matched-control reference"),
        initial_ploidy = local$initial_ploidy[[1L]],
        dose_mg = control_dose_mg,
        dose = "matched-control reference",
        treatment_status = "untreated control reference",
        day = local$day[[1L]],
        tumor_volume = NA_real_,
        tumor_volume_change = essential_tgi_calc_summary(
          local$tumor_volume_change,
          tgi_spec$control_summary
        ),
        series_type = "matched-control reference",
        reference_n = sum(is.finite(local$tumor_volume_change)),
        control_summary = tgi_spec$control_summary,
        highlight_day = highlight_day,
        is_highlight_day = local$day[[1L]] == highlight_day,
        tgi_measure = tgi_spec$measure,
        stringsAsFactors = FALSE
      )
    }
  )
  reference <- do.call(rbind, reference_rows)
  output <- rbind(individual, reference)
  output <- output[order(output$initial_ploidy, output$series_type, output$sample_id, output$day), , drop = FALSE]
  rownames(output) <- NULL
  output
}

essential_tgi_calculation_plot_data <- function(components) {
  required <- c(
    "sample_id", "initial_ploidy", "dose_mg", "dose", "treated_delta",
    "matched_control_reference_delta", "tgi_percent", "tgi_measure",
    "tgi_outcome", "tgi_control_summary", "delta_units"
  )
  missing <- setdiff(required, names(components))
  if (length(missing) > 0L) {
    stop("TGI components table is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  treated <- components
  treated$component <- "Treated-mouse growth delta"
  treated$component_value <- treated$treated_delta
  reference <- components
  reference$component <- "Matched-control reference delta"
  reference$component_value <- reference$matched_control_reference_delta
  output <- rbind(treated, reference)
  output$component <- factor(
    output$component,
    levels = c("Treated-mouse growth delta", "Matched-control reference delta")
  )
  output$tgi_label <- ifelse(
    is.finite(output$tgi_percent),
    paste0("TGI = ", formatC(output$tgi_percent, format = "f", digits = 1), "%"),
    "TGI unavailable"
  )
  output
}

essential_tgi_growth_trajectory_plot <- function(trajectory_data, tgi_spec) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required", call. = FALSE)
  individual <- trajectory_data[trajectory_data$series_type == "individual mouse", , drop = FALSE]
  reference <- trajectory_data[
    trajectory_data$series_type == "matched-control reference",
    ,
    drop = FALSE
  ]
  highlighted <- individual[individual$is_highlight_day, , drop = FALSE]
  highlighted_reference <- reference[reference$is_highlight_day, , drop = FALSE]
  highlight_day <- unique(trajectory_data$highlight_day)
  if (length(highlight_day) != 1L) stop("Trajectory data has inconsistent highlight days", call. = FALSE)
  dose_levels <- unique(individual$dose[order(individual$dose_mg)])
  individual$dose <- factor(individual$dose, levels = dose_levels)
  highlighted$dose <- factor(highlighted$dose, levels = dose_levels)
  default_colors <- c("0mg/kg" = "#777777", "30mg/kg" = "#D95F02", "120mg/kg" = "#1B9E77")
  missing_colors <- setdiff(dose_levels, names(default_colors))
  if (length(missing_colors) > 0L) {
    extra <- grDevices::hcl.colors(length(missing_colors), "Dark 3")
    names(extra) <- missing_colors
    default_colors <- c(default_colors, extra)
  }
  colors <- default_colors[dose_levels]
  summary_word <- tgi_spec$control_summary
  subtitle <- if (identical(tgi_spec$outcome, "day")) {
    paste0(
      "Day ", tgi_spec$day,
      " is the selected endpoint: TGI = 100 x (1 - treated delta / ",
      summary_word, " matched-control delta)."
    )
  } else {
    paste0(
      "AUC is integrated for each mouse before applying the matched-control ",
      summary_word, "; Day ", highlight_day,
      " is highlighted for context. The dashed curve is a pointwise control summary."
    )
  }
  ggplot2::ggplot() +
    ggplot2::geom_vline(
      xintercept = highlight_day,
      color = "#B2182B",
      linewidth = 0.9,
      alpha = 0.8
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::geom_line(
      data = individual,
      ggplot2::aes(
        x = day,
        y = tumor_volume_change,
        group = sample_id,
        color = dose
      ),
      linewidth = 0.55,
      alpha = 0.72
    ) +
    ggplot2::geom_line(
      data = reference,
      ggplot2::aes(x = day, y = tumor_volume_change, group = initial_ploidy),
      color = "black",
      linetype = "22",
      linewidth = 1.15
    ) +
    ggplot2::geom_point(
      data = highlighted,
      ggplot2::aes(x = day, y = tumor_volume_change, color = dose),
      shape = 21,
      fill = "white",
      stroke = 1,
      size = 2.5
    ) +
    ggplot2::geom_point(
      data = highlighted_reference,
      ggplot2::aes(x = day, y = tumor_volume_change),
      shape = 23,
      fill = "#B2182B",
      color = "black",
      stroke = 0.5,
      size = 3.1
    ) +
    ggplot2::facet_wrap(~initial_ploidy, nrow = 1L) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE, name = "Dose") +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(individual$day)),
      minor_breaks = NULL
    ) +
    ggplot2::labs(
      title = paste0("How ", tgi_spec$title_label, " is calculated from tumor growth"),
      subtitle = subtitle,
      x = "Study day",
      y = expression(Delta * " tumor volume from Day 0 (mm"^3 * ")"),
      caption = paste0(
        "Thin lines: individual mice. Dashed black: ", tgi_spec$reference_label,
        ". Red line and outlined points: Day ", highlight_day, "."
      )
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey70"),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.caption = ggplot2::element_text(hjust = 0)
    )
}

essential_tgi_component_plot <- function(calculation_plot_data, tgi_spec) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required", call. = FALSE)
  components <- calculation_plot_data[
    calculation_plot_data$component == "Treated-mouse growth delta",
    ,
    drop = FALSE
  ]
  sample_levels <- components$sample_id[
    order(components$initial_ploidy, components$dose_mg, components$sample_id)
  ]
  calculation_plot_data$sample_id <- factor(calculation_plot_data$sample_id, levels = sample_levels)
  components$sample_id <- factor(components$sample_id, levels = sample_levels)
  delta_label <- if (identical(tgi_spec$outcome, "auc")) {
    expression("AUC of " * Delta * " tumor volume (mm"^3 * " x day)")
  } else {
    paste0("Tumor-volume change, Day 0 to Day ", tgi_spec$day, " (mm^3)")
  }
  ggplot2::ggplot(
    calculation_plot_data,
    ggplot2::aes(x = sample_id, y = component_value)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::geom_segment(
      data = components,
      ggplot2::aes(
        x = sample_id,
        xend = sample_id,
        y = treated_delta,
        yend = matched_control_reference_delta
      ),
      inherit.aes = FALSE,
      color = "grey55",
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = component),
      shape = 21,
      color = "black",
      stroke = 0.4,
      size = 3.2
    ) +
    ggplot2::geom_text(
      data = components,
      ggplot2::aes(
        x = sample_id,
        y = treated_delta,
        label = tgi_label
      ),
      inherit.aes = FALSE,
      vjust = -0.8,
      size = 2.8
    ) +
    ggplot2::facet_grid(
      . ~ initial_ploidy,
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Treated-mouse growth delta" = "#377EB8",
        "Matched-control reference delta" = "#F4A261"
      ),
      name = NULL
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = paste0("Per-treated-mouse components of ", tgi_spec$title_label),
      subtitle = paste0(
        "Each vertical pair shows the numerator and denominator in TGI = 100 x (1 - treated delta / reference delta).\n",
        "Reference = ", tgi_spec$reference_label, "."
      ),
      x = "Treated mouse",
      y = delta_label,
      caption = "TGI labels are recalculated from the two plotted values; untreated mice contribute only to the matched-control reference."
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey70"),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.margin = ggplot2::margin(5.5, 22, 5.5, 5.5),
      plot.caption = ggplot2::element_text(hjust = 0)
    )
}

essential_tgi_calculation_visualizations <- function(
  data,
  tgi_spec,
  highlight_day = NULL,
  control_dose_mg = 0
) {
  components <- essential_tgi_calculation_components(
    data,
    tgi_spec,
    control_dose_mg = control_dose_mg
  )
  trajectory_data <- essential_tgi_trajectory_plot_data(
    data,
    tgi_spec,
    highlight_day = highlight_day,
    control_dose_mg = control_dose_mg
  )
  calculation_data <- essential_tgi_calculation_plot_data(components)
  metadata <- data.frame(
    tgi_measure = tgi_spec$measure,
    tgi_outcome = tgi_spec$outcome,
    tgi_day = if (identical(tgi_spec$outcome, "day")) tgi_spec$day else NA_real_,
    highlighted_day = unique(trajectory_data$highlight_day),
    matched_control_summary = tgi_spec$control_summary,
    matched_control_scope = "untreated mice matched by initial ploidy",
    calculation_targets = "treated mice only",
    formula = "100 * (1 - treated growth delta / matched-control reference growth delta)",
    stringsAsFactors = FALSE
  )
  list(
    trajectory_plot = essential_tgi_growth_trajectory_plot(trajectory_data, tgi_spec),
    calculation_plot = essential_tgi_component_plot(calculation_data, tgi_spec),
    trajectory_plot_data = trajectory_data,
    calculation_plot_data = calculation_data,
    calculation_components = components,
    metadata = metadata
  )
}

essential_tgi_calculation_figure_filenames <- function() {
  c(
    "TGI_calculation_growth_trajectories.pdf",
    "TGI_calculation_components_by_treated_mouse.pdf"
  )
}

essential_tgi_calculation_plot_data_filenames <- function() {
  c(
    "TGI_calculation_growth_trajectories_plot_data.csv",
    "TGI_calculation_components_by_treated_mouse_plot_data.csv"
  )
}

essential_tgi_calculation_stats_filenames <- function() {
  c(
    "TGI_calculation_components_by_treated_mouse.csv",
    "TGI_calculation_metadata.csv"
  )
}

essential_write_tgi_calculation_outputs <- function(
  cellcycle_path,
  noncellcycle_path,
  output_root,
  tgi_spec
) {
  cellcycle <- essential_read_input(cellcycle_path, "CellCycle")
  noncellcycle <- essential_read_input(noncellcycle_path, "NonCellCycle")
  selected <- essential_apply_tgi_spec(cellcycle, noncellcycle, tgi_spec)
  result <- essential_tgi_calculation_visualizations(
    selected$cellcycle,
    tgi_spec,
    highlight_day = if (identical(tgi_spec$outcome, "day")) tgi_spec$day else 17
  )

  figures_dir <- essential_ensure_dir(file.path(output_root, "Figures", "TGI_calculation"))
  stats_dir <- essential_ensure_dir(file.path(output_root, "stats", "TGI_calculation"))
  plot_data_dir <- essential_ensure_dir(file.path(output_root, "plot_data", "TGI_calculation"))

  essential_save_pdf(
    result$trajectory_plot,
    file.path(figures_dir, "TGI_calculation_growth_trajectories.pdf"),
    10,
    6.5
  )
  essential_save_pdf(
    result$calculation_plot,
    file.path(figures_dir, "TGI_calculation_components_by_treated_mouse.pdf"),
    10,
    6.5
  )
  essential_write_csv(
    result$trajectory_plot_data,
    file.path(plot_data_dir, "TGI_calculation_growth_trajectories_plot_data.csv")
  )
  essential_write_csv(
    result$calculation_plot_data,
    file.path(plot_data_dir, "TGI_calculation_components_by_treated_mouse_plot_data.csv")
  )
  essential_write_csv(
    result$calculation_components,
    file.path(stats_dir, "TGI_calculation_components_by_treated_mouse.csv")
  )
  essential_write_csv(
    result$metadata,
    file.path(stats_dir, "TGI_calculation_metadata.csv")
  )

  invisible(list(figures = 2L, stats = 2L, plot_data = 2L))
}

essential_regenerate_tgi_calculation_figures <- function(tables_root, output_root) {
  stats_dir <- file.path(tables_root, "stats", "TGI_calculation")
  plot_data_dir <- file.path(tables_root, "plot_data", "TGI_calculation")
  metadata <- essential_read_figure_table(
    file.path(stats_dir, "TGI_calculation_metadata.csv"),
    c("tgi_measure", "highlighted_day"),
    "TGI calculation metadata"
  )
  if (nrow(metadata) != 1L) {
    stop("TGI calculation metadata must contain exactly one row", call. = FALSE)
  }
  tgi_spec <- essential_tgi_spec_from_measure(metadata$tgi_measure[[1L]])
  trajectory_data <- essential_read_figure_table(
    file.path(plot_data_dir, "TGI_calculation_growth_trajectories_plot_data.csv"),
    c(
      "sample_id", "initial_ploidy", "dose_mg", "dose", "day",
      "tumor_volume_change", "series_type", "highlight_day", "is_highlight_day",
      "tgi_measure"
    ),
    "TGI calculation trajectory plot data"
  )
  calculation_data <- essential_read_figure_table(
    file.path(plot_data_dir, "TGI_calculation_components_by_treated_mouse_plot_data.csv"),
    c(
      "sample_id", "initial_ploidy", "dose_mg", "treated_delta",
      "matched_control_reference_delta", "component", "component_value",
      "tgi_percent", "tgi_label", "tgi_measure"
    ),
    "TGI calculation component plot data"
  )
  if (any(as.character(trajectory_data$tgi_measure) != tgi_spec$measure) ||
      any(as.character(calculation_data$tgi_measure) != tgi_spec$measure)) {
    stop("TGI calculation plot data do not match saved metadata", call. = FALSE)
  }
  trajectory_data$is_highlight_day <- as.character(trajectory_data$is_highlight_day) %in%
    c("TRUE", "T", "1")

  figures_dir <- essential_ensure_dir(file.path(output_root, "Figures", "TGI_calculation"))
  essential_save_pdf(
    essential_tgi_growth_trajectory_plot(trajectory_data, tgi_spec),
    file.path(figures_dir, "TGI_calculation_growth_trajectories.pdf"),
    10,
    6.5
  )
  essential_save_pdf(
    essential_tgi_component_plot(calculation_data, tgi_spec),
    file.path(figures_dir, "TGI_calculation_components_by_treated_mouse.pdf"),
    10,
    6.5
  )
  invisible(list(figures = 2L))
}
