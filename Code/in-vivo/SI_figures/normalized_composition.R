# Shared, sample-aware composition plotting for Supplementary Figures 4 and 5.
#
# Cells are first converted to within-sample cluster proportions.  Samples are
# then averaged with equal weight inside each displayed group, so sequencing
# depth cannot determine a group's contribution.  Enrichment is tested on the
# same independent-sample proportions with exact label permutations inside
# configured nuisance strata and BH correction within each panel. Panels with
# one biological sample per displayed group remain explicitly descriptive.

composition_require_columns <- function(data, columns, label = "composition data") {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(
      label, " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

composition_star <- function(q_value) {
  output <- rep("", length(q_value))
  output[is.finite(q_value) & q_value <= 0.05] <- "*"
  output[is.finite(q_value) & q_value <= 0.01] <- "**"
  output[is.finite(q_value) & q_value <= 0.001] <- "***"
  output
}

prepare_equal_sample_composition <- function(
  data,
  unit_col,
  cluster_col,
  group_col,
  cluster_levels,
  group_levels,
  x_col = group_col,
  facet_cols = character(),
  strata_cols = character()
) {
  required <- unique(c(
    unit_col,
    cluster_col,
    group_col,
    x_col,
    facet_cols,
    strata_cols
  ))
  composition_require_columns(data, required)
  if (!nrow(data)) stop("Composition data contain no cells", call. = FALSE)
  if (anyNA(data[, required, drop = FALSE])) {
    stop("Composition grouping columns contain missing values", call. = FALSE)
  }
  if (any(!nzchar(as.character(data[[unit_col]]))) ||
      any(!nzchar(as.character(data[[cluster_col]]))) ||
      any(!nzchar(as.character(data[[group_col]])))) {
    stop("Composition grouping columns contain empty values", call. = FALSE)
  }

  observed_clusters <- unique(as.character(data[[cluster_col]]))
  unexpected_clusters <- setdiff(observed_clusters, cluster_levels)
  if (length(unexpected_clusters)) {
    stop(
      "Unexpected composition clusters: ",
      paste(unexpected_clusters, collapse = ", "),
      call. = FALSE
    )
  }
  observed_groups <- unique(as.character(data[[group_col]]))
  unexpected_groups <- setdiff(observed_groups, group_levels)
  if (length(unexpected_groups)) {
    stop(
      "Unexpected composition groups: ",
      paste(unexpected_groups, collapse = ", "),
      call. = FALSE
    )
  }
  missing_groups <- setdiff(group_levels, observed_groups)
  if (length(missing_groups)) {
    stop(
      "Composition groups have no independent samples: ",
      paste(missing_groups, collapse = ", "),
      call. = FALSE
    )
  }

  unit_values <- unique(as.character(data[[unit_col]]))
  unit_factor <- factor(as.character(data[[unit_col]]), levels = unit_values)
  cluster_factor <- factor(as.character(data[[cluster_col]]), levels = cluster_levels)
  counts <- as.data.frame(
    table(unit_value = unit_factor, cluster = cluster_factor),
    stringsAsFactors = FALSE
  )
  names(counts)[names(counts) == "Freq"] <- "n_cells"
  counts$unit_value <- as.character(counts$unit_value)
  counts$cluster <- as.character(counts$cluster)
  counts$n_cells <- as.numeric(counts$n_cells)
  counts$total_cells <- ave(counts$n_cells, counts$unit_value, FUN = sum)
  if (any(!is.finite(counts$total_cells)) || any(counts$total_cells <= 0)) {
    stop("Every independent sample must contain at least one cell", call. = FALSE)
  }
  counts$sample_proportion <- counts$n_cells / counts$total_cells

  metadata <- data.frame(
    unit_value = as.character(data[[unit_col]]),
    group_value = as.character(data[[group_col]]),
    x_value = as.character(data[[x_col]]),
    stringsAsFactors = FALSE
  )
  for (column in unique(c(facet_cols, strata_cols))) {
    metadata[[column]] <- data[[column]]
  }
  metadata <- unique(metadata)
  unit_multiplicity <- table(metadata$unit_value)
  if (any(unit_multiplicity != 1L)) {
    offending <- names(unit_multiplicity)[unit_multiplicity != 1L]
    stop(
      "Each independent sample must map to one display/test group; offending: ",
      paste(offending, collapse = ", "),
      call. = FALSE
    )
  }
  counts <- merge(
    counts,
    metadata,
    by = "unit_value",
    all.x = TRUE,
    sort = FALSE
  )
  if (anyNA(counts$group_value) || anyNA(counts$x_value)) {
    stop("Independent-sample metadata failed to join", call. = FALSE)
  }
  display_columns <- c("group_value", "x_value", facet_cols)
  display_mapping <- unique(counts[, display_columns, drop = FALSE])
  group_multiplicity <- table(as.character(display_mapping$group_value))
  if (any(group_multiplicity != 1L)) {
    offending <- names(group_multiplicity)[group_multiplicity != 1L]
    stop(
      "Each composition group must map to one x/facet location; offending: ",
      paste(offending, collapse = ", "),
      call. = FALSE
    )
  }
  counts$unit_value <- factor(counts$unit_value, levels = unit_values)
  counts$cluster <- factor(counts$cluster, levels = cluster_levels)
  counts$group_value <- factor(counts$group_value, levels = group_levels)

  sample_sums <- stats::aggregate(
    sample_proportion ~ unit_value,
    data = counts,
    FUN = sum
  )
  if (any(abs(sample_sums$sample_proportion - 1) > 1e-12)) {
    stop("Within-sample cluster proportions do not sum to one", call. = FALSE)
  }
  counts
}

composition_sample_arrays <- function(
  sample_composition,
  cluster_levels,
  strata_cols = character()
) {
  composition_require_columns(
    sample_composition,
    c(
      "unit_value",
      "cluster",
      "group_value",
      "sample_proportion",
      strata_cols
    )
  )
  units <- levels(droplevels(sample_composition$unit_value))
  if (is.null(units)) units <- unique(as.character(sample_composition$unit_value))
  n_units <- length(units)
  if (n_units < 2L) {
    stop("Composition enrichment requires at least two independent samples", call. = FALSE)
  }

  proportion_matrix <- matrix(
    NA_real_,
    nrow = n_units,
    ncol = length(cluster_levels),
    dimnames = list(units, cluster_levels)
  )
  unit_index <- match(as.character(sample_composition$unit_value), units)
  cluster_index <- match(as.character(sample_composition$cluster), cluster_levels)
  proportion_matrix[cbind(unit_index, cluster_index)] <-
    as.numeric(sample_composition$sample_proportion)
  if (any(!is.finite(proportion_matrix))) {
    stop("Sample composition matrix is incomplete or nonfinite", call. = FALSE)
  }

  unit_rows <- !duplicated(as.character(sample_composition$unit_value))
  unit_metadata <- sample_composition[
    unit_rows,
    c("unit_value", "group_value", strata_cols),
    drop = FALSE
  ]
  unit_metadata <- unit_metadata[
    match(units, as.character(unit_metadata$unit_value)),
    ,
    drop = FALSE
  ]
  unit_groups <- as.character(unit_metadata$group_value)
  if (anyNA(unit_groups)) stop("Sample groups are incomplete", call. = FALSE)
  stratum <- if (length(strata_cols)) {
    do.call(
      interaction,
      c(
        unit_metadata[, strata_cols, drop = FALSE],
        list(drop = TRUE, lex.order = TRUE, sep = " | ")
      )
    )
  } else {
    factor(rep("all samples", n_units))
  }
  if (anyNA(stratum)) stop("Sample permutation strata are incomplete", call. = FALSE)

  list(
    units = units,
    n_units = n_units,
    proportion_matrix = proportion_matrix,
    unit_groups = unit_groups,
    stratum = stratum,
    strata_label = if (length(strata_cols)) {
      paste(strata_cols, collapse = "+")
    } else {
      "none"
    }
  )
}

stratified_selection_matrix <- function(is_group, stratum) {
  stratum_indices <- split(seq_along(is_group), stratum, drop = TRUE)
  selections <- lapply(stratum_indices, function(indices) {
    n_selected <- sum(is_group[indices])
    if (n_selected == 0L) return(list(integer()))
    if (n_selected == length(indices)) return(list(indices))
    utils::combn(indices, n_selected, simplify = FALSE)
  })
  n_permutations <- prod(lengths(selections))
  if (!is.finite(n_permutations) || n_permutations > 1e6) {
    stop(
      "Exact composition permutation space exceeds 1,000,000 assignments",
      call. = FALSE
    )
  }
  if (n_permutations < 2L) {
    stop(
      paste(
        "Composition enrichment is not testable within the requested strata;",
        "use descriptive mode"
      ),
      call. = FALSE
    )
  }
  selection_grid <- do.call(
    expand.grid,
    c(
      lapply(selections, seq_along),
      list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    )
  )
  output <- matrix(
    FALSE,
    nrow = length(is_group),
    ncol = nrow(selection_grid)
  )
  for (permutation_index in seq_len(nrow(selection_grid))) {
    selected <- unlist(lapply(seq_along(selections), function(stratum_index) {
      selections[[stratum_index]][[
        selection_grid[permutation_index, stratum_index]
      ]]
    }))
    output[selected, permutation_index] <- TRUE
  }
  output
}

exact_sample_enrichment <- function(
  sample_composition,
  cluster_levels,
  group_levels,
  strata_cols = character(),
  fdr_threshold = 0.05
) {
  arrays <- composition_sample_arrays(
    sample_composition,
    cluster_levels,
    strata_cols
  )

  results <- vector("list", length(group_levels) * length(cluster_levels))
  result_index <- 1L
  for (group in group_levels) {
    is_group <- arrays$unit_groups == group
    active_strata <- names(which(vapply(
      split(is_group, arrays$stratum, drop = TRUE),
      function(values) any(values) && any(!values),
      logical(1)
    )))
    active <- as.character(arrays$stratum) %in% active_strata
    n_group <- sum(is_group & active)
    n_rest <- sum(!is_group & active)
    if (n_group < 1L || n_rest < 1L) {
      stop(
        paste(
          "Every enrichment group requires at least one sample and comparator",
          "within an exchangeable stratum"
        ),
        call. = FALSE
      )
    }
    selection_matrix <- stratified_selection_matrix(is_group, arrays$stratum)
    for (cluster in cluster_levels) {
      values <- arrays$proportion_matrix[active, cluster]
      active_group <- is_group[active]
      observed <- mean(values[active_group]) - mean(values[!active_group])
      selected_sums <- as.numeric(crossprod(selection_matrix[active, , drop = FALSE], values))
      null_difference <-
        selected_sums / n_group -
        (sum(values) - selected_sums) / n_rest
      p_value <- mean(null_difference >= observed - 1e-14)
      results[[result_index]] <- data.frame(
        group_value = group,
        cluster = cluster,
        n_group_samples = n_group,
        n_rest_samples = n_rest,
        mean_group_proportion = mean(values[active_group]),
        mean_rest_proportion = mean(values[!active_group]),
        difference = observed,
        exact_permutations = length(null_difference),
        p_value = p_value,
        stringsAsFactors = FALSE
      )
      result_index <- result_index + 1L
    }
  }
  output <- do.call(rbind, results)
  output$q_value <- stats::p.adjust(output$p_value, method = "BH")
  output$testable <- TRUE
  output$enriched <- output$difference > 0 & output$q_value <= fdr_threshold
  output$significance <- ifelse(
    output$enriched,
    composition_star(output$q_value),
    ""
  )
  output$test <- paste(
    "exact equal-sample one-group-versus-rest label permutation; strata:",
    arrays$strata_label
  )
  output$strata <- arrays$strata_label
  output$adjustment <- "BH across all group-by-cluster contrasts in panel"
  output
}

descriptive_sample_enrichment <- function(
  sample_composition,
  cluster_levels,
  group_levels,
  reason
) {
  arrays <- composition_sample_arrays(sample_composition, cluster_levels)
  results <- vector("list", length(group_levels) * length(cluster_levels))
  result_index <- 1L
  for (group in group_levels) {
    is_group <- arrays$unit_groups == group
    n_group <- sum(is_group)
    n_rest <- arrays$n_units - n_group
    if (n_group < 1L || n_rest < 1L) {
      stop("Every composition group requires at least one comparator", call. = FALSE)
    }
    for (cluster in cluster_levels) {
      values <- arrays$proportion_matrix[, cluster]
      results[[result_index]] <- data.frame(
        group_value = group,
        cluster = cluster,
        n_group_samples = n_group,
        n_rest_samples = n_rest,
        mean_group_proportion = mean(values[is_group]),
        mean_rest_proportion = mean(values[!is_group]),
        difference = mean(values[is_group]) - mean(values[!is_group]),
        exact_permutations = NA_real_,
        p_value = NA_real_,
        stringsAsFactors = FALSE
      )
      result_index <- result_index + 1L
    }
  }
  output <- do.call(rbind, results)
  output$q_value <- NA_real_
  output$testable <- FALSE
  output$enriched <- FALSE
  output$significance <- ""
  output$test <- paste("descriptive only:", reason)
  output$strata <- "not applicable"
  output$adjustment <- "not applicable"
  output
}

aggregate_equal_sample_composition <- function(
  sample_composition,
  cluster_levels,
  group_levels,
  facet_cols = character()
) {
  split_key <- interaction(
    sample_composition$group_value,
    sample_composition$cluster,
    drop = TRUE,
    lex.order = TRUE
  )
  pieces <- split(sample_composition, split_key)
  output <- do.call(rbind, lapply(pieces, function(piece) {
    row <- data.frame(
      group_value = as.character(piece$group_value[[1L]]),
      cluster = as.character(piece$cluster[[1L]]),
      n_samples = length(unique(piece$unit_value)),
      sum_n_cells = sum(piece$n_cells),
      mean_proportion = mean(piece$sample_proportion),
      min_sample_proportion = min(piece$sample_proportion),
      max_sample_proportion = max(piece$sample_proportion),
      x_value = as.character(piece$x_value[[1L]]),
      stringsAsFactors = FALSE
    )
    for (facet in facet_cols) row[[facet]] <- as.character(piece[[facet]][[1L]])
    row
  }))
  rownames(output) <- NULL
  output$group_value <- factor(output$group_value, levels = group_levels)
  output$cluster <- factor(output$cluster, levels = cluster_levels)
  for (facet in facet_cols) {
    source_values <- sample_composition[[facet]]
    output[[facet]] <- if (is.factor(source_values)) {
      factor(
        output[[facet]],
        levels = levels(source_values),
        ordered = is.ordered(source_values)
      )
    } else {
      factor(output[[facet]], levels = unique(as.character(source_values)))
    }
  }
  output <- output[order(output$group_value, output$cluster), , drop = FALSE]
  group_sums <- stats::aggregate(
    mean_proportion ~ group_value,
    data = output,
    FUN = sum
  )
  if (any(abs(group_sums$mean_proportion - 1) > 1e-12)) {
    stop("Equal-sample mean cluster proportions do not sum to one", call. = FALSE)
  }
  output
}

add_composition_facets <- function(
  plot,
  facet_formula,
  facet_type,
  facet_scales,
  facet_space,
  facet_nrow
) {
  if (is.null(facet_formula)) return(plot)
  if (identical(facet_type, "grid")) {
    return(plot + ggplot2::facet_grid(
      facet_formula,
      scales = facet_scales,
      space = facet_space,
      drop = FALSE
    ))
  }
  plot + ggplot2::facet_wrap(
    facet_formula,
    scales = facet_scales,
    nrow = facet_nrow,
    drop = FALSE
  )
}

make_normalized_composition_plot <- function(
  data,
  unit_col,
  cluster_col,
  group_col,
  cluster_levels,
  group_levels,
  fill_colors,
  bar_axis = c("cluster", "group"),
  x_col = group_col,
  x_levels = NULL,
  facet_cols = character(),
  strata_cols = character(),
  facet_formula = NULL,
  facet_type = c("wrap", "grid"),
  facet_scales = "fixed",
  facet_space = "fixed",
  facet_nrow = NULL,
  title,
  subtitle = NULL,
  x_title,
  y_title,
  legend_title,
  tag,
  theme_function,
  x_text_angle = 35,
  show_n = FALSE,
  test_mode = c("permutation", "descriptive"),
  descriptive_reason = NULL,
  fdr_threshold = 0.05
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Normalized composition plots require ggplot2", call. = FALSE)
  }
  bar_axis <- match.arg(bar_axis)
  facet_type <- match.arg(facet_type)
  test_mode <- match.arg(test_mode)
  sample_data <- prepare_equal_sample_composition(
    data = data,
    unit_col = unit_col,
    cluster_col = cluster_col,
    group_col = group_col,
    cluster_levels = cluster_levels,
    group_levels = group_levels,
    x_col = x_col,
    facet_cols = facet_cols,
    strata_cols = strata_cols
  )
  plot_data <- aggregate_equal_sample_composition(
    sample_data,
    cluster_levels,
    group_levels,
    facet_cols
  )
  tests <- if (test_mode == "permutation") {
    exact_sample_enrichment(
      sample_composition = sample_data,
      cluster_levels = cluster_levels,
      group_levels = group_levels,
      strata_cols = strata_cols,
      fdr_threshold = fdr_threshold
    )
  } else {
    if (is.null(descriptive_reason) || !nzchar(descriptive_reason)) {
      stop("Descriptive composition mode requires a reason", call. = FALSE)
    }
    descriptive_sample_enrichment(
      sample_composition = sample_data,
      cluster_levels = cluster_levels,
      group_levels = group_levels,
      reason = descriptive_reason
    )
  }
  match_key <- paste(plot_data$group_value, plot_data$cluster, sep = "\r")
  test_key <- paste(tests$group_value, tests$cluster, sep = "\r")
  test_match <- match(match_key, test_key)
  if (anyNA(test_match)) {
    stop("Composition tests failed to reconcile with plot data", call. = FALSE)
  }
  plot_data$significance <- tests$significance[test_match]
  plot_data$q_value <- tests$q_value[test_match]
  plot_data$difference <- tests$difference[test_match]
  plot_data$testable <- tests$testable[test_match]
  plot_data$x_value <- factor(
    plot_data$x_value,
    levels = if (is.null(x_levels)) unique(plot_data$x_value) else x_levels
  )

  if (bar_axis == "cluster") {
    if (!all(group_levels %in% names(fill_colors))) {
      stop("Fill colors do not cover every composition group", call. = FALSE)
    }
    plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = cluster,
        y = mean_proportion,
        fill = group_value,
        group = group_value
      )
    ) +
      ggplot2::geom_col(
        position = ggplot2::position_dodge(width = 0.84),
        width = 0.78,
        color = "white",
        linewidth = 0.12
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = significance),
        position = ggplot2::position_dodge(width = 0.84),
        vjust = -0.28,
        size = 3.1,
        color = "black",
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = fill_colors[group_levels],
        drop = FALSE,
        name = legend_title
      )
  } else {
    if (!all(cluster_levels %in% names(fill_colors))) {
      stop("Fill colors do not cover every composition cluster", call. = FALSE)
    }
    plot_data$stack_midpoint <- ave(
      plot_data$mean_proportion,
      plot_data$group_value,
      FUN = function(value) 1 - cumsum(value) + value / 2
    )
    plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = x_value,
        y = mean_proportion,
        fill = cluster,
        group = cluster
      )
    ) +
      ggplot2::geom_col(
        width = 0.82,
        color = "white",
        linewidth = 0.10
      ) +
      ggplot2::geom_text(
        data = plot_data[nzchar(plot_data$significance), , drop = FALSE],
        ggplot2::aes(y = stack_midpoint, label = significance),
        size = 2.7,
        color = "black",
        fontface = "bold",
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = fill_colors[cluster_levels],
        drop = FALSE,
        name = legend_title
      )
    if (show_n) {
      n_rows <- plot_data[!duplicated(plot_data$group_value), , drop = FALSE]
      plot <- plot + ggplot2::geom_text(
        data = n_rows,
        ggplot2::aes(
          x = x_value,
          y = 0.995,
          label = paste0("n=", n_samples)
        ),
        inherit.aes = FALSE,
        size = 2.8,
        vjust = 1
      )
    }
  }

  plot <- plot +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      labels = function(values) paste0(round(100 * values), "%"),
      expand = ggplot2::expansion(mult = c(0, if (bar_axis == "cluster") 0.14 else 0))
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_title,
      y = y_title
    ) +
    theme_function(9.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = x_text_angle,
        hjust = if (x_text_angle == 0) 0.5 else 1
      )
    )
  plot <- if (bar_axis == "cluster") {
    plot + ggplot2::coord_cartesian(clip = "off")
  } else {
    plot + ggplot2::coord_cartesian(ylim = c(0, 1), clip = "off")
  }
  plot <- add_composition_facets(
    plot,
    facet_formula,
    facet_type,
    facet_scales,
    facet_space,
    facet_nrow
  )
  plot <- plot + ggplot2::labs(tag = tag)

  tests$panel_tag <- tag
  plot_data$panel_tag <- tag
  sample_data$panel_tag <- tag
  list(
    plot = plot,
    plot_data = plot_data,
    tests = tests,
    sample_data = sample_data
  )
}
