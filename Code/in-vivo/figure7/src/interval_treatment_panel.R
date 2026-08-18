# Frozen plotting-table contract and renderer for main Figure 7I.

figure7_interval_treatment_required_columns <- function() c(
  "panel_id", "more_positive_origin", "collection", "collection_label",
  "pathway", "pathway_label", "NES", "padj", "selected_direction",
  "selected_rank_within_direction", "selection_scope"
)

figure7_interval_treatment_collection_labels <- function() c(
  "H" = "Hallmark",
  "C2:CP:REACTOME" = "Reactome",
  "C5:GO:BP" = "GO BP"
)

figure7_interval_treatment_input_path <- function(repo_root, config) {
  path <- file.path(
    repo_root,
    as.character(config$interval_treatment_panel$selected_table)
  )
  if (!file.exists(path)) {
    figure7_stop("Missing frozen Figure 7I pathway table: ", path)
  }
  figure7_verify_checksum(
    path,
    as.character(config$interval_treatment_panel$selected_table_sha256),
    "frozen Figure 7I pathway table"
  )
  normalizePath(path, mustWork = TRUE)
}

figure7_validate_interval_treatment_table <- function(data, config) {
  required <- figure7_interval_treatment_required_columns()
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    figure7_stop(
      "Figure 7I pathway table is missing column(s): ",
      paste(missing, collapse = ", ")
    )
  }
  data <- data[, required, drop = FALSE]
  expected_per_direction <- as.integer(
    config$interval_treatment_panel$max_per_direction_across_collections
  )
  expected_panels <- c("negative_interaction", "positive_interaction")
  if (nrow(data) != 2L * expected_per_direction ||
      !setequal(unique(data$panel_id), expected_panels)) {
    figure7_stop(
      "Figure 7I pathway table must contain only the reviewed negative and positive interaction rows"
    )
  }
  required_nonempty <- c(
    "panel_id", "more_positive_origin", "collection", "pathway",
    "pathway_label", "selected_direction", "selection_scope"
  )
  if (any(vapply(
    data[required_nonempty],
    function(column) any(is.na(column) | !nzchar(as.character(column))),
    logical(1L)
  ))) {
    figure7_stop("Figure 7I pathway table contains a missing required value")
  }
  collection_labels <- figure7_interval_treatment_collection_labels()
  if (any(!data$collection %in% names(collection_labels))) {
    figure7_stop("Figure 7I pathway table contains an unknown collection")
  }
  data$NES <- figure7_numeric(data$NES)
  data$padj <- figure7_numeric(data$padj)
  data$selected_rank_within_direction <- suppressWarnings(as.integer(
    data$selected_rank_within_direction
  ))
  fdr_threshold <- as.numeric(config$interval_treatment_panel$fdr_threshold)
  expected_direction <- ifelse(data$NES < 0, "negative", "positive")
  expected_panel <- paste0(expected_direction, "_interaction")
  expected_origin <- ifelse(data$NES < 0, "2N", "4N")
  if (any(!is.finite(data$NES)) || any(data$NES == 0) ||
      any(!is.finite(data$padj)) || any(data$padj <= 0) ||
      any(data$padj > fdr_threshold) ||
      any(data$selected_direction != expected_direction) ||
      any(data$panel_id != expected_panel) ||
      any(data$more_positive_origin != expected_origin)) {
    figure7_stop(
      "Figure 7I pathway statistics violate the formal interaction FDR/direction contract"
    )
  }
  pathway_key <- paste(data$collection, data$pathway, sep = "\r")
  if (anyDuplicated(pathway_key)) {
    figure7_stop("Figure 7I interaction pathways are not unique")
  }
  for (direction in c("negative", "positive")) {
    local <- data[data$selected_direction == direction, , drop = FALSE]
    if (nrow(local) != expected_per_direction ||
        !identical(
          sort(local$selected_rank_within_direction),
          seq_len(expected_per_direction)
        )) {
      figure7_stop(
        "Figure 7I ", direction, " interaction selection ranks are invalid"
      )
    }
    observed_order <- local$selected_rank_within_direction[
      if (identical(direction, "negative")) {
        order(local$NES, local$padj, local$collection, local$pathway)
      } else {
        order(-local$NES, local$padj, local$collection, local$pathway)
      }
    ]
    if (!identical(observed_order, seq_len(expected_per_direction))) {
      figure7_stop(
        "Figure 7I ", direction,
        " pathways are not ranked by formal interaction NES"
      )
    }
  }
  data$collection_label <- unname(collection_labels[data$collection])
  data
}

figure7_read_interval_treatment_table <- function(path, config) {
  figure7_validate_interval_treatment_table(
    figure7_read_tsv(path, figure7_interval_treatment_required_columns()),
    config
  )
}

figure7_interval_treatment_plot <- function(data, config) {
  panel_labels <- c(
    negative_interaction = "More positive in 2N",
    positive_interaction = "More positive in 4N"
  )
  data$panel_label <- factor(
    unname(panel_labels[data$panel_id]),
    levels = unname(panel_labels)
  )
  data$collection_label <- factor(
    data$collection_label,
    levels = unname(figure7_interval_treatment_collection_labels())
  )
  data$minus_log10_fdr <- -log10(pmax(data$padj, .Machine$double.xmin))
  data$pathway_plot_key <- paste(
    data$panel_id, data$collection, data$pathway, sep = "\r"
  )
  ordered_pathway_keys <- unlist(lapply(names(panel_labels), function(panel_id) {
    local <- data[data$panel_id == panel_id, , drop = FALSE]
    local <- local[
      order(-local$selected_rank_within_direction),
      ,
      drop = FALSE
    ]
    local$pathway_plot_key
  }), use.names = FALSE)
  data$pathway_plot_key <- factor(
    data$pathway_plot_key,
    levels = ordered_pathway_keys
  )
  labels <- as.character(data$pathway_label[
    match(ordered_pathway_keys, data$pathway_plot_key)
  ])
  pathway_labels <- stats::setNames(
    vapply(
      labels,
      function(label) paste(strwrap(label, width = 43L), collapse = "\n"),
      character(1L)
    ),
    ordered_pathway_keys
  )
  symmetric_limit <- max(abs(data$NES), na.rm = TRUE) * 1.08
  interval <- config$interval_treatment_panel
  expected_per_direction <- as.integer(
    interval$max_per_direction_across_collections
  )

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = NES, y = pathway_plot_key,
      color = collection_label, size = minus_log10_fdr
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0, color = "#6B7280", linewidth = 0.4
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = NES, yend = pathway_plot_key),
      color = "#9AA4B2", linewidth = 0.7
    ) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(
      panel_label ~ ., scales = "free_y", space = "free_y"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(-symmetric_limit, symmetric_limit),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_y_discrete(labels = pathway_labels) +
    ggplot2::scale_color_manual(values = c(
      "Hallmark" = "#0072B2",
      "Reactome" = "#D55E00",
      "GO BP" = "#009E73"
    )) +
    ggplot2::scale_size_continuous(range = c(2.5, 5.6)) +
    ggplot2::labs(
      title = "Origin-dependent treatment responses",
      subtitle = sprintf(
        "Formal interaction; %.3f-%.3f pseudotime",
        as.numeric(interval$interval_start),
        as.numeric(interval$interval_end)
      ),
      x = "Treatment-by-origin interaction NES",
      y = NULL,
      color = "Collection",
      size = expression(-log[10]~"FDR"),
      caption = paste0(
        "Top ", expected_per_direction, " negative and top ",
        expected_per_direction,
        " positive formal treatment-by-origin interaction NES across collections. ",
        "All displayed pathways pass collection-wise BH FDR <= ",
        format(as.numeric(interval$fdr_threshold), trim = TRUE), "."
      )
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 8.5, color = "#374151"),
      plot.caption = ggplot2::element_text(
        size = 6.9, hjust = 0, color = "#4B5563"
      ),
      legend.position = "top",
      legend.justification = "left",
      strip.text.y = ggplot2::element_text(angle = 270, face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

figure7_build_interval_treatment_panel <- function(
  input_path,
  output_dir,
  config,
  verify_configured_checksum = TRUE,
  copy_table = TRUE
) {
  if (isTRUE(verify_configured_checksum)) {
    figure7_verify_checksum(
      input_path,
      as.character(config$interval_treatment_panel$selected_table_sha256),
      "frozen Figure 7I pathway table"
    )
  }
  data <- figure7_read_interval_treatment_table(input_path, config)
  table_output <- file.path(
    output_dir,
    "tables",
    "panel_7I_origin_comparison_selected.tsv"
  )
  if (isTRUE(copy_table)) {
    figure7_copy_file(input_path, table_output)
  }
  plot <- figure7_interval_treatment_plot(data, config)
  figure7_save_panel(
    plot,
    file.path(output_dir, "figures", config$panels$filenames[["7I"]]),
    10,
    10.5
  )
  invisible(list(
    plot = plot,
    data = data,
    input_path = normalizePath(input_path, mustWork = TRUE),
    input_sha256 = figure7_sha256(input_path),
    n_negative_interactions = sum(data$NES < 0),
    n_positive_interactions = sum(data$NES > 0)
  ))
}
