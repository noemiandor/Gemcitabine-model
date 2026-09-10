# Frozen plotting-table contract and renderer for main Figure 7I.

figure7_interval_treatment_required_columns <- function() c(
  "panel_id", "more_positive_origin", "collection", "collection_label",
  "pathway", "pathway_label", "NES", "padj", "selected_direction",
  "NES_2N", "padj_2N", "NES_4N", "padj_4N",
  "origin_effect_nes_order_concordant", "interaction_fdr_significant",
  "selected_rank_within_direction", "selection_scope"
)

figure7_interval_treatment_collection_labels <- function() c(
  "H" = "Hallmark",
  "C2:CP:REACTOME" = "Reactome",
  "C5:GO:BP" = "Gene Ontology BP"
)

figure7_interval_treatment_reader_labels <- function() c(
  GOBP_HOMOPHILIC_CELL_CELL_ADHESION =
    "Homophilic cell-cell adhesion",
  HALLMARK_INTERFERON_ALPHA_RESPONSE =
    "Interferon-alpha response",
  GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I =
    "MHC class I peptide antigen presentation",
  GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_PEPTIDE_ANTIGEN =
    "Endogenous peptide antigen presentation",
  REACTOME_INTERFERON_ALPHA_BETA_SIGNALING =
    "Interferon-alpha/beta signaling",
  GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_ANTIGEN =
    "Endogenous antigen presentation",
  GOBP_EMBRYONIC_PLACENTA_MORPHOGENESIS =
    "Embryonic placenta morphogenesis",
  GOBP_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY =
    "NK cell-mediated immunity regulation",
  GOBP_CYCLIC_NUCLEOTIDE_METABOLIC_PROCESS =
    "Cyclic nucleotide metabolism",
  GOBP_NEGATIVE_REGULATION_OF_VIRAL_GENOME_REPLICATION =
    "Negative regulation of viral replication",
  HALLMARK_MYC_TARGETS_V1 =
    "MYC targets V1",
  GOBP_RIBOSOME_BIOGENESIS =
    "Ribosome biogenesis",
  GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS =
    "Ribonucleoprotein complex biogenesis",
  GOBP_RRNA_PROCESSING =
    "rRNA processing",
  GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS =
    "Small ribosomal subunit biogenesis",
  REACTOME_RRNA_PROCESSING =
    "rRNA processing",
  GOBP_RRNA_METABOLIC_PROCESS =
    "rRNA metabolism",
  GOBP_MATURATION_OF_SSU_RRNA =
    "SSU rRNA maturation",
  REACTOME_EUKARYOTIC_TRANSLATION_INITIATION =
    "Eukaryotic translation initiation",
  GOBP_MITOCHONDRIAL_TRANSLATION =
    "Mitochondrial translation"
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
  data$NES_2N <- figure7_numeric(data$NES_2N)
  data$padj_2N <- figure7_numeric(data$padj_2N)
  data$NES_4N <- figure7_numeric(data$NES_4N)
  data$padj_4N <- figure7_numeric(data$padj_4N)
  data$selected_rank_within_direction <- suppressWarnings(as.integer(
    data$selected_rank_within_direction
  ))
  parse_logical <- function(values, label) {
    normalized <- toupper(as.character(values))
    if (any(!normalized %in% c("TRUE", "FALSE"))) {
      figure7_stop("Figure 7I pathway table has invalid ", label, " values")
    }
    normalized == "TRUE"
  }
  data$origin_effect_nes_order_concordant <- parse_logical(
    data$origin_effect_nes_order_concordant,
    "origin-effect concordance"
  )
  data$interaction_fdr_significant <- parse_logical(
    data$interaction_fdr_significant,
    "interaction-significance"
  )
  fdr_threshold <- as.numeric(config$interval_treatment_panel$fdr_threshold)
  expected_direction <- ifelse(data$NES < 0, "negative", "positive")
  expected_panel <- paste0(expected_direction, "_interaction")
  expected_origin <- ifelse(data$NES < 0, "2N", "4N")
  expected_origin_order <- ifelse(
    data$NES < 0,
    data$NES_2N > data$NES_4N,
    data$NES_4N > data$NES_2N
  )
  if (any(!is.finite(data$NES)) || any(data$NES == 0) ||
      any(!is.finite(data$padj)) || any(data$padj <= 0) ||
      any(data$padj > fdr_threshold) ||
      any(!is.finite(data$NES_2N)) || any(!is.finite(data$NES_4N)) ||
      any(!is.finite(data$padj_2N)) || any(data$padj_2N <= 0) ||
      any(data$padj_2N > 1) ||
      any(!is.finite(data$padj_4N)) || any(data$padj_4N <= 0) ||
      any(data$padj_4N > 1) ||
      any(data$selected_direction != expected_direction) ||
      any(data$panel_id != expected_panel) ||
      any(data$more_positive_origin != expected_origin) ||
      any(!data$interaction_fdr_significant) ||
      any(!data$origin_effect_nes_order_concordant) ||
      any(!expected_origin_order)) {
    figure7_stop(
      "Figure 7I pathway statistics violate the formal interaction FDR/direction/origin-effect contract"
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
    negative_interaction = "2N-more-positive\ninteraction",
    positive_interaction = "4N-more-positive\ninteraction"
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
  reader_labels <- figure7_interval_treatment_reader_labels()
  if (!setequal(names(reader_labels), as.character(data$pathway))) {
    figure7_stop(
      "Figure 7I reader labels do not match the frozen pathway selection"
    )
  }
  labels <- unname(reader_labels[as.character(data$pathway)[
    match(ordered_pathway_keys, data$pathway_plot_key)
  ]])
  pathway_labels <- stats::setNames(
    vapply(
      labels,
      function(label) paste(strwrap(label, width = 55L), collapse = "\n"),
      character(1L)
    ),
    ordered_pathway_keys
  )

  plot_data <- do.call(rbind, lapply(c("2N", "4N"), function(origin) {
    data.frame(
      panel_label = as.character(data$panel_label),
      collection_label = as.character(data$collection_label),
      pathway_plot_key = as.character(data$pathway_plot_key),
      origin = origin,
      origin_NES = data[[paste0("NES_", origin)]],
      minus_log10_fdr = data$minus_log10_fdr,
      stringsAsFactors = FALSE
    )
  }))
  plot_data$panel_label <- factor(
    plot_data$panel_label,
    levels = unname(panel_labels)
  )
  plot_data$collection_label <- factor(
    plot_data$collection_label,
    levels = unname(figure7_interval_treatment_collection_labels())
  )
  plot_data$pathway_plot_key <- factor(
    plot_data$pathway_plot_key,
    levels = ordered_pathway_keys
  )
  plot_data$origin <- factor(plot_data$origin, levels = c("2N", "4N"))

  symmetric_limit <- max(
    abs(c(data$NES_2N, data$NES_4N)),
    na.rm = TRUE
  ) * 1.08
  interval <- config$interval_treatment_panel
  expected_per_direction <- as.integer(
    interval$max_per_direction_across_collections
  )

  ggplot2::ggplot() +
    ggplot2::geom_vline(
      xintercept = 0, color = "#6B7280", linewidth = 0.4
    ) +
    ggplot2::geom_segment(
      data = data,
      ggplot2::aes(
        x = NES_2N, xend = NES_4N,
        y = pathway_plot_key, yend = pathway_plot_key
      ),
      color = "#7C8798", linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(
        x = origin_NES, y = pathway_plot_key,
        color = collection_label,
        shape = origin,
        size = minus_log10_fdr
      )
    ) +
    ggplot2::facet_grid(
      panel_label ~ ., scales = "free_y", space = "free_y"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(-symmetric_limit, symmetric_limit),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_y_discrete(
      labels = pathway_labels,
      expand = ggplot2::expansion(add = c(0.85, 1.10))
    ) +
    ggplot2::scale_color_manual(values = c(
      "Hallmark" = "#0072B2",
      "Reactome" = "#D55E00",
      "Gene Ontology BP" = "#009E73"
    )) +
    ggplot2::scale_shape_manual(
      values = c("2N" = 16, "4N" = 17),
      breaks = names(figure7_sum159_origin_labels()),
      labels = unname(figure7_sum159_origin_labels())
    ) +
    ggplot2::scale_size_continuous(range = c(2.5, 5.6)) +
    ggplot2::labs(
      title = "Origin-dependent treatment responses",
      subtitle = sprintf(
        paste0(
          "Equal-dose treatment contrast; %.3f-%.3f pseudotime\n",
          "2N: %d cells across %d mice; 4N: %d cells across %d mice"
        ),
        as.numeric(interval$interval_start),
        as.numeric(interval$interval_end),
        as.integer(interval$expected_cells_2N),
        as.integer(interval$expected_mice_per_origin),
        as.integer(interval$expected_cells_4N),
        as.integer(interval$expected_mice_per_origin)
      ),
      x = "Normalized enrichment score (treated - vehicle)",
      y = NULL,
      color = "Collection",
      shape = "Origin",
      size = expression(-log[10]~"FDR"),
      caption = paste0(
        "Pathways were selected as the top ", expected_per_direction,
        " negative and top ", expected_per_direction,
          " positive formal interaction NES across Hallmark, Reactome, and Gene Ontology BP.\n",
        "Points show joint-model origin-specific treatment NES; connectors join the 2N and 4N estimates for each pathway.\n",
        "The upper panel has a more positive response in 2N and the lower panel in 4N; more positive can also mean less depleted.\n",
        "Point size represents formal-interaction FDR; all displayed pathways pass collection-wise BH FDR <= ",
        format(as.numeric(interval$fdr_threshold), trim = TRUE), ". ",
        "GRCh38 human-tumor mouse pseudobulks; adjusted for mean within-interval pseudotime; exploratory data-selected interval."
      )
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1),
      shape = ggplot2::guide_legend(order = 2),
      size = ggplot2::guide_legend(order = 3)
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
