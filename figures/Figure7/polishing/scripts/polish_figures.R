#!/usr/bin/env Rscript

# Publication-scale Figure 7 rebuild. Scientific summaries are recomputed from
# the tracked plotting inputs; the final composite is assembled from live plot
# and grob objects, never from exported panel rasters.

parse_args <- function(args) {
  output <- list(phase = "all")
  index <- 1L
  while (index <= length(args)) {
    argument <- args[[index]]
    if (startsWith(argument, "--phase=")) {
      output$phase <- sub("^--phase=", "", argument)
    } else if (identical(argument, "--phase") && index < length(args)) {
      index <- index + 1L
      output$phase <- args[[index]]
    } else {
      stop("Unknown polishing argument: ", argument, call. = FALSE)
    }
    index <- index + 1L
  }
  if (!output$phase %in% c("subpanels", "final", "all")) {
    stop("--phase must be subpanels, final, or all", call. = FALSE)
  }
  output
}

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[[1L]])
} else {
  "figures/Figure7/polishing/scripts/polish_figures.R"
}
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."), mustWork = TRUE)
polish_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
args <- parse_args(commandArgs(trailingOnly = TRUE))

module_dir <- file.path(repo_root, "Code", "in-vivo", "figure7")
for (file in c(
  "common_io.R", "feature_species_policy.R", "tgi_data.R",
  "tgi_statistics.R", "tgi_panels.R", "context_panels.R",
  "copy_number_panel.R", "state_pathway_panel.R",
  "generated_state_pathway_reference.R",
  "state_pathway_analysis.R"
)) {
  sys.source(file.path(module_dir, "src", file), envir = .GlobalEnv)
}
sys.source(
  file.path(repo_root, "Code", "in-vivo", "SI_figures", "normalized_composition.R"),
  envir = .GlobalEnv
)
sys.source(
  file.path(repo_root, "Code", "in-vivo", "SI_figures", "shared_context_panels.R"),
  envir = .GlobalEnv
)
sys.source(
  file.path(repo_root, "Code", "in-vivo", "SI_figures", "copy_number_heatmap.R"),
  envir = .GlobalEnv
)

required_packages <- c("ggplot2", "ggrepel", "patchwork", "pheatmap", "yaml")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1L), quietly = TRUE
)]
if (length(missing_packages)) {
  stop(
    "Missing Figure 7 polishing package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

config_path <- file.path(module_dir, "figure7_config.yaml")
cellcycle_path <- file.path(
  repo_root, "Data", "in-vivo", "figure7", "processed",
  "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
)
noncellcycle_path <- file.path(
  repo_root, "Data", "in-vivo", "figure7", "processed",
  "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
)
endpoint_path <- file.path(
  repo_root, "Data", "in-vivo", "scRNAseq_Numbat", "all_ploidy.csv"
)
si_cache_dir <- file.path(repo_root, "Data", "in-vivo", "SIfigures")

build_panel_objects <- function() {
  config <- figure7_read_config(config_path)
  config <- figure7_attach_density_localization_config(
    config,
    file.path(module_dir, "density_localization_config.yaml")
  )
  endpoint <- figure7_read_endpoint_ploidy_table(endpoint_path, config)
  cellcycle <- figure7_read_cell_table(cellcycle_path, "CellCycle", config)
  noncellcycle <- figure7_read_cell_table(noncellcycle_path, "NonCellCycle", config)
  samples <- figure7_sample_table(cellcycle, noncellcycle, config, endpoint)
  cellcycle_analysis <- figure7_prepare_cellcycle(cellcycle, samples, config)
  trajectory <- figure7_growth_trajectory_data(cellcycle, samples, config)
  ecdf <- figure7_panel_b(cellcycle_analysis, samples, config)
  treated <- figure7_add_metadata(
    samples[samples$dose_mg > 0, , drop = FALSE], config
  )
  tgi_test <- figure7_add_metadata(figure7_panel_c_test(samples, config), config)
  shifts <- figure7_shift_metrics(cellcycle_analysis, samples, config)
  shift_association <- figure7_panel_d(shifts, config)
  shift_association$data <- figure7_add_metadata(shift_association$data, config)
  shift_association$test <- figure7_add_metadata(shift_association$test, config)
  endpoint_association <- figure7_panel_e(samples, config)
  endpoint_association$data <- figure7_add_metadata(endpoint_association$data, config)
  endpoint_association$test <- figure7_add_metadata(endpoint_association$test, config)

  source_plots <- list(
    A = figure7_panel_a_plot(trajectory, config),
    B = figure7_panel_b_plot(
      ecdf$data,
      ecdf$tests,
      ecdf$localization_grid,
      ecdf$localization_intervals,
      ecdf$localization_test
    ),
    C = figure7_panel_c_plot(treated, tgi_test, config),
    D = figure7_scatter_plot(
      shift_association$data,
      "shift_centered",
      "tgi_centered",
      shift_association$test,
      "CellCycle TGI association using injected-origin-matched controls",
      "Dose-centered ECDF RMSE (origin-matched untreated reference)",
      paste("Dose-centered Day", figure7_tgi_day(config), "TGI (%)"),
      "Exact within-dose permutation P",
      annotation_corner = "top-left"
    ) + ggplot2::geom_vline(
      xintercept = 0, color = "grey75", linewidth = 0.35
    ),
    E = figure7_endpoint_ploidy_plot(
      endpoint_association$data, endpoint_association$test, config
    )
  )
  context <- figure7_build_context_panels(
    si_cache_dir, repo_root, config, policy = "reviewed"
  )
  reference_root <- file.path(
    repo_root,
    as.character(config$state_pathways$reviewed_reference_root),
    as.character(config$state_pathways$reviewed_reference_id)
  )
  reference <- figure7_validate_reviewed_state_reference(reference_root, config)
  state_plot <- figure7_panel_f_plot(reference$activity, config)
  copy_number <- figure7_build_copy_number_panel(repo_root)
  mapped <- figure7_main_composite_plots(
    source_plots, context$plots, state_plot, copy_number$plot, config
  )
  list(config = config, plots = mapped, copy_number = copy_number)
}

slot_dimensions <- function() {
  spec <- figure7_publication_spec()
  content_width_in <- spec$width_in * (
    1 - spec$content_left_npc - spec$content_right_npc
  )
  widths <- c(
    A = 18 / 28, B = 10 / 28,
    C = 9 / 28, D = 10 / 28, E = 9 / 28,
    F = 10 / 28, G = 18 / 28,
    H = 1, I = 11 / 28, J = 17 / 28, K = 0.5, L = 0.5
  ) * content_width_in
  heights <- c(
    A = spec$row_heights[[1L]], B = spec$row_heights[[1L]],
    C = spec$row_heights[[2L]], D = spec$row_heights[[2L]],
    E = spec$row_heights[[2L]], F = spec$row_heights[[3L]],
    G = spec$row_heights[[3L]], H = spec$row_heights[[4L]],
    I = spec$row_heights[[5L]], J = spec$row_heights[[5L]],
    K = spec$row_heights[[6L]], L = spec$row_heights[[6L]]
  )
  data.frame(
    figure = "Figure 7",
    panel = tolower(names(widths)),
    width_in = as.numeric(widths),
    height_in = as.numeric(heights),
    stringsAsFactors = FALSE
  )
}

write_subpanels <- function(panel_objects) {
  subpanel_dir <- file.path(polish_root, "subpanels")
  layout_dir <- file.path(polish_root, "layout")
  dir.create(subpanel_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(layout_dir, recursive = TRUE, showWarnings = FALSE)
  dimensions <- slot_dimensions()
  styled <- figure7_publication_clean_plots(
    panel_objects$plots, panel_objects$config
  )
  spec <- figure7_publication_spec()
  paths <- character(nrow(dimensions))
  for (index in seq_len(nrow(dimensions))) {
    panel <- toupper(dimensions$panel[[index]])
    path <- file.path(subpanel_dir, paste0("Figure7_", panel, ".png"))
    audit_plot <- styled[[panel]]
    add_audit_tag <- function(plot) {
      plot +
        ggplot2::labs(tag = panel) +
        ggplot2::theme(
          plot.tag = ggplot2::element_text(
            family = "sans", face = "bold", size = spec$panel_tag_pt
          ),
          plot.tag.position = c(0, 1),
          plot.tag.location = "plot"
        )
    }
    panel_h_components <- attr(
      audit_plot, "figure7_panel_b_components", exact = TRUE
    )
    if (identical(panel, "H") && !is.null(panel_h_components)) {
      panel_h_components$ecdf <- add_audit_tag(panel_h_components$ecdf)
      audit_plot <- patchwork::wrap_plots(
        panel_h_components,
        ncol = 1,
        heights = c(1.35, 0.945),
        guides = "collect"
      )
    } else if (inherits(audit_plot, "ggplot")) {
      audit_plot <- add_audit_tag(audit_plot)
    }
    ggplot2::ggsave(
      path,
      plot = audit_plot,
      width = dimensions$width_in[[index]],
      height = dimensions$height_in[[index]],
      units = "in", dpi = spec$png_dpi, bg = "white", limitsize = FALSE
    )
    paths[[index]] <- file.path(
      "figures", "Figure7", "polishing", "subpanels", basename(path)
    )
  }
  dimensions$subpanel_png <- paths
  # ggsave/ragg truncate non-integer device extents to whole pixels.
  dimensions$width_px <- floor(dimensions$width_in * spec$png_dpi)
  dimensions$height_px <- floor(dimensions$height_in * spec$png_dpi)
  dimensions$dpi <- spec$png_dpi
  dimensions <- dimensions[, c(
    "figure", "panel", "subpanel_png", "width_px", "height_px",
    "width_in", "height_in", "dpi"
  )]
  utils::write.csv(
    dimensions,
    file.path(layout_dir, "subpanel_dimensions.csv"),
    row.names = FALSE,
    quote = TRUE
  )
  adopted_plan_path <- file.path(layout_dir, "layout_plan.csv")
  if (file.exists(adopted_plan_path)) {
    adopted <- utils::read.csv(adopted_plan_path, stringsAsFactors = FALSE)
    adopted$xmax <- adopted$x_in + adopted$width_in
    adopted$ymax <- adopted$y_in + adopted$height_in
    preview <- ggplot2::ggplot(adopted) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = x_in, xmax = xmax, ymin = y_in, ymax = ymax,
          fill = panel
        ),
        color = "#333333", linewidth = 0.35, show.legend = FALSE
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          x = x_in + width_in / 2,
          y = y_in + height_in / 2,
          label = toupper(panel)
        ),
        family = "sans", fontface = "bold", size = 5
      ) +
      ggplot2::coord_fixed(
        xlim = c(0, spec$width_in), ylim = c(0, spec$height_in),
        expand = FALSE, clip = "off"
      ) +
      ggplot2::scale_fill_manual(
        values = stats::setNames(
          rep(c("#DCEAF7", "#FBE5D6", "#E2F0D9", "#FFF2CC"), 3L)[seq_len(nrow(adopted))],
          adopted$panel
        )
      ) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_void(base_family = "sans") +
      ggplot2::theme(
        plot.margin = ggplot2::margin(8, 8, 8, 8, unit = "pt")
      )
    ggplot2::ggsave(
      file.path(layout_dir, "Figure_7_adopted_layout_preview.png"),
      plot = preview, width = spec$width_in, height = spec$height_in,
      units = "in", dpi = 150, bg = "white", limitsize = FALSE
    )
  }
  invisible(dimensions)
}

relative_repo_path <- function(path) {
  normalized <- normalizePath(path, mustWork = FALSE)
  prefix <- paste0(repo_root, .Platform$file.sep)
  if (!startsWith(normalized, prefix)) {
    stop("Polishing output escaped repository root: ", normalized, call. = FALSE)
  }
  substring(normalized, nchar(prefix) + 1L)
}

write_final_records <- function(output_paths) {
  relative_paths <- vapply(output_paths, relative_repo_path, character(1L))
  hashes <- vapply(output_paths, figure7_sha256, character(1L))
  sizes <- as.numeric(file.info(output_paths)$size)
  canonical_paths <- file.path(
    repo_root, "figures", "Figure7", basename(output_paths)
  )
  canonical_relative_paths <- vapply(
    canonical_paths, relative_repo_path, character(1L)
  )
  canonical_hashes <- vapply(canonical_paths, function(path) {
    if (file.exists(path)) figure7_sha256(path) else NA_character_
  }, character(1L))
  hash_status <- ifelse(
    is.na(canonical_hashes),
    "canonical_not_available",
    ifelse(hashes == canonical_hashes, "match", "mismatch")
  )
  manifest <- data.frame(
    path = relative_paths,
    size_bytes = sizes,
    sha256 = hashes,
    role = c("publication_png", "publication_vector_pdf"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    manifest, file.path(polish_root, "manifest.csv"), row.names = FALSE
  )

  direct_input_paths <- c(
    relative_repo_path(cellcycle_path),
    relative_repo_path(noncellcycle_path),
    relative_repo_path(endpoint_path),
    relative_repo_path(file.path(si_cache_dir, "manifest.tsv")),
    relative_repo_path(file.path(
      si_cache_dir, "si_figures_cell_metadata.csv"
    )),
    relative_repo_path(file.path(
      si_cache_dir, "si_figure6_endpoint_ploidy_join_audit.csv"
    )),
    "Data/in-vivo/all_ploidy.tsv",
    "Data/in-vivo/scRNAseq_Numbat/cbs_manifest.tsv",
    relative_repo_path(file.path(
      repo_root,
      "Data/in-vivo/figure7/saved_state_pathway/",
      "state_pathway_grch_human_only_initial_ploidy_day17_v3/",
      "panel_7F_pathway_activity_plot_data.tsv"
    ))
  )
  cbs_manifest <- utils::read.delim(
    file.path(repo_root, "Data/in-vivo/scRNAseq_Numbat/cbs_manifest.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  direct_input_paths <- unique(c(
    direct_input_paths,
    file.path("Data/in-vivo/scRNAseq_Numbat", cbs_manifest$filename)
  ))
  direct_input_hashes <- vapply(
    file.path(repo_root, direct_input_paths), figure7_sha256, character(1L)
  )
  dependency_paths <- c(
    relative_repo_path(config_path),
    "Code/in-vivo/figure7/density_localization_config.yaml",
    "Code/in-vivo/figure7/src/common_io.R",
    "Code/in-vivo/figure7/src/feature_species_policy.R",
    "Code/in-vivo/figure7/src/tgi_data.R",
    "Code/in-vivo/figure7/src/tgi_statistics.R",
    "Code/in-vivo/figure7/src/tgi_panels.R",
    "Code/in-vivo/figure7/src/context_panels.R",
    "Code/in-vivo/figure7/src/copy_number_panel.R",
    "Code/in-vivo/figure7/src/state_pathway_panel.R",
    "Code/in-vivo/figure7/src/generated_state_pathway_reference.R",
    "Code/in-vivo/figure7/src/state_pathway_analysis.R",
    "Code/in-vivo/SI_figures/shared_context_panels.R",
    "Code/in-vivo/SI_figures/normalized_composition.R",
    "Code/in-vivo/SI_figures/copy_number_heatmap.R",
    "figures/Figure7/polishing/layout/layout_plan.csv",
    "figures/Figure7/polishing/scripts/polish_figures.R",
    "scripts/agentRrunner.sh"
  )
  dependency_hashes <- vapply(
    file.path(repo_root, dependency_paths), figure7_sha256, character(1L)
  )
  command <- paste(
    "scripts/agentRrunner.sh",
    "figures/Figure7/polishing/scripts/polish_figures.R --phase all"
  )
  rebuild <- data.frame(
    stage = "polishing",
    manuscript_figure_name = "Figure 7",
    current_figure_name = "Figure7_reviewed_GRCh",
    source_root = "figures/Figure7/polishing",
    script_path = "figures/Figure7/polishing/scripts/polish_figures.R",
    working_directory = ".",
    command = command,
    target_image_path = canonical_relative_paths,
    target_sha256 = canonical_hashes,
    rebuild_output_path = relative_paths,
    rebuild_sha256 = hashes,
    hash_match_status = hash_status,
    direct_input_paths = paste(direct_input_paths, collapse = ";"),
    direct_input_sha256 = paste(direct_input_hashes, collapse = ";"),
    dependency_paths = paste(dependency_paths, collapse = ";"),
    dependency_sha256 = paste(dependency_hashes, collapse = ";"),
    immutable_raster_inputs = "NA",
    stringsAsFactors = FALSE
  )
  figure7_write_tsv(
    rebuild, file.path(polish_root, "figure_rebuild_manifest.tsv")
  )
  byte_report <- data.frame(
    output = relative_paths,
    sha256 = hashes,
    comparison_target = canonical_relative_paths,
    comparison_sha256 = canonical_hashes,
    status = hash_status,
    notes = c(
      paste(
        "PNG byte comparison against the independently generated canonical",
        "Manager output"
      ),
      paste(
        "PDF byte comparison is informational: cairo embeds CreationDate",
        "metadata, so independent vector exports are not byte-stable; page",
        "geometry and rendering are checked in visual_qc.md"
      )
    ),
    stringsAsFactors = FALSE
  )
  figure7_write_tsv(
    byte_report, file.path(polish_root, "figure_byte_identity_report.tsv")
  )
}

write_final <- function(panel_objects) {
  final_dir <- file.path(polish_root, "final_images")
  dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
  spec <- figure7_publication_spec()
  composite <- figure7_main_composite_object(
    panel_objects$plots, panel_objects$config
  )
  png_path <- file.path(final_dir, "Figure7_reviewed_GRCh.png")
  pdf_path <- file.path(final_dir, "Figure7_reviewed_GRCh.pdf")
  ggplot2::ggsave(
    png_path, plot = composite, device = "png", dpi = spec$png_dpi,
    width = spec$width_in, height = spec$height_in, units = "in",
    bg = "white", limitsize = FALSE
  )
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
  ggplot2::ggsave(
    pdf_path, plot = composite, device = pdf_device,
    width = spec$width_in, height = spec$height_in, units = "in",
    bg = "white", limitsize = FALSE
  )
  output_paths <- c(png_path, pdf_path)
  if (any(!file.exists(output_paths)) || any(file.info(output_paths)$size <= 0)) {
    stop("Publication-scale Figure 7 export failed", call. = FALSE)
  }
  write_final_records(output_paths)
  invisible(output_paths)
}

panel_objects <- build_panel_objects()
if (args$phase %in% c("subpanels", "all")) write_subpanels(panel_objects)
if (args$phase %in% c("final", "all")) write_final(panel_objects)
message("Completed Figure 7 polishing phase: ", args$phase)
