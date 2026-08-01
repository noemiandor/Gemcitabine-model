#!/usr/bin/env Rscript

# Assemble the reviewed Day-24/Day-31 TGI sensitivity composite from the
# tracked compact endpoint source bundle. The script intentionally consumes
# frozen plot/test tables rather than re-estimating either endpoint analysis.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[[1L]])
} else {
  "Code/in-vivo/figure7/assemble_tgi_sensitivity.R"
}
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
repo_root <- normalizePath(
  file.path(script_dir, "..", "..", ".."),
  mustWork = TRUE
)
for (file in c("common_io.R", "tgi_data.R", "tgi_panels.R")) {
  sys.source(file.path(script_dir, "src", file), envir = .GlobalEnv)
}

required_packages <- c("ggplot2", "ggrepel", "patchwork", "yaml")
missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  logical(1L),
  quietly = TRUE
)]
if (length(missing_packages)) {
  figure7_stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

args <- figure7_parse_args(commandArgs(trailingOnly = TRUE))
source_bundle_id <- "tgi_day24_day31_all_cbs_v1"
default_source_bundle_dir <- file.path(
  repo_root,
  "Data", "in-vivo", "figure7", "saved_tgi_sensitivity",
  source_bundle_id
)
source_bundle_dir <- normalizePath(
  figure7_arg(
    args,
    "source-bundle-dir",
    default_source_bundle_dir
  ),
  mustWork = TRUE
)
day24_run <- normalizePath(
  figure7_arg(
    args,
    "day24-run-dir",
    file.path(source_bundle_dir, "day24")
  ),
  mustWork = TRUE
)
day31_run <- normalizePath(
  figure7_arg(
    args,
    "day31-run-dir",
    file.path(source_bundle_dir, "day31")
  ),
  mustWork = TRUE
)
output_dir <- normalizePath(
  figure7_arg(args, "output-dir", required = TRUE),
  mustWork = FALSE
)
config_path <- normalizePath(
  figure7_arg(
    args,
    "config",
    file.path(script_dir, "figure7_config.yaml")
  ),
  mustWork = TRUE
)
output_basename <- figure7_arg(
  args,
  "output-basename",
  "Figure7_Supplement"
)
provenance_only_value <- tolower(trimws(as.character(figure7_arg(
  args,
  "provenance-only",
  "false"
))))
if (!provenance_only_value %in% c("true", "false")) {
  figure7_stop("--provenance-only must be true or false")
}
provenance_only <- identical(provenance_only_value, "true")
if (!grepl("^[A-Za-z0-9_.-]+$", output_basename)) {
  figure7_stop("--output-basename must be a simple filename stem")
}
canonical_day24_dir <- normalizePath(
  file.path(default_source_bundle_dir, "day24"),
  mustWork = TRUE
)
canonical_day31_dir <- normalizePath(
  file.path(default_source_bundle_dir, "day31"),
  mustWork = TRUE
)
if (identical(output_basename, "Figure7_Supplement") &&
    (!identical(day24_run, canonical_day24_dir) ||
      !identical(day31_run, canonical_day31_dir))) {
  figure7_stop(
    "Canonical SI8 assembly must use the tracked ", source_bundle_id,
    " source bundle"
  )
}
expected_bundle_files <- unlist(lapply(c("day24", "day31"), function(day) {
  c(
    file.path(day, "metadata", "run_config.tsv"),
    file.path(
      day,
      "tables",
      c(
        "panel_7A_plot_data.tsv",
        "panel_7C_plot_data.tsv",
        "panel_7C_test.tsv",
        "panel_7E_plot_data.tsv",
        "panel_7E_test.tsv"
      )
    )
  )
}), use.names = FALSE)
observed_bundle_files <- sort(list.files(
  source_bundle_dir,
  recursive = TRUE,
  all.files = FALSE,
  no.. = TRUE
))
if (!identical(sort(expected_bundle_files), observed_bundle_files)) {
  figure7_stop(
    "SI8 source bundle must contain exactly the 10 consumed endpoint tables ",
    "and two run configs"
  )
}

as_keyed <- function(path, label) {
  table <- figure7_read_tsv(path, c("key", "value"))
  if (anyDuplicated(table$key)) {
    figure7_stop(label, " contains duplicate keys")
  }
  stats::setNames(as.character(table$value), table$key)
}

assert_scalar <- function(observed, expected, label) {
  if (length(observed) != 1L || is.na(observed) ||
      !identical(as.character(observed), as.character(expected))) {
    figure7_stop(
      label,
      " mismatch: expected ", expected,
      ", observed ", paste(observed, collapse = ",")
    )
  }
}

read_endpoint_run <- function(run_dir, day) {
  run_config_path <- file.path(run_dir, "metadata", "run_config.tsv")
  run_config <- as_keyed(run_config_path, paste0("Day-", day, " run config"))
  config <- figure7_read_config(config_path, day)
  assert_scalar(
    run_config[["config_sha256"]],
    figure7_sha256(config_path),
    "Figure 7 config hash"
  )
  assert_scalar(run_config[["module"]], "in_vivo_figure7", "module")
  assert_scalar(run_config[["panel_set"]], "a-f", "panel set")
  assert_scalar(run_config[["tgi_outcome"]], "day", "TGI outcome")
  assert_scalar(run_config[["tgi_day"]], day, "TGI day")
  assert_scalar(
    run_config[["tgi_measure"]],
    figure7_tgi_measure(config),
    "TGI measure"
  )
  assert_scalar(
    run_config[["state_pathway_reference_id"]],
    as.character(config$state_pathways$reviewed_reference_id),
    "reviewed state-pathway reference"
  )
  assert_scalar(
    run_config[["state_pathway_reference_kind"]],
    as.character(config$state_pathways$reviewed_reference_kind),
    "reviewed state-pathway reference kind"
  )
  assert_scalar(
    tolower(run_config[["canonical_publication_allowed"]]),
    "true",
    "canonical-publication flag"
  )
  expected_endpoint_hash <- as.character(
    config$versioned_source_artifacts$panel_k_endpoint_ploidy$sha256
  )
  assert_scalar(
    run_config[["endpoint_ploidy_sha256"]],
    expected_endpoint_hash,
    "endpoint-ploidy hash"
  )
  assert_scalar(
    run_config[["endpoint_ploidy_n_cells"]],
    "14125",
    "endpoint-ploidy source cell count"
  )
  assert_scalar(
    run_config[["endpoint_ploidy_n_files"]],
    "16",
    "endpoint-ploidy source file count"
  )
  assert_scalar(
    run_config[["endpoint_ploidy_score_policy"]],
    "arithmetic_mean_of_all_finite_postprocessed_cell_ploidy_per_cbs_file",
    "endpoint-ploidy score policy"
  )
  assert_scalar(
    run_config[["endpoint_ploidy_mapping_policy"]],
    "exact_sample_growth_curve_harvest_plus_.sps.cbs",
    "endpoint-ploidy mapping policy"
  )
  endpoint_locator <- as.character(run_config[["endpoint_ploidy_input"]])
  if (!nzchar(endpoint_locator) || startsWith(endpoint_locator, "external:")) {
    figure7_stop(
      "Day-", day,
      " source run does not provide a portable endpoint-ploidy locator"
    )
  }
  endpoint_path <- normalizePath(
    file.path(repo_root, endpoint_locator),
    mustWork = TRUE
  )
  endpoint_ploidy <- figure7_read_endpoint_ploidy_table(
    endpoint_path,
    config
  )

  table_path <- function(name) file.path(run_dir, "tables", name)
  a <- figure7_read_tsv(
    table_path("panel_7A_plot_data.tsv"),
    c(
      "sample_id", "initial_ploidy", "dose", "dose_mg", "day",
      "tumor_volume", "tumor_volume_change", "series_type",
      "reference_n", "tgi_outcome", "tgi_day", "tgi_measure",
      "matched_control_summary", "matched_control_group"
    )
  )
  cdata <- figure7_read_tsv(
    table_path("panel_7C_plot_data.tsv"),
    c("sample_id", "initial_ploidy", "dose", figure7_tgi_measure(config))
  )
  ctest <- figure7_read_tsv(
    table_path("panel_7C_test.tsv"),
    c(
      "dose_adjusted_difference_high_minus_low",
      "permutation_p_two_sided", "n_group_low", "n_group_high"
    )
  )
  edata <- figure7_read_tsv(
    table_path("panel_7E_plot_data.tsv"),
    c(
      "sample_id", "initial_ploidy", "dose", "dose_mg",
      "endpoint_ploidy_file", "sample_mean_endpoint_ploidy",
      "n_endpoint_ploidy_cells", "endpoint_ploidy_source_total_cells",
      "endpoint_ploidy_source_file_count", "endpoint_ploidy_source_sha256",
      "endpoint_ploidy_score_policy", "endpoint_ploidy_mapping_policy",
      "terminal_postprocessed_cn_score",
      "terminal_cn_score_within_origin_z",
      "terminal_cn_score_nuisance_residual", "tgi_origin_dose_residual",
      "permutation_stratum", figure7_tgi_measure(config)
    )
  )
  etest <- figure7_read_tsv(
    table_path("panel_7E_test.tsv"),
    c(
      "n", "partial_correlation", "effect_per_within_origin_sd",
      "permutation_p_two_sided", "n_permutations", "permutation_mode",
      "permutation_strata", "score_variable", "score_source_sha256",
      "score_source_n_cells", "score_source_n_files",
      "treated_score_n_cells", "score_aggregation_policy",
      "sample_mapping_policy", "score_standardization",
      "adjustment_terms", "outcome_variable", "plot_x", "plot_y"
    )
  )

  for (table in list(a, cdata, ctest, edata, etest)) {
    assert_scalar(unique(table$tgi_outcome), "day", "table TGI outcome")
    assert_scalar(unique(figure7_numeric(table$tgi_day)), day, "table TGI day")
    assert_scalar(
      unique(table$tgi_measure),
      figure7_tgi_measure(config),
      "table TGI measure"
    )
    assert_scalar(
      unique(table$matched_control_summary),
      "mean",
      "matched-control summary"
    )
    assert_scalar(
      unique(table$matched_control_group),
      "initial_ploidy",
      "matched-control group"
    )
  }

  z_contract <- vapply(
    split(edata$terminal_cn_score_within_origin_z, edata$initial_ploidy),
    function(values) {
      values <- figure7_numeric(values)
      abs(mean(values)) <= 1e-10 && abs(stats::sd(values) - 1) <= 1e-10
    },
    logical(1L)
  )
  endpoint_scores <- split(
    endpoint_ploidy$cells$ploidy,
    endpoint_ploidy$cells$file
  )
  expected_files <- c(
    "2N-A2-0" = "SUM159-2N-30-0_harvest.sps.cbs",
    "2N-A2-L" = "SUM159-2N-30-L_harvest.sps.cbs",
    "2N-A4-R" = "SUM159-2N-120-R_harvest.sps.cbs",
    "2N-A4-RL" = "SUM159-2N-120-RL_harvest.sps.cbs",
    "A6-4N-O" = "SUM159-4N-30-0_harvest.sps.cbs",
    "A6-4N-RR" = "SUM159-4N-30-RR_harvest.sps.cbs",
    "4N-A8-RL" = "SUM159-4N-120-RL_harvest.sps.cbs",
    "4N-A8-RR" = "SUM159-4N-120-RR_harvest.sps.cbs"
  )
  mapped_files <- unname(expected_files[edata$sample_id])
  mapped_scores <- vapply(mapped_files, function(file) {
    mean(endpoint_scores[[file]])
  }, numeric(1L))
  mapped_counts <- vapply(mapped_files, function(file) {
    length(endpoint_scores[[file]])
  }, integer(1L))
  expected_results <- list(
    `24` = c(
      effect = -6.77351906237078,
      partial_r = -0.30605912538695,
      permutation_p = 0.5625
    ),
    `31` = c(
      effect = -2.77364130329883,
      partial_r = -0.173200987198346,
      permutation_p = 0.6875
    )
  )[[as.character(day)]]
  if (nrow(edata) != 8L || nrow(etest) != 1L || !all(z_contract) ||
      anyNA(mapped_files) ||
      any(edata$endpoint_ploidy_file != mapped_files) ||
      any(abs(
        figure7_numeric(edata$sample_mean_endpoint_ploidy) - mapped_scores
      ) > 1e-12) ||
      any(abs(
        figure7_numeric(edata$terminal_postprocessed_cn_score) - mapped_scores
      ) > 1e-12) ||
      any(figure7_numeric(edata$n_endpoint_ploidy_cells) != mapped_counts) ||
      sum(figure7_numeric(edata$n_endpoint_ploidy_cells)) != 7623L ||
      any(figure7_numeric(edata$endpoint_ploidy_source_total_cells) != 14125L) ||
      any(figure7_numeric(edata$endpoint_ploidy_source_file_count) != 16L) ||
      any(edata$endpoint_ploidy_source_sha256 != expected_endpoint_hash) ||
      any(edata$endpoint_ploidy_score_policy !=
        run_config[["endpoint_ploidy_score_policy"]]) ||
      any(edata$endpoint_ploidy_mapping_policy !=
        run_config[["endpoint_ploidy_mapping_policy"]]) ||
      any(edata$permutation_stratum != paste(
        edata$initial_ploidy,
        edata$dose_mg,
        sep = "|"
      )) ||
      !identical(
        as.character(etest$permutation_mode),
        "exact_TGI_label_enumeration_within_initial_ploidy_x_dose"
      ) ||
      !identical(
        as.character(etest$permutation_strata),
        "initial_ploidy:dose_mg"
      ) ||
      !identical(
        as.character(etest$score_variable),
        "sample_mean_all_canonical_cbs_cell_ploidy"
      ) ||
      !identical(
        as.character(etest$score_source_sha256),
        expected_endpoint_hash
      ) ||
      figure7_numeric(etest$score_source_n_cells) != 14125L ||
      figure7_numeric(etest$score_source_n_files) != 16L ||
      figure7_numeric(etest$treated_score_n_cells) != 7623L ||
      !identical(
        as.character(etest$score_aggregation_policy),
        run_config[["endpoint_ploidy_score_policy"]]
      ) ||
      !identical(
        as.character(etest$sample_mapping_policy),
        run_config[["endpoint_ploidy_mapping_policy"]]
      ) ||
      !identical(
        as.character(etest$score_standardization),
        "z_score_within_initial_ploidy"
      ) ||
      !identical(
        as.character(etest$adjustment_terms),
        "initial_ploidy+dose_mg"
      ) ||
      !identical(
        as.character(etest$plot_x),
        "terminal_cn_score_nuisance_residual"
      ) ||
      !identical(as.character(etest$plot_y), "tgi_origin_dose_residual") ||
      figure7_numeric(etest$n_permutations) != 16L ||
      abs(figure7_numeric(etest$effect_per_within_origin_sd) -
        expected_results[["effect"]]) > 1e-12 ||
      abs(figure7_numeric(etest$partial_correlation) -
        expected_results[["partial_r"]]) > 1e-12 ||
      abs(figure7_numeric(etest$permutation_p_two_sided) -
        expected_results[["permutation_p"]]) > 1e-12) {
    figure7_stop(
      "Day-", day,
      " panel K lacks the reviewed confound-safe analysis contract"
    )
  }

  list(
    run_dir = run_dir,
    run_config_path = run_config_path,
    endpoint_path = endpoint_path,
    config = config,
    a = a,
    cdata = cdata,
    ctest = ctest,
    edata = edata,
    etest = etest,
    source_tables = c(
      panel_a = table_path("panel_7A_plot_data.tsv"),
      panel_c_data = table_path("panel_7C_plot_data.tsv"),
      panel_c_test = table_path("panel_7C_test.tsv"),
      panel_e_data = table_path("panel_7E_plot_data.tsv"),
      panel_e_test = table_path("panel_7E_test.tsv")
    )
  )
}

day24 <- read_endpoint_run(day24_run, 24L)
day31 <- read_endpoint_run(day31_run, 31L)

growth_identity_columns <- c(
  "sample_id", "initial_ploidy", "dose", "dose_mg", "day",
  "tumor_volume", "tumor_volume_change", "series_type", "reference_n"
)
growth24 <- day24$a[, growth_identity_columns, drop = FALSE]
growth31 <- day31$a[, growth_identity_columns, drop = FALSE]
rownames(growth24) <- NULL
rownames(growth31) <- NULL
if (!identical(growth24, growth31)) {
  figure7_stop("Day-24 and Day-31 runs disagree on the frozen growth data")
}

figure7_sensitivity_growth_plot <- function(data) {
  individual <- data[data$series_type == "individual mouse", , drop = FALSE]
  reference <- data[
    data$series_type == "matched-control reference",
    ,
    drop = FALSE
  ]
  individual$dose <- factor(
    individual$dose,
    levels = c("0mg/kg", "30mg/kg", "120mg/kg")
  )
  endpoints <- c(24L, 31L)
  endpoint_individual <- individual[individual$day %in% endpoints, , drop = FALSE]
  endpoint_reference <- reference[reference$day %in% endpoints, , drop = FALSE]
  endpoint_individual$endpoint <- factor(
    endpoint_individual$day,
    levels = endpoints,
    labels = paste("Day", endpoints)
  )
  endpoint_reference$endpoint <- factor(
    endpoint_reference$day,
    levels = endpoints,
    labels = paste("Day", endpoints)
  )

  ggplot2::ggplot() +
    ggplot2::geom_vline(
      xintercept = endpoints[[1L]],
      color = "#B2182B",
      linetype = "22",
      linewidth = 0.9,
      alpha = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = endpoints[[2L]],
      color = "#B2182B",
      linetype = "solid",
      linewidth = 0.9,
      alpha = 0.8
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "grey75",
      linewidth = 0.35
    ) +
    ggplot2::geom_line(
      data = individual,
      ggplot2::aes(
        day, tumor_volume_change,
        group = sample_id,
        color = dose
      ),
      linewidth = 0.55,
      alpha = 0.72
    ) +
    ggplot2::geom_line(
      data = reference,
      ggplot2::aes(day, tumor_volume_change, group = initial_ploidy),
      color = "black",
      linetype = "22",
      linewidth = 1.15
    ) +
    ggplot2::geom_point(
      data = endpoint_individual,
      ggplot2::aes(
        day, tumor_volume_change,
        color = dose,
        shape = endpoint
      ),
      fill = "white",
      stroke = 1,
      size = 2.5
    ) +
    ggplot2::geom_point(
      data = endpoint_reference,
      ggplot2::aes(day, tumor_volume_change, shape = endpoint),
      fill = "#B2182B",
      color = "black",
      size = 3.1
    ) +
    ggplot2::facet_wrap(~initial_ploidy, nrow = 1) +
    ggplot2::scale_color_manual(
      values = figure7_dose_colors(),
      breaks = c("0mg/kg", "30mg/kg", "120mg/kg"),
      name = "Dose"
    ) +
    ggplot2::scale_shape_manual(
      values = c("Day 24" = 21, "Day 31" = 24),
      name = "TGI endpoint"
    ) +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(individual$day)),
      minor_breaks = NULL
    ) +
    ggplot2::labs(
      title = "How Day 24 and Day 31 TGI are calculated from tumor growth",
      subtitle = paste0(
        "At each marked endpoint: TGI = 100 x (1 - treated growth delta / ",
        "mean injected-origin-matched control growth delta)"
      ),
      x = "Days since first treatment",
      y = expression(Delta * " tumor volume from Day 0 (mm"^3 * ")"),
      caption = paste0(
        "Thin lines: individual mice. Dashed black: mean injected-origin-",
        "matched untreated controls. Red dashed: Day 24; red solid: Day 31."
      )
    ) +
    figure7_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.caption = ggplot2::element_text(hjust = 0)
    )
}

plots <- list(
  A = figure7_sensitivity_growth_plot(growth24),
  B = figure7_panel_c_plot(day24$cdata, day24$ctest, day24$config),
  C = figure7_adjusted_cn_plot(day24$edata, day24$etest, day24$config),
  D = figure7_panel_c_plot(day31$cdata, day31$ctest, day31$config),
  E = figure7_adjusted_cn_plot(day31$edata, day31$etest, day31$config)
)
design <- paste("AAAA", "BBCC", "DDEE", sep = "\n")
composite <- patchwork::wrap_plots(plots, design = design) +
  patchwork::plot_layout(heights = c(1.05, 1, 1)) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1)
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
png_path <- file.path(output_dir, paste0(output_basename, ".png"))
pdf_path <- file.path(output_dir, paste0(output_basename, ".pdf"))
panel_a_path <- file.path(
  output_dir,
  "panel_7A_day31_day_24_tgi_calculation.png"
)
if (!provenance_only) {
  ggplot2::ggsave(
    panel_a_path,
    plots$A,
    device = "png",
    dpi = 300,
    width = 15,
    height = 8,
    units = "in",
    bg = "white"
  )
  ggplot2::ggsave(
    png_path,
    composite,
    device = "png",
    dpi = 300,
    width = 16,
    height = 24,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
  ggplot2::ggsave(
    pdf_path,
    composite,
    device = grDevices::cairo_pdf,
    width = 16,
    height = 24,
    units = "in",
    limitsize = FALSE
  )
}
for (path in c(panel_a_path, png_path, pdf_path)) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    figure7_stop("Failed to write sensitivity asset: ", path)
  }
}

repo_locator <- function(path) {
  normalized <- normalizePath(path, mustWork = TRUE)
  prefix <- paste0(repo_root, .Platform$file.sep)
  if (startsWith(normalized, prefix)) {
    substring(normalized, nchar(prefix) + 1L)
  } else {
    paste0("external:", basename(normalized))
  }
}
source_files <- c(
  day24_run_config = day24$run_config_path,
  day24_endpoint_ploidy = day24$endpoint_path,
  stats::setNames(
    day24$source_tables,
    paste0("day24_", names(day24$source_tables))
  ),
  day31_run_config = day31$run_config_path,
  day31_endpoint_ploidy = day31$endpoint_path,
  stats::setNames(
    day31$source_tables,
    paste0("day31_", names(day31$source_tables))
  )
)
if (is.null(names(source_files)) || any(!nzchar(names(source_files))) ||
    anyDuplicated(names(source_files))) {
  figure7_stop("SI8 source provenance roles must be nonempty and unique")
}
provenance <- data.frame(
  key = c(
    "artifact", "assembly_script", "assembly_script_sha256",
    "source_bundle_id", "source_bundle", "source_bundle_file_count",
    "day24_source_run", "day31_source_run",
    "endpoint_ploidy_source_n_cells",
    "treated_endpoint_ploidy_n_cells",
    "endpoint_ploidy_score_policy",
    paste0("source_file:", names(source_files)),
    paste0("source_sha256:", names(source_files)),
    "day24_adjusted_slope_per_origin_sd", "day24_partial_r",
    "day24_exact_permutation_p", "day31_adjusted_slope_per_origin_sd",
    "day31_partial_r", "day31_exact_permutation_p",
    "png", "png_sha256", "pdf", "pdf_sha256"
  ),
  value = c(
    "supplementary_figure_8_tgi_sensitivity",
    repo_locator(script_path),
    figure7_sha256(script_path),
    source_bundle_id,
    repo_locator(source_bundle_dir),
    as.character(length(expected_bundle_files)),
    repo_locator(day24_run),
    repo_locator(day31_run),
    "14125",
    "7623",
    "arithmetic_mean_of_all_finite_postprocessed_cell_ploidy_per_cbs_file",
    vapply(source_files, repo_locator, character(1L)),
    vapply(source_files, figure7_sha256, character(1L)),
    as.character(day24$etest$effect_per_within_origin_sd),
    as.character(day24$etest$partial_correlation),
    as.character(day24$etest$permutation_p_two_sided),
    as.character(day31$etest$effect_per_within_origin_sd),
    as.character(day31$etest$partial_correlation),
    as.character(day31$etest$permutation_p_two_sided),
    repo_locator(png_path),
    figure7_sha256(png_path),
    repo_locator(pdf_path),
    figure7_sha256(pdf_path)
  ),
  stringsAsFactors = FALSE
)
if (anyDuplicated(provenance$key)) {
  figure7_stop("SI8 provenance keys must be unique")
}
provenance_path <- file.path(
  output_dir,
  paste0(output_basename, "_provenance.tsv")
)
figure7_write_tsv(provenance, provenance_path)

manifest_columns <- c(
  "figure", "panel", "asset_path", "source_file", "source_kind",
  "generated_by", "command", "input_data", "result_run_dir", "run_id",
  "caption_role", "asset_status", "not_regenerated_reason",
  "local_provenance_path", "citation_or_uri", "notes"
)
asset_paths <- c(pdf = pdf_path, png = png_path)
asset_locators <- vapply(asset_paths, repo_locator, character(1L))
provenance_locator <- repo_locator(provenance_path)
bundle_locator <- repo_locator(source_bundle_dir)
assembly_command <- paste(
  "Rscript Code/in-vivo/figure7/assemble_tgi_sensitivity.R",
  paste0("--source-bundle-dir=", bundle_locator),
  paste0("--output-dir=", repo_locator(output_dir))
)
si8_manifest <- data.frame(
  figure = rep("Supplementary", 2L),
  panel = c("SuppFig8", "SuppFig8_png"),
  asset_path = unname(asset_locators),
  source_file = unname(asset_locators),
  source_kind = rep("generated_panel", 2L),
  generated_by = rep(
    "Code/in-vivo/figure7/assemble_tgi_sensitivity.R",
    2L
  ),
  command = rep(assembly_command, 2L),
  input_data = unname(asset_locators),
  result_run_dir = rep(repo_locator(output_dir), 2L),
  run_id = rep(source_bundle_id, 2L),
  caption_role = c(
    "Supplementary Figure 8 final Day-24/Day-31 TGI sensitivity composite",
    "PNG derivative of Supplementary Figure 8 final composite"
  ),
  asset_status = rep("generated", 2L),
  not_regenerated_reason = rep("", 2L),
  local_provenance_path = rep(provenance_locator, 2L),
  citation_or_uri = rep("", 2L),
  notes = rep(
    paste0(
      "portable_source_bundle=", bundle_locator,
      ";canonical_endpoint_ploidy_sha256=",
      as.character(
        day24$config$versioned_source_artifacts$panel_k_endpoint_ploidy$sha256
      )
    ),
    2L
  ),
  stringsAsFactors = FALSE
)
if (!identical(names(si8_manifest), manifest_columns) ||
    nrow(si8_manifest) != 2L || anyDuplicated(si8_manifest$panel) ||
    any(si8_manifest$local_provenance_path != provenance_locator)) {
  figure7_stop("SI8 figure manifest contract is invalid")
}
figure7_write_tsv(
  si8_manifest,
  file.path(output_dir, "si8_manifest.tsv")
)
if (provenance_only) {
  message(
    "Updated reviewed Day-24/Day-31 TGI sensitivity provenance without ",
    "rewriting composite assets: ",
    file.path(output_dir, paste0(output_basename, "_provenance.tsv"))
  )
} else {
  message("Assembled reviewed Day-24/Day-31 TGI sensitivity composite: ", png_path)
}
