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
for (file in c(
  "common_io.R", "tgi_data.R", "tgi_statistics.R", "tgi_panels.R"
)) {
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
source_bundle_id <- "tgi_day24_day31_curated_cbs_v3_raw_pearson"
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
  # The compact Day-24/Day-31 source bundle records the exact analysis-time
  # config.  Main-composite presentation edits must not rewrite that frozen
  # provenance, so validate the recorded source hash directly and validate
  # every consumed scientific field against the current config below.
  assert_scalar(
    run_config[["config_sha256"]],
    "50dfdbc23c726f5a249a3ac40f4e0be877a944ec36aa512969a8b057d6b08bd2",
    "reviewed source-analysis config hash"
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
  # SI8 consumes only the Day-24/Day-31 TGI and endpoint-ploidy tables. Its
  # byte-pinned source-run config also records the panel-7F reference that was
  # current when those unrelated tables were generated, but that identity is
  # not a scientific dependency of this sensitivity analysis and therefore
  # must not be coupled to the current panel-7F reference version.
  assert_scalar(
    tolower(run_config[["canonical_publication_allowed"]]),
    "true",
    "canonical-publication flag"
  )
  expected_endpoint_hash <- as.character(
    config$versioned_source_artifacts$panel_l_endpoint_ploidy$sha256
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
    run_config[["endpoint_ploidy_score_universe_n_cells"]],
    "9832",
    "endpoint-ploidy score-universe cell count"
  )
  assert_scalar(
    run_config[["endpoint_ploidy_treated_score_n_cells"]],
    "5335",
    "treated endpoint-ploidy score cell count"
  )
  assert_scalar(
    run_config[["endpoint_ploidy_score_policy"]],
    paste0(
      "arithmetic_mean_of_finite_postprocessed_cell_ploidy_in_exact_",
      "qc_passed_cellcycle_noncellcycle_union_per_sample"
    ),
    "endpoint-ploidy score policy"
  )
  assert_scalar(
    run_config[["endpoint_ploidy_mapping_policy"]],
    paste0(
      "exact_processed_sample_barcode_to_canonical_cbs_file_cell_and_value;",
      "score_universe=reviewed_final_seurat_tumor_cells"
    ),
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
  resolve_processed_input <- function(role) {
    locator <- as.character(run_config[[paste0(role, "_input")]])
    expected_hash <- as.character(run_config[[paste0(role, "_sha256")]])
    if (!nzchar(locator) || startsWith(locator, "external:") ||
        !grepl("^[0-9a-f]{64}$", expected_hash)) {
      figure7_stop(
        "Day-", day, " source run does not provide a portable, hashed ",
        role, " locator"
      )
    }
    path <- normalizePath(file.path(repo_root, locator), mustWork = TRUE)
    figure7_verify_checksum(path, expected_hash, paste0("Day-", day, " ", role))
    path
  }
  cellcycle_path <- resolve_processed_input("cellcycle")
  noncellcycle_path <- resolve_processed_input("noncellcycle")
  expected_samples <- figure7_sample_table(
    figure7_read_cell_table(cellcycle_path, "CellCycle", config),
    figure7_read_cell_table(noncellcycle_path, "NonCellCycle", config),
    config,
    endpoint_ploidy
  )
  expected_treated <- expected_samples[
    expected_samples$dose_mg > 0,
    ,
    drop = FALSE
  ]

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
      "endpoint_ploidy_score_universe_total_cells",
      "endpoint_ploidy_source_file_count", "endpoint_ploidy_source_sha256",
      "endpoint_ploidy_score_policy", "endpoint_ploidy_mapping_policy",
      figure7_tgi_measure(config)
    )
  )
  etest <- figure7_read_tsv(
    table_path("panel_7E_test.tsv"),
    c(
      "n", "estimate", "asymptotic_p",
      "permutation_p_two_sided", "n_permutations", "permutation_mode",
      "permutation_strata", "association_type", "score_variable",
      "score_source_sha256",
      "score_source_n_cells", "score_inventory_n_cells",
      "score_source_n_files",
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

  expected_treated <- expected_treated[
    match(edata$sample_id, expected_treated$sample_id),
    ,
    drop = FALSE
  ]
  mapped_files <- expected_treated$endpoint_ploidy_file
  mapped_scores <- expected_treated$sample_mean_endpoint_ploidy
  mapped_counts <- expected_treated$n_endpoint_ploidy_cells
  expected_results <- list(
    `24` = c(
      estimate = -0.387348978978662,
      asymptotic_p = 0.343097626516948,
      permutation_p = 0.346924603174603
    ),
    `31` = c(
      estimate = -0.296948455789288,
      asymptotic_p = 0.475086351653848,
      permutation_p = 0.485714285714286
    )
  )[[as.character(day)]]
  if (nrow(edata) != 8L || nrow(etest) != 1L) {
    figure7_stop(
      "Day-", day,
      " panel L must contain exactly eight treated tumors and one test row"
    )
  }
  recomputed <- figure7_exact_cor(
    figure7_numeric(edata$sample_mean_endpoint_ploidy),
    figure7_numeric(edata[[figure7_tgi_measure(config)]])
  )
  if (anyNA(mapped_files) ||
      any(edata$endpoint_ploidy_file != mapped_files) ||
      any(abs(
        figure7_numeric(edata$sample_mean_endpoint_ploidy) - mapped_scores
      ) > 1e-12) ||
      any(figure7_numeric(edata$n_endpoint_ploidy_cells) != mapped_counts) ||
      sum(figure7_numeric(edata$n_endpoint_ploidy_cells)) != 5335L ||
      any(figure7_numeric(edata$endpoint_ploidy_source_total_cells) != 14125L) ||
      any(figure7_numeric(edata$endpoint_ploidy_score_universe_total_cells) != 9832L) ||
      any(figure7_numeric(edata$endpoint_ploidy_source_file_count) != 16L) ||
      any(edata$endpoint_ploidy_source_sha256 != expected_endpoint_hash) ||
      any(edata$endpoint_ploidy_score_policy !=
        run_config[["endpoint_ploidy_score_policy"]]) ||
      any(edata$endpoint_ploidy_mapping_policy !=
        run_config[["endpoint_ploidy_mapping_policy"]]) ||
      !identical(
        as.character(etest$permutation_mode),
        "exact_TGI_label_enumeration"
      ) ||
      !identical(
        as.character(etest$permutation_strata),
        "none"
      ) ||
      !identical(
        as.character(etest$association_type),
        "unadjusted_mouse_level_pearson"
      ) ||
      !identical(
        as.character(etest$score_variable),
        "sample_mean_qc_passed_curated_cbs_cell_ploidy"
      ) ||
      !identical(
        as.character(etest$score_source_sha256),
        expected_endpoint_hash
      ) ||
      figure7_numeric(etest$score_source_n_cells) != 9832L ||
      figure7_numeric(etest$score_inventory_n_cells) != 14125L ||
      figure7_numeric(etest$score_source_n_files) != 16L ||
      figure7_numeric(etest$treated_score_n_cells) != 5335L ||
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
        "none"
      ) ||
      !identical(
        as.character(etest$adjustment_terms),
        "none"
      ) ||
      !identical(
        as.character(etest$plot_x),
        "sample_mean_endpoint_ploidy"
      ) ||
      !identical(
        as.character(etest$plot_y),
        figure7_tgi_measure(config)
      ) ||
      figure7_numeric(etest$n_permutations) != factorial(8L) ||
      abs(figure7_numeric(etest$estimate) - recomputed$estimate) > 1e-12 ||
      abs(figure7_numeric(etest$asymptotic_p) -
        recomputed$asymptotic_p) > 1e-12 ||
      abs(figure7_numeric(etest$permutation_p_two_sided) -
        recomputed$permutation_p_two_sided) > 1e-12 ||
      abs(figure7_numeric(etest$estimate) -
        expected_results[["estimate"]]) > 1e-12 ||
      abs(figure7_numeric(etest$asymptotic_p) -
        expected_results[["asymptotic_p"]]) > 1e-12 ||
      abs(figure7_numeric(etest$permutation_p_two_sided) -
        expected_results[["permutation_p"]]) > 1e-12) {
    figure7_stop(
      "Day-", day,
      " panel L lacks the reviewed raw endpoint-ploidy association contract"
    )
  }

  list(
    run_dir = run_dir,
    run_config_path = run_config_path,
    endpoint_path = endpoint_path,
    cellcycle_path = cellcycle_path,
    noncellcycle_path = noncellcycle_path,
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
    ggplot2::facet_wrap(
      ~initial_ploidy,
      nrow = 1,
      labeller = ggplot2::as_labeller(figure7_sum159_origin_labels())
    ) +
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
      title = NULL,
      subtitle = paste0(
        "At each marked endpoint:\nTGI = 100 x (1 - treated growth delta / ",
        "mean injected-origin-matched control growth delta)"
      ),
      x = "Days since first treatment",
      y = "Tumor size change\n(mm³)",
      caption = paste0(
        "Thin lines: individual mice. Dashed black: mean injected-origin-",
        "matched untreated controls.\nRed dashed: Day 24; red solid: Day 31."
      )
    ) +
    figure7_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.caption = ggplot2::element_text(hjust = 0)
    )
}

sensitivity_tgi_group_plot <- function(endpoint) {
  figure7_panel_c_plot(endpoint$cdata, endpoint$ctest, endpoint$config) +
    ggplot2::scale_color_manual(
      values = figure7_dose_colors(),
      breaks = c("30mg/kg", "120mg/kg"),
      labels = c("30 mg/kg", "120 mg/kg"),
      name = "Dose"
    ) +
    ggplot2::scale_shape_manual(
      values = c("30mg/kg" = 16, "120mg/kg" = 17),
      breaks = c("30mg/kg", "120mg/kg"),
      labels = c("30 mg/kg", "120 mg/kg"),
      name = "Dose"
    ) +
    ggplot2::labs(
      title = paste0(
        "Day ", figure7_tgi_day(endpoint$config),
        " TGI by\ninitial ploidy"
      )
    )
}

sensitivity_endpoint_ploidy_plot <- function(endpoint) {
  day <- figure7_tgi_day(endpoint$config)
  figure7_endpoint_ploidy_plot(endpoint$edata, endpoint$etest, endpoint$config) +
    ggplot2::scale_color_manual(
      values = figure7_dose_colors(),
      labels = c("30mg/kg" = "30 mg/kg", "120mg/kg" = "120 mg/kg"),
      name = "Dose"
    ) +
    ggplot2::scale_shape_manual(
      values = c("2N" = 16, "4N" = 17),
      breaks = names(figure7_sum159_origin_labels()),
      labels = unname(figure7_sum159_origin_labels()),
      name = "Injected origin"
    ) +
    ggplot2::labs(
      title = paste0("Day ", day, " TGI vs endpoint\ntumor-cell ploidy"),
      subtitle = paste0(
        "5,335 QC-passed treated-tumor CBS cells;\n",
        "unadjusted mouse-level association"
      )
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(nrow = 1),
      color = ggplot2::guide_legend(nrow = 1)
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "left"
    )
}

plots <- list(
  A = figure7_sensitivity_growth_plot(growth24),
  B = sensitivity_tgi_group_plot(day24),
  C = sensitivity_endpoint_ploidy_plot(day24),
  D = sensitivity_tgi_group_plot(day31),
  E = sensitivity_endpoint_ploidy_plot(day31)
)
design <- paste("AAAA", "BBCC", "DDEE", sep = "\n")
composite <- patchwork::wrap_plots(plots, design = design) +
  patchwork::plot_layout(heights = c(1.05, 1, 1)) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 12),
    plot.tag.position = c(0, 1),
    axis.title = ggplot2::element_text(size = 9),
    axis.text = ggplot2::element_text(size = 7.5),
    strip.text = ggplot2::element_text(size = 8, face = "bold"),
    legend.title = ggplot2::element_text(size = 8),
    legend.text = ggplot2::element_text(size = 7.5),
    plot.title = ggplot2::element_text(size = 9.5),
    plot.subtitle = ggplot2::element_text(size = 7.5)
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
    width = 7.1,
    height = 3.8,
    units = "in",
    bg = "white"
  )
  ggplot2::ggsave(
    png_path,
    composite,
    device = "png",
    dpi = 300,
    width = 7.1,
    height = 11.2,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
  ggplot2::ggsave(
    pdf_path,
    composite,
    device = grDevices::cairo_pdf,
    width = 7.1,
    height = 11.2,
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
if (!identical(day24$cellcycle_path, day31$cellcycle_path) ||
    !identical(day24$noncellcycle_path, day31$noncellcycle_path)) {
  figure7_stop(
    "Day-24 and Day-31 sensitivity sources must bind the same exact ",
    "QC-passed CellCycle + NonCellCycle score universe"
  )
}
source_files <- c(
  cellcycle = day24$cellcycle_path,
  noncellcycle = day24$noncellcycle_path,
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
    "endpoint_ploidy_inventory_n_cells",
    "endpoint_ploidy_score_universe_n_cells",
    "treated_endpoint_ploidy_n_cells",
    "endpoint_ploidy_score_policy",
    paste0("source_file:", names(source_files)),
    paste0("source_sha256:", names(source_files)),
    "day24_pearson_r", "day24_asymptotic_p",
    "day24_exact_unrestricted_permutation_p", "day31_pearson_r",
    "day31_asymptotic_p", "day31_exact_unrestricted_permutation_p",
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
    "9832",
    "5335",
    paste0(
      "arithmetic_mean_of_finite_postprocessed_cell_ploidy_in_exact_",
      "qc_passed_cellcycle_noncellcycle_union_per_sample"
    ),
    vapply(source_files, repo_locator, character(1L)),
    vapply(source_files, figure7_sha256, character(1L)),
    as.character(day24$etest$estimate),
    as.character(day24$etest$asymptotic_p),
    as.character(day24$etest$permutation_p_two_sided),
    as.character(day31$etest$estimate),
    as.character(day31$etest$asymptotic_p),
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
        day24$config$versioned_source_artifacts$panel_l_endpoint_ploidy$sha256
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
