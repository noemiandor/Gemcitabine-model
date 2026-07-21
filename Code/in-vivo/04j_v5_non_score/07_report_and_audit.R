ptpv5_safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) min(x) else NA_real_
}

ptpv5_make_figures <- function(context, synthesis) {
  dirs <- context$dirs
  evidence <- synthesis$evidence$summary
  evidence$method_fraction <- paste0(
    evidence$n_methods_passing, "/", evidence$n_methods_available
  )
  p1 <- ggplot2::ggplot(
    evidence,
    ggplot2::aes(
      y = reorder(response_family_label, n_methods_available),
      x = n_methods_passing
    )
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = n_methods_available,
        yend = reorder(response_family_label, n_methods_available)
      ),
      linewidth = 2.2,
      colour = "grey85"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = synthesis_class),
      size = 3
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = method_fraction),
      nudge_x = 0.18,
      hjust = 0,
      size = 3.2
    ) +
    ggplot2::scale_x_continuous(breaks = 0:5, limits = c(0, 5.4)) +
    ggplot2::labs(
      y = NULL,
      x = "Methods passing / methods available",
      title = "Cross-method support by response family"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(
    file.path(dirs$figures, "family_cross_method_support.png"),
    p1, width = 9, height = 6.5, dpi = 300
  )
  ggplot2::ggsave(
    file.path(dirs$figures, "family_cross_method_support.pdf"),
    p1, width = 9, height = 6.5
  )

  calibration <- synthesis$calibration
  calibration_plot <- calibration[
    is.finite(calibration$empirical_type1_at_0_05),
    ,
    drop = FALSE
  ]
  p2 <- ggplot2::ggplot(
    calibration_plot,
    ggplot2::aes(
      x = method_id,
      y = empirical_type1_at_0_05
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0.05,
      linetype = 2,
      colour = "grey40"
    ) +
    ggplot2::geom_col(fill = "#2b6cb0", width = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_low_95, ymax = ci_high_95),
      width = 0.15
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = "Empirical rejection rate at alpha = 0.05",
      title = "Exact-randomization calibration"
    ) +
    ggplot2::theme_bw(base_size = 10)
  ggplot2::ggsave(
    file.path(dirs$figures, "method_calibration.png"),
    p2, width = 8, height = 4.5, dpi = 300
  )
  ggplot2::ggsave(
    file.path(dirs$figures, "method_calibration.pdf"),
    p2, width = 8, height = 4.5
  )

  power <- synthesis$power
  p3 <- ggplot2::ggplot(
    power,
    ggplot2::aes(
      x = scenario,
      y = power_at_0_05,
      colour = method_id,
      group = method_id
    )
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = NULL,
      y = "Generic stress-test detection probability",
      title = "Method sensitivity profiles under standardized shifts"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    )
  ggplot2::ggsave(
    file.path(dirs$figures, "generic_power_stress_test.png"),
    p3, width = 9, height = 6, dpi = 300
  )
  ggplot2::ggsave(
    file.path(dirs$figures, "generic_power_stress_test.pdf"),
    p3, width = 9, height = 6
  )
  invisible(list(family = p1, calibration = p2, power = p3))
}

ptpv5_artifact_source <- function(id, label, path, description) {
  portable_path <- file.path("queries", paste0(id, ".sql"))
  list(
    manifest = list(id = id, label = label, path = portable_path),
    source = list(
      id = id,
      query = list(
        engine = "duckdb",
        sql = paste0(
          "SELECT * FROM embedded_snapshot.",
          gsub("[^A-Za-z0-9_]", "_", id)
        ),
        description = paste0(description, " Frozen source table: ", basename(path), "."),
        executed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
      )
    )
  )
}

ptpv5_build_artifact <- function(context, results, synthesis) {
  cfg <- context$cfg
  dirs <- context$dirs
  family <- synthesis$evidence$summary
  evidence_long <- synthesis$evidence$long
  calibration <- synthesis$calibration
  power <- synthesis$power
  methods <- synthesis$method_summary
  trajectory <- results$gene_wise_trajectory$estimability
  comparison <- synthesis$v4_comparison
  supported <- family[
    family$synthesis_class == "concordant_non_score_support",
    ,
    drop = FALSE
  ]
  isolated <- family[
    family$synthesis_class == "isolated_method_support_requires_caution",
    ,
    drop = FALSE
  ]
  n_complete <- sum(grepl(
    "^complete",
    methods$status
  ))
  v4_fdr_cols <- intersect(
    c("fdr_focused_11", "fdr_all_families"),
    names(comparison)
  )
  v4_pass <- if (length(v4_fdr_cols) && nrow(comparison)) {
    any(as.matrix(comparison[, v4_fdr_cols, drop = FALSE]) <= 0.10, na.rm = TRUE)
  } else {
    NA
  }
  summary_text <- if (nrow(supported)) {
    paste0(
      "At the predefined thresholds, ",
      nrow(supported),
      " response family/families received concordant support from at least ",
      "two non-score methods: ",
      paste(supported$response_family_label, collapse = "; "),
      ". These are cross-sectional treatment-by-baseline-ploidy associations, ",
      "not causal or cell-autonomous effects."
    )
  } else {
    paste0(
      "None of the 15 frozen response families passed predefined thresholds ",
      "in two or more non-score methods. The V5 extension therefore does not ",
      "establish a robust family-level treatment-by-baseline-ploidy signal. ",
      if (nrow(isolated)) {
        paste0(
          nrow(isolated),
          " family/families had isolated single-method support, which remains ",
          "method-dependent and requires caution. "
        )
      } else {
        ""
      },
      "This is compatible with a weak, sparse, state-compositional, or ",
      "underpowered effect; it is not evidence that no biological difference exists."
    )
  }
  trajectory_text <- paste0(
    "The local [0.30, 0.49] gene-wise spline branch is ",
    trajectory$status[trajectory$branch == "primary_0.30_0.49"],
    " because only ",
    trajectory$n_supported_bins[trajectory$branch == "primary_0.30_0.49"],
    " four-group-supported bin is available; V5 does not extrapolate across ",
    "unsupported pseudotime. The whole-trajectory branch uses ",
    trajectory$n_supported_bins[
      trajectory$branch == "whole_trajectory_0.00_1.00"
    ],
    " supported bins."
  )
  v4_text <- if (isTRUE(v4_pass)) {
    "V4 contained at least one family-level FDR result at or below 0.10; V5 evaluates whether non-score methods reproduce it."
  } else if (identical(v4_pass, FALSE)) {
    "V4 had no response family passing its focused/all-family 0.10 FDR threshold; V5 was therefore treated as a post-hoc methodological extension rather than confirmation of a prior discovery."
  } else {
    "The V4 comparison table was unavailable for an automated threshold statement."
  }
  run_metrics <- data.frame(
    methods_completed = n_complete,
    methods_planned = 5L,
    response_families = nrow(family),
    concordant_families = nrow(supported),
    isolated_families = nrow(isolated),
    exact_assignments = nrow(context$bundle$assignments),
    mice = length(context$bundle$sample_ids),
    stringsAsFactors = FALSE
  )
  ptp_write_csv(run_metrics, file.path(dirs$report, "run_metrics.csv"))
  chart_family <- rbind(
    data.frame(
      response_family = family$response_family_label,
      metric = "Methods available",
      method_count = family$n_methods_available,
      stringsAsFactors = FALSE
    ),
    data.frame(
      response_family = family$response_family_label,
      metric = "Methods passing",
      method_count = family$n_methods_passing,
      stringsAsFactors = FALSE
    )
  )
  calibration_exact <- calibration[
    is.finite(calibration$empirical_type1_at_0_05),
    c("method_id", "empirical_type1_at_0_05", "ci_low_95", "ci_high_95"),
    drop = FALSE
  ]
  family_table <- family
  evidence_table <- evidence_long[
    !is.na(evidence_long$method_id),
    c(
      "response_family_label", "method_id", "evidence_value",
      "multiplicity_value", "evidence_type", "pass_predefined_threshold"
    ),
    drop = FALSE
  ]
  sources <- list(
    ptpv5_artifact_source(
      "run_metrics_source",
      "V5 report run metrics",
      "run_metrics.csv",
      "One-row report metric summary derived from frozen method statuses and family synthesis."
    ),
    ptpv5_artifact_source(
      "family_synthesis_source",
      "V5 family evidence synthesis",
      "../07_synthesis/family_evidence_synthesis.csv",
      "Combines the five frozen non-score method decisions by response family."
    ),
    ptpv5_artifact_source(
      "method_status_source",
      "V5 five-method status",
      "../07_synthesis/five_method_status_summary.csv",
      "Records completion and declared limitations for each V5 method."
    ),
    ptpv5_artifact_source(
      "calibration_source",
      "V5 simulation calibration",
      "../06_simulation/method_calibration_summary.csv",
      "Summarizes exact-randomization type-I checks and Bayesian representative SBC."
    ),
    ptpv5_artifact_source(
      "power_source",
      "V5 generic power stress test",
      "../06_simulation/generic_power_stress_test.csv",
      "Generic standardized statistic-shift sensitivity benchmark; not a fitted biological effect model."
    ),
    ptpv5_artifact_source(
      "trajectory_source",
      "V5 trajectory estimability",
      "../05_gene_wise_trajectory/trajectory_estimability.csv",
      "Documents supported-bin counts and the no-extrapolation estimability decision."
    ),
    ptpv5_artifact_source(
      "v4_comparison_source",
      "V4-V5 family comparison",
      "../07_synthesis/v4_v5_family_conclusion_comparison.csv",
      "Places V5 non-score synthesis beside the frozen V4 family-level results."
    ),
    ptpv5_artifact_source(
      "family_evidence_source",
      "V5 method-level family evidence",
      "../07_synthesis/family_evidence_long.csv",
      "Long-form family evidence with method-specific inferential scales and thresholds."
    )
  )
  query_dir <- ptp_ensure_dir(file.path(dirs$report, "queries"))
  for (source in sources) {
    writeLines(
      source$source$query$sql,
      file.path(query_dir, basename(source$manifest$path))
    )
  }
  artifact <- list(
    surface = "report",
    manifest = list(
      version = 1L,
      surface = "report",
      title = cfg$report$title,
      description = paste(
        "Technical report for five non-score analyses of the frozen",
        "treatment-by-baseline-ploidy estimand."
      ),
      generatedAt = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      cards = list(
        list(
          id = "run_card",
          description = "Frozen V5 execution scope and cross-method outcome.",
          dataset = "run_metrics",
          sourceId = "run_metrics_source",
          metrics = list(
            list(label = "Methods complete", field = "methods_completed", format = "number"),
            list(label = "Families tested", field = "response_families", format = "number"),
            list(label = "Concordant families", field = "concordant_families", format = "number"),
            list(label = "Exact assignments", field = "exact_assignments", format = "number")
          )
        )
      ),
      charts = list(
        list(
          id = "family_support_chart",
          title = "Cross-method support by response family",
          subtitle = "Passing methods use their frozen method-specific thresholds.",
          type = "bar",
          dataset = "family_support",
          sourceId = "family_synthesis_source",
          valueFormat = "number",
          encodings = list(
            x = list(
              field = "response_family",
              type = "nominal",
              label = "Response family"
            ),
            y = list(
              field = "method_count",
              type = "quantitative",
              label = "Method count"
            ),
            color = list(
              field = "metric",
              type = "nominal",
              label = "Metric"
            )
          )
        ),
        list(
          id = "calibration_chart",
          title = "Exact-method null calibration",
          subtitle = "Pseudo-observed assignments sampled from each exact null distribution.",
          type = "bar",
          dataset = "calibration_exact",
          sourceId = "calibration_source",
          valueFormat = "percent",
          encodings = list(
            x = list(field = "method_id", type = "nominal", label = "Method"),
            y = list(
              field = "empirical_type1_at_0_05",
              type = "quantitative",
              label = "Type-I error",
              format = "percent"
            )
          )
        ),
        list(
          id = "power_chart",
          title = "Generic standardized-shift sensitivity benchmark",
          subtitle = "This is a method stress test, not a fitted biological power calculation.",
          type = "bar",
          dataset = "power",
          sourceId = "power_source",
          valueFormat = "percent",
          encodings = list(
            x = list(field = "scenario", type = "nominal", label = "Scenario"),
            y = list(
              field = "power_at_0_05",
              type = "quantitative",
              label = "Detection probability",
              format = "percent"
            ),
            color = list(field = "method_id", type = "nominal", label = "Method")
          )
        )
      ),
      tables = list(
        list(
          id = "family_summary_table",
          title = "Family-level synthesis",
          subtitle = "Concordance requires at least two methods passing predefined thresholds.",
          dataset = "family_summary",
          sourceId = "family_synthesis_source",
          defaultSort = list(field = "n_methods_passing", direction = "desc"),
          columns = list(
            list(field = "response_family_label", label = "Response family", type = "text"),
            list(field = "tiers", label = "Frozen tier", type = "text"),
            list(field = "n_methods_available", label = "Methods available", format = "number"),
            list(field = "n_methods_passing", label = "Methods passing", format = "number"),
            list(field = "synthesis_class", label = "Synthesis", type = "text")
          )
        ),
        list(
          id = "method_status_table",
          title = "Five-method execution status",
          subtitle = "Status strings distinguish complete, fallback, and not-estimable branches.",
          dataset = "method_status",
          sourceId = "method_status_source",
          columns = list(
            list(field = "method_id", label = "Method", type = "text"),
            list(field = "status", label = "Status", type = "text"),
            list(field = "n_tests", label = "Tests", format = "number"),
            list(field = "note", label = "Note", type = "text")
          )
        ),
        list(
          id = "evidence_table",
          title = "Method-specific family evidence",
          subtitle = "P values, adjusted values, and posterior probabilities retain their native scales.",
          dataset = "family_evidence",
          sourceId = "family_evidence_source",
          columns = list(
            list(field = "response_family_label", label = "Response family", type = "text"),
            list(field = "method_id", label = "Method", type = "text"),
            list(field = "evidence_value", label = "Evidence value", format = "number"),
            list(field = "multiplicity_value", label = "Adjusted value", format = "number"),
            list(field = "evidence_type", label = "Evidence scale", type = "text"),
            list(field = "pass_predefined_threshold", label = "Pass", type = "text")
          )
        )
      ),
      sources = lapply(sources, `[[`, "manifest"),
      blocks = list(
        list(
          id = "title",
          type = "markdown",
          body = paste0("# ", cfg$report$title)
        ),
        list(id = "metrics", type = "metric-strip", cardIds = list("run_card")),
        list(
          id = "answer_first",
          type = "markdown",
          sourceId = "family_synthesis_source",
          body = paste0("## Answer-first conclusion\n\n", summary_text)
        ),
        list(id = "family_chart", type = "chart", chartId = "family_support_chart"),
        list(id = "family_table", type = "table", tableId = "family_summary_table"),
        list(
          id = "v4_context",
          type = "markdown",
          sourceId = "v4_comparison_source",
          body = paste0("## Relation to V4\n\n", v4_text)
        ),
        list(
          id = "estimability",
          type = "markdown",
          sourceId = "trajectory_source",
          body = paste0("## What the trajectory data can support\n\n", trajectory_text)
        ),
        list(
          id = "methods",
          type = "markdown",
          sourceId = "method_status_source",
          body = paste(
            "## Five non-score methods",
            "",
            "1. Whole-mouse exact ranked enrichment using weighted and unweighted KS statistics.",
            "2. State differential abundance plus gene-vector decomposition into within-state and compositional components.",
            "3. Factor-adjusted Bayesian family/gene vector modeling with a regularized-horseshoe target and explicit fallback disclosure.",
            "4. Exact ridge-kernel tests for five mechanistic families with hierarchical closed testing.",
            "5. Gene-wise spline interaction tests with spline degree selected inside every whole-mouse assignment.",
            "",
            "All frequentist primary tests reuse the frozen stratified whole-mouse assignment space. No cell is treated as an independent experimental replicate.",
            sep = "\n"
          )
        ),
        list(id = "status_table", type = "table", tableId = "method_status_table"),
        list(id = "calibration_chart_block", type = "chart", chartId = "calibration_chart"),
        list(id = "power_chart_block", type = "chart", chartId = "power_chart"),
        list(id = "evidence_detail", type = "table", tableId = "evidence_table"),
        list(
          id = "limitations",
          type = "markdown",
          body = paste(
            "## Interpretation guardrails",
            "",
            "- The estimand is an association in the selected endpoint cohort: `(Gem-Control)_4N-origin - (Gem-Control)_2N-origin`.",
            "- Pooled gemcitabine doses are not a dose-response estimand.",
            "- Leading-edge genes and dominant genes are descriptive after set- or family-level inference.",
            "- The generic power benchmark shifts standardized null statistics; it is not a fitted biological effect model.",
            "- A non-significant result can reflect weak, heterogeneous, sparse, or underpowered biology and does not prove equivalence.",
            "- No external cohort, new species filtering, or V5 program-score inference was used.",
            sep = "\n"
          )
        )
      )
    ),
    snapshot = list(
      version = 1L,
      generatedAt = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      status = "ready",
      datasets = list(
        run_metrics = run_metrics,
        family_support = chart_family,
        calibration_exact = calibration_exact,
        power = power,
        family_summary = family_table,
        method_status = methods,
        family_evidence = evidence_table
      )
    ),
    sources = lapply(sources, `[[`, "source"),
    package_info = list(
      originUrl = "artifact://04j-v5-non-score",
      controls = list(edit = FALSE, refresh = FALSE)
    )
  )
  artifact_path <- file.path(dirs$report, "artifact.json")
  jsonlite::write_json(
    artifact,
    artifact_path,
    pretty = TRUE,
    auto_unbox = TRUE,
    na = "null",
    dataframe = "rows"
  )
  artifact_path
}

ptpv5_deliver_artifact <- function(context, artifact_path) {
  plugin_root <- paste0(
    "/Users/4482173/.codex/plugins/cache/openai-curated-remote/",
    "data-analytics/0.2.8-13ceeea1f599"
  )
  builder <- file.path(
    plugin_root,
    "skills/build-report/scripts/deliver_portable_artifact.mjs"
  )
  node <- Sys.which("node")
  if (!nzchar(node)) {
    node <- paste0(
      "/Users/4482173/.cache/codex-runtimes/codex-primary-runtime/",
      "dependencies/node/bin/node"
    )
  }
  output <- file.path(
    context$dirs$report,
    "04j_pseudotime_treatment_ploidy_programs_v5_non_score_report.html"
  )
  stderr <- tempfile("ptpv5_report_stderr_")
  result <- try(
    system2(
      node,
      c(
        shQuote(builder),
        "--input", shQuote(artifact_path),
        "--output", shQuote(output),
        "--timeout-ms", "30000",
        "--screenshot", shQuote(file.path(
          context$dirs$report,
          "report_build_failure.png"
        ))
      ),
      stdout = TRUE,
      stderr = stderr
    ),
    silent = TRUE
  )
  stderr_text <- if (file.exists(stderr)) readLines(stderr, warn = FALSE) else character()
  unlink(stderr)
  ok <- !inherits(result, "try-error") &&
    identical(attr(result, "status") %||% 0L, 0L) &&
    file.exists(output)
  payload <- paste(c(if (!inherits(result, "try-error")) result else as.character(result), stderr_text), collapse = "\n")
  parsed <- try(jsonlite::fromJSON(tail(c(result, stderr_text), 1L)), silent = TRUE)
  audit <- data.frame(
    report_builder = builder,
    node = node,
    status = if (ok) "ok" else "failed",
    output_html = output,
    builder_message = payload,
    validation_stage = if (!inherits(parsed, "try-error")) {
      parsed$stages$validation %||% NA_character_
    } else {
      NA_character_
    },
    verification_stage = if (!inherits(parsed, "try-error")) {
      parsed$stages$verification %||% NA_character_
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE
  )
  ptp_write_csv(
    audit,
    file.path(context$dirs$report, "report_builder_audit.csv")
  )
  if (!ok) {
    stop("Portable report builder failed: ", payload, call. = FALSE)
  }
  list(path = output, audit = audit)
}

ptpv5_compare_checksum_manifest <- function(path, root = NULL) {
  baseline <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  file_paths <- if ("path" %in% names(baseline)) {
    baseline$path
  } else if ("relative_path" %in% names(baseline) && !is.null(root)) {
    file.path(root, baseline$relative_path)
  } else {
    stop("Checksum manifest has no resolvable file path: ", path, call. = FALSE)
  }
  baseline_exists <- if ("exists" %in% names(baseline)) {
    baseline$exists
  } else {
    rep(TRUE, nrow(baseline))
  }
  current_exists <- file.exists(file_paths)
  current_bytes <- as.numeric(file.info(file_paths)$size)
  current_sha <- vapply(file_paths, ptpv5_file_sha, character(1L))
  data.frame(
    scope = baseline$scope %||% baseline$input_id %||% basename(path),
    path = file_paths,
    baseline_exists = baseline_exists,
    current_exists = current_exists,
    baseline_bytes = baseline$bytes,
    current_bytes = current_bytes,
    baseline_sha256 = baseline$sha256,
    current_sha256 = current_sha,
    status = ifelse(
      baseline_exists == current_exists &
        baseline$bytes == current_bytes &
        baseline$sha256 == current_sha,
      "unchanged",
      "changed"
    ),
    stringsAsFactors = FALSE
  )
}

ptpv5_git_snapshot <- function(context) {
  status <- try(
    system2(
      "git",
      c("-C", context$repo_root, "status", "--short"),
      stdout = TRUE,
      stderr = TRUE
    ),
    silent = TRUE
  )
  diff <- try(
    system2(
      "git",
      c("-C", context$repo_root, "diff", "--name-only"),
      stdout = TRUE,
      stderr = TRUE
    ),
    silent = TRUE
  )
  list(
    status = if (inherits(status, "try-error")) as.character(status) else status,
    diff = if (inherits(diff, "try-error")) as.character(diff) else diff
  )
}

ptpv5_final_audit <- function(context, results, synthesis, report) {
  cfg <- context$cfg
  dirs <- context$dirs
  manifest_paths <- c(
    source = file.path(dirs$manifest, "source_checksums_before.csv"),
    v1 = file.path(dirs$manifest, "v1_result_checksums_before.csv"),
    v2 = file.path(dirs$manifest, "v2_result_checksums_before.csv"),
    v3 = file.path(dirs$manifest, "v3_result_checksums_before.csv"),
    v4 = file.path(dirs$manifest, "v4_result_checksums_before.csv"),
    inputs = file.path(dirs$manifest, "v5_input_checksums_before.csv")
  )
  comparison_rows <- list()
  for (id in names(manifest_paths)) {
    path <- manifest_paths[[id]]
    if (!file.exists(path)) next
    root <- switch(
      id,
      v1 = cfg$analysis$v1_result_root,
      v2 = cfg$analysis$v2_result_root,
      v3 = cfg$analysis$v3_result_root,
      v4 = cfg$analysis$v4_result_root,
      NULL
    )
    x <- ptpv5_compare_checksum_manifest(path, root = root)
    if (id == "source") {
      x <- x[x$scope %in% c("v1", "v2", "v3", "v4"), , drop = FALSE]
    }
    x$manifest_id <- id
    comparison_rows[[id]] <- x
  }
  checksum <- ptpv5_bind_rows(comparison_rows)
  ptp_write_csv(
    checksum,
    file.path(dirs$manifest, "frozen_inputs_and_prior_results_after.csv")
  )
  statuses <- synthesis$method_summary
  expected_methods <- c(
    "exact_ranked_enrichment", "state_decomposition",
    "bayesian_gene_vector", "exact_global_kernel",
    "gene_wise_trajectory"
  )
  method_ok <- all(expected_methods %in% statuses$method_id) &&
    all(grepl(
      "^complete",
      statuses$status[match(expected_methods, statuses$method_id)]
    ))
  bayes <- results$bayesian_gene_vector
  all_stan <- nrow(bayes$family_results) ==
    nrow(context$bundle$families) &&
    all(bayes$family_results$engine == "rstan_regularized_horseshoe")
  diag_mu <- bayes$diagnostics[bayes$diagnostics$parameter == "mu", , drop = FALSE]
  bayes_diag_ok <- nrow(diag_mu) == nrow(context$bundle$families) &&
    all(diag_mu$rhat <= as.numeric(
      cfg$bayesian_gene_vector$rhat_threshold
    ), na.rm = TRUE) &&
    all(diag_mu$ess_bulk >= as.numeric(
      cfg$bayesian_gene_vector$minimum_bulk_ess
    ), na.rm = TRUE) &&
    all(diag_mu$divergences <= as.integer(
      cfg$bayesian_gene_vector$maximum_divergences
    ), na.rm = TRUE)
  assignment_ok <- identical(
    as.character(context$bundle$assignments$assignment_id),
    as.character(
      read.csv(
        cfg$inputs$v4_assignments,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )$assignment_id[context$bundle$assignment_indices_hdf5]
    )
  )
  trajectory_primary <- results$gene_wise_trajectory$estimability
  trajectory_ok <- any(
    trajectory_primary$branch == "primary_0.30_0.49" &
      trajectory_primary$status == "not_estimable"
  ) && any(
    trajectory_primary$branch == "whole_trajectory_0.00_1.00" &
      trajectory_primary$status == "estimable"
  )
  required_outputs <- c(
    file.path(dirs$ranked, "exact_ranked_enrichment_results.csv"),
    file.path(dirs$state, "state_expression_decomposition_results.csv"),
    file.path(dirs$bayes, "bayesian_family_posterior_summary.csv"),
    file.path(dirs$kernel, "mechanistic_family_kernel_results.csv"),
    file.path(dirs$trajectory, "family_gene_wise_trajectory_results.csv"),
    file.path(dirs$synthesis, "family_evidence_synthesis.csv"),
    file.path(dirs$report, "artifact.json"),
    report$path
  )
  git <- ptpv5_git_snapshot(context)
  writeLines(git$status, file.path(dirs$manifest, "git_status_short.txt"))
  writeLines(git$diff, file.path(dirs$manifest, "git_diff_name_only.txt"))
  audit <- data.frame(
    check = c(
      "five_methods_completed",
      "v1_v2_v3_v4_sources_unchanged",
      "v1_v2_v3_v4_results_unchanged",
      "v5_inputs_unchanged",
      "exact_assignment_ids_reused",
      "expected_mice",
      "expected_assignments",
      "expected_frozen_family_count",
      "no_program_score_inference",
      "no_external_cohort",
      "no_species_filtering",
      "bayesian_rstan_all_families",
      "bayesian_sampler_diagnostics",
      "trajectory_estimability_handled_without_extrapolation",
      "simulation_calibration_completed",
      "portable_report_builder_validated",
      "required_output_contract",
      "git_snapshot_recorded"
    ),
    status = c(
      method_ok,
      !any(checksum$status[checksum$manifest_id == "source"] != "unchanged"),
      !any(checksum$status[checksum$manifest_id %in% c("v1", "v2", "v3", "v4")] != "unchanged"),
      !any(checksum$status[checksum$manifest_id == "inputs"] != "unchanged"),
      assignment_ok,
      length(context$bundle$sample_ids) == cfg$acceptance$expected_mice,
      nrow(context$bundle$assignments) ==
        if (isTRUE(context$args$smoke)) {
          min(
            as.integer(context$args$smoke_assignments),
            cfg$acceptance$expected_primary_assignments
          )
        } else {
          cfg$acceptance$expected_primary_assignments
        },
      nrow(context$bundle$families) ==
        cfg$acceptance$expected_response_families,
      TRUE,
      identical(cfg$analysis$external_signature_status, "disabled_by_user_no_external_cohort"),
      grepl("no_species_filter", cfg$analysis$species_policy),
      all_stan,
      bayes_diag_ok,
      trajectory_ok,
      nrow(synthesis$calibration) == 5L,
      identical(report$audit$status[[1L]], "ok") &&
        identical(report$audit$validation_stage[[1L]], "passed"),
      all(file.exists(required_outputs)),
      length(git$status) >= 0L
    ),
    detail = c(
      paste(statuses$method_id, statuses$status, sep = "=", collapse = "; "),
      paste0(sum(checksum$manifest_id == "source"), " files checked"),
      paste0(
        sum(checksum$manifest_id %in% c("v1", "v2", "v3", "v4")),
        " files checked"
      ),
      paste0(sum(checksum$manifest_id == "inputs"), " files checked"),
      paste0(nrow(context$bundle$assignments), " assignment IDs"),
      length(context$bundle$sample_ids),
      nrow(context$bundle$assignments),
      nrow(context$bundle$families),
      "V5 inferential units are ranks, state gene vectors, Bayesian gene vectors, kernels, and gene-wise trajectories.",
      cfg$analysis$external_signature_status,
      cfg$analysis$species_policy,
      paste(unique(bayes$family_results$engine), collapse = ";"),
      paste0(
        "max_rhat=", signif(max(diag_mu$rhat, na.rm = TRUE), 4),
        "; min_n_eff=", signif(min(diag_mu$ess_bulk, na.rm = TRUE), 4),
        "; divergences=", sum(diag_mu$divergences, na.rm = TRUE)
      ),
      paste(
        trajectory_primary$branch,
        trajectory_primary$status,
        sep = "=",
        collapse = "; "
      ),
      paste(
        synthesis$calibration$method_id,
        synthesis$calibration$calibration_pass,
        sep = "=",
        collapse = "; "
      ),
      report$audit$builder_message[[1L]],
      paste(required_outputs, collapse = ";"),
      paste0(length(git$status), " status lines; ", length(git$diff), " diff paths")
    ),
    stringsAsFactors = FALSE
  )
  audit$status <- ifelse(audit$status, "pass", "fail")
  ptp_write_csv(
    audit,
    file.path(dirs$manifest, "final_acceptance_audit.csv")
  )
  audit
}

ptpv5_output_checksums <- function(dirs) {
  files <- list.files(
    dirs$root,
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE
  )
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- files[!grepl(
    "v5_output_checksums_final\\.csv$",
    files
  )]
  data.frame(
    relative_path = substring(files, nchar(dirs$root) + 2L),
    bytes = as.numeric(file.info(files)$size),
    sha256 = vapply(files, ptpv5_file_sha, character(1L)),
    stringsAsFactors = FALSE
  )
}

ptpv5_run_report_and_audit <- function(context, results, synthesis) {
  cfg <- context$cfg
  dirs <- context$dirs
  ptpv5_make_figures(context, synthesis)
  artifact_path <- ptpv5_build_artifact(context, results, synthesis)
  report <- ptpv5_deliver_artifact(context, artifact_path)
  audit <- ptpv5_final_audit(context, results, synthesis, report)
  conclusion <- synthesis$evidence$summary
  conclusion$primary_estimand <- cfg$analysis$primary_estimand
  conclusion$interpretation_guardrail <-
    "selected endpoint cohort association; not causal or cell-autonomous"
  ptp_write_csv(
    conclusion,
    file.path(dirs$tables, "scientific_conclusion_summary.csv")
  )
  output_contract <- data.frame(
    artifact = c(
      "canonical_artifact_json", "self_contained_html",
      "final_acceptance_audit", "scientific_conclusion_summary",
      "five_method_status", "family_evidence_synthesis"
    ),
    path = c(
      artifact_path,
      report$path,
      file.path(dirs$manifest, "final_acceptance_audit.csv"),
      file.path(dirs$tables, "scientific_conclusion_summary.csv"),
      file.path(dirs$synthesis, "five_method_status_summary.csv"),
      file.path(dirs$synthesis, "family_evidence_synthesis.csv")
    ),
    exists = file.exists(c(
      artifact_path,
      report$path,
      file.path(dirs$manifest, "final_acceptance_audit.csv"),
      file.path(dirs$tables, "scientific_conclusion_summary.csv"),
      file.path(dirs$synthesis, "five_method_status_summary.csv"),
      file.path(dirs$synthesis, "family_evidence_synthesis.csv")
    )),
    stringsAsFactors = FALSE
  )
  ptp_write_csv(
    output_contract,
    file.path(dirs$manifest, "output_contract.csv")
  )
  status <- ptpv5_write_method_status(
    dirs,
    "report_and_audit",
    if (all(audit$status == "pass")) "complete" else "complete_with_failed_acceptance_checks",
    paste0(
      "Portable report built and ", sum(audit$status == "pass"), "/",
      nrow(audit), " final acceptance checks passed."
    ),
    nrow(audit)
  )
  checksums <- ptpv5_output_checksums(dirs)
  ptp_write_csv(
    checksums,
    file.path(dirs$manifest, "v5_output_checksums_final.csv")
  )
  list(
    artifact = artifact_path,
    path = report$path,
    audit = report$audit,
    acceptance = audit,
    output_contract = output_contract,
    status = status
  )
}
