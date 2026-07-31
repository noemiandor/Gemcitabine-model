#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else "Code/in-vivo/figure7/run_figure7.R"
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE)
for (file in c(
  "common_io.R", "feature_species_policy.R", "input_preflight.R",
  "seurat_upstream_selection.R",
  "tgi_data.R", "tgi_statistics.R",
  "tgi_panels.R", "context_panels.R", "state_pathway_panel.R",
  "generated_state_pathway_reference.R", "state_pathway_analysis.R"
)) {
  sys.source(file.path(script_dir, "src", file), envir = .GlobalEnv)
}

args <- figure7_parse_args(commandArgs(trailingOnly = TRUE))
mode <- figure7_arg(args, "mode", "standard")
allowed_modes <- c("standard", "full-analysis", "full-workflow", "render-only")
if (!mode %in% allowed_modes) figure7_stop("Unknown Figure 7 mode: ", mode)
panel_set <- figure7_arg(args, "panel-set", "a-f")
if (!panel_set %in% c("a-f", "a-e")) {
  figure7_stop("Unknown Figure 7 panel set: ", panel_set)
}
include_state_pathway <- identical(panel_set, "a-f")
required_packages <- c("ggplot2", "ggrepel", "yaml")
if (include_state_pathway) {
  required_packages <- c(required_packages, "patchwork", "pheatmap")
}
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
if (identical(mode, "full-analysis") && !include_state_pathway) {
  figure7_stop("Full analysis is incompatible with --panel-set=a-e")
}
panel_ids <- figure7_panel_ids(include_state_pathway)
config_path <- normalizePath(figure7_arg(args, "config", file.path(script_dir, "figure7_config.yaml")), mustWork = FALSE)
tgi_day_arg <- figure7_arg(args, "tgi-day", "")
config <- figure7_read_config(config_path, if (nzchar(tgi_day_arg)) tgi_day_arg else NULL)
output_dir <- normalizePath(figure7_arg(args, "output-dir", required = TRUE), mustWork = FALSE)
si_cache_dir <- ""
si_cache_policy <- "not_applicable"
si_cache_upstream_input_manifest <- ""
si_cache_source_run_config <- ""
si_cache_source_provenance <- ""
if (include_state_pathway) {
  normalized_composition_path <- file.path(
    repo_root,
    "Code",
    "in-vivo",
    "SI_figures",
    "normalized_composition.R"
  )
  shared_context_path <- file.path(
    repo_root,
    "Code",
    "in-vivo",
    "SI_figures",
    "shared_context_panels.R"
  )
  if (!file.exists(normalized_composition_path) ||
      !file.exists(shared_context_path)) {
    figure7_stop("Main Figure 7 shared SI panel helpers are unavailable")
  }
  sys.source(normalized_composition_path, envir = .GlobalEnv)
  sys.source(shared_context_path, envir = .GlobalEnv)
  si_cache_policy <- figure7_arg(args, "si-cache-policy", "reviewed")
  if (!si_cache_policy %in% c("reviewed", "generated-human-only")) {
    figure7_stop("Unknown --si-cache-policy: ", si_cache_policy)
  }
  if (identical(mode, "full-workflow") &&
      !identical(si_cache_policy, "generated-human-only")) {
    figure7_stop(
      "Full-workflow Figure 7 requires a generated-human-only SI context cache"
    )
  }
  if (!mode %in% c("full-workflow", "render-only") &&
      !identical(si_cache_policy, "reviewed")) {
    figure7_stop(
      "Standard/full-analysis Figure 7 requires the reviewed SI context cache"
    )
  }
  si_cache_dir <- normalizePath(
    figure7_arg(
      args,
      "si-table-cache-dir",
      as.character(config$si_figures$cache_root)
    ),
    mustWork = TRUE
  )
  if (identical(si_cache_policy, "generated-human-only")) {
    si_cache_upstream_input_manifest <- normalizePath(
      figure7_arg(args, "si-cache-upstream-input-manifest", required = TRUE),
      mustWork = TRUE
    )
    si_cache_source_run_config <- normalizePath(
      figure7_arg(args, "si-cache-source-run-config", required = TRUE),
      mustWork = TRUE
    )
    si_cache_source_provenance <- normalizePath(
      figure7_arg(args, "si-cache-source-provenance", required = TRUE),
      mustWork = TRUE
    )
  }
  invisible(figure7_validate_si_cache(
    si_cache_dir,
    repo_root,
    config,
    si_cache_policy,
    si_cache_upstream_input_manifest,
    si_cache_source_run_config,
    si_cache_source_provenance
  ))
}
state_pathway_results_root_arg <- figure7_arg(args, "state-pathway-results-root", "")
state_pathway_results_root <- if (nzchar(state_pathway_results_root_arg) && !identical(mode, "full-workflow")) {
  normalizePath(state_pathway_results_root_arg, mustWork = TRUE)
} else {
  ""
}

figure7_metadata_locator <- function(path) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    return("not_recorded")
  }
  if (startsWith(path, "external:") || identical(path, "not_recorded")) {
    return(path)
  }
  normalized <- normalizePath(path, mustWork = FALSE)
  prefix <- paste0(repo_root, .Platform$file.sep)
  if (startsWith(normalized, prefix)) {
    substring(normalized, nchar(prefix) + 1L)
  } else {
    paste0("external:", basename(normalized))
  }
}

write_metadata <- function(
  output_dir,
  mode,
  config,
  config_path,
  panel_ids,
  state_pathway_results_root = "",
  workflow = NULL,
  reference_identity = NULL,
  context_cache = NULL,
  cellcycle_input = "",
  cellcycle_sha256 = "",
  noncellcycle_input = "",
  noncellcycle_sha256 = ""
) {
  tgi_day <- figure7_tgi_day(config)
  tgi_measure <- figure7_tgi_measure(config)
  if (is.null(reference_identity)) {
    reference_identity <- if (!is.null(workflow)) {
      list(
        id = workflow$state_pathway_reference_id,
        kind = workflow$state_pathway_reference_kind,
        canonical_publication_allowed =
          workflow$canonical_publication_allowed
      )
    } else if (include_state_pathway) {
      list(
        id = as.character(config$state_pathways$reviewed_reference_id),
        kind = as.character(config$state_pathways$reviewed_reference_kind),
        canonical_publication_allowed = isTRUE(
          config$state_pathways$
            reviewed_reference_canonical_publication_allowed
        )
      )
    } else {
      list(
        id = "not_applicable",
        kind = "not_applicable",
        canonical_publication_allowed = FALSE
      )
    }
  }
  overall_canonical <- isTRUE(
    reference_identity$canonical_publication_allowed
  ) && (
    !include_state_pathway ||
      (!is.null(context_cache) &&
        isTRUE(context_cache$cache_canonical_publication_allowed))
  )
  run_config <- data.frame(
    key = c("module", "mode", "panel_set", "tgi_outcome", "tgi_day", "tgi_measure",
            "matched_control_summary", "matched_control_group", "etp_method", "etp_threshold",
            "state_pathway_reference_id", "state_pathway_reference_kind",
            "canonical_publication_allowed",
            "state_interval_start", "state_interval_end",
            "state_pathway_source_results_root", "config_sha256",
            "cellcycle_input", "cellcycle_sha256",
            "noncellcycle_input", "noncellcycle_sha256"),
    value = c("in_vivo_figure7", mode, panel_set, "day", as.character(tgi_day), tgi_measure, "mean", "initial_ploidy",
              config$etp$method, as.character(config$etp$threshold),
              reference_identity$id,
              reference_identity$kind,
              tolower(as.character(overall_canonical)),
              as.character(config$state_pathways$accumulated_interval$start),
              as.character(config$state_pathways$accumulated_interval$end),
              if (nzchar(state_pathway_results_root)) state_pathway_results_root else "not_recorded",
              figure7_sha256(config_path),
              figure7_metadata_locator(cellcycle_input),
              if (nzchar(cellcycle_sha256)) cellcycle_sha256 else "not_recorded",
              figure7_metadata_locator(noncellcycle_input),
              if (nzchar(noncellcycle_sha256)) noncellcycle_sha256 else "not_recorded"),
    stringsAsFactors = FALSE
  )
  if (include_state_pathway) {
    if (is.null(context_cache)) {
      figure7_stop("Full Figure 7 metadata requires an SI context cache")
    }
    mapping <- unlist(
      config$panels$main_composite$panel_order,
      use.names = TRUE
    )
    run_config <- rbind(
      run_config,
      data.frame(
        key = c(
          "main_composite_panel_set",
          "main_composite_filename",
          "main_composite_panel_order",
          "si_context_cache_policy",
          "si_context_cache_kind",
          "si_context_cache_manifest",
          "si_context_cache_manifest_sha256",
          "si_context_cache_canonical_publication_allowed",
          "si_context_upstream_input_manifest",
          "si_context_upstream_input_manifest_sha256",
          "si_context_source_run_config",
          "si_context_source_run_config_sha256",
          "si_context_source_provenance",
          "si_context_source_provenance_sha256"
        ),
        value = c(
          "a-k",
          context_cache$composite_filename,
          paste(paste(names(mapping), mapping, sep = "="), collapse = ";"),
          context_cache$cache_policy,
          context_cache$cache_kind,
          figure7_metadata_locator(context_cache$cache_manifest_path),
          context_cache$cache_manifest_sha256,
          tolower(as.character(
            context_cache$cache_canonical_publication_allowed
          )),
          figure7_metadata_locator(context_cache$upstream_input_manifest),
          context_cache$upstream_input_manifest_sha256,
          figure7_metadata_locator(context_cache$source_run_config),
          context_cache$source_run_config_sha256,
          figure7_metadata_locator(context_cache$source_provenance),
          context_cache$source_provenance_sha256
        ),
        stringsAsFactors = FALSE
      )
    )
    if (identical(context_cache$cache_policy, "reviewed")) {
      run_config <- rbind(
        run_config,
        data.frame(
          key = c(
            "reviewed_si_cache_manifest",
            "reviewed_si_cache_manifest_sha256"
          ),
          value = c(
            figure7_metadata_locator(context_cache$cache_manifest_path),
            context_cache$cache_manifest_sha256
          ),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  if (!is.null(workflow)) {
    workflow_loom_files <- if (dir.exists(workflow$loom_root)) {
      list.files(workflow$loom_root, pattern = "[.]loom$", recursive = TRUE, full.names = TRUE)
    } else {
      character()
    }
    workflow_loom_bytes <- if (length(workflow_loom_files)) sum(file.info(workflow_loom_files)$size) else 0
    run_config <- rbind(
      run_config,
      data.frame(
        key = c(
          "workflow_initial_state", "workflow_executed_stages", "workflow_scvelo_metrics",
          "workflow_scvelo_sha256", "workflow_cellcycle_input", "workflow_cellcycle_sha256",
          "workflow_noncellcycle_input", "workflow_noncellcycle_sha256",
          "raw_data_source", "raw_data_doi", "zenodo_record_id", "raw_data_dir",
          "raw_data_download_roles", "raw_data_validation_roles", "raw_download_workers",
          "raw_download_connections_per_file", "analysis_jobs",
          "raw_data_manifest_sha256", "seurat_rds_sha256",
          "workflow_seurat_rds", "loom_file_count", "loom_total_bytes",
          "workflow_loom_root",
          "workflow_state_pathway_results", "workflow_state_pathway_reference",
          "workflow_log_dir"
        ),
        value = c(
          workflow$initial_state,
          if (length(workflow$executed_stages)) paste(workflow$executed_stages, collapse = ",") else "none",
          workflow$scvelo_metrics,
          if (file.exists(workflow$scvelo_metrics)) figure7_sha256(workflow$scvelo_metrics) else "not_available",
          workflow$cellcycle,
          figure7_sha256(workflow$cellcycle),
          workflow$noncellcycle,
          figure7_sha256(workflow$noncellcycle),
          workflow$raw_data_status,
          as.character(config$raw_data$doi),
          as.character(config$raw_data$record_id),
          workflow$raw_data_dir,
          if (length(workflow$raw_download_roles)) paste(workflow$raw_download_roles, collapse = ",") else "none",
          if (length(workflow$raw_validation_roles)) paste(workflow$raw_validation_roles, collapse = ",") else "none",
          as.character(workflow$download_workers),
          as.character(workflow$download_connections_per_file),
          as.character(workflow$analysis_jobs),
          figure7_sha256(workflow$raw_manifest),
          if (nzchar(workflow$seurat_rds_sha256)) workflow$seurat_rds_sha256 else "not_available",
          workflow$seurat_rds,
          as.character(length(workflow_loom_files)),
          as.character(workflow_loom_bytes),
          workflow$loom_root,
          if (nzchar(workflow$state_pathway_root)) {
            workflow$state_pathway_root
          } else {
            "not_applicable"
          },
          if (nzchar(workflow$saved_reference)) {
            workflow$saved_reference
          } else {
            "not_applicable"
          },
          workflow$log_dir
        ),
        stringsAsFactors = FALSE
      )
    )
    reconstruction_manifest <- figure7_workflow_scalar(
      workflow,
      "seurat_reconstruction_manifest"
    )
    final_stage_manifest <- figure7_workflow_scalar(
      workflow,
      "seurat_final_stage_manifest"
    )
    lineage_config <- data.frame(
      key = c(
        "workflow_seurat_source",
        "workflow_seurat_upstream_dir",
        "workflow_seurat_reconstruction_manifest",
        "workflow_seurat_reconstruction_manifest_sha256",
        "workflow_seurat_final_stage_manifest",
        "workflow_seurat_final_stage_manifest_sha256"
      ),
      value = c(
        figure7_workflow_scalar(
          workflow,
          "seurat_source",
          "not_recorded"
        ),
        if (nzchar(reconstruction_manifest)) {
          dirname(reconstruction_manifest)
        } else {
          "not_applicable"
        },
        if (nzchar(reconstruction_manifest)) {
          reconstruction_manifest
        } else {
          "not_applicable"
        },
        figure7_workflow_scalar(
          workflow,
          "seurat_reconstruction_manifest_sha256",
          "not_available"
        ),
        if (nzchar(final_stage_manifest)) {
          final_stage_manifest
        } else {
          "not_applicable"
        },
        figure7_workflow_scalar(
          workflow,
          "seurat_final_stage_manifest_sha256",
          "not_available"
        )
      ),
      stringsAsFactors = FALSE
    )
    named_lineage_rows <- function(values, prefix) {
      if (is.null(values) || !length(values)) {
        return(data.frame(
          key = character(),
          value = character(),
          stringsAsFactors = FALSE
        ))
      }
      if (is.null(names(values)) ||
          anyNA(names(values)) ||
          any(!nzchar(names(values))) ||
          anyDuplicated(names(values))) {
        figure7_stop(
          "Workflow Seurat lineage has invalid dependency role names"
        )
      }
      data.frame(
        key = paste0(prefix, names(values)),
        value = as.character(values),
        stringsAsFactors = FALSE
      )
    }
    lineage_config <- rbind(
      lineage_config,
      named_lineage_rows(
        workflow$seurat_transitive_dependencies,
        "workflow_seurat_transitive_dependency:"
      ),
      named_lineage_rows(
        workflow$seurat_scientific_code_contracts,
        "workflow_seurat_"
      ),
      named_lineage_rows(
        workflow$seurat_upstream_common_contract,
        "workflow_seurat_upstream_common_contract:"
      )
    )
    run_config <- rbind(run_config, lineage_config)
  }
  contract <- figure7_panel_contract(
    config,
    panel_ids,
    if (include_state_pathway) {
      context_cache$composite_filename
    } else {
      figure7_main_composite_filename(config)
    }
  )
  contract$tgi_outcome <- "day"
  contract$tgi_day <- tgi_day
  contract$tgi_measure <- tgi_measure
  contract$matched_control_summary <- "mean"
  contract$matched_control_group <- "initial_ploidy"
  figure7_write_tsv(run_config, file.path(output_dir, "metadata", "run_config.tsv"))
  figure7_write_tsv(contract, file.path(output_dir, "metadata", "panel_contract.tsv"))
  session_info <- sub("[[:space:]]+$", "", utils::capture.output(sessionInfo()))
  writeLines(session_info, file.path(output_dir, "metadata", "session_info.txt"), useBytes = TRUE)
}

render_from_run <- function(
  source_dir,
  output_dir,
  config,
  config_path,
  panel_ids,
  include_state_pathway
) {
  tgi_day <- figure7_tgi_day(config)
  tgi_measure <- figure7_tgi_measure(config)
  source_dir <- normalizePath(source_dir, mustWork = TRUE)
  if (identical(source_dir, normalizePath(output_dir, mustWork = FALSE))) {
    figure7_stop("render-only output must differ from the immutable source run")
  }
  source_run_config <- figure7_read_tsv(file.path(source_dir, "metadata", "run_config.tsv"), c("key", "value"))
  if (anyDuplicated(source_run_config$key)) {
    figure7_stop("render-only source run_config.tsv contains duplicate keys")
  }
  source_value <- stats::setNames(
    as.character(source_run_config$value),
    source_run_config$key
  )
  source_scalar <- function(key, default = "") {
    value <- source_value[[key]]
    if (is.null(value) || length(value) != 1L || is.na(value)) {
      default
    } else {
      as.character(value)
    }
  }
  context_cache <- if (include_state_pathway) {
    figure7_build_context_panels(
      si_cache_dir,
      repo_root,
      config,
      policy = si_cache_policy,
      upstream_input_manifest = si_cache_upstream_input_manifest,
      source_run_config = si_cache_source_run_config,
      source_provenance = si_cache_source_provenance
    )
  } else {
    NULL
  }
  if (include_state_pathway &&
      (!identical(source_scalar("si_context_cache_policy"), si_cache_policy) ||
        !identical(
          source_scalar("si_context_cache_manifest_sha256"),
          context_cache$cache_manifest_sha256
        ))) {
    figure7_stop(
      "render-only source run does not bind the selected SI context cache"
    )
  }
  source_contract <- figure7_read_tsv(file.path(source_dir, "metadata", "panel_contract.tsv"),
                                      c("panel_id", "filename", "tgi_outcome", "tgi_day", "tgi_measure",
                                        "matched_control_summary", "matched_control_group"))
  config_hash <- source_value[["config_sha256"]]
  if (is.na(config_hash) || !identical(config_hash, figure7_sha256(config_path))) {
    figure7_stop("render-only source run was produced with a different Figure 7 config")
  }
  expected_contract <- figure7_panel_contract(
    config,
    panel_ids,
    if (include_state_pathway) {
      context_cache$composite_filename
    } else {
      figure7_main_composite_filename(config)
    }
  )
  observed_contract <- source_contract[, c("panel_id", "filename")]
  rownames(observed_contract) <- NULL
  if (!identical(observed_contract, expected_contract)) {
    figure7_stop("render-only source panel contract disagrees with the requested panel-set contract")
  }
  if (any(source_contract$tgi_outcome != "day") || any(figure7_numeric(source_contract$tgi_day) != tgi_day) ||
      any(source_contract$tgi_measure != tgi_measure) ||
      any(source_contract$matched_control_summary != "mean") ||
      any(source_contract$matched_control_group != "initial_ploidy")) {
    figure7_stop("render-only source panel contract has incompatible TGI metadata")
  }
  table_path <- function(name) file.path(source_dir, "tables", name)
  a <- figure7_read_tsv(table_path("panel_7A_plot_data.tsv"),
    c("sample_id", "initial_ploidy", "dose", "day", "tumor_volume_change", "series_type", "is_highlight_day"))
  a$is_highlight_day <- as.character(a$is_highlight_day) %in% c("TRUE", "T", "1")
  b <- figure7_read_tsv(table_path("panel_7B_plot_data.tsv"),
    c("comparison_id", "panel", "pseudotime", "mean_ecdf", "color_group", "curve_label", "line_group"))
  bt <- figure7_read_tsv(table_path("panel_7B_tests.tsv"), c("comparison_id", "panel", "annotation"))
  cdata <- figure7_read_tsv(table_path("panel_7C_plot_data.tsv"), c("sample_id", "initial_ploidy", "dose", tgi_measure))
  ct <- figure7_read_tsv(table_path("panel_7C_test.tsv"),
    c("dose_adjusted_difference_high_minus_low", "permutation_p_two_sided", "n_group_low", "n_group_high"))
  d <- figure7_read_tsv(table_path("panel_7D_plot_data.tsv"), c("sample_id", "shift_centered", "tgi_centered", "dose", "etp_group"))
  dt <- figure7_read_tsv(table_path("panel_7D_test.tsv"), c("n", "estimate", "permutation_p_two_sided"))
  e <- figure7_read_tsv(table_path("panel_7E_plot_data.tsv"), c("sample_id", "sample_mean_endpoint_ploidy", tgi_measure, "dose", "etp_group"))
  et <- figure7_read_tsv(table_path("panel_7E_test.tsv"), c("n", "estimate", "permutation_p_two_sided"))
  for (tab in list(a, b, bt, cdata, ct, d, dt, e, et)) {
    required_metadata <- c("tgi_outcome", "tgi_day", "tgi_measure", "matched_control_summary", "matched_control_group")
    missing_metadata <- setdiff(required_metadata, names(tab))
    if (length(missing_metadata) || any(tab$tgi_outcome != "day") ||
        any(figure7_numeric(tab$tgi_day) != tgi_day) ||
        any(tab$tgi_measure != tgi_measure) || any(tab$matched_control_summary != "mean") ||
        any(tab$matched_control_group != "initial_ploidy")) {
      figure7_stop("render-only source A-E tables have missing or incompatible frozen TGI metadata")
    }
  }
  f <- NULL
  provenance <- NULL
  reference_identity <- list(
    id = "not_applicable",
    kind = "not_applicable",
    canonical_publication_allowed = FALSE
  )
  if (include_state_pathway) {
    f <- figure7_read_tsv(table_path("panel_7F_pathway_activity_plot_data.tsv"),
      c("collection_id", "collection_label", "collection_display_order", "pathway_id", "pathway_label",
        "pathway_display_order", "selected_direction", "selected_rank_within_direction",
        "pseudotime", "standardized_activity"))
    pathway_meta <- unique(f[, c(
      "collection_id", "collection_display_order",
      "pathway_id", "pathway_display_order"
    )])
    pathway_keys <- paste(
      pathway_meta$collection_id,
      pathway_meta$pathway_id,
      sep = "\r"
    )
    collections <- as.character(unlist(config$state_pathways$collections))
    collection_counts <- table(factor(
      pathway_meta$collection_id,
      levels = collections
    ))
    max_per_collection <-
      as.integer(config$state_pathways$top_positive_per_collection) +
      as.integer(config$state_pathways$top_negative_per_collection)
    if (!nrow(pathway_meta) || anyDuplicated(pathway_keys) ||
        !setequal(unique(pathway_meta$collection_id), collections) ||
        any(collection_counts < 1L) ||
        any(collection_counts > max_per_collection)) {
      figure7_stop("render-only panel-7F collection contract is invalid")
    }
    f_key <- interaction(f$collection_id, f$pathway_id, drop = TRUE)
    grids <- lapply(levels(f_key), function(key) sort(figure7_numeric(f$pseudotime[f_key == key])))
    if (any(lengths(grids) != as.integer(config$state_pathways$grid_size)) ||
        !all(vapply(grids, identical, logical(1L), grids[[1L]]))) {
      figure7_stop("render-only panel-7F pathways must share the identical frozen 501-point grid")
    }
    provenance <- file.path(source_dir, "metadata", "state_pathway_provenance.tsv")
    if (!file.exists(provenance)) figure7_stop("render-only source is missing state_pathway_provenance.tsv")
    identity_keys <- c(
      "state_pathway_reference_id",
      "state_pathway_reference_kind",
      "canonical_publication_allowed"
    )
    if (!all(identity_keys %in% names(source_value))) {
      figure7_stop(
        "render-only source lacks the panel-7F publication identity contract"
      )
    }
    reference_identity <- list(
      id = source_value[["state_pathway_reference_id"]],
      kind = source_value[["state_pathway_reference_kind"]],
      canonical_publication_allowed = identical(
        tolower(source_value[["canonical_publication_allowed"]]),
        "true"
      )
    )
    historical_identity <- identical(
      unname(unlist(reference_identity)),
      unname(unlist(list(
        id = as.character(config$state_pathways$reference_id),
        kind = as.character(config$state_pathways$reference_kind),
        canonical_publication_allowed = FALSE
      )))
    )
    generated_identity <- identical(
      unname(unlist(reference_identity)),
      unname(unlist(list(
        id = as.character(config$state_pathways$generated_reference_id),
        kind = as.character(
          config$state_pathways$generated_reference_kind
        ),
        canonical_publication_allowed = FALSE
      )))
    )
    reviewed_identity <- identical(
      unname(unlist(reference_identity)),
      unname(unlist(list(
        id = as.character(config$state_pathways$reviewed_reference_id),
        kind = as.character(
          config$state_pathways$reviewed_reference_kind
        ),
        canonical_publication_allowed = TRUE
      )))
    )
    if (!historical_identity && !reviewed_identity && !generated_identity) {
      figure7_stop(
        "render-only source has an invalid panel-7F publication identity"
      )
    }
    if (historical_identity &&
        (nrow(pathway_meta) != 24L || any(collection_counts != 8L))) {
      figure7_stop(
        "render-only historical panel-7F table must contain ",
        "eight pathways in each of three collections"
      )
    }
    provenance_table <- figure7_read_tsv(
      provenance,
      c("key", "value")
    )
    if (anyDuplicated(provenance_table$key)) {
      figure7_stop(
        "render-only source state-pathway provenance contains duplicate keys"
      )
    }
    provenance_value <- stats::setNames(
      as.character(provenance_table$value),
      provenance_table$key
    )
    provenance_matches <- if (historical_identity) {
      historical_provenance_matches <- identical(
        provenance_value[["canonical_reference_id"]],
        as.character(config$state_pathways$reference_id)
      )
      if (historical_provenance_matches) {
        for (reference_file in figure7_state_required_files()) {
          source_path <- if (identical(
            reference_file,
            "state_pathway_provenance.tsv"
          )) {
            provenance
          } else {
            table_path(reference_file)
          }
          figure7_verify_checksum(
            source_path,
            config$state_pathways$expected_files[[reference_file]],
            paste("render-only historical panel-7F", reference_file)
          )
        }
      }
      historical_provenance_matches
    } else if (reviewed_identity) {
      reviewed_provenance_matches <-
        identical(
          provenance_value[["canonical_reference_id"]],
          as.character(config$state_pathways$reviewed_reference_id)
        ) &&
        identical(
          provenance_value[["reference_kind"]],
          as.character(config$state_pathways$reviewed_reference_kind)
        ) &&
        identical(
          provenance_value[["canonical_publication_allowed"]],
          "true"
        )
      if (reviewed_provenance_matches) {
        for (reference_file in figure7_state_required_files()) {
          source_path <- if (identical(
            reference_file,
            "state_pathway_provenance.tsv"
          )) {
            provenance
          } else {
            table_path(reference_file)
          }
          figure7_verify_checksum(
            source_path,
            config$state_pathways$
              reviewed_expected_files[[reference_file]],
            paste("render-only reviewed panel-7F", reference_file)
          )
        }
      }
      reviewed_provenance_matches
    } else {
      identical(
        provenance_value[["generated_reference_id"]],
        as.character(config$state_pathways$generated_reference_id)
      ) &&
        identical(
          provenance_value[["reference_kind"]],
          as.character(config$state_pathways$generated_reference_kind)
        ) &&
        identical(
          provenance_value[["canonical_publication_allowed"]],
          "false"
        )
    }
    if (!isTRUE(provenance_matches)) {
      figure7_stop(
        "render-only panel-7F provenance contradicts run metadata"
      )
    }
    if (reviewed_identity) {
      temporary_root <- tempfile("figure7_render_reviewed_reference_")
      temporary_reference <- file.path(
        temporary_root,
        as.character(config$state_pathways$reviewed_reference_id)
      )
      dir.create(temporary_reference, recursive = TRUE)
      on.exit(unlink(temporary_root, recursive = TRUE, force = TRUE), add = TRUE)
      compact_files <- setdiff(
        figure7_state_required_files(),
        "state_pathway_provenance.tsv"
      )
      compact_sources <- c(
        file.path(source_dir, "tables", compact_files),
        provenance
      )
      copied <- file.copy(
        compact_sources,
        file.path(
          temporary_reference,
          c(compact_files, "state_pathway_provenance.tsv")
        )
      )
      if (any(!copied)) {
        figure7_stop(
          "render-only reviewed panel-7F compact reference is incomplete"
        )
      }
      figure7_validate_reviewed_state_reference(
        temporary_reference,
        config
      )
    }
    if (generated_identity) {
      temporary_root <- tempfile("figure7_render_generated_reference_")
      temporary_reference <- file.path(
        temporary_root,
        as.character(config$state_pathways$generated_reference_id)
      )
      dir.create(temporary_reference, recursive = TRUE)
      on.exit(unlink(temporary_root, recursive = TRUE, force = TRUE), add = TRUE)
      compact_files <- setdiff(
        figure7_generated_state_required_files(),
        "state_pathway_provenance.tsv"
      )
      compact_sources <- c(
        file.path(source_dir, "tables", compact_files),
        provenance
      )
      copied <- file.copy(
        compact_sources,
        file.path(
          temporary_reference,
          c(compact_files, "state_pathway_provenance.tsv")
        )
      )
      if (any(!copied)) {
        figure7_stop(
          "render-only generated panel-7F compact reference is incomplete"
        )
      }
      figure7_validate_generated_state_reference(
        temporary_reference,
        config,
        config_path = config_path
      )
    }
  }
  figure7_prepare_output(output_dir)
  for (file in list.files(file.path(source_dir, "tables"), full.names = TRUE)) {
    figure7_copy_file(file, file.path(output_dir, "tables", basename(file)))
  }
  if (include_state_pathway) {
    figure7_copy_file(
      provenance,
      file.path(output_dir, "metadata", basename(provenance))
    )
  }
  filenames <- stats::setNames(
    figure7_panel_filenames(config),
    figure7_panel_ids(TRUE)
  )
  legacy_plots <- list(
    A = figure7_panel_a_plot(a, config),
    B = figure7_panel_b_plot(b, bt),
    C = figure7_panel_c_plot(cdata, ct, config),
    D = figure7_scatter_plot(
      d,
      "shift_centered",
      "tgi_centered",
      dt,
      "CellCycle TGI association after within-dose centering",
      "Dose-centered ECDF RMSE",
      paste("Dose-centered Day", tgi_day, "TGI (%)")
    ) + ggplot2::geom_vline(
      xintercept = 0,
      color = "grey75",
      linewidth = 0.35
    ),
    E = figure7_scatter_plot(
      e,
      "sample_mean_endpoint_ploidy",
      tgi_measure,
      et,
      paste(
        "Cell-cycle-associated tumor cells: Day",
        tgi_day,
        "TGI vs sample mean ETP"
      ),
      "Sample mean ETP",
      paste("Day", tgi_day, "TGI (%)")
    )
  )
  sizes <- list(
    A = c(10, 6.5), B = c(15, 5.5), C = c(6.8, 6.4),
    D = c(6.6, 6.6), E = c(6.8, 6.8)
  )
  for (panel in names(legacy_plots)) {
    dimensions <- sizes[[panel]]
    figure7_save_panel(
      legacy_plots[[panel]],
      file.path(output_dir, "figures", filenames[[paste0("7", panel)]]),
      dimensions[[1L]],
      dimensions[[2L]]
    )
  }
  if (include_state_pathway) {
    legacy_plots$F <- figure7_panel_f_plot(f, config)
    figure7_save_panel(legacy_plots$F, file.path(output_dir, "figures", filenames[["7F"]]), 9, 8)
    figure7_save_main_composite(
      figure7_main_composite_plots(
        legacy_plots[LETTERS[1:5]],
        context_cache$plots,
        legacy_plots$F,
        config
      ),
      output_dir,
      config,
      context_cache$composite_filename
    )
  }
  source_results_root <- source_run_config$value[match("state_pathway_source_results_root", source_run_config$key)]
  if (length(source_results_root) != 1L || is.na(source_results_root) || identical(source_results_root, "not_recorded")) {
    source_results_root <- ""
  }
  write_metadata(
    output_dir,
    mode,
    config,
    config_path,
    panel_ids,
    source_results_root,
    reference_identity = reference_identity,
    context_cache = context_cache,
    cellcycle_input = source_scalar("cellcycle_input"),
    cellcycle_sha256 = source_scalar("cellcycle_sha256"),
    noncellcycle_input = source_scalar("noncellcycle_input"),
    noncellcycle_sha256 = source_scalar("noncellcycle_sha256")
  )
  figure7_validate_figure_inventory(
    output_dir,
    config,
    panel_ids,
    if (include_state_pathway) {
      context_cache$composite_filename
    } else {
      figure7_main_composite_filename(config)
    }
  )
}

if (identical(mode, "render-only")) {
  render_from_run(figure7_arg(args, "source-run-dir", required = TRUE), output_dir, config, config_path,
                  panel_ids, include_state_pathway)
  message("Rendered ", length(panel_ids), " Figure 7 panels from immutable plotting tables: ", output_dir)
  quit(save = "no", status = 0L)
}

workflow <- NULL
if (identical(mode, "full-workflow")) {
  overwrite_intermediates <- figure7_flag(args, "overwrite-intermediates", FALSE)
  workflow_paths <- figure7_workflow_paths(args, repo_root, output_dir, config)
  preflight <- figure7_preflight_workflow(
    workflow_paths,
    config_path,
    config,
    include_panel_f = include_state_pathway,
    overwrite_intermediates = overwrite_intermediates
  )
  if (figure7_flag(args, "preflight-only", FALSE)) {
    cat("workflow_state\t", preflight$state, "\n", sep = "")
    cat("raw_data_status\t", preflight$raw_data_status, "\n", sep = "")
    cat("raw_download_roles\t", if (length(preflight$raw_download_roles)) paste(preflight$raw_download_roles, collapse = ",") else "none", "\n", sep = "")
    cat("raw_validation_roles\t", if (length(preflight$raw_validation_roles)) paste(preflight$raw_validation_roles, collapse = ",") else "none", "\n", sep = "")
    cat("raw_data_dir\t", workflow_paths$raw_data_dir, "\n", sep = "")
    cat("download_workers\t", workflow_paths$download_workers, "\n", sep = "")
    cat("download_connections_per_file\t", workflow_paths$download_connections_per_file, "\n", sep = "")
    cat("loom_root\t", workflow_paths$loom_root, "\n", sep = "")
    cat("seurat_rds\t", workflow_paths$seurat_rds, "\n", sep = "")
    cat("scvelo_metrics\t", workflow_paths$scvelo_metrics, "\n", sep = "")
    cat("cell_pair_source\t", preflight$cell_pair_source, "\n", sep = "")
    cat("cellcycle_input\t", preflight$paths$cellcycle, "\n", sep = "")
    cat("noncellcycle_input\t", preflight$paths$noncellcycle, "\n", sep = "")
    cat(
      "state_pathway_reference\t",
      if (preflight$frozen_reference_ready) {
        preflight$paths$frozen_reference
      } else {
        preflight$paths$saved_reference
      },
      "\n",
      sep = ""
    )
    quit(save = "no", status = 0L)
  }
  figure7_assert_empty_output(output_dir)
  workflow <- figure7_prepare_full_workflow(
    args = args,
    paths = workflow_paths,
    preflight = preflight,
    script_dir = script_dir,
    config_path = config_path,
    config = config,
    include_panel_f = include_state_pathway,
    overwrite_intermediates = overwrite_intermediates
  )
  cellcycle_path <- workflow$cellcycle
  noncellcycle_path <- workflow$noncellcycle
  state_pathway_results_root <- workflow$state_pathway_root
} else {
  cellcycle_path <- normalizePath(figure7_arg(args, "cellcycle-input", required = TRUE), mustWork = FALSE)
  noncellcycle_path <- normalizePath(figure7_arg(args, "non-cellcycle-input", required = TRUE), mustWork = FALSE)
  figure7_assert_empty_output(output_dir)
}
if (identical(mode, "full-workflow")) {
  figure7_require_workflow_file(cellcycle_path, "cellcycle-input")
  figure7_require_workflow_file(noncellcycle_path, "non-cellcycle-input")
} else {
  figure7_verify_checksum(cellcycle_path, config$inputs$cellcycle_sha256, "CellCycle processed input")
  figure7_verify_checksum(noncellcycle_path, config$inputs$noncellcycle_sha256, "NonCellCycle processed input")
}
reference <- NULL
reference_identity <- NULL
if (include_state_pathway) {
  saved_dir <- if (identical(mode, "full-workflow")) {
    workflow$saved_reference
  } else {
    normalizePath(figure7_arg(args, "saved-state-pathway-dir", required = TRUE), mustWork = FALSE)
  }
  reference <- if (identical(mode, "full-workflow")) {
    if (!identical(
          workflow$state_pathway_reference_kind,
          as.character(config$state_pathways$generated_reference_kind)
        )) {
      figure7_stop(
        "Full-workflow panel 7F must use the generated human-only ",
        "noncanonical reference profile"
      )
    }
    figure7_validate_generated_state_reference(
      saved_dir,
      config,
      expected_inputs = list(
        cellcycle = cellcycle_path,
        noncellcycle = noncellcycle_path,
        seurat_rds = workflow$seurat_rds
      ),
      config_path = config_path
    )
  } else {
    saved_id <- basename(normalizePath(saved_dir, mustWork = FALSE))
    if (identical(
          saved_id,
          as.character(config$state_pathways$reviewed_reference_id)
        )) {
      figure7_validate_reviewed_state_reference(
        saved_dir,
        config,
        verify_checksums = TRUE
      )
    } else if (identical(
          saved_id,
          as.character(config$state_pathways$reference_id)
        )) {
      figure7_validate_historical_state_reference(
        saved_dir,
        config,
        verify_checksums = TRUE
      )
    } else {
      figure7_stop(
        "Standard panel-7F reference must be the reviewed human-only v2 ",
        "or the explicit historical mixed-v1 audit reference"
      )
    }
  }
  reference_identity <- list(
    id = reference$reference_id,
    kind = reference$reference_kind,
    canonical_publication_allowed =
      reference$canonical_publication_allowed
  )
}

if (identical(mode, "full-analysis")) {
  figure7_full_analysis(
    figure7_arg(args, "seurat-rds", required = TRUE),
    figure7_arg(args, "gene-set-artifact", required = TRUE),
    reference, output_dir, config
  )
}

cellcycle <- figure7_read_cell_table(cellcycle_path, "CellCycle", config)
noncellcycle <- figure7_read_cell_table(noncellcycle_path, "NonCellCycle", config)
samples <- figure7_sample_table(cellcycle, noncellcycle, config)
data <- figure7_prepare_cellcycle(cellcycle, samples, config)
context_cache <- if (include_state_pathway) {
  figure7_build_context_panels(
    si_cache_dir,
    repo_root,
    config,
    policy = si_cache_policy,
    upstream_input_manifest = si_cache_upstream_input_manifest,
    source_run_config = si_cache_source_run_config,
    source_provenance = si_cache_source_provenance
  )
} else {
  NULL
}
figure7_prepare_output(output_dir)
ae <- figure7_build_ae(cellcycle, data, samples, output_dir, config)
if (include_state_pathway) {
  plot_f <- if (identical(mode, "full-workflow")) {
    figure7_build_generated_f(reference, output_dir, config)
  } else {
    figure7_build_f(reference, output_dir, config)
  }
  figure7_write_tsv(
    context_cache$composition_result$plot_data,
    file.path(output_dir, "tables", "main_composite_panel_F_plot_data.tsv")
  )
  figure7_write_tsv(
    context_cache$composition_result$tests,
    file.path(
      output_dir,
      "tables",
      "main_composite_panel_F_enrichment_tests.tsv"
    )
  )
  figure7_save_main_composite(
    figure7_main_composite_plots(
      ae$plots,
      context_cache$plots,
      plot_f,
      config
    ),
    output_dir,
    config,
    context_cache$composite_filename
  )
}
write_metadata(
  output_dir,
  mode,
  config,
  config_path,
  panel_ids,
  state_pathway_results_root,
  workflow,
  reference_identity,
  context_cache,
  cellcycle_path,
  figure7_sha256(cellcycle_path),
  noncellcycle_path,
  figure7_sha256(noncellcycle_path)
)
figure7_validate_figure_inventory(
  output_dir,
  config,
  panel_ids,
  if (include_state_pathway) {
    context_cache$composite_filename
  } else {
    figure7_main_composite_filename(config)
  }
)
message(
  "Generated exactly ",
  length(panel_ids),
  " Figure 7 source panels",
  if (include_state_pathway) " and the assembled A-K main composite" else "",
  ": ",
  output_dir
)
