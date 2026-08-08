# Cache-first selection of the Seurat object shared by Figure 7 and
# Supplementary Figures 4-7. Scientific reconstruction remains implemented in
# generate_final_seurat_from_cellranger.R; this file only plans and validates
# reuse without loading a multi-gigabyte Seurat object into the planner.

figure7_seurat_upstream_final_path <- function(output_root) {
  file.path(
    output_root,
    "03_final_cluster",
    "integrated_sct_cca_seurat_final_reclustered.rds"
  )
}

figure7_seurat_upstream_manifest_path <- function(output_root) {
  file.path(output_root, "reconstruction_manifest.tsv")
}

figure7_standalone_seurat_final_path <- function(output_root) {
  file.path(
    output_root,
    "03_final_cluster",
    "03_objects",
    "integrated_sct_cca_seurat_final_reclustered.rds"
  )
}

figure7_standalone_provenance_path <- function(output_root, filename) {
  file.path(output_root, "00_provenance", filename)
}

figure7_read_standalone_key_values <- function(path, label) {
  if (!file.exists(path)) figure7_stop(label, " is missing: ", path)
  table <- utils::read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    colClasses = "character"
  )
  empty_values <- !nzchar(table$value)
  allowed_empty_fields <- "private_r_library"
  if (!identical(names(table), c("field", "value")) ||
      !nrow(table) || anyNA(table) || any(!nzchar(table$field)) ||
      any(empty_values & !table$field %in% allowed_empty_fields) ||
      anyDuplicated(table$field)) {
    figure7_stop(label, " has an invalid field/value schema: ", path)
  }
  stats::setNames(table$value, table$field)
}

figure7_read_standalone_completion <- function(path) {
  if (!file.exists(path)) {
    figure7_stop("Standalone completion marker is missing: ", path)
  }
  lines <- readLines(path, warn = FALSE)
  if (!length(lines) || any(!grepl("=", lines, fixed = TRUE))) {
    figure7_stop("Standalone completion marker is malformed: ", path)
  }
  keys <- sub("=.*$", "", lines)
  values <- sub("^[^=]*=", "", lines)
  if (any(!nzchar(keys)) || any(!nzchar(values)) || anyDuplicated(keys)) {
    figure7_stop("Standalone completion marker is malformed: ", path)
  }
  stats::setNames(values, keys)
}

figure7_standalone_completion_contract_sha256 <- function(completion) {
  stable_fields <- completion[setdiff(names(completion), "completed_at")]
  figure7_contract_sha256(stable_fields)
}

figure7_validate_semantic_rds_audit <- function(output_root, final_path) {
  audit_root <- file.path(
    output_root,
    "00_validation",
    "rds_semantic_audit"
  )
  marker_path <- file.path(audit_root, "AUDIT_COMPLETE.txt")
  marker <- figure7_read_standalone_completion(marker_path)
  report_names <- c(
    "audit_summary.tsv", "metadata_comparison.tsv",
    "cluster_comparison.tsv", "graph_comparison.tsv",
    "assay_numeric_comparison.tsv", "pca_comparison.tsv",
    "umap_displacement_summary.tsv", "umap_largest_displacements.tsv",
    "command_comparison.tsv", "rds_file_identity.tsv"
  )
  report_keys <- paste0("report_sha256.", report_names)
  required_keys <- c(
    "schema_version", "status", "checks_total", "checks_passed",
    "checks_failed", "generated_rds", "generated_rds_size_bytes",
    "generated_rds_md5", "generated_rds_sha256",
    "zenodo_reference_rds", "zenodo_reference_rds_size_bytes",
    "zenodo_reference_rds_md5", "zenodo_reference_rds_sha256",
    "downstream_rds", "audit_summary", "audit_summary_sha256",
    report_keys
  )
  if (!all(required_keys %in% names(marker)) ||
      !identical(marker[["schema_version"]], "semantic_rds_audit_v1") ||
      !identical(marker[["status"]], "PASS") ||
      !identical(marker[["checks_failed"]], "0")) {
    figure7_stop("Semantic RDS audit completion marker is invalid")
  }
  final_path <- normalizePath(final_path, mustWork = TRUE)
  if (!identical(
        normalizePath(marker[["generated_rds"]], mustWork = TRUE),
        final_path
      ) ||
      !identical(
        normalizePath(marker[["downstream_rds"]], mustWork = TRUE),
        final_path
      ) ||
      !identical(
        marker[["generated_rds_size_bytes"]],
        as.character(file.info(final_path)$size)
      ) ||
      !identical(
        tolower(marker[["generated_rds_md5"]]),
        tolower(unname(tools::md5sum(final_path)))
      ) ||
      !identical(
        tolower(marker[["generated_rds_sha256"]]),
        figure7_sha256(final_path)
      )) {
    figure7_stop("Generated RDS no longer matches its semantic audit")
  }
  reference_path <- normalizePath(
    marker[["zenodo_reference_rds"]],
    mustWork = TRUE
  )
  if (!identical(
        marker[["zenodo_reference_rds_size_bytes"]],
        as.character(file.info(reference_path)$size)
      ) ||
      !identical(
        tolower(marker[["zenodo_reference_rds_md5"]]),
        tolower(unname(tools::md5sum(reference_path)))
      ) ||
      !identical(
        tolower(marker[["zenodo_reference_rds_sha256"]]),
        figure7_sha256(reference_path)
      )) {
    figure7_stop("Zenodo reference RDS no longer matches the semantic audit")
  }
  summary_path <- normalizePath(marker[["audit_summary"]], mustWork = TRUE)
  expected_summary <- normalizePath(
    file.path(audit_root, "audit_summary.tsv"),
    mustWork = TRUE
  )
  summary <- utils::read.delim(
    summary_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    colClasses = "character"
  )
  if (!identical(summary_path, expected_summary) ||
      !identical(
        names(summary),
        c(
          "group", "check", "status", "observed", "reference",
          "threshold", "details"
        )
      ) ||
      !nrow(summary) || anyNA(summary) || any(summary$status != "PASS") ||
      !identical(as.character(nrow(summary)), marker[["checks_total"]]) ||
      !identical(as.character(nrow(summary)), marker[["checks_passed"]]) ||
      !identical(
        tolower(marker[["audit_summary_sha256"]]),
        figure7_sha256(summary_path)
      )) {
    figure7_stop("Semantic RDS audit summary is inconsistent")
  }
  identity_path <- file.path(audit_root, "rds_file_identity.tsv")
  identity <- utils::read.delim(
    identity_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    colClasses = "character"
  )
  if (!identical(
        names(identity),
        c("artifact", "path", "size_bytes", "md5", "sha256", "role")
      ) ||
      !identical(identity$artifact, c("generated", "zenodo_reference")) ||
      !identical(
        normalizePath(identity$path, mustWork = TRUE),
        c(final_path, reference_path)
      ) ||
      !identical(
        tolower(identity$md5),
        tolower(c(
          marker[["generated_rds_md5"]],
          marker[["zenodo_reference_rds_md5"]]
        ))
      ) ||
      !identical(
        tolower(identity$sha256),
        tolower(c(
          marker[["generated_rds_sha256"]],
          marker[["zenodo_reference_rds_sha256"]]
        ))
      )) {
    figure7_stop("Semantic RDS audit file-identity table is inconsistent")
  }
  report_paths <- file.path(audit_root, report_names)
  if (any(!file.exists(report_paths)) || any(!vapply(
        seq_along(report_paths),
        function(index) identical(
          tolower(marker[[report_keys[[index]]]]),
          figure7_sha256(report_paths[[index]])
        ),
        logical(1L)
      ))) {
    figure7_stop("Semantic RDS audit report set is incomplete or modified")
  }
  list(
    marker = normalizePath(marker_path, mustWork = TRUE),
    marker_sha256 = figure7_sha256(marker_path),
    summary = summary_path,
    summary_sha256 = figure7_sha256(summary_path),
    identity = normalizePath(identity_path, mustWork = TRUE),
    identity_sha256 = figure7_sha256(identity_path),
    reference_rds = reference_path,
    reference_rds_sha256 = figure7_sha256(reference_path),
    generated_rds_md5 = tolower(marker[["generated_rds_md5"]]),
    reference_rds_md5 = tolower(marker[["zenodo_reference_rds_md5"]])
  )
}

figure7_validate_standalone_seurat_artifact <- function(
  output_root,
  module_dir,
  all_ploidy,
  sample_info,
  expected_rds = "",
  cellranger_root = ""
) {
  output_root <- normalizePath(output_root, mustWork = TRUE)
  final_path <- normalizePath(
    figure7_standalone_seurat_final_path(output_root),
    mustWork = TRUE
  )
  if (nzchar(expected_rds) &&
      !identical(normalizePath(expected_rds, mustWork = TRUE), final_path)) {
    figure7_stop(
      "Explicit RDS is not the standalone pipeline's recorded final object"
    )
  }

  run_manifest_path <- figure7_standalone_provenance_path(
    output_root,
    "run_manifest.tsv"
  )
  run_manifest <- utils::read.delim(
    run_manifest_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    colClasses = "character"
  )
  required_columns <- c(
    "input_type", "path", "size_bytes", "md5", "modified_time"
  )
  if (!identical(names(run_manifest), required_columns) ||
      nrow(run_manifest) != 20L || anyNA(run_manifest) ||
      sum(run_manifest$input_type == "filtered_feature_bc_matrix.h5") != 18L ||
      sum(run_manifest$input_type == "all_ploidy.tsv") != 1L ||
      sum(run_manifest$input_type == "sample_info.xlsx") != 1L ||
      anyDuplicated(run_manifest$path) ||
      any(!grepl("^[0-9a-f]{32}$", tolower(run_manifest$md5)))) {
    figure7_stop(
      "Standalone run manifest must bind exactly 18 H5 files, all_ploidy.tsv, ",
      "and sample_info.xlsx"
    )
  }
  normalized_inputs <- vapply(
    run_manifest$path,
    normalizePath,
    character(1),
    mustWork = TRUE
  )
  observed_size <- as.character(file.info(normalized_inputs)$size)
  observed_md5 <- tolower(unname(tools::md5sum(normalized_inputs)))
  if (any(run_manifest$size_bytes != observed_size) ||
      any(tolower(run_manifest$md5) != observed_md5)) {
    figure7_stop("Standalone run manifest input size/MD5 verification failed")
  }
  input_bindings <- c(
    all_ploidy.tsv = normalizePath(all_ploidy, mustWork = TRUE),
    sample_info.xlsx = normalizePath(sample_info, mustWork = TRUE)
  )
  for (input_type in names(input_bindings)) {
    if (!identical(
      unname(normalized_inputs[run_manifest$input_type == input_type]),
      unname(input_bindings[[input_type]])
    )) {
      figure7_stop("Standalone manifest binds a different ", input_type)
    }
  }
  h5_paths <- normalized_inputs[
    run_manifest$input_type == "filtered_feature_bc_matrix.h5"
  ]
  if (any(!grepl(
    "-Count-HM_filtered_feature_bc_matrix\\.h5$",
    basename(h5_paths)
  ))) {
    figure7_stop("Standalone manifest contains a noncanonical H5 filename")
  }
  if (nzchar(cellranger_root)) {
    cellranger_root <- normalizePath(cellranger_root, mustWork = TRUE)
    if (any(!startsWith(h5_paths, paste0(cellranger_root, .Platform$file.sep)))) {
      figure7_stop("Standalone H5 provenance escapes the selected Cell Ranger root")
    }
  }

  runtime_path <- figure7_standalone_provenance_path(
    output_root,
    "final_artifact_runtime.tsv"
  )
  runtime <- figure7_read_standalone_key_values(
    runtime_path,
    "Standalone final-artifact runtime marker"
  )
  required_runtime <- c(
    status = "PASS",
    phase = "post_filter",
    final_object = final_path,
    final_object_size_bytes = as.character(file.info(final_path)$size),
    final_object_md5 = tolower(unname(tools::md5sum(final_path)))
  )
  if (!all(c(
        names(required_runtime), "active_sif_md5", "R", "Seurat"
      ) %in% names(runtime)) ||
      !identical(
        unname(runtime[names(required_runtime)]),
        unname(required_runtime)
      ) || !grepl("^[0-9a-f]{32}$", tolower(runtime[["active_sif_md5"]]))) {
    figure7_stop("Standalone final-artifact runtime marker is inconsistent")
  }
  completion_path <- figure7_standalone_provenance_path(
    output_root,
    "PIPELINE_COMPLETE.txt"
  )
  completion <- figure7_read_standalone_completion(completion_path)
  if (!identical(completion[["status"]], "PASS") ||
      !identical(
        completion[["final_object"]],
        "03_final_cluster/03_objects/integrated_sct_cca_seurat_final_reclustered.rds"
      )) {
    figure7_stop("Standalone pipeline completion marker is inconsistent")
  }
  semantic_audit <- figure7_validate_semantic_rds_audit(
    output_root,
    final_path
  )

  dependency_names <- ifelse(
    run_manifest$input_type == "filtered_feature_bc_matrix.h5",
    paste0(
      "h5_",
      sub("-Count-HM_filtered_feature_bc_matrix\\.h5$", "", basename(normalized_inputs))
    ),
    sub("\\.[^.]+$", "", run_manifest$input_type)
  )
  dependencies <- stats::setNames(
    vapply(normalized_inputs, figure7_sha256, character(1)),
    dependency_names
  )
  dependencies <- c(
    dependencies,
    semantic_audit_reference_rds = semantic_audit$reference_rds_sha256
  )
  if (anyDuplicated(names(dependencies))) {
    figure7_stop("Standalone input dependency roles are not unique")
  }
  standalone_sources <- file.path(
    dirname(module_dir),
    "scRNA_Seq_analysis",
    c(
      "run_cluster_standalone.sh",
      "cluster_pipeline_standalone.R",
      "bootstrap_dependencies.R",
      "audit_generated_seurat_rds.R"
    )
  )
  standalone_sources <- normalizePath(standalone_sources, mustWork = TRUE)
  scientific_code_contracts <- stats::setNames(
    vapply(standalone_sources, figure7_sha256, character(1)),
    c(
      "scientific_code_contract:standalone_orchestrator",
      "scientific_code_contract:standalone_pipeline",
      "scientific_code_contract:standalone_bootstrap",
      "scientific_code_contract:semantic_rds_audit"
    )
  )
  common_contract <- c(
    schema_version = "standalone_cluster_v1",
    post_filter_r = runtime[["R"]],
    post_filter_seurat = runtime[["Seurat"]],
    post_filter_sif_md5 = tolower(runtime[["active_sif_md5"]]),
    completion_contract_sha256 =
      figure7_standalone_completion_contract_sha256(completion),
    semantic_audit_marker_sha256 = semantic_audit$marker_sha256,
    semantic_audit_summary_sha256 = semantic_audit$summary_sha256,
    semantic_audit_identity_sha256 = semantic_audit$identity_sha256,
    generated_rds_md5 = semantic_audit$generated_rds_md5,
    zenodo_reference_rds_md5 = semantic_audit$reference_rds_md5
  )
  list(
    valid = TRUE,
    output_root = output_root,
    rds = final_path,
    rds_sha256 = figure7_sha256(final_path),
    final_stage_manifest = runtime_path,
    final_stage_manifest_sha256 = figure7_sha256(runtime_path),
    reconstruction_manifest = run_manifest_path,
    reconstruction_manifest_sha256 = figure7_sha256(run_manifest_path),
    dependencies = dependencies,
    scientific_code_contracts = scientific_code_contracts,
    upstream_common_contract = common_contract,
    expected_dependencies = dependencies
  )
}

figure7_probe_standalone_seurat_artifact <- function(...) {
  tryCatch(
    figure7_validate_standalone_seurat_artifact(...),
    error = function(error) list(valid = FALSE, error = conditionMessage(error))
  )
}

figure7_seurat_upstream_generator_environment <- local({
  cached <- NULL
  cached_script <- ""
  function(module_dir) {
    script <- normalizePath(
      file.path(module_dir, "generate_final_seurat_from_cellranger.R"),
      mustWork = TRUE
    )
    if (is.null(cached) || !identical(script, cached_script)) {
      cached <<- new.env(parent = .GlobalEnv)
      sys.source(script, envir = cached)
      cached_script <<- script
    }
    cached
  }
})

figure7_validate_seurat_upstream_artifact <- function(
  output_root,
  module_dir,
  environment_lock,
  config,
  all_ploidy,
  sample_info,
  expected_rds = "",
  cellranger_root = ""
) {
  generator <- figure7_seurat_upstream_generator_environment(module_dir)
  generator$figure7_upstream_validate_artifact(
    output_root = output_root,
    environment_lock = environment_lock,
    config = config,
    all_ploidy = all_ploidy,
    sample_info = sample_info,
    expected_rds = expected_rds,
    cellranger_root = cellranger_root
  )
}

figure7_probe_seurat_upstream_artifact <- function(...) {
  tryCatch(
    figure7_validate_seurat_upstream_artifact(...),
    error = function(error) {
      list(valid = FALSE, error = conditionMessage(error))
    }
  )
}

figure7_probe_any_seurat_upstream_artifact <- function(
  output_root,
  module_dir,
  environment_lock,
  config,
  all_ploidy,
  sample_info,
  expected_rds = "",
  cellranger_root = ""
) {
  standalone <- figure7_probe_standalone_seurat_artifact(
    output_root = output_root,
    module_dir = module_dir,
    all_ploidy = all_ploidy,
    sample_info = sample_info,
    expected_rds = expected_rds,
    cellranger_root = cellranger_root
  )
  if (isTRUE(standalone$valid)) return(standalone)
  legacy <- figure7_probe_seurat_upstream_artifact(
    output_root = output_root,
    module_dir = module_dir,
    environment_lock = environment_lock,
    config = config,
    all_ploidy = all_ploidy,
    sample_info = sample_info,
    expected_rds = expected_rds,
    cellranger_root = cellranger_root
  )
  if (isTRUE(legacy$valid)) return(legacy)
  list(
    valid = FALSE,
    error = paste0(
      "standalone: ", standalone$error,
      "; figure7 reconstruction: ", legacy$error
    )
  )
}

figure7_validate_any_seurat_upstream_artifact <- function(...) {
  result <- figure7_probe_any_seurat_upstream_artifact(...)
  if (!isTRUE(result$valid)) {
    figure7_stop(
      "Seurat source lacks a valid standalone or Figure 7 reconstruction ",
      "chain: ", result$error
    )
  }
  result
}

figure7_probe_resumable_seurat_upstream_artifact <- function(
  output_root,
  module_dir,
  environment_lock,
  config,
  all_ploidy,
  sample_info,
  cellranger_root = ""
) {
  generator <- figure7_seurat_upstream_generator_environment(module_dir)
  generator$figure7_upstream_probe_resumable_artifact(
    output_root = output_root,
    environment_lock = environment_lock,
    config = config,
    all_ploidy = all_ploidy,
    sample_info = sample_info,
    cellranger_root = cellranger_root
  )
}

figure7_infer_seurat_upstream_root <- function(rds_path) {
  if (!nzchar(rds_path)) return("")
  expected_name <- basename(figure7_seurat_upstream_final_path("root"))
  if (!identical(basename(rds_path), expected_name)) return("")
  if (identical(basename(dirname(rds_path)), "03_objects") &&
      identical(basename(dirname(dirname(rds_path))), "03_final_cluster")) {
    return(normalizePath(dirname(dirname(dirname(rds_path))), mustWork = FALSE))
  }
  normalizePath(dirname(dirname(rds_path)), mustWork = FALSE)
}

figure7_seurat_source_dependency_values <- function(
  rds_sha256,
  transitive_dependencies = character(),
  scientific_code_contracts = character(),
  upstream_common_contract = character()
) {
  values <- c(source_seurat_rds = tolower(as.character(rds_sha256)))
  reconstructed <- length(transitive_dependencies) ||
    length(scientific_code_contracts) ||
    length(upstream_common_contract)
  if (reconstructed) {
    if (is.null(names(transitive_dependencies)) ||
        is.null(names(scientific_code_contracts)) ||
        is.null(names(upstream_common_contract)) ||
        !length(transitive_dependencies) ||
        !length(scientific_code_contracts) ||
        !length(upstream_common_contract)) {
      figure7_stop(
        "Reconstructed Seurat lineage is missing transitive scientific ",
        "dependencies"
      )
    }
    values <- c(
      values,
      stats::setNames(
        tolower(as.character(transitive_dependencies)),
        paste0(
          "seurat_transitive_dependency:",
          names(transitive_dependencies)
        )
      ),
      stats::setNames(
        tolower(as.character(scientific_code_contracts)),
        paste0(
          "seurat_",
          names(scientific_code_contracts)
        )
      ),
      seurat_upstream_scientific_common_contract =
        figure7_contract_sha256(upstream_common_contract)
    )
  }
  if (is.null(names(values)) ||
      anyNA(names(values)) ||
      any(!nzchar(names(values))) ||
      anyDuplicated(names(values)) ||
      any(!grepl("^[0-9a-f]{64}$", values))) {
    figure7_stop(
      "Seurat source lineage requires unique SHA-256 dependency roles"
    )
  }
  values
}

figure7_seurat_selection_dependency_values <- function(selection) {
  figure7_seurat_source_dependency_values(
    selection$rds_sha256,
    selection$transitive_dependencies,
    selection$scientific_code_contracts,
    selection$upstream_common_contract
  )
}

figure7_seurat_selection_audit_values <- function(selection) {
  values <- c(
    if (nzchar(selection$final_stage_manifest_sha256)) c(
      audit_seurat_final_stage_manifest =
        selection$final_stage_manifest_sha256
    ),
    if (nzchar(selection$reconstruction_manifest_sha256)) c(
      audit_seurat_reconstruction_manifest =
        selection$reconstruction_manifest_sha256
    )
  )
  value_names <- names(values)
  values <- tolower(as.character(values))
  names(values) <- value_names
  values
}

figure7_select_seurat_source <- function(
  needs_seurat,
  explicit_rds,
  deposited_rds,
  deposited_sha256,
  upstream_dir,
  cellranger_root,
  module_dir,
  environment_lock,
  config,
  all_ploidy,
  sample_info,
  required_source_kind = "auto"
) {
  empty <- list(
    required = FALSE,
    action = "not_required",
    source = "not_required",
    rds = "",
    rds_sha256 = "",
    final_stage_manifest = "",
    final_stage_manifest_sha256 = "",
    reconstruction_manifest = "",
    reconstruction_manifest_sha256 = "",
    transitive_dependencies = character(),
    scientific_code_contracts = character(),
    upstream_common_contract = character(),
    resumable_stage = "",
    resumable_stage_manifest = "",
    resumable_stage_manifest_sha256 = "",
    upstream_dir = upstream_dir,
    cellranger_root = cellranger_root
  )
  if (!isTRUE(needs_seurat)) return(empty)

  required_source_kind <- as.character(required_source_kind)
  if (!required_source_kind %in% c(
    "auto",
    "deposited",
    "manifested_reconstruction"
  )) {
    figure7_stop(
      "required_source_kind must be auto, deposited, or ",
      "manifested_reconstruction"
    )
  }
  deposited_sha256 <- tolower(as.character(deposited_sha256))
  if (!grepl("^[0-9a-f]{64}$", deposited_sha256)) {
    figure7_stop("Deposited Seurat checksum contract is invalid")
  }
  if (nzchar(cellranger_root)) {
    cellranger_root <- normalizePath(cellranger_root, mustWork = FALSE)
  }
  artifact_result <- function(validation, source, action = "reuse") {
    list(
      required = TRUE,
      action = action,
      source = source,
      rds = validation$rds,
      rds_sha256 = validation$rds_sha256,
      final_stage_manifest = validation$final_stage_manifest,
      final_stage_manifest_sha256 =
        validation$final_stage_manifest_sha256,
      reconstruction_manifest = validation$reconstruction_manifest,
      reconstruction_manifest_sha256 =
        validation$reconstruction_manifest_sha256,
      transitive_dependencies = validation$dependencies,
      scientific_code_contracts =
        validation$scientific_code_contracts,
      upstream_common_contract =
        validation$upstream_common_contract,
      resumable_stage = "",
      resumable_stage_manifest = "",
      resumable_stage_manifest_sha256 = "",
      upstream_dir = validation$output_root,
      cellranger_root = cellranger_root
    )
  }
  deposited_result <- function(path, source = NULL) {
    path <- normalizePath(path, mustWork = FALSE)
    exists <- file.exists(path)
    list(
      required = TRUE,
      action = if (exists) "validate_deposited" else "download_deposited",
      source = if (!is.null(source)) {
        source
      } else if (exists) {
        "deposited_cache"
      } else {
        "planned_deposited_download"
      },
      rds = path,
      rds_sha256 = deposited_sha256,
      final_stage_manifest = "",
      final_stage_manifest_sha256 = "",
      reconstruction_manifest = "",
      reconstruction_manifest_sha256 = "",
      transitive_dependencies = character(),
      scientific_code_contracts = character(),
      upstream_common_contract = character(),
      resumable_stage = "",
      resumable_stage_manifest = "",
      resumable_stage_manifest_sha256 = "",
      upstream_dir = upstream_dir,
      cellranger_root = ""
    )
  }

  if (nzchar(explicit_rds)) {
    explicit_rds <- normalizePath(explicit_rds, mustWork = FALSE)
    if (!file.exists(explicit_rds)) {
      figure7_stop("Explicit --seurat-rds does not exist: ", explicit_rds)
    }
    explicit_hash <- figure7_sha256(explicit_rds)
    if (identical(required_source_kind, "deposited")) {
      if (!identical(explicit_hash, deposited_sha256)) {
        figure7_stop(
          "The required deposited Seurat source has the wrong checksum"
        )
      }
      result <- deposited_result(
        explicit_rds,
        "explicit_deposited_rds"
      )
      result$action <- "reuse"
      result$rds_sha256 <- explicit_hash
      return(result)
    }
    candidate_roots <- unique(c(
      upstream_dir,
      figure7_infer_seurat_upstream_root(explicit_rds)
    ))
    candidate_roots <- candidate_roots[
      nzchar(candidate_roots) & dir.exists(candidate_roots)
    ]
    failures <- character()
    for (candidate in candidate_roots) {
      validation <- figure7_probe_any_seurat_upstream_artifact(
        output_root = candidate,
        module_dir = module_dir,
        environment_lock = environment_lock,
        config = config,
        all_ploidy = all_ploidy,
        sample_info = sample_info,
        expected_rds = explicit_rds,
        cellranger_root = if (dir.exists(cellranger_root)) {
          cellranger_root
        } else {
          ""
        }
      )
      if (isTRUE(validation$valid)) {
        return(artifact_result(
          validation,
          "explicit_reconstructed_rds"
        ))
      }
      failures <- c(failures, paste0(candidate, ": ", validation$error))
    }
    if (identical(required_source_kind, "auto") &&
        identical(explicit_hash, deposited_sha256)) {
      result <- deposited_result(
        explicit_rds,
        "explicit_deposited_rds"
      )
      result$action <- "reuse"
      result$rds_sha256 <- explicit_hash
      return(result)
    }
    figure7_stop(
      "Explicit --seurat-rds lacks the required valid reconstruction ",
      "manifest/final-stage chain",
      if (length(failures)) paste0(": ", paste(failures, collapse = " | ")) else ""
    )
  }

  if (identical(required_source_kind, "deposited")) {
    return(deposited_result(deposited_rds))
  }

  if (nzchar(upstream_dir) && dir.exists(upstream_dir)) {
    validation <- figure7_probe_any_seurat_upstream_artifact(
      output_root = upstream_dir,
      module_dir = module_dir,
      environment_lock = environment_lock,
      config = config,
      all_ploidy = all_ploidy,
      sample_info = sample_info,
      cellranger_root = if (dir.exists(cellranger_root)) {
        cellranger_root
      } else {
        ""
      }
    )
    if (isTRUE(validation$valid)) {
      return(artifact_result(validation, "reused_reconstructed_rds"))
    }
    message(
      "[seurat-source] Existing upstream cache is not reusable: ",
      validation$error
    )
    resumable <- figure7_probe_resumable_seurat_upstream_artifact(
      output_root = upstream_dir,
      module_dir = module_dir,
      environment_lock = environment_lock,
      config = config,
      all_ploidy = all_ploidy,
      sample_info = sample_info,
      cellranger_root = if (dir.exists(cellranger_root)) {
        cellranger_root
      } else {
        ""
      }
    )
    if (isTRUE(resumable$valid)) {
      return(list(
        required = TRUE,
        action = "build_upstream",
        source = paste0(
          "planned_upstream_resume_",
          resumable$stage
        ),
        rds = figure7_seurat_upstream_final_path(upstream_dir),
        rds_sha256 = "",
        final_stage_manifest = "",
        final_stage_manifest_sha256 = "",
        reconstruction_manifest =
          figure7_seurat_upstream_manifest_path(upstream_dir),
        reconstruction_manifest_sha256 = "",
        transitive_dependencies = resumable$dependencies,
        scientific_code_contracts = character(),
        upstream_common_contract = character(),
        resumable_stage = resumable$stage,
        resumable_stage_manifest = resumable$stage_manifest,
        resumable_stage_manifest_sha256 =
          resumable$stage_manifest_sha256,
        upstream_dir = upstream_dir,
        cellranger_root = cellranger_root
      ))
    }
    message(
      "[seurat-source] Existing upstream cache cannot resume: ",
      resumable$error
    )
  }

  if (nzchar(cellranger_root)) {
    if (!nzchar(upstream_dir)) {
      figure7_stop(
        "Cell Ranger reconstruction requires a Seurat-upstream output directory"
      )
    }
    if (!dir.exists(cellranger_root)) {
      figure7_stop(
        "Explicit Cell Ranger root does not exist and no reusable ",
        "upstream stage is available: ",
        cellranger_root
      )
    }
    return(list(
      required = TRUE,
      action = "build_upstream",
      source = "planned_cellranger_reconstruction",
      rds = figure7_seurat_upstream_final_path(upstream_dir),
      rds_sha256 = "",
      final_stage_manifest = "",
      final_stage_manifest_sha256 = "",
      reconstruction_manifest =
        figure7_seurat_upstream_manifest_path(upstream_dir),
      reconstruction_manifest_sha256 = "",
      transitive_dependencies = character(),
      scientific_code_contracts = character(),
      upstream_common_contract = character(),
      resumable_stage = "",
      resumable_stage_manifest = "",
      resumable_stage_manifest_sha256 = "",
      upstream_dir = upstream_dir,
      cellranger_root = cellranger_root
    ))
  }

  if (identical(
    required_source_kind,
    "manifested_reconstruction"
  )) {
    figure7_stop(
      "A manifested reconstructed Seurat source is required, but no valid ",
      "upstream cache or available Cell Ranger root can provide it"
    )
  }
  deposited_result(deposited_rds)
}

figure7_seurat_upstream_generator_args <- function(
  selection,
  environment_lock,
  config_path,
  all_ploidy,
  sample_info,
  jobs
) {
  if (!identical(selection$action, "build_upstream")) {
    figure7_stop("Seurat-upstream generator args require a build plan")
  }
  list(
    "output-dir" = selection$upstream_dir,
    "cellranger-root" = selection$cellranger_root,
    "environment-lock" = environment_lock,
    config = config_path,
    "all-ploidy" = all_ploidy,
    "sample-info" = sample_info,
    jobs = as.character(jobs),
    "force-rebuild" = "false"
  )
}
