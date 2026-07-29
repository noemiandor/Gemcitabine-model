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
      validation <- figure7_probe_seurat_upstream_artifact(
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
    validation <- figure7_probe_seurat_upstream_artifact(
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
