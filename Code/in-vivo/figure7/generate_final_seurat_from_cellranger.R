#!/usr/bin/env Rscript

# Cache-first reconstruction of the sole Seurat object consumed by Figure 7
# and Supplementary Figures 4-7.  The external raw boundary is 18 Cell Ranger
# filtered_feature_bc_matrix.h5 files.  FASTQ-to-Cell-Ranger and the creation
# of all_ploidy.tsv are not represented in the available analysis history.
#
# Tao's reviewed scripts use a sequential future plan.  The orchestration
# --jobs value is accepted and audited but never changes this single-worker
# scientific execution contract.

figure7_upstream_script_path <- function() {
  frame_files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) NA_character_ else as.character(frame$ofile)
  }, character(1L))
  command_file <- sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)
  )
  candidates <- c(
    rev(frame_files[!is.na(frame_files) & nzchar(frame_files)]),
    command_file,
    file.path(
      Sys.getenv("FIGURE7_MODULE_DIR", unset = ""),
      "generate_final_seurat_from_cellranger.R"
    ),
    "Code/in-vivo/figure7/generate_final_seurat_from_cellranger.R"
  )
  candidates <- candidates[
    basename(candidates) ==
      "generate_final_seurat_from_cellranger.R"
  ]
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) {
    stop("Cannot resolve Seurat-upstream generator path", call. = FALSE)
  }
  normalizePath(existing[[1L]], mustWork = TRUE)
}
script_path <- figure7_upstream_script_path()
script_dir <- dirname(script_path)
repo_root <- normalizePath(
  file.path(script_dir, "..", "..", ".."),
  mustWork = TRUE
)
common_path <- file.path(script_dir, "src", "common_io.R")
helper_path <- file.path(script_dir, "src", "seurat_upstream.R")
sys.source(common_path, envir = .GlobalEnv)
sys.source(helper_path, envir = .GlobalEnv)

figure7_upstream_bool <- function(value, label) {
  normalized <- tolower(trimws(as.character(value)))
  if (!normalized %in% c("true", "false", "1", "0", "yes", "no")) {
    figure7_stop(label, " must be true or false")
  }
  normalized %in% c("true", "1", "yes")
}

figure7_upstream_resolve_path <- function(
  path,
  must_work = FALSE
) {
  if (!grepl("^/", path)) path <- file.path(repo_root, path)
  normalizePath(path, mustWork = must_work)
}

figure7_upstream_legacy_sources <- c(
  `01_data.R` =
    "1148912de8feee2321e89440570a4b51b628f5421e12cf4ba8d47ba0d757c6d9",
  `01a_cell_cycle.R` =
    "558d40178fbba21fac2eac3f09972b464bdb5b90a716eed9e071cfa8a72a5c87",
  `02b_cluster_refine.R` =
    "c37da49ed1f5bf8740652cf3cd51d1b36e6d635b45104ef99192baa36e476551",
  `02d_manual_cluster_merge.R` =
    "ef033f9984e3937b90c12886a2f2b1682c411825d01bc4194049c9cbec3d56a6",
  `03_final_cluster.R` =
    "47e255458ad6054a502b15db13b60ee21c25ae48c960507cfcfd4090fc3e8b6e",
  `Utils.R` =
    "0145103a7e2a1c16bac2531037e182e0010a4ddf9cc3f944f2076f77108a2ba0"
)

figure7_upstream_stage_definitions <- function(output_root) {
  list(
    integrated = list(
      directory = file.path(output_root, "01_integrated"),
      filename = "integrated_sct_cca_seurat.rds",
      parameter_contract = paste(
        "samples=18",
        "special=2N-Cell-Culture,4N-Cell-Culture",
        "special_qc=nFeature>=200,nCount>=500,percent.mt<=20,upper_q=0.99",
        "tumor_filter=all_ploidy_cell_id",
        "SCTransform=percent.mt,return_all_genes",
        "integration=SCT-CCA,features=3000",
        "PCA=50,neighbors=1:30,clusters=0.6,UMAP=1:30",
        "future=sequential",
        "seed=1234",
        sep = ";"
      )
    ),
    cell_cycle = list(
      directory = file.path(output_root, "01a_cell_cycle"),
      filename = "integrated_sct_cca_seurat_cell_cycle.rds",
      parameter_contract = paste(
        "assay=RNA",
        "genes=Seurat::cc.genes.updated.2019",
        "candidate=score_outlier+phase_enrichment+PC_support>=2",
        "PCs=1:20,cor_threshold=0.30",
        "seed=12345",
        sep = ";"
      )
    ),
    refined = list(
      directory = file.path(output_root, "02b_cluster_refine"),
      filename = "integrated_sct_cca_seurat_cluster_refine.rds",
      parameter_contract = paste(
        "core=6,10,11,12",
        "neighbors=0,3",
        "candidates=4:4c,9:9c",
        "UMAP_k=50,region_fraction=0.60,core_fraction=0.20,hull=true",
        "component_k=5,component_mutual=false",
        "threads=1",
        "seed=12345",
        sep = ";"
      )
    ),
    merged = list(
      directory = file.path(output_root, "02d_manual_cluster_merge"),
      filename =
        "integrated_sct_cca_seurat_cluster_refine_manual_merge_test.rds",
      parameter_contract = paste(
        "normalize_RNA_if_empty=true",
        "merge_0=0,1,7",
        "merge_10=10,11,12",
        "future=sequential,threads=1",
        "seed=12345",
        sep = ";"
      )
    ),
    final = list(
      directory = file.path(output_root, "03_final_cluster"),
      filename = "integrated_sct_cca_seurat_final_reclustered.rds",
      parameter_contract = paste(
        "remove=9,4,3,9c",
        "rename=manual_merge_test:clusters",
        "derive=Ploidy,TN_from_IDs",
        "assay_preference=integrated,SCT,RNA",
        "PCA=50,UMAP=1:30,no_neighbors,no_clusters",
        "future=sequential,threads=1",
        "seed=1234",
        sep = ";"
      )
    )
  )
}

figure7_upstream_stage_output <- function(definition) {
  file.path(definition$directory, definition$filename)
}

figure7_upstream_stage_manifest <- function(definition) {
  file.path(definition$directory, "stage_manifest.tsv")
}

figure7_upstream_quarantine <- function(path, root, label) {
  if (!file.exists(path) && !dir.exists(path)) return("")
  normalized_path <- normalizePath(path, mustWork = FALSE)
  normalized_root <- normalizePath(root, mustWork = FALSE)
  if (!startsWith(
    normalized_path,
    paste0(normalized_root, .Platform$file.sep)
  )) {
    figure7_stop("Refusing to quarantine ", label, " outside output root")
  }
  target <- paste0(
    path,
    ".stale.",
    format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"),
    ".",
    Sys.getpid()
  )
  if (!file.rename(path, target)) {
    figure7_stop("Cannot preserve stale ", label, ": ", path)
  }
  message("[seurat-upstream] Preserved stale ", label, " at ", target)
  target
}

figure7_upstream_manifest_values <- function(path) {
  if (!file.exists(path)) return(NULL)
  manifest <- tryCatch(
    figure7_read_tsv(path, c("key", "value"), "upstream stage manifest"),
    error = function(error) NULL
  )
  if (is.null(manifest) ||
      anyNA(manifest$key) ||
      any(!nzchar(manifest$key)) ||
      anyDuplicated(manifest$key)) {
    return(NULL)
  }
  stats::setNames(as.character(manifest$value), manifest$key)
}

figure7_upstream_source_config_contract <- function(config) {
  source_artifacts <- config$versioned_source_artifacts
  values <- c(
    endpoint_ploidy_default_path =
      as.character(source_artifacts$endpoint_ploidy$default_path),
    endpoint_ploidy_sha256 =
      as.character(source_artifacts$endpoint_ploidy$sha256),
    sample_info_default_path =
      as.character(source_artifacts$sample_info$default_path),
    sample_info_sha256 =
      as.character(source_artifacts$sample_info$sha256)
  )
  if (anyNA(values) ||
      any(!nzchar(values)) ||
      any(!grepl("^[0-9a-f]{64}$", values[c(
        "endpoint_ploidy_sha256",
        "sample_info_sha256"
      )]))) {
    figure7_stop("Seurat-upstream source config contract is incomplete")
  }
  figure7_contract_sha256(values)
}

figure7_upstream_common_contract <- function(
  environment_lock_path,
  config,
  jobs
) {
  c(
    seurat_upstream_environment_contract_sha256 =
      figure7_environment_stage_contract_sha256(
        environment_lock_path,
        "r",
        "seurat_upstream"
      ),
    source_config_contract_sha256 =
      figure7_upstream_source_config_contract(config),
    deterministic_workers = "1"
  )
}

figure7_upstream_stage_order <- function() {
  c("integrated", "cell_cycle", "refined", "merged", "final")
}

figure7_upstream_stage_function_names <- function(stage) {
  switch(
    stage,
    integrated = c(
      "figure7_read_tsv",
      "figure7_upstream_require_namespace",
      "figure7_upstream_configure_future",
      "figure7_upstream_read_sample_info",
      "figure7_upstream_resolve_column",
      "figure7_upstream_expected_samples",
      "figure7_upstream_h5_inventory",
      "figure7_upstream_read_counts",
      "figure7_upstream_qc_special",
      "figure7_upstream_remove_doublets",
      "figure7_upstream_build_integrated"
    ),
    cell_cycle = c(
      "figure7_upstream_require_namespace",
      "figure7_upstream_require_columns",
      "figure7_upstream_normalize_feature_key",
      "figure7_upstream_feature_metadata",
      "figure7_upstream_feature_lookup",
      "figure7_upstream_resolve_gene_set",
      "figure7_upstream_safe_sd",
      "figure7_upstream_cell_cycle_candidate_table",
      "figure7_upstream_annotate_cell_cycle"
    ),
    refined = c(
      "figure7_upstream_require_namespace",
      "figure7_upstream_require_columns",
      "figure7_upstream_sort_maybe_numeric",
      "figure7_upstream_set_single_thread",
      "figure7_upstream_build_refined_levels",
      "figure7_upstream_point_in_polygon",
      "figure7_upstream_knn_indices",
      "figure7_upstream_knn_components",
      "figure7_upstream_expand_candidate_component",
      "figure7_upstream_refine_cluster_values",
      "figure7_upstream_refine_clusters"
    ),
    merged = c(
      "figure7_upstream_require_namespace",
      "figure7_upstream_require_columns",
      "figure7_upstream_sort_maybe_numeric",
      "figure7_upstream_set_single_thread",
      "figure7_upstream_configure_future",
      "figure7_upstream_manual_merge_values",
      "figure7_upstream_join_layers",
      "figure7_upstream_assay_data",
      "figure7_upstream_merge_clusters"
    ),
    final = c(
      "figure7_upstream_require_namespace",
      "figure7_upstream_require_columns",
      "figure7_upstream_sort_maybe_numeric",
      "figure7_upstream_set_single_thread",
      "figure7_upstream_configure_future",
      "figure7_upstream_join_layers",
      "figure7_upstream_assay_data",
      "figure7_upstream_derive_tn_ploidy",
      "figure7_upstream_choose_pca_assay",
      "figure7_upstream_prepare_pca_assay",
      "figure7_upstream_finalize_clusters"
    ),
    figure7_stop("Unknown Seurat-upstream stage: ", stage)
  )
}

figure7_upstream_stage_legacy_sources <- function(stage) {
  roles <- switch(
    stage,
    integrated = c("01_data.R", "Utils.R"),
    cell_cycle = "01a_cell_cycle.R",
    refined = "02b_cluster_refine.R",
    merged = "02d_manual_cluster_merge.R",
    final = "03_final_cluster.R",
    figure7_stop("Unknown Seurat-upstream stage: ", stage)
  )
  figure7_upstream_legacy_sources[roles]
}

figure7_upstream_stage_parameter_contract <- function(stage) {
  definitions <- figure7_upstream_stage_definitions(tempdir())
  if (!stage %in% names(definitions)) {
    figure7_stop("Unknown Seurat-upstream stage: ", stage)
  }
  as.character(definitions[[stage]]$parameter_contract)
}

figure7_upstream_stage_code_contract <- function(stage) {
  functions <- figure7_upstream_stage_function_names(stage)
  bodies <- vapply(
    functions,
    function(name) {
      fn <- get(name, mode = "function", inherits = TRUE)
      paste(
        deparse(fn, width.cutoff = 500L),
        collapse = "\n"
      )
    },
    character(1L)
  )
  names(bodies) <- paste0("function:", functions)
  legacy <- figure7_upstream_stage_legacy_sources(stage)
  names(legacy) <- paste0("legacy_source:", names(legacy))
  figure7_contract_sha256(c(
    bodies,
    legacy,
    parameter_contract =
      figure7_upstream_stage_parameter_contract(stage)
  ))
}

figure7_upstream_cumulative_code_contracts <- function(stage) {
  order <- figure7_upstream_stage_order()
  index <- match(stage, order)
  if (is.na(index)) {
    figure7_stop("Unknown Seurat-upstream stage: ", stage)
  }
  stages <- order[seq_len(index)]
  stats::setNames(
    vapply(
      stages,
      figure7_upstream_stage_code_contract,
      character(1L)
    ),
    paste0("scientific_code_contract:", stages)
  )
}

figure7_upstream_stage_dependency_roles <- function(stage) {
  raw_roles <- c("all_ploidy", "sample_info", "h5_inventory")
  switch(
    stage,
    integrated = raw_roles,
    cell_cycle = c("integrated_rds", raw_roles),
    refined = c("cell_cycle_rds", "integrated_rds", raw_roles),
    merged = c(
      "refined_rds", "cell_cycle_rds", "integrated_rds", raw_roles
    ),
    final = c(
      "merged_rds", "refined_rds", "cell_cycle_rds",
      "integrated_rds", raw_roles
    ),
    figure7_stop("Unknown Seurat-upstream stage: ", stage)
  )
}

figure7_upstream_manifest_dependencies <- function(manifest) {
  keys <- grep(
    "^dependency_sha256:",
    names(manifest),
    value = TRUE
  )
  stats::setNames(
    unname(as.character(manifest[keys])),
    sub("^dependency_sha256:", "", keys)
  )
}

figure7_upstream_creator_runtime <- function() {
  runtime <- figure7_r_runtime_provenance()
  stats::setNames(
    unname(as.character(runtime)),
    sub("^audit_", "creator_", names(runtime))
  )
}

figure7_upstream_h5_inventory_digest <- function(samples, hashes) {
  samples <- as.character(samples)
  hashes <- as.character(hashes)
  if (length(samples) != 18L ||
      length(hashes) != 18L ||
      anyNA(samples) ||
      any(!nzchar(samples)) ||
      anyDuplicated(samples) ||
      any(!grepl("^[0-9a-f]{64}$", hashes))) {
    figure7_stop("Cell Ranger provenance requires 18 unique H5 hashes")
  }
  order <- order(samples)
  digest::digest(
    paste(
      paste(samples[order], hashes[order], sep = "="),
      collapse = "\n"
    ),
    algo = "sha256",
    serialize = FALSE
  )
}

figure7_upstream_manifest_matches <- function(
  stage,
  definition,
  common_contract,
  expected_dependencies = character()
) {
  output <- figure7_upstream_stage_output(definition)
  manifest <- figure7_upstream_manifest_values(
    figure7_upstream_stage_manifest(definition)
  )
  if (is.null(manifest) || !file.exists(output)) return(FALSE)
  expected <- c(
    schema_version = "2",
    stage = stage,
    output_filename = definition$filename,
    parameter_contract = definition$parameter_contract,
    common_contract,
    figure7_upstream_cumulative_code_contracts(stage)
  )
  if (!all(names(expected) %in% names(manifest)) ||
      !identical(
        unname(manifest[names(expected)]),
        unname(as.character(expected))
      )) {
    return(FALSE)
  }
  expected_code_keys <- grep(
    "^scientific_code_contract:",
    names(expected),
    value = TRUE
  )
  observed_code_keys <- grep(
    "^scientific_code_contract:",
    names(manifest),
    value = TRUE
  )
  if (!setequal(expected_code_keys, observed_code_keys)) {
    return(FALSE)
  }
  dependencies <- figure7_upstream_manifest_dependencies(manifest)
  required_roles <- figure7_upstream_stage_dependency_roles(stage)
  if (!setequal(names(dependencies), required_roles) ||
      any(!grepl("^[0-9a-f]{64}$", dependencies))) {
    return(FALSE)
  }
  if (length(expected_dependencies)) {
    if (is.null(names(expected_dependencies)) ||
        any(!names(expected_dependencies) %in% required_roles) ||
        anyDuplicated(names(expected_dependencies)) ||
        any(!grepl("^[0-9a-f]{64}$", expected_dependencies)) ||
        !identical(
          unname(dependencies[names(expected_dependencies)]),
          unname(as.character(expected_dependencies))
        )) {
      return(FALSE)
    }
  }
  creator_keys <- c(
    "creator_r_runtime_version",
    "creator_r_platform",
    "creator_r_blas"
  )
  if (!all(creator_keys %in% names(manifest)) ||
      anyNA(manifest[creator_keys]) ||
      any(!nzchar(manifest[creator_keys]))) {
    return(FALSE)
  }
  audit_code_keys <- c(
    "audit_helper_sha256",
    "audit_generator_sha256"
  )
  if (!all(audit_code_keys %in% names(manifest)) ||
      any(!grepl("^[0-9a-f]{64}$", manifest[audit_code_keys]))) {
    return(FALSE)
  }
  if (identical(stage, "integrated")) {
    h5_keys <- grep("^h5_sha256:", names(manifest), value = TRUE)
    h5_samples <- sub("^h5_sha256:", "", h5_keys)
    h5_hashes <- unname(manifest[h5_keys])
    if (length(h5_keys) != 18L ||
        !identical(manifest[["h5_file_count"]], "18") ||
        anyDuplicated(h5_samples) ||
        any(!grepl("^[0-9a-f]{64}$", h5_hashes))) {
      return(FALSE)
    }
    recorded_inventory <- tryCatch(
      figure7_upstream_h5_inventory_digest(h5_samples, h5_hashes),
      error = function(error) NA_character_
    )
    if (is.na(recorded_inventory) ||
        !identical(
          dependencies[["h5_inventory"]],
          recorded_inventory
        ) ||
        !identical(
          manifest[["h5_inventory_sha256"]],
          recorded_inventory
        )) {
      return(FALSE)
    }
  }
  recorded_hash <- manifest[["output_sha256"]]
  is.character(recorded_hash) &&
    length(recorded_hash) == 1L &&
    grepl("^[0-9a-f]{64}$", recorded_hash) &&
    identical(figure7_sha256(output), recorded_hash)
}

figure7_upstream_write_manifest <- function(
  stage,
  definition,
  common_contract,
  dependencies,
  extra = character()
) {
  if (is.null(names(dependencies)) ||
      any(!nzchar(names(dependencies))) ||
      anyDuplicated(names(dependencies)) ||
      any(!grepl("^[0-9a-f]{64}$", dependencies)) ||
      !setequal(
        names(dependencies),
        figure7_upstream_stage_dependency_roles(stage)
      )) {
    figure7_stop("Upstream stage dependencies require unique SHA-256 roles")
  }
  output <- figure7_upstream_stage_output(definition)
  values <- c(
    schema_version = "2",
    stage = stage,
    output_filename = definition$filename,
    output_sha256 = figure7_sha256(output),
    parameter_contract = definition$parameter_contract,
    common_contract,
    figure7_upstream_cumulative_code_contracts(stage),
    stats::setNames(
      as.character(dependencies),
      paste0("dependency_sha256:", names(dependencies))
    ),
    figure7_upstream_creator_runtime(),
    audit_helper_sha256 = figure7_sha256(helper_path),
    audit_generator_sha256 = figure7_sha256(script_path),
    extra
  )
  if (anyDuplicated(names(values))) {
    figure7_stop("Upstream stage manifest keys are not unique")
  }
  figure7_write_tsv(
    data.frame(
      key = names(values),
      value = unname(as.character(values)),
      stringsAsFactors = FALSE
    ),
    figure7_upstream_stage_manifest(definition)
  )
}

figure7_upstream_promote_stage <- function(
  stage,
  definition,
  output_root,
  common_contract,
  dependencies,
  object,
  extra = character()
) {
  if (dir.exists(definition$directory)) {
    figure7_upstream_quarantine(
      definition$directory,
      output_root,
      paste(stage, "cache")
    )
  }
  temporary <- paste0(
    definition$directory,
    ".tmp.",
    Sys.getpid()
  )
  if (dir.exists(temporary)) {
    figure7_upstream_quarantine(
      temporary,
      output_root,
      paste(stage, "temporary cache")
    )
  }
  dir.create(temporary, recursive = TRUE, showWarnings = FALSE)
  temporary_definition <- definition
  temporary_definition$directory <- temporary
  saveRDS(
    object,
    figure7_upstream_stage_output(temporary_definition),
    version = 3
  )
  figure7_upstream_write_manifest(
    stage,
    temporary_definition,
    common_contract,
    dependencies,
    extra
  )
  if (!file.rename(temporary, definition$directory)) {
    figure7_stop("Cannot atomically promote ", stage, " Seurat cache")
  }
  invisible(figure7_upstream_stage_output(definition))
}

figure7_upstream_validate_lock <- function(path) {
  figure7_validate_r_environment(path, "seurat_upstream")
}

figure7_upstream_validate_integrated <- function(object) {
  if (!inherits(object, "Seurat")) {
    figure7_stop("Integrated cache does not contain a Seurat object")
  }
  figure7_upstream_require_columns(
    object@meta.data,
    c(
      "seurat_clusters", "sample", "Dose", "IDs", "harvest",
      "barcode_raw", "nCount_RNA", "nFeature_RNA", "percent.mt",
      "scDblFinder.class", "scDblFinder.score"
    ),
    "Integrated Seurat metadata"
  )
  if (!all(c("pca", "umap") %in% names(object@reductions))) {
    figure7_stop("Integrated Seurat cache lacks PCA/UMAP")
  }
  figure7_upstream_assert_counts(object$seurat_clusters, "integrated")
  if (ncol(object) != 42884L) {
    figure7_stop("Integrated Seurat cache must contain 42884 cells")
  }
  invisible(TRUE)
}

figure7_upstream_validate_cell_cycle <- function(object) {
  figure7_upstream_validate_integrated(object)
  figure7_upstream_require_columns(
    object@meta.data,
    c(
      "S.Score", "G2M.Score", "Phase",
      "cluster_cell_cycle_annotation"
    ),
    "Cell-cycle Seurat metadata"
  )
  candidate_clusters <- unique(as.character(object$seurat_clusters)[
    object$cluster_cell_cycle_annotation == "cell_cycle_candidate"
  ])
  if (!setequal(candidate_clusters, c("6", "10", "11", "12"))) {
    figure7_stop(
      "Cell-cycle candidate clusters differ from reviewed values 6,10,11,12"
    )
  }
  invisible(TRUE)
}

figure7_upstream_validate_refined <- function(object) {
  figure7_upstream_validate_cell_cycle(object)
  figure7_upstream_require_columns(
    object@meta.data,
    "seurat_cluster_refine",
    "Refined Seurat metadata"
  )
  figure7_upstream_assert_counts(object$seurat_cluster_refine, "refined")
  invisible(TRUE)
}

figure7_upstream_validate_merged <- function(object) {
  figure7_upstream_validate_refined(object)
  figure7_upstream_require_columns(
    object@meta.data,
    "manual_merge_test",
    "Merged Seurat metadata"
  )
  figure7_upstream_assert_counts(object$manual_merge_test, "merged")
  invisible(TRUE)
}

figure7_upstream_read_valid_stage <- function(
  stage,
  definition,
  common_contract,
  expected_dependencies,
  validator,
  force_rebuild
) {
  if (isTRUE(force_rebuild) ||
      !figure7_upstream_manifest_matches(
        stage,
        definition,
        common_contract,
        expected_dependencies
      )) {
    return(NULL)
  }
  object <- tryCatch(
    readRDS(figure7_upstream_stage_output(definition)),
    error = function(error) NULL
  )
  if (is.null(object)) return(NULL)
  valid <- tryCatch(
    {
      validator(object)
      TRUE
    },
    error = function(error) {
      message(
        "[seurat-upstream] Rejecting ",
        stage,
        " cache: ",
        conditionMessage(error)
      )
      FALSE
    }
  )
  if (!valid) return(NULL)
  message("[seurat-upstream] Reusing validated ", stage, " cache")
  object
}

figure7_upstream_previous_hash <- function(
  stage,
  definitions
) {
  figure7_sha256(figure7_upstream_stage_output(definitions[[stage]]))
}

figure7_upstream_descendant_dependencies <- function(
  parent_stage,
  definitions
) {
  parent_definition <- definitions[[parent_stage]]
  parent_manifest <- figure7_upstream_manifest_values(
    figure7_upstream_stage_manifest(parent_definition)
  )
  if (is.null(parent_manifest)) {
    figure7_stop(
      "Cannot propagate dependencies from missing ",
      parent_stage,
      " stage manifest"
    )
  }
  direct_role <- paste0(parent_stage, "_rds")
  c(
    stats::setNames(
      figure7_upstream_previous_hash(parent_stage, definitions),
      direct_role
    ),
    figure7_upstream_manifest_dependencies(parent_manifest)
  )
}

figure7_upstream_available_dependencies <- function(
  stage,
  definitions,
  raw_dependencies
) {
  roles <- figure7_upstream_stage_dependency_roles(stage)
  available <- raw_dependencies[
    intersect(names(raw_dependencies), roles)
  ]
  upstream_stages <- c(
    integrated = "integrated_rds",
    cell_cycle = "cell_cycle_rds",
    refined = "refined_rds",
    merged = "merged_rds"
  )
  for (upstream_stage in names(upstream_stages)) {
    role <- unname(upstream_stages[[upstream_stage]])
    if (!role %in% roles) next
    path <- figure7_upstream_stage_output(
      definitions[[upstream_stage]]
    )
    manifest_path <- figure7_upstream_stage_manifest(
      definitions[[upstream_stage]]
    )
    output_exists <- file.exists(path)
    manifest_exists <- file.exists(manifest_path)
    if (xor(output_exists, manifest_exists)) {
      available[[role]] <- paste(rep("0", 64L), collapse = "")
    } else if (output_exists) {
      available[[role]] <- figure7_sha256(path)
    }
  }
  available[intersect(roles, names(available))]
}

figure7_upstream_dependency_validation_summary <- function(
  stage,
  expected_dependencies
) {
  roles <- figure7_upstream_stage_dependency_roles(stage)
  verified <- intersect(roles, names(expected_dependencies))
  attested <- setdiff(roles, verified)
  c(
    dependency_validation_mode = if (length(attested)) {
      "verified_available_and_manifest_attested_absent"
    } else {
      "verified_all_current_dependencies"
    },
    dependency_roles_total = as.character(length(roles)),
    dependency_roles_verified_current = as.character(length(verified)),
    dependency_roles_manifest_attested = as.character(length(attested)),
    manifest_attested_dependency_roles = paste(attested, collapse = ",")
  )
}

figure7_upstream_write_run_manifest <- function(
  output_root,
  final_path,
  common_contract,
  requested_jobs,
  dependency_validation,
  audit_provenance = character()
) {
  final_manifest <- figure7_upstream_manifest_values(
    file.path(dirname(final_path), "stage_manifest.tsv")
  )
  if (is.null(final_manifest)) {
    figure7_stop("Cannot write reconstruction manifest without final provenance")
  }
  creator_keys <- grep("^creator_", names(final_manifest), value = TRUE)
  scientific_code_keys <- grep(
    "^scientific_code_contract:",
    names(final_manifest),
    value = TRUE
  )
  final_dependencies <- figure7_upstream_manifest_dependencies(
    final_manifest
  )
  validation_runtime <- figure7_r_runtime_provenance()
  names(validation_runtime) <- sub(
    "^audit_",
    "validation_",
    names(validation_runtime)
  )
  values <- c(
    schema_version = "2",
    module = "figure7_si4_7_seurat_upstream",
    raw_boundary =
      "external_18_cellranger_filtered_feature_bc_matrix_h5",
    versioned_source_endpoint_ploidy = "Data/in-vivo/all_ploidy.tsv",
    versioned_source_sample_info = "Data/in-vivo/sample_info.xlsx",
    fastq_to_cellranger_available = "false",
    all_ploidy_generation_available = "false",
    final_output = file.path(
      "03_final_cluster",
      basename(final_path)
    ),
    final_sha256 = figure7_sha256(final_path),
    requested_jobs = as.character(requested_jobs),
    scientific_workers = "1",
    dependency_validation,
    common_contract,
    final_manifest[scientific_code_keys],
    audit_helper_sha256 = figure7_sha256(helper_path),
    audit_generator_sha256 = figure7_sha256(script_path),
    audit_provenance,
    final_manifest[creator_keys],
    validation_runtime,
    stats::setNames(
      final_dependencies,
      paste0("final_dependency_sha256:", names(final_dependencies))
    ),
    stats::setNames(
      figure7_upstream_legacy_sources,
      paste0("legacy_source_sha256:", names(
        figure7_upstream_legacy_sources
      ))
    )
  )
  figure7_write_tsv(
    data.frame(
      key = names(values),
      value = unname(as.character(values)),
      stringsAsFactors = FALSE
    ),
    file.path(output_root, "reconstruction_manifest.tsv")
  )
}

figure7_upstream_validate_artifact <- function(
  output_root,
  environment_lock,
  config,
  all_ploidy,
  sample_info,
  expected_rds = "",
  cellranger_root = ""
) {
  output_root <- normalizePath(output_root, mustWork = TRUE)
  environment_lock <- normalizePath(environment_lock, mustWork = TRUE)
  all_ploidy <- normalizePath(all_ploidy, mustWork = TRUE)
  sample_info <- normalizePath(sample_info, mustWork = TRUE)
  figure7_verify_checksum(
    all_ploidy,
    as.character(
      config$versioned_source_artifacts$endpoint_ploidy$sha256
    ),
    "Seurat-upstream endpoint-ploidy source"
  )
  figure7_verify_checksum(
    sample_info,
    as.character(
      config$versioned_source_artifacts$sample_info$sha256
    ),
    "Seurat-upstream sample-info source"
  )
  definitions <- figure7_upstream_stage_definitions(output_root)
  final_path <- figure7_upstream_stage_output(definitions$final)
  if (!file.exists(final_path)) {
    figure7_stop("Seurat-upstream final RDS is missing: ", final_path)
  }
  raw_dependencies <- c(
    all_ploidy = figure7_sha256(all_ploidy),
    sample_info = figure7_sha256(sample_info)
  )
  if (nzchar(cellranger_root)) {
    cellranger_root <- normalizePath(cellranger_root, mustWork = TRUE)
    inventory <- figure7_upstream_h5_inventory(
      cellranger_root,
      figure7_upstream_read_sample_info(sample_info)
    )
    raw_dependencies[["h5_inventory"]] <-
      figure7_upstream_h5_inventory_digest(
        inventory$sample,
        inventory$sha256
      )
  }
  common_contract <- figure7_upstream_common_contract(
    environment_lock,
    config,
    jobs = 1L
  )
  expected_dependencies <- figure7_upstream_available_dependencies(
    "final",
    definitions,
    raw_dependencies
  )
  if (!figure7_upstream_manifest_matches(
    "final",
    definitions$final,
    common_contract,
    expected_dependencies
  )) {
    figure7_stop(
      "Seurat-upstream final stage fails its scientific dependency contract"
    )
  }
  stage_manifest <- figure7_upstream_manifest_values(
    figure7_upstream_stage_manifest(definitions$final)
  )
  reconstruction_path <- file.path(
    output_root,
    "reconstruction_manifest.tsv"
  )
  reconstruction <- figure7_upstream_manifest_values(reconstruction_path)
  if (is.null(reconstruction)) {
    figure7_stop("Seurat-upstream reconstruction manifest is missing or invalid")
  }
  final_hash <- figure7_sha256(final_path)
  required <- c(
    schema_version = "2",
    module = "figure7_si4_7_seurat_upstream",
    final_output = file.path(
      "03_final_cluster",
      basename(final_path)
    ),
    final_sha256 = final_hash,
    scientific_workers = "1",
    common_contract,
    figure7_upstream_cumulative_code_contracts("final")
  )
  if (!all(names(required) %in% names(reconstruction)) ||
      !identical(
        unname(reconstruction[names(required)]),
        unname(as.character(required))
      )) {
    figure7_stop(
      "Seurat-upstream reconstruction manifest contradicts the final stage"
    )
  }
  reconstruction_code_keys <- grep(
    "^scientific_code_contract:",
    names(reconstruction),
    value = TRUE
  )
  if (!setequal(
    reconstruction_code_keys,
    names(figure7_upstream_cumulative_code_contracts("final"))
  )) {
    figure7_stop(
      "Seurat-upstream reconstruction manifest has an incomplete ",
      "stage-code chain"
    )
  }
  audit_code_keys <- c(
    "audit_helper_sha256",
    "audit_generator_sha256"
  )
  if (!all(audit_code_keys %in% names(reconstruction)) ||
      any(!grepl(
        "^[0-9a-f]{64}$",
        reconstruction[audit_code_keys]
      ))) {
    figure7_stop(
      "Seurat-upstream reconstruction manifest lacks code audit hashes"
    )
  }
  creator_keys <- grep("^creator_", names(stage_manifest), value = TRUE)
  if (!length(creator_keys) ||
      !all(creator_keys %in% names(reconstruction)) ||
      !identical(
        unname(reconstruction[creator_keys]),
        unname(stage_manifest[creator_keys])
      )) {
    figure7_stop("Seurat-upstream creator provenance is inconsistent")
  }
  stage_dependencies <- figure7_upstream_manifest_dependencies(stage_manifest)
  recorded_dependency_keys <- paste0(
    "final_dependency_sha256:",
    names(stage_dependencies)
  )
  if (!all(recorded_dependency_keys %in% names(reconstruction)) ||
      !identical(
        unname(reconstruction[recorded_dependency_keys]),
        unname(stage_dependencies)
      )) {
    figure7_stop(
      "Seurat-upstream reconstruction manifest loses transitive dependencies"
    )
  }
  legacy_keys <- paste0(
    "legacy_source_sha256:",
    names(figure7_upstream_legacy_sources)
  )
  if (!all(legacy_keys %in% names(reconstruction)) ||
      !identical(
        unname(reconstruction[legacy_keys]),
        unname(figure7_upstream_legacy_sources)
      )) {
    figure7_stop("Seurat-upstream reviewed source provenance is inconsistent")
  }
  if (nzchar(expected_rds)) {
    expected_rds <- normalizePath(expected_rds, mustWork = TRUE)
    expected_hash <- if (identical(expected_rds, final_path)) {
      final_hash
    } else {
      figure7_sha256(expected_rds)
    }
    if (!identical(expected_hash, final_hash)) {
      figure7_stop(
        "Explicit non-deposited RDS is not bound to the reconstruction manifest"
      )
    }
  }
  list(
    valid = TRUE,
    output_root = output_root,
    rds = final_path,
    rds_sha256 = final_hash,
    final_stage_manifest =
      figure7_upstream_stage_manifest(definitions$final),
    final_stage_manifest_sha256 = figure7_sha256(
      figure7_upstream_stage_manifest(definitions$final)
    ),
    reconstruction_manifest = reconstruction_path,
    reconstruction_manifest_sha256 =
      figure7_sha256(reconstruction_path),
    dependencies = stage_dependencies,
    scientific_code_contracts =
      figure7_upstream_cumulative_code_contracts("final"),
    upstream_common_contract = common_contract,
    expected_dependencies = expected_dependencies
  )
}

figure7_upstream_validate_resumable_artifact <- function(
  output_root,
  environment_lock,
  config,
  all_ploidy,
  sample_info,
  cellranger_root = ""
) {
  output_root <- normalizePath(output_root, mustWork = TRUE)
  environment_lock <- normalizePath(environment_lock, mustWork = TRUE)
  all_ploidy <- normalizePath(all_ploidy, mustWork = TRUE)
  sample_info <- normalizePath(sample_info, mustWork = TRUE)
  figure7_verify_checksum(
    all_ploidy,
    as.character(
      config$versioned_source_artifacts$endpoint_ploidy$sha256
    ),
    "Seurat-upstream endpoint-ploidy source"
  )
  figure7_verify_checksum(
    sample_info,
    as.character(
      config$versioned_source_artifacts$sample_info$sha256
    ),
    "Seurat-upstream sample-info source"
  )
  definitions <- figure7_upstream_stage_definitions(output_root)
  raw_dependencies <- c(
    all_ploidy = figure7_sha256(all_ploidy),
    sample_info = figure7_sha256(sample_info)
  )
  if (nzchar(cellranger_root)) {
    cellranger_root <- normalizePath(cellranger_root, mustWork = TRUE)
    inventory <- figure7_upstream_h5_inventory(
      cellranger_root,
      figure7_upstream_read_sample_info(sample_info)
    )
    raw_dependencies[["h5_inventory"]] <-
      figure7_upstream_h5_inventory_digest(
        inventory$sample,
        inventory$sha256
      )
  }
  common_contract <- figure7_upstream_common_contract(
    environment_lock,
    config,
    jobs = 1L
  )
  for (stage in c("merged", "refined", "cell_cycle", "integrated")) {
    expected <- figure7_upstream_available_dependencies(
      stage,
      definitions,
      raw_dependencies
    )
    if (figure7_upstream_manifest_matches(
      stage,
      definitions[[stage]],
      common_contract,
      expected
    )) {
      manifest <- figure7_upstream_manifest_values(
        figure7_upstream_stage_manifest(definitions[[stage]])
      )
      return(list(
        valid = TRUE,
        output_root = output_root,
        stage = stage,
        stage_rds = figure7_upstream_stage_output(
          definitions[[stage]]
        ),
        stage_rds_sha256 = unname(manifest[["output_sha256"]]),
        stage_manifest = figure7_upstream_stage_manifest(
          definitions[[stage]]
        ),
        stage_manifest_sha256 = figure7_sha256(
          figure7_upstream_stage_manifest(definitions[[stage]])
        ),
        dependencies =
          figure7_upstream_manifest_dependencies(manifest),
        expected_dependencies = expected
      ))
    }
  }
  figure7_stop(
    "No scientifically manifested resumable Seurat-upstream stage was found"
  )
}

figure7_upstream_probe_resumable_artifact <- function(...) {
  tryCatch(
    figure7_upstream_validate_resumable_artifact(...),
    error = function(error) {
      list(valid = FALSE, error = conditionMessage(error))
    }
  )
}

main <- function(
  args = figure7_parse_args(commandArgs(trailingOnly = TRUE))
) {
  figure7_upstream_require_namespace("digest")
  output_root <- figure7_upstream_resolve_path(
    figure7_arg(args, "output-dir", required = TRUE),
    must_work = FALSE
  )
  if (identical(output_root, repo_root)) {
    figure7_stop("Seurat reconstruction output cannot be the repository root")
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_root)) {
    figure7_stop("Cannot create Seurat reconstruction output: ", output_root)
  }
  environment_lock <- figure7_upstream_resolve_path(
    figure7_arg(
      args,
      "environment-lock",
      "Code/in-vivo/figure7/environment_lock.tsv"
    ),
    must_work = TRUE
  )
  config_path <- figure7_upstream_resolve_path(
    figure7_arg(
      args,
      "config",
      "Code/in-vivo/figure7/figure7_config.yaml"
    ),
    must_work = TRUE
  )
  config <- figure7_read_config(config_path)
  jobs <- suppressWarnings(as.integer(
    figure7_arg(args, "jobs", "1")
  ))
  if (length(jobs) != 1L || is.na(jobs) || jobs < 1L || jobs > 16L) {
    figure7_stop("--jobs must be an integer from 1 to 16")
  }
  if (jobs != 1L) {
    message(
      "[seurat-upstream] --jobs=",
      jobs,
      " does not alter Tao's reviewed single-worker Seurat method"
    )
  }
  force_rebuild <- figure7_upstream_bool(
    figure7_arg(args, "force-rebuild", "false"),
    "--force-rebuild"
  )
  definitions <- figure7_upstream_stage_definitions(output_root)
  common_contract <- figure7_upstream_common_contract(
    environment_lock,
    config,
    jobs
  )
  audit_provenance <- c(
    audit_current_common_io_sha256 = figure7_sha256(common_path),
    audit_current_environment_lock_sha256 =
      figure7_sha256(environment_lock),
    audit_current_figure7_config_sha256 = figure7_sha256(config_path)
  )
  all_ploidy <- figure7_upstream_resolve_path(
    figure7_arg(
      args,
      "all-ploidy",
      as.character(
        config$versioned_source_artifacts$endpoint_ploidy$default_path
      )
    ),
    must_work = TRUE
  )
  sample_info <- figure7_upstream_resolve_path(
    figure7_arg(
      args,
      "sample-info",
      as.character(
        config$versioned_source_artifacts$sample_info$default_path
      )
    ),
    must_work = TRUE
  )
  figure7_verify_checksum(
    all_ploidy,
    as.character(
      config$versioned_source_artifacts$endpoint_ploidy$sha256
    ),
    "versioned endpoint-ploidy source"
  )
  figure7_verify_checksum(
    sample_info,
    as.character(
      config$versioned_source_artifacts$sample_info$sha256
    ),
    "versioned sample-info source"
  )
  raw_dependencies <- c(
    all_ploidy = figure7_sha256(all_ploidy),
    sample_info = figure7_sha256(sample_info)
  )
  cellranger_root_argument <- figure7_arg(
    args,
    "cellranger-root",
    ""
  )
  cellranger_root <- NULL
  current_h5_inventory <- NULL
  if (nzchar(cellranger_root_argument)) {
    candidate_cellranger_root <- figure7_upstream_resolve_path(
      cellranger_root_argument,
      must_work = FALSE
    )
    if (dir.exists(candidate_cellranger_root)) {
      cellranger_root <- normalizePath(
        candidate_cellranger_root,
        mustWork = TRUE
      )
      current_h5_inventory <- figure7_upstream_h5_inventory(
        cellranger_root,
        figure7_upstream_read_sample_info(sample_info)
      )
      raw_dependencies[["h5_inventory"]] <-
        figure7_upstream_h5_inventory_digest(
          current_h5_inventory$sample,
          current_h5_inventory$sha256
        )
    } else {
      message(
        "[seurat-upstream] Cell Ranger root is currently unavailable; ",
        "manifested downstream stages will be tried before it is required"
      )
    }
  }

  final_expected <- figure7_upstream_available_dependencies(
    "final",
    definitions,
    raw_dependencies
  )
  final <- figure7_upstream_read_valid_stage(
    "final",
    definitions$final,
    common_contract,
    final_expected,
    function(object) figure7_upstream_validate_final(object, TRUE),
    force_rebuild
  )
  if (!is.null(final)) {
    final_path <- figure7_upstream_stage_output(definitions$final)
    figure7_upstream_write_run_manifest(
      output_root,
      final_path,
      common_contract,
      jobs,
      figure7_upstream_dependency_validation_summary(
        "final",
        final_expected
      ),
      audit_provenance
    )
    cat(final_path, "\n", sep = "")
    return(invisible(final_path))
  }

  merged_expected <- figure7_upstream_available_dependencies(
    "merged",
    definitions,
    raw_dependencies
  )
  merged <- figure7_upstream_read_valid_stage(
    "merged",
    definitions$merged,
    common_contract,
    merged_expected,
    figure7_upstream_validate_merged,
    force_rebuild
  )
  if (is.null(merged)) {
    refined_expected <- figure7_upstream_available_dependencies(
      "refined",
      definitions,
      raw_dependencies
    )
    refined <- figure7_upstream_read_valid_stage(
      "refined",
      definitions$refined,
      common_contract,
      refined_expected,
      figure7_upstream_validate_refined,
      force_rebuild
    )
    if (is.null(refined)) {
      cell_cycle_expected <- figure7_upstream_available_dependencies(
        "cell_cycle",
        definitions,
        raw_dependencies
      )
      annotated <- figure7_upstream_read_valid_stage(
        "cell_cycle",
        definitions$cell_cycle,
        common_contract,
        cell_cycle_expected,
        figure7_upstream_validate_cell_cycle,
        force_rebuild
      )
      if (is.null(annotated)) {
        integrated_expected <- figure7_upstream_available_dependencies(
          "integrated",
          definitions,
          raw_dependencies
        )
        integrated <- figure7_upstream_read_valid_stage(
          "integrated",
          definitions$integrated,
          common_contract,
          integrated_expected,
          figure7_upstream_validate_integrated,
          force_rebuild
        )
        if (is.null(integrated)) {
          if (is.null(cellranger_root)) {
            if (nzchar(cellranger_root_argument)) {
              figure7_stop(
                "Cell Ranger root is required to rebuild the missing ",
                "integrated stage but does not exist: ",
                candidate_cellranger_root
              )
            }
            figure7_stop("Missing required argument --cellranger-root")
          }
          figure7_upstream_validate_lock(environment_lock)
          built <- figure7_upstream_build_integrated(
            cellranger_root,
            all_ploidy,
            sample_info,
            jobs
          )
          integrated <- built$object
          figure7_upstream_validate_integrated(integrated)
          inventory_digest <- figure7_upstream_h5_inventory_digest(
            built$h5_inventory$sample,
            built$h5_inventory$sha256
          )
          extra <- c(
            raw_boundary =
              "external_18_cellranger_filtered_feature_bc_matrix_h5",
            h5_inventory_sha256 = inventory_digest,
            h5_file_count = as.character(nrow(built$h5_inventory)),
            stats::setNames(
              built$h5_inventory$sha256,
              paste0("h5_sha256:", built$h5_inventory$sample)
            )
          )
          dependencies <- c(
            all_ploidy = figure7_sha256(all_ploidy),
            sample_info = figure7_sha256(sample_info),
            h5_inventory = inventory_digest
          )
          figure7_upstream_promote_stage(
            "integrated",
            definitions$integrated,
            output_root,
            common_contract,
            dependencies,
            integrated,
            extra
          )
          figure7_write_tsv(
            built$sample_summary,
            file.path(
              definitions$integrated$directory,
              "sample_filter_summary.tsv"
            )
          )
        }
        figure7_upstream_validate_lock(environment_lock)
        annotation <- figure7_upstream_annotate_cell_cycle(integrated)
        annotated <- annotation$object
        figure7_upstream_validate_cell_cycle(annotated)
        figure7_upstream_promote_stage(
          "cell_cycle",
          definitions$cell_cycle,
          output_root,
          common_contract,
          figure7_upstream_descendant_dependencies(
            "integrated",
            definitions
          ),
          annotated,
          c(
            candidate_clusters = paste(
              annotation$candidate_table$cluster[
                annotation$candidate_table$final_candidate_flag
              ],
              collapse = ","
            )
          )
        )
      }
      figure7_upstream_validate_lock(environment_lock)
      refinement <- figure7_upstream_refine_clusters(annotated)
      refined <- refinement$object
      figure7_upstream_validate_refined(refined)
      figure7_upstream_promote_stage(
        "refined",
        definitions$refined,
        output_root,
        common_contract,
        figure7_upstream_descendant_dependencies(
          "cell_cycle",
          definitions
        ),
        refined,
        c(
          selected_4c = as.character(sum(
            as.character(refined$seurat_cluster_refine) == "4c"
          )),
          selected_9c = as.character(sum(
            as.character(refined$seurat_cluster_refine) == "9c"
          ))
        )
      )
    }
    figure7_upstream_validate_lock(environment_lock)
    merged <- figure7_upstream_merge_clusters(refined, jobs)
    figure7_upstream_validate_merged(merged)
    figure7_upstream_promote_stage(
      "merged",
      definitions$merged,
      output_root,
      common_contract,
      figure7_upstream_descendant_dependencies(
        "refined",
        definitions
      ),
      merged
    )
  }
  figure7_upstream_validate_lock(environment_lock)
  set.seed(1234)
  final <- figure7_upstream_finalize_clusters(
    merged,
    jobs,
    rerun_reductions = TRUE
  )
  figure7_upstream_validate_final(final, TRUE)
  figure7_upstream_promote_stage(
    "final",
    definitions$final,
    output_root,
    common_contract,
    figure7_upstream_descendant_dependencies(
      "merged",
      definitions
    ),
    final
  )
  final_path <- figure7_upstream_stage_output(definitions$final)
  figure7_upstream_write_run_manifest(
    output_root,
    final_path,
    common_contract,
    jobs,
    figure7_upstream_dependency_validation_summary(
      "final",
      figure7_upstream_available_dependencies(
        "final",
        definitions,
        raw_dependencies
      )
    ),
    audit_provenance
  )
  cat(final_path, "\n", sep = "")
  invisible(final_path)
}

if (sys.nframe() == 0L) main()
