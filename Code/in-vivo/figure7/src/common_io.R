# Shared I/O and output-contract helpers for the narrow Figure 7 module.

figure7_stop <- function(...) stop(..., call. = FALSE)

figure7_parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) figure7_stop("Unexpected positional argument: ", arg)
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L || !nzchar(pieces[[1L]])) {
      figure7_stop("Arguments must use --name=value syntax: ", arg)
    }
    out[[gsub("-", "_", pieces[[1L]])]] <- pieces[[2L]]
  }
  out
}

figure7_arg <- function(args, name, default = NULL, required = FALSE) {
  value <- args[[gsub("-", "_", name)]]
  if (is.null(value) || !nzchar(value)) value <- default
  if (isTRUE(required) && (is.null(value) || !nzchar(value))) {
    figure7_stop("Missing required argument --", gsub("_", "-", name))
  }
  value
}

figure7_sha256 <- function(path) {
  if (!file.exists(path)) figure7_stop("Cannot checksum missing file: ", path)
  if (!requireNamespace("digest", quietly = TRUE)) figure7_stop("R package 'digest' is required")
  unname(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

figure7_required_environment_pins <- function() {
  python_versions <- c(
    python = "3.10.13",
    numpy = "1.26.4",
    pandas = "2.2.3",
    anndata = "0.10.9",
    scanpy = "1.10.4",
    scvelo = "0.3.4",
    scipy = "1.13.1",
    `scikit-learn` = "1.5.2",
    numba = "0.60.0",
    loompy = "3.0.8",
    h5py = "3.11.0",
    igraph = "0.11.8",
    leidenalg = "0.10.2",
    `umap-learn` = "0.5.6",
    pynndescent = "0.6.0",
    llvmlite = "0.43.0",
    joblib = "1.5.3",
    threadpoolctl = "3.6.0",
    statsmodels = "0.14.6",
    networkx = "3.4.2",
    `numpy-groupies` = "0.11.3",
    matplotlib = "3.8.4"
  )
  r_versions <- c(
    R = "4.5.0",
    digest = "0.6.39",
    future = "1.70.0",
    `future.apply` = "1.20.2",
    readr = "2.2.0",
    readxl = "1.4.5",
    yaml = "2.3.12",
    Seurat = "4.4.0",
    SeuratObject = "4.1.4",
    Matrix = "1.7-3",
    scDblFinder = "1.22.0",
    SingleCellExperiment = "1.32.0",
    hdf5r = "1.3.12",
    sctransform = "0.4.3",
    RANN = "2.6.2",
    edgeR = "4.8.2",
    limma = "3.66.0",
    fgsea = "1.36.2",
    msigdbr = "26.1.0"
  )
  r_stage_packages <- list(
    scvelo = c(
      "R", "digest", "yaml", "Seurat", "SeuratObject", "readr"
    ),
    download = c("R", "digest"),
    celllevel = c("R", "digest", "readr", "readxl"),
    state = c(
      "R", "digest", "yaml", "Seurat", "SeuratObject", "Matrix",
      "edgeR", "limma", "fgsea", "msigdbr", "readr"
    ),
    state_export = c("R", "digest", "readr"),
    si_raw = c(
      "R", "digest", "yaml", "Seurat", "SeuratObject", "Matrix",
      "future", "future.apply", "fgsea", "msigdbr"
    ),
    seurat_upstream = c(
      "R", "digest", "future", "scDblFinder",
      "SingleCellExperiment", "Seurat", "SeuratObject", "Matrix",
      "hdf5r", "sctransform", "RANN", "readxl", "yaml"
    )
  )
  python <- data.frame(
    ecosystem = "python",
    stage = "scvelo",
    package = names(python_versions),
    version = unname(python_versions),
    stringsAsFactors = FALSE
  )
  r <- do.call(rbind, lapply(names(r_stage_packages), function(stage) {
    packages <- r_stage_packages[[stage]]
    data.frame(
      ecosystem = "r",
      stage = stage,
      package = packages,
      version = unname(r_versions[packages]),
      stringsAsFactors = FALSE
    )
  }))
  if (anyNA(r$version)) {
    figure7_stop("Internal Figure 7 R pin inventory is incomplete")
  }
  rbind(python, r)
}

figure7_environment_pin_key <- function(data) {
  paste(data$ecosystem, data$stage, data$package, sep = "\r")
}

figure7_read_environment_lock <- function(path) {
  if (!file.exists(path)) {
    figure7_stop("Missing Figure 7 environment lock: ", path)
  }
  lock <- utils::read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    colClasses = "character"
  )
  required_columns <- c("ecosystem", "stage", "package", "version")
  if (!identical(names(lock), required_columns) || !nrow(lock) ||
      anyNA(lock) ||
      any(!nzchar(trimws(as.matrix(lock))))) {
    figure7_stop(
      "Figure 7 environment lock must contain exactly nonempty columns: ",
      paste(required_columns, collapse = ", ")
    )
  }
  lock[] <- lapply(lock, trimws)
  if (anyDuplicated(figure7_environment_pin_key(lock))) {
    figure7_stop("Figure 7 environment lock contains duplicate pins")
  }

  required <- figure7_required_environment_pins()
  observed_key <- figure7_environment_pin_key(lock)
  required_key <- figure7_environment_pin_key(required)
  index <- match(required_key, observed_key)
  if (anyNA(index)) {
    missing <- paste(
      required$ecosystem[is.na(index)],
      required$stage[is.na(index)],
      required$package[is.na(index)],
      sep = "/"
    )
    figure7_stop(
      "Figure 7 environment lock is missing required exact pin(s): ",
      paste(missing, collapse = ", ")
    )
  }
  observed_version <- lock$version[index]
  mismatch <- observed_version != required$version
  if (any(mismatch)) {
    detail <- paste0(
      required$ecosystem[mismatch], "/",
      required$stage[mismatch], "/",
      required$package[mismatch], "=",
      observed_version[mismatch], " (expected ",
      required$version[mismatch], ")"
    )
    figure7_stop(
      "Figure 7 environment lock changes a reviewed exact pin: ",
      paste(detail, collapse = ", ")
    )
  }
  lock
}

figure7_normalize_r_version <- function(value) {
  tryCatch(
    as.character(utils::package_version(as.character(value))),
    error = function(error) as.character(value)
  )
}

figure7_validate_r_environment <- function(
  lock_path,
  stage,
  packages = NULL
) {
  lock <- figure7_read_environment_lock(lock_path)
  stage <- as.character(stage)
  if (length(stage) != 1L || !nzchar(stage)) {
    figure7_stop("R environment validation requires one explicit stage")
  }
  r_lock <- lock[
    lock$ecosystem == "r" & lock$stage == stage,
    ,
    drop = FALSE
  ]
  if (!nrow(r_lock)) {
    figure7_stop("Environment lock has no R pins for stage: ", stage)
  }
  if (is.null(packages)) {
    packages <- r_lock$package
  } else {
    packages <- unique(as.character(packages))
    packages <- packages[nzchar(packages)]
  }
  unpinned <- setdiff(packages, r_lock$package)
  if (length(unpinned)) {
    figure7_stop(
      "Environment lock stage ",
      stage,
      " does not explicitly pin R package(s): ",
      paste(unpinned, collapse = ", ")
    )
  }
  expected <- vapply(
    packages,
    function(package) {
      rows <- r_lock[r_lock$package == package, , drop = FALSE]
      if (nrow(rows) != 1L) {
        figure7_stop(
          "Environment lock must contain one stage-specific R pin for ",
          stage,
          "/",
          package
        )
      }
      rows$version[[1L]]
    },
    character(1L)
  )
  observed <- vapply(
    packages,
    function(package) {
      if (identical(package, "R")) {
        paste(R.version$major, sub("[[:space:]].*$", "", R.version$minor),
              sep = ".")
      } else {
        if (!requireNamespace(package, quietly = TRUE)) {
          figure7_stop(
            "Pinned R package is unavailable for raw computation: ",
            package
          )
        }
        as.character(utils::packageVersion(package))
      }
    },
    character(1L)
  )
  mismatch <- vapply(
    packages,
    function(package) {
      !identical(
        figure7_normalize_r_version(observed[[package]]),
        figure7_normalize_r_version(expected[[package]])
      )
    },
    logical(1L)
  )
  if (any(mismatch)) {
    detail <- paste0(
      packages[mismatch], "=", observed[mismatch],
      " (expected ", expected[mismatch], ")"
    )
    figure7_stop(
      "R runtime does not match the exact raw-computation lock: ",
      paste(detail, collapse = ", ")
    )
  }
  stats::setNames(observed, packages)
}

figure7_environment_stage_contract <- function(
  lock_path,
  ecosystem,
  stage
) {
  lock <- figure7_read_environment_lock(lock_path)
  pins <- lock[
    lock$ecosystem == ecosystem & lock$stage == stage,
    ,
    drop = FALSE
  ]
  if (!nrow(pins)) {
    figure7_stop(
      "Environment lock has no ",
      ecosystem,
      " pins for stage: ",
      stage
    )
  }
  pins <- pins[order(pins$package), , drop = FALSE]
  stats::setNames(
    as.character(pins$version),
    paste(pins$ecosystem, pins$stage, pins$package, sep = "/")
  )
}

figure7_contract_sha256 <- function(values) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    figure7_stop("R package 'digest' is required")
  }
  value_names <- names(values)
  values <- as.character(values)
  if (!length(values) || is.null(value_names) ||
      anyNA(value_names) || any(!nzchar(value_names)) ||
      anyDuplicated(value_names) || anyNA(values)) {
    figure7_stop("Contract values require unique names and nonmissing values")
  }
  names(values) <- value_names
  values <- values[order(names(values))]
  unname(digest::digest(
    paste(names(values), values, sep = "=", collapse = "\n"),
    algo = "sha256",
    serialize = FALSE
  ))
}

figure7_environment_stage_contract_sha256 <- function(
  lock_path,
  ecosystem,
  stage
) {
  figure7_contract_sha256(
    figure7_environment_stage_contract(lock_path, ecosystem, stage)
  )
}

figure7_validate_python_environment <- function(
  lock_path,
  python,
  stage = "scvelo"
) {
  lock <- figure7_read_environment_lock(lock_path)
  pins <- lock[
    lock$ecosystem == "python" & lock$stage == stage,
    ,
    drop = FALSE
  ]
  if (!nrow(pins)) {
    figure7_stop("Environment lock has no Python pins for stage: ", stage)
  }
  if (!nzchar(python) || !file.exists(python)) {
    figure7_stop("Missing pinned Python executable for raw computation")
  }
  distributions <- pins$package[pins$package != "python"]
  code <- paste0(
    "import importlib.metadata, platform\n",
    "print('python=' + platform.python_version())\n",
    paste0(
      "print(",
      vapply(distributions, function(package) {
        sprintf(
          "%s + '=' + importlib.metadata.version(%s)",
          dQuote(package),
          dQuote(package)
        )
      }, character(1L)),
      ")",
      collapse = "\n"
    )
  )
  output <- suppressWarnings(system2(
    python,
    c("-c", shQuote(code)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (!is.null(status) && as.integer(status) != 0L) {
    figure7_stop(
      "Pinned Python dependency preflight failed: ",
      paste(output, collapse = "\n")
    )
  }
  pieces <- strsplit(output[grepl("=", output, fixed = TRUE)], "=", fixed = TRUE)
  valid <- lengths(pieces) == 2L
  if (!all(valid)) {
    figure7_stop("Pinned Python dependency preflight returned malformed output")
  }
  observed <- stats::setNames(
    vapply(pieces, `[[`, character(1L), 2L),
    vapply(pieces, `[[`, character(1L), 1L)
  )
  expected <- stats::setNames(pins$version, pins$package)
  if (!all(names(expected) %in% names(observed)) ||
      !identical(
        unname(observed[names(expected)]),
        unname(expected)
      )) {
    mismatched <- names(expected)[
      is.na(observed[names(expected)]) |
        observed[names(expected)] != expected
    ]
    detail <- paste0(
      mismatched, "=",
      observed[mismatched],
      " (expected ", expected[mismatched], ")"
    )
    figure7_stop(
      "Python runtime does not match the exact raw-computation lock: ",
      paste(detail, collapse = ", ")
    )
  }
  observed[names(expected)]
}

figure7_r_runtime_provenance <- function() {
  session <- utils::sessionInfo()
  blas <- as.character(session$BLAS)
  if (!length(blas) || is.na(blas) || !nzchar(blas)) {
    blas <- "unavailable"
  } else {
    blas <- basename(blas)
  }
  c(
    audit_r_runtime_version = paste(
      R.version$major,
      sub("[[:space:]].*$", "", R.version$minor),
      sep = "."
    ),
    audit_r_platform = as.character(R.version$platform),
    audit_r_blas = blas
  )
}

figure7_feature_species_contract_values <- function(config) {
  policy <- config$feature_species
  required <- c(
    "policy_id", "human_prefix", "mouse_prefix",
    "unknown_feature_policy"
  )
  missing <- setdiff(required, names(policy))
  if (length(missing)) {
    figure7_stop(
      "Feature-species config is missing field(s): ",
      paste(missing, collapse = ", ")
    )
  }
  values <- c(
    policy_id = as.character(policy$policy_id),
    human_prefix = as.character(policy$human_prefix),
    mouse_prefix = as.character(policy$mouse_prefix),
    unknown_feature_policy =
      as.character(policy$unknown_feature_policy)
  )
  expected <- c(
    policy_id = "grch_human_tumor_only_v2",
    human_prefix = "GRCh38-",
    mouse_prefix = "GRCm39-",
    unknown_feature_policy = "reject"
  )
  if (!identical(values, expected)) {
    figure7_stop(
      "Figure 7/SI7 feature-species contract must use exact GRCh38-/",
      "GRCm39- prefixes and reject unrecognized features"
    )
  }
  values
}

figure7_read_config <- function(path, tgi_day = NULL) {
  if (!file.exists(path)) figure7_stop("Missing Figure 7 config: ", path)
  if (!requireNamespace("yaml", quietly = TRUE)) figure7_stop("R package 'yaml' is required")
  config <- yaml::read_yaml(path)
  required <- c(
    "schema_version", "module", "raw_data", "versioned_source_artifacts",
    "feature_species", "inputs",
    "si_figures", "tgi", "statistics", "gene_sets", "etp", "intervals",
    "state_pathways", "panels"
  )
  missing <- setdiff(required, names(config))
  if (length(missing)) figure7_stop("Config is missing section(s): ", paste(missing, collapse = ", "))
  if (!identical(as.character(config$module), "in_vivo_figure7")) figure7_stop("Unexpected config module")
  figure7_feature_species_contract_values(config)
  if (!identical(as.character(config$gene_sets$provider), "msigdbr") ||
      !identical(as.character(config$gene_sets$package_version), "26.1.0") ||
      !identical(
        as.character(config$gene_sets$database_release),
        "2026.1.Hs"
      ) ||
      !identical(as.character(config$gene_sets$species), "Homo sapiens")) {
    figure7_stop("Figure 7 legacy raw-rebuild gene-set contract is invalid")
  }
  if (!identical(as.character(config$raw_data$doi), "10.5281/zenodo.21463392") ||
      !identical(as.integer(config$raw_data$required_loom_files), 18L) ||
      !identical(
        as.character(config$raw_data$seurat_rds_sha256),
        "727b8a5e5868da911c3b0873838fb1b0023377ed21ea5498dc6493acbbef6d98"
      )) {
    figure7_stop("Figure 7 raw-data contract is not the reviewed Zenodo record")
  }
  source_roles <- c("endpoint_ploidy", "sample_info", "growth_curve")
  sources <- config$versioned_source_artifacts
  if (!all(source_roles %in% names(sources)) ||
      !identical(
        as.character(sources$source_revision),
        "9d7e994554dc16f48dd1ca8f6001b7011589381d"
      ) ||
      !identical(
        as.character(sources$endpoint_ploidy$default_path),
        "Data/in-vivo/all_ploidy.tsv"
      ) ||
      !identical(
        as.character(sources$sample_info$default_path),
        "Data/in-vivo/sample_info.xlsx"
      ) ||
      any(!vapply(
        source_roles,
        function(role) {
          value <- as.character(sources[[role]]$sha256)
          length(value) == 1L && grepl("^[0-9a-f]{64}$", value)
        },
        logical(1L)
      ))) {
    figure7_stop(
      "Figure 7 versioned source-artifact checksum contract is incomplete"
    )
  }
  configured_day <- suppressWarnings(as.integer(config$tgi$day))
  if (length(configured_day) != 1L || !is.finite(configured_day) || configured_day < 0L ||
      !identical(as.character(config$tgi$outcome), "day") ||
      !identical(as.character(config$tgi$matched_control_summary), "mean") ||
      !identical(as.character(config$tgi$matched_control_group), "initial_ploidy")) {
    figure7_stop("Figure 7 requires endpoint-day TGI and the mean initial-ploidy-matched control reference")
  }
  if (!is.null(tgi_day)) {
    selected_day <- suppressWarnings(as.integer(tgi_day))
    if (length(selected_day) != 1L || !is.finite(selected_day) || selected_day < 0L) {
      figure7_stop("--tgi-day must be one non-negative integer")
    }
    config$tgi$day <- selected_day
  }
  if (!isTRUE(all.equal(as.numeric(config$etp$threshold), 2.24, tolerance = 0))) {
    figure7_stop("Figure 7 requires the reference-balanced ETP threshold 2.24")
  }
  if (!identical(as.integer(unlist(config$panels$selected_direct_comparison_ids)), c(1L, 8L, 9L))) {
    figure7_stop("Figure 7 panel 7B requires comparison IDs 1, 8, and 9")
  }
  si <- config$si_figures
  si_required <- c(
    "cluster_order", "ploidy_levels", "dose_levels", "plot_shuffle_seed",
    "raw_rebuild_analysis_seed",
    "raw_rebuild_fgsea_nperm_simple",
    "umap_reduction", "pca_reduction", "cluster_id_field",
    "base_cluster_field", "clustering_resolution",
    "cluster_annotation_field", "sample_field", "dose_field",
    "ploidy_field", "context_field", "cellcycle_mapping", "inclusion",
    "source_qc_fields", "source_qc_policy", "si7_feature_species_policy",
    "si7_gene_set_database", "raw_rebuild_species_policy",
    "raw_rebuild_gene_set_database"
  )
  si_missing <- setdiff(si_required, names(si))
  if (length(si_missing)) {
    figure7_stop(
      "Figure 7 supplementary config is missing field(s): ",
      paste(si_missing, collapse = ", ")
    )
  }
  state <- config$state_pathways
  state_identity <- c(
    reference_id = as.character(state$reference_id),
    reference_kind = as.character(state$reference_kind),
    reference_canonical_publication_allowed =
      tolower(as.character(state$reference_canonical_publication_allowed)),
    generated_reference_kind =
      as.character(state$generated_reference_kind),
    generated_reference_id =
      as.character(state$generated_reference_id),
    generated_canonical_publication_allowed =
      tolower(as.character(state$generated_canonical_publication_allowed))
  )
  expected_state_identity <- c(
    reference_id = "taoli_04i_etp2_24_day17_v1",
    reference_kind = "historical_mixed_frozen",
    reference_canonical_publication_allowed = "false",
    generated_reference_kind = "generated_human_only",
    generated_reference_id =
      "runtime_state_pathway_grch_human_only_v2",
    generated_canonical_publication_allowed = "false"
  )
  if (!identical(state_identity, expected_state_identity)) {
    figure7_stop(
      "Panel-7F historical/generated publication identity is invalid"
    )
  }
  config
}

figure7_si_raw_config_contract_values <- function(config) {
  si <- config$si_figures
  mapping <- unlist(si$cellcycle_mapping, use.names = TRUE)
  inclusion <- si$inclusion
  indexed <- function(prefix, values) {
    values <- as.character(unlist(values, use.names = FALSE))
    stats::setNames(
      values,
      paste0(prefix, "[", seq_along(values), "]")
    )
  }
  named <- function(prefix, values) {
    values <- unlist(values, use.names = TRUE)
    if (is.null(names(values)) || any(!nzchar(names(values)))) {
      figure7_stop(prefix, " requires named config values")
    }
    stats::setNames(
      as.character(values),
      paste0(prefix, ".", names(values))
    )
  }
  c(
    stats::setNames(
      figure7_feature_species_contract_values(config),
      paste0(
        "feature_species.",
        names(figure7_feature_species_contract_values(config))
      )
    ),
    `versioned_source_artifacts.endpoint_ploidy.sha256` =
      as.character(config$versioned_source_artifacts$endpoint_ploidy$sha256),
    `gene_sets.provider` = as.character(config$gene_sets$provider),
    `gene_sets.package_version` =
      as.character(config$gene_sets$package_version),
    `gene_sets.database_release` =
      as.character(config$gene_sets$database_release),
    `gene_sets.species` = as.character(config$gene_sets$species),
    `si_figures.raw_rebuild_analysis_seed` =
      as.character(si$raw_rebuild_analysis_seed),
    `si_figures.raw_rebuild_fgsea_nperm_simple` =
      as.character(figure7_si_raw_fgsea_nperm_simple(config)),
    `si_figures.umap_reduction` = as.character(si$umap_reduction),
    `si_figures.cluster_id_field` = as.character(si$cluster_id_field),
    `si_figures.cluster_annotation_field` =
      as.character(si$cluster_annotation_field),
    `si_figures.ploidy_field` = as.character(si$ploidy_field),
    `si_figures.context_field` = as.character(si$context_field),
    `si_figures.dose_field` = as.character(si$dose_field),
    `si_figures.sample_field` = as.character(si$sample_field),
    indexed("si_figures.cluster_order", si$cluster_order),
    named("si_figures.cellcycle_mapping", mapping),
    `si_figures.inclusion.context` =
      as.character(inclusion$context),
    indexed(
      "si_figures.inclusion.ploidy_levels",
      inclusion$ploidy_levels
    ),
    indexed(
      "si_figures.inclusion.dose_levels",
      inclusion$dose_levels
    ),
    indexed(
      "si_figures.inclusion.required_nonmissing",
      inclusion$required_nonmissing
    )
  )
}

figure7_si_raw_config_contract_sha256 <- function(config) {
  figure7_contract_sha256(figure7_si_raw_config_contract_values(config))
}

figure7_si_raw_fgsea_nperm_simple <- function(
  config,
  supplied = NULL
) {
  reviewed <- 10000L
  configured <- suppressWarnings(as.integer(
    config$si_figures$raw_rebuild_fgsea_nperm_simple
  ))
  if (length(configured) != 1L || is.na(configured) ||
      configured != reviewed) {
    figure7_stop(
      "si_figures.raw_rebuild_fgsea_nperm_simple must equal the reviewed ",
      reviewed
    )
  }
  if (!is.null(supplied)) {
    supplied <- suppressWarnings(as.integer(supplied))
    if (length(supplied) != 1L || is.na(supplied) ||
        supplied != reviewed) {
      figure7_stop(
        "--fgsea-nperm-simple must equal the reviewed value ",
        reviewed
      )
    }
  }
  reviewed
}

figure7_si_raw_dependency_values <- function(
  si_raw_scientific_code_contract,
  si_raw_config_contract,
  si_raw_runtime_contract,
  versioned_endpoint_ploidy,
  seurat_source_dependencies
) {
  values <- c(
    si_raw_scientific_code_contract =
      si_raw_scientific_code_contract,
    si_raw_config_contract = si_raw_config_contract,
    si_raw_runtime_contract = si_raw_runtime_contract,
    versioned_endpoint_ploidy = versioned_endpoint_ploidy,
    seurat_source_dependencies
  )
  source_roles <- names(seurat_source_dependencies)
  if (is.null(source_roles) ||
      !"source_seurat_rds" %in% source_roles ||
      anyNA(source_roles) ||
      any(!nzchar(source_roles)) ||
      anyDuplicated(source_roles)) {
    figure7_stop(
      "SI raw Seurat lineage requires a unique source_seurat_rds role"
    )
  }
  if (any(!grepl("^[0-9a-f]{64}$", values))) {
    figure7_stop(
      "SI raw computational dependencies require exact SHA-256 values"
    )
  }
  values
}

figure7_si_raw_work_fingerprint <- function(dependencies) {
  figure7_contract_sha256(dependencies)
}

figure7_state_config_contract_sha256 <- function(config) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    figure7_stop("R package 'digest' is required")
  }
  state_fields <- config$state_pathways[
    setdiff(
      names(config$state_pathways),
      c(
        "reference_id", "reference_root", "expected_files",
        "comparison_tolerances", "generated_reference_id",
        "generated_reference_kind",
        "generated_canonical_publication_allowed"
      )
    )
  ]
  values <- c(
    unlist(config$feature_species, use.names = TRUE),
    config$etp$method,
    config$etp$threshold,
    unlist(config$intervals, use.names = TRUE),
    unlist(config$gene_sets, use.names = TRUE),
    unlist(state_fields, use.names = TRUE)
  )
  unname(digest::digest(
    paste(as.character(values), collapse = "\n"),
    algo = "sha256",
    serialize = FALSE
  ))
}

figure7_tgi_day <- function(config) as.integer(config$tgi$day)

figure7_tgi_measure <- function(config) {
  paste0("TGI_percent_Day_", figure7_tgi_day(config))
}

figure7_tgi_delta_measure <- function(config) {
  paste0("tumor_volume_delta_Day_", figure7_tgi_day(config))
}

figure7_tgi_label <- function(config) {
  paste("Day", figure7_tgi_day(config))
}

figure7_add_tgi_metadata <- function(data, config) {
  data$tgi_outcome <- "day"
  data$tgi_day <- figure7_tgi_day(config)
  data$tgi_measure <- figure7_tgi_measure(config)
  data$matched_control_summary <- "mean"
  data$matched_control_group <- "initial_ploidy"
  data
}

figure7_verify_checksum <- function(path, expected, label = basename(path)) {
  if (!file.exists(path)) figure7_stop("Missing ", label, ": ", path)
  expected <- as.character(expected)
  if (!grepl("^[0-9a-f]{64}$", expected)) {
    figure7_stop("Canonical SHA-256 has not been frozen for ", label, ": ", expected)
  }
  observed <- figure7_sha256(path)
  if (!identical(observed, expected)) {
    figure7_stop("SHA-256 mismatch for ", label, ": expected ", expected, "; observed ", observed)
  }
  invisible(observed)
}

figure7_verify_named_checksums <- function(paths, expected) {
  if (is.null(names(paths)) || is.null(names(expected)) ||
      !identical(sort(names(paths)), sort(names(expected)))) {
    figure7_stop("Named checksum paths and expectations must contain identical keys")
  }
  observed <- vapply(names(expected), function(key) {
    figure7_verify_checksum(paths[[key]], expected[[key]], key)
  }, character(1L))
  invisible(observed)
}

figure7_assert_empty_output <- function(path) {
  if (dir.exists(path)) {
    existing <- list.files(path, recursive = TRUE, all.files = TRUE, no.. = TRUE,
                           include.dirs = FALSE)
    allowed_manager_skeleton <- c("logs/stdout.log", "logs/stderr.log")
    unexpected <- setdiff(existing, allowed_manager_skeleton)
    if (length(unexpected)) {
      figure7_stop("Output directory contains pre-existing analysis files: ", paste(unexpected, collapse = ", "))
    }
  }
  invisible(TRUE)
}

figure7_prepare_output <- function(path) {
  figure7_assert_empty_output(path)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  for (subdir in c("figures", "tables", "metadata", "logs")) {
    dir.create(file.path(path, subdir), recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

figure7_read_tsv <- function(path, required = character(), label = basename(path)) {
  if (!file.exists(path)) figure7_stop("Missing ", label, ": ", path)
  data <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  character_columns <- vapply(data, is.character, logical(1L))
  data[character_columns] <- lapply(data[character_columns], function(x) {
    x <- gsub("\\n", "\n", x, fixed = TRUE)
    x <- gsub("\\r", "\r", x, fixed = TRUE)
    gsub("\\t", "\t", x, fixed = TRUE)
  })
  missing <- setdiff(required, names(data))
  if (length(missing)) figure7_stop(label, " is missing column(s): ", paste(missing, collapse = ", "))
  data
}

figure7_write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  encoded <- as.data.frame(data, stringsAsFactors = FALSE)
  character_columns <- vapply(encoded, is.character, logical(1L))
  encoded[character_columns] <- lapply(encoded[character_columns], function(x) {
    x <- gsub("\t", "\\t", x, fixed = TRUE)
    x <- gsub("\r", "\\r", x, fixed = TRUE)
    gsub("\n", "\\n", x, fixed = TRUE)
  })
  utils::write.table(encoded, path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "NA")
  invisible(path)
}

figure7_save_panel <- function(plot, pdf_path, width, height, png_dpi = 300) {
  png_path <- sub("[.]pdf$", ".png", pdf_path, ignore.case = TRUE)
  if (identical(png_path, pdf_path)) figure7_stop("Figure 7 panel path must end in .pdf: ", pdf_path)
  ggplot2::ggsave(pdf_path, plot, device = grDevices::cairo_pdf,
                  width = width, height = height, units = "in")
  ggplot2::ggsave(png_path, plot, device = "png", dpi = png_dpi,
                  width = width, height = height, units = "in", bg = "white")
  paths <- c(pdf_path, png_path)
  if (any(!file.exists(paths)) || any(file.info(paths)$size <= 0)) {
    figure7_stop("Failed to write PDF/PNG panel pair: ", pdf_path)
  }
  invisible(paths)
}

figure7_panel_ids <- function(include_panel_f = TRUE) {
  if (isTRUE(include_panel_f)) c("7A", "7B", "7C", "7D", "7E", "7F") else c("7A", "7B", "7C", "7D", "7E")
}

figure7_panel_filenames <- function(config, panel_ids = figure7_panel_ids(TRUE)) {
  filenames <- unname(unlist(config$panels$filenames[panel_ids], use.names = FALSE))
  sub("day17", paste0("day", figure7_tgi_day(config)), filenames, ignore.case = TRUE)
}

figure7_panel_asset_filenames <- function(config, panel_ids = figure7_panel_ids(TRUE)) {
  pdfs <- figure7_panel_filenames(config, panel_ids)
  c(pdfs, sub("[.]pdf$", ".png", pdfs, ignore.case = TRUE))
}

figure7_validate_figure_inventory <- function(output_dir, config, panel_ids = figure7_panel_ids(TRUE)) {
  figures_dir <- file.path(output_dir, "figures")
  all_figures <- list.files(
    output_dir, pattern = "[.](pdf|png|svg|tiff?|jpg|jpeg)$",
    recursive = TRUE, full.names = TRUE, ignore.case = TRUE
  )
  observed <- sort(basename(all_figures))
  expected <- sort(figure7_panel_asset_filenames(config, panel_ids))
  if (!identical(observed, expected) ||
      any(dirname(normalizePath(all_figures)) != normalizePath(figures_dir))) {
    figure7_stop("Figure inventory mismatch. Expected: ", paste(expected, collapse = ", "),
                 "; observed: ", paste(observed, collapse = ", "))
  }
  if (any(file.info(all_figures)$size <= 0)) figure7_stop("Every Figure 7 asset must be nonempty")
  invisible(TRUE)
}

figure7_copy_file <- function(from, to) {
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(from, to, overwrite = FALSE, copy.mode = TRUE)) figure7_stop("Could not copy ", from, " to ", to)
  invisible(to)
}

figure7_workflow_scalar <- function(workflow, name, default = "") {
  if (!is.list(workflow) || length(name) != 1L || is.na(name) ||
      !nzchar(name) || length(default) != 1L || is.na(default)) {
    figure7_stop("Invalid workflow metadata scalar request")
  }
  value <- workflow[[name]]
  if (is.null(value) || length(value) != 1L) {
    return(as.character(default))
  }
  value <- as.character(value)
  if (is.na(value) || !nzchar(value)) as.character(default) else value
}
