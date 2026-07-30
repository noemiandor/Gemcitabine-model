#!/usr/bin/env Rscript

# Cache-first orchestrator for Supplementary Figures 4-7.
# A valid reviewed 11-table cache renders without raw-analysis dependencies.
# The fallback prefers a validated reconstructed Seurat object, resumes its
# upstream stages when possible, and uses the deposited Seurat RDS only when no
# reconstruction source is available.

parse_args <- function(tokens) {
  out <- list()
  i <- 1L
  while (i <= length(tokens)) {
    token <- tokens[[i]]
    if (grepl("^--[^=]+=", token)) {
      key <- sub("^--([^=]+)=.*$", "\\1", token)
      out[[key]] <- sub("^--[^=]+=", "", token)
      i <- i + 1L
    } else if (startsWith(token, "--")) {
      key <- sub("^--", "", token)
      if (i < length(tokens) && !startsWith(tokens[[i + 1L]], "--")) {
        out[[key]] <- tokens[[i + 1L]]
        i <- i + 2L
      } else {
        out[[key]] <- "true"
        i <- i + 1L
      }
    } else {
      stop("Unexpected positional argument: ", token, call. = FALSE)
    }
  }
  out
}

arg <- function(args, name, default = NULL) {
  value <- args[[name]]
  if (!is.null(value) && length(value) == 1L && nzchar(value)) value else default
}

flag <- function(args, name, default = FALSE) {
  value <- arg(args, name, if (default) "true" else "false")
  normalized <- tolower(trimws(value))
  if (!normalized %in% c("true", "false", "1", "0", "yes", "no")) {
    stop("--", name, " must be true or false", call. = FALSE)
  }
  normalized %in% c("true", "1", "yes")
}

script_path <- function() {
  hit <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(hit)) return(NA_character_)
  normalizePath(sub("^--file=", "", hit[[1L]]), mustWork = TRUE)
}

resolve_path <- function(path, root, must_work = FALSE) {
  candidate <- if (grepl("^/", path)) path else file.path(root, path)
  normalizePath(candidate, mustWork = must_work)
}

sha256 <- function(path) {
  if (!file.exists(path)) stop("Cannot checksum missing file: ", path, call. = FALSE)
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("R package 'digest' is required", call. = FALSE)
  }
  unname(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

run_stage <- function(label, script, values, log_path) {
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  values <- values[!vapply(values, is.null, logical(1L))]
  cli <- paste0(
    "--", names(values), "=",
    vapply(values, as.character, character(1L))
  )
  message("[si-figures] ", label)
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(script), vapply(cli, shQuote, character(1L))),
    stdout = log_path,
    stderr = log_path
  )
  if (!identical(status, 0L)) {
    detail <- if (file.exists(log_path)) {
      paste(tail(readLines(log_path, warn = FALSE), 80L), collapse = "\n")
    } else {
      "No stage log was written."
    }
    stop(label, " failed; log: ", log_path, "\n", detail, call. = FALSE)
  }
  invisible(log_path)
}

validate_cache <- function(
  validator,
  cache_dir,
  policy = "corrected-human-only"
) {
  cli <- c(shQuote(validator), "--cache-dir", shQuote(cache_dir))
  if (!identical(policy, "corrected-human-only")) {
    cli <- c(cli, "--si7-policy", shQuote(policy))
  }
  identical(suppressWarnings(system2("python3", cli)), 0L)
}

read_manifest <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    utils::read.delim(
      path, check.names = FALSE, stringsAsFactors = FALSE,
      quote = "", comment.char = ""
    ),
    error = function(error) NULL
  )
}

portable_locator <- function(path, repo_root) {
  normalized <- normalizePath(path, mustWork = TRUE)
  normalized_root <- normalizePath(repo_root, mustWork = TRUE)
  prefix <- paste0(normalized_root, .Platform$file.sep)
  if (startsWith(normalized, prefix)) {
    substring(normalized, nchar(prefix) + 1L)
  } else {
    paste0("external:", basename(normalized))
  }
}

augment_manifest_with_scvelo_audit <- function(
  source_manifest,
  scvelo_metrics,
  destination,
  repo_root
) {
  manifest <- read_manifest(source_manifest)
  if (is.null(manifest) ||
      !identical(
        names(manifest),
        c("role", "locator", "sha256", "bytes")
      )) {
    stop(
      "Cannot augment malformed generated SI input manifest",
      call. = FALSE
    )
  }
  manifest <- manifest[
    manifest$role != "audit_optional_scvelo",
    ,
    drop = FALSE
  ]
  manifest <- rbind(
    manifest,
    data.frame(
      role = "audit_optional_scvelo",
      locator = portable_locator(scvelo_metrics, repo_root),
      sha256 = sha256(scvelo_metrics),
      bytes = as.numeric(file.info(scvelo_metrics)$size),
      stringsAsFactors = FALSE
    )
  )
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    manifest,
    destination,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = ""
  )
  invisible(destination)
}

manifest_matches <- function(manifest, expected) {
  if (is.null(manifest) ||
      !all(c("role", "sha256") %in% names(manifest)) ||
      anyNA(manifest$role) ||
      any(!nzchar(manifest$role)) ||
      anyDuplicated(manifest$role)) {
    return(FALSE)
  }
  observed <- stats::setNames(
    tolower(as.character(manifest$sha256)),
    as.character(manifest$role)
  )
  extras <- setdiff(names(observed), names(expected))
  all(startsWith(extras, "audit_")) &&
    all(names(expected) %in% names(observed)) &&
    identical(
      unname(observed[names(expected)]),
      unname(tolower(as.character(expected)))
    )
}

is_within <- function(path, root) {
  path <- normalizePath(path, mustWork = FALSE)
  root <- normalizePath(root, mustWork = FALSE)
  identical(path, root) ||
    startsWith(path, paste0(root, .Platform$file.sep))
}

quarantine_generated <- function(path, root) {
  if (!is_within(path, root)) {
    stop("Refusing to modify SI cache outside --intermediate-dir: ", path,
         call. = FALSE)
  }
  if (!dir.exists(path)) return(invisible(""))
  target <- paste0(
    path, ".stale.",
    format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"),
    ".", Sys.getpid()
  )
  if (!file.rename(path, target)) {
    stop("Cannot preserve stale generated SI cache: ", path, call. = FALSE)
  }
  invisible(target)
}

render_cache <- function(
  renderer,
  cache,
  config,
  output,
  generated,
  log_name = "00_render.log",
  upstream_input_manifest = NULL
) {
  values <- list(
    "table-cache-dir" = cache,
    config = config,
    "output-dir" = output
  )
  if (generated) {
    values[["allow-generated-human-only-si7"]] <- "true"
  }
  if (!is.null(upstream_input_manifest) &&
      nzchar(upstream_input_manifest)) {
    values[["upstream-input-manifest"]] <- upstream_input_manifest
  }
  run_stage(
    "Render Supplementary Figure 4-7 composites",
    renderer,
    values,
    file.path(output, "logs", log_name)
  )
}

computational_dependencies <- function(
  si_raw_scientific_code_contract,
  si_raw_config_contract,
  si_raw_runtime_contract,
  versioned_endpoint_ploidy,
  seurat_source_dependencies,
  audit_optional_scvelo = NULL
) {
  # audit_optional_scvelo is accepted so callers/tests can state the complete
  # input set, but it is intentionally absent from the computational lineage.
  c(
    si_raw_scientific_code_contract =
      si_raw_scientific_code_contract,
    si_raw_config_contract = si_raw_config_contract,
    si_raw_runtime_contract = si_raw_runtime_contract,
    versioned_endpoint_ploidy = versioned_endpoint_ploidy,
    seurat_source_dependencies
  )
}

computational_fingerprint <- function(dependency_hashes) {
  if (!length(dependency_hashes) || is.null(names(dependency_hashes)) ||
      anyNA(names(dependency_hashes)) || any(!nzchar(names(dependency_hashes)))) {
    stop("Computational dependency hashes must be a named vector", call. = FALSE)
  }
  if (anyDuplicated(names(dependency_hashes))) {
    stop("Computational fingerprint roles must be unique", call. = FALSE)
  }
  dependency_names <- names(dependency_hashes)
  dependency_hashes <- tolower(as.character(dependency_hashes))
  names(dependency_hashes) <- dependency_names
  if (any(!grepl("^[0-9a-f]{64}$", dependency_hashes))) {
    stop("Computational dependencies require SHA-256 values", call. = FALSE)
  }
  dependency_hashes <- dependency_hashes[order(names(dependency_hashes))]
  digest::digest(
    paste(
      paste(names(dependency_hashes), dependency_hashes, sep = "="),
      collapse = "\n"
    ),
    algo = "sha256",
    serialize = FALSE
  )
}

manifest_computational_dependencies <- function(manifest) {
  required <- c("role", "locator", "sha256", "bytes")
  if (is.null(manifest) ||
      !identical(names(manifest), required) ||
      anyNA(manifest$role) ||
      any(!nzchar(as.character(manifest$role))) ||
      anyDuplicated(manifest$role) ||
      anyNA(manifest$locator) ||
      any(!nzchar(as.character(manifest$locator)))) {
    return(NULL)
  }
  bytes <- suppressWarnings(as.numeric(manifest$bytes))
  if (anyNA(bytes) || any(!is.finite(bytes)) || any(bytes < 0)) {
    return(NULL)
  }
  roles <- as.character(manifest$role)
  hashes <- tolower(as.character(manifest$sha256))
  names(hashes) <- roles
  if (any(!grepl("^[0-9a-f]{64}$", hashes))) return(NULL)
  hashes[!startsWith(roles, "audit_")]
}

read_key_values <- function(path) {
  table <- read_manifest(path)
  if (is.null(table) ||
      !all(c("key", "value") %in% names(table)) ||
      anyNA(table$key) ||
      any(!nzchar(as.character(table$key))) ||
      anyDuplicated(table$key)) {
    return(NULL)
  }
  stats::setNames(as.character(table$value), as.character(table$key))
}

si_reconstruction_contract <- function(
  generator,
  environment_lock,
  config,
  all_ploidy,
  sample_info
) {
  transitive_roles <- generator$figure7_upstream_stage_dependency_roles(
    "final"
  )
  scientific <- generator$figure7_upstream_cumulative_code_contracts(
    "final"
  )
  common <- generator$figure7_upstream_common_contract(
    environment_lock,
    config,
    jobs = 1L
  )
  fixed <- c(
    stats::setNames(
      c(sha256(all_ploidy), sha256(sample_info)),
      paste0(
        "seurat_transitive_dependency:",
        c("all_ploidy", "sample_info")
      )
    ),
    stats::setNames(
      unname(scientific),
      paste0("seurat_", names(scientific))
    ),
    seurat_upstream_scientific_common_contract =
      figure7_contract_sha256(common)
  )
  list(
    transitive_roles = paste0(
      "seurat_transitive_dependency:",
      transitive_roles
    ),
    scientific_roles = paste0("seurat_", names(scientific)),
    fixed = fixed
  )
}

si_current_h5_inventory_hash <- function(
  generator,
  cellranger_root,
  sample_info
) {
  if (!nzchar(cellranger_root) || !dir.exists(cellranger_root)) return("")
  inventory <- generator$figure7_upstream_h5_inventory(
    normalizePath(cellranger_root, mustWork = TRUE),
    generator$figure7_upstream_read_sample_info(sample_info)
  )
  generator$figure7_upstream_h5_inventory_digest(
    inventory$sample,
    inventory$sha256
  )
}

si_present_reconstruction_hashes <- function(generator, upstream_dir) {
  if (!nzchar(upstream_dir) || !dir.exists(upstream_dir)) {
    return(character())
  }
  definitions <- generator$figure7_upstream_stage_definitions(upstream_dir)
  stage_roles <- c(
    integrated = "integrated_rds",
    cell_cycle = "cell_cycle_rds",
    refined = "refined_rds",
    merged = "merged_rds"
  )
  paths <- vapply(
    names(stage_roles),
    function(stage) {
      generator$figure7_upstream_stage_output(definitions[[stage]])
    },
    character(1L)
  )
  present <- file.exists(paths)
  if (!any(present)) return(character())
  stats::setNames(
    vapply(paths[present], sha256, character(1L)),
    paste0(
      "seurat_transitive_dependency:",
      unname(stage_roles[names(paths)[present]])
    )
  )
}

si_cache_candidate_dependencies <- function(
  build_dir,
  base_dependencies,
  deposited_sha256,
  reconstruction_contract,
  generator,
  explicit_rds = "",
  deposited_rds = "",
  upstream_dir = "",
  cellranger_root = "",
  sample_info = "",
  require_fingerprint_name = TRUE
) {
  manifest <- read_manifest(
    file.path(build_dir, "metadata", "input_manifest.tsv")
  )
  observed <- manifest_computational_dependencies(manifest)
  if (is.null(observed)) return(NULL)
  locators <- as.character(manifest$locator)
  if (any(grepl("^/|^[A-Za-z]:[/\\\\]", locators)) ||
      any(grepl("(^|/)[.][.](/|$)", locators))) {
    return(NULL)
  }
  source_role <- "source_seurat_rds"
  if (!source_role %in% names(observed)) return(NULL)
  source_hash <- unname(observed[[source_role]])
  reconstructed_roles <- c(
    reconstruction_contract$transitive_roles,
    reconstruction_contract$scientific_roles,
    "seurat_upstream_scientific_common_contract"
  )
  present_reconstructed_roles <- intersect(
    names(observed),
    reconstructed_roles
  )
  reconstructed <- if (!length(present_reconstructed_roles)) {
    FALSE
  } else if (setequal(
    present_reconstructed_roles,
    reconstructed_roles
  )) {
    TRUE
  } else {
    return(NULL)
  }
  if (!reconstructed &&
      !identical(source_hash, tolower(deposited_sha256))) {
    return(NULL)
  }
  source_roles <- c(
    source_role,
    if (reconstructed) reconstructed_roles
  )
  required_audit_roles <- c(
    "audit_si_table_builder",
    "audit_environment_validator",
    "audit_environment_lock",
    "audit_figure7_config",
    "audit_seurat_selection",
    if (reconstructed) {
      c(
        "audit_seurat_upstream_generator",
        "audit_seurat_upstream_helper",
        "audit_seurat_final_stage_manifest",
        "audit_seurat_reconstruction_manifest"
      )
    }
  )
  if (!all(required_audit_roles %in% as.character(manifest$role))) {
    return(NULL)
  }
  expected_roles <- c(names(base_dependencies), source_roles)
  if (!setequal(names(observed), expected_roles)) return(NULL)
  if (!identical(
    unname(observed[names(base_dependencies)]),
    unname(tolower(as.character(base_dependencies)))
  )) {
    return(NULL)
  }
  if (reconstructed && !identical(
    unname(observed[names(reconstruction_contract$fixed)]),
    unname(tolower(as.character(reconstruction_contract$fixed)))
  )) {
    return(NULL)
  }

  present_sources <- character()
  if (nzchar(explicit_rds) && file.exists(explicit_rds)) {
    present_sources <- c(present_sources, sha256(explicit_rds))
  } else {
    upstream_final <- if (nzchar(upstream_dir)) {
      figure7_seurat_upstream_final_path(upstream_dir)
    } else {
      ""
    }
    if (reconstructed && nzchar(upstream_final) &&
        file.exists(upstream_final)) {
      present_sources <- c(present_sources, sha256(upstream_final))
    }
    if (!reconstructed && nzchar(deposited_rds) &&
        file.exists(deposited_rds)) {
      present_sources <- c(present_sources, sha256(deposited_rds))
    }
  }
  if (length(present_sources) &&
      any(tolower(present_sources) != source_hash)) {
    return(NULL)
  }

  if (reconstructed) {
    relevant_upstream_dir <- upstream_dir
    if (nzchar(explicit_rds) && file.exists(explicit_rds)) {
      inferred_upstream_dir <- figure7_infer_seurat_upstream_root(
        explicit_rds
      )
      if (nzchar(inferred_upstream_dir) &&
          dir.exists(inferred_upstream_dir)) {
        relevant_upstream_dir <- inferred_upstream_dir
      }
    }
    present_stage_hashes <- si_present_reconstruction_hashes(
      generator,
      relevant_upstream_dir
    )
    if (length(present_stage_hashes) && !identical(
      unname(observed[names(present_stage_hashes)]),
      unname(tolower(present_stage_hashes))
    )) {
      return(NULL)
    }
    current_h5 <- si_current_h5_inventory_hash(
      generator,
      cellranger_root,
      sample_info
    )
    if (nzchar(current_h5) && !identical(
      unname(observed[[
        "seurat_transitive_dependency:h5_inventory"
      ]]),
      tolower(current_h5)
    )) {
      return(NULL)
    }
  }

  fingerprint <- computational_fingerprint(observed)
  expected_name <- paste0(
    "generated_human_only_",
    substr(fingerprint, 1L, 20L)
  )
  if (isTRUE(require_fingerprint_name) &&
      !identical(basename(build_dir), expected_name)) {
    return(NULL)
  }
  run_config <- read_key_values(
    file.path(build_dir, "metadata", "run_config.tsv")
  )
  expected_run_config <- c(
    schema_version = "1",
    module = "si_figures_raw_table_build",
    output_table_count = "11",
    figures_supported = "4,5,6,7",
    si7_canonical_publication_allowed = "false",
    si7_species_policy_id = "grch_human_tumor_only_v2",
    si7_source_data_layer_reused = "false",
    seurat_source_kind = if (reconstructed) {
      "manifested_reconstruction"
    } else {
      "deposited"
    },
    si_raw_scientific_code_contract_sha256 =
      unname(base_dependencies[["si_raw_scientific_code_contract"]]),
    si_raw_config_contract_sha256 =
      unname(base_dependencies[["si_raw_config_contract"]]),
    si_raw_runtime_contract_sha256 =
      unname(base_dependencies[["si_raw_runtime_contract"]]),
    work_dependency_fingerprint = fingerprint
  )
  if (is.null(run_config) ||
      !all(names(expected_run_config) %in% names(run_config)) ||
      !identical(
        unname(tolower(run_config[names(expected_run_config)])),
        unname(tolower(expected_run_config))
      )) {
    return(NULL)
  }
  list(
    build_dir = build_dir,
    dependencies = observed,
    fingerprint = fingerprint,
    reconstructed = reconstructed,
    source_hash = source_hash
  )
}

find_reusable_generated_cache <- function(
  intermediate_root,
  validator,
  base_dependencies,
  deposited_sha256,
  reconstruction_contract,
  generator,
  explicit_rds = "",
  deposited_rds = "",
  upstream_dir = "",
  cellranger_root = "",
  sample_info = "",
  generated_cache_dir = ""
) {
  candidates <- if (nzchar(generated_cache_dir)) {
    generated_cache_dir
  } else if (dir.exists(intermediate_root)) {
    list.dirs(intermediate_root, recursive = FALSE, full.names = TRUE)
  } else {
    character()
  }
  candidates <- candidates[
    grepl("^generated_human_only_[0-9a-f]{20}$", basename(candidates))
  ]
  reusable <- lapply(candidates, function(build_dir) {
    if (!validate_cache(
      validator,
      file.path(build_dir, "tables"),
      policy = "generated-human-only"
    )) {
      return(NULL)
    }
    si_cache_candidate_dependencies(
      build_dir,
      base_dependencies,
      deposited_sha256,
      reconstruction_contract,
      generator,
      explicit_rds,
      deposited_rds,
      upstream_dir,
      cellranger_root,
      sample_info
    )
  })
  reusable <- Filter(Negate(is.null), reusable)
  if (length(reusable) > 1L) {
    stop(
      "Multiple generated SI caches match the current scientific lineage; ",
      "select one with --generated-cache-dir: ",
      paste(vapply(reusable, function(x) x$build_dir, character(1L)),
            collapse = ", "),
      call. = FALSE
    )
  }
  if (nzchar(generated_cache_dir) && !length(reusable)) {
    stop(
      "Explicit --generated-cache-dir is not a reusable generated SI cache: ",
      generated_cache_dir,
      call. = FALSE
    )
  }
  if (length(reusable)) reusable[[1L]] else NULL
}

raw_builder_args <- function(
  all_ploidy,
  seurat_rds,
  config,
  output,
  seurat_upstream_dir = NULL,
  sample_info = NULL,
  cellranger_root = NULL,
  scvelo_metrics = NULL,
  work_cache_dir = NULL,
  work_dependency_fingerprint = NULL,
  workers = NULL
) {
  values <- list(
    "all-ploidy" = all_ploidy,
    "seurat-rds" = seurat_rds,
    config = config,
    "output-dir" = output
  )
  if (!is.null(scvelo_metrics) && nzchar(scvelo_metrics)) {
    values[["scvelo-metrics"]] <- scvelo_metrics
  }
  if (!is.null(seurat_upstream_dir) && nzchar(seurat_upstream_dir)) {
    values[["seurat-upstream-dir"]] <- seurat_upstream_dir
  }
  if (!is.null(sample_info) && nzchar(sample_info)) {
    values[["sample-info"]] <- sample_info
  }
  if (!is.null(cellranger_root) && nzchar(cellranger_root)) {
    values[["cellranger-root"]] <- cellranger_root
  }
  if (!is.null(work_cache_dir) && nzchar(work_cache_dir)) {
    values[["work-cache-dir"]] <- work_cache_dir
  }
  if (!is.null(work_dependency_fingerprint) &&
      nzchar(work_dependency_fingerprint)) {
    values[["work-dependency-fingerprint"]] <-
      work_dependency_fingerprint
  }
  if (!is.null(workers)) values[["workers"]] <- as.integer(workers)
  values
}

raw_audit_args <- function(
  all_ploidy,
  seurat_rds,
  config,
  scvelo_metrics,
  seurat_upstream_dir = NULL,
  sample_info = NULL,
  cellranger_root = NULL
) {
  values <- list(
    "all-ploidy" = all_ploidy,
    "seurat-rds" = seurat_rds,
    config = config,
    "scvelo-metrics" = scvelo_metrics,
    "audit-only" = "true"
  )
  if (!is.null(seurat_upstream_dir) && nzchar(seurat_upstream_dir)) {
    values[["seurat-upstream-dir"]] <- seurat_upstream_dir
  }
  if (!is.null(sample_info) && nzchar(sample_info)) {
    values[["sample-info"]] <- sample_info
  }
  if (!is.null(cellranger_root) && nzchar(cellranger_root)) {
    values[["cellranger-root"]] <- cellranger_root
  }
  values
}

raw_rds_acquisition <- function(
  rds_arg,
  raw_data_dir,
  seurat_filename,
  allow_download
) {
  explicit <- nzchar(rds_arg)
  list(
    path = if (explicit) {
      rds_arg
    } else {
      file.path(raw_data_dir, seurat_filename)
    },
    explicit = explicit,
    verify_with_downloader = !explicit,
    allow_download = isTRUE(allow_download)
  )
}

raw_rds_needed <- function(build_reusable, scvelo_metrics = NULL) {
  !isTRUE(build_reusable) ||
    (!is.null(scvelo_metrics) && nzchar(scvelo_metrics))
}

materialize_seurat_selection <- function(
  selection,
  script_dir,
  environment_lock,
  config,
  config_path,
  all_ploidy,
  sample_info,
  jobs,
  raw_data_dir,
  raw_manifest,
  allow_download,
  log_dir
) {
  if (identical(selection$action, "build_upstream")) {
    run_stage(
      if (nzchar(selection$resumable_stage)) {
        paste0(
          "Resume Seurat reconstruction from validated ",
          selection$resumable_stage,
          " cache"
        )
      } else {
        "Reconstruct final Seurat object from Cell Ranger matrices"
      },
      file.path(
        script_dir,
        "..",
        "figure7",
        "generate_final_seurat_from_cellranger.R"
      ),
      figure7_seurat_upstream_generator_args(
        selection,
        environment_lock,
        config_path,
        all_ploidy,
        sample_info,
        jobs
      ),
      file.path(log_dir, "00_seurat_upstream.log")
    )
    validation <- figure7_validate_seurat_upstream_artifact(
      output_root = selection$upstream_dir,
      module_dir = file.path(script_dir, "..", "figure7"),
      environment_lock = environment_lock,
      config = config,
      all_ploidy = all_ploidy,
      sample_info = sample_info,
      expected_rds = selection$rds,
      cellranger_root = if (dir.exists(selection$cellranger_root)) {
        selection$cellranger_root
      } else {
        ""
      }
    )
    selection$source <- if (nzchar(selection$resumable_stage)) {
      paste0(
        "resumed_reconstructed_rds_from_",
        selection$resumable_stage
      )
    } else {
      "reconstructed_rds_from_cellranger"
    }
    selection$action <- "reuse"
    selection$rds <- validation$rds
    selection$rds_sha256 <- validation$rds_sha256
    selection$final_stage_manifest <- validation$final_stage_manifest
    selection$final_stage_manifest_sha256 <-
      validation$final_stage_manifest_sha256
    selection$reconstruction_manifest <-
      validation$reconstruction_manifest
    selection$reconstruction_manifest_sha256 <-
      validation$reconstruction_manifest_sha256
    selection$transitive_dependencies <- validation$dependencies
    selection$scientific_code_contracts <-
      validation$scientific_code_contracts
    selection$upstream_common_contract <-
      validation$upstream_common_contract
  }

  if (selection$action %in% c(
    "validate_deposited",
    "download_deposited"
  )) {
    if (identical(selection$action, "download_deposited") &&
        !isTRUE(allow_download)) {
      stop(
        "Deposited Seurat RDS is missing and automatic download is disabled",
        call. = FALSE
      )
    }
    run_stage(
      "Download/validate deposited Seurat RDS only",
      file.path(
        script_dir,
        "..",
        "figure7",
        "download_figure7_raw_data.R"
      ),
      list(
        "raw-data-dir" = raw_data_dir,
        manifest = raw_manifest,
        roles = "seurat_rds",
        "allow-download" = if (isTRUE(allow_download)) "true" else "false"
      ),
      file.path(log_dir, "00_download_seurat_rds.log")
    )
    selection$action <- "reuse"
    selection$source <- if (file.exists(selection$rds)) {
      "validated_deposited_rds"
    } else {
      selection$source
    }
  }

  if (!file.exists(selection$rds)) {
    stop(
      "Selected Seurat RDS is unavailable after input preparation: ",
      selection$rds,
      call. = FALSE
    )
  }
  observed <- sha256(selection$rds)
  if (nzchar(selection$rds_sha256) &&
      !identical(observed, tolower(selection$rds_sha256))) {
    stop("Selected Seurat RDS changed before SI computation", call. = FALSE)
  }
  selection$rds_sha256 <- observed
  selection
}

main <- function(args = parse_args(commandArgs(trailingOnly = TRUE))) {
  script <- script_path()
  if (is.na(script)) stop("Cannot resolve SI orchestrator path", call. = FALSE)
  script_dir <- dirname(script)
  repo_root <- normalizePath(
    dirname(dirname(dirname(script_dir))),
    mustWork = TRUE
  )
  if (!requireNamespace("yaml", quietly = TRUE) ||
      !requireNamespace("digest", quietly = TRUE)) {
    stop("R packages yaml and digest are required", call. = FALSE)
  }

  config_path <- resolve_path(
    arg(args, "config", "Code/in-vivo/figure7/figure7_config.yaml"),
    repo_root,
    must_work = TRUE
  )
  config <- yaml::read_yaml(config_path)
  mode <- tolower(arg(args, "mode", "auto"))
  if (!mode %in% c("auto", "plot-only", "full-workflow")) {
    stop(
      "--mode must be auto, plot-only, or full-workflow",
      call. = FALSE
    )
  }
  output_arg <- arg(args, "output-dir", NULL)
  if (is.null(output_arg)) {
    stop("Missing required --output-dir", call. = FALSE)
  }
  output_dir <- resolve_path(output_arg, repo_root, must_work = FALSE)
  canonical_cache <- resolve_path(
    arg(
      args,
      "table-cache-dir",
      as.character(config$si_figures$cache_root)
    ),
    repo_root,
    must_work = FALSE
  )
  renderer <- file.path(script_dir, "generate_supplementary_figures.R")
  builder <- file.path(script_dir, "build_raw_supplementary_tables.R")
  validator <- file.path(
    repo_root,
    "Code",
    "tools",
    "validate_si_figures_table_cache.py"
  )

  # This branch deliberately precedes every raw-source, environment, and
  # reconstruction check.
  if (!identical(mode, "full-workflow") &&
      validate_cache(validator, canonical_cache)) {
    render_cache(renderer, canonical_cache, config_path, output_dir, FALSE)
    return(invisible(output_dir))
  }
  if (identical(mode, "plot-only")) {
    stop(
      "The reviewed SI table cache is invalid or missing; ",
      "plot-only mode will not acquire raw data. Use --mode=full-workflow ",
      "or --mode=auto to enable raw reconstruction.",
      call. = FALSE
    )
  }

  intermediate_root <- resolve_path(
    arg(
      args,
      "intermediate-dir",
      "Results/in-vivo/SI_figures/intermediates"
    ),
    repo_root,
    must_work = FALSE
  )
  dir.create(intermediate_root, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(intermediate_root)) {
    stop(
      "Cannot create SI intermediate root: ",
      intermediate_root,
      call. = FALSE
    )
  }

  figure7_dir <- file.path(repo_root, "Code", "in-vivo", "figure7")
  environment_validator <- file.path(
    figure7_dir,
    "src",
    "common_io.R"
  )
  seurat_selection_script <- file.path(
    figure7_dir,
    "src",
    "seurat_upstream_selection.R"
  )
  runtime_environment <- environment(si_reconstruction_contract)
  sys.source(environment_validator, envir = runtime_environment)
  sys.source(seurat_selection_script, envir = runtime_environment)
  config <- figure7_read_config(config_path)

  all_ploidy <- resolve_path(
    arg(
      args,
      "all-ploidy",
      as.character(
        config$versioned_source_artifacts$endpoint_ploidy$default_path
      )
    ),
    repo_root,
    must_work = TRUE
  )
  all_ploidy_sha256 <- sha256(all_ploidy)
  if (!identical(
    all_ploidy_sha256,
    as.character(
      config$versioned_source_artifacts$endpoint_ploidy$sha256
    )
  )) {
    stop(
      "Endpoint-ploidy input differs from its frozen checksum",
      call. = FALSE
    )
  }
  sample_info <- resolve_path(
    arg(
      args,
      "sample-info",
      as.character(
        config$versioned_source_artifacts$sample_info$default_path
      )
    ),
    repo_root,
    must_work = TRUE
  )
  if (!identical(
    sha256(sample_info),
    as.character(config$versioned_source_artifacts$sample_info$sha256)
  )) {
    stop("Sample-info input differs from its frozen checksum", call. = FALSE)
  }
  environment_lock <- resolve_path(
    as.character(config$raw_data$environment_lock),
    repo_root,
    must_work = TRUE
  )

  builder_contract_environment <- new.env(parent = runtime_environment)
  sys.source(builder, envir = builder_contract_environment)
  si_raw_scientific_code_contract_sha256 <-
    builder_contract_environment$si_raw_scientific_code_contract(builder)
  si_raw_config_contract_sha256 <-
    figure7_si_raw_config_contract_sha256(config)
  si_raw_runtime_contract_sha256 <-
    figure7_environment_stage_contract_sha256(
      environment_lock,
      "r",
      "si_raw"
    )
  base_dependencies <- computational_dependencies(
    si_raw_scientific_code_contract_sha256,
    si_raw_config_contract_sha256,
    si_raw_runtime_contract_sha256,
    all_ploidy_sha256,
    character()
  )

  generator <- figure7_seurat_upstream_generator_environment(figure7_dir)
  reconstruction_contract <- si_reconstruction_contract(
    generator,
    environment_lock,
    config,
    all_ploidy,
    sample_info
  )

  raw_data_dir <- resolve_path(
    arg(
      args,
      "raw-data-dir",
      as.character(config$raw_data$default_root)
    ),
    repo_root,
    must_work = FALSE
  )
  deposited_rds <- file.path(
    raw_data_dir,
    as.character(config$raw_data$seurat_filename)
  )
  explicit_rds_arg <- arg(args, "seurat-rds", "")
  explicit_rds <- if (nzchar(explicit_rds_arg)) {
    resolve_path(explicit_rds_arg, repo_root, must_work = FALSE)
  } else {
    ""
  }
  upstream_dir <- resolve_path(
    arg(
      args,
      "seurat-upstream-dir",
      file.path(intermediate_root, "seurat_upstream")
    ),
    repo_root,
    must_work = FALSE
  )
  cellranger_arg <- arg(args, "cellranger-root", "")
  cellranger_root <- if (nzchar(cellranger_arg)) {
    resolve_path(cellranger_arg, repo_root, must_work = FALSE)
  } else {
    ""
  }
  metrics_arg <- arg(args, "scvelo-metrics", "")
  metrics_path <- if (nzchar(metrics_arg)) {
    resolve_path(metrics_arg, repo_root, must_work = TRUE)
  } else {
    NULL
  }
  generated_cache_arg <- arg(args, "generated-cache-dir", "")
  generated_cache_dir <- if (nzchar(generated_cache_arg)) {
    resolve_path(generated_cache_arg, repo_root, must_work = FALSE)
  } else {
    ""
  }

  reusable <- find_reusable_generated_cache(
    intermediate_root,
    validator,
    base_dependencies,
    as.character(config$raw_data$seurat_rds_sha256),
    reconstruction_contract,
    generator,
    explicit_rds,
    deposited_rds,
    upstream_dir,
    cellranger_root,
    sample_info,
    generated_cache_dir
  )
  if (!is.null(reusable) && is.null(metrics_path)) {
    message(
      "[si-figures] Reusing validated generated cache: ",
      reusable$build_dir
    )
    render_cache(
      renderer,
      file.path(reusable$build_dir, "tables"),
      config_path,
      output_dir,
      TRUE,
      "01_render.log",
      file.path(
        reusable$build_dir,
        "metadata",
        "input_manifest.tsv"
      )
    )
    return(invisible(output_dir))
  }

  workers <- suppressWarnings(as.numeric(arg(args, "workers", "1")))
  if (length(workers) != 1L || !is.finite(workers) ||
      workers < 1L || workers > 64L || workers != floor(workers)) {
    stop("--workers must be one integer from 1 to 64", call. = FALSE)
  }
  workers <- as.integer(workers)
  raw_manifest <- resolve_path(
    as.character(config$raw_data$required_manifest),
    repo_root,
    must_work = FALSE
  )
  selection <- figure7_select_seurat_source(
    needs_seurat = TRUE,
    explicit_rds = explicit_rds,
    deposited_rds = deposited_rds,
    deposited_sha256 = as.character(config$raw_data$seurat_rds_sha256),
    upstream_dir = upstream_dir,
    cellranger_root = cellranger_root,
    module_dir = figure7_dir,
    environment_lock = environment_lock,
    config = config,
    all_ploidy = all_ploidy,
    sample_info = sample_info,
    required_source_kind = if (is.null(reusable)) {
      "auto"
    } else if (isTRUE(reusable$reconstructed)) {
      "manifested_reconstruction"
    } else {
      "deposited"
    }
  )
  selection <- materialize_seurat_selection(
    selection,
    script_dir,
    environment_lock,
    config,
    config_path,
    all_ploidy,
    sample_info,
    workers,
    raw_data_dir,
    raw_manifest,
    flag(args, "download-missing-raw", TRUE),
    intermediate_root
  )
  source_dependencies <- figure7_seurat_selection_dependency_values(
    selection
  )
  expected <- computational_dependencies(
    si_raw_scientific_code_contract_sha256,
    si_raw_config_contract_sha256,
    si_raw_runtime_contract_sha256,
    all_ploidy_sha256,
    source_dependencies,
    audit_optional_scvelo = metrics_path
  )

  reconstructed <- nzchar(selection$reconstruction_manifest)
  forwarded_upstream <- if (reconstructed) selection$upstream_dir else NULL
  forwarded_sample_info <- if (reconstructed) sample_info else NULL
  forwarded_cellranger <- if (
    reconstructed &&
      nzchar(selection$cellranger_root)
  ) {
    selection$cellranger_root
  } else {
    NULL
  }

  if (!is.null(reusable)) {
    if (!identical(
      unname(reusable$dependencies[names(expected)]),
      unname(tolower(as.character(expected)))
    ) || !setequal(names(reusable$dependencies), names(expected))) {
      stop(
        "Optional scVelo audit cannot recreate the Seurat source recorded ",
        "by the reusable SI cache",
        call. = FALSE
      )
    }
    run_stage(
      "Audit optional scVelo metrics against selected Seurat metadata",
      builder,
      raw_audit_args(
        all_ploidy = all_ploidy,
        seurat_rds = selection$rds,
        config = config_path,
        scvelo_metrics = metrics_path,
        seurat_upstream_dir = forwarded_upstream,
        sample_info = forwarded_sample_info,
        cellranger_root = forwarded_cellranger
      ),
      file.path(
        intermediate_root,
        paste0(basename(reusable$build_dir), ".audit.log")
      )
    )
    audited_manifest <- file.path(
      intermediate_root,
      "render_manifests",
      paste0(
        basename(reusable$build_dir),
        ".scvelo_",
        substr(sha256(metrics_path), 1L, 20L),
        ".tsv"
      )
    )
    augment_manifest_with_scvelo_audit(
      file.path(
        reusable$build_dir,
        "metadata",
        "input_manifest.tsv"
      ),
      metrics_path,
      audited_manifest,
      repo_root
    )
    render_cache(
      renderer,
      file.path(reusable$build_dir, "tables"),
      config_path,
      output_dir,
      TRUE,
      "01_render.log",
      audited_manifest
    )
    return(invisible(output_dir))
  }

  fingerprint <- computational_fingerprint(expected)
  build_dir <- file.path(
    intermediate_root,
    paste0("generated_human_only_", substr(fingerprint, 1L, 20L))
  )
  work_cache_dir <- paste0(build_dir, ".work")
  if (dir.exists(build_dir)) {
    quarantine_generated(build_dir, intermediate_root)
  }
  staging <- paste0(build_dir, ".tmp.", Sys.getpid())
  if (dir.exists(staging)) {
    quarantine_generated(staging, intermediate_root)
  }
  run_stage(
    "Build raw Supplementary Figure plot tables",
    builder,
    raw_builder_args(
      all_ploidy = all_ploidy,
      seurat_rds = selection$rds,
      config = config_path,
      output = staging,
      seurat_upstream_dir = forwarded_upstream,
      sample_info = forwarded_sample_info,
      cellranger_root = forwarded_cellranger,
      scvelo_metrics = metrics_path,
      work_cache_dir = work_cache_dir,
      work_dependency_fingerprint = fingerprint,
      workers = workers
    ),
    file.path(
      intermediate_root,
      paste0(basename(build_dir), ".log")
    )
  )
  if (!validate_cache(
    validator,
    file.path(staging, "tables"),
    policy = "generated-human-only"
  )) {
    stop("Generated SI plot-facing cache failed validation", call. = FALSE)
  }
  staged_manifest <- read_manifest(
    file.path(staging, "metadata", "input_manifest.tsv")
  )
  staged_dependencies <- manifest_computational_dependencies(
    staged_manifest
  )
  if (is.null(staged_dependencies) ||
      !setequal(names(staged_dependencies), names(expected)) ||
      !identical(
        unname(staged_dependencies[names(expected)]),
        unname(tolower(as.character(expected)))
      )) {
    stop(
      "Generated SI build has the wrong computational lineage",
      call. = FALSE
    )
  }
  staged_candidate <- si_cache_candidate_dependencies(
    staging,
    base_dependencies,
    as.character(config$raw_data$seurat_rds_sha256),
    reconstruction_contract,
    generator,
    explicit_rds = selection$rds,
    deposited_rds = deposited_rds,
    upstream_dir = if (reconstructed) selection$upstream_dir else "",
    cellranger_root = if (reconstructed) {
      selection$cellranger_root
    } else {
      ""
    },
    sample_info = sample_info,
    require_fingerprint_name = FALSE
  )
  if (is.null(staged_candidate) ||
      !identical(staged_candidate$fingerprint, fingerprint)) {
    stop(
      "Generated SI build has contradictory or incomplete provenance",
      call. = FALSE
    )
  }
  if (!file.rename(staging, build_dir)) {
    stop("Cannot atomically promote generated SI cache", call. = FALSE)
  }

  render_cache(
    renderer,
    file.path(build_dir, "tables"),
    config_path,
    output_dir,
    TRUE,
    "01_render.log",
    file.path(build_dir, "metadata", "input_manifest.tsv")
  )
  invisible(output_dir)
}

if (sys.nframe() == 0L) main()
