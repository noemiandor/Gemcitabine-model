.ptpv4_source_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(e) ""
)
.ptpv4_script_dir <- if (nzchar(.ptpv4_source_file)) {
  dirname(.ptpv4_source_file)
} else {
  file.path(getwd(), "Code/in-vivo")
}

source(file.path(
  .ptpv4_script_dir,
  "04j_pseudotime_treatment_ploidy_programs_v3_util.R"
))

ptpv4_thread_env <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
}

ptpv4_now <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
}

ptpv4_message <- function(..., log_file = NULL) {
  msg <- paste0("[", ptpv4_now(), "] ", paste0(..., collapse = ""))
  message(msg)
  if (!is.null(log_file)) {
    cat(msg, "\n", file = log_file, append = TRUE)
  }
  invisible(msg)
}

ptpv4_bind_rows <- function(rows) {
  rows <- rows[!vapply(
    rows,
    function(x) is.null(x) || nrow(as.data.frame(x)) == 0L,
    logical(1L)
  )]
  if (length(rows) == 0L) {
    return(data.frame())
  }
  rows <- lapply(rows, as.data.frame, stringsAsFactors = FALSE)
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    for (nm in setdiff(all_names, names(x))) {
      x[[nm]] <- NA
    }
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, rows)
}

ptpv4_merge_lists <- function(x, y) {
  if (!is.list(x) || !is.list(y)) {
    return(y)
  }
  out <- x
  for (nm in names(y)) {
    if (
      nm %in% names(out) &&
      is.list(out[[nm]]) &&
      is.list(y[[nm]]) &&
      !is.null(names(out[[nm]])) &&
      !is.null(names(y[[nm]]))
    ) {
      out[[nm]] <- ptpv4_merge_lists(out[[nm]], y[[nm]])
    } else {
      out[[nm]] <- y[[nm]]
    }
  }
  out
}

ptpv4_read_config <- function(path) {
  v4 <- yaml::read_yaml(path)
  base_path <- v4$analysis$base_v3_config
  if (is.null(base_path) || !file.exists(base_path)) {
    stop("V4 config requires an existing analysis.base_v3_config.", call. = FALSE)
  }
  base <- ptpv3_read_config(base_path)
  cfg <- ptpv4_merge_lists(base, v4)
  cfg$v4 <- v4
  cfg$analysis <- ptpv4_merge_lists(base$analysis, v4$analysis)
  cfg$simulation <- v4$simulation
  cfg$v2_methods$self_contained$mroast_rotations_observed <-
    as.integer(v4$self_contained$primary_mroast_rotations)
  cfg$v2_methods$self_contained$mroast_seed <-
    as.integer(v4$reproducibility$rotation_seed)
  cfg
}

ptpv4_validate_config <- function(cfg) {
  required <- list(
    analysis = c(
      "name", "version", "schema_version", "method_version",
      "v1_result_root", "v2_result_root", "v3_result_root", "v4_result_root"
    ),
    reproducibility = c(
      "master_seed", "permutation_seed", "simulation_seed",
      "crossfit_seed", "coordinate_seed", "dose_seed"
    ),
    exact_permutation = c(
      "treated_per_ploidy", "expected_unique_assignments"
    ),
    simulation = c("null_replicates", "power_replicates", "scenarios"),
    multiscale = c("n_bins", "window_widths"),
    trajectory_global = c("initial_bins", "spline_df_candidates"),
    dose = c("expected_assignments", "candidates")
  )
  rows <- list()
  for (section in names(required)) {
    sec <- cfg[[section]]
    for (field in required[[section]]) {
      ok <- !is.null(sec) && !is.null(sec[[field]])
      rows[[length(rows) + 1L]] <- data.frame(
        section = section,
        field = field,
        status = if (ok) "ok" else "missing",
        stringsAsFactors = FALSE
      )
    }
  }
  out <- ptpv4_bind_rows(rows)
  if (any(out$status != "ok")) {
    stop(
      "V4 configuration schema is incomplete: ",
      paste(
        paste0(out$section[out$status != "ok"], ".", out$field[out$status != "ok"]),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  out
}

ptpv4_dirs <- function(output_root) {
  d <- function(...) ptp_ensure_dir(file.path(output_root, ...))
  primary_root <- d("primary_initial_ploidy")
  list(
    root = ptp_ensure_dir(output_root),
    manifest = d("00_manifest"),
    frozen = d("frozen_definitions"),
    qc = d("qc"),
    simulation = d("simulation"),
    primary = list(
      root = primary_root,
      score_models = d("primary_initial_ploidy", "score_models"),
      fry_mroast = d("primary_initial_ploidy", "fry_mroast"),
      camera = d("primary_initial_ploidy", "camera_diagnostic"),
      mouse_gene = d(
        "primary_initial_ploidy",
        "mouse_level_gene_statistics"
      ),
      adaptive = d("primary_initial_ploidy", "adaptive_gene_set"),
      family = d(
        "primary_initial_ploidy",
        "response_family_hierarchical"
      ),
      multiscale = d("primary_initial_ploidy", "multiscale_pseudotime"),
      trajectory = d("primary_initial_ploidy", "trajectory_global"),
      robustness = d("primary_initial_ploidy", "robustness")
    ),
    marginal = d("secondary_marginal_total_effect"),
    dose = d("secondary_dose"),
    etp = d("secondary_etp"),
    end_time = d("secondary_end_time_ploidy"),
    coordinates = d("coordinate_sensitivity"),
    crossfit = d("crossfit_04i"),
    figures = d("figures"),
    tables = d("tables"),
    report = d("report"),
    logs = d("logs"),
    checkpoints = d("checkpoints")
  )
}

ptpv4_prepare_output <- function(path, overwrite = FALSE) {
  preserved_names <- c(
    "source_checksums_before.csv",
    "v1_result_checksums_before.csv",
    "v2_result_checksums_before.csv",
    "v3_result_checksums_before.csv"
  )
  preserved <- list()
  manifest <- file.path(path, "00_manifest")
  if (dir.exists(manifest)) {
    for (nm in preserved_names) {
      p <- file.path(manifest, nm)
      if (file.exists(p)) {
        preserved[[nm]] <- readBin(p, what = "raw", n = file.info(p)$size)
      }
    }
  }
  checkpoint_dir <- file.path(path, "checkpoints")
  resumable <- dir.exists(checkpoint_dir) &&
    length(list.files(
      checkpoint_dir,
      pattern = "\\.rds$",
      full.names = TRUE
    )) > 0L
  if (dir.exists(path) && isTRUE(overwrite) && !resumable) {
    unlink(path, recursive = TRUE, force = TRUE)
  } else if (dir.exists(path) && isTRUE(overwrite) && resumable) {
    message(
      "Preserving existing V4 checkpoint tree for schema-validated resume: ",
      checkpoint_dir
    )
  } else if (dir.exists(path)) {
    existing <- list.files(path, all.files = TRUE, no.. = TRUE)
    allowed <- identical(sort(existing), "00_manifest") &&
      all(list.files(manifest, all.files = TRUE, no.. = TRUE) %in% preserved_names)
    if (length(existing) > 0L && !allowed) {
      stop(
        "V4 output root is not empty. Use --overwrite=TRUE: ",
        path,
        call. = FALSE
      )
    }
  }
  dir.create(file.path(path, "00_manifest"), recursive = TRUE, showWarnings = FALSE)
  for (nm in names(preserved)) {
    writeBin(preserved[[nm]], file.path(path, "00_manifest", nm))
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

ptpv4_file_sha <- function(path) {
  if (!file.exists(path) || dir.exists(path)) {
    return(NA_character_)
  }
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

ptpv4_object_sha <- function(x) {
  digest::digest(x, algo = "sha256")
}

ptpv4_config_sha <- function(args, cfg) {
  source_dir <- dirname(args$config)
  ptpv4_object_sha(list(
    config_file_sha = ptpv4_file_sha(args$config),
    entry_script_sha = ptpv4_file_sha(file.path(
      source_dir,
      "04j_pseudotime_treatment_ploidy_programs_v4.R"
    )),
    utility_script_sha = ptpv4_file_sha(file.path(
      source_dir,
      "04j_pseudotime_treatment_ploidy_programs_v4_util.R"
    )),
    config = cfg$v4,
    method_version = cfg$analysis$method_version,
    schema_version = cfg$analysis$schema_version
  ))
}

ptpv4_checkpoint_path <- function(dirs, id) {
  file.path(dirs$checkpoints, paste0(ptp_safe_name(id), ".rds"))
}

ptpv4_checkpoint_read <- function(
  dirs,
  id,
  cfg_sha,
  cfg,
  input_sha,
  required_fields = character(0L)
) {
  path <- ptpv4_checkpoint_path(dirs, id)
  if (!file.exists(path)) {
    return(NULL)
  }
  x <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(x) || !is.list(x)) {
    return(NULL)
  }
  ok <- identical(x$config_checksum, cfg_sha) &&
    identical(x$method_version, cfg$analysis$method_version) &&
    identical(x$schema_version, cfg$checkpoint$schema_version) &&
    identical(x$input_checksum, input_sha) &&
    identical(x$status, "complete") &&
    all(required_fields %in% names(x$data))
  if (!ok) {
    return(NULL)
  }
  x$data
}

ptpv4_checkpoint_write <- function(dirs, id, data, cfg_sha, cfg, input_sha) {
  path <- ptpv4_checkpoint_path(dirs, id)
  tmp <- paste0(path, ".tmp")
  saveRDS(
    list(
      checkpoint_id = id,
      created_at = ptpv4_now(),
      config_checksum = cfg_sha,
      method_version = cfg$analysis$method_version,
      schema_version = cfg$checkpoint$schema_version,
      input_checksum = input_sha,
      status = "complete",
      data = data
    ),
    tmp,
    compress = "xz"
  )
  if (!file.rename(tmp, path)) {
    stop("Could not finalize checkpoint: ", path, call. = FALSE)
  }
  invisible(path)
}

ptpv4_checkpoint_compute <- function(
  dirs,
  id,
  cfg_sha,
  cfg,
  input_sha,
  compute,
  required_fields = character(0L)
) {
  checkpoint <- ptpv4_checkpoint_read(
    dirs,
    id,
    cfg_sha,
    cfg,
    input_sha,
    required_fields = required_fields
  )
  if (!is.null(checkpoint)) {
    return(checkpoint)
  }
  value <- compute()
  ptpv4_checkpoint_write(
    dirs,
    id,
    value,
    cfg_sha,
    cfg,
    input_sha
  )
  value
}

ptpv4_checksum_tree <- function(path, scope, phase) {
  if (!dir.exists(path)) {
    return(data.frame(
      scope = scope,
      phase = phase,
      path = path,
      relative_path = "",
      bytes = NA_real_,
      sha256 = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  files <- list.files(
    path,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE
  )
  files <- sort(files[!file.info(files)$isdir])
  root_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  file_norm <- normalizePath(files, winslash = "/", mustWork = FALSE)
  data.frame(
    scope = scope,
    phase = phase,
    path = file_norm,
    relative_path = substring(file_norm, nchar(root_norm) + 2L),
    bytes = file.info(files)$size,
    sha256 = vapply(files, ptpv4_file_sha, character(1L)),
    stringsAsFactors = FALSE
  )
}

ptpv4_compare_checksum_tables <- function(before, after) {
  b <- before[, c("scope", "relative_path", "sha256"), drop = FALSE]
  names(b)[3L] <- "sha256_before"
  a <- after[, c("scope", "relative_path", "sha256"), drop = FALSE]
  names(a)[3L] <- "sha256_after"
  out <- merge(
    b,
    a,
    by = c("scope", "relative_path"),
    all = TRUE,
    sort = FALSE
  )
  out$status <- ifelse(
    is.na(out$sha256_before),
    "added",
    ifelse(
      is.na(out$sha256_after),
      "removed",
      ifelse(out$sha256_before == out$sha256_after, "unchanged", "changed")
    )
  )
  out
}

ptpv4_record_input_manifest <- function(args, cfg, repo_root, dirs, seurat_rds) {
  inputs <- data.frame(
    input = c(
      "cell_metadata",
      "noncell_metadata",
      "sample_info",
      "seurat_rds",
      "v4_config",
      "v3_config",
      "plan"
    ),
    path = c(
      args$cell_metadata,
      args$noncell_metadata,
      args$sample_info,
      seurat_rds,
      args$config,
      cfg$analysis$base_v3_config,
      cfg$analysis$plan_path
    ),
    stringsAsFactors = FALSE
  )
  inputs$path <- normalizePath(
    inputs$path,
    winslash = "/",
    mustWork = FALSE
  )
  inputs$exists <- file.exists(inputs$path)
  inputs$bytes <- suppressWarnings(file.info(inputs$path)$size)
  inputs$sha256 <- vapply(inputs$path, ptpv4_file_sha, character(1L))
  v3_path <- cfg$analysis$v3_manifest_path
  v3 <- if (file.exists(v3_path)) {
    read.csv(v3_path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  inputs$v3_manifest_sha256 <- NA_character_
  inputs$v3_manifest_input <- NA_character_
  if (nrow(v3) > 0L) {
    for (i in seq_len(nrow(inputs))) {
      hit <- which(
        normalizePath(v3$path, winslash = "/", mustWork = FALSE) ==
          inputs$path[[i]]
      )
      if (length(hit) == 1L) {
        inputs$v3_manifest_sha256[[i]] <- v3$sha256[[hit]]
        inputs$v3_manifest_input[[i]] <- v3$input[[hit]]
      }
    }
  }
  inputs$v3_comparison_status <- ifelse(
    !inputs$exists,
    "missing_current_input",
    ifelse(
      is.na(inputs$v3_manifest_sha256),
      "not_recorded_in_v3_manifest",
      ifelse(
        inputs$sha256 == inputs$v3_manifest_sha256,
        "matched_v3_manifest",
        "checksum_mismatch_v3_manifest"
      )
    )
  )
  ptp_write_csv(inputs, file.path(dirs$manifest, "input_manifest.csv"))
  ptp_write_csv(
    inputs[, c(
      "input", "path", "sha256", "v3_manifest_input",
      "v3_manifest_sha256", "v3_comparison_status"
    )],
    file.path(dirs$manifest, "input_checksum_v3_comparison.csv")
  )
  inputs
}

ptpv4_program_role_table <- function(programs) {
  role <- ptpv3_program_role_table(programs)
  role$analysis_role[role$tier == "negative_control"] <- "matched_null_control"
  role$directional_eligibility <-
    role$analysis_role == "mechanistic_directional_signature"
  role$directional_claim_eligible <- role$directional_eligibility
  role
}

ptpv4_membership_sha <- function(genes) {
  genes <- sort(unique(ptp_clean_gene_symbols(genes)))
  digest::digest(
    paste(unname(genes), collapse = "\n"),
    algo = "sha256",
    serialize = FALSE
  )
}

ptpv4_membership_objects <- function(programs, role_table) {
  hashes <- vapply(programs, function(p) {
    ptpv4_membership_sha(p$observed_genes)
  }, character(1L))
  alias_rows <- lapply(names(programs), function(id) {
    r <- role_table[match(id, role_table$program_id), , drop = FALSE]
    data.frame(
      program_id = id,
      membership_id = hashes[[id]],
      program_label = r$program_label,
      tier = r$tier,
      tier_label = r$tier_label,
      response_family = r$response_family,
      response_family_label = r$response_family_label,
      analysis_role = r$analysis_role,
      directional_eligibility = r$directional_eligibility,
      control_status = r$control_status,
      n_genes = length(unique(programs[[id]]$observed_genes)),
      stringsAsFactors = FALSE
    )
  })
  alias <- ptpv4_bind_rows(alias_rows)
  split_ids <- split(alias$program_id, alias$membership_id)
  memberships <- lapply(names(split_ids), function(mid) {
    ids <- sort(split_ids[[mid]])
    genes <- sort(unique(ptp_clean_gene_symbols(programs[[ids[[1L]]]]$observed_genes)))
    directional_ids <- ids[
      role_table$directional_eligibility[
        match(ids, role_table$program_id)
      ] %in% TRUE
    ]
    directional_id <- if (length(directional_ids) == 1L) {
      directional_ids[[1L]]
    } else {
      NA_character_
    }
    weights <- NULL
    if (!is.na(directional_id)) {
      w <- programs[[directional_id]]$observed_gene_weights
      w$gene <- ptp_clean_gene_symbols(w$gene)
      w <- w[match(genes, w$gene), c("gene", "signed_weight"), drop = FALSE]
      weights <- w
    }
    list(
      id = mid,
      aliases = ids,
      genes = genes,
      directional_program_id = directional_id,
      directional_weights = weights
    )
  })
  names(memberships) <- names(split_ids)
  membership_table <- ptpv4_bind_rows(lapply(memberships, function(m) {
    a <- alias[alias$membership_id == m$id, , drop = FALSE]
    data.frame(
      membership_id = m$id,
      membership_short_id = paste0("m_", substr(m$id, 1L, 16L)),
      n_genes = length(m$genes),
      n_aliases = nrow(a),
      aliases = paste(a$program_id, collapse = ";"),
      tiers = paste(sort(unique(a$tier)), collapse = ";"),
      response_families = paste(
        sort(unique(a$response_family)),
        collapse = ";"
      ),
      directional_eligible = !is.na(m$directional_program_id),
      directional_program_id = m$directional_program_id,
      genes = paste(m$genes, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  pair_rows <- list()
  mids <- names(memberships)
  for (i in seq_along(mids)) {
    if (i == length(mids)) {
      next
    }
    for (j in seq.int(i + 1L, length(mids))) {
      a <- memberships[[mids[[i]]]]$genes
      b <- memberships[[mids[[j]]]]$genes
      ni <- length(intersect(a, b))
      nu <- length(union(a, b))
      pair_rows[[length(pair_rows) + 1L]] <- data.frame(
        membership_id_a = mids[[i]],
        membership_id_b = mids[[j]],
        n_intersection = ni,
        n_union = nu,
        jaccard = if (nu > 0L) ni / nu else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  all_genes <- unlist(lapply(memberships, `[[`, "genes"), use.names = FALSE)
  multiplicity <- sort(table(all_genes), decreasing = TRUE)
  list(
    memberships = memberships,
    alias = alias,
    membership_table = membership_table,
    pairwise_jaccard = ptpv4_bind_rows(pair_rows),
    gene_multiplicity = data.frame(
      gene_symbol = names(multiplicity),
      n_unique_memberships = as.integer(multiplicity),
      stringsAsFactors = FALSE
    )
  )
}

ptpv4_write_frozen_universe <- function(
  programs,
  role_table,
  membership,
  cfg,
  dirs
) {
  ptp_write_csv(
    role_table,
    file.path(dirs$frozen, "program_universe_role_eligibility.csv")
  )
  ptp_write_csv(
    ptp_program_weight_table(programs),
    file.path(dirs$frozen, "program_score_gene_weights.csv")
  )
  ptp_write_csv(
    membership$alias,
    file.path(dirs$frozen, "program_label_to_membership_alias.csv")
  )
  ptp_write_csv(
    membership$membership_table,
    file.path(dirs$frozen, "unique_memberships.csv")
  )
  ptp_write_csv(
    membership$pairwise_jaccard,
    file.path(dirs$frozen, "unique_membership_pairwise_jaccard.csv")
  )
  ptp_write_csv(
    membership$gene_multiplicity,
    file.path(dirs$frozen, "unique_membership_gene_multiplicity.csv")
  )
  family_rows <- lapply(split(
    membership$alias,
    membership$alias$response_family
  ), function(x) {
    mids <- unique(x$membership_id)
    genes <- sort(unique(unlist(
      lapply(membership$memberships[mids], `[[`, "genes"),
      use.names = FALSE
    )))
    data.frame(
      response_family = x$response_family[[1L]],
      response_family_label = x$response_family_label[[1L]],
      tiers = paste(sort(unique(x$tier)), collapse = ";"),
      n_labels = nrow(x),
      n_unique_memberships = length(mids),
      n_union_genes = length(genes),
      union_membership_id = ptpv4_membership_sha(genes),
      membership_ids = paste(sort(mids), collapse = ";"),
      union_genes = paste(genes, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  family <- ptpv4_bind_rows(family_rows)
  ptp_write_csv(
    family,
    file.path(dirs$frozen, "response_family_definitions.csv")
  )
  ptp_write_csv(
    data.frame(
      branch = "external_signature",
      status = cfg$analysis$external_signature_status,
      network_used = FALSE,
      external_cohort_used = FALSE,
      stringsAsFactors = FALSE
    ),
    file.path(dirs$frozen, "external_signature_status.csv")
  )
  ptp_write_csv(
    data.frame(
      policy = "species_handling",
      status = cfg$analysis$species_policy,
      species_filter_added = FALSE,
      species_split_added = FALSE,
      species_ratio_qc_added = FALSE,
      stringsAsFactors = FALSE
    ),
    file.path(dirs$frozen, "species_policy_audit.csv")
  )
  family
}

ptpv4_compare_program_universe_v3 <- function(programs, membership, cfg, dirs) {
  v3_weights_path <- file.path(
    cfg$analysis$v3_result_root,
    "frozen_definitions",
    "program_score_gene_weights.csv"
  )
  if (!file.exists(v3_weights_path)) {
    out <- data.frame(
      program_id = names(programs),
      current_membership_id = membership$alias$membership_id[
        match(names(programs), membership$alias$program_id)
      ],
      v3_membership_id = NA_character_,
      status = "v3_weight_file_missing",
      stringsAsFactors = FALSE
    )
    ptp_write_csv(
      out,
      file.path(dirs$manifest, "program_universe_v3_comparison.csv")
    )
    return(out)
  }
  v3 <- read.csv(v3_weights_path, check.names = FALSE, stringsAsFactors = FALSE)
  if ("observed" %in% names(v3)) {
    observed <- tolower(trimws(as.character(v3$observed))) %in%
      c("true", "t", "1")
    v3 <- v3[observed, , drop = FALSE]
  }
  v3$gene <- ptp_clean_gene_symbols(v3$gene)
  v3_hash <- vapply(split(v3$gene, v3$program_id), function(g) {
    ptpv4_membership_sha(g)
  }, character(1L))
  current <- membership$alias[, c("program_id", "membership_id"), drop = FALSE]
  names(current)[2L] <- "current_membership_id"
  current$v3_membership_id <- unname(v3_hash[current$program_id])
  current$status <- ifelse(
    is.na(current$v3_membership_id),
    "missing_in_v3",
    ifelse(
      current$current_membership_id == current$v3_membership_id,
      "unchanged",
      "changed"
    )
  )
  ptp_write_csv(
    current,
    file.path(dirs$manifest, "program_universe_v3_comparison.csv")
  )
  current
}

ptpv4_build_assignment_space <- function(mouse_meta, treated_per_ploidy = 4L) {
  mouse_meta <- mouse_meta[order(
    mouse_meta$initial_ploidy,
    mouse_meta$sample_id
  ), , drop = FALSE]
  strata <- split(mouse_meta$sample_id, mouse_meta$initial_ploidy)
  if (length(strata) != 2L || any(lengths(strata) != 8L)) {
    stop("Exact assignment space requires two strata of eight mice.", call. = FALSE)
  }
  strata <- lapply(strata, sort)
  combo <- lapply(strata, function(ids) {
    m <- utils::combn(ids, treated_per_ploidy)
    out <- matrix(FALSE, nrow = length(ids), ncol = ncol(m))
    rownames(out) <- ids
    for (j in seq_len(ncol(m))) {
      out[, j] <- ids %in% m[, j]
    }
    out
  })
  names(combo) <- names(strata)
  n1 <- ncol(combo[[1L]])
  n2 <- ncol(combo[[2L]])
  rows <- vector("list", n1 * n2)
  k <- 0L
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      k <- k + 1L
      treated <- sort(c(
        rownames(combo[[1L]])[combo[[1L]][, i]],
        rownames(combo[[2L]])[combo[[2L]][, j]]
      ))
      rows[[k]] <- data.frame(
        assignment_index = k,
        assignment_id = sprintf("perm_%04d", k),
        stratum1_combination = i,
        stratum2_combination = j,
        treated_samples = paste(treated, collapse = ";"),
        assignment_checksum = ptpv3_assignment_checksum(
          mouse_meta$sample_id,
          treated
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  table <- do.call(rbind, rows)
  observed_treated <- sort(unique(
    mouse_meta$sample_id[mouse_meta$treatment != "control"]
  ))
  obs_sha <- ptpv3_assignment_checksum(mouse_meta$sample_id, observed_treated)
  table$is_observed_assignment <- table$assignment_checksum == obs_sha
  list(
    table = table,
    strata = strata,
    combinations = combo,
    observed_index = which(table$is_observed_assignment),
    observed_treated = observed_treated
  )
}

ptpv4_build_mouse_standardized_expression <- function(
  primary_fit,
  eligible_subbins,
  dirs
) {
  meta <- primary_fit$metadata
  E <- ptp_logcpm_matrix(primary_fit)
  eligible_subbins <- as.character(eligible_subbins)
  sample_ids <- sort(unique(meta$sample_id))
  expected <- expand.grid(
    sample_id = sample_ids,
    subbin_id = eligible_subbins,
    stringsAsFactors = FALSE
  )
  expected$base_weight <- 1 / length(eligible_subbins)
  key <- paste(meta$sample_id, meta$subbin_id, sep = "__")
  expected$key <- paste(expected$sample_id, expected$subbin_id, sep = "__")
  expected$observation_index <- match(expected$key, key)
  expected$available <- is.finite(expected$observation_index)
  expected$availability_reason <- ifelse(
    expected$available,
    "retained_original_minimum_cell_rule",
    "missing_or_failed_original_minimum_cell_rule"
  )
  available_weight <- aggregate(
    base_weight ~ sample_id,
    data = expected[expected$available, , drop = FALSE],
    FUN = sum
  )
  names(available_weight)[2L] <- "available_base_weight_sum"
  expected <- merge(
    expected,
    available_weight,
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )
  expected$final_weight <- ifelse(
    expected$available,
    expected$base_weight / expected$available_base_weight_sum,
    0
  )
  expected <- expected[order(
    match(expected$sample_id, sample_ids),
    match(expected$subbin_id, eligible_subbins)
  ), , drop = FALSE]
  operator <- matrix(
    0,
    nrow = length(sample_ids),
    ncol = nrow(meta),
    dimnames = list(sample_ids, meta$sample_subbin_id)
  )
  for (i in seq_len(nrow(expected))) {
    if (expected$available[[i]]) {
      operator[
        expected$sample_id[[i]],
        expected$observation_index[[i]]
      ] <- expected$final_weight[[i]]
    }
  }
  subbin_means <- vapply(eligible_subbins, function(sb) {
    idx <- which(meta$subbin_id == sb)
    if (!length(idx)) {
      return(rep(NA_real_, nrow(E)))
    }
    rowMeans(E[, idx, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(E)))
  colnames(subbin_means) <- eligible_subbins
  residual <- E
  for (j in seq_len(ncol(E))) {
    residual[, j] <- E[, j] - subbin_means[, meta$subbin_id[[j]]]
  }
  mouse_expression <- residual %*% t(operator)
  colnames(mouse_expression) <- sample_ids
  mouse_meta <- unique(meta[, c(
    "sample_id",
    "initial_ploidy",
    "treatment",
    "dose_mg"
  ), drop = FALSE])
  mouse_meta <- mouse_meta[match(sample_ids, mouse_meta$sample_id), , drop = FALSE]
  operator_long <- expected[, c(
    "sample_id",
    "subbin_id",
    "base_weight",
    "available",
    "availability_reason",
    "available_base_weight_sum",
    "final_weight",
    "observation_index"
  )]
  operator_long$standardization_rule <-
    "equal_base_weights_then_support_only_renormalization"
  operator_long$operator_checksum <- ptpv4_object_sha(operator)
  ptp_write_csv(
    operator_long,
    file.path(dirs$primary$mouse_gene, "mouse_subbin_weights.csv")
  )
  ptp_write_csv(
    data.frame(
      sample_id = rownames(operator),
      operator,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    file.path(
      dirs$primary$mouse_gene,
      "final_standardization_operator.csv"
    )
  )
  ptp_write_csv(
    mouse_meta,
    file.path(dirs$primary$mouse_gene, "mouse_metadata.csv")
  )
  ptp_write_csv(
    data.frame(
      gene = rownames(subbin_means),
      gene_symbol = ptp_clean_gene_symbols(rownames(subbin_means)),
      subbin_means,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    file.path(
      dirs$primary$mouse_gene,
      "treatment_blind_subbin_nuisance_means.csv"
    )
  )
  list(
    expression = mouse_expression,
    metadata = mouse_meta,
    operator = operator,
    operator_long = operator_long,
    nuisance_means = subbin_means,
    residual_subbin_expression = residual
  )
}

ptpv4_stratum_stats <- function(Y, A) {
  n_t <- colSums(A)
  n_c <- nrow(A) - n_t
  mean_t <- sweep(Y %*% A, 2L, n_t, "/")
  mean_c <- sweep(Y %*% (!A), 2L, n_c, "/")
  ss_t <- Y^2 %*% A - sweep(mean_t^2, 2L, n_t, "*")
  ss_c <- Y^2 %*% (!A) - sweep(mean_c^2, 2L, n_c, "*")
  var_t <- sweep(ss_t, 2L, pmax(n_t - 1L, 1L), "/")
  var_c <- sweep(ss_c, 2L, pmax(n_c - 1L, 1L), "/")
  var_contrast <- sweep(var_t, 2L, n_t, "/") +
    sweep(var_c, 2L, n_c, "/")
  list(
    difference = mean_t - mean_c,
    variance = pmax(var_contrast, 0),
    mean_treated = mean_t,
    mean_control = mean_c
  )
}

ptpv4_interaction_stats <- function(Y, assignments) {
  s1 <- assignments$strata[[1L]]
  s2 <- assignments$strata[[2L]]
  Y1 <- Y[, s1, drop = FALSE]
  Y2 <- Y[, s2, drop = FALSE]
  z1 <- ptpv4_stratum_stats(Y1, assignments$combinations[[1L]])
  z2 <- ptpv4_stratum_stats(Y2, assignments$combinations[[2L]])
  B <- nrow(assignments$table)
  estimate <- matrix(NA_real_, nrow = nrow(Y), ncol = B)
  se <- matrix(NA_real_, nrow = nrow(Y), ncol = B)
  k <- 0L
  for (i in seq_len(ncol(z1$difference))) {
    cols <- k + seq_len(ncol(z2$difference))
    estimate[, cols] <- sweep(
      z2$difference,
      1L,
      z1$difference[, i],
      "-"
    )
    se[, cols] <- sqrt(sweep(
      z2$variance,
      1L,
      z1$variance[, i],
      "+"
    ))
    k <- max(cols)
  }
  tstat <- estimate / se
  tstat[!is.finite(tstat)] <- NA_real_
  list(estimate = estimate, se = se, t = tstat)
}

ptpv4_reduced_residuals <- function(mouse_expression, mouse_meta) {
  dat <- mouse_meta
  dat$initial_ploidy <- factor(dat$initial_ploidy, levels = c("2N", "4N"))
  dat$treatment <- factor(dat$treatment, levels = c("control", "gemcitabine"))
  design <- stats::model.matrix(~ initial_ploidy + treatment, data = dat)
  q <- qr(design)
  coef <- qr.coef(q, t(mouse_expression))
  fitted <- t(design %*% coef)
  residual <- mouse_expression - fitted
  list(
    design = design,
    coefficients = coef,
    fitted = fitted,
    residual = residual,
    design_columns = colnames(design)
  )
}

ptpv4_init_gene_hdf5 <- function(path, n_genes, n_assignments, cfg) {
  if (file.exists(path)) {
    unlink(path)
  }
  rhdf5::h5createFile(path)
  chunk <- c(
    min(as.integer(cfg$mouse_standardization$hdf5_gene_chunk_size), n_genes),
    min(
      as.integer(cfg$mouse_standardization$hdf5_assignment_chunk_size),
      n_assignments
    )
  )
  for (dataset in c(
    "estimate",
    "standard_error",
    "t_statistic",
    "freedman_lane_t_statistic"
  )) {
    rhdf5::h5createDataset(
      path,
      dataset,
      dims = c(n_genes, n_assignments),
      chunk = chunk,
      level = as.integer(
        cfg$mouse_standardization$hdf5_compression_level
      ),
      storage.mode = "double"
    )
  }
  rhdf5::h5createDataset(
    path,
    "gene_ids",
    dims = n_genes,
    storage.mode = "character"
  )
  rhdf5::h5createDataset(
    path,
    "assignment_ids",
    dims = n_assignments,
    storage.mode = "character"
  )
  invisible(path)
}

ptpv4_compute_gene_statistics <- function(
  mouse_bundle,
  assignments,
  cfg,
  args,
  dirs,
  cfg_sha,
  input_sha
) {
  checkpoint <- ptpv4_checkpoint_read(
    dirs,
    "primary_gene_statistics",
    cfg_sha,
    cfg,
    input_sha,
    required_fields = c("observed", "freedman_lane", "hdf5_path", "index")
  )
  if (!is.null(checkpoint) && file.exists(checkpoint$hdf5_path)) {
    return(checkpoint)
  }
  Y <- mouse_bundle$expression
  B_full <- nrow(assignments$table)
  B <- if (isTRUE(args$smoke)) {
    min(as.integer(args$smoke_permutations), B_full)
  } else {
    B_full
  }
  keep_assignments <- seq_len(B)
  if (!(assignments$observed_index %in% keep_assignments)) {
    keep_assignments[[B]] <- assignments$observed_index
  }
  assign_small <- assignments
  assign_small$table <- assignments$table[keep_assignments, , drop = FALSE]
  if (B < B_full) {
    # For smoke, construct a direct assignment list so arbitrary selected rows
    # remain valid without changing the full-run enumeration order.
    table_smoke <- assign_small$table
    combos1 <- assignments$combinations[[1L]]
    combos2 <- assignments$combinations[[2L]]
    assign_small$combinations <- list(
      combos1[, table_smoke$stratum1_combination, drop = FALSE],
      combos2[, table_smoke$stratum2_combination, drop = FALSE]
    )
    names(assign_small$combinations) <- names(assignments$combinations)
    assign_small$direct_pairing <- TRUE
  } else {
    assign_small$direct_pairing <- FALSE
  }
  h5_path <- file.path(
    dirs$primary$mouse_gene,
    "permuted_gene_statistics.h5"
  )
  reduced <- ptpv4_reduced_residuals(Y, mouse_bundle$metadata)
  chunk_size <- if (isTRUE(args$smoke)) {
    min(64L, nrow(Y))
  } else {
    as.integer(cfg$mouse_standardization$hdf5_gene_chunk_size)
  }
  chunks <- split(seq_len(nrow(Y)), ceiling(seq_len(nrow(Y)) / chunk_size))
  progress <- ptpv4_checkpoint_read(
    dirs,
    "primary_gene_statistics_progress",
    cfg_sha,
    cfg,
    input_sha,
    required_fields = c(
      "completed_chunks",
      "hdf5_path",
      "observed_rows",
      "fl_rows",
      "index_rows"
    )
  )
  if (!is.null(progress) && file.exists(progress$hdf5_path)) {
    h5_path <- progress$hdf5_path
    completed_chunks <- as.integer(progress$completed_chunks)
    observed_rows <- progress$observed_rows
    fl_rows <- progress$fl_rows
    index_rows <- progress$index_rows
    ptpv4_message(
      "Resuming primary gene statistics after ",
      length(completed_chunks),
      " completed HDF5 chunks.",
      log_file = file.path(dirs$logs, "run.log")
    )
  } else {
    ptpv4_init_gene_hdf5(h5_path, nrow(Y), B, cfg)
    completed_chunks <- integer(0L)
    observed_rows <- vector("list", length(chunks))
    fl_rows <- vector("list", length(chunks))
    index_rows <- vector("list", length(chunks))
  }
  remaining_chunks <- setdiff(seq_along(chunks), completed_chunks)
  for (ci in remaining_chunks) {
    idx <- chunks[[ci]]
    if (isTRUE(assign_small$direct_pairing)) {
      # Direct smoke computation is intentionally simple and small.
      stat <- lapply(seq_len(B), function(b) {
        treated <- strsplit(
          assign_small$table$treated_samples[[b]],
          ";",
          fixed = TRUE
        )[[1L]]
        labels <- mouse_bundle$metadata$sample_id %in% treated
        one <- ptpv4_single_assignment_stats(
          Y[idx, , drop = FALSE],
          mouse_bundle$metadata,
          labels
        )
        one
      })
      estimate <- do.call(cbind, lapply(stat, `[[`, "estimate"))
      se <- do.call(cbind, lapply(stat, `[[`, "se"))
      tstat <- estimate / se
      fl_stat <- lapply(seq_len(B), function(b) {
        treated <- strsplit(
          assign_small$table$treated_samples[[b]],
          ";",
          fixed = TRUE
        )[[1L]]
        labels <- mouse_bundle$metadata$sample_id %in% treated
        ptpv4_single_assignment_stats(
          reduced$residual[idx, , drop = FALSE],
          mouse_bundle$metadata,
          labels
        )$t
      })
      fl_t <- do.call(cbind, fl_stat)
    } else {
      stat <- ptpv4_interaction_stats(Y[idx, , drop = FALSE], assign_small)
      estimate <- stat$estimate
      se <- stat$se
      tstat <- stat$t
      fl_t <- ptpv4_interaction_stats(
        reduced$residual[idx, , drop = FALSE],
        assign_small
      )$t
    }
    row_index <- min(idx):max(idx)
    rhdf5::h5write(
      estimate,
      h5_path,
      "estimate",
      index = list(row_index, seq_len(B))
    )
    rhdf5::h5write(
      se,
      h5_path,
      "standard_error",
      index = list(row_index, seq_len(B))
    )
    rhdf5::h5write(
      tstat,
      h5_path,
      "t_statistic",
      index = list(row_index, seq_len(B))
    )
    rhdf5::h5write(
      fl_t,
      h5_path,
      "freedman_lane_t_statistic",
      index = list(row_index, seq_len(B))
    )
    obs_col <- which(assign_small$table$is_observed_assignment)
    if (!length(obs_col)) {
      obs_col <- B
    }
    p <- vapply(seq_len(nrow(tstat)), function(i) {
      (1 + sum(abs(tstat[i, ]) >= abs(tstat[i, obs_col]), na.rm = TRUE)) /
        (1 + sum(is.finite(tstat[i, ])))
    }, numeric(1L))
    p_fl <- vapply(seq_len(nrow(fl_t)), function(i) {
      (1 + sum(abs(fl_t[i, ]) >= abs(fl_t[i, obs_col]), na.rm = TRUE)) /
        (1 + sum(is.finite(fl_t[i, ])))
    }, numeric(1L))
    observed_rows[[ci]] <- data.frame(
      gene = rownames(Y)[idx],
      gene_symbol = ptp_clean_gene_symbols(rownames(Y)[idx]),
      estimate = estimate[, obs_col],
      standard_error = se[, obs_col],
      t_statistic = tstat[, obs_col],
      exact_empirical_p = p,
      stringsAsFactors = FALSE
    )
    fl_rows[[ci]] <- data.frame(
      gene = rownames(Y)[idx],
      gene_symbol = ptp_clean_gene_symbols(rownames(Y)[idx]),
      reduced_residual_observed_t = fl_t[, obs_col],
      freedman_lane_empirical_p = p_fl,
      stringsAsFactors = FALSE
    )
    index_rows[[ci]] <- data.frame(
      chunk_id = ci,
      row_start = min(idx),
      row_end = max(idx),
      n_rows = length(idx),
      hdf5_path = h5_path,
      stringsAsFactors = FALSE
    )
    ptpv4_message(
      "Gene-statistics chunk ",
      ci,
      "/",
      length(chunks),
      " written.",
      log_file = file.path(dirs$logs, "run.log")
    )
    completed_chunks <- sort(unique(c(completed_chunks, ci)))
    ptpv4_checkpoint_write(
      dirs,
      "primary_gene_statistics_progress",
      list(
        completed_chunks = completed_chunks,
        hdf5_path = h5_path,
        observed_rows = observed_rows,
        fl_rows = fl_rows,
        index_rows = index_rows
      ),
      cfg_sha,
      cfg,
      input_sha
    )
  }
  rhdf5::h5write(
    rownames(Y),
    h5_path,
    "gene_ids"
  )
  rhdf5::h5write(
    assign_small$table$assignment_id,
    h5_path,
    "assignment_ids"
  )
  observed <- ptpv4_bind_rows(observed_rows)
  observed$fdr_all_genes <- stats::p.adjust(
    observed$exact_empirical_p,
    method = "BH"
  )
  freedman_lane <- ptpv4_bind_rows(fl_rows)
  freedman_lane$fdr_all_genes <- stats::p.adjust(
    freedman_lane$freedman_lane_empirical_p,
    method = "BH"
  )
  index <- ptpv4_bind_rows(index_rows)
  ptp_write_csv(
    assign_small$table,
    file.path(dirs$primary$mouse_gene, "primary_4900_assignments.csv")
  )
  ptp_write_csv(
    observed,
    file.path(
      dirs$primary$mouse_gene,
      "observed_mouse_level_gene_statistics.csv"
    )
  )
  ptp_write_csv(
    freedman_lane,
    file.path(
      dirs$primary$mouse_gene,
      "freedman_lane_gene_statistics.csv"
    )
  )
  ptp_write_csv(
    index,
    file.path(
      dirs$primary$mouse_gene,
      "permuted_gene_statistics_hdf5_index.csv"
    )
  )
  h5_audit <- data.frame(
    artifact = h5_path,
    sha256 = ptpv4_file_sha(h5_path),
    n_genes = nrow(Y),
    n_assignments = B,
    datasets = paste(
      c(
        "estimate",
        "standard_error",
        "t_statistic",
        "freedman_lane_t_statistic",
        "gene_ids",
        "assignment_ids"
      ),
      collapse = ";"
    ),
    studentization_recomputed_each_assignment = TRUE,
    stringsAsFactors = FALSE
  )
  ptp_write_csv(
    h5_audit,
    file.path(
      dirs$primary$mouse_gene,
      "finite_sample_hdf5_audit.csv"
    )
  )
  out <- list(
    observed = observed,
    freedman_lane = freedman_lane,
    hdf5_path = h5_path,
    index = index,
    assignments = assign_small$table,
    reduced_model = reduced,
    hdf5_audit = h5_audit
  )
  ptpv4_checkpoint_write(
    dirs,
    "primary_gene_statistics",
    out,
    cfg_sha,
    cfg,
    input_sha
  )
  out
}

ptpv4_single_assignment_stats <- function(Y, meta, treated) {
  ploidies <- c("2N", "4N")
  diff <- var <- vector("list", 2L)
  for (k in seq_along(ploidies)) {
    idx <- which(meta$initial_ploidy == ploidies[[k]])
    tr <- idx[treated[idx]]
    co <- idx[!treated[idx]]
    mt <- rowMeans(Y[, tr, drop = FALSE])
    mc <- rowMeans(Y[, co, drop = FALSE])
    vt <- apply(Y[, tr, drop = FALSE], 1L, stats::var)
    vc <- apply(Y[, co, drop = FALSE], 1L, stats::var)
    diff[[k]] <- mt - mc
    var[[k]] <- vt / length(tr) + vc / length(co)
  }
  estimate <- diff[[2L]] - diff[[1L]]
  se <- sqrt(var[[2L]] + var[[1L]])
  list(estimate = estimate, se = se, t = estimate / se)
}

ptpv4_load_primary_data <- function(args, cfg, dirs, seurat_rds) {
  meta <- ptp_read_cell_metadata(args$cell_metadata)
  loaded <- ptp_load_counts_and_qc(
    seurat_rds,
    args$assay,
    args$counts_layer,
    meta
  )
  counts <- loaded$counts
  meta <- loaded$metadata
  ptp_write_csv(
    attr(counts, "seurat_audit"),
    file.path(dirs$qc, "expression_source_audit.csv")
  )
  ptp_write_csv(
    ptp_audit_counts(counts),
    file.path(dirs$qc, "count_matrix_audit.csv")
  )
  matched <- ptp_match_cells(
    meta,
    counts,
    cfg$qc$min_metadata_to_counts_match_rate,
    dirs$qc
  )
  meta <- matched$meta
  counts <- matched$counts
  limited <- ptp_limit_for_smoke(
    meta,
    counts,
    args$max_cells,
    args$max_genes,
    args$seed
  )
  meta <- limited$meta
  counts <- limited$counts
  ptp_write_csv(
    ptp_assignment_audit(meta),
    file.path(dirs$qc, "sample_assignment_audit.csv")
  )
  ptp_write_csv(
    ptp_sample_metadata(meta),
    file.path(dirs$qc, "sample_level_qc.csv")
  )
  expanded <- ptpv2_expand_program_config(cfg)
  cfg2 <- expanded$cfg
  primary_pb <- ptp_choose_primary_grid(meta, counts, cfg2)
  if (!isTRUE(primary_pb$estimable)) {
    stop("Primary matched-state common support is not estimable.", call. = FALSE)
  }
  primary_fit <- ptp_fit_cellmeans(
    primary_pb$counts,
    primary_pb$metadata,
    primary_pb$eligible_subbins,
    "ptp_group",
    cfg2$qc
  )
  observed_symbols <- unique(ptp_clean_gene_symbols(primary_fit$retained_genes))
  programs <- ptp_load_programs(
    cfg2,
    cfg2$qc$min_program_genes_observed,
    observed_symbols
  )
  role_table <- ptpv4_program_role_table(programs)
  ptp_write_csv(
    primary_pb$support,
    file.path(dirs$qc, "primary_common_support.csv")
  )
  ptp_write_csv(
    primary_fit$metadata,
    file.path(dirs$qc, "primary_model_mouse_subbin_metadata.csv")
  )
  list(
    metadata = meta,
    counts = counts,
    config = cfg2,
    primary_pb = primary_pb,
    primary_fit = primary_fit,
    programs = programs,
    role_table = role_table
  )
}

ptpv4_empirical_p_upper <- function(statistic) {
  x <- as.numeric(statistic)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) {
    return(out)
  }
  n <- sum(ok)
  out[ok] <- (1 + rank(-x[ok], ties.method = "max")) / (1 + n)
  out
}

ptpv4_empirical_p_lower <- function(statistic) {
  x <- as.numeric(statistic)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) {
    return(out)
  }
  n <- sum(ok)
  out[ok] <- (1 + rank(x[ok], ties.method = "max")) / (1 + n)
  out
}

ptpv4_true_maxmean <- function(T) {
  positive <- colMeans(pmax(T, 0), na.rm = TRUE)
  negative <- colMeans(pmax(-T, 0), na.rm = TRUE)
  pmax(positive, negative)
}

ptpv4_higher_criticism <- function(T) {
  p <- 2 * stats::pnorm(-abs(T))
  thresholds <- c(0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50)
  n <- nrow(T)
  out <- matrix(NA_real_, nrow = length(thresholds), ncol = ncol(T))
  for (i in seq_along(thresholds)) {
    t0 <- thresholds[[i]]
    count <- colSums(p <= t0, na.rm = TRUE)
    out[i, ] <- (count - n * t0) / sqrt(pmax(n * t0 * (1 - t0), 1e-8))
  }
  apply(out, 2L, max, na.rm = TRUE)
}

ptpv4_covariance_quadratic <- function(T, mouse_expression, lambda = 0.20) {
  Z <- t(scale(t(mouse_expression)))
  Z[!is.finite(Z)] <- 0
  g <- nrow(Z)
  if (g < 2L) {
    return(rep(NA_real_, ncol(T)))
  }
  rank_max <- min(ncol(Z) - 1L, g)
  sv <- tryCatch(
    svd(Z / sqrt(max(ncol(Z) - 1L, 1L)), nu = rank_max, nv = 0L),
    error = function(e) NULL
  )
  if (is.null(sv) || length(sv$d) == 0L) {
    return(colSums(T^2, na.rm = TRUE))
  }
  eig <- head(sv$d^2, ncol(sv$u))
  keep <- is.finite(eig) & eig > 1e-10
  if (!any(keep)) {
    return(colSums(T^2, na.rm = TRUE) / lambda)
  }
  U <- sv$u[, keep, drop = FALSE]
  eig <- eig[keep]
  projection <- crossprod(U, T)
  adjustment <- 1 / (lambda + (1 - lambda) * eig) - 1 / lambda
  colSums(T^2, na.rm = TRUE) / lambda +
    colSums(sweep(projection^2, 1L, adjustment, "*"), na.rm = TRUE)
}

ptpv4_membership_component_statistics <- function(
  membership,
  T,
  mouse_expression,
  cfg
) {
  stats <- list()
  if (!is.na(membership$directional_program_id)) {
    w <- membership$directional_weights$signed_weight
    w[!is.finite(w)] <- 0
    if (sum(abs(w)) > 0L) {
      stats$legal_directional_signed_mean <-
        as.numeric(crossprod(w / sum(abs(w)), T))
    }
  }
  stats$true_maxmean <- ptpv4_true_maxmean(T)
  stats$mean_square <- colMeans(T^2, na.rm = TRUE)
  stats$max_absolute <- apply(abs(T), 2L, max, na.rm = TRUE)
  stats$aspu_gamma_1 <- abs(colMeans(T, na.rm = TRUE))
  stats$aspu_gamma_2 <- sqrt(colMeans(T^2, na.rm = TRUE))
  stats$aspu_gamma_4 <- colMeans(T^4, na.rm = TRUE)^(1 / 4)
  stats$aspu_gamma_8 <- colMeans(T^8, na.rm = TRUE)^(1 / 8)
  stats$aspu_gamma_infinity <- stats$max_absolute
  stats$higher_criticism <- ptpv4_higher_criticism(T)
  stats$shrinkage_covariance_quadratic <- ptpv4_covariance_quadratic(
    T,
    mouse_expression,
    lambda = as.numeric(cfg$adaptive_gene_set$covariance_ridge_lambda)
  )
  stats
}

ptpv4_nested_adaptive_one <- function(
  membership,
  T,
  mouse_expression,
  observed_index,
  cfg
) {
  components <- ptpv4_membership_component_statistics(
    membership,
    T,
    mouse_expression,
    cfg
  )
  p_by_assignment <- list()
  component_rows <- list()
  for (nm in names(components)) {
    s <- components[[nm]]
    s_test <- if (identical(nm, "legal_directional_signed_mean")) abs(s) else s
    p_all <- ptpv4_empirical_p_upper(s_test)
    p_by_assignment[[nm]] <- p_all
    component_rows[[length(component_rows) + 1L]] <- data.frame(
      membership_id = membership$id,
      component = nm,
      observed_statistic = s[[observed_index]],
      empirical_p = p_all[[observed_index]],
      n_assignments = sum(is.finite(s)),
      stringsAsFactors = FALSE
    )
  }
  p_matrix <- do.call(rbind, p_by_assignment)
  minp <- apply(p_matrix, 2L, min, na.rm = TRUE)
  minp[!is.finite(minp)] <- NA_real_
  adaptive_p_by_assignment <- ptpv4_empirical_p_lower(minp)
  list(
    membership_id = membership$id,
    component_rows = ptpv4_bind_rows(component_rows),
    component_statistics = components,
    component_p_by_assignment = p_by_assignment,
    minp_by_assignment = minp,
    adaptive_p_by_assignment = adaptive_p_by_assignment,
    observed_min_component_p = minp[[observed_index]],
    adaptive_empirical_p = adaptive_p_by_assignment[[observed_index]]
  )
}

ptpv4_read_program_gene_matrices <- function(
  gene_stats,
  mouse_bundle,
  membership
) {
  gene_symbols <- ptp_clean_gene_symbols(rownames(mouse_bundle$expression))
  union_genes <- sort(unique(unlist(
    lapply(membership$memberships, `[[`, "genes"),
    use.names = FALSE
  )))
  idx <- match(union_genes, gene_symbols)
  keep <- is.finite(idx)
  union_genes <- union_genes[keep]
  idx <- idx[keep]
  T <- rhdf5::h5read(
    gene_stats$hdf5_path,
    "t_statistic",
    index = list(idx, NULL)
  )
  if (is.vector(T)) {
    T <- matrix(T, nrow = length(idx))
  }
  rownames(T) <- union_genes
  mouse <- mouse_bundle$expression[idx, , drop = FALSE]
  rownames(mouse) <- union_genes
  list(
    t = T,
    mouse_expression = mouse,
    row_index = idx,
    genes = union_genes
  )
}

ptpv4_apply_membership_fdr <- function(
  results,
  membership,
  p_column = "adaptive_empirical_p"
) {
  if (!(p_column %in% names(results))) {
    stop(
      "Membership FDR input lacks P-value column: ",
      p_column,
      call. = FALSE
    )
  }
  alias <- membership$alias
  meta <- ptpv4_bind_rows(lapply(results$membership_id, function(mid) {
    a <- alias[alias$membership_id == mid, , drop = FALSE]
    data.frame(
      membership_id = mid,
      focused = any(a$tier %in% c("tier1_core", "tier2_upgraded_focused")),
      tier3 = any(a$tier == "tier3_exploratory_only"),
      control = any(a$tier == "negative_control" | a$control_status %in% TRUE),
      tiers = paste(sort(unique(a$tier)), collapse = ";"),
      response_families = paste(
        sort(unique(a$response_family)),
        collapse = ";"
      ),
      stringsAsFactors = FALSE
    )
  }))
  out <- merge(results, meta, by = "membership_id", all.x = TRUE, sort = FALSE)
  p_values <- out[[p_column]]
  out$fdr_all_memberships <- stats::p.adjust(
    p_values,
    "BH"
  )
  out$fdr_focused_tier1_tier2 <- NA_real_
  out$fdr_tier3 <- NA_real_
  out$fdr_controls <- NA_real_
  idx <- which(out$focused %in% TRUE & is.finite(p_values))
  if (length(idx)) {
    out$fdr_focused_tier1_tier2[idx] <- stats::p.adjust(
      p_values[idx],
      "BH"
    )
  }
  idx <- which(out$tier3 %in% TRUE & is.finite(p_values))
  if (length(idx)) {
    out$fdr_tier3[idx] <- stats::p.adjust(
      p_values[idx],
      "BH"
    )
  }
  idx <- which(out$control %in% TRUE & is.finite(p_values))
  if (length(idx)) {
    out$fdr_controls[idx] <- stats::p.adjust(
      p_values[idx],
      "BH"
    )
  }
  out
}

ptpv4_run_adaptive_memberships <- function(
  gene_stats,
  mouse_bundle,
  membership,
  cfg,
  args,
  dirs,
  cfg_sha,
  input_sha,
  output_dir = NULL,
  checkpoint_id = "primary_adaptive_memberships"
) {
  if (is.null(output_dir)) {
    output_dir <- dirs$primary$adaptive
  }
  checkpoint <- ptpv4_checkpoint_read(
    dirs,
    checkpoint_id,
    cfg_sha,
    cfg,
    input_sha,
    required_fields = c(
      "membership_results",
      "component_results",
      "minp_matrix",
      "gene_matrices"
    )
  )
  if (!is.null(checkpoint)) {
    return(checkpoint)
  }
  gene_matrices <- ptpv4_read_program_gene_matrices(
    gene_stats,
    mouse_bundle,
    membership
  )
  observed_index <- which(gene_stats$assignments$is_observed_assignment)
  if (!length(observed_index)) {
    observed_index <- nrow(gene_stats$assignments)
  }
  mids <- names(membership$memberships)
  if (isTRUE(args$smoke)) {
    mids <- head(mids, min(8L, length(mids)))
  }
  worker <- function(mid) {
    m <- membership$memberships[[mid]]
    idx <- match(m$genes, rownames(gene_matrices$t))
    idx <- idx[is.finite(idx)]
    if (length(idx) < 2L) {
      return(list(
        status = "not_estimable",
        reason = "fewer_than_two_genes_in_mouse_expression",
        membership_id = mid
      ))
    }
    m$genes <- rownames(gene_matrices$t)[idx]
    if (!is.null(m$directional_weights)) {
      m$directional_weights <- m$directional_weights[
        match(m$genes, m$directional_weights$gene),
        ,
        drop = FALSE
      ]
    }
    out <- ptpv4_nested_adaptive_one(
      m,
      gene_matrices$t[idx, , drop = FALSE],
      gene_matrices$mouse_expression[idx, , drop = FALSE],
      observed_index,
      cfg
    )
    out$status <- "ok"
    out
  }
  workers <- ptpv2_workers(args, length(mids), "adaptive_workers")
  parallel_out <- ptp_parallel_lapply(
    mids,
    worker,
    workers = workers,
    task_label = checkpoint_id
  )
  raw <- parallel_out$results
  ok <- raw[vapply(raw, function(x) identical(x$status, "ok"), logical(1L))]
  failed <- raw[!vapply(raw, function(x) identical(x$status, "ok"), logical(1L))]
  component_results <- ptpv4_bind_rows(lapply(ok, `[[`, "component_rows"))
  membership_results <- ptpv4_bind_rows(lapply(ok, function(x) {
    data.frame(
      membership_id = x$membership_id,
      method_id = "nested_exact_adaptive_gene_set",
      observed_min_component_p = x$observed_min_component_p,
      adaptive_empirical_p = x$adaptive_empirical_p,
      n_assignments = nrow(gene_stats$assignments),
      status = "ok",
      stringsAsFactors = FALSE
    )
  }))
  if (length(failed)) {
    membership_results <- ptpv4_bind_rows(c(
      list(membership_results),
      lapply(failed, function(x) {
        data.frame(
          membership_id = x$membership_id,
          method_id = "nested_exact_adaptive_gene_set",
          observed_min_component_p = NA_real_,
          adaptive_empirical_p = NA_real_,
          n_assignments = nrow(gene_stats$assignments),
          status = x$status,
          reason = x$reason,
          stringsAsFactors = FALSE
        )
      })
    ))
  }
  minp_matrix <- do.call(cbind, lapply(ok, `[[`, "minp_by_assignment"))
  colnames(minp_matrix) <- vapply(ok, `[[`, character(1L), "membership_id")
  adaptive_assignment_p <- do.call(cbind, lapply(
    ok,
    `[[`,
    "adaptive_p_by_assignment"
  ))
  colnames(adaptive_assignment_p) <- colnames(minp_matrix)
  maxT <- apply(
    -log10(pmax(minp_matrix, .Machine$double.xmin)),
    1L,
    max,
    na.rm = TRUE
  )
  membership_results$westfall_young_maxT_p <- vapply(
    membership_results$membership_id,
    function(mid) {
      if (!(mid %in% colnames(minp_matrix))) {
        return(NA_real_)
      }
      score <- -log10(pmax(
        minp_matrix[observed_index, mid],
        .Machine$double.xmin
      ))
      (1 + sum(maxT >= score, na.rm = TRUE)) /
        (1 + sum(is.finite(maxT)))
    },
    numeric(1L)
  )
  membership_results <- ptpv4_apply_membership_fdr(
    membership_results,
    membership
  )
  alias_results <- merge(
    membership$alias,
    membership_results,
    by = "membership_id",
    all.x = TRUE,
    sort = FALSE
  )
  ptp_write_csv(
    component_results,
    file.path(output_dir, "adaptive_component_empirical_pvalues.csv")
  )
  ptp_write_csv(
    membership_results,
    file.path(output_dir, "unique_membership_adaptive_results.csv")
  )
  ptp_write_csv(
    alias_results,
    file.path(output_dir, "program_label_adaptive_alias_results.csv")
  )
  saveRDS(
    list(
      minp_by_assignment = minp_matrix,
      adaptive_p_by_assignment = adaptive_assignment_p,
      assignment_ids = gene_stats$assignments$assignment_id
    ),
    file.path(output_dir, "adaptive_nested_permutation_distributions.rds"),
    compress = "xz"
  )
  out <- list(
    membership_results = membership_results,
    alias_results = alias_results,
    component_results = component_results,
    minp_matrix = minp_matrix,
    adaptive_assignment_p = adaptive_assignment_p,
    raw = raw,
    gene_matrices = gene_matrices,
    observed_index = observed_index,
    parallel_log = parallel_out$log
  )
  ptpv4_checkpoint_write(
    dirs,
    checkpoint_id,
    out,
    cfg_sha,
    cfg,
    input_sha
  )
  out
}

ptpv4_family_definitions <- function(membership) {
  split_alias <- split(
    membership$alias,
    membership$alias$response_family
  )
  lapply(names(split_alias), function(fid) {
    a <- split_alias[[fid]]
    mids <- sort(unique(a$membership_id))
    genes <- sort(unique(unlist(
      lapply(membership$memberships[mids], `[[`, "genes"),
      use.names = FALSE
    )))
    list(
      id = fid,
      label = a$response_family_label[[1L]],
      tiers = sort(unique(a$tier)),
      membership_ids = mids,
      genes = genes,
      union_membership_id = ptpv4_membership_sha(genes)
    )
  })
}

ptpv4_run_family_hierarchy <- function(
  adaptive,
  membership,
  mouse_bundle,
  cfg,
  args,
  dirs,
  output_dir = NULL
) {
  if (is.null(output_dir)) {
    output_dir <- dirs$primary$family
  }
  families <- ptpv4_family_definitions(membership)
  if (isTRUE(args$smoke)) {
    families <- head(families, min(4L, length(families)))
  }
  T_all <- adaptive$gene_matrices$t
  mouse_all <- adaptive$gene_matrices$mouse_expression
  obs <- adaptive$observed_index
  family_rows <- list()
  union_component_rows <- list()
  route_distributions <- list()
  for (f in families) {
    idx <- match(f$genes, rownames(T_all))
    idx <- idx[is.finite(idx)]
    pseudo <- list(
      id = f$union_membership_id,
      genes = rownames(T_all)[idx],
      directional_program_id = NA_character_,
      directional_weights = NULL
    )
    union <- ptpv4_nested_adaptive_one(
      pseudo,
      T_all[idx, , drop = FALSE],
      mouse_all[idx, , drop = FALSE],
      obs,
      cfg
    )
    mids <- intersect(f$membership_ids, colnames(adaptive$adaptive_assignment_p))
    if (length(mids)) {
      member_min <- apply(
        adaptive$adaptive_assignment_p[, mids, drop = FALSE],
        1L,
        min,
        na.rm = TRUE
      )
      member_route_p <- ptpv4_empirical_p_lower(member_min)
    } else {
      member_min <- rep(NA_real_, ncol(T_all))
      member_route_p <- rep(NA_real_, ncol(T_all))
    }
    union_route_p <- union$adaptive_p_by_assignment
    route_min <- pmin(union_route_p, member_route_p, na.rm = TRUE)
    route_min[!is.finite(route_min)] <- NA_real_
    omnibus_by_assignment <- ptpv4_empirical_p_lower(route_min)
    family_rows[[length(family_rows) + 1L]] <- data.frame(
      response_family = f$id,
      response_family_label = f$label,
      tiers = paste(f$tiers, collapse = ";"),
      n_unique_memberships = length(f$membership_ids),
      n_union_genes = length(f$genes),
      union_membership_id = f$union_membership_id,
      union_adaptive_p = union$adaptive_empirical_p,
      member_minP = member_route_p[[obs]],
      family_omnibus_p = omnibus_by_assignment[[obs]],
      stringsAsFactors = FALSE
    )
    u <- union$component_rows
    u$response_family <- f$id
    union_component_rows[[length(union_component_rows) + 1L]] <- u
    route_distributions[[f$id]] <- data.frame(
      assignment_id = seq_along(route_min),
      union_route_p = union_route_p,
      member_route_p = member_route_p,
      selected_route_min_p = route_min,
      family_omnibus_p_by_assignment = omnibus_by_assignment,
      stringsAsFactors = FALSE
    )
  }
  family_results <- ptpv4_bind_rows(family_rows)
  family_results$fdr_all_families <- stats::p.adjust(
    family_results$family_omnibus_p,
    "BH"
  )
  family_results$fdr_focused_11 <- NA_real_
  family_results$fdr_tier3_3 <- NA_real_
  family_results$fdr_controls <- NA_real_
  focused <- grepl("tier1_core|tier2_upgraded_focused", family_results$tiers)
  tier3 <- grepl("tier3_exploratory_only", family_results$tiers)
  controls <- grepl("negative_control", family_results$tiers)
  if (any(focused)) {
    family_results$fdr_focused_11[focused] <- stats::p.adjust(
      family_results$family_omnibus_p[focused],
      "BH"
    )
  }
  if (any(tier3)) {
    family_results$fdr_tier3_3[tier3] <- stats::p.adjust(
      family_results$family_omnibus_p[tier3],
      "BH"
    )
  }
  if (any(controls)) {
    family_results$fdr_controls[controls] <- stats::p.adjust(
      family_results$family_omnibus_p[controls],
      "BH"
    )
  }
  within_rows <- list()
  for (f in families) {
    mids <- intersect(f$membership_ids, colnames(adaptive$minp_matrix))
    if (!length(mids)) {
      next
    }
    scores <- -log10(pmax(
      adaptive$minp_matrix[, mids, drop = FALSE],
      .Machine$double.xmin
    ))
    max_family <- apply(scores, 1L, max, na.rm = TRUE)
    for (mid in mids) {
      observed_score <- scores[obs, mid]
      within_rows[[length(within_rows) + 1L]] <- data.frame(
        response_family = f$id,
        membership_id = mid,
        within_family_westfall_young_p =
          (1 + sum(max_family >= observed_score, na.rm = TRUE)) /
          (1 + sum(is.finite(max_family))),
        stringsAsFactors = FALSE
      )
    }
  }
  within_family <- ptpv4_bind_rows(within_rows)
  alias <- merge(
    membership$alias,
    adaptive$membership_results,
    by = "membership_id",
    all.x = TRUE,
    sort = FALSE
  )
  alias <- merge(
    alias,
    family_results[, c(
      "response_family",
      "family_omnibus_p",
      "fdr_all_families",
      "fdr_focused_11",
      "fdr_tier3_3",
      "fdr_controls"
    )],
    by = "response_family",
    all.x = TRUE,
    sort = FALSE
  )
  alias <- merge(
    alias,
    within_family,
    by = c("response_family", "membership_id"),
    all.x = TRUE,
    sort = FALSE
  )
  alias$focused_claim_gate <- with(
    alias,
    tier %in% c("tier1_core", "tier2_upgraded_focused") &
      is.finite(fdr_focused_11) &
      fdr_focused_11 < 0.10 &
      is.finite(within_family_westfall_young_p) &
      within_family_westfall_young_p < 0.10
  )
  ptp_write_csv(
    family_results,
    file.path(output_dir, "response_family_hierarchical_results.csv")
  )
  ptp_write_csv(
    ptpv4_bind_rows(union_component_rows),
    file.path(output_dir, "family_union_component_results.csv")
  )
  ptp_write_csv(
    within_family,
    file.path(output_dir, "within_family_westfall_young.csv")
  )
  ptp_write_csv(
    alias,
    file.path(output_dir, "program_family_claim_gate.csv")
  )
  saveRDS(
    route_distributions,
    file.path(output_dir, "family_nested_permutation_distributions.rds"),
    compress = "xz"
  )
  list(
    family_results = family_results,
    within_family = within_family,
    alias_claims = alias,
    route_distributions = route_distributions,
    union_components = ptpv4_bind_rows(union_component_rows)
  )
}

ptpv4_membership_scalar_scores <- function(mouse_bundle, membership) {
  genes <- ptp_clean_gene_symbols(rownames(mouse_bundle$expression))
  Z <- t(scale(t(mouse_bundle$expression)))
  Z[!is.finite(Z)] <- 0
  scores <- matrix(
    NA_real_,
    nrow = length(membership$memberships),
    ncol = ncol(Z),
    dimnames = list(names(membership$memberships), colnames(Z))
  )
  for (mid in names(membership$memberships)) {
    m <- membership$memberships[[mid]]
    idx <- match(m$genes, genes)
    keep <- is.finite(idx)
    idx <- idx[keep]
    if (length(idx) < 2L) {
      next
    }
    if (!is.na(m$directional_program_id)) {
      w <- m$directional_weights$signed_weight[keep]
      w[!is.finite(w)] <- 0
      if (sum(abs(w)) > 0L) {
        scores[mid, ] <- as.numeric(crossprod(
          w / sum(abs(w)),
          Z[idx, , drop = FALSE]
        ))
      }
    } else {
      scores[mid, ] <- colMeans(Z[idx, , drop = FALSE], na.rm = TRUE)
    }
  }
  scores
}

ptpv4_run_influence <- function(
  adaptive,
  family,
  mouse_bundle,
  membership,
  cfg,
  dirs
) {
  T <- adaptive$gene_matrices$t
  obs <- adaptive$observed_index
  influence_rows <- list()
  dominant_rows <- list()
  for (mid in names(membership$memberships)) {
    if (!(mid %in% colnames(adaptive$minp_matrix))) {
      next
    }
    m <- membership$memberships[[mid]]
    idx <- match(m$genes, rownames(T))
    idx <- idx[is.finite(idx)]
    z <- T[idx, obs]
    genes <- rownames(T)[idx]
    contribution <- abs(z) / sum(abs(z), na.rm = TRUE)
    ord <- order(contribution, decreasing = TRUE)
    dominant_rows[[length(dominant_rows) + 1L]] <- data.frame(
      membership_id = mid,
      gene_symbol = genes,
      observed_gene_t = z,
      absolute_contribution = contribution,
      contribution_rank = rank(-contribution, ties.method = "first"),
      dominant_top = seq_along(genes) %in% ord[
        seq_len(min(as.integer(cfg$adaptive_gene_set$influence_top_n), length(ord)))
      ],
      stringsAsFactors = FALSE
    )
    full_stat <- ptpv4_true_maxmean(matrix(z, ncol = 1L))[[1L]]
    for (g in seq_along(genes)) {
      z_drop <- z[-g]
      drop_stat <- if (length(z_drop) >= 2L) {
        ptpv4_true_maxmean(matrix(z_drop, ncol = 1L))[[1L]]
      } else {
        NA_real_
      }
      influence_rows[[length(influence_rows) + 1L]] <- data.frame(
        membership_id = mid,
        gene_symbol = genes[[g]],
        n_genes = length(genes),
        set_size_class = if (
          length(genes) <= as.integer(cfg$adaptive_gene_set$small_set_max_genes)
        ) {
          "small"
        } else if (
          length(genes) <= as.integer(cfg$adaptive_gene_set$medium_set_max_genes)
        ) {
          "medium"
        } else {
          "large"
        },
        full_true_maxmean = full_stat,
        leave_one_gene_true_maxmean = drop_stat,
        absolute_change = abs(drop_stat - full_stat),
        stringsAsFactors = FALSE
      )
    }
  }
  dominant <- ptpv4_bind_rows(dominant_rows)
  gene_influence <- ptpv4_bind_rows(influence_rows)
  scalar <- ptpv4_membership_scalar_scores(mouse_bundle, membership)
  meta <- mouse_bundle$metadata
  treated <- meta$treatment != "control"
  full_stats <- ptpv4_single_assignment_stats(scalar, meta, treated)
  loo_rows <- list()
  for (i in seq_len(nrow(meta))) {
    keep <- seq_len(nrow(meta)) != i
    one <- ptpv4_single_assignment_stats(
      scalar[, keep, drop = FALSE],
      meta[keep, , drop = FALSE],
      treated[keep]
    )
    loo_rows[[i]] <- data.frame(
      membership_id = rownames(scalar),
      left_out_mouse = meta$sample_id[[i]],
      full_interaction = full_stats$estimate,
      leave_one_mouse_interaction = one$estimate,
      full_se = full_stats$se,
      leave_one_mouse_se = one$se,
      sign_concordant = sign(one$estimate) == sign(full_stats$estimate),
      standardized_change = abs(one$estimate - full_stats$estimate) /
        pmax(full_stats$se, 1e-8),
      stringsAsFactors = FALSE
    )
  }
  loo <- ptpv4_bind_rows(loo_rows)
  sig_mids <- adaptive$membership_results$membership_id[
    (adaptive$membership_results$fdr_all_memberships < 0.10) |
      (adaptive$membership_results$fdr_focused_tier1_tier2 < 0.10) |
      (adaptive$membership_results$fdr_tier3 < 0.10) |
      (adaptive$membership_results$fdr_controls < 0.10)
  ]
  sig_mids <- unique(sig_mids[!is.na(sig_mids)])
  permutation_loo_gene <- data.frame(
    membership_id = character(0L),
    gene_symbol = character(0L),
    robustness_empirical_p = numeric(0L),
    status = character(0L),
    stringsAsFactors = FALSE
  )
  if (length(sig_mids)) {
    rows <- list()
    for (mid in sig_mids) {
      m <- membership$memberships[[mid]]
      idx <- match(m$genes, rownames(T))
      idx <- idx[is.finite(idx)]
      for (g in seq_along(idx)) {
        keep <- idx[-g]
        stat <- ptpv4_true_maxmean(T[keep, , drop = FALSE])
        p_all <- ptpv4_empirical_p_upper(stat)
        rows[[length(rows) + 1L]] <- data.frame(
          membership_id = mid,
          gene_symbol = rownames(T)[idx[[g]]],
          robustness_empirical_p = p_all[[obs]],
          status = "permutation_calibrated_true_maxmean_leave_one_gene",
          stringsAsFactors = FALSE
        )
      }
    }
    permutation_loo_gene <- ptpv4_bind_rows(rows)
  }
  alias <- membership$alias
  ppp_memberships <- names(membership$memberships)[vapply(
    membership$memberships,
    function(m) "PPP4R2" %in% m$genes,
    logical(1L)
  )]
  hdr_alias <- alias[
    grepl(
      "homolog|repair|fork|hdr",
      paste(alias$program_id, alias$program_label, alias$response_family),
      ignore.case = TRUE
    ),
    ,
    drop = FALSE
  ]
  predeclared <- data.frame(
    audit_id = c("PPP4R2_membership", "HDR_repair_programs"),
    status = c(
      if (length(ppp_memberships)) "present_audited_no_membership_change" else "not_present",
      if (nrow(hdr_alias)) "present_audited_no_membership_change" else "not_present"
    ),
    detail = c(
      paste(ppp_memberships, collapse = ";"),
      paste(unique(hdr_alias$program_id), collapse = ";")
    ),
    stringsAsFactors = FALSE
  )
  ptp_write_csv(
    dominant,
    file.path(dirs$primary$robustness, "dominant_gene_contributions.csv")
  )
  ptp_write_csv(
    gene_influence,
    file.path(dirs$primary$robustness, "leave_one_gene_influence.csv")
  )
  ptp_write_csv(
    loo,
    file.path(dirs$primary$robustness, "leave_one_mouse_interactions.csv")
  )
  ptp_write_csv(
    permutation_loo_gene,
    file.path(
      dirs$primary$robustness,
      "permutation_calibrated_leave_one_gene_robustness.csv"
    )
  )
  ptp_write_csv(
    predeclared,
    file.path(dirs$primary$robustness, "PPP4R2_HDR_predeclared_audit.csv")
  )
  list(
    dominant = dominant,
    gene_influence = gene_influence,
    loo_mouse = loo,
    permutation_loo_gene = permutation_loo_gene,
    predeclared_audit = predeclared
  )
}

ptpv4_run_score_models <- function(primary, cfg, args, dirs) {
  fit <- primary$primary_fit
  programs <- primary$programs
  role <- primary$role_table
  logcpm <- ptp_logcpm_matrix(fit)
  score_defs <- ptp_program_score_methods_from_logcpm(
    logcpm,
    programs,
    cfg,
    pb_counts = fit$dge$counts
  )
  contrasts <- ptp_primary_contrasts(
    primary$primary_pb$eligible_subbins,
    colnames(fit$design)
  )
  result_rows <- list()
  qc_rows <- list()
  score_objects <- list()
  for (method_id in c(
    "signed_weighted_mean_zscore",
    "signed_pca_eigengene",
    "signed_rank_mean",
    "trimmed_signed_mean_zscore"
  )) {
    if (!(method_id %in% names(score_defs))) {
      next
    }
    sm <- score_defs[[method_id]]
    pfit <- ptp_fit_program_scores(sm$score_matrix, fit)
    res <- ptp_apply_program_contrasts(
      pfit,
      contrasts,
      "primary_matched_state",
      programs
    )
    res$score_method_id <- method_id
    res$score_method_label <- sm$method$label %||% method_id
    res <- ptpv3_apply_fdr_family(
      res,
      "p_value",
      "score_method_id",
      role
    )
    result_rows[[method_id]] <- res
    q <- sm$qc_summary
    q$score_method_id <- method_id
    qc_rows[[method_id]] <- q
    score_objects[[method_id]] <- sm$score_matrix
  }
  residual <- ptpv2_within_subbin_residual_scores(fit, programs)
  residual_fit <- ptpv2_fit_program_scores_weighted(
    residual$score_matrix,
    fit,
    residual$observation_weights
  )
  residual_res <- ptp_apply_program_contrasts(
    residual_fit,
    contrasts,
    "primary_matched_state",
    programs
  )
  residual_res$score_method_id <- "within_subbin_residual_zscore"
  residual_res$score_method_label <-
    "Within-subbin residual standardized z-score"
  residual_res <- ptpv3_apply_fdr_family(
    residual_res,
    "p_value",
    "score_method_id",
    role
  )
  result_rows$within_subbin_residual_zscore <- residual_res
  residual$qc_summary$score_method_id <- "within_subbin_residual_zscore"
  qc_rows$within_subbin_residual_zscore <- residual$qc_summary
  score_objects$within_subbin_residual_zscore <- residual$score_matrix
  results <- ptpv4_bind_rows(result_rows)
  qc <- ptpv4_bind_rows(qc_rows)
  precision <- ptpv4_bind_rows(lapply(result_rows, ptp_precision_table))
  score_long <- ptpv4_bind_rows(lapply(names(score_objects), function(mid) {
    x <- score_objects[[mid]]
    data.frame(
      score_method_id = mid,
      program_id = rep(rownames(x), times = ncol(x)),
      sample_subbin_id = rep(colnames(x), each = nrow(x)),
      program_score = as.numeric(x),
      stringsAsFactors = FALSE
    )
  }))
  ptp_write_csv(
    results,
    file.path(dirs$primary$score_models, "all_five_score_model_results.csv")
  )
  ptp_write_csv(
    qc,
    file.path(dirs$primary$score_models, "all_five_score_model_qc.csv")
  )
  ptp_write_csv(
    precision,
    file.path(dirs$primary$score_models, "score_model_precision_and_mde.csv")
  )
  ptp_write_csv(
    score_long,
    file.path(dirs$primary$score_models, "mouse_subbin_program_scores.csv")
  )
  ptp_write_csv(
    residual$gene_scale,
    file.path(dirs$primary$score_models, "within_subbin_gene_scale.csv")
  )
  list(
    results = results,
    qc = qc,
    precision = precision,
    score_matrices = score_objects,
    residual = residual
  )
}

ptpv4_select_mixed_p <- function(df, method, role_table) {
  out <- ptpv3_join_program_meta(df, role_table)
  directional <- out$directional_eligibility %in% TRUE
  if (identical(method, "mroast_msq")) {
    selected <- out$PValue.Mixed
    source <- rep("PValue.Mixed", nrow(out))
  } else {
    selected <- ifelse(directional, out$PValue, out$PValue.Mixed)
    source <- ifelse(directional, "PValue", "PValue.Mixed")
  }
  out$selected_p_value <- selected
  out$selected_p_source <- source
  out$directional_p_legal <- directional
  out$method_id <- method
  ptpv3_apply_fdr_family(
    out,
    "selected_p_value",
    "method_id",
    role_table
  )
}

ptpv4_run_fry_mroast <- function(
  primary_fit,
  programs,
  role_table,
  cfg,
  args,
  dirs,
  output_dir = NULL,
  coordinate_mode = FALSE
) {
  if (is.null(output_dir)) {
    output_dir <- dirs$primary$fry_mroast
  }
  cfg_local <- cfg
  cfg_local$v2_methods$self_contained$mroast_rotations_observed <- if (
    isTRUE(args$smoke)
  ) {
    19L
  } else if (isTRUE(coordinate_mode)) {
    as.integer(cfg$self_contained$coordinate_mroast_rotations)
  } else {
    as.integer(cfg$self_contained$primary_mroast_rotations)
  }
  cfg_local$v2_methods$self_contained$mroast_seed <-
    as.integer(cfg$reproducibility$rotation_seed)
  contrasts <- ptp_primary_contrasts(
    unique(primary_fit$metadata$subbin_id),
    colnames(primary_fit$design)
  )
  raw <- ptpv2_self_contained_tests(
    primary_fit,
    contrasts$treatment_by_initial_ploidy_interaction,
    programs,
    role_table,
    cfg_local,
    args
  )
  fry <- ptpv4_select_mixed_p(raw$fry, "fry", role_table)
  mroast_mean <- ptpv4_select_mixed_p(
    raw$mroast_mean,
    "mroast_mean",
    role_table
  )
  mroast_msq <- ptpv4_select_mixed_p(
    raw$mroast_msq,
    "mroast_msq",
    role_table
  )
  ptp_write_csv(
    fry,
    file.path(output_dir, "fry_mixed_p_corrected.csv")
  )
  ptp_write_csv(
    mroast_mean,
    file.path(output_dir, "mroast_mean_mixed_p_corrected.csv")
  )
  ptp_write_csv(
    mroast_msq,
    file.path(output_dir, "mroast_msq_mixed_p_corrected.csv")
  )
  audit <- ptpv4_bind_rows(lapply(
    list(fry = fry, mroast_mean = mroast_mean, mroast_msq = mroast_msq),
    function(x) {
      aggregate(
        program_id ~ selected_p_source + directional_p_legal,
        data = x,
        FUN = length
      )
    }
  ))
  ptp_write_csv(
    audit,
    file.path(output_dir, "mixed_p_multiplicity_audit.csv")
  )
  list(fry = fry, mroast_mean = mroast_mean, mroast_msq = mroast_msq)
}

ptpv4_run_camera <- function(primary, cfg, dirs) {
  gene <- ptp_apply_contrasts(
    primary$primary_fit,
    ptp_primary_contrasts(
      primary$primary_pb$eligible_subbins,
      colnames(primary$primary_fit$design)
    ),
    "primary_matched_state"
  )
  grid <- ptpv2_camera_grid(
    gene,
    primary$programs,
    ptpv3_branch_role_table(primary$role_table),
    cfg
  )
  grid <- ptpv3_apply_fdr_family(
    grid,
    "PValue",
    "method_id",
    primary$role_table
  )
  residual_corr <- ptpv2_residual_corr(
    primary$primary_fit,
    primary$programs,
    ptpv3_branch_role_table(primary$role_table)
  )
  x <- gene[
    gene$contrast_id == "treatment_by_initial_ploidy_interaction",
    ,
    drop = FALSE
  ]
  statistics <- x$t_statistic
  names(statistics) <- x$gene_symbol
  empirical_rows <- list()
  for (id in names(primary$programs)) {
    idx <- which(names(statistics) %in% primary$programs[[id]]$observed_genes)
    rho <- residual_corr$empirical_inter_gene_correlation[
      match(id, residual_corr$program_id)
    ]
    rho <- pmin(pmax(rho, -0.20), 0.99)
    if (length(idx) < 2L || !is.finite(rho)) {
      empirical_rows[[length(empirical_rows) + 1L]] <- data.frame(
        program_id = id,
        method_id = "cameraPR_residual_empirical_rho",
        PValue = NA_real_,
        Direction = NA_character_,
        empirical_rho = rho,
        status = "not_estimable",
        stringsAsFactors = FALSE
      )
      next
    }
    cam <- tryCatch(
      limma::cameraPR(
        statistic = statistics,
        index = list(program = idx),
        inter.gene.cor = rho,
        use.ranks = FALSE,
        sort = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(cam)) {
      empirical_rows[[length(empirical_rows) + 1L]] <- data.frame(
        program_id = id,
        method_id = "cameraPR_residual_empirical_rho",
        PValue = NA_real_,
        Direction = NA_character_,
        empirical_rho = rho,
        status = "failed",
        stringsAsFactors = FALSE
      )
    } else {
      cam <- as.data.frame(cam, stringsAsFactors = FALSE)
      cam$program_id <- id
      cam$method_id <- "cameraPR_residual_empirical_rho"
      cam$empirical_rho <- rho
      cam$status <- "ok"
      empirical_rows[[length(empirical_rows) + 1L]] <- cam
    }
  }
  empirical <- ptpv3_apply_fdr_family(
    ptpv4_bind_rows(empirical_rows),
    "PValue",
    "method_id",
    primary$role_table
  )
  all <- ptpv4_bind_rows(list(grid, empirical))
  controls <- all[
    all$tier == "negative_control" & is.finite(all$PValue),
    ,
    drop = FALSE
  ]
  negative_fail <- any(controls$fdr_controls < 0.10, na.rm = TRUE) ||
    sum(controls$PValue < 0.05, na.rm = TRUE) > 1L
  all$interpretation_status <- if (negative_fail) {
    "diagnostic_only"
  } else {
    "calibrated_diagnostic"
  }
  calibration <- data.frame(
    method_family = "camera",
    n_control_tests = nrow(controls),
    n_control_nominal_0_05 = sum(controls$PValue < 0.05, na.rm = TRUE),
    n_control_fdr_0_10 = sum(controls$fdr_controls < 0.10, na.rm = TRUE),
    status = if (negative_fail) "diagnostic_only" else "calibrated_diagnostic",
    rule = cfg$camera$negative_control_fail_rule,
    stringsAsFactors = FALSE
  )
  ptp_write_csv(
    grid,
    file.path(dirs$primary$camera, "camera_original_rho_grid.csv")
  )
  ptp_write_csv(
    residual_corr,
    file.path(dirs$primary$camera, "residual_empirical_correlation.csv")
  )
  ptp_write_csv(
    empirical,
    file.path(dirs$primary$camera, "camera_residual_empirical_rho.csv")
  )
  ptp_write_csv(
    calibration,
    file.path(dirs$primary$camera, "camera_negative_control_calibration.csv")
  )
  list(
    results = all,
    grid = grid,
    empirical = empirical,
    residual_correlation = residual_corr,
    calibration = calibration
  )
}

ptpv4_raw_program_scores <- function(logcpm, programs) {
  row_symbols <- ptp_clean_gene_symbols(rownames(logcpm))
  out <- matrix(
    NA_real_,
    nrow = length(programs),
    ncol = ncol(logcpm),
    dimnames = list(names(programs), colnames(logcpm))
  )
  for (id in names(programs)) {
    p <- programs[[id]]
    idx <- match(p$observed_genes, row_symbols)
    keep <- is.finite(idx)
    idx <- idx[keep]
    w <- p$observed_gene_weights$signed_weight[
      match(row_symbols[idx], p$observed_gene_weights$gene)
    ]
    ok <- is.finite(w) & w != 0
    if (sum(ok) < 2L) {
      next
    }
    out[id, ] <- as.numeric(crossprod(
      w[ok] / sum(abs(w[ok])),
      logcpm[idx[ok], , drop = FALSE]
    ))
  }
  out
}

ptpv4_control_reference_assignment <- function(
  raw_scores,
  meta,
  operator,
  treated_samples,
  prior = 4
) {
  designated_control <- !(meta$sample_id %in% treated_samples)
  global_sd <- apply(raw_scores, 1L, stats::sd, na.rm = TRUE)
  R <- raw_scores * NA_real_
  reference_n <- integer(ncol(raw_scores))
  reference_ids <- character(ncol(raw_scores))
  self_inclusion <- logical(ncol(raw_scores))
  for (j in seq_len(ncol(raw_scores))) {
    candidates <- which(
      meta$subbin_id == meta$subbin_id[[j]] &
        meta$initial_ploidy == meta$initial_ploidy[[j]] &
        designated_control
    )
    if (designated_control[[j]]) {
      candidates <- candidates[meta$sample_id[candidates] != meta$sample_id[[j]]]
    }
    reference_mice <- unique(meta$sample_id[candidates])
    reference_n[[j]] <- length(reference_mice)
    reference_ids[[j]] <- paste(sort(reference_mice), collapse = ";")
    self_inclusion[[j]] <- meta$sample_id[[j]] %in% reference_mice
    if (length(reference_mice) < 2L || length(candidates) < 2L) {
      next
    }
    mu <- rowMeans(raw_scores[, candidates, drop = FALSE], na.rm = TRUE)
    sd0 <- apply(raw_scores[, candidates, drop = FALSE], 1L, stats::sd)
    shrunk <- sqrt(
      ((length(candidates) - 1L) * sd0^2 + prior * global_sd^2) /
        ((length(candidates) - 1L) + prior)
    )
    shrunk[!is.finite(shrunk) | shrunk <= 0] <-
      global_sd[!is.finite(shrunk) | shrunk <= 0]
    R[, j] <- (raw_scores[, j] - mu) / pmax(shrunk, 1e-8)
  }
  mouse_scores <- R %*% t(operator)
  list(
    observation_scores = R,
    mouse_scores = mouse_scores,
    reference_n = reference_n,
    reference_ids = reference_ids,
    self_inclusion = self_inclusion
  )
}

ptpv4_run_control_reference <- function(
  primary,
  mouse_bundle,
  assignments,
  cfg,
  args,
  dirs
) {
  raw_scores <- ptpv4_raw_program_scores(
    ptp_logcpm_matrix(primary$primary_fit),
    primary$programs
  )
  meta <- primary$primary_fit$metadata
  B <- if (isTRUE(args$smoke)) {
    nrow(assignments$table[
      seq_len(min(args$smoke_permutations, nrow(assignments$table))),
      ,
      drop = FALSE
    ])
  } else {
    nrow(assignments$table)
  }
  keep <- seq_len(B)
  if (!(assignments$observed_index %in% keep)) {
    keep[[B]] <- assignments$observed_index
  }
  table <- assignments$table[keep, , drop = FALSE]
  t_matrix <- matrix(
    NA_real_,
    nrow = nrow(raw_scores),
    ncol = nrow(table),
    dimnames = list(rownames(raw_scores), table$assignment_id)
  )
  estimate_matrix <- t_matrix
  se_matrix <- t_matrix
  observed_audit <- NULL
  for (b in seq_len(nrow(table))) {
    treated <- strsplit(table$treated_samples[[b]], ";", fixed = TRUE)[[1L]]
    ref <- ptpv4_control_reference_assignment(
      raw_scores,
      meta,
      mouse_bundle$operator,
      treated,
      prior = as.numeric(cfg$control_reference$sd_shrinkage_prior_weight)
    )
    label <- mouse_bundle$metadata$sample_id %in% treated
    stat <- ptpv4_single_assignment_stats(
      ref$mouse_scores,
      mouse_bundle$metadata,
      label
    )
    t_matrix[, b] <- stat$t
    estimate_matrix[, b] <- stat$estimate
    se_matrix[, b] <- stat$se
    if (isTRUE(table$is_observed_assignment[[b]])) {
      observed_audit <- data.frame(
        sample_subbin_id = meta$sample_subbin_id,
        sample_id = meta$sample_id,
        subbin_id = meta$subbin_id,
        initial_ploidy = meta$initial_ploidy,
        observed_designated_control =
          !(meta$sample_id %in% treated),
        reference_mouse_count = ref$reference_n,
        reference_mouse_ids = ref$reference_ids,
        control_self_inclusion = ref$self_inclusion,
        stringsAsFactors = FALSE
      )
    }
    if (b %% 500L == 0L) {
      ptpv4_message(
        "Control-reference assignments ",
        b,
        "/",
        nrow(table),
        log_file = file.path(dirs$logs, "run.log")
      )
    }
  }
  obs <- which(table$is_observed_assignment)
  if (!length(obs)) {
    obs <- nrow(table)
  }
  p <- vapply(seq_len(nrow(t_matrix)), function(i) {
    (1 + sum(abs(t_matrix[i, ]) >= abs(t_matrix[i, obs]), na.rm = TRUE)) /
      (1 + sum(is.finite(t_matrix[i, ])))
  }, numeric(1L))
  crit <- stats::qt(0.975, df = 12)
  results <- data.frame(
    program_id = rownames(t_matrix),
    method_id = "control_reference_residual_score",
    estimate = estimate_matrix[, obs],
    se = se_matrix[, obs],
    ci_low = estimate_matrix[, obs] - crit * se_matrix[, obs],
    ci_high = estimate_matrix[, obs] + crit * se_matrix[, obs],
    p_value = p,
    n_assignments = nrow(table),
    reference_recomputed_each_assignment = TRUE,
    control_self_inclusion_allowed = FALSE,
    stringsAsFactors = FALSE
  )
  results <- ptpv3_apply_fdr_family(
    results,
    "p_value",
    "method_id",
    primary$role_table
  )
  ptp_write_csv(
    results,
    file.path(
      dirs$primary$score_models,
      "control_reference_residual_score.csv"
    )
  )
  ptp_write_csv(
    observed_audit,
    file.path(
      dirs$primary$score_models,
      "control_reference_observed_training_audit.csv"
    )
  )
  ptp_write_csv(
    data.frame(
      n_assignments = nrow(table),
      reference_recomputed_each_assignment = TRUE,
      reference_scale = "program_score_level",
      treated_reference = "assignment_controls",
      control_reference = "leave_one_control_mouse_out",
      any_observed_self_inclusion = any(
        observed_audit$control_self_inclusion,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    ),
    file.path(
      dirs$primary$score_models,
      "control_reference_permutation_audit.csv"
    )
  )
  list(
    results = results,
    observed_audit = observed_audit,
    t_matrix = t_matrix,
    assignments = table
  )
}

ptpv4_run_covariance_whitened <- function(
  primary,
  mouse_bundle,
  assignments,
  cfg,
  dirs
) {
  row_symbols <- ptp_clean_gene_symbols(rownames(mouse_bundle$expression))
  meta <- mouse_bundle$metadata
  controls <- which(meta$treatment == "control")
  result_rows <- list()
  audit_rows <- list()
  dominant_rows <- list()
  for (id in primary$role_table$program_id[
    primary$role_table$directional_eligibility %in% TRUE
  ]) {
    p <- primary$programs[[id]]
    idx <- match(p$observed_genes, row_symbols)
    keep <- is.finite(idx)
    idx <- idx[keep]
    genes <- row_symbols[idx]
    w <- p$observed_gene_weights$signed_weight[
      match(genes, p$observed_gene_weights$gene)
    ]
    valid <- is.finite(w) & w != 0
    idx <- idx[valid]
    genes <- genes[valid]
    w <- w[valid]
    if (length(idx) > as.integer(cfg$covariance_whitened$max_genes)) {
      ord <- order(abs(w), decreasing = TRUE)[
        seq_len(as.integer(cfg$covariance_whitened$max_genes))
      ]
      idx <- idx[ord]
      genes <- genes[ord]
      w <- w[ord]
    }
    X <- t(scale(t(mouse_bundle$expression[idx, , drop = FALSE])))
    X[!is.finite(X)] <- 0
    S <- stats::cov(t(X[, controls, drop = FALSE]))
    S[!is.finite(S)] <- 0
    lambda <- as.numeric(cfg$covariance_whitened$ridge_lambda)
    ridge <- lambda * stats::median(diag(S)[diag(S) > 0], na.rm = TRUE)
    if (!is.finite(ridge) || ridge <= 0) {
      ridge <- lambda
    }
    Sr <- S + diag(ridge, nrow(S))
    eig <- eigen(Sr, symmetric = TRUE, only.values = TRUE)$values
    condition <- max(eig) / pmax(min(eig), 1e-12)
    fallback <- ""
    inverse <- tryCatch(solve(Sr), error = function(e) NULL)
    if (
      is.null(inverse) ||
      !is.finite(condition) ||
      condition > as.numeric(cfg$covariance_whitened$condition_number_fallback)
    ) {
      inverse <- diag(1 / pmax(diag(Sr), 1e-8), nrow(Sr))
      fallback <- "diagonal_inverse"
    }
    a <- matrix(w, ncol = 1L)
    direction <- inverse %*% a
    denom <- sqrt(as.numeric(crossprod(a, direction)))
    score <- as.numeric(crossprod(direction / denom, X))
    stat <- ptpv4_interaction_stats(
      matrix(score, nrow = 1L, dimnames = list(id, colnames(X))),
      assignments
    )
    obs <- assignments$observed_index
    pvalue <- (1 + sum(
      abs(stat$t[1L, ]) >= abs(stat$t[1L, obs]),
      na.rm = TRUE
    )) / (1 + sum(is.finite(stat$t[1L, ])))
    crit <- stats::qt(0.975, df = 12)
    result_rows[[length(result_rows) + 1L]] <- data.frame(
      program_id = id,
      method_id = "covariance_whitened_score",
      estimate = stat$estimate[1L, obs],
      se = stat$se[1L, obs],
      ci_low = stat$estimate[1L, obs] - crit * stat$se[1L, obs],
      ci_high = stat$estimate[1L, obs] + crit * stat$se[1L, obs],
      p_value = pvalue,
      condition_number = condition,
      effective_rank = sum(eig)^2 / sum(eig^2),
      ridge_value = ridge,
      fallback = fallback,
      n_genes = length(idx),
      stringsAsFactors = FALSE
    )
    contribution <- abs(as.numeric(direction)) /
      sum(abs(as.numeric(direction)))
    dominant_rows[[length(dominant_rows) + 1L]] <- data.frame(
      program_id = id,
      gene_symbol = genes,
      whitened_direction_weight = as.numeric(direction),
      absolute_contribution = contribution,
      contribution_rank = rank(-contribution, ties.method = "first"),
      stringsAsFactors = FALSE
    )
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      program_id = id,
      status = "ok",
      condition_number = condition,
      effective_rank = sum(eig)^2 / sum(eig^2),
      fallback = fallback,
      metadata_complete = TRUE,
      stringsAsFactors = FALSE
    )
  }
  results <- ptpv3_apply_fdr_family(
    ptpv4_bind_rows(result_rows),
    "p_value",
    "method_id",
    primary$role_table
  )
  audit <- ptpv4_bind_rows(audit_rows)
  dominant <- ptpv4_bind_rows(dominant_rows)
  ptp_write_csv(
    results,
    file.path(
      dirs$primary$score_models,
      "covariance_whitened_directional_results.csv"
    )
  )
  ptp_write_csv(
    audit,
    file.path(
      dirs$primary$score_models,
      "covariance_whitened_metadata_audit.csv"
    )
  )
  ptp_write_csv(
    dominant,
    file.path(
      dirs$primary$score_models,
      "covariance_whitened_dominant_genes.csv"
    )
  )
  list(results = results, audit = audit, dominant = dominant)
}

ptpv4_run_multidimensional <- function(
  mouse_bundle,
  membership,
  assignments,
  cfg,
  args,
  dirs
) {
  symbols <- ptp_clean_gene_symbols(rownames(mouse_bundle$expression))
  result_rows <- list()
  loading_rows <- list()
  audit_rows <- list()
  mids <- names(membership$memberships)
  if (isTRUE(args$smoke)) {
    mids <- head(mids, min(8L, length(mids)))
  }
  for (mid in mids) {
    m <- membership$memberships[[mid]]
    idx <- match(m$genes, symbols)
    idx <- idx[is.finite(idx)]
    if (length(idx) < 3L) {
      result_rows[[length(result_rows) + 1L]] <- data.frame(
        membership_id = mid,
        method_id = "treatment_blind_multidimensional_pca",
        status = "not_estimable",
        reason = "fewer_than_three_genes",
        p_value = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    X <- t(scale(t(mouse_bundle$expression[idx, , drop = FALSE])))
    X[!is.finite(X)] <- 0
    pc <- tryCatch(
      stats::prcomp(t(X), center = TRUE, scale. = FALSE),
      error = function(e) NULL
    )
    if (is.null(pc) || ncol(pc$x) == 0L) {
      result_rows[[length(result_rows) + 1L]] <- data.frame(
        membership_id = mid,
        method_id = "treatment_blind_multidimensional_pca",
        status = "not_estimable",
        reason = "PCA_failed",
        p_value = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    variance <- pc$sdev^2 / sum(pc$sdev^2)
    k <- which(cumsum(variance) >=
      as.numeric(cfg$multidimensional$variance_explained_target))[[1L]]
    k <- min(k, as.integer(cfg$multidimensional$max_pcs), ncol(pc$x))
    score <- t(pc$x[, seq_len(k), drop = FALSE])
    stat <- ptpv4_interaction_stats(score, assignments)
    joint <- colSums(stat$t^2, na.rm = TRUE)
    obs <- assignments$observed_index
    pvalue <- (1 + sum(joint >= joint[[obs]], na.rm = TRUE)) /
      (1 + sum(is.finite(joint)))
    basis <- pc$rotation[, seq_len(k), drop = FALSE]
    basis_sha <- ptpv4_object_sha(list(
      genes = symbols[idx],
      loadings = basis,
      variance = variance[seq_len(k)]
    ))
    result_rows[[length(result_rows) + 1L]] <- data.frame(
      membership_id = mid,
      method_id = "treatment_blind_multidimensional_pca",
      status = "ok",
      n_genes = length(idx),
      n_pcs = k,
      variance_explained = sum(variance[seq_len(k)]),
      observed_joint_statistic = joint[[obs]],
      p_value = pvalue,
      n_assignments = nrow(assignments$table),
      basis_checksum = basis_sha,
      stringsAsFactors = FALSE
    )
    load <- data.frame(
      membership_id = mid,
      gene_symbol = rep(symbols[idx], times = k),
      pc = rep(seq_len(k), each = length(idx)),
      loading = as.numeric(basis),
      basis_checksum = basis_sha,
      stringsAsFactors = FALSE
    )
    loading_rows[[length(loading_rows) + 1L]] <- load
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      membership_id = mid,
      basis_training = cfg$multidimensional$basis_training,
      treatment_labels_used_for_basis = FALSE,
      basis_retrained_per_assignment = FALSE,
      basis_invariant_across_assignments = TRUE,
      basis_checksum = basis_sha,
      n_assignment_checks = nrow(assignments$table),
      status = "ok",
      stringsAsFactors = FALSE
    )
  }
  results <- ptpv4_apply_membership_fdr(
    ptpv4_bind_rows(result_rows),
    membership,
    p_column = "p_value"
  )
  alias <- merge(
    membership$alias,
    results,
    by = "membership_id",
    all.x = TRUE,
    sort = FALSE
  )
  loadings <- ptpv4_bind_rows(loading_rows)
  audit <- ptpv4_bind_rows(audit_rows)
  ptp_write_csv(
    results,
    file.path(
      dirs$primary$adaptive,
      "multidimensional_unique_membership_results.csv"
    )
  )
  ptp_write_csv(
    alias,
    file.path(
      dirs$primary$adaptive,
      "multidimensional_program_alias_results.csv"
    )
  )
  ptp_write_csv(
    loadings,
    file.path(dirs$primary$adaptive, "multidimensional_pc_loadings.csv")
  )
  ptp_write_csv(
    audit,
    file.path(dirs$primary$adaptive, "multidimensional_basis_audit.csv")
  )
  list(results = results, alias = alias, loadings = loadings, audit = audit)
}

ptpv4_assignment_columns <- function(assignments, args) {
  n <- nrow(assignments$table)
  if (!isTRUE(args$smoke) || n <= as.integer(args$smoke_permutations)) {
    return(seq_len(n))
  }
  keep <- seq_len(min(as.integer(args$smoke_permutations), n))
  if (!(assignments$observed_index %in% keep)) {
    keep[[length(keep)]] <- assignments$observed_index
  }
  unique(keep)
}

ptpv4_stratum_stats_missing <- function(Y, A) {
  valid <- is.finite(Y)
  Y0 <- Y
  Y0[!valid] <- 0
  V <- valid * 1
  n_t <- V %*% A
  n_c <- V %*% (!A)
  sum_t <- Y0 %*% A
  sum_c <- Y0 %*% (!A)
  sq_t <- (Y0^2) %*% A
  sq_c <- (Y0^2) %*% (!A)
  mean_t <- sum_t / n_t
  mean_c <- sum_c / n_c
  ss_t <- sq_t - n_t * mean_t^2
  ss_c <- sq_c - n_c * mean_c^2
  var_t <- ss_t / pmax(n_t - 1, 1)
  var_c <- ss_c / pmax(n_c - 1, 1)
  var_contrast <- var_t / n_t + var_c / n_c
  invalid <- n_t < 2 | n_c < 2
  difference <- mean_t - mean_c
  difference[invalid] <- NA_real_
  var_contrast[invalid] <- NA_real_
  list(
    difference = difference,
    variance = pmax(var_contrast, 0),
    n_treated = n_t,
    n_control = n_c
  )
}

ptpv4_interaction_stats_missing <- function(Y, mouse_meta, assignments) {
  s1 <- assignments$strata[[1L]]
  s2 <- assignments$strata[[2L]]
  Y1 <- Y[, match(s1, mouse_meta$sample_id), drop = FALSE]
  Y2 <- Y[, match(s2, mouse_meta$sample_id), drop = FALSE]
  z1 <- ptpv4_stratum_stats_missing(Y1, assignments$combinations[[1L]])
  z2 <- ptpv4_stratum_stats_missing(Y2, assignments$combinations[[2L]])
  n1 <- ncol(z1$difference)
  n2 <- ncol(z2$difference)
  B <- n1 * n2
  estimate <- matrix(NA_real_, nrow = nrow(Y), ncol = B)
  se <- matrix(NA_real_, nrow = nrow(Y), ncol = B)
  k <- 0L
  for (i in seq_len(n1)) {
    cols <- k + seq_len(n2)
    estimate[, cols] <- sweep(
      z2$difference,
      1L,
      z1$difference[, i],
      "-"
    )
    se[, cols] <- sqrt(sweep(
      z2$variance,
      1L,
      z1$variance[, i],
      "+"
    ))
    k <- max(cols)
  }
  tstat <- estimate / se
  tstat[!is.finite(tstat)] <- NA_real_
  rownames(estimate) <- rownames(Y)
  rownames(se) <- rownames(Y)
  rownames(tstat) <- rownames(Y)
  colnames(estimate) <- assignments$table$assignment_id
  colnames(se) <- assignments$table$assignment_id
  colnames(tstat) <- assignments$table$assignment_id
  list(estimate = estimate, se = se, t = tstat)
}

ptpv4_grid_regions <- function(n_bins, widths = NULL, prefix = "bin") {
  if (is.null(widths)) {
    widths <- 1L
  }
  regions <- list()
  rows <- list()
  for (w in as.integer(widths)) {
    for (i in seq_len(n_bins - w + 1L)) {
      j <- i + w - 1L
      id <- if (length(widths) == 1L && w == 1L) {
        sprintf("%s%02d", prefix, i)
      } else {
        sprintf("w%02d_b%02d_%02d", w, i, j)
      }
      start <- (i - 1L) / n_bins
      end <- j / n_bins
      regions[[id]] <- list(
        id = id,
        label = sprintf("[%.2f, %.2f%s", start, end, if (j == n_bins) "]" else ")"),
        start = start,
        end = end,
        include_start = TRUE,
        include_end = j == n_bins
      )
      rows[[length(rows) + 1L]] <- data.frame(
        region_id = id,
        start_bin = i,
        end_bin = j,
        width_bins = w,
        start = start,
        end = end,
        midpoint = (start + end) / 2,
        stringsAsFactors = FALSE
      )
    }
  }
  list(regions = regions, table = ptpv4_bind_rows(rows))
}

ptpv4_region_support <- function(pb_meta, min_mice = 2L) {
  keys <- expand.grid(
    region_id = unique(pb_meta$region_id),
    initial_ploidy = c("2N", "4N"),
    treatment = c("control", "gemcitabine"),
    stringsAsFactors = FALSE
  )
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, , drop = FALSE]
    x <- pb_meta[
      pb_meta$region_id == k$region_id &
        pb_meta$initial_ploidy == k$initial_ploidy &
        pb_meta$treatment == k$treatment,
      ,
      drop = FALSE
    ]
    data.frame(
      region_id = k$region_id,
      initial_ploidy = k$initial_ploidy,
      treatment = k$treatment,
      n_mice_total = length(unique(x$sample_id)),
      n_mice_qc = length(unique(x$sample_id[x$retained_for_model])),
      total_cells = sum(x$cell_count, na.rm = TRUE),
      minimum_cells = if (nrow(x)) min(x$cell_count) else NA_real_,
      support_pass = length(unique(x$sample_id[x$retained_for_model])) >=
        min_mice,
      stringsAsFactors = FALSE
    )
  })
  out <- ptpv4_bind_rows(rows)
  eligible <- aggregate(support_pass ~ region_id, out, all)
  names(eligible)[2L] <- "eligible_support_only"
  merge(out, eligible, by = "region_id", all.x = TRUE, sort = FALSE)
}

ptpv4_region_expression <- function(pb, region_id, sample_ids, cfg) {
  meta <- pb$metadata[
    pb$metadata$region_id == region_id & pb$metadata$retained_for_model,
    ,
    drop = FALSE
  ]
  if (nrow(meta) < 8L) {
    return(NULL)
  }
  counts <- pb$counts[, meta$sample_region_id, drop = FALSE]
  keep_gene <- Matrix::rowSums(counts) > 0
  if (sum(keep_gene) < 10L) {
    return(NULL)
  }
  dge <- edgeR::DGEList(counts = counts[keep_gene, , drop = FALSE])
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  E0 <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  E <- matrix(
    NA_real_,
    nrow = nrow(E0),
    ncol = length(sample_ids),
    dimnames = list(rownames(E0), sample_ids)
  )
  E[, match(meta$sample_id, sample_ids)] <- E0
  list(expression = E, metadata = meta, dge = dge)
}

ptpv4_entity_definitions <- function(membership) {
  entities <- lapply(membership$memberships, function(m) {
    list(
      id = m$id,
      entity_type = "unique_membership",
      genes = m$genes,
      directional_program_id = m$directional_program_id,
      directional_weights = m$directional_weights
    )
  })
  families <- ptpv4_family_definitions(membership)
  for (f in families) {
    id <- paste0("family__", f$id)
    entities[[id]] <- list(
      id = id,
      entity_type = "response_family_union",
      response_family = f$id,
      genes = f$genes,
      directional_program_id = NA_character_,
      directional_weights = NULL
    )
  }
  entities
}

ptpv4_entity_scores <- function(expression, entities) {
  symbols <- ptp_clean_gene_symbols(rownames(expression))
  Z <- t(scale(t(expression)))
  Z[!is.finite(Z)] <- NA_real_
  out <- matrix(
    NA_real_,
    nrow = length(entities),
    ncol = ncol(expression),
    dimnames = list(names(entities), colnames(expression))
  )
  for (id in names(entities)) {
    e <- entities[[id]]
    idx <- match(e$genes, symbols)
    keep <- is.finite(idx)
    idx <- idx[keep]
    if (length(idx) < 2L) {
      next
    }
    if (!is.na(e$directional_program_id)) {
      w <- e$directional_weights$signed_weight[keep]
      ok <- is.finite(w) & w != 0
      if (sum(ok) >= 2L) {
        out[id, ] <- as.numeric(crossprod(
          w[ok] / sum(abs(w[ok])),
          Z[idx[ok], , drop = FALSE]
        ))
      }
    } else {
      out[id, ] <- colMeans(Z[idx, , drop = FALSE], na.rm = TRUE)
    }
  }
  out
}

ptpv4_entity_metadata <- function(entities, membership) {
  ptpv4_bind_rows(lapply(entities, function(e) {
    if (identical(e$entity_type, "unique_membership")) {
      a <- membership$alias[
        membership$alias$membership_id == e$id,
        ,
        drop = FALSE
      ]
      data.frame(
        entity_id = e$id,
        entity_type = e$entity_type,
        response_family = paste(sort(unique(a$response_family)), collapse = ";"),
        tiers = paste(sort(unique(a$tier)), collapse = ";"),
        n_genes = length(e$genes),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        entity_id = e$id,
        entity_type = e$entity_type,
        response_family = e$response_family,
        tiers = paste(
          sort(unique(membership$alias$tier[
            membership$alias$response_family == e$response_family
          ])),
          collapse = ";"
        ),
        n_genes = length(e$genes),
        stringsAsFactors = FALSE
      )
    }
  }))
}

ptpv4_run_multiscale <- function(
  meta,
  counts,
  mouse_meta,
  assignments,
  membership,
  cfg,
  args,
  dirs,
  output_dir = NULL,
  prefix = "primary"
) {
  if (is.null(output_dir)) {
    output_dir <- dirs$primary$multiscale
  }
  grid <- ptpv4_grid_regions(
    as.integer(cfg$multiscale$n_bins),
    unlist(cfg$multiscale$window_widths),
    prefix = "bin"
  )
  pb <- ptp_construct_region_pseudobulk(
    meta,
    counts,
    grid$regions,
    cfg$qc,
    ptp_endpoint_ploidy_specs(cfg)
  )
  support <- ptpv4_region_support(
    pb$metadata,
    as.integer(cfg$multiscale$support_min_mice_per_group)
  )
  support <- merge(support, grid$table, by = "region_id", all.x = TRUE)
  ptp_write_csv(
    grid$table,
    file.path(output_dir, "frozen_contiguous_windows.csv")
  )
  ptp_write_csv(
    support,
    file.path(output_dir, "window_support_audit.csv")
  )
  eligible <- unique(support$region_id[support$eligible_support_only])
  entities <- ptpv4_entity_definitions(membership)
  entity_meta <- ptpv4_entity_metadata(entities, membership)
  assignment_cols <- ptpv4_assignment_columns(assignments, args)
  assignment_ids <- assignments$table$assignment_id[assignment_cols]
  obs_full <- assignments$observed_index
  obs <- match(obs_full, assignment_cols)
  if (!is.finite(obs)) {
    stop("Observed assignment was lost from multiscale subset.", call. = FALSE)
  }
  max_abs <- matrix(
    -Inf,
    nrow = length(entities),
    ncol = length(assignment_cols),
    dimnames = list(names(entities), assignment_ids)
  )
  best_window <- rep(NA_character_, length(entities))
  best_observed_abs <- rep(-Inf, length(entities))
  names(best_window) <- names(best_observed_abs) <- names(entities)
  observed_rows <- list()
  for (wi in seq_along(eligible)) {
    region_id <- eligible[[wi]]
    rb <- ptpv4_region_expression(
      pb,
      region_id,
      mouse_meta$sample_id,
      cfg
    )
    if (is.null(rb)) {
      next
    }
    score <- ptpv4_entity_scores(rb$expression, entities)
    stat_full <- ptpv4_interaction_stats_missing(score, mouse_meta, assignments)
    tstat <- stat_full$t[, assignment_cols, drop = FALSE]
    estimate <- stat_full$estimate[, assignment_cols, drop = FALSE]
    se <- stat_full$se[, assignment_cols, drop = FALSE]
    local_abs <- abs(tstat)
    replace <- is.finite(local_abs) &
      (!is.finite(max_abs) | local_abs > max_abs)
    max_abs[replace] <- local_abs[replace]
    obs_abs <- local_abs[, obs]
    better <- is.finite(obs_abs) & obs_abs > best_observed_abs
    best_window[better] <- region_id
    best_observed_abs[better] <- obs_abs[better]
    info <- grid$table[grid$table$region_id == region_id, , drop = FALSE]
    observed_rows[[length(observed_rows) + 1L]] <- data.frame(
      entity_id = rownames(score),
      entity_type = vapply(entities, `[[`, character(1L), "entity_type"),
      region_id = region_id,
      start_bin = info$start_bin,
      end_bin = info$end_bin,
      width_bins = info$width_bins,
      start = info$start,
      end = info$end,
      observed_estimate = estimate[, obs],
      observed_se = se[, obs],
      observed_t = tstat[, obs],
      unadjusted_empirical_p = vapply(seq_len(nrow(tstat)), function(i) {
        (1 + sum(local_abs[i, ] >= local_abs[i, obs], na.rm = TRUE)) /
          (1 + sum(is.finite(local_abs[i, ])))
      }, numeric(1L)),
      inference_role = "location_only_until_scan_max_calibration",
      stringsAsFactors = FALSE
    )
    if (wi %% 10L == 0L) {
      ptpv4_message(
        prefix,
        " multiscale windows ",
        wi,
        "/",
        length(eligible),
        log_file = file.path(dirs$logs, "run.log")
      )
    }
  }
  max_abs[!is.finite(max_abs)] <- NA_real_
  scan_results <- data.frame(
    entity_id = rownames(max_abs),
    best_window = best_window,
    observed_scan_max_abs_t = max_abs[, obs],
    scan_adjusted_empirical_p = vapply(seq_len(nrow(max_abs)), function(i) {
      (1 + sum(max_abs[i, ] >= max_abs[i, obs], na.rm = TRUE)) /
        (1 + sum(is.finite(max_abs[i, ])))
    }, numeric(1L)),
    n_eligible_windows = length(eligible),
    n_assignments = ncol(max_abs),
    best_window_inference_role = "localization_only",
    stringsAsFactors = FALSE
  )
  scan_results <- merge(
    scan_results,
    entity_meta,
    by.x = "entity_id",
    by.y = "entity_id",
    all.x = TRUE,
    sort = FALSE
  )
  scan_results$fdr_all_entities <- stats::p.adjust(
    scan_results$scan_adjusted_empirical_p,
    "BH"
  )
  scan_results$fdr_unique_memberships <- NA_real_
  scan_results$fdr_response_families <- NA_real_
  ii <- which(scan_results$entity_type == "unique_membership")
  ff <- which(scan_results$entity_type == "response_family_union")
  if (length(ii)) {
    scan_results$fdr_unique_memberships[ii] <- stats::p.adjust(
      scan_results$scan_adjusted_empirical_p[ii],
      "BH"
    )
  }
  if (length(ff)) {
    scan_results$fdr_response_families[ff] <- stats::p.adjust(
      scan_results$scan_adjusted_empirical_p[ff],
      "BH"
    )
  }
  observed_windows <- ptpv4_bind_rows(observed_rows)
  permutation_max <- data.frame(
    entity_id = rep(rownames(max_abs), times = ncol(max_abs)),
    assignment_id = rep(colnames(max_abs), each = nrow(max_abs)),
    permutation_scan_max_abs_t = as.numeric(max_abs),
    stringsAsFactors = FALSE
  )
  ptp_write_csv(
    observed_windows,
    file.path(output_dir, "all_window_observed_statistics.csv")
  )
  ptp_write_csv(
    permutation_max,
    file.path(output_dir, "permutation_window_scan_max_statistics.csv")
  )
  ptp_write_csv(
    scan_results,
    file.path(output_dir, "multiscale_scan_adjusted_results.csv")
  )
  saveRDS(
    list(
      max_abs_t = max_abs,
      assignment_ids = assignment_ids,
      observed_index = obs,
      eligible_windows = eligible
    ),
    file.path(output_dir, "multiscale_permutation_distributions.rds"),
    compress = "xz"
  )
  list(
    grid = grid$table,
    support = support,
    eligible_windows = eligible,
    observed_windows = observed_windows,
    permutation_max = max_abs,
    scan_results = scan_results
  )
}

ptpv4_count_support_interval <- function(meta, start, end, include_end, cfg) {
  hit <- meta$pseudotime >= start &
    if (include_end) meta$pseudotime <= end else meta$pseudotime < end
  x <- meta[hit, , drop = FALSE]
  mice <- unique(meta[, c("sample_id", "initial_ploidy", "treatment")])
  tab <- as.data.frame(table(
    factor(x$sample_id, levels = mice$sample_id)
  ), stringsAsFactors = FALSE)
  names(tab) <- c("sample_id", "cell_count")
  tab$cell_count <- as.integer(tab$cell_count)
  tab <- merge(mice, tab, by = "sample_id", all.x = TRUE, sort = FALSE)
  tab$retained <- tab$cell_count >= as.integer(cfg$qc$min_cells_per_mouse_subbin)
  group <- aggregate(
    retained ~ initial_ploidy + treatment,
    tab,
    sum
  )
  grid <- expand.grid(
    initial_ploidy = c("2N", "4N"),
    treatment = c("control", "gemcitabine"),
    stringsAsFactors = FALSE
  )
  group <- merge(grid, group, by = c("initial_ploidy", "treatment"), all.x = TRUE)
  group$retained[is.na(group$retained)] <- 0
  list(mouse = tab, group = group, min_support = min(group$retained))
}

ptpv4_merge_trajectory_bins <- function(meta, cfg) {
  n <- as.integer(cfg$trajectory_global$initial_bins)
  bins <- data.frame(
    original_start_bin = seq_len(n),
    original_end_bin = seq_len(n),
    start = (seq_len(n) - 1) / n,
    end = seq_len(n) / n,
    stringsAsFactors = FALSE
  )
  audit <- list()
  iteration <- 0L
  min_required <- as.integer(cfg$trajectory_global$support_min_mice_per_group)
  repeat {
    support <- lapply(seq_len(nrow(bins)), function(i) {
      ptpv4_count_support_interval(
        meta,
        bins$start[[i]],
        bins$end[[i]],
        bins$original_end_bin[[i]] == n,
        cfg
      )
    })
    min_support <- vapply(support, `[[`, numeric(1L), "min_support")
    if (all(min_support >= min_required) || nrow(bins) == 1L) {
      break
    }
    failing <- which(min_support < min_required)[[1L]]
    candidates <- integer(0L)
    if (failing > 1L) {
      candidates <- c(candidates, failing - 1L)
    }
    if (failing < nrow(bins)) {
      candidates <- c(candidates, failing + 1L)
    }
    candidate_rows <- lapply(candidates, function(j) {
      lo <- min(failing, j)
      hi <- max(failing, j)
      s <- ptpv4_count_support_interval(
        meta,
        bins$start[[lo]],
        bins$end[[hi]],
        bins$original_end_bin[[hi]] == n,
        cfg
      )
      data.frame(
        neighbor = j,
        side = if (j < failing) "left" else "right",
        min_support = s$min_support,
        total_support = sum(s$group$retained),
        stringsAsFactors = FALSE
      )
    })
    candidate_table <- ptpv4_bind_rows(candidate_rows)
    candidate_table$side_priority <- ifelse(candidate_table$side == "left", 0, 1)
    candidate_table <- candidate_table[order(
      -candidate_table$min_support,
      -candidate_table$total_support,
      candidate_table$side_priority
    ), , drop = FALSE]
    chosen <- candidate_table$neighbor[[1L]]
    lo <- min(failing, chosen)
    hi <- max(failing, chosen)
    iteration <- iteration + 1L
    audit[[iteration]] <- data.frame(
      iteration = iteration,
      failing_bin_before = paste0(
        bins$original_start_bin[[failing]],
        "-",
        bins$original_end_bin[[failing]]
      ),
      failing_min_group_support = min_support[[failing]],
      chosen_side = candidate_table$side[[1L]],
      chosen_neighbor_before = paste0(
        bins$original_start_bin[[chosen]],
        "-",
        bins$original_end_bin[[chosen]]
      ),
      resulting_min_group_support = candidate_table$min_support[[1L]],
      rule = cfg$trajectory_global$merge_rule,
      stringsAsFactors = FALSE
    )
    merged <- data.frame(
      original_start_bin = bins$original_start_bin[[lo]],
      original_end_bin = bins$original_end_bin[[hi]],
      start = bins$start[[lo]],
      end = bins$end[[hi]],
      stringsAsFactors = FALSE
    )
    bins <- rbind(
      if (lo > 1L) bins[seq_len(lo - 1L), , drop = FALSE] else NULL,
      merged,
      if (hi < nrow(bins)) bins[seq.int(hi + 1L, nrow(bins)), , drop = FALSE] else NULL
    )
    rownames(bins) <- NULL
  }
  bins$region_id <- sprintf("merged_bin_%02d", seq_len(nrow(bins)))
  bins$midpoint <- (bins$start + bins$end) / 2
  bins$width <- bins$end - bins$start
  bins$include_end <- bins$original_end_bin == n
  list(bins = bins, audit = ptpv4_bind_rows(audit))
}

ptpv4_run_trajectory_global <- function(
  meta,
  counts,
  mouse_meta,
  assignments,
  membership,
  cfg,
  args,
  dirs,
  output_dir = NULL,
  prefix = "primary"
) {
  if (is.null(output_dir)) {
    output_dir <- dirs$primary$trajectory
  }
  merged <- ptpv4_merge_trajectory_bins(meta, cfg)
  regions <- lapply(seq_len(nrow(merged$bins)), function(i) {
    x <- merged$bins[i, , drop = FALSE]
    list(
      id = x$region_id,
      label = sprintf("[%.2f, %.2f%s", x$start, x$end, if (x$include_end) "]" else ")"),
      start = x$start,
      end = x$end,
      include_start = TRUE,
      include_end = x$include_end
    )
  })
  names(regions) <- merged$bins$region_id
  pb <- ptp_construct_region_pseudobulk(
    meta,
    counts,
    regions,
    cfg$qc,
    ptp_endpoint_ploidy_specs(cfg)
  )
  support <- ptpv4_region_support(
    pb$metadata,
    as.integer(cfg$trajectory_global$support_min_mice_per_group)
  )
  support <- merge(
    support,
    merged$bins,
    by = "region_id",
    all.x = TRUE,
    sort = FALSE
  )
  ptp_write_csv(
    merged$bins,
    file.path(output_dir, "support_only_merged_bin_definitions.csv")
  )
  ptp_write_csv(
    merged$audit,
    file.path(output_dir, "support_only_bin_merge_audit.csv")
  )
  ptp_write_csv(
    support,
    file.path(output_dir, "trajectory_mouse_support.csv")
  )
  eligible <- unique(support$region_id[support$eligible_support_only])
  eligible <- merged$bins$region_id[merged$bins$region_id %in% eligible]
  entities <- ptpv4_entity_definitions(membership)
  entity_meta <- ptpv4_entity_metadata(entities, membership)
  assignment_cols <- ptpv4_assignment_columns(assignments, args)
  obs <- match(assignments$observed_index, assignment_cols)
  B <- length(assignment_cols)
  stat_by_bin <- list()
  for (region_id in eligible) {
    rb <- ptpv4_region_expression(
      pb,
      region_id,
      mouse_meta$sample_id,
      cfg
    )
    if (is.null(rb)) {
      next
    }
    score <- ptpv4_entity_scores(rb$expression, entities)
    full <- ptpv4_interaction_stats_missing(score, mouse_meta, assignments)
    stat_by_bin[[region_id]] <- list(
      estimate = full$estimate[, assignment_cols, drop = FALSE],
      se = full$se[, assignment_cols, drop = FALSE],
      t = full$t[, assignment_cols, drop = FALSE]
    )
  }
  eligible <- intersect(eligible, names(stat_by_bin))
  if (length(eligible) < 2L) {
    status <- data.frame(
      method_id = "trajectory_wide_global_interaction",
      status = "not_estimable",
      reason = "fewer_than_two_support_eligible_merged_bins",
      stringsAsFactors = FALSE
    )
    ptp_write_csv(status, file.path(output_dir, "trajectory_global_results.csv"))
    return(list(status = status, support = support, merge = merged))
  }
  mids <- merged$bins$midpoint[match(eligible, merged$bins$region_id)]
  widths <- merged$bins$width[match(eligible, merged$bins$region_id)]
  curve_rows <- list()
  global_rows <- list()
  integrated_rows <- list()
  distribution <- list()
  basis_rows <- list()
  interval <- as.numeric(unlist(cfg$trajectory_global$integration_interval))
  for (entity_id in names(entities)) {
    T <- do.call(rbind, lapply(stat_by_bin, function(x) x$t[entity_id, ]))
    Est <- do.call(rbind, lapply(stat_by_bin, function(x) x$estimate[entity_id, ]))
    Se <- do.call(rbind, lapply(stat_by_bin, function(x) x$se[entity_id, ]))
    rownames(T) <- rownames(Est) <- rownames(Se) <- eligible
    max_abs <- apply(abs(T), 2L, max, na.rm = TRUE)
    crit <- suppressWarnings(stats::quantile(
      max_abs[-obs],
      probs = 0.95,
      na.rm = TRUE,
      names = FALSE
    ))
    curve_rows[[length(curve_rows) + 1L]] <- data.frame(
      entity_id = entity_id,
      entity_type = entities[[entity_id]]$entity_type,
      region_id = eligible,
      midpoint = mids,
      start = merged$bins$start[match(eligible, merged$bins$region_id)],
      end = merged$bins$end[match(eligible, merged$bins$region_id)],
      observed_estimate = Est[, obs],
      observed_se = Se[, obs],
      observed_t = T[, obs],
      simultaneous_95_low = Est[, obs] - crit * Se[, obs],
      simultaneous_95_high = Est[, obs] + crit * Se[, obs],
      simultaneous_maxT_critical = crit,
      supported = TRUE,
      stringsAsFactors = FALSE
    )
    df_stats <- list()
    for (df in as.integer(unlist(cfg$trajectory_global$spline_df_candidates))) {
      if (nrow(T) <= df) {
        next
      }
      basis <- splines::ns(mids, df = df)
      q <- qr.Q(qr(basis))
      projected <- crossprod(q, T)
      s <- colSums(projected^2, na.rm = TRUE)
      df_stats[[paste0("df", df)]] <- s
      basis_rows[[length(basis_rows) + 1L]] <- data.frame(
        entity_id = entity_id,
        spline_df = df,
        basis_checksum = ptpv4_object_sha(basis),
        treatment_labels_used = FALSE,
        selected_inside_permutation = TRUE,
        stringsAsFactors = FALSE
      )
    }
    if (!length(df_stats)) {
      next
    }
    df_matrix <- do.call(rbind, df_stats)
    selected <- apply(df_matrix, 2L, max, na.rm = TRUE)
    selected_df <- rownames(df_matrix)[max.col(
      t(df_matrix),
      ties.method = "first"
    )]
    p <- (1 + sum(selected >= selected[[obs]], na.rm = TRUE)) /
      (1 + sum(is.finite(selected)))
    global_rows[[length(global_rows) + 1L]] <- data.frame(
      entity_id = entity_id,
      entity_type = entities[[entity_id]]$entity_type,
      observed_global_statistic = selected[[obs]],
      observed_selected_spline_df = selected_df[[obs]],
      exact_whole_mouse_permutation_p = p,
      n_assignments = B,
      df_selection_calibrated_inside_permutation = TRUE,
      stringsAsFactors = FALSE
    )
    overlap <- pmax(
      0,
      pmin(
        merged$bins$end[match(eligible, merged$bins$region_id)],
        interval[[2L]]
      ) - pmax(
        merged$bins$start[match(eligible, merged$bins$region_id)],
        interval[[1L]]
      )
    )
    w <- overlap / sum(overlap)
    if (sum(overlap) > 0) {
      int_est <- colSums(Est * w, na.rm = TRUE)
      int_se <- sqrt(colSums((Se * w)^2, na.rm = TRUE))
      int_t <- int_est / int_se
      int_p <- (1 + sum(
        abs(int_t) >= abs(int_t[[obs]]),
        na.rm = TRUE
      )) / (1 + sum(is.finite(int_t)))
      integrated_rows[[length(integrated_rows) + 1L]] <- data.frame(
        entity_id = entity_id,
        entity_type = entities[[entity_id]]$entity_type,
        interval_start = interval[[1L]],
        interval_end = interval[[2L]],
        observed_integrated_interaction = int_est[[obs]],
        observed_integrated_se = int_se[[obs]],
        observed_integrated_t = int_t[[obs]],
        exact_empirical_p = int_p,
        reference_distribution = "same_ploidy_stratified_4900_assignments",
        stringsAsFactors = FALSE
      )
    }
    distribution[[entity_id]] <- list(
      spline_df_statistics = df_stats,
      selected_global_statistic = selected,
      selected_df = selected_df,
      max_abs_bin_t = max_abs
    )
  }
  global <- merge(
    ptpv4_bind_rows(global_rows),
    entity_meta,
    by = "entity_id",
    all.x = TRUE,
    sort = FALSE
  )
  global$fdr_all_entities <- stats::p.adjust(
    global$exact_whole_mouse_permutation_p,
    "BH"
  )
  ii <- global$entity_type.x %in% "unique_membership" |
    global$entity_type.y %in% "unique_membership"
  ff <- global$entity_type.x %in% "response_family_union" |
    global$entity_type.y %in% "response_family_union"
  global$fdr_unique_memberships <- NA_real_
  global$fdr_response_families <- NA_real_
  if (any(ii)) {
    global$fdr_unique_memberships[ii] <- stats::p.adjust(
      global$exact_whole_mouse_permutation_p[ii],
      "BH"
    )
  }
  if (any(ff)) {
    global$fdr_response_families[ff] <- stats::p.adjust(
      global$exact_whole_mouse_permutation_p[ff],
      "BH"
    )
  }
  curve <- ptpv4_bind_rows(curve_rows)
  integrated <- ptpv4_bind_rows(integrated_rows)
  ptp_write_csv(
    curve,
    file.path(output_dir, "trajectory_interaction_curve_simultaneous_band.csv")
  )
  ptp_write_csv(
    global,
    file.path(output_dir, "trajectory_global_results.csv")
  )
  ptp_write_csv(
    integrated,
    file.path(output_dir, "standardized_integrated_interaction_030_049.csv")
  )
  ptp_write_csv(
    ptpv4_bind_rows(basis_rows),
    file.path(output_dir, "trajectory_spline_basis_audit.csv")
  )
  saveRDS(
    list(
      distributions = distribution,
      assignment_ids = assignments$table$assignment_id[assignment_cols],
      observed_index = obs,
      eligible_bins = eligible
    ),
    file.path(output_dir, "trajectory_global_permutation_distributions.rds"),
    compress = "xz"
  )
  list(
    support = support,
    merge = merged,
    curve = curve,
    global = global,
    integrated = integrated,
    distributions = distribution
  )
}

ptpv4_run_adaptive_from_matrices <- function(
  T,
  mouse_expression,
  assignments,
  membership,
  cfg,
  args,
  output_dir,
  prefix = "matrix"
) {
  obs <- assignments$observed_index
  mids <- names(membership$memberships)
  if (isTRUE(args$smoke)) {
    mids <- head(mids, min(8L, length(mids)))
  }
  raw <- lapply(mids, function(mid) {
    m <- membership$memberships[[mid]]
    idx <- match(m$genes, rownames(T))
    idx <- idx[is.finite(idx)]
    if (length(idx) < 2L) {
      return(list(
        status = "not_estimable",
        reason = "fewer_than_two_genes",
        membership_id = mid
      ))
    }
    m$genes <- rownames(T)[idx]
    if (!is.null(m$directional_weights)) {
      m$directional_weights <- m$directional_weights[
        match(m$genes, m$directional_weights$gene),
        ,
        drop = FALSE
      ]
    }
    x <- ptpv4_nested_adaptive_one(
      m,
      T[idx, , drop = FALSE],
      mouse_expression[idx, , drop = FALSE],
      obs,
      cfg
    )
    x$status <- "ok"
    x
  })
  ok <- raw[vapply(raw, function(x) identical(x$status, "ok"), logical(1L))]
  failed <- raw[!vapply(raw, function(x) identical(x$status, "ok"), logical(1L))]
  component <- ptpv4_bind_rows(lapply(ok, `[[`, "component_rows"))
  results <- ptpv4_bind_rows(lapply(ok, function(x) {
    data.frame(
      membership_id = x$membership_id,
      method_id = paste0(prefix, "_nested_exact_adaptive_gene_set"),
      observed_min_component_p = x$observed_min_component_p,
      adaptive_empirical_p = x$adaptive_empirical_p,
      n_assignments = nrow(assignments$table),
      status = "ok",
      stringsAsFactors = FALSE
    )
  }))
  if (length(failed)) {
    results <- ptpv4_bind_rows(c(
      list(results),
      lapply(failed, function(x) {
        data.frame(
          membership_id = x$membership_id,
          method_id = paste0(prefix, "_nested_exact_adaptive_gene_set"),
          observed_min_component_p = NA_real_,
          adaptive_empirical_p = NA_real_,
          n_assignments = nrow(assignments$table),
          status = x$status,
          reason = x$reason,
          stringsAsFactors = FALSE
        )
      })
    ))
  }
  minp <- do.call(cbind, lapply(ok, `[[`, "minp_by_assignment"))
  adaptive_p <- do.call(cbind, lapply(ok, `[[`, "adaptive_p_by_assignment"))
  colnames(minp) <- colnames(adaptive_p) <- vapply(
    ok,
    `[[`,
    character(1L),
    "membership_id"
  )
  maxT <- apply(
    -log10(pmax(minp, .Machine$double.xmin)),
    1L,
    max,
    na.rm = TRUE
  )
  results$westfall_young_maxT_p <- vapply(results$membership_id, function(mid) {
    if (!(mid %in% colnames(minp))) {
      return(NA_real_)
    }
    score <- -log10(pmax(minp[obs, mid], .Machine$double.xmin))
    (1 + sum(maxT >= score, na.rm = TRUE)) /
      (1 + sum(is.finite(maxT)))
  }, numeric(1L))
  results <- ptpv4_apply_membership_fdr(results, membership)
  alias <- merge(
    membership$alias,
    results,
    by = "membership_id",
    all.x = TRUE,
    sort = FALSE
  )
  ptp_write_csv(
    component,
    file.path(output_dir, paste0(prefix, "_adaptive_component_pvalues.csv"))
  )
  ptp_write_csv(
    results,
    file.path(output_dir, paste0(prefix, "_unique_membership_adaptive.csv"))
  )
  ptp_write_csv(
    alias,
    file.path(output_dir, paste0(prefix, "_program_alias_adaptive.csv"))
  )
  saveRDS(
    list(
      minp_by_assignment = minp,
      adaptive_p_by_assignment = adaptive_p,
      assignment_ids = assignments$table$assignment_id
    ),
    file.path(output_dir, paste0(prefix, "_adaptive_distributions.rds")),
    compress = "xz"
  )
  list(
    membership_results = results,
    alias_results = alias,
    component_results = component,
    minp_matrix = minp,
    adaptive_assignment_p = adaptive_p,
    raw = raw,
    gene_matrices = list(t = T, mouse_expression = mouse_expression),
    observed_index = obs
  )
}

ptpv4_aggregate_mouse_counts <- function(meta, counts) {
  sample_meta <- ptp_sample_metadata(meta)
  ids <- sample_meta$sample_id
  group <- factor(meta$sample_id, levels = ids)
  design <- Matrix::sparse.model.matrix(~ 0 + group)
  colnames(design) <- ids
  pb <- counts[, meta$cell_id, drop = FALSE] %*% design
  rownames(pb) <- rownames(counts)
  sample_meta <- sample_meta[match(ids, sample_meta$sample_id), , drop = FALSE]
  list(counts = pb, metadata = sample_meta)
}

ptpv4_fit_independent_mouse_model <- function(pb, cfg) {
  meta <- pb$metadata
  counts <- pb$counts
  y0 <- edgeR::DGEList(counts = counts)
  keep <- rowSums(edgeR::cpm(y0) > cfg$qc$gene_filter_cpm) >=
    cfg$qc$gene_filter_min_observations
  y <- edgeR::calcNormFactors(
    edgeR::DGEList(counts = counts[keep, , drop = FALSE]),
    method = "TMM"
  )
  meta$initial_ploidy <- factor(meta$initial_ploidy, levels = c("2N", "4N"))
  meta$treatment <- factor(meta$treatment, levels = c("control", "gemcitabine"))
  design <- stats::model.matrix(~ initial_ploidy * treatment, data = meta)
  v <- limma::voomWithQualityWeights(y, design, plot = FALSE)
  fit <- limma::eBayes(limma::lmFit(v, design), robust = TRUE)
  list(
    fit = fit,
    voom = v,
    dge = y,
    design = design,
    metadata = meta,
    retained_genes = rownames(y$counts),
    group_col = "independent_mouse"
  )
}

ptpv4_independent_contrasts <- function(design_columns) {
  make <- function(values) {
    out <- rep(0, length(design_columns))
    names(out) <- design_columns
    out[names(values)] <- values
    out
  }
  list(
    delta_2N_origin = make(c(treatmentgemcitabine = 1)),
    delta_4N_origin = make(c(
      treatmentgemcitabine = 1,
      `initial_ploidy4N:treatmentgemcitabine` = 1
    )),
    treatment_by_initial_ploidy_interaction = make(c(
      `initial_ploidy4N:treatmentgemcitabine` = 1
    ))
  )
}

ptpv4_fit_independent_scores <- function(score_matrix, model) {
  keep <- rowSums(is.finite(score_matrix)) == ncol(score_matrix)
  score_matrix <- score_matrix[keep, , drop = FALSE]
  fit <- limma::eBayes(
    limma::lmFit(score_matrix, model$design),
    robust = FALSE
  )
  list(
    fit = fit,
    design = model$design,
    metadata = model$metadata,
    score_matrix = score_matrix
  )
}

ptpv4_marginal_residual_scores <- function(model, programs) {
  E <- ptp_logcpm_matrix(model)
  meta <- model$metadata
  reduced <- stats::model.matrix(~ initial_ploidy + treatment, data = meta)
  fit <- limma::eBayes(limma::lmFit(E, reduced), robust = TRUE)
  residual <- E - fit$coefficients %*% t(reduced)
  scale0 <- sqrt(fit$s2.post)
  Z <- sweep(residual, 1L, pmax(scale0, 1e-8), "/")
  symbols <- ptp_clean_gene_symbols(rownames(Z))
  score <- matrix(
    NA_real_,
    nrow = length(programs),
    ncol = ncol(Z),
    dimnames = list(names(programs), colnames(Z))
  )
  for (id in names(programs)) {
    p <- programs[[id]]
    idx <- match(p$observed_genes, symbols)
    keep <- is.finite(idx)
    idx <- idx[keep]
    w <- p$observed_gene_weights$signed_weight[keep]
    ok <- is.finite(w) & w != 0
    if (sum(ok) >= 2L) {
      score[id, ] <- as.numeric(crossprod(
        w[ok] / sum(abs(w[ok])),
        Z[idx[ok], , drop = FALSE]
      ))
    }
  }
  list(score_matrix = score, gene_scale = scale0, reduced_design = reduced)
}

ptpv4_run_marginal_fry_mroast <- function(
  model,
  programs,
  role_table,
  cfg,
  args,
  output_dir
) {
  symbols <- ptp_clean_gene_symbols(rownames(model$voom$E))
  index <- lapply(programs, function(p) {
    which(symbols %in% p$observed_genes)
  })
  index <- index[lengths(index) >= 2L]
  contrast <- ptpv4_independent_contrasts(
    colnames(model$design)
  )[["treatment_by_initial_ploidy_interaction"]]
  nrot <- if (isTRUE(args$smoke)) {
    19L
  } else {
    as.integer(cfg$self_contained$primary_mroast_rotations)
  }
  fry_raw <- as.data.frame(limma::fry(
    model$voom,
    index = index,
    design = model$design,
    contrast = contrast,
    sort = "none"
  ), stringsAsFactors = FALSE)
  fry_raw$program_id <- rownames(fry_raw)
  set.seed(as.integer(cfg$reproducibility$rotation_seed))
  mean_raw <- as.data.frame(limma::mroast(
    model$voom,
    index = index,
    design = model$design,
    contrast = contrast,
    nrot = nrot,
    set.statistic = "mean",
    sort = "none"
  ), stringsAsFactors = FALSE)
  mean_raw$program_id <- rownames(mean_raw)
  set.seed(as.integer(cfg$reproducibility$rotation_seed) + 1L)
  msq_raw <- as.data.frame(limma::mroast(
    model$voom,
    index = index,
    design = model$design,
    contrast = contrast,
    nrot = nrot,
    set.statistic = "msq",
    sort = "none"
  ), stringsAsFactors = FALSE)
  msq_raw$program_id <- rownames(msq_raw)
  fry <- ptpv4_select_mixed_p(fry_raw, "marginal_fry", role_table)
  mean <- ptpv4_select_mixed_p(mean_raw, "marginal_mroast_mean", role_table)
  msq <- ptpv4_select_mixed_p(msq_raw, "mroast_msq", role_table)
  ptp_write_csv(fry, file.path(output_dir, "marginal_fry_mixed_p.csv"))
  ptp_write_csv(mean, file.path(output_dir, "marginal_mroast_mean_mixed_p.csv"))
  ptp_write_csv(msq, file.path(output_dir, "marginal_mroast_msq_mixed_p.csv"))
  list(fry = fry, mroast_mean = mean, mroast_msq = msq)
}

ptpv4_run_marginal <- function(
  primary,
  membership,
  assignments,
  cfg,
  args,
  dirs
) {
  pb <- ptpv4_aggregate_mouse_counts(primary$metadata, primary$counts)
  model <- ptpv4_fit_independent_mouse_model(pb, cfg)
  logcpm <- ptp_logcpm_matrix(model)
  defs <- ptp_program_score_methods_from_logcpm(
    logcpm,
    primary$programs,
    cfg,
    pb_counts = model$dge$counts
  )
  contrasts <- ptpv4_independent_contrasts(colnames(model$design))
  rows <- list()
  score_matrices <- list()
  for (method_id in c(
    "signed_weighted_mean_zscore",
    "signed_pca_eigengene",
    "signed_rank_mean",
    "trimmed_signed_mean_zscore"
  )) {
    if (!(method_id %in% names(defs))) {
      next
    }
    pfit <- ptpv4_fit_independent_scores(defs[[method_id]]$score_matrix, model)
    x <- ptp_apply_program_contrasts(
      pfit,
      contrasts,
      "marginal_total_effect",
      primary$programs
    )
    x$score_method_id <- method_id
    x <- ptpv3_apply_fdr_family(
      x,
      "p_value",
      "score_method_id",
      primary$role_table
    )
    rows[[method_id]] <- x
    score_matrices[[method_id]] <- defs[[method_id]]$score_matrix
  }
  residual <- ptpv4_marginal_residual_scores(model, primary$programs)
  pfit <- ptpv4_fit_independent_scores(residual$score_matrix, model)
  r <- ptp_apply_program_contrasts(
    pfit,
    contrasts,
    "marginal_total_effect",
    primary$programs
  )
  r$score_method_id <- "within_subbin_residual_zscore"
  r$score_method_definition <-
    "marginal reduced-model residual analogue; no subbin conditioning"
  r <- ptpv3_apply_fdr_family(
    r,
    "p_value",
    "score_method_id",
    primary$role_table
  )
  rows$within_subbin_residual_zscore <- r
  score_matrices$within_subbin_residual_zscore <- residual$score_matrix
  score_results <- ptpv4_bind_rows(rows)
  ptp_write_csv(
    pb$metadata,
    file.path(dirs$marginal, "marginal_mouse_metadata.csv")
  )
  ptp_write_csv(
    data.frame(
      sample_id = colnames(pb$counts),
      library_size = Matrix::colSums(pb$counts),
      stringsAsFactors = FALSE
    ),
    file.path(dirs$marginal, "marginal_mouse_pseudobulk_qc.csv")
  )
  ptp_write_csv(
    score_results,
    file.path(dirs$marginal, "marginal_five_score_model_results.csv")
  )
  fry <- ptpv4_run_marginal_fry_mroast(
    model,
    primary$programs,
    primary$role_table,
    cfg,
    args,
    dirs$marginal
  )
  union_genes <- sort(unique(unlist(
    lapply(membership$memberships, `[[`, "genes"),
    use.names = FALSE
  )))
  symbols <- ptp_clean_gene_symbols(rownames(logcpm))
  idx <- match(union_genes, symbols)
  keep <- is.finite(idx)
  idx <- idx[keep]
  union_genes <- union_genes[keep]
  E <- logcpm[idx, match(
    assignments$strata |> unlist(use.names = FALSE),
    model$metadata$sample_id
  ), drop = FALSE]
  colnames(E) <- assignments$strata |> unlist(use.names = FALSE)
  rownames(E) <- union_genes
  mouse_meta <- model$metadata[
    match(colnames(E), model$metadata$sample_id),
    ,
    drop = FALSE
  ]
  T <- ptpv4_interaction_stats(E, assignments)$t
  rownames(T) <- union_genes
  adaptive <- ptpv4_run_adaptive_from_matrices(
    T,
    E,
    assignments,
    membership,
    cfg,
    args,
    dirs$marginal,
    prefix = "marginal"
  )
  family <- ptpv4_run_family_hierarchy(
    adaptive,
    membership,
    list(expression = E, metadata = mouse_meta),
    cfg,
    args,
    dirs,
    output_dir = dirs$marginal
  )
  occupancy <- ptp_occupancy_models(primary$metadata, cfg$regions)
  ptp_write_csv(
    occupancy,
    file.path(dirs$marginal, "pseudotime_region_occupancy_interactions.csv")
  )
  adaptive$gene_matrices <- NULL
  list(
    pseudobulk = pb,
    model = model,
    score_results = score_results,
    score_matrices = score_matrices,
    fry_mroast = fry,
    adaptive = adaptive,
    family = family,
    occupancy = occupancy
  )
}

ptpv4_run_endpoint_and_standard_dose <- function(primary, cfg, args, dirs) {
  pb <- primary$primary_pb
  programs <- primary$programs
  dose_fit <- ptp_fit_cellmeans(
    pb$counts,
    pb$metadata,
    pb$eligible_subbins,
    "dose_group",
    cfg$qc
  )
  dose_contrasts <- ptp_dose_contrasts(
    pb$eligible_subbins,
    colnames(dose_fit$design)
  )
  dose_scores <- ptp_program_scores_from_logcpm(
    ptp_logcpm_matrix(dose_fit),
    programs
  )
  dose_program <- ptp_apply_program_contrasts(
    ptp_fit_program_scores(dose_scores, dose_fit),
    dose_contrasts,
    "dose_specific_secondary",
    programs
  )
  dose_program <- ptpv3_apply_fdr_family(
    dose_program,
    "p_value",
    "contrast_id",
    primary$role_table
  )
  ptp_write_csv(
    dose_program,
    file.path(dirs$dose, "standard_dose_specific_score_results.csv")
  )
  endpoint <- ptp_run_endpoint_ploidy_models(
    pb,
    programs,
    cfg,
    args
  )
  if (nrow(endpoint$continuous_program)) {
    endpoint$continuous_program <- ptpv3_apply_fdr_family(
      endpoint$continuous_program,
      "p_value",
      "contrast_id",
      primary$role_table
    )
  }
  if (nrow(endpoint$group_program)) {
    endpoint$group_program <- ptpv3_apply_fdr_family(
      endpoint$group_program,
      "p_value",
      "endpoint_ploidy_method",
      primary$role_table
    )
  }
  ptp_write_csv(
    endpoint$status,
    file.path(dirs$end_time, "endpoint_ploidy_status.csv")
  )
  ptp_write_csv(
    endpoint$continuous_program,
    file.path(dirs$etp, "continuous_etp_program_interactions.csv")
  )
  ptp_write_csv(
    endpoint$group_program,
    file.path(dirs$end_time, "threshold_end_time_ploidy_interactions.csv")
  )
  list(dose = dose_program, endpoint = endpoint)
}

ptpv4_build_dose_stratum_space <- function(ids) {
  ids <- sort(ids)
  controls <- utils::combn(ids, 4L, simplify = FALSE)
  labels <- list()
  rows <- list()
  k <- 0L
  for (control in controls) {
    remaining <- setdiff(ids, control)
    dose30 <- utils::combn(remaining, 2L, simplify = FALSE)
    for (d30 in dose30) {
      k <- k + 1L
      d120 <- setdiff(remaining, d30)
      lab <- rep("dose120", length(ids))
      names(lab) <- ids
      lab[control] <- "control"
      lab[d30] <- "dose30"
      labels[[k]] <- lab
      rows[[k]] <- data.frame(
        stratum_assignment = k,
        control_samples = paste(sort(control), collapse = ";"),
        dose30_samples = paste(sort(d30), collapse = ";"),
        dose120_samples = paste(sort(d120), collapse = ";"),
        checksum = ptpv4_object_sha(lab),
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    table = ptpv4_bind_rows(rows),
    labels = do.call(cbind, labels),
    ids = ids
  )
}

ptpv4_build_dose_assignment_space <- function(mouse_meta) {
  strata <- split(mouse_meta$sample_id, mouse_meta$initial_ploidy)
  space <- lapply(strata, ptpv4_build_dose_stratum_space)
  names(space) <- names(strata)
  n1 <- ncol(space[[1L]]$labels)
  n2 <- ncol(space[[2L]]$labels)
  grid <- expand.grid(
    stratum1_assignment = seq_len(n1),
    stratum2_assignment = seq_len(n2),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid <- grid[order(grid$stratum1_assignment, grid$stratum2_assignment), ]
  grid$assignment_index <- seq_len(nrow(grid))
  grid$assignment_id <- sprintf("dose_perm_%06d", grid$assignment_index)
  observed <- setNames(
    ifelse(
      mouse_meta$dose_mg == 0,
      "control",
      ifelse(mouse_meta$dose_mg == 30, "dose30", "dose120")
    ),
    mouse_meta$sample_id
  )
  obs1 <- which(vapply(seq_len(n1), function(i) {
    identical(
      unname(space[[1L]]$labels[, i]),
      unname(observed[space[[1L]]$ids])
    )
  }, logical(1L)))
  obs2 <- which(vapply(seq_len(n2), function(i) {
    identical(
      unname(space[[2L]]$labels[, i]),
      unname(observed[space[[2L]]$ids])
    )
  }, logical(1L)))
  grid$is_observed_assignment <-
    grid$stratum1_assignment == obs1 &
    grid$stratum2_assignment == obs2
  grid$assignment_checksum <- vapply(seq_len(nrow(grid)), function(i) {
    ptpv4_object_sha(c(
      space[[1L]]$labels[, grid$stratum1_assignment[[i]]],
      space[[2L]]$labels[, grid$stratum2_assignment[[i]]]
    ))
  }, character(1L))
  list(
    table = grid,
    strata = space,
    observed_index = which(grid$is_observed_assignment),
    observed_labels = observed
  )
}

ptpv4_group_mean_variance <- function(Y, membership_matrix) {
  n <- colSums(membership_matrix)
  mean <- sweep(Y %*% membership_matrix, 2L, n, "/")
  ss <- Y^2 %*% membership_matrix - sweep(mean^2, 2L, n, "*")
  variance <- sweep(ss, 2L, pmax(n - 1, 1), "/")
  list(mean = mean, variance = variance, n = n)
}

ptpv4_dose_stratum_statistics <- function(Y, label_matrix) {
  control <- label_matrix == "control"
  d30 <- label_matrix == "dose30"
  d120 <- label_matrix == "dose120"
  anydose <- d30 | d120
  c0 <- ptpv4_group_mean_variance(Y, control)
  c30 <- ptpv4_group_mean_variance(Y, d30)
  c120 <- ptpv4_group_mean_variance(Y, d120)
  ca <- ptpv4_group_mean_variance(Y, anydose)
  contrast <- function(a, b) {
    list(
      estimate = a$mean - b$mean,
      variance = sweep(a$variance, 2L, a$n, "/") +
        sweep(b$variance, 2L, b$n, "/")
    )
  }
  pooled <- contrast(ca, c0)
  d30c <- contrast(c30, c0)
  d120c <- contrast(c120, c0)
  nonlinear <- contrast(c120, c30)
  X <- (d30 * 1) + (d120 * 2)
  xbar <- colMeans(X)
  XC <- sweep(X, 2L, xbar, "-")
  sxx <- colSums(XC^2)
  beta <- sweep(Y %*% XC, 2L, sxx, "/")
  ybar <- rowMeans(Y)
  sst <- rowSums((Y - ybar)^2)
  sse <- matrix(sst, nrow = nrow(Y), ncol = ncol(XC)) -
    sweep(beta^2, 2L, sxx, "*")
  slope_var <- sweep(pmax(sse / 6, 0), 2L, sxx, "/")
  list(
    pooled_any_dose_interaction = pooled,
    linear_dose_trend_interaction = list(
      estimate = beta,
      variance = slope_var
    ),
    dose30_vs_control_interaction = d30c,
    dose120_vs_control_interaction = d120c,
    dose120_vs_30_nonlinear_interaction = nonlinear
  )
}

ptpv4_dose_one_entity <- function(
  y,
  dose_space,
  observed_index,
  assignment_keep = NULL
) {
  z <- lapply(seq_along(dose_space$strata), function(k) {
    ids <- dose_space$strata[[k]]$ids
    ptpv4_dose_stratum_statistics(
      matrix(y[ids], nrow = 1L, dimnames = list(NULL, ids)),
      dose_space$strata[[k]]$labels
    )
  })
  candidates <- names(z[[1L]])
  statistics <- list()
  for (candidate in candidates) {
    e1 <- as.numeric(z[[1L]][[candidate]]$estimate)
    e2 <- as.numeric(z[[2L]][[candidate]]$estimate)
    v1 <- as.numeric(z[[1L]][[candidate]]$variance)
    v2 <- as.numeric(z[[2L]][[candidate]]$variance)
    estimate <- outer(e1, e2, function(a, b) b - a)
    se <- sqrt(outer(v1, v2, "+"))
    tstat <- as.vector(t(estimate / se))
    tstat[!is.finite(tstat)] <- NA_real_
    statistics[[candidate]] <- abs(tstat)
  }
  joint <- sqrt(Reduce(`+`, lapply(statistics, `^`, 2)))
  statistics$joint_dose_by_ploidy_statistic <- joint
  if (!is.null(assignment_keep)) {
    statistics <- lapply(statistics, `[`, assignment_keep)
    observed_index <- match(observed_index, assignment_keep)
  }
  p_by_assignment <- lapply(statistics, ptpv4_empirical_p_upper)
  minp <- apply(do.call(rbind, p_by_assignment), 2L, min, na.rm = TRUE)
  adaptive <- ptpv4_empirical_p_lower(minp)
  component <- data.frame(
    candidate = names(statistics),
    observed_statistic = vapply(statistics, `[[`, numeric(1L), observed_index),
    empirical_p = vapply(
      p_by_assignment,
      `[[`,
      numeric(1L),
      observed_index
    ),
    stringsAsFactors = FALSE
  )
  list(
    component = component,
    minp = minp,
    adaptive_by_assignment = adaptive,
    observed_adaptive_p = adaptive[[observed_index]]
  )
}

ptpv4_run_exact_dose <- function(
  mouse_bundle,
  membership,
  cfg,
  args,
  dirs,
  cfg_sha,
  input_sha
) {
  dose_space <- ptpv4_build_dose_assignment_space(mouse_bundle$metadata)
  if (
    nrow(dose_space$table) != as.integer(cfg$dose$expected_assignments) ||
    length(dose_space$observed_index) != 1L
  ) {
    stop("Dose assignment space failed its frozen cardinality audit.", call. = FALSE)
  }
  ptp_write_csv(
    dose_space$table,
    file.path(dirs$dose, "dose_176400_assignments.csv")
  )
  ptp_write_csv(
    ptpv4_bind_rows(lapply(names(dose_space$strata), function(s) {
      x <- dose_space$strata[[s]]$table
      x$initial_ploidy <- s
      x
    })),
    file.path(dirs$dose, "dose_stratum_420_assignments.csv")
  )
  score <- ptpv4_membership_scalar_scores(mouse_bundle, membership)
  mids <- rownames(score)
  if (isTRUE(args$smoke)) {
    mids <- head(mids, min(8L, length(mids)))
  }
  assignment_keep <- NULL
  if (isTRUE(args$smoke)) {
    assignment_keep <- seq_len(min(
      as.integer(args$smoke_dose_assignments),
      nrow(dose_space$table)
    ))
    if (!(dose_space$observed_index %in% assignment_keep)) {
      assignment_keep[[length(assignment_keep)]] <- dose_space$observed_index
    }
  }
  rows <- list()
  components <- list()
  distributions <- list()
  for (i in seq_along(mids)) {
    mid <- mids[[i]]
    one <- ptpv4_dose_one_entity(
      score[mid, ],
      dose_space,
      dose_space$observed_index,
      assignment_keep
    )
    c0 <- one$component
    c0$membership_id <- mid
    components[[i]] <- c0
    rows[[i]] <- data.frame(
      membership_id = mid,
      method_id = "dose_nested_exact_adaptive",
      adaptive_empirical_p = one$observed_adaptive_p,
      n_assignments = if (is.null(assignment_keep)) {
        nrow(dose_space$table)
      } else {
        length(assignment_keep)
      },
      complete_frozen_assignment_space = is.null(assignment_keep),
      stringsAsFactors = FALSE
    )
    distributions[[mid]] <- list(
      minp = one$minp,
      adaptive_by_assignment = one$adaptive_by_assignment
    )
    if (i %% as.integer(cfg$dose$checkpoint_every_memberships) == 0L) {
      ptpv4_checkpoint_write(
        dirs,
        paste0("dose_memberships_", i),
        list(
          completed = mids[seq_len(i)],
          rows = ptpv4_bind_rows(rows),
          components = ptpv4_bind_rows(components)
        ),
        cfg_sha,
        cfg,
        input_sha
      )
      ptpv4_message(
        "Dose memberships ",
        i,
        "/",
        length(mids),
        log_file = file.path(dirs$logs, "run.log")
      )
    }
  }
  results <- ptpv4_apply_membership_fdr(
    ptpv4_bind_rows(rows),
    membership
  )
  component <- ptpv4_bind_rows(components)
  alias <- merge(
    membership$alias,
    results,
    by = "membership_id",
    all.x = TRUE,
    sort = FALSE
  )
  audit <- data.frame(
    expected_assignments = as.integer(cfg$dose$expected_assignments),
    generated_assignments = nrow(dose_space$table),
    unique_assignment_checksums = length(unique(
      dose_space$table$assignment_checksum
    )),
    observed_assignment_count = sum(dose_space$table$is_observed_assignment),
    enumerated_complete_space = TRUE,
    inference_assignment_count = if (is.null(assignment_keep)) {
      nrow(dose_space$table)
    } else {
      length(assignment_keep)
    },
    smoke_subset = isTRUE(args$smoke),
    stringsAsFactors = FALSE
  )
  ptp_write_csv(
    component,
    file.path(dirs$dose, "dose_adaptive_component_results.csv")
  )
  ptp_write_csv(
    results,
    file.path(dirs$dose, "dose_unique_membership_adaptive_results.csv")
  )
  ptp_write_csv(
    alias,
    file.path(dirs$dose, "dose_program_alias_results.csv")
  )
  ptp_write_csv(
    audit,
    file.path(dirs$dose, "dose_assignment_space_audit.csv")
  )
  saveRDS(
    distributions,
    file.path(dirs$dose, "dose_nested_permutation_distributions.rds"),
    compress = "xz"
  )
  list(
    space = dose_space,
    results = results,
    alias = alias,
    component = component,
    distributions = distributions,
    audit = audit
  )
}

ptpv4_simulation_programs <- function(programs, role_table, cfg) {
  tab <- role_table[match(names(programs), role_table$program_id), , drop = FALSE]
  tab$size_class <- ifelse(
    tab$n_genes_observed <= 50,
    "small",
    ifelse(tab$n_genes_observed <= 250, "medium", "large")
  )
  selected <- character(0L)
  for (s in c("small", "medium", "large")) {
    ids <- sort(tab$program_id[
      tab$size_class == s & tab$tier != "negative_control"
    ])
    selected <- c(selected, head(ids, 3L))
  }
  controls <- sort(tab$program_id[tab$tier == "negative_control"])
  if (length(controls) && length(selected)) {
    selected[[length(selected)]] <- controls[[1L]]
  }
  n_target <- as.integer(cfg$simulation$representative_program_count)
  selected <- head(unique(selected), n_target)
  list(
    programs = programs[selected],
    table = tab[match(selected, tab$program_id), , drop = FALSE]
  )
}

ptpv4_simulation_reduced_model <- function(primary_fit, gene_idx) {
  meta <- primary_fit$metadata
  meta$subbin_id <- factor(meta$subbin_id)
  meta$initial_ploidy <- factor(meta$initial_ploidy, levels = c("2N", "4N"))
  meta$treatment <- factor(meta$treatment, levels = c("control", "gemcitabine"))
  design <- stats::model.matrix(
    ~ subbin_id + initial_ploidy + treatment,
    data = meta
  )
  E <- ptp_logcpm_matrix(primary_fit)[gene_idx, , drop = FALSE]
  q <- qr(design)
  coef <- qr.coef(q, t(E))
  fitted <- t(design %*% coef)
  residual <- E - fitted
  list(
    expression = E,
    fitted = fitted,
    residual = residual,
    design = design,
    metadata = meta,
    design_columns = colnames(design)
  )
}

ptpv4_simulation_mouse_expression <- function(E, meta, operator) {
  subbins <- unique(meta$subbin_id)
  nuisance <- vapply(subbins, function(sb) {
    rowMeans(E[, meta$subbin_id == sb, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(E)))
  colnames(nuisance) <- subbins
  R <- E
  for (j in seq_len(ncol(R))) {
    R[, j] <- E[, j] - nuisance[, as.character(meta$subbin_id[[j]])]
  }
  R %*% t(operator)
}

ptpv4_fast_program_scores <- function(
  mouse_expression,
  programs,
  residual_gene_scale = NULL
) {
  symbols <- ptp_clean_gene_symbols(rownames(mouse_expression))
  Z <- t(scale(t(mouse_expression)))
  Z[!is.finite(Z)] <- 0
  if (is.null(residual_gene_scale)) {
    residual_gene_scale <- apply(mouse_expression, 1L, stats::sd)
  }
  residual_gene_scale[!is.finite(residual_gene_scale) |
    residual_gene_scale <= 0] <- 1
  RZ <- sweep(mouse_expression, 1L, residual_gene_scale, "/")
  methods <- c(
    "signed_weighted_mean_zscore",
    "signed_pca_eigengene",
    "signed_rank_mean",
    "trimmed_signed_mean_zscore",
    "within_subbin_residual_zscore"
  )
  out <- lapply(methods, function(x) {
    matrix(
      NA_real_,
      nrow = length(programs),
      ncol = ncol(mouse_expression),
      dimnames = list(names(programs), colnames(mouse_expression))
    )
  })
  names(out) <- methods
  for (id in names(programs)) {
    p <- programs[[id]]
    idx <- match(p$observed_genes, symbols)
    keep <- is.finite(idx)
    idx <- idx[keep]
    w <- p$observed_gene_weights$signed_weight[keep]
    ok <- is.finite(w) & w != 0
    idx <- idx[ok]
    w <- w[ok]
    if (length(idx) < 2L) {
      next
    }
    wn <- w / sum(abs(w))
    primary <- as.numeric(crossprod(wn, Z[idx, , drop = FALSE]))
    out$signed_weighted_mean_zscore[id, ] <- primary
    signed <- sweep(Z[idx, , drop = FALSE], 1L, sign(w), "*")
    gram <- crossprod(signed)
    eig <- tryCatch(eigen(gram, symmetric = TRUE), error = function(e) NULL)
    if (!is.null(eig) && ncol(eig$vectors)) {
      pc <- eig$vectors[, 1L] * sqrt(pmax(eig$values[[1L]], 0))
      cor0 <- suppressWarnings(stats::cor(pc, primary))
      if (is.finite(cor0) && cor0 < 0) {
        pc <- -pc
      }
      out$signed_pca_eigengene[id, ] <- as.numeric(scale(pc))
    }
    rank_score <- vapply(seq_len(ncol(Z)), function(j) {
      x <- Z[idx, j] * sign(w)
      r <- rank(x, ties.method = "average")
      mean((r - (length(r) + 1) / 2) / pmax((length(r) - 1) / 2, 1))
    }, numeric(1L))
    out$signed_rank_mean[id, ] <- rank_score
    trimmed <- vapply(seq_len(ncol(Z)), function(j) {
      x <- Z[idx, j] * sign(w)
      mean(x, trim = 0.10, na.rm = TRUE)
    }, numeric(1L))
    out$trimmed_signed_mean_zscore[id, ] <- trimmed
    out$within_subbin_residual_zscore[id, ] <- as.numeric(crossprod(
      wn,
      RZ[idx, , drop = FALSE]
    ))
  }
  out
}

ptpv4_simulation_score_statistics <- function(
  scores,
  mouse_meta,
  true_scores = NULL
) {
  treated <- mouse_meta$treatment != "control"
  rows <- list()
  for (method in names(scores)) {
    stat <- ptpv4_single_assignment_stats(scores[[method]], mouse_meta, treated)
    true_effect <- rep(0, nrow(scores[[method]]))
    if (!is.null(true_scores)) {
      true_stat <- ptpv4_single_assignment_stats(
        true_scores[[method]],
        mouse_meta,
        treated
      )
      true_effect <- true_stat$estimate
    }
    df <- 12
    p <- 2 * stats::pt(-abs(stat$t), df = df)
    crit <- stats::qt(0.975, df = df)
    rows[[method]] <- data.frame(
      score_method_id = method,
      program_id = rownames(scores[[method]]),
      estimate = stat$estimate,
      se = stat$se,
      ci_low = stat$estimate - crit * stat$se,
      ci_high = stat$estimate + crit * stat$se,
      p_value = p,
      true_effect = true_effect,
      stringsAsFactors = FALSE
    )
  }
  ptpv4_bind_rows(rows)
}

ptpv4_binomial_ci <- function(k, n, conf = 0.95) {
  if (!is.finite(k) || !is.finite(n) || n <= 0) {
    return(c(NA_real_, NA_real_))
  }
  as.numeric(stats::binom.test(k, n, conf.level = conf)$conf.int)
}

ptpv4_run_simulation <- function(
  primary,
  mouse_bundle,
  cfg,
  args,
  dirs,
  cfg_sha,
  input_sha
) {
  selected <- ptpv4_simulation_programs(
    primary$programs,
    primary$role_table,
    cfg
  )
  programs <- selected$programs
  all_genes <- sort(unique(unlist(
    lapply(programs, `[[`, "observed_genes"),
    use.names = FALSE
  )))
  symbols <- ptp_clean_gene_symbols(rownames(primary$primary_fit$dge$counts))
  idx <- match(all_genes, symbols)
  keep <- is.finite(idx)
  idx <- idx[keep]
  all_genes <- all_genes[keep]
  reduced <- ptpv4_simulation_reduced_model(primary$primary_fit, idx)
  rownames(reduced$expression) <- rownames(reduced$fitted) <-
    rownames(reduced$residual) <- all_genes
  residual_scale <- apply(reduced$residual, 1L, stats::sd)
  scenarios <- ptpv4_bind_rows(lapply(
    cfg$simulation$scenarios,
    as.data.frame,
    stringsAsFactors = FALSE
  ))
  if (isTRUE(args$smoke)) {
    scenarios$replicates <- ifelse(
      scenarios$standardized_effect_size == 0,
      as.integer(args$smoke_null_replicates),
      as.integer(args$smoke_power_replicates)
    )
  }
  ptp_write_csv(
    scenarios,
    file.path(dirs$simulation, "frozen_simulation_scenarios.csv")
  )
  ptp_write_csv(
    selected$table,
    file.path(dirs$simulation, "representative_programs.csv")
  )
  mouse_ids <- unique(reduced$metadata$sample_id)
  obs_to_mouse <- match(reduced$metadata$sample_id, mouse_ids)
  base_mouse <- ptpv4_simulation_mouse_expression(
    reduced$fitted,
    reduced$metadata,
    mouse_bundle$operator
  )
  base_scores <- ptpv4_fast_program_scores(
    base_mouse,
    programs,
    residual_scale
  )
  scenario_worker <- function(si) {
    sc <- scenarios[si, , drop = FALSE]
    checkpoint_id <- paste0("simulation_", sc$id)
    old <- ptpv4_checkpoint_read(
      dirs,
      checkpoint_id,
      cfg_sha,
      cfg,
      input_sha,
      required_fields = c("raw", "fdp")
    )
    if (!is.null(old)) {
      return(old)
    }
    nrep <- as.integer(sc$replicates)
    raw_rows <- list()
    fdp_rows <- list()
    size_class <- as.character(sc$program_size)
    candidate <- selected$table$program_id[
      selected$table$size_class == size_class
    ]
    if (!length(candidate)) {
      candidate <- selected$table$program_id
    }
    for (r in seq_len(nrep)) {
      seed <- as.integer(cfg$reproducibility$simulation_seed) +
        si * 100000L + r
      set.seed(seed)
      signs <- sample(c(-1, 1), length(mouse_ids), replace = TRUE)
      simE <- reduced$fitted +
        sweep(reduced$residual, 2L, signs[obs_to_mouse], "*")
      target <- if (sc$standardized_effect_size == 0) {
        selected$table$program_id[[1L]]
      } else {
        candidate[((r - 1L) %% length(candidate)) + 1L]
      }
      signal <- matrix(
        0,
        nrow = nrow(simE),
        ncol = ncol(simE),
        dimnames = dimnames(simE)
      )
      active <- character(0L)
      if (as.numeric(sc$active_gene_fraction) > 0) {
        genes <- intersect(programs[[target]]$observed_genes, rownames(simE))
        n_active <- max(
          1L,
          min(
            length(genes),
            round(length(genes) * as.numeric(sc$active_gene_fraction))
          )
        )
        active <- sample(genes, n_active, replace = FALSE)
        direction <- if (identical(as.character(sc$direction_pattern), "mixed")) {
          sample(c(-1, 1), n_active, replace = TRUE)
        } else {
          rep(1, n_active)
        }
        if (identical(as.character(sc$direction_pattern), "sparse")) {
          direction <- sample(c(-1, 1), n_active, replace = TRUE)
        }
        active_idx <- match(active, rownames(simE))
        interaction_obs <- reduced$metadata$initial_ploidy == "4N" &
          reduced$metadata$treatment == "gemcitabine"
        effect <- as.numeric(sc$standardized_effect_size) *
          pmax(residual_scale[active_idx], 0.05) * direction
        signal[active_idx, interaction_obs] <- effect
        simE <- simE + signal
      }
      mouseE <- ptpv4_simulation_mouse_expression(
        simE,
        reduced$metadata,
        mouse_bundle$operator
      )
      scores <- ptpv4_fast_program_scores(mouseE, programs, residual_scale)
      noiseless_mouse <- ptpv4_simulation_mouse_expression(
        reduced$fitted + signal,
        reduced$metadata,
        mouse_bundle$operator
      )
      noiseless_scores <- ptpv4_fast_program_scores(
        noiseless_mouse,
        programs,
        residual_scale
      )
      true_scores <- Map(`-`, noiseless_scores, base_scores)
      stat <- ptpv4_simulation_score_statistics(
        scores,
        mouse_bundle$metadata,
        true_scores
      )
      stat$scenario_id <- sc$id
      stat$replicate <- r
      stat$replicate_seed <- seed
      stat$target_program_id <- target
      stat$is_target_program <- stat$program_id == target &
        as.numeric(sc$standardized_effect_size) != 0
      stat$is_null_program <- !stat$is_target_program
      stat$active_genes <- paste(sort(active), collapse = ";")
      stat$n_active_genes <- length(active)
      stat$program_size_class <- size_class
      stat$active_gene_fraction <- as.numeric(sc$active_gene_fraction)
      stat$direction_pattern <- as.character(sc$direction_pattern)
      stat$standardized_effect_size <- as.numeric(
        sc$standardized_effect_size
      )
      stat$reject_0_05 <- is.finite(stat$p_value) & stat$p_value < 0.05
      stat$covered <- is.finite(stat$ci_low) &
        stat$ci_low <= stat$true_effect &
        stat$ci_high >= stat$true_effect
      stat$direction_correct <- ifelse(
        stat$is_target_program & stat$true_effect != 0,
        sign(stat$estimate) == sign(stat$true_effect),
        NA
      )
      raw_rows[[r]] <- stat
      fdp_method <- lapply(split(stat, stat$score_method_id), function(x) {
        q <- stats::p.adjust(x$p_value, "BH")
        rejected <- is.finite(q) & q < as.numeric(cfg$simulation$q)
        false <- rejected & x$is_null_program
        data.frame(
          scenario_id = sc$id,
          replicate = r,
          score_method_id = x$score_method_id[[1L]],
          discoveries_q_0_10 = sum(rejected),
          false_discoveries_q_0_10 = sum(false),
          replicate_FDP = sum(false) / max(sum(rejected), 1L),
          target_rejected_q_0_10 = any(rejected & x$is_target_program),
          replicate_failure = any(!is.finite(x$p_value)),
          stringsAsFactors = FALSE
        )
      })
      fdp_rows[[r]] <- ptpv4_bind_rows(fdp_method)
    }
    out <- list(
      raw = ptpv4_bind_rows(raw_rows),
      fdp = ptpv4_bind_rows(fdp_rows)
    )
    ptpv4_checkpoint_write(
      dirs,
      checkpoint_id,
      out,
      cfg_sha,
      cfg,
      input_sha
    )
    out
  }
  workers <- ptpv2_workers(
    args,
    nrow(scenarios),
    "simulation_workers"
  )
  parallel <- ptp_parallel_lapply(
    seq_len(nrow(scenarios)),
    scenario_worker,
    workers = workers,
    task_label = "v4_simulation"
  )
  raw <- ptpv4_bind_rows(lapply(parallel$results, `[[`, "raw"))
  fdp <- ptpv4_bind_rows(lapply(parallel$results, `[[`, "fdp"))
  sentinel <- selected$table$program_id[[1L]]
  summary_rows <- list()
  groups <- split(
    raw,
    interaction(raw$scenario_id, raw$score_method_id, drop = TRUE)
  )
  for (x in groups) {
    scenario_id <- x$scenario_id[[1L]]
    method <- x$score_method_id[[1L]]
    sc <- scenarios[match(scenario_id, scenarios$id), , drop = FALSE]
    is_null <- as.numeric(sc$standardized_effect_size) == 0
    replicate_rows <- if (is_null) {
      x[x$program_id == sentinel, , drop = FALSE]
    } else {
      x[x$is_target_program, , drop = FALSE]
    }
    n <- nrow(replicate_rows)
    reject_rate <- mean(replicate_rows$reject_0_05, na.rm = TRUE)
    ci <- ptpv4_binomial_ci(
      sum(replicate_rows$reject_0_05, na.rm = TRUE),
      sum(is.finite(replicate_rows$p_value))
    )
    f <- fdp[
      fdp$scenario_id == scenario_id &
        fdp$score_method_id == method,
      ,
      drop = FALSE
    ]
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      scenario_id = scenario_id,
      score_method_id = method,
      replicates = n,
      type_I_error = if (is_null) reject_rate else NA_real_,
      type_I_mc_ci_low = if (is_null) ci[[1L]] else NA_real_,
      type_I_mc_ci_high = if (is_null) ci[[2L]] else NA_real_,
      type_I_ci_covers_0_05 = if (is_null) {
        ci[[1L]] <= 0.05 && ci[[2L]] >= 0.05
      } else {
        NA
      },
      type_I_in_0035_0065_reference = if (is_null) {
        reject_rate >= 0.035 && reject_rate <= 0.065
      } else {
        NA
      },
      empirical_FDR_q_0_10 = mean(f$replicate_FDP, na.rm = TRUE),
      global_null_fdr_le_0_12 = if (is_null) {
        mean(f$replicate_FDP, na.rm = TRUE) <= 0.12
      } else {
        NA
      },
      power_0_05 = if (!is_null) reject_rate else NA_real_,
      power_q_0_10 = if (!is_null) {
        mean(f$target_rejected_q_0_10, na.rm = TRUE)
      } else {
        NA_real_
      },
      bias = mean(
        replicate_rows$estimate - replicate_rows$true_effect,
        na.rm = TRUE
      ),
      ci_coverage = mean(replicate_rows$covered, na.rm = TRUE),
      coverage_in_0925_0975 = mean(
        replicate_rows$covered,
        na.rm = TRUE
      ) >= 0.925 && mean(
        replicate_rows$covered,
        na.rm = TRUE
      ) <= 0.975,
      direction_accuracy = if (!is_null) {
        mean(replicate_rows$direction_correct, na.rm = TRUE)
      } else {
        NA_real_
      },
      failure_rate = mean(f$replicate_failure, na.rm = TRUE),
      mcse_reject = sqrt(
        reject_rate * (1 - reject_rate) /
          max(sum(is.finite(replicate_rows$p_value)), 1L)
      ),
      score_specific_true_effect_mean = mean(
        replicate_rows$true_effect,
        na.rm = TRUE
      ),
      score_specific_true_effect_sd = stats::sd(
        replicate_rows$true_effect,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }
  summary <- ptpv4_bind_rows(summary_rows)
  null_summary <- summary[is.finite(summary$type_I_error), , drop = FALSE]
  method_status <- ptpv4_bind_rows(lapply(
    unique(summary$score_method_id),
    function(method) {
      n0 <- null_summary[null_summary$score_method_id == method, , drop = FALSE]
      power_rows <- summary[
        summary$score_method_id == method &
          is.finite(summary$power_0_05),
        ,
        drop = FALSE
      ]
      calibrated <- nrow(n0) == 1L &&
        isTRUE(n0$type_I_ci_covers_0_05) &&
        isTRUE(n0$global_null_fdr_le_0_12) &&
        all(power_rows$coverage_in_0925_0975 %in% TRUE) &&
        all(power_rows$failure_rate <= 0.05, na.rm = TRUE)
      data.frame(
        method_id = method,
        status = if (calibrated) "calibrated" else "uncalibrated",
        type_I_ci_covers_0_05 = n0$type_I_ci_covers_0_05 %||% FALSE,
        global_null_fdr_le_0_12 = n0$global_null_fdr_le_0_12 %||% FALSE,
        all_power_coverage_in_0925_0975 =
          all(power_rows$coverage_in_0925_0975 %in% TRUE),
        maximum_failure_rate = max(
          summary$failure_rate[summary$score_method_id == method],
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  ))
  audit <- data.frame(
    reduced_model_terms = paste(colnames(reduced$design), collapse = ";"),
    interaction_in_null_fitted_mean = FALSE,
    residual_resampling_unit = "whole_mouse",
    residual_resampling_method = "wild_sign_flip",
    preserves_within_mouse_subbin_correlation = TRUE,
    preserves_missing_subbin_pattern = TRUE,
    preserves_gene_covariance = TRUE,
    preserves_residual_scale = TRUE,
    one_target_program_per_power_replicate = TRUE,
    active_genes_random_seeded_without_replacement = TRUE,
    score_true_effect_recomputed_from_noiseless_signal = TRUE,
    empirical_fdr_replicate_unit = TRUE,
    stringsAsFactors = FALSE
  )
  ptp_write_csv(
    raw,
    file.path(dirs$simulation, "simulation_replicate_program_results.csv")
  )
  ptp_write_csv(
    fdp,
    file.path(dirs$simulation, "simulation_replicate_fdp.csv")
  )
  ptp_write_csv(
    summary,
    file.path(dirs$simulation, "simulation_method_calibration_summary.csv")
  )
  ptp_write_csv(
    method_status,
    file.path(dirs$simulation, "method_calibration_status.csv")
  )
  ptp_write_csv(
    audit,
    file.path(dirs$simulation, "simulation_generation_audit.csv")
  )
  list(
    raw = raw,
    fdp = fdp,
    summary = summary,
    method_status = method_status,
    audit = audit,
    representative_programs = selected$table,
    parallel_log = parallel$log
  )
}

ptpv4_training_tmm_reference <- function(counts) {
  counts <- as.matrix(counts)
  lib <- colSums(counts)
  calc_quantile <- get(".calcFactorQuantile", asNamespace("edgeR"))
  calc_tmm <- get(".calcFactorTMM", asNamespace("edgeR"))
  f75 <- suppressWarnings(calc_quantile(counts, lib.size = lib, p = 0.75))
  ref <- if (stats::median(f75) < 1e-20) {
    which.max(colSums(sqrt(counts)))
  } else {
    which.min(abs(f75 - mean(f75)))
  }
  raw <- vapply(seq_len(ncol(counts)), function(i) {
    calc_tmm(
      obs = counts[, i],
      ref = counts[, ref],
      libsize.obs = lib[[i]],
      libsize.ref = lib[[ref]],
      logratioTrim = 0.3,
      sumTrim = 0.05,
      doWeighting = TRUE,
      Acutoff = -1e10
    )
  }, numeric(1L))
  scale_factor <- exp(mean(log(raw[is.finite(raw) & raw > 0])))
  factors <- raw / scale_factor
  names(factors) <- colnames(counts)
  list(
    reference_column = ref,
    reference_id = colnames(counts)[[ref]],
    reference_counts = counts[, ref],
    reference_library_size = lib[[ref]],
    scale_factor = scale_factor,
    training_factors = factors
  )
}

ptpv4_project_tmm_factor <- function(counts, reference) {
  calc_tmm <- get(".calcFactorTMM", asNamespace("edgeR"))
  lib <- sum(counts)
  raw <- calc_tmm(
    obs = as.numeric(counts),
    ref = as.numeric(reference$reference_counts),
    libsize.obs = lib,
    libsize.ref = reference$reference_library_size,
    logratioTrim = 0.3,
    sumTrim = 0.05,
    doWeighting = TRUE,
    Acutoff = -1e10
  )
  raw / reference$scale_factor
}

ptpv4_logcpm_with_factors <- function(counts, factors) {
  counts <- as.matrix(counts)
  lib <- colSums(counts) * factors
  log2(sweep(counts + 0.5, 2L, lib + 1, "/") * 1e6)
}

ptpv4_run_crossfit_04i <- function(
  primary,
  cfg,
  args,
  repo_root,
  dirs,
  cfg_sha,
  input_sha
) {
  checkpoint <- ptpv4_checkpoint_read(
    dirs,
    "crossfit_04i",
    cfg_sha,
    cfg,
    input_sha,
    required_fields = c("status", "fold_audit", "weights", "scores", "results")
  )
  if (
    !is.null(checkpoint) &&
      nrow(checkpoint$status) > 0L &&
      all(checkpoint$status$status == "ok")
  ) {
    return(checkpoint)
  }
  meta <- primary$metadata
  counts <- primary$counts
  gene_universe <- ptp_clean_gene_symbols(rownames(counts))
  gene_universe_sha <- ptpv4_object_sha(sort(unique(gene_universe)))
  pst_cfg <- pst_read_config(file.path(
    repo_root,
    "Code/in-vivo/04i_pseudotime_state_pathways_config.yaml"
  ))
  model_spec <- list(
    model_id = "primary_initial_ploidy",
    covariate_terms = c("initial_ploidy_factor"),
    covariate_mode = "initial_ploidy",
    etp_method = NA_character_,
    etp_threshold = NA_real_
  )
  sample_ids <- sort(unique(meta$sample_id))
  run_sample_ids <- if (isTRUE(args$smoke)) {
    head(sample_ids, min(2L, length(sample_ids)))
  } else {
    sample_ids
  }
  families <- sort(unique(primary$role_table$response_family[
    primary$role_table$tier %in% c(
      "tier1_core",
      "tier2_upgraded_focused",
      "tier3_exploratory_only"
    )
  ]))
  primary_meta <- primary$primary_fit$metadata
  primary_pb_counts <- primary$primary_pb$counts[
    ,
    primary_meta$sample_subbin_id,
    drop = FALSE
  ]
  score_mat <- matrix(
    NA_real_,
    nrow = length(families),
    ncol = nrow(primary_meta),
    dimnames = list(
      paste0("crossfit_04i_", ptp_safe_name(families)),
      primary_meta$sample_subbin_id
    )
  )
  worker <- function(held) {
    held_cols <- which(primary_meta$sample_id == held)
    train_meta <- meta[meta$sample_id != held, , drop = FALSE]
    train_counts <- counts[, train_meta$cell_id, drop = FALSE]
    step <- "construct_training_04i_pseudobulk"
    tryCatch({
      pb04i <- pst_construct_pseudobulk(
        train_meta,
        train_counts,
        20L,
        5L
      )
      step <- "evaluate_training_04i_coverage"
      coverage <- pst_region_coverage(
        pb04i$cell_metadata,
        pb04i$metadata,
        pst_cfg
      )
      coverage_check <- pst_check_primary_coverage(coverage)
      step <- "fit_training_04i_model"
      mf <- pst_fit_pseudotime_model(
        pb04i$counts,
        pb04i$metadata,
        5L,
        coverage_check$n_contributing_mice[[1L]],
        include_dose = TRUE,
        model_spec = model_spec
      )
      step <- "derive_training_04i_contrast"
      contrast <- pst_contrast_vector(
        mf,
        pst_cfg,
        "primary_adjacent_state",
        501L
      )
      step <- "extract_training_04i_gene_statistics"
      genes <- pst_contrast_table(
        mf,
        contrast,
        "primary_adjacent_state"
      )
      genes$gene_symbol <- ptp_clean_gene_symbols(genes$gene_symbol)
      train_cols <- which(primary_meta$sample_id != held)
      held_raw <- primary_pb_counts[, held_cols, drop = FALSE]
      train_raw <- primary_pb_counts[, train_cols, drop = FALSE]
      step <- "filter_training_only_expression_genes"
      keep_gene <- Matrix::rowSums(train_raw) > 0
      train_raw <- train_raw[keep_gene, , drop = FALSE]
      held_raw <- held_raw[keep_gene, , drop = FALSE]
      step <- "fit_training_only_tmm_reference"
      training_reference <- ptpv4_training_tmm_reference(train_raw)
      step <- "normalize_training_only_expression"
      train_log <- ptpv4_logcpm_with_factors(
        train_raw,
        training_reference$training_factors
      )
      step <- "project_held_out_tmm_factors"
      held_factors <- vapply(seq_len(ncol(held_raw)), function(j) {
        ptpv4_project_tmm_factor(
          held_raw[, j],
          training_reference
        )
      }, numeric(1L))
      step <- "normalize_held_out_expression"
      held_log <- ptpv4_logcpm_with_factors(held_raw, held_factors)
      step <- "fit_training_only_center_and_scale"
      training_mean <- rowMeans(train_log, na.rm = TRUE)
      training_sd <- apply(train_log, 1L, stats::sd)
      training_sd[!is.finite(training_sd) | training_sd <= 0] <- NA_real_
      held_z <- sweep(
        sweep(held_log, 1L, training_mean, "-"),
        1L,
        training_sd,
        "/"
      )
      held_z[!is.finite(held_z)] <- NA_real_
      symbols <- ptp_clean_gene_symbols(rownames(held_z))
      step <- "fit_and_project_family_weights"
      fold_scores <- matrix(
        NA_real_,
        nrow = length(families),
        ncol = length(held_cols),
        dimnames = list(
          paste0("crossfit_04i_", ptp_safe_name(families)),
          primary_meta$sample_subbin_id[held_cols]
        )
      )
      fold_weights <- list()
      for (family in families) {
        fam_programs <- primary$role_table$program_id[
          primary$role_table$response_family == family
        ]
        fam_genes <- sort(unique(unlist(
          lapply(primary$programs[fam_programs], `[[`, "observed_genes"),
          use.names = FALSE
        )))
        g <- genes[
          genes$gene_symbol %in% fam_genes &
            is.finite(genes$t_statistic),
          ,
          drop = FALSE
        ]
        g <- g[!duplicated(g$gene_symbol), , drop = FALSE]
        idx <- match(g$gene_symbol, symbols)
        keep <- is.finite(idx) & is.finite(g$t_statistic)
        idx <- idx[keep]
        g <- g[keep, , drop = FALSE]
        if (length(idx) < 2L) {
          next
        }
        w <- g$t_statistic
        winsor <- stats::quantile(
          w,
          probs = c(0.01, 0.99),
          na.rm = TRUE,
          names = FALSE
        )
        w <- pmin(pmax(w, winsor[[1L]]), winsor[[2L]])
        w <- w / sum(abs(w), na.rm = TRUE)
        score_id <- paste0("crossfit_04i_", ptp_safe_name(family))
        fold_scores[score_id, ] <- as.numeric(crossprod(
          w,
          held_z[idx, , drop = FALSE]
        ))
        fold_weights[[length(fold_weights) + 1L]] <- data.frame(
          held_out_mouse_id = held,
          response_family = family,
          gene_symbol = symbols[idx],
          signed_weight = w,
          source_training_contrast = "04i_primary_adjacent_state",
          stringsAsFactors = FALSE
        )
      }
      weights <- ptpv4_bind_rows(fold_weights)
      normalization_audit <- data.frame(
        held_out_mouse_id = held,
        held_out_sample_subbin_id = primary_meta$sample_subbin_id[held_cols],
        projected_tmm_factor = held_factors,
        training_reference_sample_subbin =
          training_reference$reference_id,
        training_reference_checksum = ptpv4_object_sha(
          training_reference$reference_counts
        ),
        training_mean_checksum = ptpv4_object_sha(training_mean),
        training_sd_checksum = ptpv4_object_sha(training_sd),
        stringsAsFactors = FALSE
      )
      list(
        status = data.frame(
          held_out_mouse_id = held,
          training_mouse_ids = paste(
            sort(unique(train_meta$sample_id)),
            collapse = ";"
          ),
          training_contains_heldout = held %in% unique(train_meta$sample_id),
          training_input_checksum = ptpv4_object_sha(sort(train_meta$cell_id)),
          training_gene_universe = "original_04i_full_target_cell_gene_universe",
          training_gene_universe_n = nrow(counts),
          training_gene_universe_checksum = gene_universe_sha,
          training_04i_retained_genes = nrow(genes),
          training_04i_gene_filter_checksum = ptpv4_object_sha(
            sort(genes$gene_symbol)
          ),
          fitted_weight_checksum = ptpv4_object_sha(weights),
          held_out_projection_checksum = ptpv4_object_sha(fold_scores),
          normalization_training_only = TRUE,
          gene_filtering_training_only = TRUE,
          centering_training_only = TRUE,
          scaling_training_only = TRUE,
          weight_fitting_training_only = TRUE,
          reused_full_data_04i_statistics = FALSE,
          reused_full_data_04j_gene_z = FALSE,
          status = "ok",
          error = "",
          stringsAsFactors = FALSE
        ),
        weights = weights,
        scores = fold_scores,
        normalization = normalization_audit
      )
    }, error = function(e) {
      list(
        status = data.frame(
          held_out_mouse_id = held,
          training_mouse_ids = paste(
            sort(setdiff(sample_ids, held)),
            collapse = ";"
          ),
          training_contains_heldout = FALSE,
          training_input_checksum = "",
          training_gene_universe = "original_04i_full_target_cell_gene_universe",
          training_gene_universe_n = nrow(counts),
          training_gene_universe_checksum = gene_universe_sha,
          training_04i_retained_genes = NA_integer_,
          training_04i_gene_filter_checksum = "",
          fitted_weight_checksum = "",
          held_out_projection_checksum = "",
          normalization_training_only = TRUE,
          gene_filtering_training_only = TRUE,
          centering_training_only = TRUE,
          scaling_training_only = TRUE,
          weight_fitting_training_only = TRUE,
          reused_full_data_04i_statistics = FALSE,
          reused_full_data_04j_gene_z = FALSE,
          status = "failed",
          error = paste(step, conditionMessage(e), sep = ": "),
          stringsAsFactors = FALSE
        ),
        weights = data.frame(),
        scores = matrix(
          NA_real_,
          nrow = length(families),
          ncol = length(held_cols),
          dimnames = list(
            paste0("crossfit_04i_", ptp_safe_name(families)),
            primary_meta$sample_subbin_id[held_cols]
          )
        ),
        normalization = data.frame()
      )
    })
  }
  workers <- ptpv2_workers(args, length(run_sample_ids), "crossfit_workers")
  parallel <- ptp_parallel_lapply(
    run_sample_ids,
    worker,
    workers = workers,
    task_label = "v4_crossfit_04i"
  )
  for (x in parallel$results) {
    score_mat[rownames(x$scores), colnames(x$scores)] <- x$scores
  }
  fold_audit <- ptpv4_bind_rows(lapply(parallel$results, `[[`, "status"))
  weights <- ptpv4_bind_rows(lapply(parallel$results, `[[`, "weights"))
  normalization <- ptpv4_bind_rows(lapply(
    parallel$results,
    `[[`,
    "normalization"
  ))
  results <- data.frame()
  complete <- rowSums(is.finite(score_mat)) == ncol(score_mat)
  if (any(complete)) {
    projection_programs <- lapply(rownames(score_mat)[complete], function(id) {
      list(
        id = id,
        label = id,
        response_family = id,
        response_family_label = id,
        tier = "projection",
        tier_label = "Internal cross-fitted exploratory bridge",
        family = id,
        upgrade_status = "exploratory_bridge",
        expected_direction = "exploratory",
        score_mode = "crossfit_projection",
        score_direction_label = "training-fold 04i state contrast",
        score_direction_source = "04i refit excluding held-out mouse",
        score_component_policy = "winsorized_l1",
        source = "04i_crossfit",
        set_names = character(0L),
        missing_set_names = character(0L),
        genes = character(0L),
        estimable = TRUE,
        component_gene_weights = data.frame(),
        observed_genes = character(0L),
        observed_gene_weights = data.frame(
          gene = character(0L),
          signed_weight = numeric(0L)
        )
      )
    })
    names(projection_programs) <- rownames(score_mat)[complete]
    pfit <- ptp_fit_program_scores(
      score_mat[complete, , drop = FALSE],
      primary$primary_fit
    )
    results <- ptp_apply_program_contrasts(
      pfit,
      ptp_primary_contrasts(
        primary$primary_pb$eligible_subbins,
        colnames(primary$primary_fit$design)
      ),
      "cross_fitted_04i_state_projection",
      projection_programs
    )
    results$validation_role <-
      "internal_exploratory_bridge_not_independent_validation"
  }
  status <- data.frame(
    method_id = "cross_fitted_04i_state_projection",
    status = if (
      all(fold_audit$status == "ok") &&
      all(!fold_audit$training_contains_heldout) &&
      all(fold_audit$normalization_training_only) &&
      all(!fold_audit$reused_full_data_04i_statistics) &&
      all(!fold_audit$reused_full_data_04j_gene_z)
    ) {
      "ok"
    } else {
      "failed"
    },
    n_folds = nrow(fold_audit),
    n_failed = sum(fold_audit$status != "ok"),
    validation_role = "not_independent_validation_or_positive_control",
    stringsAsFactors = FALSE
  )
  score_long <- data.frame(
    program_id = rep(rownames(score_mat), times = ncol(score_mat)),
    sample_subbin_id = rep(colnames(score_mat), each = nrow(score_mat)),
    score = as.numeric(score_mat),
    stringsAsFactors = FALSE
  )
  ptp_write_csv(
    status,
    file.path(dirs$crossfit, "cross_fitted_04i_projection_status.csv")
  )
  ptp_write_csv(
    fold_audit,
    file.path(dirs$crossfit, "cross_fitted_04i_fold_audit.csv")
  )
  ptp_write_csv(
    normalization,
    file.path(dirs$crossfit, "cross_fitted_04i_normalization_audit.csv")
  )
  ptp_write_csv(
    weights,
    file.path(dirs$crossfit, "cross_fitted_04i_fold_weights.csv")
  )
  ptp_write_csv(
    score_long,
    file.path(dirs$crossfit, "cross_fitted_04i_out_of_fold_scores.csv")
  )
  ptp_write_csv(
    results,
    file.path(dirs$crossfit, "cross_fitted_04i_projection_model_results.csv")
  )
  out <- list(
    status = status,
    fold_audit = fold_audit,
    normalization = normalization,
    weights = weights,
    scores = score_long,
    score_matrix = score_mat,
    results = results,
    parallel_log = parallel$log
  )
  ptpv4_checkpoint_write(
    dirs,
    "crossfit_04i",
    out,
    cfg_sha,
    cfg,
    input_sha
  )
  out
}

ptpv4_log_normalize_cells <- function(counts) {
  lib <- Matrix::colSums(counts)
  lib[!is.finite(lib) | lib <= 0] <- 1
  Matrix::t(Matrix::t(counts) / lib) * 1e4
}

ptpv4_select_variable_genes <- function(
  counts,
  n_features = 1500L,
  excluded_symbols = character(0L)
) {
  symbols <- ptp_clean_gene_symbols(rownames(counts))
  keep <- !(symbols %in% excluded_symbols) & Matrix::rowSums(counts) > 0
  x <- ptpv4_log_normalize_cells(counts[keep, , drop = FALSE])
  x@x <- log1p(x@x)
  mu <- Matrix::rowMeans(x)
  second <- Matrix::rowMeans(x^2)
  variance <- pmax(second - mu^2, 0)
  ord <- order(variance, decreasing = TRUE, na.last = NA)
  ord <- head(ord, min(as.integer(n_features), length(ord)))
  data.frame(
    gene = rownames(x)[ord],
    gene_symbol = symbols[keep][ord],
    variance = variance[ord],
    stringsAsFactors = FALSE
  )
}

ptpv4_scale_01 <- function(x) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  names(out) <- names(x)
  if (sum(ok) < 2L || diff(range(x[ok])) <= 0) {
    return(out)
  }
  out[ok] <- (x[ok] - min(x[ok])) / diff(range(x[ok]))
  out
}

ptpv4_learn_monocle3_training <- function(
  counts,
  meta,
  variable_genes,
  root_cluster = "6",
  dimensions = 30L
) {
  ids <- intersect(meta$cell_id, colnames(counts))
  meta <- meta[match(ids, meta$cell_id), , drop = FALSE]
  genes <- intersect(variable_genes$gene, rownames(counts))
  x <- counts[genes, ids, drop = FALSE]
  cell_data <- meta
  rownames(cell_data) <- cell_data$cell_id
  gene_data <- data.frame(
    gene_short_name = ptp_clean_gene_symbols(rownames(x)),
    row.names = rownames(x),
    stringsAsFactors = FALSE
  )
  cds <- monocle3::new_cell_data_set(
    x,
    cell_metadata = cell_data,
    gene_metadata = gene_data
  )
  cds <- monocle3::preprocess_cds(
    cds,
    num_dim = min(as.integer(dimensions), nrow(x) - 1L, ncol(x) - 1L),
    norm_method = "log",
    method = "PCA"
  )
  cds <- monocle3::reduce_dimension(
    cds,
    reduction_method = "UMAP",
    preprocess_method = "PCA"
  )
  cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP")
  cds <- monocle3::learn_graph(
    cds,
    use_partition = FALSE,
    learn_graph_control = list(
      minimal_branch_len = 10,
      geodesic_distance_ratio = 0.5
    )
  )
  root_cells <- meta$cell_id[as.character(meta$cluster) == root_cluster]
  root_cells <- intersect(root_cells, colnames(cds))
  if (!length(root_cells)) {
    stop("No root-cluster cells are present in the training set.", call. = FALSE)
  }
  cds <- monocle3::order_cells(
    cds,
    reduction_method = "UMAP",
    root_cells = root_cells
  )
  pt <- monocle3::pseudotime(cds)
  pt <- pt[colnames(cds)]
  scaled <- ptpv4_scale_01(pt)
  if (mean(is.finite(scaled)) < 0.80) {
    stop(
      "Monocle3 yielded finite pseudotime for fewer than 80% of training cells.",
      call. = FALSE
    )
  }
  list(
    cds = cds,
    pseudotime = scaled,
    raw_pseudotime = pt,
    root_cells = root_cells,
    variable_genes = genes
  )
}

ptpv4_project_training_pseudotime <- function(
  counts,
  training_ids,
  all_ids,
  genes,
  training_pseudotime,
  dimensions = 30L,
  k = 15L
) {
  x <- counts[genes, all_ids, drop = FALSE]
  x <- ptpv4_log_normalize_cells(x)
  x@x <- log1p(x@x)
  dense_train <- t(as.matrix(x[, match(training_ids, all_ids), drop = FALSE]))
  pc <- stats::prcomp(
    dense_train,
    center = TRUE,
    scale. = FALSE,
    rank. = min(
      as.integer(dimensions),
      nrow(dense_train) - 1L,
      ncol(dense_train) - 1L
    )
  )
  dense_all <- t(as.matrix(x))
  centered <- sweep(dense_all, 2L, pc$center, "-")
  projected <- centered %*% pc$rotation
  train_pc <- projected[match(training_ids, all_ids), , drop = FALSE]
  knn <- FNN::get.knnx(
    train_pc,
    projected,
    k = min(as.integer(k), nrow(train_pc))
  )
  distance <- pmax(knn$nn.dist, 1e-8)
  weight <- 1 / distance
  projected_pt <- rowSums(
    weight * matrix(
      training_pseudotime[knn$nn.index],
      nrow = nrow(knn$nn.index)
    ),
    na.rm = TRUE
  ) / rowSums(weight)
  names(projected_pt) <- all_ids
  projected_pt[training_ids] <- training_pseudotime[training_ids]
  list(
    pseudotime = ptpv4_scale_01(projected_pt),
    pca = pc,
    projected_pcs = projected,
    knn_index = knn$nn.index,
    knn_distance = knn$nn.dist
  )
}

ptpv4_coordinate_diagnostics <- function(
  coordinate_id,
  coordinate,
  meta,
  training_ids,
  excluded_genes,
  root_cells,
  projection_method,
  projection_checksum,
  treatment_label_audit
) {
  ids <- meta$cell_id
  x <- data.frame(
    coordinate_id = coordinate_id,
    cell_id = ids,
    sample_id = meta$sample_id,
    cluster = meta$cluster,
    initial_ploidy = meta$initial_ploidy,
    treatment = meta$treatment,
    original_pseudotime = meta$pseudotime,
    coordinate = coordinate[ids],
    training_cell = ids %in% training_ids,
    root_cell = ids %in% root_cells,
    stringsAsFactors = FALSE
  )
  grid <- seq(0, 1, length.out = 201L)
  ecdf_rows <- ptpv4_bind_rows(lapply(
    split(x, interaction(x$initial_ploidy, x$treatment, drop = TRUE)),
    function(g) {
      data.frame(
        coordinate_id = coordinate_id,
        group = paste(g$initial_ploidy[[1L]], g$treatment[[1L]], sep = "__"),
        grid = grid,
        coordinate_ecdf = stats::ecdf(g$coordinate[is.finite(g$coordinate)])(grid),
        original_ecdf = stats::ecdf(
          g$original_pseudotime[is.finite(g$original_pseudotime)]
        )(grid),
        stringsAsFactors = FALSE
      )
    }
  ))
  density_rows <- ptpv4_bind_rows(lapply(
    split(x, interaction(x$initial_ploidy, x$treatment, drop = TRUE)),
    function(g) {
      values <- g$coordinate[is.finite(g$coordinate)]
      if (length(values) < 2L || stats::sd(values) == 0) {
        return(data.frame())
      }
      d <- stats::density(values, from = 0, to = 1, n = 201L)
      data.frame(
        coordinate_id = coordinate_id,
        group = paste(g$initial_ploidy[[1L]], g$treatment[[1L]], sep = "__"),
        coordinate = d$x,
        density = d$y,
        stringsAsFactors = FALSE
      )
    }
  ))
  summary <- data.frame(
    coordinate_id = coordinate_id,
    status = if (mean(is.finite(x$coordinate)) >= 0.80) "ok" else "not_estimable",
    n_cells = nrow(x),
    n_training_cells = length(training_ids),
    n_training_mice = length(unique(meta$sample_id[ids %in% training_ids])),
    finite_coordinate_fraction = mean(is.finite(x$coordinate)),
    correlation_with_original = suppressWarnings(stats::cor(
      x$coordinate,
      x$original_pseudotime,
      use = "complete.obs"
    )),
    excluded_gene_count = length(excluded_genes),
    excluded_gene_checksum = ptpv4_object_sha(sort(excluded_genes)),
    root_cell_count = length(root_cells),
    root_cell_checksum = ptpv4_object_sha(sort(root_cells)),
    orientation = "root_cluster_6_to_later_pseudotime",
    projection_method = projection_method,
    projection_checksum = projection_checksum,
    treatment_labels_used_for_training_subset =
      treatment_label_audit$used_for_training_subset,
    treatment_labels_used_for_orientation =
      treatment_label_audit$used_for_orientation,
    treatment_response_used_for_direction_selection = FALSE,
    S_Score_status = if ("S.Score" %in% names(meta)) "available" else "not_available",
    G2M_Score_status = if ("G2M.Score" %in% names(meta)) "available" else "not_available",
    phase_status = if ("Phase" %in% names(meta)) "available" else "not_available",
    stringsAsFactors = FALSE
  )
  list(cells = x, summary = summary, ecdf = ecdf_rows, density = density_rows)
}

ptpv4_coordinate_dirs <- function(root) {
  d <- function(...) ptp_ensure_dir(file.path(root, ...))
  list(
    root = ptp_ensure_dir(root),
    score_models = d("score_models"),
    fry_mroast = d("fry_mroast"),
    mouse_gene = d("mouse_level_gene_statistics"),
    adaptive = d("adaptive_gene_set"),
    family = d("response_family_hierarchical"),
    multiscale = d("multiscale_pseudotime"),
    trajectory = d("trajectory_global"),
    qc = d("qc")
  )
}

ptpv4_run_coordinate_branch <- function(
  coordinate_id,
  coordinate,
  primary,
  membership,
  assignments,
  cfg,
  args,
  dirs
) {
  root <- file.path(dirs$coordinates, coordinate_id)
  cd <- ptpv4_coordinate_dirs(root)
  meta <- primary$metadata
  meta$pseudotime <- coordinate[meta$cell_id]
  finite <- is.finite(meta$pseudotime)
  meta <- meta[finite, , drop = FALSE]
  counts <- primary$counts[, meta$cell_id, drop = FALSE]
  attempted <- list()
  out <- tryCatch({
    attempted[[1L]] <- data.frame(
      step = "choose_primary_grid",
      status = "attempted",
      stringsAsFactors = FALSE
    )
    pb <- ptp_choose_primary_grid(meta, counts, cfg)
    if (!isTRUE(pb$estimable)) {
      stop("coordinate_common_support_not_estimable", call. = FALSE)
    }
    fit <- ptp_fit_cellmeans(
      pb$counts,
      pb$metadata,
      pb$eligible_subbins,
      "ptp_group",
      cfg$qc
    )
    local <- list(
      metadata = meta,
      counts = counts,
      config = cfg,
      primary_pb = pb,
      primary_fit = fit,
      programs = primary$programs,
      role_table = primary$role_table
    )
    fake <- dirs
    fake$primary$score_models <- cd$score_models
    fake$primary$fry_mroast <- cd$fry_mroast
    fake$primary$mouse_gene <- cd$mouse_gene
    fake$primary$adaptive <- cd$adaptive
    fake$primary$family <- cd$family
    fake$primary$multiscale <- cd$multiscale
    fake$primary$trajectory <- cd$trajectory
    score <- ptpv4_run_score_models(local, cfg, args, fake)
    fry <- ptpv4_run_fry_mroast(
      fit,
      primary$programs,
      primary$role_table,
      cfg,
      args,
      fake,
      output_dir = cd$fry_mroast,
      coordinate_mode = TRUE
    )
    mouse <- ptpv4_build_mouse_standardized_expression(
      fit,
      pb$eligible_subbins,
      fake
    )
    union_genes <- sort(unique(unlist(
      lapply(membership$memberships, `[[`, "genes"),
      use.names = FALSE
    )))
    symbols <- ptp_clean_gene_symbols(rownames(mouse$expression))
    idx <- match(union_genes, symbols)
    keep <- is.finite(idx)
    idx <- idx[keep]
    union_genes <- union_genes[keep]
    E <- mouse$expression[idx, , drop = FALSE]
    rownames(E) <- union_genes
    T <- ptpv4_interaction_stats(E, assignments)$t
    rownames(T) <- union_genes
    adaptive <- ptpv4_run_adaptive_from_matrices(
      T,
      E,
      assignments,
      membership,
      cfg,
      args,
      cd$adaptive,
      prefix = coordinate_id
    )
    family <- ptpv4_run_family_hierarchy(
      adaptive,
      membership,
      mouse,
      cfg,
      args,
      dirs,
      output_dir = cd$family
    )
    multiscale <- ptpv4_run_multiscale(
      meta,
      counts,
      mouse$metadata,
      assignments,
      membership,
      cfg,
      args,
      dirs,
      output_dir = cd$multiscale,
      prefix = coordinate_id
    )
    trajectory <- ptpv4_run_trajectory_global(
      meta,
      counts,
      mouse$metadata,
      assignments,
      membership,
      cfg,
      args,
      dirs,
      output_dir = cd$trajectory,
      prefix = coordinate_id
    )
    attempted[[2L]] <- data.frame(
      step = "all_required_coordinate_methods",
      status = "ok",
      stringsAsFactors = FALSE
    )
    adaptive$gene_matrices <- NULL
    list(
      status = data.frame(
        coordinate_id = coordinate_id,
        status = "ok",
        reason = "",
        eligible_subbins = paste(pb$eligible_subbins, collapse = ";"),
        n_mouse_subbin_observations = nrow(fit$metadata),
        stringsAsFactors = FALSE
      ),
      score = score,
      fry_mroast = fry,
      adaptive = adaptive,
      family = family,
      multiscale = multiscale,
      trajectory = trajectory
    )
  }, error = function(e) {
    attempted[[length(attempted) + 1L]] <- data.frame(
      step = "coordinate_branch",
      status = "failed",
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
    list(
      status = data.frame(
        coordinate_id = coordinate_id,
        status = "not_estimable",
        reason = conditionMessage(e),
        eligible_subbins = "",
        n_mouse_subbin_observations = NA_integer_,
        stringsAsFactors = FALSE
      )
    )
  })
  ptp_write_csv(
    out$status,
    file.path(root, "coordinate_branch_status.csv")
  )
  ptp_write_csv(
    ptpv4_bind_rows(attempted),
    file.path(root, "coordinate_branch_attempt_audit.csv")
  )
  out
}

ptpv4_coordinate_status_row <- function(diagnostics, branch_status) {
  row <- diagnostics$summary
  row$coordinate_learning_status <- row$status
  row$analysis_status <- branch_status$status[[1L]]
  row$analysis_reason <- branch_status$reason[[1L]]
  row$status <- row$analysis_status
  row$reason <- row$analysis_reason
  row$eligible_subbins <- branch_status$eligible_subbins[[1L]]
  row$n_mouse_subbin_observations <-
    branch_status$n_mouse_subbin_observations[[1L]]
  row
}

ptpv4_run_coordinates <- function(
  primary,
  membership,
  assignments,
  cfg,
  args,
  dirs,
  primary_results
) {
  universe <- sort(unique(unlist(
    lapply(membership$memberships, `[[`, "genes"),
    use.names = FALSE
  )))
  ptp_write_csv(
    data.frame(
      gene_symbol = universe,
      freeze_status = "frozen_before_coordinate_results",
      exclusion_list_checksum = ptpv4_object_sha(universe),
      stringsAsFactors = FALSE
    ),
    file.path(
      dirs$coordinates,
      "target_genes_excluded_frozen_gene_list.csv"
    )
  )
  meta <- primary$metadata
  counts <- primary$counts
  original <- setNames(meta$pseudotime, meta$cell_id)
  original_diag <- ptpv4_coordinate_diagnostics(
    "original_coordinate",
    original,
    meta,
    meta$cell_id,
    character(0L),
    meta$cell_id[as.character(meta$cluster) == "6"],
    "supplied_frozen_pseudotime",
    ptpv4_object_sha(original),
    list(
      used_for_training_subset = FALSE,
      used_for_orientation = FALSE
    )
  )
  original_root <- file.path(dirs$coordinates, "original_coordinate")
  ptp_ensure_dir(original_root)
  original_summary <- original_diag$summary
  original_summary$coordinate_learning_status <- "supplied_frozen_coordinate"
  original_summary$analysis_status <- "ok"
  original_summary$analysis_reason <- ""
  original_summary$status <- "ok"
  original_summary$reason <- ""
  original_summary$eligible_subbins <- paste(
    primary$primary_pb$eligible_subbins,
    collapse = ";"
  )
  original_summary$n_mouse_subbin_observations <-
    nrow(primary$primary_fit$metadata)
  ptp_write_csv(
    original_diag$cells,
    file.path(original_root, "coordinate_cell_values.csv")
  )
  ptp_write_csv(
    original_summary,
    file.path(original_root, "coordinate_summary.csv")
  )
  ptp_write_csv(
    data.frame(
      method = c(
        "five_score_models",
        "adaptive_exact_gene_set",
        "family_hierarchical",
        "multiscale_scan",
        "trajectory_global",
        "fry_mroast"
      ),
      status = "completed_in_primary_initial_ploidy_identical_coordinate",
      source_path = c(
        dirs$primary$score_models,
        dirs$primary$adaptive,
        dirs$primary$family,
        dirs$primary$multiscale,
        dirs$primary$trajectory,
        dirs$primary$fry_mroast
      ),
      source_checksum = vapply(c(
        dirs$primary$score_models,
        dirs$primary$adaptive,
        dirs$primary$family,
        dirs$primary$multiscale,
        dirs$primary$trajectory,
        dirs$primary$fry_mroast
      ), ptpv4_object_sha, character(1L)),
      stringsAsFactors = FALSE
    ),
    file.path(original_root, "original_coordinate_method_map.csv")
  )
  learning <- list()
  diagnostics <- list(original_coordinate = original_diag)
  status_rows <- list(original_summary)
  control_ids <- meta$cell_id[meta$treatment == "control"]
  untreated <- tryCatch({
    vars <- ptpv4_select_variable_genes(
      counts[, control_ids, drop = FALSE],
      n_features = as.integer(
        cfg$coordinates$untreated_learned_coordinate$variable_genes
      )
    )
    learned <- ptpv4_learn_monocle3_training(
      counts,
      meta[match(control_ids, meta$cell_id), , drop = FALSE],
      vars,
      root_cluster = as.character(
        cfg$coordinates$untreated_learned_coordinate$root_cluster
      ),
      dimensions = as.integer(
        cfg$coordinates$untreated_learned_coordinate$pca_dimensions
      )
    )
    projected <- ptpv4_project_training_pseudotime(
      counts,
      training_ids = control_ids,
      all_ids = meta$cell_id,
      genes = learned$variable_genes,
      training_pseudotime = learned$pseudotime,
      dimensions = as.integer(
        cfg$coordinates$untreated_learned_coordinate$pca_dimensions
      ),
      k = as.integer(
        cfg$coordinates$untreated_learned_coordinate$projection_k
      )
    )
    checksum <- ptpv4_object_sha(list(
      genes = learned$variable_genes,
      pca_center = projected$pca$center,
      pca_rotation = projected$pca$rotation,
      training_pseudotime = learned$pseudotime,
      k = cfg$coordinates$untreated_learned_coordinate$projection_k
    ))
    diag <- ptpv4_coordinate_diagnostics(
      "untreated_learned_coordinate",
      projected$pseudotime,
      meta,
      control_ids,
      character(0L),
      learned$root_cells,
      "control_monocle3_training_plus_training_PCA_kNN_projection",
      checksum,
      list(
        used_for_training_subset = TRUE,
        used_for_orientation = FALSE
      )
    )
    list(
      status = "ok",
      coordinate = projected$pseudotime,
      diagnostics = diag,
      variable_genes = vars,
      root_cells = learned$root_cells,
      projection_checksum = checksum
    )
  }, error = function(e) {
    list(status = "not_estimable", reason = conditionMessage(e))
  })
  if (identical(untreated$status, "ok")) {
    root <- file.path(dirs$coordinates, "untreated_learned_coordinate")
    ptp_ensure_dir(root)
    ptp_write_csv(
      untreated$diagnostics$cells,
      file.path(root, "coordinate_cell_values.csv")
    )
    ptp_write_csv(
      untreated$diagnostics$summary,
      file.path(root, "coordinate_summary.csv")
    )
    ptp_write_csv(
      untreated$diagnostics$ecdf,
      file.path(root, "coordinate_ecdf.csv")
    )
    ptp_write_csv(
      untreated$diagnostics$density,
      file.path(root, "coordinate_density.csv")
    )
    ptp_write_csv(
      untreated$variable_genes,
      file.path(root, "training_variable_genes.csv")
    )
    branch <- ptpv4_run_coordinate_branch(
      "untreated_learned_coordinate",
      untreated$coordinate,
      primary,
      membership,
      assignments,
      cfg,
      args,
      dirs
    )
    learning$untreated_learned_coordinate <- branch
    diagnostics$untreated_learned_coordinate <- untreated$diagnostics
    summary <- ptpv4_coordinate_status_row(
      untreated$diagnostics,
      branch$status
    )
    ptp_write_csv(
      summary,
      file.path(root, "coordinate_summary.csv")
    )
    status_rows[[length(status_rows) + 1L]] <- summary
  } else {
    row <- data.frame(
      coordinate_id = "untreated_learned_coordinate",
      status = "not_estimable",
      reason = untreated$reason,
      stringsAsFactors = FALSE
    )
    ptp_write_csv(
      row,
      file.path(
        ptp_ensure_dir(file.path(
          dirs$coordinates,
          "untreated_learned_coordinate"
        )),
        "coordinate_summary.csv"
      )
    )
    status_rows[[length(status_rows) + 1L]] <- row
  }
  excluded <- tryCatch({
    vars <- ptpv4_select_variable_genes(
      counts,
      n_features = as.integer(
        cfg$coordinates$target_genes_excluded_coordinate$variable_genes
      ),
      excluded_symbols = universe
    )
    learned <- ptpv4_learn_monocle3_training(
      counts,
      meta,
      vars,
      root_cluster = as.character(
        cfg$coordinates$target_genes_excluded_coordinate$root_cluster
      ),
      dimensions = as.integer(
        cfg$coordinates$target_genes_excluded_coordinate$pca_dimensions
      )
    )
    checksum <- ptpv4_object_sha(list(
      excluded_genes = universe,
      variable_genes = learned$variable_genes,
      raw_pseudotime = learned$raw_pseudotime,
      root_cells = learned$root_cells
    ))
    diag <- ptpv4_coordinate_diagnostics(
      "target_genes_excluded_coordinate",
      learned$pseudotime,
      meta,
      meta$cell_id,
      universe,
      learned$root_cells,
      "all_target_cells_monocle3_after_frozen_tested_gene_exclusion",
      checksum,
      list(
        used_for_training_subset = FALSE,
        used_for_orientation = FALSE
      )
    )
    list(
      status = "ok",
      coordinate = learned$pseudotime,
      diagnostics = diag,
      variable_genes = vars,
      root_cells = learned$root_cells,
      projection_checksum = checksum
    )
  }, error = function(e) {
    list(status = "not_estimable", reason = conditionMessage(e))
  })
  if (identical(excluded$status, "ok")) {
    root <- file.path(dirs$coordinates, "target_genes_excluded_coordinate")
    ptp_ensure_dir(root)
    ptp_write_csv(
      excluded$diagnostics$cells,
      file.path(root, "coordinate_cell_values.csv")
    )
    ptp_write_csv(
      excluded$diagnostics$summary,
      file.path(root, "coordinate_summary.csv")
    )
    ptp_write_csv(
      excluded$diagnostics$ecdf,
      file.path(root, "coordinate_ecdf.csv")
    )
    ptp_write_csv(
      excluded$diagnostics$density,
      file.path(root, "coordinate_density.csv")
    )
    ptp_write_csv(
      excluded$variable_genes,
      file.path(root, "training_variable_genes.csv")
    )
    branch <- ptpv4_run_coordinate_branch(
      "target_genes_excluded_coordinate",
      excluded$coordinate,
      primary,
      membership,
      assignments,
      cfg,
      args,
      dirs
    )
    learning$target_genes_excluded_coordinate <- branch
    diagnostics$target_genes_excluded_coordinate <- excluded$diagnostics
    summary <- ptpv4_coordinate_status_row(
      excluded$diagnostics,
      branch$status
    )
    ptp_write_csv(
      summary,
      file.path(root, "coordinate_summary.csv")
    )
    status_rows[[length(status_rows) + 1L]] <- summary
  } else {
    row <- data.frame(
      coordinate_id = "target_genes_excluded_coordinate",
      status = "not_estimable",
      reason = excluded$reason,
      stringsAsFactors = FALSE
    )
    ptp_write_csv(
      row,
      file.path(
        ptp_ensure_dir(file.path(
          dirs$coordinates,
          "target_genes_excluded_coordinate"
        )),
        "coordinate_summary.csv"
      )
    )
    status_rows[[length(status_rows) + 1L]] <- row
  }
  status <- ptpv4_bind_rows(status_rows)
  ptp_write_csv(
    status,
    file.path(dirs$coordinates, "coordinate_sensitivity_status.csv")
  )
  list(
    status = status,
    diagnostics = diagnostics,
    branches = learning
  )
}

ptpv4_effective_rank_table <- function(adaptive, membership) {
  X <- adaptive$gene_matrices$mouse_expression
  rows <- lapply(names(membership$memberships), function(mid) {
    genes <- membership$memberships[[mid]]$genes
    idx <- match(genes, rownames(X))
    idx <- idx[is.finite(idx)]
    if (length(idx) < 2L) {
      er <- NA_real_
    } else {
      Z <- t(scale(t(X[idx, , drop = FALSE])))
      Z[!is.finite(Z)] <- 0
      singular <- tryCatch(
        svd(Z, nu = 0L, nv = 0L)$d^2,
        error = function(e) numeric(0L)
      )
      er <- if (length(singular) && sum(singular^2) > 0) {
        sum(singular)^2 / sum(singular^2)
      } else {
        NA_real_
      }
    }
    data.frame(
      membership_id = mid,
      n_genes = length(idx),
      effective_rank = er,
      stringsAsFactors = FALSE
    )
  })
  ptpv4_bind_rows(rows)
}

ptpv4_result_classification <- function(results, membership) {
  rows <- list()
  add_program <- function(df, method, p_col, fdr_cols) {
    if (is.null(df) || !nrow(df) || !(p_col %in% names(df))) {
      return()
    }
    x <- df
    x$method <- method
    x$p_value_normalized <- x[[p_col]]
    x$fdr_normalized <- NA_real_
    for (fdr_col in fdr_cols) {
      if (fdr_col %in% names(x)) {
        idx <- is.na(x$fdr_normalized) & is.finite(x[[fdr_col]])
        x$fdr_normalized[idx] <- x[[fdr_col]][idx]
      }
    }
    keep <- intersect(c(
      "method", "program_id", "membership_id", "program_label",
      "tier", "response_family", "p_value_normalized", "fdr_normalized"
    ), names(x))
    rows[[length(rows) + 1L]] <<- x[, keep, drop = FALSE]
  }
  primary_score <- results$score$results[
    results$score$results$contrast_id ==
      "treatment_by_initial_ploidy_interaction",
    ,
    drop = FALSE
  ]
  add_program(
    primary_score,
    "five_score_models",
    "p_value",
    c(
      "fdr_focused_tier1_tier2",
      "fdr_tier3_exploratory",
      "fdr_negative_control",
      "fdr_all_universe"
    )
  )
  for (nm in c("fry", "mroast_mean", "mroast_msq")) {
    add_program(
      results$fry_mroast[[nm]],
      nm,
      "selected_p_value",
      c(
        "fdr_focused_tier1_tier2",
        "fdr_tier3_exploratory",
        "fdr_negative_control",
        "fdr_all_universe"
      )
    )
  }
  add_program(
    results$camera$results,
    "camera",
    "PValue",
    c(
      "fdr_focused_tier1_tier2",
      "fdr_tier3_exploratory",
      "fdr_controls",
      "fdr_all_universe"
    )
  )
  add_program(
    results$control_reference$results,
    "control_reference",
    "p_value",
    c(
      "fdr_focused_tier1_tier2",
      "fdr_tier3_exploratory",
      "fdr_negative_control",
      "fdr_all_universe"
    )
  )
  add_program(
    results$covariance_whitened$results,
    "covariance_whitened",
    "p_value",
    c(
      "fdr_focused_tier1_tier2",
      "fdr_tier3_exploratory",
      "fdr_negative_control",
      "fdr_all_universe"
    )
  )
  adaptive <- results$adaptive$membership_results
  adaptive <- merge(
    adaptive,
    membership$alias[!duplicated(membership$alias$membership_id), c(
      "membership_id", "program_id", "program_label", "tier",
      "response_family"
    )],
    by = "membership_id",
    all.x = TRUE,
    sort = FALSE
  )
  add_program(
    adaptive,
    "nested_exact_adaptive",
    "adaptive_empirical_p",
    c(
      "fdr_focused_tier1_tier2",
      "fdr_tier3",
      "fdr_controls",
      "fdr_all_memberships"
    )
  )
  out <- ptpv4_bind_rows(rows)
  if (!nrow(out)) {
    return(out)
  }
  out$classification <- ifelse(
    is.finite(out$fdr_normalized) & out$fdr_normalized < 0.10,
    "passes_relevant_FDR",
    ifelse(
      is.finite(out$p_value_normalized) & out$p_value_normalized < 0.05,
      "nominal_only",
      "not_significant"
    )
  )
  out
}

ptpv4_evidence_matrix <- function(classification) {
  if (is.null(classification) || !nrow(classification)) {
    return(data.frame())
  }
  split_rows <- split(
    classification,
    interaction(classification$method, classification$tier, drop = TRUE)
  )
  ptpv4_bind_rows(lapply(split_rows, function(x) {
    data.frame(
      method = x$method[[1L]],
      tier = x$tier[[1L]],
      n_tests = nrow(x),
      n_nominal_0_05 = sum(x$p_value_normalized < 0.05, na.rm = TRUE),
      n_relevant_fdr_0_10 = sum(x$fdr_normalized < 0.10, na.rm = TRUE),
      minimum_p = min(x$p_value_normalized, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

ptpv4_method_calibration_table <- function(
  simulation,
  adaptive,
  assignments,
  camera,
  membership
) {
  sim <- simulation$method_status
  if (nrow(sim)) {
    sim$calibration_basis <- "reduced-null_whole-mouse_simulation"
    sim$scientific_use <- ifelse(
      sim$status == "calibrated",
      "eligible_subject_to_other_qualification",
      "cannot_support_positive_conclusion"
    )
  }
  control_mids <- unique(membership$alias$membership_id[
    membership$alias$tier == "negative_control"
  ])
  control <- adaptive$membership_results[
    adaptive$membership_results$membership_id %in% control_mids,
    ,
    drop = FALSE
  ]
  controls_ok <- nrow(control) > 0L &&
    !any(control$fdr_controls < 0.10, na.rm = TRUE) &&
    sum(control$adaptive_empirical_p < 0.05, na.rm = TRUE) <= 1L
  assignments_ok <- nrow(assignments$table) == 4900L &&
    length(unique(assignments$table$assignment_checksum)) == 4900L &&
    sum(assignments$table$is_observed_assignment) == 1L
  exact_status <- if (assignments_ok && controls_ok) {
    "calibrated_exact_randomization"
  } else {
    "uncalibrated_negative_control_or_assignment_failure"
  }
  exact <- data.frame(
    method_id = c(
      "nested_exact_adaptive_gene_set",
      "response_family_nested_omnibus",
      "multiscale_window_scan_maxT",
      "trajectory_global_exact_permutation",
      "dose_nested_exact_permutation"
    ),
    status = exact_status,
    calibration_basis =
      "complete_stratified_randomization_space_and_negative_controls",
    scientific_use = if (grepl("^calibrated", exact_status)) {
      "eligible_subject_to_other_qualification"
    } else {
      "cannot_support_positive_conclusion"
    },
    stringsAsFactors = FALSE
  )
  rotation <- data.frame(
    method_id = c(
      "fry", "mroast_mean", "mroast_msq",
      "control_reference_residual_score",
      "covariance_whitened_score",
      "treatment_blind_multidimensional_pca"
    ),
    status = "uncalibrated_no_dedicated_null_simulation",
    calibration_basis = "not_covered_by_five-score_simulation",
    scientific_use = "exploratory_or_supporting_only",
    stringsAsFactors = FALSE
  )
  camera_row <- data.frame(
    method_id = "camera",
    status = camera$calibration$status[[1L]],
    calibration_basis = "negative_control_camera_grid",
    scientific_use = if (camera$calibration$status[[1L]] == "diagnostic_only") {
      "diagnostic_only"
    } else {
      "exploratory_or_supporting_only"
    },
    stringsAsFactors = FALSE
  )
  ptpv4_bind_rows(list(sim, exact, rotation, camera_row))
}

ptpv4_scientific_qualification <- function(
  family,
  adaptive,
  influence,
  calibration,
  membership,
  coordinates
) {
  alias <- family$alias_claims
  exact_calibrated <- any(
    calibration$method_id == "nested_exact_adaptive_gene_set" &
      grepl("^calibrated", calibration$status)
  )
  loo <- influence$loo_mouse
  loo_summary <- if (nrow(loo)) {
    aggregate(
      sign_concordant ~ membership_id,
      loo,
      function(x) mean(x, na.rm = TRUE)
    )
  } else {
    data.frame(membership_id = character(0L), sign_concordant = numeric(0L))
  }
  dominant <- influence$dominant
  dominant_summary <- if (nrow(dominant)) {
    aggregate(
      absolute_contribution ~ membership_id,
      dominant,
      max,
      na.rm = TRUE
    )
  } else {
    data.frame(
      membership_id = character(0L),
      absolute_contribution = numeric(0L)
    )
  }
  out <- merge(alias, loo_summary, by = "membership_id", all.x = TRUE)
  out <- merge(out, dominant_summary, by = "membership_id", all.x = TRUE)
  out$exact_method_calibrated <- exact_calibrated
  out$mouse_stability_pass <- is.finite(out$sign_concordant) &
    out$sign_concordant >= 0.75
  out$gene_stability_pass <- is.finite(out$absolute_contribution) &
    out$absolute_contribution < 0.50
  out$coordinate_status <- if (
    all(coordinates$status$status %in% c("ok", "not_estimable"))
  ) {
    "no_silent_coordinate_failure"
  } else {
    "coordinate_audit_failure"
  }
  out$confirmatory_qualified <- out$focused_claim_gate &
    out$exact_method_calibrated &
    out$mouse_stability_pass &
    out$gene_stability_pass &
    out$coordinate_status == "no_silent_coordinate_failure"
  out$qualification <- ifelse(
    out$confirmatory_qualified,
    "focused_confirmatory",
    ifelse(
      out$focused_claim_gate,
      "statistical_gate_only_not_fully_qualified",
      ifelse(
        is.finite(out$adaptive_empirical_p) &
          out$adaptive_empirical_p < 0.05,
        "suggestive_or_exploratory",
        "not_significant"
      )
    )
  )
  out
}

ptpv4_tier_summary <- function(family_results, qualification) {
  rows <- lapply(
    list(
      Tier1 = "tier1_core",
      Tier2 = "tier2_upgraded_focused",
      Tier3 = "tier3_exploratory_only",
      controls = "negative_control"
    ),
    function(pattern) {
      fam <- family_results[grepl(pattern, family_results$tiers), , drop = FALSE]
      q <- qualification[qualification$tier == pattern, , drop = FALSE]
      relevant_fdr <- if (pattern %in% c(
        "tier1_core",
        "tier2_upgraded_focused"
      )) {
        fam$fdr_focused_11
      } else if (pattern == "tier3_exploratory_only") {
        fam$fdr_tier3_3
      } else {
        fam$fdr_controls
      }
      data.frame(
        tier = pattern,
        n_families = nrow(fam),
        n_family_fdr_0_10 = sum(relevant_fdr < 0.10, na.rm = TRUE),
        n_confirmatory_qualified_program_labels =
          sum(q$confirmatory_qualified, na.rm = TRUE),
        minimum_family_p = if (nrow(fam)) {
          min(fam$family_omnibus_p, na.rm = TRUE)
        } else {
          NA_real_
        },
        conclusion = if (sum(q$confirmatory_qualified, na.rm = TRUE) > 0L) {
          "confirmatory evidence present"
        } else if (sum(relevant_fdr < 0.10, na.rm = TRUE) > 0L) {
          "family statistical signal did not pass all scientific qualification gates"
        } else {
          "no family passes its pre-specified FDR family"
        },
        stringsAsFactors = FALSE
      )
    }
  )
  out <- ptpv4_bind_rows(rows)
  out$display_tier <- names(rows)
  out
}

ptpv4_safe_plot <- function(path, title, builder) {
  plot <- tryCatch(builder(), error = function(e) NULL)
  if (is.null(plot)) {
    plot <- ggplot2::ggplot(
      data.frame(x = 0, y = 0, label = "No estimable rows for this figure."),
      ggplot2::aes(x, y)
    ) +
      ggplot2::geom_text(ggplot2::aes(label = label), size = 4) +
      ggplot2::xlim(-1, 1) +
      ggplot2::ylim(-1, 1) +
      ggplot2::theme_void() +
      ggplot2::labs(title = title)
  }
  ggplot2::ggsave(
    path,
    plot,
    width = 9,
    height = 6,
    dpi = 160,
    bg = "white"
  )
  invisible(path)
}

ptpv4_make_figures <- function(
  results,
  primary,
  membership,
  calibration,
  classification,
  evidence,
  dirs
) {
  catalog <- list()
  add <- function(section, title, stub, legend, interpretation, builder) {
    n <- length(catalog) + 1L
    filename <- sprintf("fig%02d_%s.png", n, stub)
    path <- file.path(dirs$figures, filename)
    ptpv4_safe_plot(path, paste0("Figure ", n, ". ", title), builder)
    catalog[[n]] <<- data.frame(
      figure_number = n,
      section = section,
      title = title,
      filename = filename,
      absolute_path = path,
      legend = legend,
      interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  }
  add(
    3,
    "Common-support and mouse coverage on the frozen primary grid",
    "common_support_mouse_coverage",
    "Bars show retained target-cell counts for each mouse and eligible primary pseudotime subbin; facets separate initial-ploidy strata and fill indicates treatment.",
    "The plot documents the support actually entering the matched-state model. Missing or low-support mouse-subbins are not interpolated.",
    function() {
      x <- primary$primary_fit$metadata
      ggplot2::ggplot(
        x,
        ggplot2::aes(
          x = subbin_id,
          y = cell_count,
          fill = treatment,
          group = sample_id
        )
      ) +
        ggplot2::geom_col(position = "dodge") +
        ggplot2::facet_wrap(~ initial_ploidy) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
        ) +
        ggplot2::labs(x = "Frozen common-support subbin", y = "Target cells")
    }
  )
  add(
    2,
    "Exact membership size and effective expression rank",
    "program_size_effective_rank",
    "Each point is one exact unique gene membership. Effective rank is calculated from its 16-mouse standardized expression covariance; color denotes the represented tier set.",
    "Large gene sets can have much smaller effective rank because genes are correlated. Exact aliases are represented once.",
    function() {
      x <- ptpv4_effective_rank_table(results$adaptive, membership)
      x <- merge(
        x,
        membership$membership_table[, c("membership_id", "tiers")],
        by = "membership_id"
      )
      ggplot2::ggplot(
        x,
        ggplot2::aes(n_genes, effective_rank, color = tiers)
      ) +
        ggplot2::geom_point(alpha = 0.75, size = 2) +
        ggplot2::scale_x_log10() +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::labs(
          x = "Genes in exact membership (log scale)",
          y = "Effective rank",
          color = "Tier set"
        )
    }
  )
  add(
    9,
    "Five-score treatment effects and the primary interaction forest",
    "score_interaction_forest",
    "Points and horizontal intervals are estimates and 95% confidence intervals from the signed weighted mean z-score model. Panels show the Gem-Control effect in 2N-origin mice, the effect in 4N-origin mice, and their difference.",
    "This is a true effect-size forest plot. Intervals crossing zero do not provide two-sided evidence for that contrast.",
    function() {
      x <- results$score$results[
        results$score$results$score_method_id ==
          "signed_weighted_mean_zscore",
        ,
        drop = FALSE
      ]
      reps <- primary$role_table[
        order(primary$role_table$response_family, primary$role_table$program_id),
        ,
        drop = FALSE
      ]
      reps <- reps[!duplicated(reps$response_family), ]
      x <- x[x$program_id %in% reps$program_id, , drop = FALSE]
      x$program_label <- factor(
        x$program_label,
        levels = rev(unique(x$program_label))
      )
      ggplot2::ggplot(
        x,
        ggplot2::aes(
          x = estimate,
          y = program_label,
          xmin = ci_low,
          xmax = ci_high,
          color = tier_label
        )
      ) +
        ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
        ggplot2::geom_errorbar(orientation = "y", width = 0.15) +
        ggplot2::geom_point(size = 1.8) +
        ggplot2::facet_wrap(~ contrast_id, scales = "free_x") +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::labs(
          x = "Standardized score effect",
          y = NULL,
          color = "Tier"
        )
    }
  )
  add(
    4,
    "Reduced-null simulation calibration and power",
    "simulation_calibration",
    "Panels show empirical type-I error, 95% CI coverage, and power at alpha 0.05 for each frozen score definition and simulation scenario. Dashed reference lines mark 0.05 for type-I error and 0.95 for coverage.",
    "Methods outside the frozen calibration bands are marked uncalibrated rather than being repaired after seeing results.",
    function() {
      x <- results$simulation$summary
      metric <- rbind(
        data.frame(x, metric = "Type-I error", value = x$type_I_error),
        data.frame(x, metric = "95% CI coverage", value = x$ci_coverage),
        data.frame(x, metric = "Power", value = x$power_0_05)
      )
      metric <- metric[is.finite(metric$value), ]
      ggplot2::ggplot(
        metric,
        ggplot2::aes(
          x = scenario_id,
          y = value,
          color = score_method_id,
          group = score_method_id
        )
      ) +
        ggplot2::geom_hline(
          data = data.frame(
            metric = c("Type-I error", "95% CI coverage"),
            y = c(0.05, 0.95)
          ),
          ggplot2::aes(yintercept = y),
          linetype = 2
        ) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::facet_wrap(~ metric, scales = "free_x") +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
        ) +
        ggplot2::labs(x = "Frozen scenario", y = "Probability", color = "Score")
    }
  )
  add(
    10,
    "Adaptive component tests and nested omnibus calibration",
    "adaptive_components_omnibus",
    "Each point is an empirical P value for a component statistic or the permutation-calibrated adaptive omnibus, summarized across exact unique memberships.",
    "The omnibus P values include statistic selection inside the nested permutation; uncalibrated minimum component P values are not used as final inference.",
    function() {
      c0 <- results$adaptive$component_results[, c(
        "membership_id", "component", "empirical_p"
      )]
      o <- results$adaptive$membership_results[, c(
        "membership_id", "adaptive_empirical_p"
      )]
      names(o)[2L] <- "empirical_p"
      o$component <- "nested_adaptive_omnibus"
      x <- rbind(c0, o[, names(c0)])
      ggplot2::ggplot(
        x,
        ggplot2::aes(
          x = component,
          y = -log10(pmax(empirical_p, 1 / 4901))
        )
      ) +
        ggplot2::geom_boxplot(outlier.shape = NA, fill = "#dbeafe") +
        ggplot2::geom_jitter(width = 0.18, alpha = 0.25, size = 0.7) +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        ) +
        ggplot2::labs(x = "Component", y = "-log10 empirical P")
    }
  )
  add(
    11,
    "Response-family hierarchical omnibus evidence",
    "family_hierarchical_heatmap",
    "Rows are the 15 frozen response families. Color shows the permutation-calibrated family omnibus evidence; labels report the relevant family-level FDR.",
    "Tier 1+2, Tier 3, and controls use separate FDR families; the transparent all-family FDR is also retained in the tables.",
    function() {
      x <- results$family$family_results
      x$relevant_fdr <- ifelse(
        grepl("tier1_core|tier2_upgraded_focused", x$tiers),
        x$fdr_focused_11,
        ifelse(
          grepl("tier3_exploratory_only", x$tiers),
          x$fdr_tier3_3,
          x$fdr_controls
        )
      )
      x$label <- paste0(x$response_family_label, "\nFDR=", signif(x$relevant_fdr, 3))
      ggplot2::ggplot(
        x,
        ggplot2::aes(
          x = 1,
          y = reorder(response_family_label, -family_omnibus_p),
          fill = -log10(pmax(family_omnibus_p, 1 / 4901))
        )
      ) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = signif(relevant_fdr, 3))) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        ) +
        ggplot2::labs(
          x = NULL,
          y = NULL,
          fill = "-log10 family P"
        )
    }
  )
  add(
    8,
    "Negative-control calibration across exact and camera methods",
    "negative_control_calibration",
    "Points show negative-control P values from the exact adaptive test and camera variants. The horizontal line marks nominal P=0.05.",
    "A systematic excess of small control P values disqualifies the corresponding method; camera is fixed as diagnostic-only when its control rule fails.",
    function() {
      control_mids <- unique(membership$alias$membership_id[
        membership$alias$tier == "negative_control"
      ])
      a <- results$adaptive$membership_results[
        results$adaptive$membership_results$membership_id %in% control_mids,
        ,
        drop = FALSE
      ]
      a$method <- "nested_exact_adaptive"
      a$p <- a$adaptive_empirical_p
      c0 <- results$camera$results[
        results$camera$results$tier == "negative_control",
        ,
        drop = FALSE
      ]
      c0$method <- c0$method_id
      c0$p <- c0$PValue
      x <- rbind(
        a[, c("method", "p")],
        c0[, c("method", "p")]
      )
      ggplot2::ggplot(x, ggplot2::aes(method, p, color = method)) +
        ggplot2::geom_hline(yintercept = 0.05, linetype = 2) +
        ggplot2::geom_jitter(width = 0.15, height = 0, size = 2) +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
          legend.position = "none"
        ) +
        ggplot2::labs(x = "Method", y = "Negative-control P value")
    }
  )
  add(
    12,
    "Permutation-adjusted multiscale pseudotime window scan",
    "multiscale_window_scan",
    "Curves show observed studentized interactions across all support-eligible contiguous windows for response-family gene unions; color denotes window width.",
    "Best windows localize effects only. Inference uses the across-window maximum statistic and its exact permutation-adjusted P value.",
    function() {
      x <- results$multiscale$observed_windows
      x <- x[x$entity_type == "response_family_union", , drop = FALSE]
      x$response_family <- sub("^family__", "", x$entity_id)
      highlight <- unique(membership$alias$response_family[
        grepl(
          "nucleotide|repair|interferon|rna|gemcitabine",
          membership$alias$response_family,
          ignore.case = TRUE
        )
      ])
      x <- x[x$response_family %in% head(highlight, 6L), , drop = FALSE]
      x$midpoint <- (x$start + x$end) / 2
      ggplot2::ggplot(
        x,
        ggplot2::aes(
          midpoint,
          observed_t,
          color = factor(width_bins),
          group = interaction(entity_id, width_bins)
        )
      ) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2) +
        ggplot2::geom_line(alpha = 0.75) +
        ggplot2::facet_wrap(~ response_family, scales = "free_y") +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::labs(
          x = "Window midpoint",
          y = "Observed interaction t",
          color = "Width (bins)"
        )
    }
  )
  add(
    13,
    "Trajectory-wide interaction curves, simultaneous bands, and support",
    "trajectory_interaction_curve_support",
    "Lines are response-family interaction estimates over support-only merged bins; ribbons are 95% permutation maxT simultaneous bands and point size shows minimum group mouse support.",
    "Unsupported regions are absent rather than interpolated. Global P values select spline df inside every whole-mouse permutation.",
    function() {
      x <- results$trajectory$curve
      x <- x[x$entity_type == "response_family_union", , drop = FALSE]
      x$response_family <- sub("^family__", "", x$entity_id)
      highlight <- unique(membership$alias$response_family[
        grepl(
          "nucleotide|repair|interferon|rna|gemcitabine",
          membership$alias$response_family,
          ignore.case = TRUE
        )
      ])
      x <- x[x$response_family %in% head(highlight, 6L), , drop = FALSE]
      support <- aggregate(
        n_mice_qc ~ region_id,
        results$trajectory$support,
        min
      )
      x <- merge(x, support, by = "region_id", all.x = TRUE)
      ggplot2::ggplot(
        x,
        ggplot2::aes(midpoint, observed_estimate, group = response_family)
      ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(
            ymin = simultaneous_95_low,
            ymax = simultaneous_95_high
          ),
          fill = "#bfdbfe",
          alpha = 0.5
        ) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2) +
        ggplot2::geom_line(color = "#1d4ed8") +
        ggplot2::geom_point(ggplot2::aes(size = n_mice_qc)) +
        ggplot2::facet_wrap(~ response_family, scales = "free_y") +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::labs(
          x = "Pseudotime",
          y = "Interaction estimate",
          size = "Minimum mice/group"
        )
    }
  )
  add(
    14,
    "Marginal total-effect versus matched-state interaction",
    "marginal_vs_matched",
    "Each point is one frozen program under the signed weighted mean score. The x-axis is the matched-state interaction and the y-axis is the marginal total-expression interaction; bars show marginal 95% confidence intervals.",
    "Concordance indicates similar marginal and state-conditioned evidence. Discordance can reflect state occupancy or conditioning and is not formal causal mediation.",
    function() {
      a <- results$score$results[
        results$score$results$score_method_id ==
          "signed_weighted_mean_zscore" &
          results$score$results$contrast_id ==
          "treatment_by_initial_ploidy_interaction",
        c("program_id", "estimate", "tier_label"),
        drop = FALSE
      ]
      names(a)[2L] <- "matched_estimate"
      b <- results$marginal$score_results[
        results$marginal$score_results$score_method_id ==
          "signed_weighted_mean_zscore" &
          results$marginal$score_results$contrast_id ==
          "treatment_by_initial_ploidy_interaction",
        c("program_id", "estimate", "ci_low", "ci_high"),
        drop = FALSE
      ]
      names(b)[2:4] <- c("marginal_estimate", "marginal_low", "marginal_high")
      x <- merge(a, b, by = "program_id")
      ggplot2::ggplot(
        x,
        ggplot2::aes(
          matched_estimate,
          marginal_estimate,
          color = tier_label
        )
      ) +
        ggplot2::geom_hline(yintercept = 0, linetype = 3) +
        ggplot2::geom_vline(xintercept = 0, linetype = 3) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = marginal_low, ymax = marginal_high),
          width = 0
        ) +
        ggplot2::geom_point() +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::labs(
          x = "Matched-state interaction",
          y = "Marginal interaction",
          color = "Tier"
        )
    }
  )
  add(
    15,
    "Dose-adaptive interaction evidence",
    "dose_interaction",
    "Points show the nested exact dose-adaptive empirical P value for each exact membership, grouped by its first response-family alias. The line marks nominal P=0.05.",
    "The six candidate dose statistics are selected inside the complete 176,400-assignment permutation calibration; dose results remain secondary.",
    function() {
      x <- merge(
        results$dose$results,
        membership$alias[!duplicated(membership$alias$membership_id), c(
          "membership_id", "response_family", "tier"
        )],
        by = "membership_id",
        all.x = TRUE
      )
      ggplot2::ggplot(
        x,
        ggplot2::aes(
          reorder(response_family, adaptive_empirical_p),
          adaptive_empirical_p,
          color = tier
        )
      ) +
        ggplot2::geom_hline(yintercept = 0.05, linetype = 2) +
        ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.7) +
        ggplot2::coord_flip() +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::labs(
          x = "Response family",
          y = "Dose-adaptive exact P",
          color = "Tier"
        )
    }
  )
  add(
    16,
    "Alternative-coordinate concordance with original pseudotime",
    "coordinate_concordance",
    "Cells are plotted by original pseudotime versus each estimable alternative coordinate. Color indicates treatment and panels identify the coordinate.",
    "Coordinates were learned with their frozen training/exclusion rules; treatment response was not used to orient or select a favorable direction.",
    function() {
      x <- ptpv4_bind_rows(lapply(
        results$coordinates$diagnostics,
        `[[`,
        "cells"
      ))
      x <- x[x$coordinate_id != "original_coordinate", , drop = FALSE]
      if (nrow(x) > 1500L) {
        set.seed(74167)
        x <- x[sample(seq_len(nrow(x)), 1500L), , drop = FALSE]
      }
      ggplot2::ggplot(
        x,
        ggplot2::aes(
          original_pseudotime,
          coordinate,
          color = treatment
        )
      ) +
        ggplot2::geom_point(alpha = 0.35, size = 0.8) +
        ggplot2::facet_wrap(~ coordinate_id) +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::labs(
          x = "Original pseudotime",
          y = "Alternative coordinate",
          color = "Treatment"
        )
    }
  )
  add(
    18,
    "Leave-one-mouse interaction influence",
    "leave_one_mouse_influence",
    "For each mouse, points summarize the maximum standardized change across exact memberships after leaving that mouse out; color denotes initial ploidy and treatment.",
    "Large values identify influential mice. Confirmatory qualification requires direction stability across the leave-one-mouse analysis.",
    function() {
      x <- results$influence$loo_mouse
      x <- merge(
        x,
        primary$primary_fit$metadata[!duplicated(
          primary$primary_fit$metadata$sample_id
        ), c("sample_id", "initial_ploidy", "treatment")],
        by.x = "left_out_mouse",
        by.y = "sample_id",
        all.x = TRUE
      )
      x <- aggregate(
        standardized_change ~ left_out_mouse + initial_ploidy + treatment,
        x,
        max,
        na.rm = TRUE
      )
      ggplot2::ggplot(
        x,
        ggplot2::aes(
          reorder(left_out_mouse, standardized_change),
          standardized_change,
          color = interaction(initial_ploidy, treatment)
        )
      ) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::coord_flip() +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::labs(
          x = "Left-out mouse",
          y = "Maximum standardized interaction change",
          color = "Group"
        )
    }
  )
  add(
    18,
    "Dominant-gene and leave-one-gene influence",
    "dominant_gene_leave_one_gene",
    "Points are the genes with the largest absolute contribution to their membership's observed exact statistic; x-position is contribution and color is the leave-one-gene change in true maxmean.",
    "A single dominant gene is reported as sparse or gene-driven evidence and is not described as broad pathway activation.",
    function() {
      d <- results$influence$dominant
      g <- results$influence$gene_influence
      x <- merge(
        d[d$dominant_top, ],
        g[, c("membership_id", "gene_symbol", "absolute_change")],
        by = c("membership_id", "gene_symbol")
      )
      x <- x[order(x$absolute_contribution, decreasing = TRUE), ]
      x <- head(x, 40L)
      ggplot2::ggplot(
        x,
        ggplot2::aes(
          absolute_contribution,
          reorder(gene_symbol, absolute_contribution),
          color = absolute_change
        )
      ) +
        ggplot2::geom_point(size = 2.3) +
        ggplot2::scale_color_viridis_c() +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::labs(
          x = "Absolute contribution",
          y = "Gene",
          color = "Leave-one-gene change"
        )
    }
  )
  add(
    19,
    "Method-by-tier evidence matrix",
    "method_by_tier_evidence_matrix",
    "Tiles show the number of results passing the relevant FDR<0.10 within each method and tier. Printed labels are counts.",
    "Counts summarize concordance but do not turn exact aliases or correlated methods into independent evidence.",
    function() {
      ggplot2::ggplot(
        evidence,
        ggplot2::aes(method, tier, fill = n_relevant_fdr_0_10)
      ) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(
          ggplot2::aes(label = n_relevant_fdr_0_10),
          size = 3
        ) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
        ) +
        ggplot2::labs(
          x = "Method",
          y = "Tier",
          fill = "FDR<0.10"
        )
    }
  )
  out <- ptpv4_bind_rows(catalog)
  ptp_write_csv(out, file.path(dirs$figures, "figure_catalog.csv"))
  out
}

ptpv4_html_table <- function(df, n = 30L) {
  ptpv3_html_table(df, n = n)
}

ptpv4_figure_html <- function(catalog, section) {
  x <- catalog[catalog$section == section, , drop = FALSE]
  if (!nrow(x)) {
    return("")
  }
  paste(vapply(seq_len(nrow(x)), function(i) {
    r <- x[i, , drop = FALSE]
    src <- paste0(
      "data:image/png;base64,",
      base64enc::base64encode(r$absolute_path, linewidth = 0L)
    )
    paste0(
      "<figure id='figure-", r$figure_number, "'>",
      "<h3>Figure ", r$figure_number, ". ",
      ptp_html_escape(r$title), "</h3>",
      "<img src='", src, "' alt='Figure ", r$figure_number, ". ",
      ptp_html_escape(r$title), "'/>",
      "<figcaption><b>Legend.</b> ", ptp_html_escape(r$legend),
      "<br/><b>Result interpretation.</b> ",
      ptp_html_escape(r$interpretation), "</figcaption></figure>"
    )
  }, character(1L)), collapse = "\n")
}

ptpv4_write_html_report <- function(
  results,
  primary,
  membership,
  calibration,
  qualification,
  tier_summary,
  acceptance,
  catalog,
  dirs
) {
  family <- results$family$family_results
  tier1 <- family[grepl("tier1_core", family$tiers), , drop = FALSE]
  tier2 <- family[grepl("tier2_upgraded_focused", family$tiers), , drop = FALSE]
  tier3 <- family[grepl("tier3_exploratory_only", family$tiers), , drop = FALSE]
  controls <- family[grepl("negative_control", family$tiers), , drop = FALSE]
  confirmed <- qualification[qualification$confirmatory_qualified, , drop = FALSE]
  conclusion <- if (nrow(confirmed)) {
    paste0(
      nrow(confirmed),
      " program-label rows satisfy the complete focused confirmatory qualification. ",
      "Exact aliases remain one membership-level evidence unit."
    )
  } else {
    paste0(
      "No focused program satisfies the full confirmatory qualification. ",
      "Failure to reject does not establish ploidy independence; the 16-mouse ",
      "design can miss localized, sparse, or unstable interactions."
    )
  }
  sections <- c(
    paste0(
      "<section><h2>1. Scientific question and estimands</h2>",
      "<p>The primary estimand is ",
      ptp_html_escape(cfg_text <- "(Gem-Control)_4N-origin - (Gem-Control)_2N-origin"),
      " within frozen matched pseudotime support, using initial ploidy as the primary ploidy variable.</p>",
      "</section>"
    ),
    paste0(
      "<section><h2>2. Frozen definitions and post-hoc status</h2>",
      "<p>V4 is a post-hoc methodological repair inspired by V3 review. ",
      "The method, memberships, windows, weights, seeds, assignment spaces, ",
      "coordinates and FDR families were frozen before V4 result interpretation. ",
      "The external branch is disabled_by_user_no_external_cohort.</p>",
      ptpv4_html_table(membership$membership_table, 15L),
      ptpv4_figure_html(catalog, 2), "</section>"
    ),
    paste0(
      "<section><h2>3. Input/cohort/common-support QC</h2>",
      ptpv4_html_table(primary$primary_fit$metadata, 20L),
      ptpv4_figure_html(catalog, 3), "</section>"
    ),
    paste0(
      "<section><h2>4. Simulation and method calibration</h2>",
      ptpv4_html_table(results$simulation$summary, 50L),
      ptpv4_html_table(calibration, 50L),
      ptpv4_figure_html(catalog, 4), "</section>"
    ),
    paste0(
      "<section><h2>5. Tier 1 family results</h2>",
      ptpv4_html_table(tier1, 30L), "</section>"
    ),
    paste0(
      "<section><h2>6. Tier 2 family results</h2>",
      ptpv4_html_table(tier2, 30L), "</section>"
    ),
    paste0(
      "<section><h2>7. Tier 3 exploratory family results</h2>",
      ptpv4_html_table(tier3, 30L), "</section>"
    ),
    paste0(
      "<section><h2>8. Negative controls</h2>",
      ptpv4_html_table(controls, 20L),
      ptpv4_html_table(results$camera$calibration, 10L),
      ptpv4_figure_html(catalog, 8), "</section>"
    ),
    paste0(
      "<section><h2>9. Score-model comparison</h2>",
      ptpv4_html_table(
        results$score$results[
          results$score$results$contrast_id ==
            "treatment_by_initial_ploidy_interaction",
          ,
          drop = FALSE
        ],
        50L
      ),
      ptpv4_figure_html(catalog, 9), "</section>"
    ),
    paste0(
      "<section><h2>10. Adaptive gene-set results</h2>",
      ptpv4_html_table(results$adaptive$membership_results, 50L),
      ptpv4_figure_html(catalog, 10), "</section>"
    ),
    paste0(
      "<section><h2>11. Family hierarchical results</h2>",
      ptpv4_html_table(results$family$alias_claims, 50L),
      ptpv4_figure_html(catalog, 11), "</section>"
    ),
    paste0(
      "<section><h2>12. Multiscale pseudotime scan</h2>",
      ptpv4_html_table(results$multiscale$scan_results, 50L),
      ptpv4_figure_html(catalog, 12), "</section>"
    ),
    paste0(
      "<section><h2>13. Trajectory-wide interaction</h2>",
      ptpv4_html_table(results$trajectory$global, 50L),
      ptpv4_html_table(results$trajectory$integrated, 30L),
      ptpv4_figure_html(catalog, 13), "</section>"
    ),
    paste0(
      "<section><h2>14. Marginal versus matched-state comparison</h2>",
      "<p>These branches distinguish marginal total-expression evidence from ",
      "state-conditioned evidence; they are not formal causal mediation.</p>",
      ptpv4_html_table(results$marginal$score_results, 40L),
      ptpv4_figure_html(catalog, 14), "</section>"
    ),
    paste0(
      "<section><h2>15. Dose/ETP/end-time analyses</h2>",
      ptpv4_html_table(results$dose$results, 40L),
      ptpv4_html_table(results$endpoint$endpoint$continuous_program, 20L),
      ptpv4_html_table(results$endpoint$endpoint$group_program, 20L),
      ptpv4_figure_html(catalog, 15), "</section>"
    ),
    paste0(
      "<section><h2>16. Coordinate sensitivity</h2>",
      ptpv4_html_table(results$coordinates$status, 20L),
      ptpv4_figure_html(catalog, 16), "</section>"
    ),
    paste0(
      "<section><h2>17. 04i cross-fit bridge</h2>",
      "<p>04i is an internal cross-fitted exploratory bridge, not independent ",
      "validation or a positive control.</p>",
      ptpv4_html_table(results$crossfit$status, 10L),
      ptpv4_html_table(results$crossfit$fold_audit, 20L), "</section>"
    ),
    paste0(
      "<section><h2>18. Leave-one-mouse/gene influence</h2>",
      ptpv4_html_table(qualification, 50L),
      ptpv4_figure_html(catalog, 18), "</section>"
    ),
    paste0(
      "<section><h2>19. Cross-method concordance</h2>",
      ptpv4_html_table(results$classification, 60L),
      ptpv4_figure_html(catalog, 19), "</section>"
    ),
    paste0(
      "<section><h2>20. Integrated conclusion</h2><p>",
      ptp_html_escape(conclusion), "</p>",
      ptpv4_html_table(tier_summary, 10L),
      if (nrow(confirmed)) ptpv4_html_table(confirmed, 30L) else "",
      "</section>"
    ),
    paste0(
      "<section><h2>21. Limitations and final audit</h2>",
      "<p>Nominal P&lt;0.05, camera-only signals, unadjusted best-window results, ",
      "single-score positives, Tier-3-only signals, sensitivity-only signals, ",
      "and methods failing calibration are exploratory or suggestive only. ",
      "Failure to reject is not evidence of ploidy independence. Sparse or ",
      "gene-driven effects are labeled as such rather than broad activation.</p>",
      ptpv4_html_table(acceptance, 100L), "</section>"
    )
  )
  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'/>",
    "<title>04j V4 pseudotime treatment-by-ploidy report</title>",
    "<style>",
    "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;",
    "margin:28px;line-height:1.5;color:#1f2937}",
    "table{border-collapse:collapse;font-size:11px;margin:12px 0;",
    "max-width:100%;display:block;overflow:auto}",
    "td,th{border:1px solid #d8dee9;padding:4px 6px}",
    "th{background:#edf2f7;position:sticky;top:0}",
    "img{max-width:100%;border:1px solid #d8dee9}",
    "figure{margin:24px 0}figcaption{font-size:13px;margin-top:8px}",
    "h1,h2,h3{color:#111827}section{margin:28px 0}",
    "</style></head><body>",
    "<h1>04j V4 matched-state treatment-by-initial-ploidy analysis</h1>",
    paste(sections, collapse = "\n"),
    "</body></html>"
  )
  path <- file.path(
    dirs$report,
    "04j_pseudotime_treatment_ploidy_programs_v4_report.html"
  )
  writeLines(html, path, useBytes = TRUE)
  path
}

ptpv4_finalize_checksums <- function(cfg, dirs) {
  before_path <- file.path(dirs$manifest, "source_checksums_before.csv")
  before <- read.csv(before_path, check.names = FALSE, stringsAsFactors = FALSE)
  after <- before
  after$phase <- "after"
  after$bytes <- suppressWarnings(file.info(after$path)$size)
  after$sha256 <- vapply(after$path, ptpv4_file_sha, character(1L))
  ptp_write_csv(
    after,
    file.path(dirs$manifest, "source_checksums_after.csv")
  )
  source_comp <- ptpv4_compare_checksum_tables(before, after)
  source_comp$expected <- ifelse(
    source_comp$scope == "plan_document",
    "changed_by_v4_amendment",
    "unchanged"
  )
  source_comp$acceptance_status <- ifelse(
    source_comp$scope == "plan_document",
    ifelse(source_comp$status == "changed", "ok", "failed"),
    ifelse(source_comp$status == "unchanged", "ok", "failed")
  )
  ptp_write_csv(
    source_comp,
    file.path(dirs$manifest, "source_checksum_comparison.csv")
  )
  result_comp <- list()
  roots <- list(
    v1_result_root = cfg$analysis$v1_result_root,
    v2_result_root = cfg$analysis$v2_result_root,
    v3_result_root = cfg$analysis$v3_result_root
  )
  for (scope in names(roots)) {
    before_file <- file.path(
      dirs$manifest,
      sub("_root$", "", paste0(scope, "_checksums_before.csv"))
    )
    before_file <- file.path(
      dirs$manifest,
      paste0(sub("_root$", "", scope), "_checksums_before.csv")
    )
    if (!file.exists(before_file)) {
      before_file <- file.path(
        dirs$manifest,
        paste0(sub("_result_root$", "_result", scope), "_checksums_before.csv")
      )
    }
    b <- read.csv(before_file, check.names = FALSE, stringsAsFactors = FALSE)
    a <- ptpv4_checksum_tree(roots[[scope]], scope, "after")
    out_name <- paste0(sub("_root$", "", scope), "_checksums_after.csv")
    ptp_write_csv(a, file.path(dirs$manifest, out_name))
    comp <- ptpv4_compare_checksum_tables(b, a)
    result_comp[[scope]] <- comp
    ptp_write_csv(
      comp,
      file.path(
        dirs$manifest,
        paste0(sub("_root$", "", scope), "_checksum_comparison.csv")
      )
    )
  }
  comparison <- ptpv4_bind_rows(result_comp)
  ptp_write_csv(
    comparison,
    file.path(dirs$manifest, "v1_v2_v3_result_checksum_comparison.csv")
  )
  list(source = source_comp, results = comparison)
}

ptpv4_build_acceptance <- function(
  results,
  primary,
  membership,
  assignments,
  inputs,
  checksums,
  calibration,
  catalog,
  report_path,
  dirs
) {
  report_text <- if (file.exists(report_path)) {
    paste(readLines(report_path, warn = FALSE), collapse = "\n")
  } else {
    ""
  }
  input_required <- inputs$input %in% c(
    "cell_metadata", "noncell_metadata", "sample_info", "seurat_rds"
  )
  input_ok <- inputs$v3_comparison_status[input_required] ==
    "matched_v3_manifest"
  input_detail <- paste(
    inputs$input[input_required],
    inputs$v3_comparison_status[input_required],
    sep = "=",
    collapse = ";"
  )
  role_complete <- nrow(primary$role_table) == 106L &&
    all(nzchar(primary$role_table$tier)) &&
    all(nzchar(primary$role_table$response_family)) &&
    all(nzchar(primary$role_table$analysis_role))
  alias_complete <- nrow(membership$alias) == 106L &&
    all(membership$alias$membership_id %in%
      membership$membership_table$membership_id)
  dose_audit <- results$dose$audit
  mixed_ok <- all(
    results$fry_mroast$fry$selected_p_source[
      !results$fry_mroast$fry$directional_p_legal
    ] == "PValue.Mixed"
  ) &&
    all(
      results$fry_mroast$mroast_msq$selected_p_source == "PValue.Mixed"
    )
  crossfit_ok <- nrow(results$crossfit$fold_audit) == 16L &&
    all(results$crossfit$fold_audit$status == "ok") &&
    all(!results$crossfit$fold_audit$training_contains_heldout) &&
    all(results$crossfit$fold_audit$normalization_training_only) &&
    all(!results$crossfit$fold_audit$reused_full_data_04j_gene_z)
  pca_ok <- nrow(results$multidimensional$audit) > 0L &&
    all(!results$multidimensional$audit$treatment_labels_used_for_basis) &&
    all(results$multidimensional$audit$basis_invariant_across_assignments)
  no_program_one <- !any(
    primary$role_table$program_id == "1" |
      membership$alias$program_id == "1"
  )
  simulation_numeric <- nrow(results$simulation$summary) > 0L &&
    all(c(
      "type_I_error", "empirical_FDR_q_0_10", "ci_coverage",
      "mcse_reject", "score_specific_true_effect_mean"
    ) %in% names(results$simulation$summary))
  figures_embedded <- nrow(catalog) >= 15L &&
    lengths(regmatches(
      report_text,
      gregexpr("data:image/png;base64", report_text, fixed = TRUE)
    ))[[1L]] >= nrow(catalog)
  rows <- list(
    data.frame(
      check = "v1_v2_v3_code_unchanged",
      status = if (
        all(checksums$source$acceptance_status == "ok")
      ) "ok" else "failed",
      detail = paste(table(checksums$source$status), collapse = ";"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "v1_v2_v3_result_roots_unchanged",
      status = if (
        all(checksums$results$status == "unchanged")
      ) "ok" else "failed",
      detail = paste(table(checksums$results$status), collapse = ";"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "input_checksums_match_v3_manifest",
      status = if (all(input_ok)) "ok" else "failed",
      detail = input_detail,
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "program_universe_matches_v3",
      status = if (
        all(results$program_comparison$status == "unchanged")
      ) "ok" else "failed",
      detail = paste(table(results$program_comparison$status), collapse = ";"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "106_labels_have_tier_family_role",
      status = if (role_complete) "ok" else "failed",
      detail = paste(nrow(primary$role_table), "labels"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "unique_membership_alias_mapping_complete",
      status = if (alias_complete) "ok" else "failed",
      detail = paste(
        nrow(membership$alias),
        "labels;",
        nrow(membership$membership_table),
        "unique memberships"
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "primary_4900_assignments_unique_observed",
      status = if (
        nrow(assignments$table) == 4900L &&
          length(unique(assignments$table$assignment_checksum)) == 4900L &&
          sum(assignments$table$is_observed_assignment) == 1L
      ) "ok" else "failed",
      detail = paste(
        nrow(assignments$table),
        "assignments;",
        sum(assignments$table$is_observed_assignment),
        "observed"
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "dose_assignment_space",
      status = if (
        dose_audit$generated_assignments == 176400L &&
          dose_audit$unique_assignment_checksums == 176400L &&
          dose_audit$observed_assignment_count == 1L
      ) "ok" else "failed",
      detail = paste0(
        dose_audit$generated_assignments,
        " generated; inference=",
        dose_audit$inference_assignment_count
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "fry_mroast_mixed_p_fdr",
      status = if (mixed_ok) "ok" else "failed",
      detail = "broad/mixed fry and mroast use PValue.Mixed; msq always Mixed",
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "crossfit_heldout_no_training_leakage",
      status = if (crossfit_ok) "ok" else "failed",
      detail = paste(
        nrow(results$crossfit$fold_audit),
        "folds;",
        sum(results$crossfit$fold_audit$status != "ok"),
        "failed"
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "PCA_permutation_no_label_leakage",
      status = if (pca_ok) "ok" else "failed",
      detail = "basis trained on all 16 mice without treatment labels",
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "trajectory_support_only",
      status = if (
        nrow(results$trajectory$support) > 0L &&
          all(results$trajectory$support$eligible_support_only[
            results$trajectory$support$region_id %in%
              unique(results$trajectory$curve$region_id)
          ])
      ) "ok" else "failed",
      detail = paste(length(unique(results$trajectory$curve$region_id)), "bins"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "eligible_methods_result_or_not_estimable",
      status = if (
        nrow(results$coordinates$status) == 3L &&
          all(results$coordinates$status$status %in% c("ok", "not_estimable"))
      ) "ok" else "failed",
      detail = paste(
        results$coordinates$status$coordinate_id,
        results$coordinates$status$status,
        collapse = ";"
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "no_program_id_1_or_missing_metadata",
      status = if (no_program_one && role_complete) "ok" else "failed",
      detail = "checked frozen labels and membership aliases",
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "simulation_numeric_calibration",
      status = if (simulation_numeric) "ok" else "failed",
      detail = paste(
        nrow(results$simulation$summary),
        "scenario-method rows"
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "negative_controls_in_calibration",
      status = if (
        any(calibration$calibration_basis ==
          "complete_stratified_randomization_space_and_negative_controls")
      ) "ok" else "failed",
      detail = paste(
        sum(membership$alias$tier == "negative_control"),
        "control labels"
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "report_figures_all_embedded",
      status = if (figures_embedded) "ok" else "failed",
      detail = paste(nrow(catalog), "catalogued PNG figures"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "figure_metadata_complete",
      status = if (
        all(nzchar(catalog$title)) &&
          all(nzchar(catalog$legend)) &&
          all(nzchar(catalog$interpretation))
      ) "ok" else "failed",
      detail = "continuous figure number, title, legend, interpretation",
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "no_proxy_forest",
      status = if (
        !any(grepl("forest_proxy", list.files(
          dirs$root,
          recursive = TRUE
        ), fixed = TRUE)) &&
          !grepl("forest_proxy", report_text, fixed = TRUE)
      ) "ok" else "failed",
      detail = "true estimate-and-CI forest only",
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "no_external_cohort",
      status = if (
        identical(
          primary$config$analysis$external_signature_status,
          "disabled_by_user_no_external_cohort"
        )
      ) "ok" else "failed",
      detail = primary$config$analysis$external_signature_status,
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "no_species_filter_or_split_added",
      status = if (
        grepl(
          "no_species_filter_or_split",
          primary$config$analysis$species_policy,
          fixed = TRUE
        )
      ) "ok" else "failed",
      detail = primary$config$analysis$species_policy,
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "git_diff_check",
      status = "pending_external_command",
      detail = "updated after git diff --check",
      stringsAsFactors = FALSE
    )
  )
  ptpv4_bind_rows(rows)
}

ptpv4_mde_table <- function(score_results, current_total_mice = 16L) {
  x <- score_results[
    score_results$score_method_id == "signed_weighted_mean_zscore" &
      score_results$contrast_id ==
        "treatment_by_initial_ploidy_interaction",
    ,
    drop = FALSE
  ]
  if (!nrow(x)) {
    return(data.frame())
  }
  multiplier <- stats::qnorm(0.975) + stats::qnorm(0.80)
  x$mde_80_power_alpha_0_05 <- multiplier * x$se
  target <- pmax(abs(x$estimate), 0.25 * x$mde_80_power_alpha_0_05)
  x$estimated_total_mice_for_observed_effect <- ceiling(
    current_total_mice *
      (x$mde_80_power_alpha_0_05 / target)^2
  )
  x$estimated_total_mice_for_observed_effect <-
    pmax(x$estimated_total_mice_for_observed_effect, current_total_mice)
  x$estimated_additional_mice <- pmax(
    0,
    x$estimated_total_mice_for_observed_effect - current_total_mice
  )
  x[, c(
    "program_id", "program_label", "tier", "response_family",
    "estimate", "se", "mde_80_power_alpha_0_05",
    "estimated_total_mice_for_observed_effect",
    "estimated_additional_mice"
  ), drop = FALSE]
}

ptpv4_matched_marginal_comparison <- function(primary_score, marginal_score) {
  a <- primary_score[
    primary_score$contrast_id ==
      "treatment_by_initial_ploidy_interaction",
    c(
      "program_id", "score_method_id", "estimate", "se", "ci_low",
      "ci_high", "p_value", "tier", "response_family"
    ),
    drop = FALSE
  ]
  b <- marginal_score[
    marginal_score$contrast_id ==
      "treatment_by_initial_ploidy_interaction",
    c(
      "program_id", "score_method_id", "estimate", "se", "ci_low",
      "ci_high", "p_value"
    ),
    drop = FALSE
  ]
  names(a)[3:7] <- paste0("matched_", names(a)[3:7])
  names(b)[3:7] <- paste0("marginal_", names(b)[3:7])
  out <- merge(
    a,
    b,
    by = c("program_id", "score_method_id"),
    all = TRUE,
    sort = FALSE
  )
  out$direction_concordant <- sign(out$matched_estimate) ==
    sign(out$marginal_estimate)
  out
}

ptpv4_copy_or_create_baselines <- function(cfg, dirs, repo_root) {
  required <- c(
    "source_checksums_before.csv",
    "v1_result_checksums_before.csv",
    "v2_result_checksums_before.csv",
    "v3_result_checksums_before.csv"
  )
  missing <- required[!file.exists(file.path(dirs$manifest, required))]
  canonical <- file.path(cfg$analysis$v4_result_root, "00_manifest")
  if (length(missing) && dir.exists(canonical)) {
    for (nm in missing) {
      src <- file.path(canonical, nm)
      if (file.exists(src)) {
        file.copy(src, file.path(dirs$manifest, nm), overwrite = FALSE)
      }
    }
  }
  if (!file.exists(file.path(dirs$manifest, "source_checksums_before.csv"))) {
    legacy <- file.path(repo_root, c(
      "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs.R",
      "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_util.R",
      "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_config.yaml",
      "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_v2.R",
      "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_v2_util.R",
      "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_v2_config.yaml",
      "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_v3.R",
      "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_v3_util.R",
      "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_v3_config.yaml",
      "Data/in-vivo/Plan/pseudotime_treatment_by_ploidy_matched_state_plan.md"
    ))
    scope <- ifelse(
      grepl("/Plan/", legacy, fixed = TRUE),
      "plan_document",
      "source_code"
    )
    x <- data.frame(
      scope = scope,
      phase = "before",
      path = legacy,
      relative_path = sub(paste0("^", repo_root, "/"), "", legacy),
      bytes = file.info(legacy)$size,
      sha256 = vapply(legacy, ptpv4_file_sha, character(1L)),
      stringsAsFactors = FALSE
    )
    ptp_write_csv(
      x,
      file.path(dirs$manifest, "source_checksums_before.csv")
    )
  }
  roots <- list(
    v1_result = cfg$analysis$v1_result_root,
    v2_result = cfg$analysis$v2_result_root,
    v3_result = cfg$analysis$v3_result_root
  )
  for (nm in names(roots)) {
    path <- file.path(dirs$manifest, paste0(nm, "_checksums_before.csv"))
    if (!file.exists(path)) {
      x <- ptpv4_checksum_tree(
        roots[[nm]],
        paste0(nm, "_root"),
        "before"
      )
      ptp_write_csv(x, path)
    }
  }
  invisible(required)
}

ptpv4_whitespace_diff_check <- function(repo_root, files) {
  tracked <- system2(
    "git",
    c("-C", repo_root, "diff", "--check"),
    stdout = TRUE,
    stderr = TRUE
  )
  tracked_status <- attr(tracked, "status") %||% 0L
  source_rows <- list()
  for (path in files) {
    lines <- readLines(path, warn = FALSE)
    bad <- grep("[[:blank:]]+$", lines)
    source_rows[[length(source_rows) + 1L]] <- data.frame(
      path = path,
      trailing_whitespace_lines = paste(bad, collapse = ";"),
      status = if (length(bad)) "failed" else "ok",
      stringsAsFactors = FALSE
    )
  }
  source <- ptpv4_bind_rows(source_rows)
  list(
    status = if (tracked_status == 0L && all(source$status == "ok")) {
      "ok"
    } else {
      "failed"
    },
    detail = paste(
      c(tracked, paste(source$path, source$status, sep = "=")),
      collapse = ";"
    ),
    source = source
  )
}

ptpv4_run_workflow <- function(args, repo_root, script_dir) {
  started <- Sys.time()
  ptpv4_thread_env()
  set.seed(as.integer(args$seed))
  worker_fields <- c(
    "workers", "permutation_workers", "simulation_workers",
    "crossfit_workers", "adaptive_workers", "multiscale_workers",
    "trajectory_workers", "dose_workers", "coordinate_workers",
    "rotation_test_workers", "model_workers", "within_model_workers"
  )
  for (field in worker_fields) {
    args[[field]] <- max(1L, min(8L, as.integer(args[[field]])))
  }
  ptp_check_packages()
  pst_check_packages()
  required <- c(
    "rhdf5", "base64enc", "FNN", "monocle3", "edgeR", "limma",
    "ggplot2", "digest", "yaml", "Matrix"
  )
  missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required installed packages: ", paste(missing, collapse = ", "))
  }
  cfg <- ptpv4_read_config(args$config)
  schema <- ptpv4_validate_config(cfg)
  output_root <- ptpv4_prepare_output(args$output_root, args$overwrite)
  dirs <- ptpv4_dirs(output_root)
  ptpv4_copy_or_create_baselines(cfg, dirs, repo_root)
  preexisting_checkpoints <- list.files(
    dirs$checkpoints,
    pattern = "\\.rds$",
    full.names = TRUE
  )
  log_file <- file.path(dirs$logs, "run.log")
  ptpv4_message(
    "Starting 04j V4 workflow at ",
    output_root,
    "; smoke=",
    isTRUE(args$smoke),
    "; preexisting checkpoints=",
    length(preexisting_checkpoints),
    log_file = log_file
  )
  ptp_write_csv(
    data.frame(
      parameter = names(args),
      value = vapply(args, as.character, character(1L)),
      stringsAsFactors = FALSE
    ),
    file.path(dirs$manifest, "analysis_parameters.csv")
  )
  ptp_write_csv(schema, file.path(dirs$manifest, "config_schema_validation.csv"))
  ptp_write_csv(
    ptp_package_versions(),
    file.path(dirs$manifest, "package_versions.csv")
  )
  file.copy(
    args$config,
    file.path(dirs$frozen, basename(args$config)),
    overwrite = TRUE
  )
  file.copy(
    cfg$analysis$plan_path,
    file.path(dirs$frozen, basename(cfg$analysis$plan_path)),
    overwrite = TRUE
  )
  v4_sources <- file.path(script_dir, c(
    "04j_pseudotime_treatment_ploidy_programs_v4.R",
    "04j_pseudotime_treatment_ploidy_programs_v4_util.R",
    "04j_pseudotime_treatment_ploidy_programs_v4_config.yaml"
  ))
  ptp_write_csv(
    data.frame(
      path = v4_sources,
      bytes = file.info(v4_sources)$size,
      sha256 = vapply(v4_sources, ptpv4_file_sha, character(1L)),
      stringsAsFactors = FALSE
    ),
    file.path(dirs$manifest, "v4_source_manifest.csv")
  )
  seurat_rds <- ptp_resolve_seurat_rds(
    args$seurat_rds,
    cfg,
    args$results_root
  )
  inputs <- ptpv4_record_input_manifest(
    args,
    cfg,
    repo_root,
    dirs,
    seurat_rds
  )
  input_sha <- ptpv4_object_sha(inputs[, c("input", "path", "sha256")])
  cfg_sha <- ptpv4_config_sha(args, cfg)
  previous_contract_path <- file.path(
    dirs$manifest,
    "checkpoint_contract.csv"
  )
  previous_contract <- if (file.exists(previous_contract_path)) {
    tryCatch(
      read.csv(
        previous_contract_path,
        check.names = FALSE,
        stringsAsFactors = FALSE
      ),
      error = function(e) data.frame()
    )
  } else {
    data.frame()
  }
  eligible_preexisting_checkpoint_count <- if (
    nrow(previous_contract) == 1L &&
      identical(previous_contract$config_checksum[[1L]], cfg_sha) &&
      identical(previous_contract$input_checksum[[1L]], input_sha) &&
      identical(
        previous_contract$method_version[[1L]],
        cfg$analysis$method_version
      ) &&
      identical(
        as.character(previous_contract$checkpoint_schema[[1L]]),
        as.character(cfg$checkpoint$schema_version)
      )
  ) {
    length(preexisting_checkpoints)
  } else {
    0L
  }
  ptpv4_message(
    "Schema-eligible preexisting checkpoints=",
    eligible_preexisting_checkpoint_count,
    "/",
    length(preexisting_checkpoints),
    log_file = log_file
  )
  ptp_write_csv(
    data.frame(
      config_checksum = cfg_sha,
      input_checksum = input_sha,
      method_version = cfg$analysis$method_version,
      checkpoint_schema = cfg$checkpoint$schema_version,
      stringsAsFactors = FALSE
    ),
    file.path(dirs$manifest, "checkpoint_contract.csv")
  )
  primary <- ptpv4_load_primary_data(args, cfg, dirs, seurat_rds)
  cfg <- primary$config
  if (file.exists(args$sample_info)) {
    crosswalk <- ptp_processing_batch_crosswalk(
      primary$metadata,
      args$sample_info,
      ptp_endpoint_ploidy_specs(cfg)
    )
    ptp_write_csv(
      crosswalk,
      file.path(dirs$qc, "processing_batch_crosswalk.csv")
    )
    ptp_write_csv(
      ptp_batch_group_separability(crosswalk),
      file.path(dirs$qc, "processing_batch_separability_audit.csv")
    )
  }
  membership <- ptpv4_membership_objects(
    primary$programs,
    primary$role_table
  )
  family_defs <- ptpv4_write_frozen_universe(
    primary$programs,
    primary$role_table,
    membership,
    cfg,
    dirs
  )
  program_comparison <- ptpv4_compare_program_universe_v3(
    primary$programs,
    membership,
    cfg,
    dirs
  )
  ptp_write_csv(
    ptp_subbins_table(cfg$primary_subbins, "primary_four_bin"),
    file.path(dirs$frozen, "primary_subbin_definitions.csv")
  )
  ptp_write_csv(
    ptp_subbins_table(cfg$fallback_subbins, "fallback_two_bin"),
    file.path(dirs$frozen, "fallback_subbin_definitions.csv")
  )
  ptp_write_csv(
    ptp_regions_table(cfg),
    file.path(dirs$frozen, "pseudotime_region_definitions.csv")
  )
  ptp_write_csv(
    data.frame(
      seed_name = names(cfg$reproducibility),
      value = unlist(cfg$reproducibility),
      stringsAsFactors = FALSE
    ),
    file.path(dirs$frozen, "random_seeds_and_threads.csv")
  )
  mouse_meta <- unique(primary$primary_fit$metadata[, c(
    "sample_id", "initial_ploidy", "treatment", "dose_mg"
  )])
  assignments <- ptpv4_build_assignment_space(
    mouse_meta,
    as.integer(cfg$exact_permutation$treated_per_ploidy)
  )
  assignment_audit <- data.frame(
    check = c(
      "assignment_count", "unique_checksums", "observed_assignment_count"
    ),
    observed = c(
      nrow(assignments$table),
      length(unique(assignments$table$assignment_checksum)),
      sum(assignments$table$is_observed_assignment)
    ),
    expected = c(4900L, 4900L, 1L),
    status = ifelse(
      c(
        nrow(assignments$table),
        length(unique(assignments$table$assignment_checksum)),
        sum(assignments$table$is_observed_assignment)
      ) == c(4900L, 4900L, 1L),
      "ok",
      "failed"
    ),
    stringsAsFactors = FALSE
  )
  if (any(assignment_audit$status != "ok")) {
    stop("Primary exact assignment audit failed.", call. = FALSE)
  }
  ptp_write_csv(
    assignments$table,
    file.path(dirs$frozen, "primary_4900_assignments.csv")
  )
  ptp_write_csv(
    assignment_audit,
    file.path(dirs$manifest, "primary_assignment_audit.csv")
  )
  ptpv4_message("Running five primary score models.", log_file = log_file)
  score <- ptpv4_checkpoint_compute(
    dirs,
    "primary_five_score_models",
    cfg_sha,
    cfg,
    input_sha,
    function() ptpv4_run_score_models(primary, cfg, args, dirs)
  )
  ptpv4_message("Running corrected fry/mroast and camera.", log_file = log_file)
  fry_mroast <- ptpv4_checkpoint_compute(
    dirs,
    "primary_fry_mroast",
    cfg_sha,
    cfg,
    input_sha,
    function() ptpv4_run_fry_mroast(
      primary$primary_fit,
      primary$programs,
      primary$role_table,
      cfg,
      args,
      dirs
    )
  )
  camera <- ptpv4_checkpoint_compute(
    dirs,
    "primary_camera_grid",
    cfg_sha,
    cfg,
    input_sha,
    function() ptpv4_run_camera(primary, cfg, dirs)
  )
  mouse_bundle <- ptpv4_checkpoint_compute(
    dirs,
    "primary_mouse_standardized_expression",
    cfg_sha,
    cfg,
    input_sha,
    function() ptpv4_build_mouse_standardized_expression(
      primary$primary_fit,
      primary$primary_pb$eligible_subbins,
      dirs
    )
  )
  if (ncol(mouse_bundle$expression) != 16L) {
    stop("Mouse-standardized expression did not yield 16 mice.", call. = FALSE)
  }
  gene_stats <- ptpv4_compute_gene_statistics(
    mouse_bundle,
    assignments,
    cfg,
    args,
    dirs,
    cfg_sha,
    input_sha
  )
  adaptive <- ptpv4_run_adaptive_memberships(
    gene_stats,
    mouse_bundle,
    membership,
    cfg,
    args,
    dirs,
    cfg_sha,
    input_sha
  )
  family <- ptpv4_checkpoint_compute(
    dirs,
    "primary_response_family_hierarchy",
    cfg_sha,
    cfg,
    input_sha,
    function() ptpv4_run_family_hierarchy(
      adaptive,
      membership,
      mouse_bundle,
      cfg,
      args,
      dirs
    )
  )
  influence <- ptpv4_checkpoint_compute(
    dirs,
    "primary_influence",
    cfg_sha,
    cfg,
    input_sha,
    function() ptpv4_run_influence(
      adaptive,
      family,
      mouse_bundle,
      membership,
      cfg,
      dirs
    )
  )
  control_reference <- ptpv4_checkpoint_compute(
    dirs,
    "primary_control_reference",
    cfg_sha,
    cfg,
    input_sha,
    function() ptpv4_run_control_reference(
      primary,
      mouse_bundle,
      assignments,
      cfg,
      args,
      dirs
    )
  )
  covariance_whitened <- ptpv4_checkpoint_compute(
    dirs,
    "primary_covariance_whitened",
    cfg_sha,
    cfg,
    input_sha,
    function() ptpv4_run_covariance_whitened(
      primary,
      mouse_bundle,
      assignments,
      cfg,
      dirs
    )
  )
  multidimensional <- ptpv4_checkpoint_compute(
    dirs,
    "primary_multidimensional_pca",
    cfg_sha,
    cfg,
    input_sha,
    function() ptpv4_run_multidimensional(
      mouse_bundle,
      membership,
      assignments,
      cfg,
      args,
      dirs
    )
  )
  ptpv4_message("Running exact multiscale scan.", log_file = log_file)
  multiscale <- ptpv4_checkpoint_compute(
    dirs,
    "primary_multiscale_scan",
    cfg_sha,
    cfg,
    input_sha,
    function() ptpv4_run_multiscale(
      primary$metadata,
      primary$counts,
      mouse_bundle$metadata,
      assignments,
      membership,
      cfg,
      args,
      dirs
    )
  )
  ptpv4_message("Running exact trajectory-global branch.", log_file = log_file)
  trajectory <- ptpv4_checkpoint_compute(
    dirs,
    "primary_trajectory_global",
    cfg_sha,
    cfg,
    input_sha,
    function() ptpv4_run_trajectory_global(
      primary$metadata,
      primary$counts,
      mouse_bundle$metadata,
      assignments,
      membership,
      cfg,
      args,
      dirs
    )
  )
  ptpv4_message("Running marginal total-effect branch.", log_file = log_file)
  marginal <- ptpv4_checkpoint_compute(
    dirs,
    "secondary_marginal_total_effect",
    cfg_sha,
    cfg,
    input_sha,
    function() ptpv4_run_marginal(
      primary,
      membership,
      assignments,
      cfg,
      args,
      dirs
    )
  )
  endpoint <- ptpv4_checkpoint_compute(
    dirs,
    "secondary_etp_endtime_standard_dose",
    cfg_sha,
    cfg,
    input_sha,
    function() ptpv4_run_endpoint_and_standard_dose(
      primary,
      cfg,
      args,
      dirs
    )
  )
  ptpv4_message("Running 176400-assignment dose branch.", log_file = log_file)
  dose <- ptpv4_checkpoint_compute(
    dirs,
    "secondary_exact_dose_complete",
    cfg_sha,
    cfg,
    input_sha,
    function() ptpv4_run_exact_dose(
      mouse_bundle,
      membership,
      cfg,
      args,
      dirs,
      cfg_sha,
      input_sha
    )
  )
  ptpv4_message("Running reduced-null simulation.", log_file = log_file)
  simulation <- ptpv4_run_simulation(
    primary,
    mouse_bundle,
    cfg,
    args,
    dirs,
    cfg_sha,
    input_sha
  )
  ptpv4_message("Running corrected cross-fit 04i.", log_file = log_file)
  crossfit <- ptpv4_run_crossfit_04i(
    primary,
    cfg,
    args,
    repo_root,
    dirs,
    cfg_sha,
    input_sha
  )
  core_results <- list(
    score = score,
    fry_mroast = fry_mroast,
    camera = camera,
    adaptive = adaptive,
    family = family,
    influence = influence,
    control_reference = control_reference,
    covariance_whitened = covariance_whitened,
    multidimensional = multidimensional,
    multiscale = multiscale,
    trajectory = trajectory,
    marginal = marginal,
    endpoint = endpoint,
    dose = dose,
    simulation = simulation,
    crossfit = crossfit,
    program_comparison = program_comparison
  )
  ptpv4_message("Running alternative-coordinate branches.", log_file = log_file)
  coordinates <- ptpv4_checkpoint_compute(
    dirs,
    "coordinate_sensitivity_complete",
    cfg_sha,
    cfg,
    input_sha,
    function() ptpv4_run_coordinates(
      primary,
      membership,
      assignments,
      cfg,
      args,
      dirs,
      core_results
    )
  )
  results <- c(core_results, list(coordinates = coordinates))
  classification <- ptpv4_result_classification(results, membership)
  evidence <- ptpv4_evidence_matrix(classification)
  results$classification <- classification
  results$evidence <- evidence
  ptp_write_csv(
    classification,
    file.path(dirs$tables, "result_classification_by_method.csv")
  )
  ptp_write_csv(
    evidence,
    file.path(dirs$tables, "method_by_tier_evidence_matrix.csv")
  )
  calibration <- ptpv4_method_calibration_table(
    simulation,
    adaptive,
    assignments,
    camera,
    membership
  )
  ptp_write_csv(
    calibration,
    file.path(dirs$manifest, "method_calibration_status.csv")
  )
  qualification <- ptpv4_scientific_qualification(
    family,
    adaptive,
    influence,
    calibration,
    membership,
    coordinates
  )
  tier_summary <- ptpv4_tier_summary(
    family$family_results,
    qualification
  )
  mde <- ptpv4_mde_table(score$results)
  matched_marginal <- ptpv4_matched_marginal_comparison(
    score$results,
    marginal$score_results
  )
  ptp_write_csv(
    qualification,
    file.path(dirs$tables, "scientific_qualification_by_program.csv")
  )
  ptp_write_csv(
    tier_summary,
    file.path(dirs$tables, "tier_scientific_summary.csv")
  )
  ptp_write_csv(
    mde,
    file.path(dirs$tables, "minimum_detectable_effect_and_mouse_guidance.csv")
  )
  ptp_write_csv(
    matched_marginal,
    file.path(dirs$tables, "matched_vs_marginal_score_comparison.csv")
  )
  catalog <- ptpv4_make_figures(
    results,
    primary,
    membership,
    calibration,
    classification,
    evidence,
    dirs
  )
  checksums <- ptpv4_finalize_checksums(cfg, dirs)
  placeholder_report <- file.path(
    dirs$report,
    "04j_pseudotime_treatment_ploidy_programs_v4_report.html"
  )
  acceptance <- ptpv4_build_acceptance(
    results,
    primary,
    membership,
    assignments,
    inputs,
    checksums,
    calibration,
    catalog,
    placeholder_report,
    dirs
  )
  report <- ptpv4_write_html_report(
    results,
    primary,
    membership,
    calibration,
    qualification,
    tier_summary,
    acceptance,
    catalog,
    dirs
  )
  acceptance <- ptpv4_build_acceptance(
    results,
    primary,
    membership,
    assignments,
    inputs,
    checksums,
    calibration,
    catalog,
    report,
    dirs
  )
  diff <- ptpv4_whitespace_diff_check(repo_root, c(
    file.path(script_dir, "04j_pseudotime_treatment_ploidy_programs_v4.R"),
    file.path(script_dir, "04j_pseudotime_treatment_ploidy_programs_v4_util.R"),
    file.path(script_dir, "04j_pseudotime_treatment_ploidy_programs_v4_config.yaml"),
    cfg$analysis$plan_path
  ))
  acceptance$status[acceptance$check == "git_diff_check"] <- diff$status
  acceptance$detail[acceptance$check == "git_diff_check"] <- diff$detail
  ptp_write_csv(
    diff$source,
    file.path(dirs$manifest, "source_whitespace_audit.csv")
  )
  ptp_write_csv(
    acceptance,
    file.path(dirs$manifest, "acceptance_checks.csv")
  )
  final_audit <- acceptance
  final_audit$acceptance_class <- ifelse(
    final_audit$status == "ok",
    "passed",
    ifelse(
      final_audit$check %in% c(
        "input_checksums_match_v3_manifest",
        "eligible_methods_result_or_not_estimable"
      ),
      "audited_limitation",
      "failed"
    )
  )
  ptp_write_csv(
    final_audit,
    file.path(dirs$manifest, "final_acceptance_audit.csv")
  )
  report <- ptpv4_write_html_report(
    results,
    primary,
    membership,
    calibration,
    qualification,
    tier_summary,
    final_audit,
    catalog,
    dirs
  )
  finished <- Sys.time()
  timing <- data.frame(
    started_at = format(started, "%Y-%m-%d %H:%M:%S %Z"),
    finished_at = format(finished, "%Y-%m-%d %H:%M:%S %Z"),
    elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
    elapsed_hours = as.numeric(difftime(finished, started, units = "hours")),
    smoke = isTRUE(args$smoke),
    preexisting_checkpoint_count = length(preexisting_checkpoints),
    eligible_preexisting_checkpoint_count =
      eligible_preexisting_checkpoint_count,
    checkpoint_resume_used = eligible_preexisting_checkpoint_count > 0L,
    workers = args$workers,
    stringsAsFactors = FALSE
  )
  ptp_write_csv(timing, file.path(dirs$manifest, "run_timing.csv"))
  writeLines(
    utils::capture.output(sessionInfo()),
    file.path(dirs$manifest, "sessionInfo.txt"),
    useBytes = TRUE
  )
  ptpv4_message(
    "04j V4 workflow complete in ",
    signif(timing$elapsed_hours, 4),
    " hours. Report: ",
    report,
    log_file = log_file
  )
  invisible(output_root)
}
