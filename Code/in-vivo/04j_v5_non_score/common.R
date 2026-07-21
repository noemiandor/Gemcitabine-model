.ptpv5_source_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(e) ""
)
.ptpv5_module_dir <- if (nzchar(.ptpv5_source_file)) {
  dirname(.ptpv5_source_file)
} else {
  file.path(getwd(), "Code/in-vivo/04j_v5_non_score")
}
.ptpv5_script_dir <- dirname(.ptpv5_module_dir)

source(file.path(
  .ptpv5_script_dir,
  "04j_pseudotime_treatment_ploidy_programs_util.R"
))

ptpv5_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

ptpv5_thread_env <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
}

ptpv5_message <- function(..., log_file = NULL) {
  text <- paste0("[", ptpv5_now(), "] ", paste0(..., collapse = ""))
  message(text)
  if (!is.null(log_file)) cat(text, "\n", file = log_file, append = TRUE)
  invisible(text)
}

ptpv5_bind_rows <- function(rows) {
  rows <- rows[!vapply(
    rows,
    function(x) is.null(x) || nrow(as.data.frame(x)) == 0L,
    logical(1L)
  )]
  if (!length(rows)) return(data.frame())
  rows <- lapply(rows, as.data.frame, stringsAsFactors = FALSE)
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    for (nm in setdiff(all_names, names(x))) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, rows)
}

ptpv5_read_config <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The yaml package is required.", call. = FALSE)
  }
  yaml::read_yaml(path)
}

ptpv5_validate_config <- function(cfg) {
  required <- list(
    analysis = c(
      "name", "version", "schema_version", "method_version",
      "primary_estimand", "v1_result_root", "v2_result_root",
      "v3_result_root", "v4_result_root", "v5_result_root"
    ),
    inputs = c(
      "cell_metadata", "seurat_rds", "v4_hdf5",
      "v4_mouse_expression_checkpoint", "v4_gene_statistics_checkpoint",
      "v4_trajectory_checkpoint", "v4_marginal_checkpoint",
      "v4_unique_memberships", "v4_program_alias",
      "v4_program_roles", "v4_response_families",
      "v4_gene_weights", "v4_assignments", "v4_trajectory_support"
    ),
    exact_randomization = c("expected_mice", "expected_assignments"),
    acceptance = c(
      "expected_program_labels", "expected_unique_memberships",
      "expected_response_families", "expected_primary_assignments",
      "expected_mice", "expected_genes", "expected_methods"
    )
  )
  rows <- list()
  for (section in names(required)) {
    for (field in required[[section]]) {
      ok <- !is.null(cfg[[section]]) && !is.null(cfg[[section]][[field]])
      rows[[length(rows) + 1L]] <- data.frame(
        section = section,
        field = field,
        status = if (ok) "ok" else "missing",
        stringsAsFactors = FALSE
      )
    }
  }
  out <- ptpv5_bind_rows(rows)
  if (any(out$status != "ok")) {
    stop(
      "V5 configuration schema is incomplete: ",
      paste0(
        out$section[out$status != "ok"], ".",
        out$field[out$status != "ok"],
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  out
}

ptpv5_dirs <- function(output_root) {
  d <- function(...) ptp_ensure_dir(file.path(output_root, ...))
  list(
    root = ptp_ensure_dir(output_root),
    manifest = d("00_manifest"),
    frozen = d("frozen_definitions"),
    qc = d("qc"),
    ranked = d("01_exact_ranked_enrichment"),
    state = d("02_state_decomposition"),
    bayes = d("03_bayesian_gene_vector"),
    kernel = d("04_exact_global_kernel"),
    trajectory = d("05_gene_wise_trajectory"),
    simulation = d("06_simulation"),
    synthesis = d("07_synthesis"),
    figures = d("figures"),
    tables = d("tables"),
    report = d("report"),
    checkpoints = d("checkpoints"),
    logs = d("logs")
  )
}

ptpv5_prepare_output <- function(path, overwrite = FALSE) {
  ptp_ensure_dir(path)
  existing <- list.files(path, all.files = TRUE, no.. = TRUE)
  allowed_initial <- length(existing) == 0L ||
    identical(sort(existing), "00_manifest")
  checkpoint_dir <- file.path(path, "checkpoints")
  has_checkpoints <- dir.exists(checkpoint_dir) &&
    length(list.files(checkpoint_dir, pattern = "\\.rds$", full.names = TRUE)) > 0L
  if (!isTRUE(overwrite) && !allowed_initial && !has_checkpoints) {
    stop(
      "V5 output root is not empty. Use --overwrite=TRUE to resume or replace V5 outputs: ",
      path,
      call. = FALSE
    )
  }
  if (isTRUE(overwrite)) {
    message(
      "V5 overwrite is non-destructive: baseline manifests and valid checkpoints ",
      "are preserved; method outputs are replaced as they are regenerated."
    )
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

ptpv5_file_sha <- function(path) {
  if (!file.exists(path) || dir.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

ptpv5_object_sha <- function(x) digest::digest(x, algo = "sha256")

ptpv5_config_sha <- function(args, cfg, module_dir) {
  module_files <- sort(list.files(module_dir, pattern = "\\.R$", full.names = TRUE))
  ptpv5_object_sha(list(
    config_sha = ptpv5_file_sha(args$config),
    entry_sha = ptpv5_file_sha(file.path(
      dirname(args$config),
      "04j_pseudotime_treatment_ploidy_programs_v5.R"
    )),
    module_sha = setNames(vapply(
      module_files,
      ptpv5_file_sha,
      character(1L)
    ), basename(module_files)),
    config = cfg,
    smoke = isTRUE(args$smoke),
    method_version = cfg$analysis$method_version
  ))
}

ptpv5_checkpoint_path <- function(dirs, id) {
  file.path(dirs$checkpoints, paste0(id, ".rds"))
}

ptpv5_save_checkpoint <- function(
  dirs,
  id,
  data,
  cfg,
  config_sha,
  input_sha
) {
  payload <- list(
    checkpoint_id = id,
    created_at = ptpv5_now(),
    config_checksum = config_sha,
    method_version = cfg$analysis$method_version,
    schema_version = cfg$checkpoint$schema_version,
    input_checksum = input_sha,
    status = "complete",
    data = data
  )
  saveRDS(payload, ptpv5_checkpoint_path(dirs, id), compress = "xz")
  invisible(payload)
}

ptpv5_load_checkpoint <- function(
  dirs,
  id,
  cfg,
  config_sha,
  input_sha,
  quiet = FALSE
) {
  path <- ptpv5_checkpoint_path(dirs, id)
  if (!file.exists(path)) return(NULL)
  x <- try(readRDS(path), silent = TRUE)
  if (inherits(x, "try-error") || !is.list(x)) return(NULL)
  checks <- c(
    identical(x$status, "complete"),
    identical(x$schema_version, cfg$checkpoint$schema_version),
    identical(x$method_version, cfg$analysis$method_version),
    identical(x$config_checksum, config_sha),
    identical(x$input_checksum, input_sha)
  )
  if (!all(checks)) {
    if (!quiet) message("Ignoring stale V5 checkpoint: ", basename(path))
    return(NULL)
  }
  if (!quiet) message("Resuming V5 checkpoint: ", basename(path))
  x$data
}

ptpv5_manifest_for_paths <- function(paths, ids = names(paths)) {
  if (is.null(ids) || any(!nzchar(ids))) ids <- basename(paths)
  info <- file.info(unname(paths))
  data.frame(
    input_id = ids,
    path = unname(paths),
    exists = file.exists(unname(paths)),
    bytes = as.numeric(info$size),
    sha256 = vapply(unname(paths), ptpv5_file_sha, character(1L)),
    stringsAsFactors = FALSE
  )
}

ptpv5_split_genes <- function(x) {
  x <- as.character(x %||% "")
  unique(ptp_clean_gene_symbols(unlist(strsplit(x, ";", fixed = TRUE))))
}

ptpv5_finite_p <- function(observed, null, two_sided = FALSE) {
  null <- as.numeric(null)
  null <- null[is.finite(null)]
  if (!is.finite(observed) || !length(null)) return(NA_real_)
  if (two_sided) {
    observed <- abs(observed)
    null <- abs(null)
  }
  (1 + sum(null >= observed)) / (1 + length(null))
}

ptpv5_empirical_rank_p_matrix <- function(stat_matrix, larger = TRUE) {
  stat_matrix <- as.matrix(stat_matrix)
  B <- nrow(stat_matrix)
  out <- matrix(NA_real_, nrow = B, ncol = ncol(stat_matrix))
  for (j in seq_len(ncol(stat_matrix))) {
    x <- stat_matrix[, j]
    ok <- is.finite(x)
    if (!any(ok)) next
    if (larger) {
      out[ok, j] <- (1 + vapply(
        x[ok],
        function(v) sum(x[ok] >= v),
        integer(1L)
      )) / (1 + sum(ok))
    } else {
      out[ok, j] <- (1 + vapply(
        x[ok],
        function(v) sum(x[ok] <= v),
        integer(1L)
      )) / (1 + sum(ok))
    }
  }
  out
}

ptpv5_true_maxmean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(mean(pmax(x, 0)), mean(pmax(-x, 0)))
}

ptpv5_vector_components <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(true_maxmean = NA_real_, mean_square = NA_real_, max_absolute = NA_real_))
  }
  c(
    true_maxmean = ptpv5_true_maxmean(x),
    mean_square = mean(x^2),
    max_absolute = max(abs(x))
  )
}

ptpv5_nested_component_test <- function(stat_matrix, observed_index) {
  stat_matrix <- as.matrix(stat_matrix)
  component_p <- ptpv5_empirical_rank_p_matrix(stat_matrix, larger = TRUE)
  min_p <- apply(component_p, 1L, function(x) {
    x <- x[is.finite(x)]
    if (length(x)) min(x) else NA_real_
  })
  adaptive_p <- ptpv5_empirical_rank_p_matrix(matrix(-min_p, ncol = 1L))[, 1L]
  list(
    component_p = component_p,
    min_component_p = min_p,
    adaptive_assignment_p = adaptive_p,
    observed_adaptive_p = adaptive_p[[observed_index]]
  )
}

ptpv5_assignment_matrix <- function(assignments, sample_ids) {
  out <- matrix(
    FALSE,
    nrow = length(sample_ids),
    ncol = nrow(assignments),
    dimnames = list(sample_ids, assignments$assignment_id)
  )
  for (j in seq_len(nrow(assignments))) {
    treated <- strsplit(
      assignments$treated_samples[[j]],
      ";",
      fixed = TRUE
    )[[1L]]
    out[, j] <- sample_ids %in% treated
  }
  out
}

ptpv5_subset_assignment_bundle <- function(bundle, args) {
  if (!isTRUE(args$smoke)) return(bundle)
  B <- nrow(bundle$assignments)
  obs <- bundle$observed_index
  target <- max(4L, min(as.integer(args$smoke_assignments), B))
  keep <- unique(c(obs, seq_len(target)))
  keep <- keep[seq_len(min(length(keep), target))]
  if (!obs %in% keep) keep[[length(keep)]] <- obs
  keep <- sort(unique(keep))
  bundle$assignments <- bundle$assignments[keep, , drop = FALSE]
  bundle$observed_index <- match(obs, keep)
  bundle$assignment_indices_hdf5 <- keep
  bundle$is_smoke <- TRUE
  bundle
}

ptpv5_load_bundle <- function(cfg, args, dirs) {
  input_paths <- unlist(cfg$inputs, use.names = TRUE)
  missing <- input_paths[!file.exists(input_paths)]
  if (length(missing)) {
    stop(
      "Required V5 inputs are missing: ",
      paste(names(missing), missing, sep = "=", collapse = "; "),
      call. = FALSE
    )
  }
  input_manifest <- ptpv5_manifest_for_paths(input_paths, names(input_paths))
  ptp_write_csv(input_manifest, file.path(dirs$manifest, "v5_input_checksums_runtime.csv"))
  input_sha <- ptpv5_object_sha(input_manifest[, c("input_id", "sha256")])

  gene_cp <- readRDS(cfg$inputs$v4_gene_statistics_checkpoint)$data
  mouse_cp <- readRDS(cfg$inputs$v4_mouse_expression_checkpoint)$data
  memberships <- read.csv(
    cfg$inputs$v4_unique_memberships,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  alias <- read.csv(
    cfg$inputs$v4_program_alias,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  roles <- read.csv(
    cfg$inputs$v4_program_roles,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  families <- read.csv(
    cfg$inputs$v4_response_families,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  weights <- read.csv(
    cfg$inputs$v4_gene_weights,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  assignments <- read.csv(
    cfg$inputs$v4_assignments,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  observed_index <- which(
    tolower(as.character(assignments$is_observed_assignment)) %in%
      c("true", "t", "1")
  )
  if (length(observed_index) != 1L) {
    stop("V5 requires exactly one observed assignment.", call. = FALSE)
  }
  if (!identical(as.character(assignments$assignment_id), as.character(gene_cp$assignments$assignment_id))) {
    stop("V4 assignment table and gene checkpoint assignment order differ.", call. = FALSE)
  }
  gene_ids <- as.character(gene_cp$observed$gene)
  gene_symbols <- ptp_clean_gene_symbols(gene_cp$observed$gene_symbol)
  symbol_to_row <- split(seq_along(gene_symbols), gene_symbols)
  symbol_to_row <- vapply(
    symbol_to_row,
    function(i) i[[1L]],
    integer(1L)
  )
  sample_ids <- as.character(mouse_cp$metadata$sample_id)
  if (!identical(sample_ids, colnames(mouse_cp$expression))) {
    mouse_cp$expression <- mouse_cp$expression[, sample_ids, drop = FALSE]
  }
  checks <- data.frame(
    check = c(
      "mice", "assignments", "genes", "unique_memberships",
      "program_labels", "response_families", "observed_assignment"
    ),
    observed = c(
      length(sample_ids), nrow(assignments), length(gene_ids),
      nrow(memberships), nrow(alias), nrow(families), length(observed_index)
    ),
    expected = c(
      cfg$acceptance$expected_mice,
      cfg$acceptance$expected_primary_assignments,
      cfg$acceptance$expected_genes,
      cfg$acceptance$expected_unique_memberships,
      cfg$acceptance$expected_program_labels,
      cfg$acceptance$expected_response_families,
      1L
    ),
    stringsAsFactors = FALSE
  )
  checks$status <- ifelse(checks$observed == checks$expected, "ok", "failed")
  ptp_write_csv(checks, file.path(dirs$qc, "frozen_input_dimension_audit.csv"))
  if (any(checks$status != "ok")) {
    stop("Frozen V5 input dimensions failed validation.", call. = FALSE)
  }
  bundle <- list(
    gene_checkpoint = gene_cp,
    mouse = mouse_cp,
    memberships = memberships,
    alias = alias,
    roles = roles,
    families = families,
    weights = weights,
    assignments = assignments,
    observed_index = observed_index,
    assignment_indices_hdf5 = seq_len(nrow(assignments)),
    sample_ids = sample_ids,
    gene_ids = gene_ids,
    gene_symbols = gene_symbols,
    symbol_to_row = symbol_to_row,
    input_manifest = input_manifest,
    input_sha = input_sha,
    is_smoke = FALSE
  )
  ptpv5_subset_assignment_bundle(bundle, args)
}

ptpv5_set_indices <- function(genes, bundle) {
  genes <- ptp_clean_gene_symbols(genes)
  idx <- unname(bundle$symbol_to_row[genes])
  sort(unique(idx[is.finite(idx)]))
}

ptpv5_frozen_sets <- function(bundle, args) {
  memberships <- bundle$memberships
  families <- bundle$families
  if (isTRUE(args$smoke)) {
    memberships <- head(memberships, as.integer(args$smoke_memberships))
    families <- head(families, as.integer(args$smoke_families))
  }
  membership_sets <- lapply(memberships$genes, function(x) {
    ptpv5_set_indices(ptpv5_split_genes(x), bundle)
  })
  names(membership_sets) <- memberships$membership_id
  family_sets <- lapply(families$union_genes, function(x) {
    ptpv5_set_indices(ptpv5_split_genes(x), bundle)
  })
  names(family_sets) <- families$response_family
  list(
    memberships = memberships,
    families = families,
    membership_sets = membership_sets,
    family_sets = family_sets
  )
}

ptpv5_add_fdr <- function(results, p_col, prefix = "fdr") {
  results[[paste0(prefix, "_all")]] <- stats::p.adjust(results[[p_col]], "BH")
  if ("tiers" %in% names(results)) {
    focused <- grepl("tier1|tier2", results$tiers)
    tier3 <- grepl("tier3", results$tiers) & !focused
    control <- grepl("negative_control", results$response_family %||% "")
    results[[paste0(prefix, "_focused")]] <- NA_real_
    results[[paste0(prefix, "_tier3")]] <- NA_real_
    results[[paste0(prefix, "_controls")]] <- NA_real_
    results[[paste0(prefix, "_focused")]][focused] <-
      stats::p.adjust(results[[p_col]][focused], "BH")
    results[[paste0(prefix, "_tier3")]][tier3] <-
      stats::p.adjust(results[[p_col]][tier3], "BH")
    results[[paste0(prefix, "_controls")]][control] <-
      stats::p.adjust(results[[p_col]][control], "BH")
  }
  results
}

ptpv5_write_method_status <- function(
  dirs,
  method_id,
  status,
  note,
  n_tests = NA_integer_,
  path = NULL
) {
  row <- data.frame(
    method_id = method_id,
    status = status,
    note = note,
    n_tests = n_tests,
    completed_at = ptpv5_now(),
    stringsAsFactors = FALSE
  )
  if (is.null(path)) {
    path <- file.path(dirs$qc, paste0(method_id, "_status.csv"))
  }
  ptp_write_csv(row, path)
  row
}

ptpv5_read_method_statuses <- function(dirs) {
  files <- list.files(
    dirs$qc,
    pattern = "_status\\.csv$",
    full.names = TRUE
  )
  ptpv5_bind_rows(lapply(files, function(path) {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }))
}

ptpv5_write_session_info <- function(dirs, args, cfg, config_sha, input_sha) {
  capture.output(sessionInfo(), file = file.path(dirs$manifest, "sessionInfo.txt"))
  ptp_write_csv(
    data.frame(
      key = c(
        "created_at", "config", "output_root", "config_checksum",
        "input_checksum", "method_version", "schema_version",
        "smoke", "workers", "master_seed"
      ),
      value = c(
        ptpv5_now(), args$config, args$output_root, config_sha,
        input_sha, cfg$analysis$method_version, cfg$analysis$schema_version,
        isTRUE(args$smoke), args$workers, args$seed
      ),
      stringsAsFactors = FALSE
    ),
    file.path(dirs$manifest, "run_manifest.csv")
  )
}

ptpv5_copy_frozen_definitions <- function(bundle, cfg, dirs) {
  map <- c(
    v4_unique_memberships = "unique_memberships.csv",
    v4_program_alias = "program_label_to_membership_alias.csv",
    v4_program_roles = "program_universe_role_eligibility.csv",
    v4_response_families = "response_family_definitions.csv",
    v4_assignments = "primary_4900_assignments.csv"
  )
  for (id in names(map)) {
    file.copy(cfg$inputs[[id]], file.path(dirs$frozen, map[[id]]), overwrite = TRUE)
  }
  ptp_write_csv(
    data.frame(
      policy = c("program_score_inference", "external_cohort", "species_filtering"),
      status = c("forbidden", "disabled", "not_added"),
      stringsAsFactors = FALSE
    ),
    file.path(dirs$frozen, "v5_analysis_policy.csv")
  )
}

ptpv5_run_workflow <- function(args, repo_root, script_dir, module_dir) {
  ptpv5_thread_env()
  options(stringsAsFactors = FALSE, warn = 1)
  args$output_root <- normalizePath(
    ptpv5_prepare_output(args$output_root, args$overwrite),
    winslash = "/",
    mustWork = FALSE
  )
  dirs <- ptpv5_dirs(args$output_root)
  log_file <- file.path(dirs$logs, "v5_workflow.log")
  cfg <- ptpv5_read_config(args$config)
  config_audit <- ptpv5_validate_config(cfg)
  ptp_write_csv(config_audit, file.path(dirs$qc, "config_schema_audit.csv"))
  config_sha <- ptpv5_config_sha(args, cfg, module_dir)
  set.seed(as.integer(args$seed))
  ptpv5_message("Loading and validating frozen V4 inputs.", log_file = log_file)
  bundle <- ptpv5_load_bundle(cfg, args, dirs)
  ptpv5_write_session_info(
    dirs, args, cfg, config_sha, bundle$input_sha
  )
  ptpv5_copy_frozen_definitions(bundle, cfg, dirs)
  context <- list(
    args = args,
    cfg = cfg,
    dirs = dirs,
    bundle = bundle,
    repo_root = repo_root,
    script_dir = script_dir,
    module_dir = module_dir,
    log_file = log_file,
    config_sha = config_sha,
    input_sha = bundle$input_sha
  )

  runners <- list(
    exact_ranked_enrichment = ptpv5_run_exact_ranked_enrichment,
    state_decomposition = ptpv5_run_state_decomposition,
    bayesian_gene_vector = ptpv5_run_bayesian_gene_vector,
    exact_global_kernel = ptpv5_run_exact_global_kernel,
    gene_wise_trajectory = ptpv5_run_gene_wise_trajectory
  )
  results <- list()
  for (id in names(runners)) {
    ptpv5_message("Starting method: ", id, log_file = log_file)
    results[[id]] <- tryCatch(
      runners[[id]](context),
      error = function(e) {
        ptpv5_write_method_status(
          dirs, id, "failed", conditionMessage(e)
        )
        saveRDS(
          list(error = conditionMessage(e), calls = sys.calls()),
          file.path(dirs$checkpoints, paste0(id, "_failure.rds"))
        )
        stop("V5 method failed [", id, "]: ", conditionMessage(e), call. = FALSE)
      }
    )
  }
  ptpv5_message("Running calibration and cross-method synthesis.", log_file = log_file)
  synthesis <- ptpv5_run_simulation_and_synthesis(context, results)
  ptpv5_message("Building report and final audit.", log_file = log_file)
  report <- ptpv5_run_report_and_audit(context, results, synthesis)
  ptpv5_message("V5 workflow complete.", log_file = log_file)
  invisible(list(results = results, synthesis = synthesis, report = report))
}
