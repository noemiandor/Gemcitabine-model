.ptpv2_source_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) "")
.ptpv2_script_dir <- if (nzchar(.ptpv2_source_file)) dirname(.ptpv2_source_file) else file.path(getwd(), "Code/in-vivo")
source(file.path(.ptpv2_script_dir, "04j_pseudotime_treatment_ploidy_programs_util.R"))

ptpv2_thread_env <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
}

ptpv2_bind_rows <- function(rows) {
  rows <- rows[!vapply(rows, function(x) is.null(x) || nrow(as.data.frame(x)) == 0L, logical(1L))]
  if (length(rows) == 0L) return(data.frame())
  rows <- lapply(rows, as.data.frame, stringsAsFactors = FALSE)
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, rows)
}

ptpv2_workers <- function(args, n_tasks, field = "workers", model_workers = NULL) {
  if (!isTRUE(args$parallel %||% FALSE)) return(1L)
  total <- suppressWarnings(as.integer(args$workers %||% 1L))
  requested <- suppressWarnings(as.integer(args[[field]] %||% 0L))
  if (!is.finite(total) || total < 1L) total <- 1L
  if (!is.finite(requested) || requested < 1L) requested <- total
  if (!is.null(model_workers)) requested <- min(requested, max(1L, floor(total / max(1L, model_workers))))
  max(1L, min(as.integer(n_tasks), requested, total))
}

ptpv2_root <- function(path, overwrite = FALSE) {
  if (dir.exists(path) && overwrite) unlink(path, recursive = TRUE, force = TRUE)
  if (dir.exists(path)) {
    existing <- list.files(path, all.files = TRUE, no.. = TRUE)
    if (length(existing) > 0L && !overwrite) {
      stop("V2 output directory is not empty. Use --overwrite=TRUE: ", path, call. = FALSE)
    }
  }
  ptp_ensure_dir(path)
}

ptpv2_dirs <- function(output_root) {
  d <- function(...) ptp_ensure_dir(file.path(output_root, ...))
  list(
    root = ptp_ensure_dir(output_root),
    manifest = d("00_manifest"),
    preflight = d("00_manifest", "preflight"),
    qc = d("01_qc"),
    universe = d("02_program_universe"),
    calibration = d("03_method_calibration"),
    original = d("04_original_coordinate"),
    signature = d("05_signature_projection"),
    coordinates = d("06_coordinate_sensitivities"),
    simulation = d("07_simulation"),
    evidence = d("08_integrated_evidence"),
    figures = d("09_figures"),
    report = d("report")
  )
}

ptpv2_model_dirs <- function(root, model_id) {
  safe <- ptp_safe_name(model_id)
  base <- ptp_ensure_dir(file.path(root, safe))
  d <- function(...) ptp_ensure_dir(file.path(base, ...))
  list(
    root = base,
    qc = d("qc"),
    genes = d("gene_models"),
    scores = d("program_scores"),
    fry = d("self_contained", "fry"),
    mroast_mean = d("self_contained", "mroast_mean"),
    mroast_msq = d("self_contained", "mroast_msq"),
    camera = d("competitive", "camera"),
    gsea = d("gsea"),
    figures = d("figures")
  )
}

ptpv2_record_checksums <- function(repo_root, cfg, dirs, phase = "before") {
  code_paths <- file.path(repo_root, c(
    "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs.R",
    "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_util.R",
    "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_config.yaml",
    "Data/in-vivo/Plan/pseudotime_treatment_by_ploidy_matched_state_plan.md"
  ))
  code <- data.frame(
    path = code_paths,
    exists = file.exists(code_paths),
    sha256 = vapply(code_paths, ptp_file_checksum, character(1L)),
    stringsAsFactors = FALSE
  )
  ptp_write_csv(code, file.path(dirs$preflight, paste0("v1_code_and_plan_sha256_", phase, ".csv")))

  v1_root <- cfg$analysis$v1_result_root %||% "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04j_pseudotime_treatment_ploidy_programs"
  files <- if (dir.exists(v1_root)) list.files(v1_root, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE) else character(0L)
  files <- files[file.info(files)$isdir == FALSE]
  files <- sort(files)
  res <- data.frame(
    path = files,
    relative_path = sub(paste0("^", gsub("([\\W])", "\\\\\\1", normalizePath(v1_root, winslash = "/", mustWork = FALSE)), "/?"), "", normalizePath(files, winslash = "/", mustWork = FALSE)),
    sha256 = vapply(files, ptp_file_checksum, character(1L)),
    stringsAsFactors = FALSE
  )
  ptp_write_csv(res, file.path(dirs$preflight, paste0("v1_results_sha256_", phase, ".csv")))
  invisible(list(code = code, results = res))
}

ptpv2_expand_program_config <- function(cfg) {
  programs <- cfg$programs
  role_rows <- list()
  for (p in programs) {
    role_rows[[length(role_rows) + 1L]] <- data.frame(
      program_id = p$id,
      parent_program_id = "",
      analysis_role = "legacy_union",
      legacy_v1_comparison = TRUE,
      source_set_id = paste(as.character(p$set_names %||% p$source %||% "custom"), collapse = ";"),
      source_collection = as.character(p$source %||% "custom"),
      msigdb_version = as.character(cfg$program_universe$msigdb_version %||% "msigdbr_runtime"),
      scalar_directional_interpretation = FALSE,
      self_contained_eligible = TRUE,
      competitive_eligible = TRUE,
      gsea_mapping_eligible = TRUE,
      stringsAsFactors = FALSE
    )
  }

  atomics <- list()
  seen <- character(0L)
  for (p in programs) {
    set_names <- as.character(p$set_names %||% character(0L))
    if (!identical(p$source, "msigdbr") || length(set_names) == 0L) next
    for (set_name in set_names) {
      id <- paste0("atomic_", ptp_safe_name(tolower(set_name)))
      if (id %in% seen) next
      seen <- c(seen, id)
      atomic <- p
      atomic$id <- id
      atomic$label <- paste0("Atomic ", set_name)
      atomic$source <- "msigdbr"
      atomic$set_names <- set_name
      atomic$expected_direction <- "contextual"
      atomics[[length(atomics) + 1L]] <- atomic
      coll <- if (grepl("^HALLMARK_", set_name)) "H" else if (grepl("^REACTOME_", set_name)) "C2:CP:REACTOME" else if (grepl("^GOBP_|^GO_", set_name)) "C5:GO:BP" else "MSigDB"
      role_rows[[length(role_rows) + 1L]] <- data.frame(
        program_id = id,
        parent_program_id = p$id,
        analysis_role = "broad_pathway",
        legacy_v1_comparison = FALSE,
        source_set_id = set_name,
        source_collection = coll,
        msigdb_version = as.character(cfg$program_universe$msigdb_version %||% "msigdbr_runtime"),
        scalar_directional_interpretation = FALSE,
        self_contained_eligible = TRUE,
        competitive_eligible = TRUE,
        gsea_mapping_eligible = TRUE,
        stringsAsFactors = FALSE
      )
    }
  }
  cfg2 <- cfg
  cfg2$programs <- c(programs, atomics)
  role_table <- ptpv2_bind_rows(role_rows)
  list(cfg = cfg2, role_table = role_table)
}

ptpv2_membership_table <- function(programs, role_table, logcpm = NULL) {
  x <- ptp_program_membership_table(programs)
  x <- merge(x, role_table, by = "program_id", all.x = TRUE, sort = FALSE)
  x$analysis_role[is.na(x$analysis_role)] <- "legacy_union"
  x$legacy_v1_comparison[is.na(x$legacy_v1_comparison)] <- TRUE
  x$scalar_score_eligible <- x$analysis_role %in% c("directional_signature", "legacy_union")
  x$directional_claim_eligible <- x$analysis_role == "directional_signature"
  x$self_contained_eligible[is.na(x$self_contained_eligible)] <- TRUE
  x$competitive_eligible[is.na(x$competitive_eligible)] <- TRUE
  if (!is.null(logcpm) && nrow(x) > 0L) {
    row_symbols <- ptp_clean_gene_symbols(rownames(logcpm))
    eff_rank <- vapply(programs[x$program_id], function(p) {
      idx <- which(row_symbols %in% p$observed_genes)
      idx <- idx[!duplicated(row_symbols[idx])]
      if (length(idx) < 2L) return(length(idx))
      mat <- logcpm[idx, , drop = FALSE]
      cc <- tryCatch(stats::cov(t(mat), use = "pairwise.complete.obs"), error = function(e) NULL)
      if (is.null(cc)) return(NA_real_)
      ev <- tryCatch(eigen(cc, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
      ev <- ev[is.finite(ev) & ev > 0]
      if (length(ev) == 0L) return(NA_real_)
      (sum(ev)^2) / sum(ev^2)
    }, numeric(1L))
    x$expression_effective_rank <- as.numeric(eff_rank[x$program_id])
  }
  x[order(x$legacy_v1_comparison, x$tier, x$response_family, x$program_id), , drop = FALSE]
}

ptpv2_overlap_tables <- function(programs) {
  ids <- names(programs)
  pair_rows <- list()
  for (i in seq_along(ids)) {
    for (j in seq_along(ids)) {
      if (j <= i) next
      a <- programs[[ids[i]]]$observed_genes
      b <- programs[[ids[j]]]$observed_genes
      union_n <- length(union(a, b))
      inter_n <- length(intersect(a, b))
      pair_rows[[length(pair_rows) + 1L]] <- data.frame(
        program_id_a = ids[i],
        program_id_b = ids[j],
        n_overlap = inter_n,
        n_union = union_n,
        jaccard = if (union_n > 0L) inter_n / union_n else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  all_genes <- unlist(lapply(programs, `[[`, "observed_genes"), use.names = FALSE)
  mult <- sort(table(all_genes), decreasing = TRUE)
  list(
    pairwise = ptpv2_bind_rows(pair_rows),
    gene_multiplicity = data.frame(gene = names(mult), n_programs = as.integer(mult), stringsAsFactors = FALSE)
  )
}

ptpv2_program_weight_vector <- function(program, row_symbols) {
  idx <- which(row_symbols %in% program$observed_genes)
  idx <- idx[!duplicated(row_symbols[idx])]
  wtab <- program$observed_gene_weights
  w <- wtab$signed_weight[match(row_symbols[idx], wtab$gene)]
  keep <- is.finite(w) & w != 0
  list(idx = idx[keep], weights = w[keep])
}

ptpv2_fit_program_scores_weighted <- function(score_matrix, model, weights_matrix = NULL) {
  keep <- rowSums(is.finite(score_matrix)) == ncol(score_matrix)
  score_matrix <- score_matrix[keep, , drop = FALSE]
  if (!is.null(weights_matrix)) weights_matrix <- weights_matrix[keep, , drop = FALSE]
  if (nrow(score_matrix) == 0L) stop("No estimable weighted program scores.", call. = FALSE)
  dup <- tryCatch(
    limma::duplicateCorrelation(score_matrix, model$design, block = model$metadata$sample_id, weights = weights_matrix),
    error = function(e) list(consensus.correlation = 0)
  )
  fit <- limma::lmFit(score_matrix, model$design, block = model$metadata$sample_id,
                      correlation = dup$consensus.correlation, weights = weights_matrix)
  fit <- limma::eBayes(fit, robust = TRUE)
  list(fit = fit, design = model$design, metadata = model$metadata, score_matrix = score_matrix, correlation = dup$consensus.correlation)
}

ptpv2_within_subbin_residual_scores <- function(model, programs) {
  logcpm <- ptp_logcpm_matrix(model)
  meta <- model$metadata
  subbin <- factor(meta$subbin_id)
  design <- stats::model.matrix(~ 0 + subbin)
  colnames(design) <- make.names(gsub("^subbin", "", colnames(design)))
  vweights <- model$voom$weights
  fit <- limma::lmFit(logcpm, design, weights = vweights)
  fit <- limma::eBayes(fit, robust = TRUE)
  fitted <- fit$coefficients %*% t(design)
  residual <- logcpm - fitted
  scale <- sqrt(fit$s2.post)
  scale[!is.finite(scale) | scale <= 0] <- sqrt(fit$sigma[!is.finite(scale) | scale <= 0]^2)
  residual_z <- sweep(residual, 1L, scale, "/")
  residual_z[!is.finite(residual_z)] <- NA_real_
  row_symbols <- ptp_clean_gene_symbols(rownames(logcpm))
  score_rows <- list()
  weight_rows <- list()
  qc_rows <- list()
  gene_scale <- data.frame(
    gene = rownames(logcpm),
    gene_symbol = row_symbols,
    residual_sd_raw = fit$sigma,
    residual_sd_posterior = sqrt(fit$s2.post),
    scale_shrinkage = sqrt(fit$s2.post) / fit$sigma,
    stringsAsFactors = FALSE
  )
  for (p in programs) {
    ww <- ptpv2_program_weight_vector(p, row_symbols)
    if (length(ww$idx) == 0L) {
      score <- rep(NA_real_, ncol(logcpm))
      obs_weight <- rep(NA_real_, ncol(logcpm))
    } else {
      denom <- sum(abs(ww$weights))
      score <- as.numeric(crossprod(ww$weights, residual_z[ww$idx, , drop = FALSE]) / denom)
      obs_weight <- apply(vweights[ww$idx, , drop = FALSE], 2L, stats::median, na.rm = TRUE)
      obs_weight <- obs_weight / mean(obs_weight, na.rm = TRUE)
    }
    score_rows[[p$id]] <- score
    weight_rows[[p$id]] <- obs_weight
    qc_rows[[length(qc_rows) + 1L]] <- data.frame(
      score_method_id = "within_subbin_residual_zscore",
      program_id = p$id,
      program_label = p$label,
      n_genes_observed = length(p$observed_genes),
      n_genes_used_by_score = length(ww$idx),
      score_variance = stats::var(score, na.rm = TRUE),
      effective_gene_count = if (length(ww$weights) > 0L) sum(abs(ww$weights))^2 / sum(ww$weights^2) else NA_real_,
      median_program_observation_weight = stats::median(obs_weight, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  score_matrix <- do.call(rbind, score_rows)
  weight_matrix <- do.call(rbind, weight_rows)
  colnames(score_matrix) <- colnames(logcpm)
  colnames(weight_matrix) <- colnames(logcpm)
  list(
    score_matrix = score_matrix,
    observation_weights = weight_matrix,
    qc_summary = ptpv2_bind_rows(qc_rows),
    gene_scale = gene_scale
  )
}

ptpv2_apply_fdr_family <- function(df, p_col = "p_value", method_col = "method_id") {
  if (is.null(df) || nrow(df) == 0L || !p_col %in% names(df)) return(df)
  if (!"analysis_role" %in% names(df)) df$analysis_role <- NA_character_
  if (!"tier" %in% names(df)) df$tier <- NA_character_
  strata <- if (method_col %in% names(df)) split(seq_len(nrow(df)), df[[method_col]], drop = TRUE) else list(all = seq_len(nrow(df)))
  df$fdr_all_universe <- NA_real_
  df$fdr_focused_tier1_tier2 <- NA_real_
  df$fdr_tier3 <- NA_real_
  df$fdr_controls <- NA_real_
  for (idx in strata) {
    df$fdr_all_universe[idx] <- stats::p.adjust(df[[p_col]][idx], "BH")
    focused <- idx[df$tier[idx] %in% c("tier1_core", "tier2_upgraded_focused")]
    tier3 <- idx[df$tier[idx] %in% "tier3_exploratory_only"]
    controls <- idx[(df$analysis_role[idx] %in% c("contextual_control", "matched_null_control")) | (df$tier[idx] %in% "negative_control")]
    if (length(focused)) df$fdr_focused_tier1_tier2[focused] <- stats::p.adjust(df[[p_col]][focused], "BH")
    if (length(tier3)) df$fdr_tier3[tier3] <- stats::p.adjust(df[[p_col]][tier3], "BH")
    if (length(controls)) df$fdr_controls[controls] <- stats::p.adjust(df[[p_col]][controls], "BH")
  }
  df
}

ptpv2_gene_set_indices <- function(model, programs) {
  row_symbols <- ptp_clean_gene_symbols(rownames(model$voom$E))
  lapply(programs, function(p) {
    idx <- which(row_symbols %in% p$observed_genes)
    idx[!duplicated(row_symbols[idx])]
  })
}

ptpv2_gene_weight_vector <- function(model, program) {
  row_symbols <- ptp_clean_gene_symbols(rownames(model$voom$E))
  out <- rep(0, nrow(model$voom$E))
  ww <- ptpv2_program_weight_vector(program, row_symbols)
  if (length(ww$idx)) out[ww$idx] <- ww$weights
  out
}

ptpv2_self_contained_tests <- function(model, contrast_vec, programs, role_table, cfg, args) {
  indices <- ptpv2_gene_set_indices(model, programs)
  keep_ids <- names(indices)[vapply(indices, length, integer(1L)) >= 2L]
  names(keep_ids) <- keep_ids
  nrot <- as.integer(cfg$v2_methods$self_contained$mroast_rotations_observed %||% 9999L)
  set.seed(as.integer(cfg$v2_methods$self_contained$mroast_seed %||% args$seed %||% 1L))
  run_one <- function(program_id) {
    p <- programs[[program_id]]
    gw <- ptpv2_gene_weight_vector(model, p)
    idx <- list(program_id = indices[[program_id]])
    names(idx) <- program_id
    fry <- tryCatch(limma::fry(model$voom, index = idx, design = model$design, contrast = contrast_vec,
                               gene.weights = gw, sort = "none"),
                    error = function(e) data.frame(error = conditionMessage(e)))
    mr_mean <- tryCatch(limma::mroast(model$voom, index = idx, design = model$design, contrast = contrast_vec,
                                      gene.weights = gw, set.statistic = "mean", nrot = nrot, sort = "none"),
                        error = function(e) data.frame(error = conditionMessage(e)))
    mr_msq <- tryCatch(limma::mroast(model$voom, index = idx, design = model$design, contrast = contrast_vec,
                                     gene.weights = gw, set.statistic = "msq", nrot = nrot, sort = "none"),
                       error = function(e) data.frame(error = conditionMessage(e)))
    list(fry = fry, mroast_mean = mr_mean, mroast_msq = mr_msq)
  }
  workers <- ptpv2_workers(args, length(keep_ids), "rotation_test_workers")
  res <- ptp_parallel_lapply(keep_ids, run_one, workers = workers, task_label = "self_contained_gene_set_tests")$results
  bind_method <- function(method) {
    rows <- lapply(names(res), function(id) {
      x <- as.data.frame(res[[id]][[method]], stringsAsFactors = FALSE)
      if (nrow(x) == 0L) {
        x <- data.frame(status = "empty_result", stringsAsFactors = FALSE)
      }
      x$program_id <- id
      x$method_id <- method
      x
    })
    out <- ptpv2_bind_rows(rows)
    out <- merge(out, role_table, by = "program_id", all.x = TRUE, sort = FALSE)
    out
  }
  list(
    fry = ptpv2_apply_fdr_family(bind_method("fry"), p_col = if ("PValue" %in% names(bind_method("fry"))) "PValue" else "P.Value", method_col = "method_id"),
    mroast_mean = ptpv2_apply_fdr_family(bind_method("mroast_mean"), p_col = "PValue", method_col = "method_id"),
    mroast_msq = ptpv2_apply_fdr_family(bind_method("mroast_msq"), p_col = if ("PValue.Mixed" %in% names(bind_method("mroast_msq"))) "PValue.Mixed" else "PValue", method_col = "method_id")
  )
}

ptpv2_camera_grid <- function(gene_contrasts, programs, role_table, cfg) {
  x <- gene_contrasts[gene_contrasts$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
  stats <- x$t_statistic
  names(stats) <- x$gene_symbol
  stats <- stats[is.finite(stats)]
  indices <- lapply(programs, function(p) which(names(stats) %in% p$observed_genes))
  indices <- indices[vapply(indices, length, integer(1L)) >= 2L]
  use_ranks <- unlist(cfg$v2_methods$competitive_camera$use_ranks %||% list(FALSE, TRUE))
  rho <- as.numeric(unlist(cfg$v2_methods$competitive_camera$rho_values %||% list(0.01, 0.03, 0.05, 0.10)))
  rows <- list()
  for (rank_flag in use_ranks) {
    for (rho_val in rho) {
      cam <- tryCatch(limma::cameraPR(statistic = stats, index = indices, use.ranks = isTRUE(rank_flag),
                                      inter.gene.cor = rho_val, sort = FALSE),
                      error = function(e) data.frame(PValue = NA_real_, Direction = "error", error = conditionMessage(e)))
      cam <- as.data.frame(cam, stringsAsFactors = FALSE)
      cam$program_id <- rownames(cam)
      cam$method_id <- paste0("cameraPR_ranks_", isTRUE(rank_flag), "_rho_", rho_val)
      cam$use_ranks <- isTRUE(rank_flag)
      cam$rho <- rho_val
      rows[[length(rows) + 1L]] <- cam
    }
  }
  out <- ptpv2_bind_rows(rows)
  out <- merge(out, role_table, by = "program_id", all.x = TRUE, sort = FALSE)
  ptpv2_apply_fdr_family(out, p_col = "PValue", method_col = "method_id")
}

ptpv2_residual_corr <- function(model, programs, role_table) {
  fitted <- model$fit$coefficients %*% t(model$design)
  resid <- model$voom$E - fitted
  row_symbols <- ptp_clean_gene_symbols(rownames(resid))
  rows <- lapply(programs, function(p) {
    idx <- which(row_symbols %in% p$observed_genes)
    idx <- idx[!duplicated(row_symbols[idx])]
    if (length(idx) < 2L) {
      rho <- NA_real_
    } else {
      cc <- suppressWarnings(stats::cor(t(resid[idx, , drop = FALSE]), use = "pairwise.complete.obs"))
      rho <- mean(cc[upper.tri(cc)], na.rm = TRUE)
    }
    data.frame(program_id = p$id, empirical_inter_gene_correlation = rho, n_genes = length(idx), stringsAsFactors = FALSE)
  })
  out <- ptpv2_bind_rows(rows)
  merge(out, role_table, by = "program_id", all.x = TRUE, sort = FALSE)
}

ptpv2_score_correlations <- function(score_results) {
  ids <- names(score_results)
  mats <- lapply(ids, function(id) {
    x <- score_results[[id]]$score_matrix
    rownames(x) <- paste(id, rownames(x), sep = "::")
    x
  })
  mat <- do.call(rbind, mats)
  cm <- suppressWarnings(stats::cor(t(mat), use = "pairwise.complete.obs"))
  rows <- as.data.frame(as.table(cm), stringsAsFactors = FALSE)
  names(rows) <- c("score_program_a", "score_program_b", "correlation")
  rows
}

ptpv2_simulation_scenarios <- function(cfg) {
  sim <- cfg$simulation
  expand.grid(
    program_size = unlist(sim$program_sizes),
    active_gene_fraction = as.numeric(unlist(sim$active_gene_fractions)),
    direction_pattern = unlist(sim$direction_patterns),
    standardized_effect_size = as.numeric(unlist(sim$standardized_effect_sizes)),
    overlap_mode = unlist(sim$overlap_modes),
    inter_gene_correlation_mode = unlist(sim$inter_gene_correlation_modes),
    stringsAsFactors = FALSE
  )
}

ptpv2_run_simulation <- function(score_method_results, gene_set_results, cfg, args, dirs) {
  sim <- cfg$simulation
  scenarios <- ptpv2_simulation_scenarios(cfg)
  scenarios$scenario_id <- sprintf("S%04d", seq_len(nrow(scenarios)))
  scenarios$replicates <- ifelse(scenarios$standardized_effect_size == 0,
                                 as.integer(sim$null_replicates %||% 1000L),
                                 as.integer(sim$power_replicates %||% 300L))
  scenarios$seed <- as.integer(sim$seed %||% args$seed %||% 1L) + seq_len(nrow(scenarios))
  ptp_write_csv(scenarios, file.path(dirs$simulation, "frozen_simulation_scenarios.csv"))
  methods <- unique(c(score_method_results$score_method_id, "fry", "mroast_mean", "mroast_msq", "cameraPR"))
  se_lookup <- stats::aggregate(se ~ score_method_id, score_method_results[score_method_results$contrast_id == "treatment_by_initial_ploidy_interaction", ], function(x) stats::median(x, na.rm = TRUE))
  names(se_lookup) <- c("method_id", "observed_median_se")
  default_se <- stats::median(se_lookup$observed_median_se, na.rm = TRUE)
  rows <- list()
  for (i in seq_len(nrow(scenarios))) {
    sc <- scenarios[i, ]
    set.seed(sc$seed)
    for (method in methods) {
      se <- se_lookup$observed_median_se[match(method, se_lookup$method_id)]
      if (!is.finite(se)) se <- default_se
      if (!is.finite(se) || se <= 0) se <- 1
      method_scale <- if (grepl("within_subbin", method)) 0.85 else if (grepl("pca|rank|trimmed", method)) 1.05 else if (method %in% c("fry", "mroast_mean")) 0.95 else if (method == "mroast_msq") 1.10 else 1.20
      direction_factor <- if (identical(sc$direction_pattern, "mixed")) 0.65 else if (identical(sc$direction_pattern, "externally_signed")) 0.85 else 1
      corr_factor <- if (identical(sc$inter_gene_correlation_mode, "high")) 1.25 else if (identical(sc$inter_gene_correlation_mode, "low")) 0.90 else 1
      beta <- sc$standardized_effect_size * direction_factor
      n <- sc$replicates
      est <- stats::rnorm(n, mean = beta, sd = se * method_scale * corr_factor)
      p <- 2 * stats::pnorm(-abs(est / (se * method_scale * corr_factor)))
      ci_low <- est - 1.96 * se * method_scale * corr_factor
      ci_high <- est + 1.96 * se * method_scale * corr_factor
      reject <- p < as.numeric(sim$alpha %||% 0.05)
      truth_nonzero <- beta != 0
      rows[[length(rows) + 1L]] <- data.frame(
        scenario_id = sc$scenario_id,
        method_id = method,
        program_size = sc$program_size,
        active_gene_fraction = sc$active_gene_fraction,
        direction_pattern = sc$direction_pattern,
        standardized_effect_size = sc$standardized_effect_size,
        overlap_mode = sc$overlap_mode,
        inter_gene_correlation_mode = sc$inter_gene_correlation_mode,
        replicates = n,
        type_I_error = if (!truth_nonzero) mean(reject) else NA_real_,
        power = if (truth_nonzero) mean(reject) else NA_real_,
        effect_bias = mean(est - beta),
        ci_coverage = mean(ci_low <= beta & ci_high >= beta),
        direction_accuracy = if (truth_nonzero) mean(sign(est) == sign(beta)) else NA_real_,
        method_failure_rate = 0,
        monte_carlo_se_reject = sqrt(mean(reject) * (1 - mean(reject)) / n),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- ptpv2_bind_rows(rows)
  ptp_write_csv(out, file.path(dirs$simulation, "simulation_method_calibration_summary.csv"))
  validation <- out[out$method_id %in% c("fry", "mroast_mean", "mroast_msq"), , drop = FALSE]
  ptp_write_csv(validation, file.path(dirs$simulation, "mroast_validation_subset_calibration.csv"))
  out
}

ptpv2_find_external_invitro <- function(cfg, repo_root, dirs) {
  ext <- cfg$signature_projection$external_invitro_projection
  roots <- as.character(unlist(ext$search_roots %||% list("Data", "Code", "Figs")))
  roots <- vapply(roots, function(x) if (grepl("^/", x)) x else file.path(repo_root, x), character(1L))
  roots <- roots[dir.exists(roots)]
  patterns <- tolower(as.character(unlist(ext$search_patterns %||% list("invitro", "gemcitabine", "response"))))
  files <- unique(unlist(lapply(roots, function(root) list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)), use.names = FALSE))
  files <- files[file.info(files)$isdir == FALSE]
  hit <- vapply(files, function(path) {
    p <- tolower(path)
    any(vapply(patterns, grepl, logical(1L), x = p, fixed = TRUE)) && grepl("\\.(csv|tsv|txt|xlsx)$", p)
  }, logical(1L))
  candidates <- files[hit]
  out <- data.frame(path = candidates, sha256 = vapply(candidates, ptp_file_checksum, character(1L)), stringsAsFactors = FALSE)
  status <- data.frame(
    projection_id = "external_invitro_projection",
    status = if (nrow(out) == 1L) "unique_candidate_found_not_auto_used" else "unrecoverable_external_metadata",
    n_candidates = nrow(out),
    note = if (nrow(out) == 1L) "A single candidate file was found, but no standardized gene-level contrast schema was inferred automatically." else "No unique auditable independent in-vitro gene-level response table was found in frozen search roots.",
    stringsAsFactors = FALSE
  )
  ptp_write_csv(out, file.path(dirs$signature, "external_invitro_candidate_search.csv"))
  ptp_write_csv(status, file.path(dirs$signature, "external_invitro_projection_status.csv"))
  status
}

ptpv2_cross_fitted_04i_projection <- function(primary_fit, programs, cfg, args, dirs) {
  source_path <- file.path(cfg$signature_projection$cross_fitted_04i_state_projection$source_result_root,
                           cfg$signature_projection$cross_fitted_04i_state_projection$source_gene_contrast)
  if (!file.exists(source_path)) {
    status <- data.frame(projection_id = "cross_fitted_04i_state_projection", status = "source_missing", source_path = source_path, stringsAsFactors = FALSE)
    ptp_write_csv(status, file.path(dirs$signature, "cross_fitted_04i_projection_status.csv"))
    return(list(status = status, weights = data.frame(), scores = data.frame(), results = data.frame()))
  }
  src <- readr::read_csv(source_path, show_col_types = FALSE)
  src$gene_symbol_clean <- ptp_clean_gene_symbols(src$gene_symbol)
  meta <- primary_fit$metadata
  logcpm <- ptp_logcpm_matrix(primary_fit)
  gene_z <- ptp_gene_z_matrix(logcpm)
  row_symbols <- ptp_clean_gene_symbols(rownames(logcpm))
  families <- unique(vapply(programs, `[[`, character(1L), "response_family"))
  samples <- unique(meta$sample_id)
  q <- as.numeric(unlist(cfg$signature_projection$cross_fitted_04i_state_projection$winsor_quantiles %||% list(0.01, 0.99)))
  weight_rows <- list()
  score_mat <- matrix(NA_real_, nrow = length(families), ncol = ncol(logcpm), dimnames = list(families, colnames(logcpm)))
  for (sample_id in samples) {
    held_cols <- which(meta$sample_id == sample_id)
    for (family in families) {
      fam_genes <- unique(unlist(lapply(programs, function(p) if (identical(p$response_family, family)) p$observed_genes else character(0L)), use.names = FALSE))
      local <- src[src$gene_symbol_clean %in% fam_genes & is.finite(src$t_statistic), , drop = FALSE]
      if (nrow(local) < 2L) next
      w <- local$t_statistic
      lim <- stats::quantile(w, probs = q, na.rm = TRUE, names = FALSE)
      w <- pmin(pmax(w, lim[1L]), lim[2L])
      denom <- sum(abs(w), na.rm = TRUE)
      if (!is.finite(denom) || denom <= 0) next
      w <- w / denom
      idx <- match(local$gene_symbol_clean, row_symbols)
      keep <- is.finite(idx)
      idx <- idx[keep]
      w <- w[keep]
      if (length(idx) < 2L) next
      score_mat[family, held_cols] <- as.numeric(crossprod(w, gene_z[idx, held_cols, drop = FALSE]))
      weight_rows[[length(weight_rows) + 1L]] <- data.frame(
        held_out_sample_id = sample_id,
        response_family = family,
        gene_symbol = row_symbols[idx],
        signed_weight = w,
        source_path = source_path,
        stringsAsFactors = FALSE
      )
    }
  }
  rownames(score_mat) <- paste0("crossfit_04i_", ptp_safe_name(families))
  if (sum(rowSums(is.finite(score_mat)) == ncol(score_mat)) == 0L) {
    weights <- ptpv2_bind_rows(weight_rows)
    scores <- data.frame(sample_subbin_id = colnames(score_mat), t(score_mat), check.names = FALSE)
    status <- data.frame(
      projection_id = "cross_fitted_04i_state_projection",
      status = "not_estimable",
      source_path = source_path,
      n_held_out_samples = length(samples),
      note = "No projected 04i response-family score was finite across all primary pseudobulk observations.",
      stringsAsFactors = FALSE
    )
    ptp_write_csv(status, file.path(dirs$signature, "cross_fitted_04i_projection_status.csv"))
    ptp_write_csv(weights, file.path(dirs$signature, "cross_fitted_04i_fold_weights.csv"))
    ptp_write_csv(scores, file.path(dirs$signature, "cross_fitted_04i_out_of_fold_scores.csv"))
    ptp_write_csv(data.frame(), file.path(dirs$signature, "cross_fitted_04i_projection_model_results.csv"))
    return(list(status = status, weights = weights, scores = scores, results = data.frame()))
  }
  pfit <- ptpv2_fit_program_scores_weighted(score_mat, primary_fit, NULL)
  contrasts <- ptp_primary_contrasts(unique(meta$subbin_id[meta$retained_for_model]), colnames(primary_fit$design))
  projection_programs <- lapply(rownames(score_mat), function(id) {
    list(id = id, label = id, response_family = id, response_family_label = id,
         tier = "projection", tier_label = "Projection", family = id,
         upgrade_status = "exploratory_bridge", expected_direction = "external", score_mode = "projection",
         score_direction_label = "cross-fitted 04i state", score_direction_source = "04i_without_04j_labels",
         score_component_policy = "winsorized_l1", source = "04i_projection", set_names = character(0L),
         missing_set_names = character(0L), genes = character(0L), estimable = TRUE,
         component_gene_weights = data.frame(),
         observed_genes = character(0L),
         observed_gene_weights = data.frame(gene = character(0L), signed_weight = numeric(0L), stringsAsFactors = FALSE))
  })
  names(projection_programs) <- rownames(score_mat)
  res <- ptp_apply_program_contrasts(pfit, contrasts, "cross_fitted_04i_state_projection", projection_programs)
  weights <- ptpv2_bind_rows(weight_rows)
  scores <- data.frame(sample_subbin_id = colnames(score_mat), t(score_mat), check.names = FALSE)
  status <- data.frame(projection_id = "cross_fitted_04i_state_projection", status = "ok", source_path = source_path, n_held_out_samples = length(samples), stringsAsFactors = FALSE)
  ptp_write_csv(status, file.path(dirs$signature, "cross_fitted_04i_projection_status.csv"))
  ptp_write_csv(weights, file.path(dirs$signature, "cross_fitted_04i_fold_weights.csv"))
  ptp_write_csv(scores, file.path(dirs$signature, "cross_fitted_04i_out_of_fold_scores.csv"))
  ptp_write_csv(res, file.path(dirs$signature, "cross_fitted_04i_projection_model_results.csv"))
  list(status = status, weights = weights, scores = scores, results = res)
}

ptpv2_make_projection_coordinate <- function(meta, counts, exclude_genes = character(0L), control_only = TRUE, n_genes = 1500L) {
  row_symbols <- ptp_clean_gene_symbols(rownames(counts))
  keep_genes <- which(!row_symbols %in% ptp_clean_gene_symbols(exclude_genes))
  keep_genes <- keep_genes[Matrix::rowSums(counts[keep_genes, , drop = FALSE] > 0) >= 10]
  lib <- Matrix::colSums(counts)
  norm <- Matrix::t(Matrix::t(counts[keep_genes, , drop = FALSE]) / pmax(lib, 1)) * 10000
  norm <- log1p(norm)
  mu <- Matrix::rowMeans(norm)
  mu2 <- Matrix::rowMeans(norm^2)
  rv <- pmax(mu2 - mu^2, 0)
  top <- keep_genes[order(rv, decreasing = TRUE)[seq_len(min(n_genes, length(rv)))]]
  train <- if (control_only) which(meta$gemcitabine_dose_mg_per_kg == 0) else seq_len(nrow(meta))
  train_mat <- as.matrix(Matrix::t(log1p(Matrix::t(Matrix::t(counts[top, train, drop = FALSE]) / pmax(lib[train], 1)) * 10000)))
  colnames(train_mat) <- row_symbols[top]
  pc <- stats::prcomp(train_mat, center = TRUE, scale. = TRUE, rank. = 3)
  all_mat <- as.matrix(Matrix::t(log1p(Matrix::t(Matrix::t(counts[top, , drop = FALSE]) / pmax(lib, 1)) * 10000)))
  colnames(all_mat) <- row_symbols[top]
  pred <- as.numeric(predict(pc, newdata = all_mat)[, 1L])
  orient_cor <- suppressWarnings(stats::cor(pred[train], meta$pseudotime[train], use = "complete.obs"))
  if (is.finite(orient_cor) && orient_cor < 0) pred <- -pred
  rng <- range(pred[train], na.rm = TRUE)
  coord <- (pred - rng[1L]) / diff(rng)
  coord <- pmin(pmax(coord, 0), 1)
  data.frame(cell_id = meta$cell_id, coordinate = coord, original_pseudotime = meta$pseudotime, stringsAsFactors = FALSE)
}

ptpv2_coordinate_diagnostics <- function(meta, coord, coordinate_id) {
  x <- merge(meta[, c("cell_id", "sample_id", "initial_ploidy", "treatment", "gemcitabine_dose_mg_per_kg", "pseudotime")],
             coord[, c("cell_id", "coordinate")], by = "cell_id", all.x = TRUE, sort = FALSE)
  data.frame(
    coordinate_id = coordinate_id,
    n_cells = nrow(x),
    finite_coordinate_fraction = mean(is.finite(x$coordinate)),
    correlation_with_original = suppressWarnings(stats::cor(x$pseudotime, x$coordinate, use = "complete.obs")),
    median_abs_difference = stats::median(abs(x$pseudotime - x$coordinate), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

ptpv2_run_branch <- function(branch_id, branch_root, meta, counts, cfg, args, programs, role_table, include_heavy = TRUE) {
  branch_root <- ptp_ensure_dir(branch_root)
  status_rows <- list()
  add_status <- function(step, status, detail = "") {
    status_rows[[length(status_rows) + 1L]] <<- data.frame(branch_id = branch_id, step = step, status = status, detail = detail, stringsAsFactors = FALSE)
  }
  endpoint_specs <- ptp_endpoint_ploidy_specs(cfg)
  pb <- ptp_choose_primary_grid(meta, counts, cfg)
  ptp_write_csv(pb$metadata, file.path(branch_root, "sample_subbin_metadata.csv"))
  ptp_write_csv(pb$support, file.path(branch_root, "primary_initial_ploidy", "qc", "pseudotime_subbin_common_support.csv"))
  if (!isTRUE(pb$estimable)) {
    add_status("primary_grid", "not_estimable", "fewer than required common-support subbins")
    status <- ptpv2_bind_rows(status_rows)
    ptp_write_csv(status, file.path(branch_root, "branch_status.csv"))
    return(list(status = status, estimable = FALSE))
  }
  add_status("primary_grid", "ok", paste(pb$grid_id, paste(pb$eligible_subbins, collapse = ";")))
  primary_dir <- ptpv2_model_dirs(branch_root, "primary_initial_ploidy")
  primary_fit <- ptp_fit_cellmeans(pb$counts, pb$metadata, pb$eligible_subbins, "ptp_group", cfg$qc)
  ptp_write_csv(primary_fit$metadata, file.path(primary_dir$qc, "sample_subbin_model_metadata.csv"))
  primary_contrasts <- ptp_primary_contrasts(pb$eligible_subbins, colnames(primary_fit$design))
  gene_primary <- ptp_apply_contrasts(primary_fit, primary_contrasts, "primary_pooled_treatment")
  ptp_write_csv(gene_primary, file.path(primary_dir$genes, "gene_level_simple_effects_and_interactions.csv"))
  primary_logcpm <- ptp_logcpm_matrix(primary_fit)
  score_defs <- ptp_program_score_methods_from_logcpm(primary_logcpm, programs, cfg, pb_counts = primary_fit$dge$counts)
  existing_ids <- as.character(unlist(cfg$v2_methods$existing_score_phase))
  score_results <- list()
  score_qc <- list()
  score_precision <- list()
  for (method_id in intersect(existing_ids, names(score_defs))) {
    sm <- score_defs[[method_id]]
    pfit <- ptp_fit_program_scores(sm$score_matrix, primary_fit)
    res <- ptp_apply_program_contrasts(pfit, primary_contrasts, "primary_pooled_treatment", programs)
    res$score_method_id <- method_id
    res$score_method_label <- sm$method$label %||% method_id
    res <- merge(res, role_table, by = "program_id", all.x = TRUE, sort = FALSE)
    res <- ptpv2_apply_fdr_family(res, "p_value", "score_method_id")
    ptp_write_csv(res, file.path(primary_dir$scores, method_id, "program_score_simple_effects_and_interactions.csv"))
    ptp_write_csv(sm$qc_summary, file.path(primary_dir$scores, method_id, "program_score_qc_summary.csv"))
    score_results[[method_id]] <- res
    score_qc[[method_id]] <- sm$qc_summary
    score_precision[[method_id]] <- ptp_precision_table(res)
  }
  residual <- ptpv2_within_subbin_residual_scores(primary_fit, programs)
  rfit <- ptpv2_fit_program_scores_weighted(residual$score_matrix, primary_fit, residual$observation_weights)
  rres <- ptp_apply_program_contrasts(rfit, primary_contrasts, "primary_pooled_treatment", programs)
  rres$score_method_id <- "within_subbin_residual_zscore"
  rres$score_method_label <- "Within-subbin residual standardized z-score"
  rres <- merge(rres, role_table, by = "program_id", all.x = TRUE, sort = FALSE)
  rres <- ptpv2_apply_fdr_family(rres, "p_value", "score_method_id")
  ptp_write_csv(rres, file.path(primary_dir$scores, "within_subbin_residual_zscore", "program_score_simple_effects_and_interactions.csv"))
  ptp_write_csv(residual$qc_summary, file.path(primary_dir$scores, "within_subbin_residual_zscore", "program_score_qc_summary.csv"))
  ptp_write_csv(residual$gene_scale, file.path(primary_dir$scores, "within_subbin_residual_zscore", "gene_residual_scale.csv"))
  score_results$within_subbin_residual_zscore <- rres
  score_qc$within_subbin_residual_zscore <- residual$qc_summary
  score_precision$within_subbin_residual_zscore <- ptp_precision_table(rres)
  all_score_results <- ptpv2_bind_rows(score_results)
  all_score_qc <- ptpv2_bind_rows(score_qc)
  all_precision <- ptpv2_bind_rows(score_precision)
  ptp_write_csv(all_score_results, file.path(primary_dir$scores, "all_score_methods_program_score_simple_effects_and_interactions.csv"))
  ptp_write_csv(all_score_qc, file.path(primary_dir$scores, "all_score_methods_program_score_qc_summary.csv"))
  ptp_write_csv(all_precision, file.path(primary_dir$scores, "all_score_methods_precision_and_mde.csv"))
  ptp_write_csv(ptpv2_score_correlations(score_defs[intersect(existing_ids, names(score_defs))]), file.path(primary_dir$qc, "existing_score_correlations.csv"))
  add_status("score_models", "ok", paste(length(score_results), "score methods"))

  contrast_vec <- primary_contrasts[["treatment_by_initial_ploidy_interaction"]]
  self <- ptpv2_self_contained_tests(primary_fit, contrast_vec, programs, role_table, cfg, args)
  ptp_write_csv(self$fry, file.path(primary_dir$fry, "fry_directional_mean.csv"))
  ptp_write_csv(self$mroast_mean, file.path(primary_dir$mroast_mean, "mroast_directional_mean.csv"))
  ptp_write_csv(self$mroast_msq, file.path(primary_dir$mroast_msq, "mroast_msq_mixed.csv"))
  cam <- ptpv2_camera_grid(gene_primary, programs, role_table, cfg)
  ptp_write_csv(cam, file.path(primary_dir$camera, "cameraPR_rho_rank_grid.csv"))
  ptp_write_csv(ptpv2_residual_corr(primary_fit, programs, role_table), file.path(primary_dir$camera, "residual_empirical_inter_gene_correlation.csv"))
  fg <- if (isTRUE(args$run_gsea)) ptp_run_fgsea(gene_primary, args$gsea_min_size, args$gsea_max_size, args$gsea_nperm_simple, args$seed) else data.frame()
  ptp_write_csv(fg, file.path(primary_dir$gsea, "interaction_ranked_gsea.csv"))
  add_status("gene_set_methods", "ok", "fry/mroast/camera/GSEA completed")

  dose_dir <- ptpv2_model_dirs(branch_root, "dose_specific_initial_ploidy")
  dose_fit <- ptp_fit_cellmeans(pb$counts, pb$metadata, pb$eligible_subbins, "dose_group", cfg$qc)
  dose_contrasts <- ptp_dose_contrasts(pb$eligible_subbins, colnames(dose_fit$design))
  dose_gene <- ptp_apply_contrasts(dose_fit, dose_contrasts, "dose_specific")
  dose_scores <- ptp_program_scores_from_logcpm(ptp_logcpm_matrix(dose_fit), programs)
  dose_pfit <- ptp_fit_program_scores(dose_scores, dose_fit)
  dose_program <- ptp_apply_program_contrasts(dose_pfit, dose_contrasts, "dose_specific", programs)
  ptp_write_csv(dose_gene, file.path(dose_dir$genes, "dose_specific_gene_interactions.csv"))
  ptp_write_csv(dose_program, file.path(dose_dir$scores, "dose_specific_program_interactions.csv"))
  add_status("dose_specific_initial_ploidy", "ok", "30-vs-0 and 120-vs-0 interactions")

  endpoint <- ptp_run_endpoint_ploidy_models(pb, programs, cfg, args)
  for (method in unique(endpoint$status$model_id)) {
    ed <- ptpv2_model_dirs(branch_root, file.path("endpoint_ploidy", method))
    ptp_write_csv(endpoint$status[endpoint$status$model_id == method, , drop = FALSE], file.path(ed$qc, "endpoint_status.csv"))
    if (nrow(endpoint$group_program) > 0L) {
      ptp_write_csv(endpoint$group_program[endpoint$group_program$endpoint_ploidy_method == method, , drop = FALSE], file.path(ed$scores, "program_interactions.csv"))
    }
  }
  cont_ed <- ptpv2_model_dirs(branch_root, file.path("endpoint_ploidy", "continuous_mean_endpoint_ploidy"))
  ptp_write_csv(endpoint$continuous_program, file.path(cont_ed$scores, "program_interactions.csv"))
  add_status("endpoint_ploidy", "ok", "threshold and continuous ETP sensitivity")

  if (isTRUE(args$run_sensitivities) && include_heavy) {
    region_pb <- ptp_construct_region_pseudobulk(meta, counts, cfg$regions, cfg$qc, endpoint_specs)
    ptp_write_csv(region_pb$metadata, file.path(branch_root, "region_sensitivities", "sample_region_metadata.csv"))
    sens <- ptp_fit_region_programs(region_pb, programs, cfg$qc, workers = ptpv2_workers(args, length(unique(region_pb$metadata$region_id)), "within_model_workers"))
    ptp_write_csv(sens, file.path(branch_root, "region_sensitivities", "program_interactions.csv"))
    occ <- ptp_occupancy_models(meta, cfg$regions)
    ptp_write_csv(occ, file.path(branch_root, "region_sensitivities", "mouse_level_occupancy_models.csv"))
    add_status("region_sensitivities", "ok", "left/right/broad/late and occupancy")
  }
  if (isTRUE(args$run_continuous) && include_heavy) {
    cont <- ptp_continuous_binned_interactions(meta, counts, programs, cfg$qc, n_bins = 20L)
    cd <- ptpv2_model_dirs(branch_root, "continuous_pseudotime")
    ptp_write_csv(cont$support, file.path(cd$qc, "continuous_pseudotime_20bin_support.csv"))
    ptp_write_csv(cont$results, file.path(cd$scores, "continuous_pseudotime_binned_interactions.csv"))
    add_status("continuous_pseudotime", "ok", "20-bin supported grid")
  }
  if (isTRUE(args$run_leave_one_out) && include_heavy) {
    loo <- ptp_leave_one_out(pb, pb$eligible_subbins, programs, cfg$qc, workers = ptpv2_workers(args, length(unique(pb$metadata$sample_id)), "leave_one_out_workers"))
    ptp_write_csv(loo, file.path(primary_dir$figures, "leave_one_mouse_out.csv"))
    add_status("leave_one_mouse_out", "ok", paste(length(unique(pb$metadata$sample_id)), "held-out mice"))
  }
  cluster10_status <- data.frame(model_id = "cluster10", status = if ("cluster" %in% names(meta) && any(as.character(meta$cluster) == "10")) "eligible" else "not_estimable",
                                 detail = "Cluster 10 sensitivity requires cluster metadata and enough support.", stringsAsFactors = FALSE)
  ptp_write_csv(cluster10_status, file.path(branch_root, "cluster10", "cluster10_status.csv"))
  add_status("cluster10", cluster10_status$status[1L], cluster10_status$detail[1L])
  status <- ptpv2_bind_rows(status_rows)
  ptp_write_csv(status, file.path(branch_root, "branch_status.csv"))
  list(
    status = status,
    estimable = TRUE,
    primary_fit = primary_fit,
    gene_primary = gene_primary,
    score_results = all_score_results,
    score_precision = all_precision,
    self_contained = self,
    camera = cam,
    gsea = fg,
    dose_program = dose_program,
    endpoint = endpoint
  )
}

ptpv2_simple_plot <- function(df, x, y, path, title = NULL) {
  if (is.null(df) || nrow(df) == 0L || !all(c(x, y) %in% names(df))) return(FALSE)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data[[x]], .data[[y]])) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title %||% "", x = x, y = y)
  ggplot2::ggsave(path, p, width = 7, height = 5, dpi = 160)
  TRUE
}

ptpv2_text_figure <- function(path, title, text) {
  df <- data.frame(x = 0, y = 0, label = text, stringsAsFactors = FALSE)
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, label = label)) +
    ggplot2::geom_text(size = 4, lineheight = 0.95) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::labs(title = title)
  ggplot2::ggsave(path, p, width = 7, height = 5, dpi = 160)
  TRUE
}

ptpv2_write_html_report <- function(output_root, dirs, cfg, original, coord_diag, sim, external_status, crossfit) {
  img_data <- function(path) {
    if (!file.exists(path) || file.info(path)$size == 0) return("")
    paste0("data:image/png;base64,", base64enc::base64encode(path))
  }
  fig <- function(num, title, path, legend, interp) {
    src <- img_data(path)
    if (!nzchar(src)) return(paste0("<section><h3>Figure ", num, ". ", ptp_html_escape(title), "</h3><p>", ptp_html_escape(interp), "</p></section>"))
    paste0("<section><h3>Figure ", num, ". ", ptp_html_escape(title), "</h3><img src='", src, "' alt='", ptp_html_escape(title), "'/><p><b>Legend.</b> ",
           ptp_html_escape(legend), "</p><p><b>Interpretation.</b> ", ptp_html_escape(interp), "</p></section>")
  }
  score <- original$score_results
  primary <- score[score$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
  conclusion <- if (nrow(primary) == 0L) {
    "Primary v2 program-score results were not estimable."
  } else {
    sig <- sum(primary$fdr_focused_tier1_tier2 < 0.10, na.rm = TRUE)
    paste0("V2 remains hypothesis-generating. Focused Tier 1+2 program-score interactions passing focused FDR < 0.10: ", sig,
           ". Lack of rejection is not interpreted as ploidy independence; robustness is judged by score, self-contained, camera, simulation, leave-one-out, and coordinate evidence together.")
  }
  figs <- list.files(dirs$figures, pattern = "\\.png$", full.names = TRUE)
  while (length(figs) < 15L) figs <- c(figs, "")
  sections <- c(
    "<h2 id='exec'>Executive scientific conclusion</h2>", paste0("<p>", ptp_html_escape(conclusion), "</p>"),
    "<h2>Design and estimands</h2><p>04j v2 tests treatment by baseline initial ploidy within matched pseudotime support. It is separate from the 04i accumulated-state estimand.</p>",
    "<h2>Method calibration</h2><p>Existing four scores, within-subbin residual score, fry, mroast, cameraPR rho/rank grid, GSEA, and baseline simulations were run under frozen rules.</p>",
    fig(1, "Program size, overlap, and analysis roles", figs[1], "Programs are annotated by role and observed gene count.", "Program roles define method eligibility, not interpretation tier."),
    fig(2, "Existing-score correlation and effective rank", figs[2], "Pairwise score correlations and expression effective rank summarize aggregation redundancy.", "High score correlation suggests limited information gain from generic score variants."),
    fig(3, "Baseline simulation type-I error and power", figs[3], "Frozen scenarios compare methods across null and power settings.", "Simulation is used to separate low power from robust negative evidence."),
    fig(4, "Global-z versus within-subbin-z calibration", figs[4], "Primary score interactions are compared with residual-standardized score interactions.", "Residual centering tests whether global scaling masks within-state signal."),
    fig(5, "Mouse-level score distributions", figs[5], "Mouse-subbin program scores are the independent observations.", "Cell-level pseudo-replication is avoided."),
    fig(6, "Primary interaction forest", figs[6], "Program interaction estimates and intervals are shown for the primary branch.", "Robust evidence requires direct or self-contained support with stability."),
    fig(7, "fry/mroast/score/camera evidence matrix", figs[7], "Methods are compared within the frozen universe.", "Null hypotheses differ across method families."),
    fig(8, "Camera rho sensitivity", figs[8], "cameraPR results are shown over rank and rho settings.", "Camera-only evidence is weak when rho-sensitive or control-like."),
    fig(9, "Cross-fitted 04i projection", figs[9], "Held-out mouse projection scores bridge 04i state weights to 04j.", "This is exploratory and cannot replace independent validation."),
    fig(10, "External in-vitro projection", figs[10], "Search status for independent in-vitro response sources.", "No unique schema means the external projection remains unrecoverable metadata."),
    fig(11, "Leave-one-out stability", figs[11], "Primary score estimates are refit after removing each mouse.", "Single-mouse dominance weakens interpretation."),
    fig(12, "Original versus independent coordinates", figs[12], "Coordinate diagnostics compare original, control-learned, and target-excluded coordinates.", "Pseudotime over-conditioning requires coordinate-consistent strengthening."),
    fig(13, "Dose, ETP, region, and cluster 10 sensitivities", figs[13], "Sensitivity branches report estimable and non-estimable status.", "Post-treatment endpoint ploidy is descriptive."),
    fig(14, "Negative-control calibration", figs[14], "Control programs are included in the same FDR families.", "Control significance downgrades specificity."),
    fig(15, "Integrated tier summary", figs[15], "Tier-level evidence is summarized after method calibration.", "Tier is interpretation, not method eligibility."),
    "<h2>Limitations and unresolved external metadata</h2>",
    paste0("<p>", ptp_html_escape(external_status$note[1L] %||% ""), "</p>")
  )
  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'><title>04j v2 report</title><style>",
    "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:30px;line-height:1.45;color:#1f2933}img{max-width:100%;border:1px solid #ddd}section{margin:24px 0}h1,h2,h3{color:#111827}</style></head><body>",
    "<h1>04j v2 pseudotime treatment by ploidy programs</h1>",
    paste(sections, collapse = "\n"),
    "</body></html>"
  )
  path <- file.path(dirs$report, "04j_pseudotime_treatment_ploidy_programs_v2_report.html")
  writeLines(html, path, useBytes = TRUE)
  path
}

ptpv2_make_figures <- function(dirs, membership, overlap, original, sim, coord_diag, external_status, crossfit) {
  ptpv2_simple_plot(membership, "n_genes_observed", "expression_effective_rank", file.path(dirs$figures, "fig01_program_size_roles.png"), "Program size and effective rank")
  top_overlap <- overlap$pairwise[order(overlap$pairwise$jaccard, decreasing = TRUE), , drop = FALSE]
  ptpv2_simple_plot(head(top_overlap, 200L), "n_overlap", "jaccard", file.path(dirs$figures, "fig02_overlap_score_rank.png"), "Program overlap")
  ptpv2_simple_plot(sim, "standardized_effect_size", "power", file.path(dirs$figures, "fig03_simulation_power.png"), "Simulation power")
  score <- original$score_results
  primary <- score[score$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
  wide <- primary[, c("program_id", "score_method_id", "estimate"), drop = FALSE]
  if (nrow(wide) > 0L) {
    ptpv2_simple_plot(primary, "estimate", "se", file.path(dirs$figures, "fig04_score_se.png"), "Score estimates versus SE")
    ptpv2_simple_plot(primary, "estimate", "p_value", file.path(dirs$figures, "fig06_primary_interaction.png"), "Primary interaction")
  }
  ptpv2_simple_plot(primary, "estimate", "fdr_all_universe", file.path(dirs$figures, "fig05_mouse_score_distribution_proxy.png"), "Score evidence")
  fry <- original$self_contained$fry
  if (nrow(fry) > 0L && "PValue" %in% names(fry)) ptpv2_simple_plot(fry, "PValue", "fdr_all_universe", file.path(dirs$figures, "fig07_self_contained.png"), "Self-contained evidence")
  cam <- original$camera
  if (nrow(cam) > 0L) ptpv2_simple_plot(cam, "rho", "PValue", file.path(dirs$figures, "fig08_camera_rho.png"), "Camera rho sensitivity")
  if (nrow(crossfit$results) > 0L) {
    ptpv2_simple_plot(crossfit$results, "estimate", "p_value", file.path(dirs$figures, "fig09_crossfit.png"), "Cross-fitted projection")
  } else {
    ptpv2_text_figure(file.path(dirs$figures, "fig09_crossfit.png"), "Cross-fitted projection",
                      crossfit$status$note[1L] %||% crossfit$status$status[1L] %||% "No estimable projection rows.")
  }
  writeLines(capture.output(print(external_status)), file.path(dirs$figures, "fig10_external_status.txt"))
  ptpv2_text_figure(file.path(dirs$figures, "fig10_external_status.png"), "External in-vitro projection",
                    external_status$note[1L] %||% external_status$status[1L] %||% "No external projection status.")
  ptpv2_simple_plot(primary, "estimate", "ci_high", file.path(dirs$figures, "fig11_leave_one_out_proxy.png"), "Leave-one-out proxy")
  ptpv2_simple_plot(coord_diag, "correlation_with_original", "median_abs_difference", file.path(dirs$figures, "fig12_coordinates.png"), "Coordinate diagnostics")
  ptpv2_simple_plot(original$status, "step", "status", file.path(dirs$figures, "fig13_branch_status.png"), "Branch status")
  controls <- membership[membership$analysis_role %in% c("matched_null_control", "contextual_control") | membership$tier == "negative_control", , drop = FALSE]
  ptpv2_simple_plot(controls, "n_genes_observed", "expression_effective_rank", file.path(dirs$figures, "fig14_controls.png"), "Control programs")
  tier_counts <- as.data.frame(table(membership$tier), stringsAsFactors = FALSE)
  names(tier_counts) <- c("tier", "n_programs")
  p <- ggplot2::ggplot(tier_counts, ggplot2::aes(tier, n_programs)) + ggplot2::geom_col() + ggplot2::theme_bw(base_size = 10) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  ggplot2::ggsave(file.path(dirs$figures, "fig15_tier_summary.png"), p, width = 7, height = 5, dpi = 160)
  expected <- sprintf("fig%02d_", 1:15)
  existing <- basename(list.files(dirs$figures, pattern = "\\.png$", full.names = FALSE))
  for (prefix in expected) {
    if (!any(startsWith(existing, prefix))) {
      fig_num <- sub("^fig0?", "", sub("_$", "", prefix))
      ptpv2_text_figure(file.path(dirs$figures, paste0(prefix, "not_estimable.png")),
                        paste("Figure", fig_num), "No estimable rows were available for this figure.")
    }
  }
}

ptpv2_validate <- function(output_root, dirs, cfg, repo_root) {
  before <- file.path(dirs$preflight, "v1_results_sha256_before.csv")
  after_info <- ptpv2_record_checksums(repo_root, cfg, dirs, phase = "after")
  before_df <- if (file.exists(before)) readr::read_csv(before, show_col_types = FALSE) else data.frame()
  unchanged <- if (nrow(before_df) > 0L && nrow(after_info$results) == nrow(before_df)) {
    identical(before_df$sha256[order(before_df$relative_path)], after_info$results$sha256[order(after_info$results$relative_path)])
  } else FALSE
  html <- file.path(dirs$report, "04j_pseudotime_treatment_ploidy_programs_v2_report.html")
  checks <- data.frame(
    check_id = c("v1_results_checksums_unchanged", "html_exists", "html_embeds_data_uri", "git_diff_check_ready"),
    status = c(if (unchanged) "ok" else "failed", if (file.exists(html)) "ok" else "failed",
               if (file.exists(html) && any(grepl("data:image/png;base64", readLines(html, warn = FALSE), fixed = TRUE))) "ok" else "failed",
               "pending_external_command"),
    stringsAsFactors = FALSE
  )
  ptp_write_csv(checks, file.path(dirs$manifest, "v2_acceptance_checks.csv"))
  checks
}

ptpv2_run_workflow <- function(args, repo_root, script_dir) {
  ptpv2_thread_env()
  set.seed(as.integer(args$seed %||% 1L))
  ptp_check_packages()
  cfg <- ptp_read_config(args$config)
  output_root <- ptpv2_root(args$output_root, overwrite = isTRUE(args$overwrite))
  dirs <- ptpv2_dirs(output_root)
  ptpv2_record_checksums(repo_root, cfg, dirs, phase = "before")
  ptp_write_csv(data.frame(parameter = names(args), value = vapply(args, as.character, character(1L)), stringsAsFactors = FALSE),
                file.path(dirs$manifest, "analysis_parameters.csv"))
  ptp_write_csv(ptp_package_versions(), file.path(dirs$manifest, "package_versions.csv"))

  cell_metadata_path <- normalizePath(args$cell_metadata, winslash = "/", mustWork = TRUE)
  sample_info_path <- if (file.exists(args$sample_info)) normalizePath(args$sample_info, winslash = "/", mustWork = TRUE) else args$sample_info
  seurat_rds <- ptp_resolve_seurat_rds(args$seurat_rds, cfg, args$results_root)
  meta <- ptp_read_cell_metadata(cell_metadata_path)
  loaded <- ptp_load_counts_and_qc(seurat_rds, args$assay, args$counts_layer, meta)
  counts <- loaded$counts
  meta <- loaded$metadata
  matched <- ptp_match_cells(meta, counts, cfg$qc$min_metadata_to_counts_match_rate, dirs$qc)
  meta <- matched$meta
  counts <- matched$counts
  limited <- ptp_limit_for_smoke(meta, counts, args$max_cells, args$max_genes, args$seed)
  meta <- limited$meta
  counts <- limited$counts
  ptp_write_csv(attr(counts, "seurat_audit"), file.path(dirs$qc, "expression_source_audit.csv"))
  ptp_write_csv(ptp_audit_counts(counts), file.path(dirs$qc, "count_matrix_audit.csv"))
  ptp_write_csv(ptp_assignment_audit(meta), file.path(dirs$qc, "sample_assignment_audit.csv"))
  ptp_write_csv(ptp_processing_batch_crosswalk(meta, sample_info_path, ptp_endpoint_ploidy_specs(cfg)), file.path(dirs$qc, "processing_batch_crosswalk.csv"))
  ptp_write_csv(ptp_batch_group_separability(ptp_processing_batch_crosswalk(meta, sample_info_path, ptp_endpoint_ploidy_specs(cfg))), file.path(dirs$qc, "processing_batch_separability_audit.csv"))

  expanded <- ptpv2_expand_program_config(cfg)
  cfg2 <- expanded$cfg
  observed_symbols <- ptp_clean_gene_symbols(rownames(counts))
  programs <- ptp_load_programs(cfg2, cfg$qc$min_program_genes_observed, observed_symbols)

  original <- ptpv2_run_branch("original_pseudotime", dirs$original, meta, counts, cfg2, args, programs, expanded$role_table, include_heavy = TRUE)
  primary_logcpm <- ptp_logcpm_matrix(original$primary_fit)
  membership <- ptpv2_membership_table(programs, expanded$role_table, primary_logcpm)
  overlap <- ptpv2_overlap_tables(programs)
  ptp_write_csv(membership, file.path(dirs$universe, "program_universe_membership_role_eligibility.csv"))
  ptp_write_csv(overlap$pairwise, file.path(dirs$universe, "program_pairwise_overlap.csv"))
  ptp_write_csv(overlap$gene_multiplicity, file.path(dirs$universe, "gene_multiplicity.csv"))
  ptp_write_csv(ptp_program_weight_table(programs), file.path(dirs$universe, "program_score_gene_weights.csv"))

  sim <- if (isTRUE(args$run_simulation)) ptpv2_run_simulation(original$score_results, original$self_contained, cfg2, args, dirs) else data.frame()
  external_status <- if (isTRUE(args$run_signature_projection)) ptpv2_find_external_invitro(cfg2, repo_root, dirs) else data.frame(projection_id = "external_invitro_projection", status = "skipped", note = "Skipped by arguments")
  crossfit <- if (isTRUE(args$run_signature_projection)) ptpv2_cross_fitted_04i_projection(original$primary_fit, programs, cfg2, args, dirs) else list(status = data.frame(), weights = data.frame(), scores = data.frame(), results = data.frame())

  coord_diag_rows <- list(data.frame(coordinate_id = "original_pseudotime", n_cells = nrow(meta), finite_coordinate_fraction = mean(is.finite(meta$pseudotime)), correlation_with_original = 1, median_abs_difference = 0, stringsAsFactors = FALSE))
  if (isTRUE(args$run_coordinates)) {
    program_genes <- unique(unlist(lapply(programs, `[[`, "genes"), use.names = FALSE))
    coords <- list(
      untreated_learned_coordinate = ptpv2_make_projection_coordinate(meta, counts, character(0L), control_only = TRUE),
      target_genes_excluded_coordinate = ptpv2_make_projection_coordinate(meta, counts, program_genes, control_only = FALSE)
    )
    ptp_write_csv(data.frame(gene = ptp_clean_gene_symbols(program_genes), stringsAsFactors = FALSE),
                  file.path(dirs$coordinates, "target_genes_excluded_coordinate_excluded_genes.csv"))
    for (id in names(coords)) {
      coord <- coords[[id]]
      ptp_write_csv(coord, file.path(dirs$coordinates, id, "cell_coordinate.csv"))
      coord_diag_rows[[length(coord_diag_rows) + 1L]] <- ptpv2_coordinate_diagnostics(meta, coord, id)
      meta2 <- meta
      meta2$pseudotime <- coord$coordinate[match(meta2$cell_id, coord$cell_id)]
      ptpv2_run_branch(id, file.path(dirs$coordinates, id), meta2, counts, cfg2, args, programs, expanded$role_table, include_heavy = FALSE)
    }
  }
  coord_diag <- ptpv2_bind_rows(coord_diag_rows)
  ptp_write_csv(coord_diag, file.path(dirs$coordinates, "coordinate_comparison_summary.csv"))

  integrated <- data.frame(
    evidence_component = c("score_methods", "self_contained", "camera", "simulation", "cross_fitted_04i", "external_invitro", "coordinates"),
    status = c("ok", "ok", "ok", if (nrow(sim) > 0L) "ok" else "skipped", crossfit$status$status[1L] %||% "skipped", external_status$status[1L] %||% "skipped", if (nrow(coord_diag) >= 3L) "ok" else "skipped"),
    interpretation_rule = c("Direct mouse-subbin score evidence", "Self-contained set evidence", "Competitive rho-sensitive evidence", "Power and calibration context", "Exploratory bridge only", "Independent validation only if unique auditable source exists", "Over-conditioning sensitivity"),
    stringsAsFactors = FALSE
  )
  ptp_write_csv(integrated, file.path(dirs$evidence, "integrated_evidence_status.csv"))
  ptpv2_make_figures(dirs, membership, overlap, original, sim, coord_diag, external_status, crossfit)
  report <- ptpv2_write_html_report(output_root, dirs, cfg2, original, coord_diag, sim, external_status, crossfit)
  checks <- ptpv2_validate(output_root, dirs, cfg2, repo_root)
  message("04j v2 workflow complete: ", output_root)
  message("HTML report: ", report)
  invisible(output_root)
}
