.ptpv3_source_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) "")
.ptpv3_script_dir <- if (nzchar(.ptpv3_source_file)) dirname(.ptpv3_source_file) else file.path(getwd(), "Code/in-vivo")

source(file.path(.ptpv3_script_dir, "04j_pseudotime_treatment_ploidy_programs_v2_util.R"))
source(file.path(.ptpv3_script_dir, "04i_pseudotime_state_pathways_util.R"))

ptpv3_thread_env <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
}

ptpv3_bind_rows <- function(rows) {
  rows <- rows[!vapply(rows, function(x) is.null(x) || nrow(as.data.frame(x)) == 0L, logical(1L))]
  if (length(rows) == 0L) return(data.frame())
  rows <- lapply(rows, as.data.frame, stringsAsFactors = FALSE)
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    for (nm in setdiff(all_names, names(x))) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, rows)
}

ptpv3_read_config <- function(path) {
  v3 <- yaml::read_yaml(path)
  base_path <- v3$analysis$base_v2_config
  if (is.null(base_path) || !file.exists(base_path)) stop("V3 config requires existing analysis.base_v2_config.", call. = FALSE)
  cfg <- ptp_read_config(base_path)
  cfg$analysis$name <- v3$analysis$name
  cfg$analysis$version <- v3$analysis$version
  cfg$analysis$frozen_date <- v3$analysis$frozen_date
  cfg$analysis$post_hoc_status <- v3$analysis$post_hoc_status
  cfg$analysis$plan_path <- v3$analysis$plan_path
  cfg$analysis$v1_result_root <- v3$analysis$v1_result_root
  cfg$analysis$v2_result_root <- v3$analysis$v2_result_root
  cfg$analysis$v3_result_root <- v3$analysis$v3_result_root
  cfg$analysis$guardrail <- v3$analysis$guardrail
  cfg$v3 <- v3
  if (!is.null(v3$simulation)) cfg$simulation <- v3$simulation
  if (!is.null(v3$self_contained_gene_set_tests$mroast_rotations)) {
    cfg$v2_methods$self_contained$mroast_rotations_observed <- as.integer(v3$self_contained_gene_set_tests$mroast_rotations)
  }
  cfg$v2_methods$self_contained$mroast_seed <- as.integer(v3$reproducibility$rotation_seed %||% 41043L)
  cfg
}

ptpv3_dirs <- function(output_root) {
  d <- function(...) ptp_ensure_dir(file.path(output_root, ...))
  list(
    root = ptp_ensure_dir(output_root),
    manifest = d("manifest"),
    frozen = d("frozen_definitions"),
    qc = d("qc"),
    primary = list(
      root = d("primary_initial_ploidy"),
      score_models = d("primary_initial_ploidy", "score_models"),
      exact = d("primary_initial_ploidy", "exact_permutation"),
      fry_mroast = d("primary_initial_ploidy", "fry_mroast"),
      camera = d("primary_initial_ploidy", "camera"),
      multidim = d("primary_initial_ploidy", "multidimensional"),
      trajectory = d("primary_initial_ploidy", "trajectory_global")
    ),
    secondary_etp = d("secondary_etp"),
    secondary_end_time_ploidy = d("secondary_end_time_ploidy"),
    secondary_dose = d("secondary_dose"),
    crossfit = d("crossfit_04i"),
    directional = d("directional_signatures"),
    simulation = d("simulation"),
    controls = d("controls"),
    figures = d("figures"),
    tables = d("tables"),
    report = d("report"),
    logs = d("logs"),
    checkpoints = d("checkpoints")
  )
}

ptpv3_checksum_tree <- function(path, label, phase, dirs) {
  if (!file.exists(path)) {
    out <- data.frame(label = label, phase = phase, path = path, relative_path = "", exists = FALSE,
                      is_dir = FALSE, bytes = NA_real_, sha256 = NA_character_, stringsAsFactors = FALSE)
    ptp_write_csv(out, file.path(dirs$manifest, paste0(label, "_sha256_", phase, ".csv")))
    return(out)
  }
  if (dir.exists(path)) {
    files <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    files <- files[file.info(files)$isdir == FALSE]
  } else {
    files <- path
  }
  files <- sort(files)
  root_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  file_norm <- normalizePath(files, winslash = "/", mustWork = FALSE)
  out <- data.frame(
    label = label,
    phase = phase,
    path = file_norm,
    relative_path = if (dir.exists(path)) sub(paste0("^", gsub("([\\W])", "\\\\\\1", root_norm), "/?"), "", file_norm) else basename(file_norm),
    exists = file.exists(files),
    is_dir = FALSE,
    bytes = suppressWarnings(file.info(files)$size),
    sha256 = vapply(files, ptp_file_checksum, character(1L)),
    stringsAsFactors = FALSE
  )
  ptp_write_csv(out, file.path(dirs$manifest, paste0(label, "_sha256_", phase, ".csv")))
  out
}

ptpv3_compare_checksums <- function(before, after, label, dirs) {
  b <- before[, c("relative_path", "sha256"), drop = FALSE]
  names(b)[2L] <- "sha256_before"
  a <- after[, c("relative_path", "sha256"), drop = FALSE]
  names(a)[2L] <- "sha256_after"
  x <- merge(b, a, by = "relative_path", all = TRUE, sort = FALSE)
  x$status <- ifelse(is.na(x$sha256_before), "added",
                     ifelse(is.na(x$sha256_after), "removed",
                            ifelse(x$sha256_before == x$sha256_after, "unchanged", "changed")))
  x$label <- label
  x <- x[, c("label", setdiff(names(x), "label")), drop = FALSE]
  ptp_write_csv(x, file.path(dirs$manifest, paste0(label, "_sha256_comparison.csv")))
  x
}

ptpv3_program_role_table <- function(programs) {
  x <- ptp_program_membership_table(programs)
  x$analysis_role <- ifelse(x$tier == "negative_control", "matched_null_control",
                            ifelse(x$program_id == "gemcitabine_handling_custom", "mechanistic_directional_signature",
                                   ifelse(grepl("^atomic_", x$program_id), "broad_pathway", "legacy_non_directional")))
  x$analysis_role[x$source == "custom" & x$tier != "negative_control" & x$program_id != "gemcitabine_handling_custom"] <- "legacy_non_directional"
  x$eligibility <- isTRUE(TRUE) & x$estimable & x$n_genes_observed >= 2L
  x$directional_eligibility <- x$analysis_role == "mechanistic_directional_signature"
  x$control_status <- x$tier == "negative_control" | x$analysis_role %in% c("matched_null_control", "contextual_control")
  x$scalar_score_eligible <- x$eligibility
  x$self_contained_eligible <- x$eligibility
  x$competitive_eligible <- x$eligibility
  x$gsea_mapping_eligible <- x$eligibility
  x$legacy_v1_comparison <- x$analysis_role == "legacy_non_directional"
  x$directional_claim_eligible <- x$directional_eligibility
  x$source_set_id <- x$set_names
  x$source_collection <- x$source
  x$msigdb_version <- "msigdbr_runtime"
  x
}

ptpv3_branch_role_table <- function(role_table) {
  base_cols <- c("program_label", "tier", "tier_label", "response_family", "response_family_label",
                 "family", "upgrade_status", "expected_direction", "score_mode",
                 "score_direction_label", "score_direction_source", "score_component_policy",
                 "n_genes_observed", "n_positive_weight_genes", "n_negative_weight_genes",
                 "sum_abs_observed_weights")
  role_table[, setdiff(names(role_table), base_cols), drop = FALSE]
}

ptpv3_join_program_meta <- function(df, role_table) {
  if (is.null(df) || nrow(df) == 0L || !"program_id" %in% names(df)) return(as.data.frame(df))
  meta_cols <- c("program_id", "program_label", "tier", "tier_label", "response_family",
                 "response_family_label", "family", "upgrade_status", "expected_direction",
                 "analysis_role", "eligibility", "directional_eligibility", "control_status",
                 "scalar_score_eligible", "self_contained_eligible", "competitive_eligible",
                 "gsea_mapping_eligible", "legacy_v1_comparison", "directional_claim_eligible",
                 "source_set_id", "source_collection", "msigdb_version")
  meta <- role_table[, intersect(meta_cols, names(role_table)), drop = FALSE]
  drop_cols <- setdiff(intersect(names(df), names(meta)), "program_id")
  df <- df[, setdiff(names(df), drop_cols), drop = FALSE]
  merge(df, meta, by = "program_id", all.x = TRUE, sort = FALSE)
}

ptpv3_apply_fdr_family <- function(df, p_col = "p_value", method_col = "method_id", role_table = NULL) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!is.null(role_table)) df <- ptpv3_join_program_meta(df, role_table)
  if (nrow(df) == 0L || !p_col %in% names(df)) return(df)
  if (!method_col %in% names(df)) df[[method_col]] <- "method"
  df$fdr_all_universe <- NA_real_
  df$fdr_focused_tier1_tier2 <- NA_real_
  df$fdr_tier3 <- NA_real_
  df$fdr_controls <- NA_real_
  for (method in unique(df[[method_col]])) {
    idx <- which(df[[method_col]] == method & is.finite(df[[p_col]]) & isTRUE(TRUE))
    if (!length(idx)) next
    df$fdr_all_universe[idx] <- stats::p.adjust(df[[p_col]][idx], "BH")
    focused <- idx[df$tier[idx] %in% c("tier1_core", "tier2_upgraded_focused")]
    tier3 <- idx[df$tier[idx] %in% "tier3_exploratory_only"]
    controls <- idx[(df$control_status[idx] %in% TRUE) | (df$tier[idx] %in% "negative_control")]
    if (length(focused)) df$fdr_focused_tier1_tier2[focused] <- stats::p.adjust(df[[p_col]][focused], "BH")
    if (length(tier3)) df$fdr_tier3[tier3] <- stats::p.adjust(df[[p_col]][tier3], "BH")
    if (length(controls)) df$fdr_controls[controls] <- stats::p.adjust(df[[p_col]][controls], "BH")
  }
  df
}

ptpv3_assert_complete_program_meta <- function(df, role_table, artifact, dirs) {
  if (is.null(df) || nrow(df) == 0L || !"program_id" %in% names(df)) return(TRUE)
  x <- ptpv3_join_program_meta(df[, "program_id", drop = FALSE], role_table)
  bad <- x[x$eligibility %in% TRUE & (!nzchar(x$tier %||% "") | is.na(x$tier) |
                                      !nzchar(x$response_family %||% "") | is.na(x$response_family)), , drop = FALSE]
  status <- data.frame(
    artifact = artifact,
    status = if (nrow(bad) == 0L) "ok" else "failed",
    n_rows = nrow(df),
    n_bad_eligible_programs = nrow(bad),
    stringsAsFactors = FALSE
  )
  ptp_write_csv(status, file.path(dirs$qc, paste0("metadata_acceptance_", ptp_safe_name(artifact), ".csv")))
  if (nrow(bad) > 0L) stop("Eligible programs missing tier/response metadata in ", artifact, call. = FALSE)
  TRUE
}

ptpv2_self_contained_tests <- function(model, contrast_vec, programs, role_table, cfg, args) {
  indices <- ptpv2_gene_set_indices(model, programs)
  keep_ids <- names(indices)[vapply(indices, length, integer(1L)) >= 2L]
  names(keep_ids) <- keep_ids
  nrot <- as.integer(cfg$v2_methods$self_contained$mroast_rotations_observed %||% 9999L)
  master_seed <- as.integer(cfg$v2_methods$self_contained$mroast_seed %||% args$seed %||% 1L)
  corr <- model$duplicate_correlation_second %||% model$duplicate_correlation %||% 0
  block <- model$metadata$sample_id
  run_one <- function(program_id) {
    i <- match(program_id, keep_ids)
    set.seed(master_seed + i)
    p <- programs[[program_id]]
    gw <- ptpv2_gene_weight_vector(model, p)
    idx <- list(indices[[program_id]])
    names(idx) <- program_id
    common_args <- list(y = model$voom, index = idx, design = model$design, contrast = contrast_vec,
                        gene.weights = gw, block = block, correlation = corr, sort = "none")
    fry <- tryCatch(do.call(limma::fry, common_args), error = function(e) data.frame(error = conditionMessage(e)))
    mean_args <- c(common_args, list(set.statistic = "mean", nrot = nrot))
    mr_mean <- tryCatch(do.call(limma::mroast, mean_args), error = function(e) data.frame(error = conditionMessage(e)))
    msq_args <- c(common_args, list(set.statistic = "msq", nrot = nrot))
    mr_msq <- tryCatch(do.call(limma::mroast, msq_args), error = function(e) data.frame(error = conditionMessage(e)))
    add_audit <- function(x, method) {
      x <- as.data.frame(x, stringsAsFactors = FALSE)
      x$block_column <- "sample_id"
      x$duplicate_correlation_second_pass <- corr
      x$n_rotations <- if (grepl("mroast", method)) nrot else NA_integer_
      x$rotation_seed <- if (grepl("mroast", method)) master_seed + i else NA_integer_
      x$program_id <- program_id
      x$method_id <- method
      x
    }
    list(fry = add_audit(fry, "fry"),
         mroast_mean = add_audit(mr_mean, "mroast_mean"),
         mroast_msq = add_audit(mr_msq, "mroast_msq"))
  }
  workers <- ptpv2_workers(args, length(keep_ids), "rotation_test_workers")
  res <- ptp_parallel_lapply(keep_ids, run_one, workers = workers, task_label = "v3_self_contained_gene_set_tests")$results
  bind_method <- function(method) ptpv3_bind_rows(lapply(res, `[[`, method))
  list(
    fry = ptpv3_apply_fdr_family(bind_method("fry"), p_col = if ("PValue" %in% names(bind_method("fry"))) "PValue" else "P.Value", method_col = "method_id", role_table = role_table),
    mroast_mean = ptpv3_apply_fdr_family(bind_method("mroast_mean"), p_col = "PValue", method_col = "method_id", role_table = role_table),
    mroast_msq = ptpv3_apply_fdr_family(bind_method("mroast_msq"), p_col = if ("PValue.Mixed" %in% names(bind_method("mroast_msq"))) "PValue.Mixed" else "PValue", method_col = "method_id", role_table = role_table)
  )
}

ptpv3_assignment_table <- function(model_meta, treated_per_ploidy = 4L) {
  sample_tbl <- unique(model_meta[, c("sample_id", "initial_ploidy", "treatment"), drop = FALSE])
  sample_tbl <- sample_tbl[order(sample_tbl$initial_ploidy, sample_tbl$sample_id), , drop = FALSE]
  ploidies <- sort(unique(sample_tbl$initial_ploidy))
  if (length(ploidies) != 2L) stop("Exact permutation expects two initial-ploidy strata.", call. = FALSE)
  combs <- lapply(ploidies, function(pl) {
    ids <- sort(sample_tbl$sample_id[sample_tbl$initial_ploidy == pl])
    if (length(ids) != 8L) stop("Exact permutation expects 8 mice in ploidy stratum ", pl, call. = FALSE)
    utils::combn(ids, treated_per_ploidy, simplify = FALSE)
  })
  names(combs) <- ploidies
  rows <- list()
  k <- 0L
  for (a in combs[[1L]]) {
    for (b in combs[[2L]]) {
      k <- k + 1L
      treated <- sort(c(a, b))
      rows[[k]] <- data.frame(
        permutation_id = sprintf("perm_%04d", k),
        treated_samples = paste(treated, collapse = ";"),
        assignment_checksum = ptpv3_assignment_checksum(sample_tbl$sample_id, treated),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out
}

ptpv3_assignment_checksum <- function(sample_ids, treated_samples) {
  sample_ids <- sort(unique(as.character(sample_ids)))
  labels <- data.frame(sample_id = sample_ids,
                       permuted_treatment = ifelse(sample_ids %in% treated_samples, "treated", "control"),
                       stringsAsFactors = FALSE)
  row.names(labels) <- NULL
  digest::digest(labels, algo = "sha256")
}

ptpv3_observed_treated_samples <- function(model_meta) {
  if ("dose_mg" %in% names(model_meta)) {
    dose <- suppressWarnings(as.numeric(model_meta$dose_mg))
    treated <- model_meta$sample_id[is.finite(dose) & dose > 0]
    if (length(treated)) return(sort(unique(treated)))
  }
  treatment <- tolower(trimws(as.character(model_meta$treatment %||% "")))
  control_labels <- c("control", "vehicle", "untreated", "none", "0", "na", "")
  sort(unique(model_meta$sample_id[!(treatment %in% control_labels)]))
}

ptpv3_contrast_weight_vector <- function(model_meta, treated_samples, eligible_subbins) {
  meta <- model_meta
  meta$perm_treatment <- ifelse(meta$sample_id %in% treated_samples, "treated", "control")
  w <- rep(0, nrow(meta))
  n_subbins <- length(eligible_subbins)
  for (sb in eligible_subbins) {
    for (pl in unique(meta$initial_ploidy)) {
      is4 <- grepl("4", pl)
      for (tr in c("treated", "control")) {
        idx <- which(meta$subbin_id == sb & meta$initial_ploidy == pl & meta$perm_treatment == tr)
        if (!length(idx)) next
        sign <- if (is4 && tr == "treated") 1 else if (is4 && tr == "control") -1 else if (!is4 && tr == "treated") -1 else 1
        w[idx] <- sign / (length(idx) * n_subbins)
      }
    }
  }
  w
}

ptpv3_program_indices <- function(model, programs) {
  row_symbols <- ptp_clean_gene_symbols(rownames(model$voom$E))
  lapply(programs, function(p) {
    idx <- which(row_symbols %in% p$observed_genes)
    idx[!duplicated(row_symbols[idx])]
  })
}

ptpv3_exact_permutation <- function(primary_fit, programs, role_table, cfg, args, dirs, eligible_subbins) {
  meta <- primary_fit$metadata
  perms <- ptpv3_assignment_table(meta, treated_per_ploidy = as.integer(cfg$v3$exact_permutation$treatment_labels_per_ploidy %||% 4L))
  observed_treated <- ptpv3_observed_treated_samples(meta)
  obs_checksum <- ptpv3_assignment_checksum(meta$sample_id, observed_treated)
  perms$is_observed_assignment <- perms$assignment_checksum == obs_checksum
  weights_obs <- ptpv3_contrast_weight_vector(meta, observed_treated, eligible_subbins)
  E <- ptp_logcpm_matrix(primary_fit)
  Ez <- t(scale(t(E)))
  Ez[!is.finite(Ez)] <- 0
  indices <- ptpv3_program_indices(primary_fit, programs)
  program_ids <- names(indices)[vapply(indices, length, integer(1L)) >= 2L]
  comp <- c("directional_signed_mean", "maxmean", "msq", "maxabs")
  B <- nrow(perms)
  P <- length(program_ids)
  C <- length(comp)
  obs <- matrix(NA_real_, nrow = P, ncol = C, dimnames = list(program_ids, comp))
  perm_stats <- array(NA_real_, dim = c(B, P, C), dimnames = list(perms$permutation_id, program_ids, comp))
  stat_components <- function(gstat, p, directional) {
    if (length(gstat) < 2L) return(rep(NA_real_, C))
    w <- p$observed_gene_weights$signed_weight[match(ptp_clean_gene_symbols(rownames(Ez))[indices[[p$id]]], p$observed_gene_weights$gene)]
    w[!is.finite(w) | w == 0] <- 1
    signed <- if (directional) sum(w * gstat, na.rm = TRUE) / sum(abs(w), na.rm = TRUE) else NA_real_
    c(directional_signed_mean = signed,
      maxmean = abs(mean(gstat, na.rm = TRUE)),
      msq = mean(gstat^2, na.rm = TRUE),
      maxabs = max(abs(gstat), na.rm = TRUE))
  }
  gobs <- as.numeric(Ez %*% weights_obs / sqrt(sum(weights_obs^2)))
  names(gobs) <- rownames(Ez)
  for (j in seq_along(program_ids)) {
    id <- program_ids[[j]]
    obs[j, ] <- stat_components(gobs[indices[[id]]], programs[[id]], role_table$directional_eligibility[match(id, role_table$program_id)] %in% TRUE)
  }
  chunk_size <- 200L
  chunks <- split(seq_len(B), ceiling(seq_len(B) / chunk_size))
  for (chunk in chunks) {
    W <- vapply(chunk, function(i) {
      treated <- strsplit(perms$treated_samples[[i]], ";", fixed = TRUE)[[1L]]
      ptpv3_contrast_weight_vector(meta, treated, eligible_subbins)
    }, numeric(nrow(meta)))
    W <- sweep(W, 2L, sqrt(colSums(W^2)), "/")
    G <- Ez %*% W
    for (j in seq_along(program_ids)) {
      id <- program_ids[[j]]
      idx <- indices[[id]]
      directional <- role_table$directional_eligibility[match(id, role_table$program_id)] %in% TRUE
      w <- programs[[id]]$observed_gene_weights$signed_weight[match(ptp_clean_gene_symbols(rownames(Ez))[idx], programs[[id]]$observed_gene_weights$gene)]
      w[!is.finite(w) | w == 0] <- 1
      vals <- G[idx, , drop = FALSE]
      if (directional) perm_stats[chunk, j, "directional_signed_mean"] <- as.numeric(crossprod(w / sum(abs(w)), vals))
      perm_stats[chunk, j, "maxmean"] <- abs(colMeans(vals, na.rm = TRUE))
      perm_stats[chunk, j, "msq"] <- colMeans(vals^2, na.rm = TRUE)
      perm_stats[chunk, j, "maxabs"] <- apply(abs(vals), 2L, max, na.rm = TRUE)
    }
  }
  component_rows <- list()
  adaptive_rows <- list()
  minp_perm <- matrix(NA_real_, nrow = B, ncol = P, dimnames = list(perms$permutation_id, program_ids))
  for (j in seq_along(program_ids)) {
    id <- program_ids[[j]]
    comp_p <- rep(NA_real_, C)
    for (k in seq_along(comp)) {
      vals <- perm_stats[, j, k]
      o <- obs[j, k]
      if (!is.finite(o) || all(!is.finite(vals))) next
      extreme_obs <- if (comp[[k]] == "directional_signed_mean") abs(vals) >= abs(o) else vals >= o
      comp_p[[k]] <- (1 + sum(extreme_obs, na.rm = TRUE)) / (1 + sum(is.finite(vals)))
      ranks <- vapply(vals, function(v) {
        if (!is.finite(v)) return(NA_real_)
        ext <- if (comp[[k]] == "directional_signed_mean") abs(vals) >= abs(v) else vals >= v
        (1 + sum(ext, na.rm = TRUE)) / (1 + sum(is.finite(vals)))
      }, numeric(1L))
      if (all(is.na(minp_perm[, j]))) minp_perm[, j] <- ranks else minp_perm[, j] <- pmin(minp_perm[, j], ranks, na.rm = TRUE)
      component_rows[[length(component_rows) + 1L]] <- data.frame(
        program_id = id,
        statistic_component = comp[[k]],
        observed_statistic = o,
        empirical_p = comp_p[[k]],
        observed_rank = 1 + sum(if (comp[[k]] == "directional_signed_mean") abs(vals) > abs(o) else vals > o, na.rm = TRUE),
        finite_sample_formula = cfg$v3$reproducibility$finite_sample_p_formula %||% "(1 + extreme) / (1 + permutations)",
        status = if (comp[[k]] == "directional_signed_mean" && !(role_table$directional_eligibility[match(id, role_table$program_id)] %in% TRUE)) "not_estimable_no_legal_direction" else "ok",
        stringsAsFactors = FALSE
      )
    }
    obs_minp <- min(comp_p, na.rm = TRUE)
    adaptive_p <- (1 + sum(minp_perm[, j] <= obs_minp, na.rm = TRUE)) / (1 + sum(is.finite(minp_perm[, j])))
    adaptive_rows[[length(adaptive_rows) + 1L]] <- data.frame(program_id = id, method_id = "exact_permutation_adaptive_omnibus",
                                                              observed_min_component_p = obs_minp, p_value = adaptive_p,
                                                              n_permutations = B, stringsAsFactors = FALSE)
  }
  maxT <- apply(-log10(pmax(minp_perm, .Machine$double.xmin)), 1L, max, na.rm = TRUE)
  adaptive <- ptpv3_bind_rows(adaptive_rows)
  adaptive$westfall_young_maxT_p <- vapply(seq_len(nrow(adaptive)), function(i) {
    obs_score <- -log10(pmax(adaptive$observed_min_component_p[[i]], .Machine$double.xmin))
    (1 + sum(maxT >= obs_score, na.rm = TRUE)) / (1 + sum(is.finite(maxT)))
  }, numeric(1L))
  adaptive <- ptpv3_apply_fdr_family(adaptive, "p_value", "method_id", role_table)
  components <- ptpv3_apply_fdr_family(ptpv3_bind_rows(component_rows), "empirical_p", "statistic_component", role_table)
  validation <- data.frame(
    check = c("unique_assignments", "expected_assignments", "observed_assignment_present", "whole_mouse_unit", "ploidy_strata_preserved"),
    status = c(if (length(unique(perms$assignment_checksum)) == B) "ok" else "failed",
               if (B == as.integer(cfg$v3$exact_permutation$expected_unique_assignments %||% 4900L)) "ok" else "failed",
               if (any(perms$is_observed_assignment)) "ok" else "failed",
               "ok", "ok"),
    value = c(length(unique(perms$assignment_checksum)), B, sum(perms$is_observed_assignment), "sample_id", paste(sort(unique(meta$initial_ploidy)), collapse = ";")),
    stringsAsFactors = FALSE
  )
  ptp_write_csv(perms, file.path(dirs$primary$exact, "permutation_assignments.csv"))
  ptp_write_csv(validation, file.path(dirs$primary$exact, "permutation_validation.csv"))
  ptp_write_csv(components, file.path(dirs$primary$exact, "exact_permutation_component_pvalues.csv"))
  ptp_write_csv(adaptive, file.path(dirs$primary$exact, "exact_permutation_adaptive_omnibus.csv"))
  saveRDS(list(permutations = perms, component_stats = perm_stats, minp_perm = minp_perm),
          file.path(dirs$primary$exact, "exact_permutation_null_distributions.rds"))
  list(permutations = perms, validation = validation, components = components, adaptive = adaptive, minp_perm = minp_perm)
}

ptpv3_external_signature_status <- function(cfg, repo_root, dirs) {
  roots <- c("Data", "Code", "Figs", "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels")
  roots <- vapply(roots, function(x) if (grepl("^/", x)) x else file.path(repo_root, x), character(1L))
  files <- unlist(lapply(roots[file.exists(roots)], function(root) {
    list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  }), use.names = FALSE)
  hit <- files[grepl("in[-_ ]?vitro|gemcitabine.*response|response.*gemcitabine", basename(files), ignore.case = TRUE)]
  hit <- hit[!grepl("04j_pseudotime_treatment_ploidy_programs_v[23]/05_signature_projection", hit, fixed = FALSE)]
  candidates <- data.frame(path = normalizePath(hit, winslash = "/", mustWork = FALSE),
                           sha256 = vapply(hit, ptp_file_checksum, character(1L)),
                           stringsAsFactors = FALSE)
  status <- data.frame(method_id = "independent_compact_directional_score",
                       status = if (nrow(candidates) == 1L) "candidate_requires_manual_schema_review" else "not_estimable",
                       reason = if (nrow(candidates) == 1L) "one filename candidate found but no frozen parser/schema is defined" else "no unique auditable external in-vitro gemcitabine response source found",
                       n_candidates = nrow(candidates),
                       stringsAsFactors = FALSE)
  ptp_write_csv(candidates, file.path(dirs$directional, "external_invitro_candidate_search.csv"))
  ptp_write_csv(status, file.path(dirs$directional, "independent_compact_directional_score_status.csv"))
  status
}

ptpv3_control_reference_residual <- function(primary_fit, programs, role_table, cfg, dirs) {
  meta <- primary_fit$metadata
  E <- ptp_logcpm_matrix(primary_fit)
  groups <- paste(meta$subbin_id, meta$initial_ploidy, sep = "__")
  is_control <- meta$treatment == "control"
  R <- E * NA_real_
  audit_rows <- list()
  prior <- as.numeric(cfg$v3$control_reference_residual$sd_shrinkage_prior_weight %||% 4)
  global_sd <- apply(E, 1L, stats::sd, na.rm = TRUE)
  for (g in unique(groups)) {
    train <- which(groups == g & is_control)
    test <- which(groups == g)
    mu <- rowMeans(E[, train, drop = FALSE], na.rm = TRUE)
    sd0 <- apply(E[, train, drop = FALSE], 1L, stats::sd, na.rm = TRUE)
    sd_shrunk <- sqrt(((pmax(length(train) - 1L, 1L) * sd0^2) + prior * global_sd^2) / (pmax(length(train) - 1L, 1L) + prior))
    sd_shrunk[!is.finite(sd_shrunk) | sd_shrunk <= 0] <- global_sd[!is.finite(sd_shrunk) | sd_shrunk <= 0]
    R[, test] <- (E[, test, drop = FALSE] - mu) / pmax(sd_shrunk, 1e-6)
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(reference_group = g, n_control_train_observations = length(train),
                                                        n_projected_observations = length(test),
                                                        median_sd_shrinkage = stats::median(sd_shrunk / pmax(sd0, 1e-6), na.rm = TRUE),
                                                        stringsAsFactors = FALSE)
  }
  scores <- ptp_program_scores_from_logcpm(R, programs)
  pfit <- ptp_fit_program_scores(scores, primary_fit)
  contrasts <- ptp_primary_contrasts(unique(meta$subbin_id[meta$retained_for_model]), colnames(primary_fit$design))
  res <- ptp_apply_program_contrasts(pfit, contrasts, "control_reference_residual_score", programs)
  res$score_method_id <- "control_reference_residual_score"
  res$score_method_label <- "Control-reference residual score"
  res <- ptpv3_apply_fdr_family(res, "p_value", "score_method_id", role_table)
  ptp_write_csv(ptpv3_bind_rows(audit_rows), file.path(dirs$controls, "control_reference_training_audit.csv"))
  ptp_write_csv(res, file.path(dirs$primary$score_models, "control_reference_residual_score_interactions.csv"))
  res
}

ptpv3_covariance_whitened <- function(primary_fit, programs, role_table, cfg, dirs) {
  meta <- primary_fit$metadata
  E <- t(scale(t(ptp_logcpm_matrix(primary_fit))))
  E[!is.finite(E)] <- 0
  row_symbols <- ptp_clean_gene_symbols(rownames(E))
  lambda <- as.numeric(cfg$v3$covariance_whitened_score$ridge_lambda %||% 0.1)
  max_genes <- as.integer(cfg$v3$covariance_whitened_score$max_genes %||% 80L)
  scores <- matrix(NA_real_, nrow = length(programs), ncol = ncol(E), dimnames = list(names(programs), colnames(E)))
  audit <- list()
  for (id in names(programs)) {
    p <- programs[[id]]
    directional <- role_table$directional_eligibility[match(id, role_table$program_id)] %in% TRUE
    if (!directional) {
      audit[[length(audit) + 1L]] <- data.frame(program_id = id, status = "not_estimable", reason = "no legal frozen directional weight vector", stringsAsFactors = FALSE)
      next
    }
    ww <- ptpv2_program_weight_vector(p, row_symbols)
    if (length(ww$idx) < 2L) {
      audit[[length(audit) + 1L]] <- data.frame(program_id = id, status = "not_estimable", reason = "fewer than two observed weighted genes", stringsAsFactors = FALSE)
      next
    }
    if (length(ww$idx) > max_genes) {
      ord <- order(abs(ww$weights), decreasing = TRUE)[seq_len(max_genes)]
      ww$idx <- ww$idx[ord]
      ww$weights <- ww$weights[ord]
    }
    train <- which(meta$treatment == "control")
    X <- t(E[ww$idx, train, drop = FALSE])
    S <- stats::cov(X, use = "pairwise.complete.obs")
    S[!is.finite(S)] <- 0
    diagS <- diag(S)
    ridge <- lambda * stats::median(diagS[is.finite(diagS) & diagS > 0], na.rm = TRUE)
    if (!is.finite(ridge) || ridge <= 0) ridge <- lambda
    Sr <- S + diag(ridge, nrow(S))
    ev <- tryCatch(eigen(Sr, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
    condition <- if (all(is.finite(ev)) && min(ev) > 0) max(ev) / min(ev) else Inf
    inv <- tryCatch(solve(Sr), error = function(e) NULL)
    if (is.null(inv) || !is.finite(condition) || condition > as.numeric(cfg$v3$covariance_whitened_score$condition_number_fallback %||% 1e6)) {
      audit[[length(audit) + 1L]] <- data.frame(program_id = id, status = "fallback_visible", reason = "covariance inversion unstable", ridge_lambda = lambda,
                                                condition_number = condition, n_genes = length(ww$idx), stringsAsFactors = FALSE)
      next
    }
    a <- matrix(ww$weights, ncol = 1L)
    denom <- sqrt(as.numeric(t(a) %*% inv %*% a))
    scores[id, ] <- as.numeric(t(a) %*% inv %*% E[ww$idx, , drop = FALSE] / denom)
    audit[[length(audit) + 1L]] <- data.frame(program_id = id, status = "ok", reason = "", ridge_lambda = lambda,
                                              condition_number = condition, effective_rank = (sum(ev)^2) / sum(ev^2),
                                              n_genes = length(ww$idx), stringsAsFactors = FALSE)
  }
  keep <- rowSums(is.finite(scores)) == ncol(scores)
  res <- data.frame()
  if (any(keep)) {
    pfit <- ptp_fit_program_scores(scores[keep, , drop = FALSE], primary_fit)
    contrasts <- ptp_primary_contrasts(unique(meta$subbin_id[meta$retained_for_model]), colnames(primary_fit$design))
    res <- ptp_apply_program_contrasts(pfit, contrasts, "covariance_whitened_score", programs[rownames(scores)[keep]])
    res$score_method_id <- "covariance_whitened_score"
    res <- ptpv3_apply_fdr_family(res, "p_value", "score_method_id", role_table)
  }
  ptp_write_csv(ptpv3_bind_rows(audit), file.path(dirs$directional, "covariance_whitened_score_audit.csv"))
  ptp_write_csv(res, file.path(dirs$primary$score_models, "covariance_whitened_score_interactions.csv"))
  res
}

ptpv3_crossfit_04i <- function(primary_fit, programs, role_table, cfg, args, repo_root, dirs) {
  if (!isTRUE(args$run_crossfit_04i)) return(list(status = data.frame(status = "skipped"), results = data.frame()))
  meta <- pst_read_cell_metadata(args$cell_metadata)
  counts <- pst_load_counts(ptp_resolve_seurat_rds(args$seurat_rds, cfg, args$results_root), args$assay, args$counts_layer)
  matched <- pst_match_cells(meta, counts, 0.99, dirs$crossfit)
  meta <- matched$meta
  counts <- matched$counts
  program_gene_universe <- unique(unlist(lapply(programs, `[[`, "observed_genes"), use.names = FALSE))
  keep_gene <- ptp_clean_gene_symbols(rownames(counts)) %in% program_gene_universe
  counts <- counts[keep_gene, , drop = FALSE]
  gene_universe_checksum <- digest::digest(sort(ptp_clean_gene_symbols(rownames(counts))), algo = "sha256")
  pst_cfg <- pst_read_config(file.path(repo_root, "Code/in-vivo/04i_pseudotime_state_pathways_config.yaml"))
  model_spec <- list(model_id = "primary_initial_ploidy", covariate_terms = c("initial_ploidy_factor"),
                     covariate_mode = "initial_ploidy", etp_method = NA_character_, etp_threshold = NA_real_)
  sample_ids <- sort(unique(meta$sample_id))
  primary_logcpm <- ptp_logcpm_matrix(primary_fit)
  gene_z <- t(scale(t(primary_logcpm)))
  gene_z[!is.finite(gene_z)] <- 0
  row_symbols <- ptp_clean_gene_symbols(rownames(gene_z))
  families <- sort(unique(role_table$response_family[role_table$tier %in% c("tier1_core", "tier2_upgraded_focused", "tier3_exploratory_only")]))
  score_mat <- matrix(NA_real_, nrow = length(families), ncol = ncol(gene_z),
                      dimnames = list(paste0("crossfit_04i_", ptp_safe_name(families)), colnames(gene_z)))
  fold_rows <- list()
  weight_rows <- list()
  worker_fun <- function(held) {
    train_meta <- meta[meta$sample_id != held, , drop = FALSE]
    train_counts <- counts[, train_meta$cell_id, drop = FALSE]
    held_cols <- which(primary_fit$metadata$sample_id == held)
    out <- tryCatch({
      pb <- pst_construct_pseudobulk(train_meta, train_counts, args$n_pseudotime_bins %||% 20L, args$min_cells_per_sample_bin %||% 5L)
      coverage <- pst_region_coverage(pb$cell_metadata, pb$metadata, pst_cfg)
      coverage_check <- pst_check_primary_coverage(coverage)
      mf <- pst_fit_pseudotime_model(pb$counts, pb$metadata, args$spline_df %||% 5L,
                                     coverage_check$n_contributing_mice[[1L]],
                                     include_dose = TRUE, model_spec = model_spec)
      contrast <- pst_contrast_vector(mf, pst_cfg, "primary_adjacent_state", args$grid_size %||% 501L)
      genes <- pst_contrast_table(mf, contrast, "primary_adjacent_state")
      genes$gene_symbol <- ptp_clean_gene_symbols(genes$gene_symbol)
      fold_weights <- list()
      fold_scores <- matrix(NA_real_, nrow = length(families), ncol = length(held_cols),
                            dimnames = list(paste0("crossfit_04i_", ptp_safe_name(families)), colnames(gene_z)[held_cols]))
      for (family in families) {
        fam_programs <- role_table$program_id[role_table$response_family == family]
        fam_genes <- unique(unlist(lapply(programs[fam_programs], `[[`, "observed_genes"), use.names = FALSE))
        g <- genes[genes$gene_symbol %in% fam_genes & is.finite(genes$t_statistic), , drop = FALSE]
        if (nrow(g) < 2L) next
        w <- g$t_statistic
        q <- stats::quantile(w, probs = c(0.01, 0.99), na.rm = TRUE)
        w <- pmin(pmax(w, q[[1L]]), q[[2L]])
        w <- w / sum(abs(w), na.rm = TRUE)
        idx <- match(g$gene_symbol, row_symbols)
        keep <- is.finite(idx)
        idx <- idx[keep]
        w <- w[keep]
        if (length(idx) < 2L) next
        fold_scores[paste0("crossfit_04i_", ptp_safe_name(family)), ] <- as.numeric(crossprod(w, gene_z[idx, held_cols, drop = FALSE]))
        fold_weights[[length(fold_weights) + 1L]] <- data.frame(held_out_mouse_id = held, response_family = family,
                                                                 gene_symbol = row_symbols[idx], weight = w,
                                                                 stringsAsFactors = FALSE)
      }
      list(status = data.frame(held_out_mouse_id = held,
                               training_mouse_ids = paste(sort(unique(train_meta$sample_id)), collapse = ";"),
                               training_contains_heldout = held %in% unique(train_meta$sample_id),
                               training_input_checksum = digest::digest(sort(train_meta$cell_id), algo = "sha256"),
                               training_gene_universe = "frozen_04j_program_genes_only",
                               training_gene_universe_checksum = gene_universe_checksum,
                               fitted_weight_checksum = digest::digest(ptpv3_bind_rows(fold_weights), algo = "sha256"),
                               held_out_projection_checksum = digest::digest(fold_scores, algo = "sha256"),
                               nonzero_genes = sum(abs(ptpv3_bind_rows(fold_weights)$weight) > 0, na.rm = TRUE),
                               status = "ok", error = "", stringsAsFactors = FALSE),
           weights = ptpv3_bind_rows(fold_weights), scores = fold_scores)
    }, error = function(e) {
      list(status = data.frame(held_out_mouse_id = held, training_mouse_ids = "",
                               training_contains_heldout = NA, training_input_checksum = "",
                               training_gene_universe = "frozen_04j_program_genes_only",
                               training_gene_universe_checksum = gene_universe_checksum,
                               fitted_weight_checksum = "", held_out_projection_checksum = "",
                               nonzero_genes = 0L, status = "failed", error = conditionMessage(e),
                               stringsAsFactors = FALSE),
           weights = data.frame(), scores = matrix(NA_real_, nrow = length(families), ncol = length(held_cols),
                                                   dimnames = list(paste0("crossfit_04i_", ptp_safe_name(families)), colnames(gene_z)[held_cols])))
    })
    out
  }
  workers <- ptpv2_workers(args, length(sample_ids), "crossfit_workers")
  res <- ptp_parallel_lapply(sample_ids, worker_fun, workers = workers, task_label = "v3_crossfit_04i")$results
  for (r in res) {
    fold_rows[[length(fold_rows) + 1L]] <- r$status
    if (nrow(r$weights) > 0L) weight_rows[[length(weight_rows) + 1L]] <- r$weights
    score_mat[rownames(r$scores), colnames(r$scores)] <- r$scores
  }
  fold_audit <- ptpv3_bind_rows(fold_rows)
  weights <- ptpv3_bind_rows(weight_rows)
  score_long <- data.frame(sample_subbin_id = colnames(score_mat), t(score_mat), check.names = FALSE)
  res_model <- data.frame()
  if (all(rowSums(is.finite(score_mat)) == ncol(score_mat))) {
    projection_programs <- lapply(rownames(score_mat), function(id) {
      list(id = id, label = id, response_family = id, response_family_label = id,
           tier = "projection", tier_label = "Projection", family = id, upgrade_status = "exploratory_bridge",
           expected_direction = "external", score_mode = "projection", score_direction_label = "cross-fitted 04i state",
           score_direction_source = "04i refit without held-out mouse", score_component_policy = "winsorized_l1",
           source = "04i_crossfit", set_names = character(0L), missing_set_names = character(0L), genes = character(0L),
           estimable = TRUE, component_gene_weights = data.frame(), observed_genes = character(0L),
           observed_gene_weights = data.frame(gene = character(0L), signed_weight = numeric(0L), stringsAsFactors = FALSE))
    })
    names(projection_programs) <- rownames(score_mat)
    pfit <- ptp_fit_program_scores(score_mat, primary_fit)
    contrasts <- ptp_primary_contrasts(unique(primary_fit$metadata$subbin_id[primary_fit$metadata$retained_for_model]), colnames(primary_fit$design))
    res_model <- ptp_apply_program_contrasts(pfit, contrasts, "cross_fitted_04i_state_projection", projection_programs)
  }
  status <- data.frame(method_id = "cross_fitted_04i_state_projection",
                       status = if (all(fold_audit$status == "ok") && all(!fold_audit$training_contains_heldout)) "ok" else "failed",
                       n_folds = nrow(fold_audit),
                       n_failed = sum(fold_audit$status != "ok"),
                       stringsAsFactors = FALSE)
  ptp_write_csv(status, file.path(dirs$crossfit, "cross_fitted_04i_projection_status.csv"))
  ptp_write_csv(fold_audit, file.path(dirs$crossfit, "cross_fitted_04i_fold_audit.csv"))
  ptp_write_csv(weights, file.path(dirs$crossfit, "cross_fitted_04i_fold_weights.csv"))
  ptp_write_csv(score_long, file.path(dirs$crossfit, "cross_fitted_04i_out_of_fold_scores.csv"))
  ptp_write_csv(res_model, file.path(dirs$crossfit, "cross_fitted_04i_projection_model_results.csv"))
  list(status = status, fold_audit = fold_audit, weights = weights, scores = score_long, results = res_model)
}

ptpv3_multidimensional <- function(primary_fit, programs, role_table, exact, cfg, dirs, eligible_subbins) {
  meta <- primary_fit$metadata
  E <- t(scale(t(ptp_logcpm_matrix(primary_fit))))
  E[!is.finite(E)] <- 0
  row_symbols <- ptp_clean_gene_symbols(rownames(E))
  weights_obs <- ptpv3_contrast_weight_vector(meta, sort(unique(meta$sample_id[meta$treatment == "treated"])), eligible_subbins)
  perms <- exact$permutations
  rows <- list()
  loading_rows <- list()
  for (id in names(programs)) {
    idx <- which(row_symbols %in% programs[[id]]$observed_genes)
    idx <- idx[!duplicated(row_symbols[idx])]
    if (length(idx) < 3L) {
      rows[[length(rows) + 1L]] <- data.frame(program_id = id, method_id = "multidimensional_pathway_interaction",
                                              status = "not_estimable", reason = "fewer than three observed genes", p_value = NA_real_,
                                              stringsAsFactors = FALSE)
      next
    }
    train <- which(meta$treatment == "control")
    pc <- tryCatch(stats::prcomp(t(E[idx, train, drop = FALSE]), center = TRUE, scale. = FALSE), error = function(e) NULL)
    if (is.null(pc) || ncol(pc$rotation) == 0L) {
      rows[[length(rows) + 1L]] <- data.frame(program_id = id, method_id = "multidimensional_pathway_interaction",
                                              status = "failed", reason = "PCA failed", p_value = NA_real_, stringsAsFactors = FALSE)
      next
    }
    varex <- pc$sdev^2 / sum(pc$sdev^2)
    k <- min(which(cumsum(varex) >= as.numeric(cfg$v3$multidimensional_pathway$variance_explained_target %||% 0.8))[1L],
             as.integer(cfg$v3$multidimensional_pathway$max_pcs %||% 3L), ncol(pc$rotation))
    basis <- pc$rotation[, seq_len(k), drop = FALSE]
    all_scores <- t(E[idx, , drop = FALSE]) %*% basis
    obs_stat <- sum(as.numeric(crossprod(weights_obs, all_scores))^2)
    perm_stat <- vapply(seq_len(nrow(perms)), function(i) {
      treated <- strsplit(perms$treated_samples[[i]], ";", fixed = TRUE)[[1L]]
      w <- ptpv3_contrast_weight_vector(meta, treated, eligible_subbins)
      sum(as.numeric(crossprod(w, all_scores))^2)
    }, numeric(1L))
    pval <- (1 + sum(perm_stat >= obs_stat, na.rm = TRUE)) / (1 + sum(is.finite(perm_stat)))
    rows[[length(rows) + 1L]] <- data.frame(program_id = id, method_id = "multidimensional_pathway_interaction",
                                            status = "ok", reason = "", n_pcs = k,
                                            variance_explained = sum(varex[seq_len(k)]),
                                            joint_statistic = obs_stat, p_value = pval,
                                            n_permutations = nrow(perms), stringsAsFactors = FALSE)
    for (pcid in seq_len(k)) {
      ord <- order(abs(basis[, pcid]), decreasing = TRUE)[seq_len(min(20L, nrow(basis)))]
      loading_rows[[length(loading_rows) + 1L]] <- data.frame(program_id = id, pc = pcid,
                                                              gene_symbol = row_symbols[idx][ord],
                                                              loading = basis[ord, pcid],
                                                              stringsAsFactors = FALSE)
    }
  }
  out <- ptpv3_apply_fdr_family(ptpv3_bind_rows(rows), "p_value", "method_id", role_table)
  ptp_write_csv(out, file.path(dirs$primary$multidim, "multidimensional_pathway_interactions.csv"))
  ptp_write_csv(ptpv3_bind_rows(loading_rows), file.path(dirs$primary$multidim, "multidimensional_pc_loadings_top_genes.csv"))
  out
}

ptpv3_trajectory_global <- function(meta, counts, programs, role_table, cfg, args, dirs) {
  bins <- lapply(seq_len(as.integer(cfg$v3$trajectory_global$n_bins %||% 20L)), function(i) {
    start <- (i - 1L) / as.integer(cfg$v3$trajectory_global$n_bins %||% 20L)
    end <- i / as.integer(cfg$v3$trajectory_global$n_bins %||% 20L)
    list(id = sprintf("bin%02d", i), label = sprintf("[%.2f,%.2f%s", start, end, if (i == as.integer(cfg$v3$trajectory_global$n_bins %||% 20L)) "]" else ")"),
         start = start, end = end, include_start = TRUE, include_end = i == as.integer(cfg$v3$trajectory_global$n_bins %||% 20L))
  })
  pb <- ptp_construct_subbin_pseudobulk(meta, counts, bins, cfg$qc, ptp_endpoint_ploidy_specs(cfg))
  keep_obs <- pb$metadata$retained_for_model & pb$metadata$library_size > 0
  md <- pb$metadata[keep_obs, , drop = FALSE]
  cts <- pb$counts[, md$sample_subbin_id, drop = FALSE]
  y <- edgeR::DGEList(counts = cts)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 1)
  scores <- ptp_program_scores_from_logcpm(logcpm, programs)
  bin_table <- data.frame(subbin_id = vapply(bins, `[[`, character(1L), "id"),
                          bin_midpoint = vapply(bins, function(x) (as.numeric(x$start) + as.numeric(x$end)) / 2, numeric(1L)),
                          stringsAsFactors = FALSE)
  md$bin_midpoint <- bin_table$bin_midpoint[match(md$subbin_id, bin_table$subbin_id)]
  keep_design <- stats::complete.cases(md[, c("sample_subbin_id", "bin_midpoint", "initial_ploidy", "treatment", "sample_id"), drop = FALSE])
  md <- md[keep_design, , drop = FALSE]
  logcpm <- logcpm[, md$sample_subbin_id, drop = FALSE]
  scores <- ptp_program_scores_from_logcpm(logcpm, programs)
  md$initial_ploidy <- factor(md$initial_ploidy)
  md$treatment <- ifelse(suppressWarnings(as.numeric(md$dose_mg)) == 0, "control", "treated")
  md$treatment <- factor(md$treatment, levels = c("control", "treated"))
  design <- stats::model.matrix(~ splines::ns(bin_midpoint, df = as.integer(cfg$v3$trajectory_global$spline_df %||% 3L)) * initial_ploidy * treatment, data = md)
  if (nrow(design) != nrow(md)) {
    kept <- suppressWarnings(as.integer(rownames(design)))
    kept <- kept[is.finite(kept) & kept >= 1L & kept <= nrow(md)]
    md <- md[kept, , drop = FALSE]
    logcpm <- logcpm[, md$sample_subbin_id, drop = FALSE]
    scores <- ptp_program_scores_from_logcpm(logcpm, programs)
  }
  scores <- as.matrix(scores)
  colnames(scores) <- colnames(logcpm)
  rownames(design) <- md$sample_subbin_id
  scores <- scores[, md$sample_subbin_id, drop = FALSE]
  keep_score <- rowSums(is.finite(scores)) == ncol(scores)
  scores <- scores[keep_score, , drop = FALSE]
  if (nrow(scores) == 0L || ncol(scores) != nrow(design)) {
    status <- data.frame(method_id = "trajectory_wide_global_interaction", status = "failed",
                         reason = paste0("score/design dimension mismatch: score_cols=", ncol(scores), ", design_rows=", nrow(design)),
                         stringsAsFactors = FALSE)
    ptp_write_csv(status, file.path(dirs$primary$trajectory, "trajectory_wide_global_interactions.csv"))
    return(status)
  }
  dup <- limma::duplicateCorrelation(scores, design, block = md$sample_id)
  fit <- limma::lmFit(scores, design, block = md$sample_id, correlation = dup$consensus.correlation)
  fit <- limma::eBayes(fit, robust = TRUE)
  coef_idx <- grep("bin_midpoint.*initial_ploidy.*treatment|bin_midpoint.*treatment.*initial_ploidy", colnames(design))
  if (!length(coef_idx)) coef_idx <- grep(":", colnames(design))
  tt <- limma::topTable(fit, coef = coef_idx, number = Inf, sort.by = "none")
  tt$program_id <- rownames(tt)
  names(tt)[names(tt) == "P.Value"] <- "p_value"
  out <- data.frame(program_id = tt$program_id, method_id = "trajectory_wide_global_interaction",
                    global_F = tt$F, p_value = tt$p_value, spline_df = as.integer(cfg$v3$trajectory_global$spline_df %||% 3L),
                    duplicate_correlation = dup$consensus.correlation, stringsAsFactors = FALSE)
  out <- ptpv3_apply_fdr_family(out, "p_value", "method_id", role_table)
  ptp_write_csv(md, file.path(dirs$primary$trajectory, "trajectory_20bin_pseudobulk_metadata.csv"))
  ptp_write_csv(out, file.path(dirs$primary$trajectory, "trajectory_wide_global_interactions.csv"))
  out
}

ptpv3_select_sim_programs <- function(programs, role_table, n = 9L) {
  x <- role_table[role_table$eligibility %in% TRUE, , drop = FALSE]
  x$size <- cut(x$n_genes_observed, breaks = c(-Inf, 20, 100, Inf), labels = c("small", "medium", "large"))
  ids <- unlist(lapply(levels(x$size), function(sz) head(x$program_id[x$size == sz][order(x$n_genes_observed[x$size == sz])], ceiling(n / 3))), use.names = FALSE)
  ids <- unique(ids)
  programs[intersect(ids, names(programs))]
}

ptpv3_run_simulation <- function(primary_fit, programs, role_table, cfg, args, dirs) {
  if (!isTRUE(args$run_simulation)) return(data.frame())
  sim_programs <- ptpv3_select_sim_programs(programs, role_table, as.integer(cfg$simulation$representative_program_count %||% 9L))
  scenarios <- ptpv3_bind_rows(lapply(cfg$simulation$scenarios, as.data.frame, stringsAsFactors = FALSE))
  ptp_write_csv(scenarios, file.path(dirs$simulation, "frozen_simulation_scenarios.csv"))
  E <- ptp_logcpm_matrix(primary_fit)
  fit_mean <- primary_fit$fit$coefficients %*% t(primary_fit$design)
  resid <- E - fit_mean
  meta <- primary_fit$metadata
  contrast_w <- ptpv3_contrast_weight_vector(meta, sort(unique(meta$sample_id[meta$treatment == "treated"])), unique(meta$subbin_id[meta$retained_for_model]))
  score_methods <- intersect(c("signed_weighted_mean_zscore", "signed_pca_eigengene", "signed_rank_mean", "trimmed_signed_mean_zscore"),
                             vapply(ptp_score_method_specs(cfg), function(x) x$id, character(1L)))
  worker <- function(sid) {
    sc <- scenarios[sid, , drop = FALSE]
    set.seed(as.integer(cfg$v3$reproducibility$simulation_seed %||% args$seed) + sid)
    reps <- as.integer(sc$replicates[[1L]])
    rows <- list()
    for (r in seq_len(reps)) {
      boot_cols <- sample(seq_len(ncol(resid)), ncol(resid), replace = TRUE)
      simE <- fit_mean + resid[, boot_cols, drop = FALSE]
      beta <- as.numeric(sc$standardized_effect_size[[1L]])
      if (is.finite(beta) && beta != 0) {
        for (p in sim_programs) {
          idx <- which(ptp_clean_gene_symbols(rownames(simE)) %in% p$observed_genes)
          if (!length(idx)) next
          active_n <- max(1L, floor(length(idx) * as.numeric(sc$active_gene_fraction[[1L]])))
          active <- idx[seq_len(min(active_n, length(idx)))]
          direction <- rep(1, length(active))
          if (identical(sc$direction_pattern[[1L]], "mixed") && length(direction) >= 2L) direction[seq(2, length(direction), by = 2)] <- -1
          if (identical(sc$direction_pattern[[1L]], "sparse") && length(direction) >= 2L) direction[-seq_len(max(1L, floor(length(direction) / 3)))] <- 0
          simE[active, ] <- simE[active, , drop = FALSE] + beta * outer(direction, contrast_w / max(abs(contrast_w)))
        }
      }
      defs <- ptp_program_score_methods_from_logcpm(simE, sim_programs, cfg)
      for (mid in intersect(score_methods, names(defs))) {
        pfit <- tryCatch(ptp_fit_program_scores(defs[[mid]]$score_matrix, primary_fit), error = function(e) NULL)
        if (is.null(pfit)) next
        cons <- ptp_primary_contrasts(unique(meta$subbin_id[meta$retained_for_model]), colnames(primary_fit$design))
        res <- ptp_apply_program_contrasts(pfit, cons["treatment_by_initial_ploidy_interaction"], "simulation", sim_programs)
        if (!nrow(res)) next
        res$scenario_id <- sc$id[[1L]]
        res$replicate <- r
        res$score_method_id <- mid
        res$true_effect <- beta
        rows[[length(rows) + 1L]] <- res[, c("scenario_id", "replicate", "score_method_id", "program_id", "estimate", "se", "ci_low", "ci_high", "p_value", "true_effect"), drop = FALSE]
      }
    }
    out <- ptpv3_bind_rows(rows)
    ptp_write_csv(out, file.path(dirs$simulation, "checkpoints", paste0("simulation_", ptp_safe_name(sc$id[[1L]]), ".csv")))
    out
  }
  workers <- ptpv2_workers(args, nrow(scenarios), "simulation_workers")
  raw <- ptpv3_bind_rows(ptp_parallel_lapply(seq_len(nrow(scenarios)), worker, workers = workers, task_label = "v3_simulation")$results)
  raw <- ptpv3_join_program_meta(raw, role_table)
  raw$reject_0_05 <- raw$p_value < 0.05
  raw$covered <- raw$ci_low <= raw$true_effect & raw$ci_high >= raw$true_effect
  summary <- ptpv3_bind_rows(lapply(split(raw, interaction(raw$scenario_id, raw$score_method_id, drop = TRUE)), function(df) {
    data.frame(scenario_id = df$scenario_id[[1L]], score_method_id = df$score_method_id[[1L]],
               n_tests = nrow(df), true_effect = unique(df$true_effect)[1L],
               type_I_error = if (unique(df$true_effect)[1L] == 0) mean(df$reject_0_05, na.rm = TRUE) else NA_real_,
               power = if (unique(df$true_effect)[1L] != 0) mean(df$reject_0_05, na.rm = TRUE) else NA_real_,
               bias = mean(df$estimate - df$true_effect, na.rm = TRUE),
               interval_coverage = mean(df$covered, na.rm = TRUE),
               direction_accuracy = mean(sign(df$estimate) == sign(df$true_effect), na.rm = TRUE),
               failure_rate = mean(!is.finite(df$p_value)),
               mcse_reject = sqrt(mean(df$reject_0_05, na.rm = TRUE) * (1 - mean(df$reject_0_05, na.rm = TRUE)) / nrow(df)),
               stringsAsFactors = FALSE)
  }))
  vectorized_status <- data.frame(layer = "vectorized_expansion", status = "not_used",
                                  reason = "All frozen scenarios were run through expression-scale score recomputation; no vectorized expansion was needed.",
                                  stringsAsFactors = FALSE)
  ptp_write_csv(raw, file.path(dirs$simulation, "simulation_replicate_results.csv"))
  ptp_write_csv(summary, file.path(dirs$simulation, "simulation_method_calibration_summary.csv"))
  ptp_write_csv(vectorized_status, file.path(dirs$simulation, "vectorized_expansion_status.csv"))
  summary
}

ptpv3_plot_or_skip <- function(df, path, x, y, color = NULL, title = "") {
  if (is.null(df) || nrow(df) == 0L || !all(c(x, y) %in% names(df))) return(FALSE)
  aes <- if (!is.null(color) && color %in% names(df)) ggplot2::aes(.data[[x]], .data[[y]], color = .data[[color]]) else ggplot2::aes(.data[[x]], .data[[y]])
  p <- ggplot2::ggplot(df, aes) + ggplot2::geom_point(alpha = 0.75, size = 1.8) +
    ggplot2::theme_bw(base_size = 10) + ggplot2::labs(title = title, x = x, y = y, color = color %||% "")
  ggplot2::ggsave(path, p, width = 7.5, height = 5.2, dpi = 160)
  TRUE
}

ptpv3_make_figures <- function(dirs, role_table, original, exact, sim, crossfit, control_ref, cov_white, multidim, trajectory) {
  coverage <- original$primary_fit$metadata
  cov_plot <- aggregate(cell_count ~ sample_id + initial_ploidy + treatment + subbin_id, coverage, sum)
  p <- ggplot2::ggplot(cov_plot, ggplot2::aes(subbin_id, cell_count, fill = treatment)) +
    ggplot2::geom_col(position = "dodge") + ggplot2::facet_wrap(~ initial_ploidy) +
    ggplot2::theme_bw(base_size = 10) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
    ggplot2::labs(x = "Pseudotime subbin", y = "Cells", title = "Mouse/subbin/cell coverage")
  ggplot2::ggsave(file.path(dirs$figures, "fig01_coverage_common_support.png"), p, width = 8, height = 5, dpi = 160)

  ptpv3_plot_or_skip(role_table, file.path(dirs$figures, "fig02_program_size_effective_rank.png"),
                     "n_genes_observed", "n_positive_weight_genes", "tier_label", "Program size and direction metadata")

  score <- original$score_results
  primary <- score[score$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
  ptpv3_plot_or_skip(primary, file.path(dirs$figures, "fig03_score_model_forest_proxy.png"),
                     "estimate", "p_value", "tier_label", "Score model interaction estimates")

  ptpv3_plot_or_skip(sim, file.path(dirs$figures, "fig04_simulation_calibration.png"),
                     "true_effect", "power", "score_method_id", "Simulation calibration")

  exact_ad <- exact$adaptive
  ptpv3_plot_or_skip(exact_ad, file.path(dirs$figures, "fig05_exact_adaptive_pvalues.png"),
                     "observed_min_component_p", "p_value", "tier_label", "Exact permutation adaptive P values")

  controls <- primary[primary$tier == "negative_control", , drop = FALSE]
  ptpv3_plot_or_skip(controls, file.path(dirs$figures, "fig06_negative_control_distribution.png"),
                     "estimate", "p_value", "program_label", "Negative-control score evidence")

  if (!is.null(crossfit$fold_audit) && nrow(crossfit$fold_audit) > 0L) {
    cf <- crossfit$fold_audit
    cf$fold_index <- seq_len(nrow(cf))
    ptpv3_plot_or_skip(cf, file.path(dirs$figures, "fig07_crossfit_fold_nonzero_genes.png"),
                       "fold_index", "nonzero_genes", "status", "Cross-fit fold audit")
  }

  ptpv3_plot_or_skip(multidim, file.path(dirs$figures, "fig08_multidimensional_pvalues.png"),
                     "joint_statistic", "p_value", "tier_label", "Multidimensional pathway interactions")
  ptpv3_plot_or_skip(trajectory, file.path(dirs$figures, "fig09_trajectory_global_pvalues.png"),
                     "global_F", "p_value", "tier_label", "Trajectory-wide global interactions")

  evidence_row <- function(method, df) {
    if (is.null(df) || nrow(df) == 0L || !"p_value" %in% names(df) || !"tier_label" %in% names(df)) return(data.frame())
    data.frame(method = method, p_value = df$p_value, tier_label = df$tier_label, stringsAsFactors = FALSE)
  }
  method_rows <- list(
    evidence_row("score_models", primary),
    evidence_row("exact_permutation", exact_ad),
    evidence_row("control_reference", control_ref),
    evidence_row("covariance_whitened", cov_white),
    evidence_row("multidimensional", multidim),
    evidence_row("trajectory_global", trajectory)
  )
  mat <- ptpv3_bind_rows(method_rows)
  if (nrow(mat) > 0L) {
    mat$sig <- mat$p_value < 0.05
    agg <- aggregate(sig ~ method + tier_label, mat, function(x) sum(x, na.rm = TRUE))
    p <- ggplot2::ggplot(agg, ggplot2::aes(method, tier_label, fill = sig)) +
      ggplot2::geom_tile(color = "white") + ggplot2::geom_text(ggplot2::aes(label = sig), size = 3) +
      ggplot2::theme_bw(base_size = 10) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
      ggplot2::labs(x = NULL, y = NULL, fill = "nominal P<0.05", title = "Evidence matrix by method and tier")
    ggplot2::ggsave(file.path(dirs$figures, "fig10_evidence_matrix.png"), p, width = 8.5, height = 5.5, dpi = 160)
  }
}

ptpv3_html_table <- function(df, n = 20L) {
  if (is.null(df) || nrow(df) == 0L || ncol(df) == 0L) return("<p>No rows.</p>")
  df <- head(as.data.frame(df, stringsAsFactors = FALSE), n)
  df[] <- lapply(df, function(x) {
    if (is.numeric(x)) ifelse(is.finite(x), signif(x, 4), NA) else as.character(x)
  })
  header <- paste0("<tr>", paste0("<th>", ptp_html_escape(names(df)), "</th>", collapse = ""), "</tr>")
  rows <- apply(df, 1L, function(r) paste0("<tr>", paste0("<td>", ptp_html_escape(r), "</td>", collapse = ""), "</tr>"))
  paste0("<table>", header, paste(rows, collapse = "\n"), "</table>")
}

ptpv3_img <- function(path, title, legend, interpretation) {
  if (!file.exists(path)) return(paste0("<section><h3>", ptp_html_escape(title), "</h3><p>", ptp_html_escape(interpretation), "</p></section>"))
  src <- paste0("data:image/png;base64,", base64enc::base64encode(path, linewidth = 0L))
  paste0("<section><h3>", ptp_html_escape(title), "</h3><img src=\"", src, "\" alt=\"", ptp_html_escape(title),
         "\"/><p><b>Legend.</b> ", ptp_html_escape(legend), "</p><p><b>Interpretation.</b> ",
         ptp_html_escape(interpretation), "</p></section>")
}

ptpv3_result_classification <- function(result_list) {
  rows <- list()
  for (nm in names(result_list)) {
    df <- result_list[[nm]]
    if (is.null(df) || nrow(df) == 0L || !"p_value" %in% names(df)) next
    p <- df
    p$method <- nm
    fdr_col <- if ("fdr_focused_tier1_tier2" %in% names(p)) "fdr_focused_tier1_tier2" else "fdr_all_universe"
    p$class <- ifelse(is.finite(p[[fdr_col]]) & p[[fdr_col]] < 0.10, "significant_predefined_FDR",
                      ifelse(is.finite(p$p_value) & p$p_value < 0.05, "nominal_only", "non_significant"))
    rows[[length(rows) + 1L]] <- p[, intersect(c("method", "program_id", "program_label", "tier", "response_family", "p_value", fdr_col, "class"), names(p)), drop = FALSE]
  }
  ptpv3_bind_rows(rows)
}

ptpv3_write_html_report <- function(dirs, cfg, original, exact, sim, external_status, crossfit, control_ref, cov_white, multidim, trajectory, acceptance, result_classes) {
  primary <- original$score_results[original$score_results$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
  focused_sig <- sum(primary$fdr_focused_tier1_tier2 < 0.10, na.rm = TRUE)
  control_nominal <- sum(primary$tier == "negative_control" & primary$p_value < 0.05, na.rm = TRUE)
  conclusion <- if (focused_sig > 0L) {
    paste0(focused_sig, " focused score-model interactions pass focused FDR < 0.10; interpret only with method concordance and controls.")
  } else {
    "No focused primary score-model interaction passes focused FDR < 0.10; failure to reject is not evidence of ploidy independence."
  }
  figs <- list.files(dirs$figures, pattern = "\\.png$", full.names = TRUE)
  fig_html <- paste(vapply(seq_along(figs), function(i) {
    ptpv3_img(figs[[i]], paste0("Figure ", i, ". ", tools::file_path_sans_ext(basename(figs[[i]]))),
              "Generated from V3 result CSV/RDS objects in the same result tree.",
              "This figure is descriptive unless the corresponding section identifies a pre-defined FDR family or exact empirical test.")
  }, character(1L)), collapse = "\n")
  sections <- c(
    "<h2>1. Scientific question and estimands</h2><p>V3 tests whether the gemcitabine-associated transcriptional program response differs by baseline initial ploidy at matched pseudotime/state support. The primary estimand remains treatment by initial ploidy.</p>",
    "<h2>2. Frozen definitions and post-hoc status</h2><p>This is a post-hoc methodological analysis. Program universe, tiers, roles, eligibility, seeds, permutation space, FDR families, score methods and simulation scenarios were exported before V3 real-result interpretation.</p>",
    "<h2>3. Cohort/design/common-support QC</h2>", ptpv3_html_table(original$primary_fit$metadata, 12L),
    "<h2>4. Program universe, response families and tiers</h2>", ptpv3_html_table(readr::read_csv(file.path(dirs$frozen, "program_universe_role_eligibility.csv"), show_col_types = FALSE), 30L),
    "<h2>5. Method calibration and simulation</h2>", ptpv3_html_table(sim, 30L),
    "<h2>6. Tier 1 results by response family</h2>", ptpv3_html_table(primary[primary$tier == "tier1_core", ], 30L),
    "<h2>7. Tier 2 results by response family</h2>", ptpv3_html_table(primary[primary$tier == "tier2_upgraded_focused", ], 30L),
    "<h2>8. Tier 3 exploratory results by response family</h2>", ptpv3_html_table(primary[primary$tier == "tier3_exploratory_only", ], 30L),
    "<h2>9. Negative/context controls</h2><p>Nominal negative-control score-model hits: ", control_nominal, ".</p>", ptpv3_html_table(primary[primary$tier == "negative_control", ], 20L),
    "<h2>10. Score model comparison</h2>", ptpv3_html_table(original$score_results, 30L),
    "<h2>11. Exact permutation adaptive tests</h2>", ptpv3_html_table(exact$adaptive, 30L),
    "<h2>12. fry/mroast/camera/GSEA gene-set evidence</h2>", ptpv3_html_table(original$self_contained$fry, 20L), ptpv3_html_table(original$camera, 20L),
    "<h2>13. Cross-fitted 04i and external directional signatures</h2>", ptpv3_html_table(crossfit$fold_audit, 20L), ptpv3_html_table(external_status, 5L),
    "<h2>14. Multidimensional pathway results</h2>", ptpv3_html_table(multidim, 30L),
    "<h2>15. Trajectory-wide global interactions</h2>", ptpv3_html_table(trajectory, 30L),
    "<h2>16. ETP/end-time ploidy/dose sensitivity analyses</h2>", ptpv3_html_table(original$dose_program, 20L), ptpv3_html_table(original$endpoint$continuous_program, 20L),
    "<h2>17. Cross-method concordance and discordance</h2>", ptpv3_html_table(result_classes, 40L),
    "<h2>18. Integrated scientific conclusion</h2><p>", ptp_html_escape(conclusion), " Controls and simulation calibration are considered before any positive interpretation; nominal-only or single-method results are not confirmatory.</p>",
    "<h2>19. Limitations, not-estimable methods and audit information</h2>", ptpv3_html_table(acceptance, 50L),
    fig_html
  )
  html <- paste0("<!doctype html><html><head><meta charset='utf-8'><title>04j V3 report</title><style>",
                 "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:28px;line-height:1.45;color:#1f2937}table{border-collapse:collapse;font-size:12px;margin:12px 0;max-width:100%;display:block;overflow:auto}td,th{border:1px solid #d8dee9;padding:4px 6px}th{background:#edf2f7}img{max-width:100%;border:1px solid #d8dee9}h1,h2,h3{color:#111827}section{margin:20px 0}</style></head><body>",
                 "<h1>04j V3 pseudotime treatment-by-ploidy methodological analysis</h1>",
                 paste(sections, collapse = "\n"), "</body></html>")
  path <- file.path(dirs$report, "04j_pseudotime_treatment_ploidy_programs_v3_report.html")
  writeLines(html, path, useBytes = TRUE)
  path
}

ptpv3_acceptance_table <- function(cfg, before_after, role_table, exact, crossfit, external_status, sim, dirs) {
  rows <- list(
    data.frame(check = "v1_v2_checksum_unchanged", status = if (all(before_after$status == "unchanged")) "ok" else "failed",
               detail = paste(table(before_after$status), collapse = ";"), stringsAsFactors = FALSE),
    data.frame(check = "program_universe_frozen", status = if (file.exists(file.path(dirs$frozen, "program_universe_role_eligibility.csv"))) "ok" else "failed",
               detail = ptp_file_checksum(file.path(dirs$frozen, "program_universe_role_eligibility.csv")), stringsAsFactors = FALSE),
    data.frame(check = "exact_4900_unique", status = if (all(exact$validation$status == "ok")) "ok" else "failed",
               detail = paste(exact$validation$check, exact$validation$status, collapse = ";"), stringsAsFactors = FALSE),
    data.frame(check = "eligible_program_tier_complete", status = if (all(nzchar(role_table$tier[role_table$eligibility %in% TRUE]))) "ok" else "failed",
               detail = paste(sum(role_table$eligibility %in% TRUE), "eligible programs"), stringsAsFactors = FALSE),
    data.frame(check = "crossfit_excludes_heldout", status = if (!is.null(crossfit$fold_audit) && nrow(crossfit$fold_audit) > 0L && all(!crossfit$fold_audit$training_contains_heldout)) "ok" else "failed",
               detail = paste0("folds=", nrow(crossfit$fold_audit %||% data.frame())), stringsAsFactors = FALSE),
    data.frame(check = "broad_pathways_no_directional_claim", status = if (all(!role_table$directional_claim_eligible[role_table$analysis_role == "broad_pathway"])) "ok" else "failed",
               detail = paste(sum(role_table$analysis_role == "broad_pathway"), "broad pathways"), stringsAsFactors = FALSE),
    data.frame(check = "simulation_completed", status = if (nrow(sim) > 0L) "ok" else "failed", detail = paste(nrow(sim), "summary rows"), stringsAsFactors = FALSE),
    data.frame(check = "external_directional_signature", status = external_status$status[[1L]], detail = external_status$reason[[1L]], stringsAsFactors = FALSE),
    data.frame(check = "html_embedded_images", status = "pending", detail = "checked after report write", stringsAsFactors = FALSE)
  )
  ptpv3_bind_rows(rows)
}

ptpv3_copy_branch_outputs <- function(original, dirs) {
  ptp_write_csv(original$score_results, file.path(dirs$primary$score_models, "all_existing_score_methods_interactions.csv"))
  ptp_write_csv(original$score_precision, file.path(dirs$primary$score_models, "all_existing_score_methods_precision.csv"))
  ptp_write_csv(original$self_contained$fry, file.path(dirs$primary$fry_mroast, "fry_block_duplicateCorrelation.csv"))
  ptp_write_csv(original$self_contained$mroast_mean, file.path(dirs$primary$fry_mroast, "mroast_mean_block_duplicateCorrelation.csv"))
  ptp_write_csv(original$self_contained$mroast_msq, file.path(dirs$primary$fry_mroast, "mroast_msq_block_duplicateCorrelation.csv"))
  ptp_write_csv(original$camera, file.path(dirs$primary$camera, "cameraPR_rho_rank_grid.csv"))
  ptp_write_csv(original$dose_program, file.path(dirs$secondary_dose, "dose_specific_program_interactions.csv"))
  ptp_write_csv(original$endpoint$continuous_program, file.path(dirs$secondary_etp, "continuous_mean_endpoint_ploidy_program_interactions.csv"))
  ptp_write_csv(original$endpoint$group_program, file.path(dirs$secondary_end_time_ploidy, "threshold_endpoint_ploidy_program_interactions.csv"))
}

ptpv3_run_workflow <- function(args, repo_root, script_dir) {
  ptpv3_thread_env()
  set.seed(args$seed)
  ptp_check_packages()
  pst_check_packages()
  cfg <- ptpv3_read_config(args$config)
  output_root <- ptp_clean_dir(args$output_root, overwrite = args$overwrite)
  dirs <- ptpv3_dirs(output_root)

  ptp_write_csv(data.frame(parameter = names(args), value = vapply(args, as.character, character(1L)), stringsAsFactors = FALSE),
                file.path(dirs$manifest, "analysis_parameters.csv"))
  ptp_write_csv(ptp_package_versions(), file.path(dirs$manifest, "package_versions.csv"))
  file.copy(args$config, file.path(dirs$frozen, basename(args$config)), overwrite = TRUE)

  legacy_files <- file.path(repo_root, c(
    "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs.R",
    "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_util.R",
    "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_config.yaml",
    "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_v2.R",
    "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_v2_util.R",
    "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_v2_config.yaml",
    "Data/in-vivo/Plan/pseudotime_treatment_by_ploidy_matched_state_plan.md"
  ))
  code_before <- data.frame(label = "v1_v2_code_plan", phase = "before", path = legacy_files,
                            relative_path = sub(paste0("^", repo_root, "/"), "", legacy_files),
                            exists = file.exists(legacy_files), is_dir = FALSE,
                            bytes = suppressWarnings(file.info(legacy_files)$size),
                            sha256 = vapply(legacy_files, ptp_file_checksum, character(1L)), stringsAsFactors = FALSE)
  ptp_write_csv(code_before, file.path(dirs$manifest, "v1_v2_code_plan_sha256_before.csv"))
  v1_before <- ptpv3_checksum_tree(cfg$analysis$v1_result_root, "v1_results", "before", dirs)
  v2_before <- ptpv3_checksum_tree(cfg$analysis$v2_result_root, "v2_results", "before", dirs)

  cell_metadata_path <- normalizePath(args$cell_metadata, winslash = "/", mustWork = TRUE)
  sample_info_path <- if (file.exists(args$sample_info)) normalizePath(args$sample_info, winslash = "/", mustWork = TRUE) else args$sample_info
  seurat_rds <- ptp_resolve_seurat_rds(args$seurat_rds, cfg, args$results_root)
  ptp_write_csv(data.frame(input = c("cell_metadata", "sample_info", "seurat_rds", "config", "base_v2_config", "plan"),
                           path = c(cell_metadata_path, sample_info_path, seurat_rds, args$config, cfg$v3$analysis$base_v2_config, cfg$analysis$plan_path),
                           sha256 = c(ptp_file_checksum(cell_metadata_path), ptp_file_checksum(sample_info_path), ptp_file_checksum(seurat_rds), ptp_file_checksum(args$config), ptp_file_checksum(cfg$v3$analysis$base_v2_config), ptp_file_checksum(cfg$analysis$plan_path)),
                           stringsAsFactors = FALSE), file.path(dirs$manifest, "input_manifest.csv"))

  meta <- ptp_read_cell_metadata(cell_metadata_path)
  loaded <- ptp_load_counts_and_qc(seurat_rds, args$assay, args$counts_layer, meta)
  counts <- loaded$counts
  meta <- loaded$metadata
  ptp_write_csv(attr(counts, "seurat_audit"), file.path(dirs$qc, "expression_source_audit.csv"))
  ptp_write_csv(ptp_audit_counts(counts), file.path(dirs$qc, "count_matrix_audit.csv"))
  matched <- ptp_match_cells(meta, counts, cfg$qc$min_metadata_to_counts_match_rate, dirs$qc)
  meta <- matched$meta
  counts <- matched$counts
  limited <- ptp_limit_for_smoke(meta, counts, args$max_cells, args$max_genes, args$seed)
  meta <- limited$meta
  counts <- limited$counts
  ptp_write_csv(ptp_assignment_audit(meta), file.path(dirs$qc, "sample_assignment_audit.csv"))
  ptp_write_csv(ptp_sample_metadata(meta), file.path(dirs$qc, "mouse_region_coverage.csv"))
  if (file.exists(sample_info_path)) {
    crosswalk <- ptp_processing_batch_crosswalk(meta, sample_info_path, ptp_endpoint_ploidy_specs(cfg))
    ptp_write_csv(crosswalk, file.path(dirs$qc, "processing_batch_crosswalk.csv"))
    ptp_write_csv(ptp_batch_group_separability(crosswalk), file.path(dirs$qc, "processing_batch_separability_audit.csv"))
  }

  expanded <- ptpv2_expand_program_config(cfg)
  cfg2 <- expanded$cfg
  pb0 <- ptp_choose_primary_grid(meta, counts, cfg2)
  if (!isTRUE(pb0$estimable)) stop("Primary matched pseudotime grid is not estimable.", call. = FALSE)
  fit0 <- ptp_fit_cellmeans(pb0$counts, pb0$metadata, pb0$eligible_subbins, "ptp_group", cfg2$qc)
  observed_symbols <- unique(ptp_clean_gene_symbols(fit0$retained_genes))
  programs <- ptp_load_programs(cfg2, cfg2$qc$min_program_genes_observed, observed_symbols)
  role_table <- ptpv3_program_role_table(programs)
  ptp_write_csv(role_table, file.path(dirs$frozen, "program_universe_role_eligibility.csv"))
  ptp_write_csv(ptp_program_weight_table(programs), file.path(dirs$frozen, "program_score_gene_weights.csv"))
  ptp_write_csv(ptp_subbins_table(cfg2$primary_subbins, "primary_four_bin"), file.path(dirs$frozen, "primary_subbin_definitions.csv"))
  ptp_write_csv(ptp_subbins_table(cfg2$fallback_subbins, "fallback_two_bin"), file.path(dirs$frozen, "fallback_subbin_definitions.csv"))
  ptp_write_csv(ptp_regions_table(cfg2), file.path(dirs$frozen, "pseudotime_region_definitions.csv"))
  ptp_write_csv(ptp_endpoint_ploidy_specs_table(ptp_endpoint_ploidy_specs(cfg2)), file.path(dirs$frozen, "endpoint_ploidy_definitions.csv"))
  ptp_write_csv(data.frame(seed_name = names(cfg$v3$reproducibility), value = unlist(cfg$v3$reproducibility), stringsAsFactors = FALSE),
                file.path(dirs$frozen, "random_seeds_and_threads.csv"))
  ptpv3_assert_complete_program_meta(role_table, role_table, "frozen_program_universe", dirs)

  message("Running V3 primary branch checkpoint.")
  branch_role_table <- ptpv3_branch_role_table(role_table)
  original <- ptpv2_run_branch("original_coordinate", file.path(dirs$checkpoints, "original_coordinate_branch"),
                               meta, counts, cfg2, args, programs, branch_role_table, include_heavy = TRUE)
  original$score_results <- ptpv3_apply_fdr_family(original$score_results, "p_value", "score_method_id", role_table)
  original$self_contained$fry <- ptpv3_apply_fdr_family(original$self_contained$fry, if ("PValue" %in% names(original$self_contained$fry)) "PValue" else "P.Value", "method_id", role_table)
  original$self_contained$mroast_mean <- ptpv3_apply_fdr_family(original$self_contained$mroast_mean, "PValue", "method_id", role_table)
  original$self_contained$mroast_msq <- ptpv3_apply_fdr_family(original$self_contained$mroast_msq, if ("PValue.Mixed" %in% names(original$self_contained$mroast_msq)) "PValue.Mixed" else "PValue", "method_id", role_table)
  original$camera <- ptpv3_apply_fdr_family(original$camera, "PValue", "method_id", role_table)
  original$dose_program <- ptpv3_apply_fdr_family(original$dose_program, "p_value", "contrast_id", role_table)
  original$endpoint$group_program <- ptpv3_apply_fdr_family(original$endpoint$group_program, "p_value", "endpoint_ploidy_method", role_table)
  original$endpoint$continuous_program <- ptpv3_apply_fdr_family(original$endpoint$continuous_program, "p_value", "contrast_id", role_table)
  ptpv3_copy_branch_outputs(original, dirs)
  ptpv3_assert_complete_program_meta(original$score_results, role_table, "score_models", dirs)
  ptpv3_assert_complete_program_meta(original$self_contained$fry, role_table, "fry", dirs)
  ptpv3_assert_complete_program_meta(original$self_contained$mroast_msq, role_table, "mroast_msq", dirs)
  ptpv3_assert_complete_program_meta(original$camera, role_table, "camera", dirs)

  external_status <- ptpv3_external_signature_status(cfg2, repo_root, dirs)
  exact <- if (isTRUE(args$run_exact_permutation)) ptpv3_exact_permutation(original$primary_fit, programs, role_table, cfg2, args, dirs, pb0$eligible_subbins) else list(adaptive = data.frame(), validation = data.frame())
  control_ref <- if (isTRUE(args$run_control_reference)) ptpv3_control_reference_residual(original$primary_fit, programs, role_table, cfg2, dirs) else data.frame()
  cov_white <- if (isTRUE(args$run_covariance_whitened)) ptpv3_covariance_whitened(original$primary_fit, programs, role_table, cfg2, dirs) else data.frame()
  crossfit <- if (isTRUE(args$run_crossfit_04i)) ptpv3_crossfit_04i(original$primary_fit, programs, role_table, cfg2, args, repo_root, dirs) else list(status = data.frame(status = "skipped"), fold_audit = data.frame(), results = data.frame())
  multidim <- if (isTRUE(args$run_multidimensional) && nrow(exact$permutations) > 0L) ptpv3_multidimensional(original$primary_fit, programs, role_table, exact, cfg2, dirs, pb0$eligible_subbins) else data.frame()
  trajectory <- if (isTRUE(args$run_trajectory_global)) ptpv3_trajectory_global(meta, counts, programs, role_table, cfg2, args, dirs) else data.frame()
  sim <- if (isTRUE(args$run_simulation)) ptpv3_run_simulation(original$primary_fit, programs, role_table, cfg2, args, dirs) else data.frame()

  result_classes <- ptpv3_result_classification(list(
    score_models = original$score_results[original$score_results$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE],
    exact_permutation = exact$adaptive,
    control_reference = control_ref[control_ref$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE],
    covariance_whitened = cov_white[cov_white$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE],
    multidimensional = multidim,
    trajectory_global = trajectory
  ))
  ptp_write_csv(result_classes, file.path(dirs$tables, "result_classification_by_method.csv"))

  code_after <- data.frame(label = "v1_v2_code_plan", phase = "after", path = legacy_files,
                           relative_path = sub(paste0("^", repo_root, "/"), "", legacy_files),
                           exists = file.exists(legacy_files), is_dir = FALSE,
                           bytes = suppressWarnings(file.info(legacy_files)$size),
                           sha256 = vapply(legacy_files, ptp_file_checksum, character(1L)), stringsAsFactors = FALSE)
  ptp_write_csv(code_after, file.path(dirs$manifest, "v1_v2_code_plan_sha256_after.csv"))
  v1_after <- ptpv3_checksum_tree(cfg$analysis$v1_result_root, "v1_results", "after", dirs)
  v2_after <- ptpv3_checksum_tree(cfg$analysis$v2_result_root, "v2_results", "after", dirs)
  comp_code <- ptpv3_compare_checksums(code_before, code_after, "v1_v2_code_plan", dirs)
  comp_v1 <- ptpv3_compare_checksums(v1_before, v1_after, "v1_results", dirs)
  comp_v2 <- ptpv3_compare_checksums(v2_before, v2_after, "v2_results", dirs)
  before_after <- rbind(comp_code, comp_v1, comp_v2)

  acceptance <- ptpv3_acceptance_table(cfg2, before_after, role_table, exact, crossfit, external_status, sim, dirs)
  ptp_write_csv(acceptance, file.path(dirs$manifest, "v3_acceptance_checks.csv"))
  ptpv3_make_figures(dirs, role_table, original, exact, sim, crossfit, control_ref, cov_white, multidim, trajectory)
  report <- ptpv3_write_html_report(dirs, cfg2, original, exact, sim, external_status, crossfit, control_ref, cov_white, multidim, trajectory, acceptance, result_classes)
  html <- readLines(report, warn = FALSE)
  acceptance$status[acceptance$check == "html_embedded_images"] <- if (any(grepl("data:image/png;base64", html, fixed = TRUE))) "ok" else "failed"
  ptp_write_csv(acceptance, file.path(dirs$manifest, "v3_acceptance_checks.csv"))

  session <- utils::capture.output(sessionInfo())
  writeLines(session, file.path(dirs$manifest, "sessionInfo.txt"), useBytes = TRUE)
  message("04j V3 workflow complete: ", output_root)
  invisible(output_root)
}
