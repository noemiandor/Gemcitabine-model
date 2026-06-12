run_enrichment_or_stop <- function(values,
                                   annotations,
                                   cancer,
                                   direction,
                                   permute_n,
                                   pvalue_cutoff,
                                   enrichment_fn = enrichment) {
  tryCatch(
    enrichment_fn(
      values,
      annotations,
      permute.n = permute_n,
      normalize = FALSE,
      pvalue.cutoff = pvalue_cutoff
    )$pvalue,
    error = function(e) {
      stop(
        sprintf(
          "Enrichment failed for cancer=%s direction=%s drugs=%d groups=%d: %s",
          cancer,
          direction,
          length(values),
          length(unique(annotations$group)),
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )
}

plot_barplot_or_stop <- function(plot_values,
                                 colors,
                                 cancer,
                                 xlab,
                                 plotting_fn = barplot) {
  tryCatch(
    plotting_fn(
      plot_values,
      col = colors,
      main = cancer,
      horiz = TRUE,
      las = 2,
      cex.lab = 0.7,
      cex.names = 0.35,
      xlab = xlab
    ),
    error = function(e) {
      stop(
        sprintf(
          "Plotting failed for cancer=%s plotted_drugs=%d: %s",
          cancer,
          length(plot_values),
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )
}

normalize_cell_line_name <- function(x) {
  x <- toupper(as.character(x))
  gsub("[^A-Z0-9]", "", x)
}

audit_normalized_key_collisions <- function(dt, raw_col, key_col, source_name) {
  missing_cols <- setdiff(c(raw_col, key_col), colnames(dt))
  if (length(missing_cols) > 0) {
    stop("Collision audit input is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  tmp <- unique(as.data.frame(dt)[, c(raw_col, key_col), drop = FALSE])
  names(tmp) <- c("raw_name", "normalized_key")
  tmp <- tmp[!is.na(tmp$normalized_key) & nzchar(tmp$normalized_key), , drop = FALSE]

  out <- data.table::as.data.table(tmp)
  out <- out[, .(
    n_raw_names = data.table::uniqueN(raw_name),
    raw_names = paste(sort(unique(raw_name)), collapse = " | ")
  ), by = normalized_key]
  out <- out[n_raw_names > 1]
  out[, source := source_name]
  out[]
}

validate_ploidy_key_collisions <- function(appCL,
                                           key_col = "CELL_LINE_KEY",
                                           raw_col = "Cell iname",
                                           ploidy_col = "ploidy",
                                           tolerance = 0.05,
                                           tables_dir = NULL) {
  missing_cols <- setdiff(c(key_col, raw_col, ploidy_col), colnames(appCL))
  if (length(missing_cols) > 0) {
    stop("Ploidy collision validation is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  appCL_dt <- data.table::as.data.table(appCL)
  summary_dt <- appCL_dt[, .(
    n_rows = .N,
    n_raw_names = data.table::uniqueN(get(raw_col)),
    ploidy_min = min(get(ploidy_col), na.rm = TRUE),
    ploidy_max = max(get(ploidy_col), na.rm = TRUE),
    ploidy_range = max(get(ploidy_col), na.rm = TRUE) - min(get(ploidy_col), na.rm = TRUE),
    raw_names = paste(sort(unique(get(raw_col))), collapse = " | ")
  ), by = key_col]
  data.table::setnames(summary_dt, key_col, "normalized_key")

  collapsed <- summary_dt[n_raw_names > 1 & ploidy_range <= tolerance]
  failures <- summary_dt[n_raw_names > 1 & ploidy_range > tolerance]

  if (!is.null(tables_dir)) {
    write_tsv(as.data.frame(summary_dt), file.path(tables_dir, "cell_line_key_ploidy_summary.tsv"))
    write_tsv(as.data.frame(collapsed), file.path(tables_dir, "cell_line_key_ploidy_near_duplicate_collapses.tsv"))
    if (nrow(failures) > 0) {
      write_tsv(as.data.frame(failures), file.path(tables_dir, "cell_line_key_collisions_ploidy_FAIL.tsv"))
    }
  }

  if (nrow(failures) > 0) {
    stop(
      "Normalized ploidy cell-line keys are not unique and have discordant ploidy values. Review cell_line_key_collisions_ploidy_FAIL.tsv.",
      call. = FALSE
    )
  }

  invisible(list(summary = summary_dt, collapsed = collapsed, failures = failures))
}

write_cell_line_matching_delta <- function(dr, appCL, path) {
  raw_names <- if ("CELL_LINE_NAME_RAW" %in% colnames(dr)) dr$CELL_LINE_NAME_RAW else dr$CELL_LINE_NAME
  gdsc <- unique(data.frame(
    gdsc_cell_line_name = dr$CELL_LINE_NAME,
    gdsc_cell_line_name_raw = raw_names,
    normalized_key = dr$CELL_LINE_KEY,
    stringsAsFactors = FALSE
  ))

  ploidy <- unique(data.frame(
    ploidy_cell_line_name = appCL$`Cell iname`,
    normalized_key = appCL$CELL_LINE_KEY,
    stringsAsFactors = FALSE
  ))

  ploidy_by_key <- data.table::as.data.table(ploidy)
  ploidy_by_key <- ploidy_by_key[, .(
    normalized_match_count = .N,
    normalized_match_raw_names = paste(sort(unique(ploidy_cell_line_name)), collapse = " | ")
  ), by = normalized_key]

  delta <- data.table::as.data.table(gdsc)
  delta[, canonical_raw_match := gdsc_cell_line_name %in% appCL$`Cell iname`]
  delta <- merge(delta, ploidy_by_key, by = "normalized_key", all.x = TRUE, sort = FALSE)
  delta[, normalized_key_match := !is.na(normalized_match_count)]
  delta[is.na(normalized_match_count), normalized_match_count := 0L]
  delta[is.na(normalized_match_raw_names), normalized_match_raw_names := ""]
  delta[, match_status := data.table::fifelse(
    canonical_raw_match & normalized_key_match,
    "matched_by_both",
    data.table::fifelse(
      canonical_raw_match & !normalized_key_match,
      "raw_only",
      data.table::fifelse(!canonical_raw_match & normalized_key_match, "normalized_only", "unmatched")
    )
  )]
  delta <- delta[order(match_status, gdsc_cell_line_name)]
  write_tsv(as.data.frame(delta), path)
  invisible(delta)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

validate_required_groups <- function(groups, required_groups = c("SIGNALING", "CYTOTOXIC")) {
  missing_groups <- setdiff(required_groups, unique(groups))
  if (length(missing_groups) > 0) {
    stop(
      "Missing required annotation groups after filtering: ",
      paste(missing_groups, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

resolve_duplicate_drug_cell_lines <- function(drug_table,
                                             metric,
                                             strategy = "lowest_rmse",
                                             qc_dir = NULL) {
  required_cols <- c("DRUG_NAME", "CELL_LINE_NAME", metric)
  missing_cols <- setdiff(required_cols, colnames(drug_table))
  if (length(missing_cols) > 0) {
    stop("Drug table is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (!strategy %in% c("lowest_rmse", "first_input_order")) {
    stop("Unsupported duplicate resolution strategy: ", strategy, call. = FALSE)
  }

  key_cols <- c("DRUG_NAME", "CELL_LINE_NAME")
  duplicate_mask <- duplicated(drug_table[, key_cols]) | duplicated(drug_table[, key_cols], fromLast = TRUE)
  duplicate_rows <- drug_table[duplicate_mask, , drop = FALSE]

  if (!is.null(qc_dir)) {
    dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
    write_tsv(duplicate_rows, file.path(qc_dir, "duplicate_drug_cell_line_records.tsv"))
  }

  if (strategy == "first_input_order") {
    resolved <- drug_table[!duplicated(drug_table[, key_cols]), , drop = FALSE]
  } else {
    order_cols <- list(drug_table$DRUG_NAME, drug_table$CELL_LINE_NAME)
    if ("RMSE" %in% colnames(drug_table)) {
      order_cols <- c(order_cols, list(is.na(drug_table$RMSE), drug_table$RMSE))
    }
    for (col in c("NLME_RESULT_ID", "NLME_CURVE_ID", "DRUG_ID", "COSMIC_ID")) {
      if (col %in% colnames(drug_table)) {
        order_cols <- c(order_cols, list(drug_table[[col]]))
      }
    }
    ordered <- drug_table[do.call(order, order_cols), , drop = FALSE]
    resolved <- ordered[!duplicated(ordered[, key_cols]), , drop = FALSE]
  }

  remaining_duplicates <- duplicated(resolved[, key_cols]) | duplicated(resolved[, key_cols], fromLast = TRUE)
  if (any(remaining_duplicates)) {
    stop("Duplicate resolution failed to produce one row per DRUG_NAME and CELL_LINE_NAME.", call. = FALSE)
  }

  summary <- data.frame(
    metric = metric,
    strategy = strategy,
    input_rows = nrow(drug_table),
    duplicate_rows = nrow(duplicate_rows),
    duplicate_keys = if (nrow(duplicate_rows) == 0) 0 else nrow(unique(duplicate_rows[, key_cols])),
    resolved_rows = nrow(resolved),
    removed_rows = nrow(drug_table) - nrow(resolved),
    stringsAsFactors = FALSE
  )
  if (!is.null(qc_dir)) {
    write_tsv(summary, file.path(qc_dir, "duplicate_resolution_summary.tsv"))
  }

  attr(resolved, "duplicate_resolution_summary") <- summary
  resolved
}
