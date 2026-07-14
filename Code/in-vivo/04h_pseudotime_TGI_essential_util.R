essential_parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (grepl("^--[^=]+=", token)) {
      key <- sub("^--([^=]+)=.*$", "\\1", token)
      out[[key]] <- sub("^--[^=]+=", "", token)
      i <- i + 1L
    } else if (grepl("^--", token)) {
      key <- sub("^--", "", token)
      if (i < length(args) && !grepl("^--", args[[i + 1L]])) {
        out[[key]] <- args[[i + 1L]]
        i <- i + 2L
      } else {
        out[[key]] <- "TRUE"
        i <- i + 1L
      }
    } else {
      i <- i + 1L
    }
  }
  out
}

essential_arg <- function(args, name, default = NULL) {
  value <- args[[name]]
  if (!is.null(value) && length(value) == 1L && nzchar(value)) value else default
}

essential_bool <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

essential_safe_numeric <- function(x) suppressWarnings(as.numeric(x))

essential_ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

essential_write_csv <- function(x, path) {
  essential_ensure_dir(dirname(path))
  write.csv(x, path, row.names = FALSE, na = "NA")
  invisible(path)
}

essential_clean_character <- function(x) {
  out <- trimws(as.character(x))
  out[out %in% c("", "NA", "NaN", "NULL", "None")] <- NA_character_
  out
}

essential_first_existing_column <- function(data, candidates, label) {
  matches <- candidates[candidates %in% names(data)]
  if (length(matches) == 0L) {
    stop("Cannot find ", label, ". Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  matches[[1L]]
}

essential_resolve_column_case_insensitive <- function(data, candidates) {
  lower_names <- tolower(names(data))
  for (candidate in candidates) {
    index <- match(tolower(candidate), lower_names)
    if (!is.na(index)) return(names(data)[[index]])
  }
  NA_character_
}

essential_standardize_ploidy <- function(x) {
  out <- toupper(essential_clean_character(x))
  out <- gsub("\\s+", "", out)
  out[out %in% c("2", "2.0", "2N")] <- "2N"
  out[out %in% c("4", "4.0", "4N")] <- "4N"
  out[!(out %in% c("2N", "4N"))] <- NA_character_
  out
}

essential_standardize_dose <- function(x) {
  out <- essential_clean_character(x)
  out[out %in% c("0", "0mg", "0 mg/kg", "0mg/kg")] <- "0mg/kg"
  out[out %in% c("30", "30mg", "30 mg/kg", "30mg/kg")] <- "30mg/kg"
  out[out %in% c("120", "120mg", "120 mg/kg", "120mg/kg")] <- "120mg/kg"
  out
}

essential_dose_panel <- function(x) {
  out <- sub("mg/kg$", "", essential_standardize_dose(x))
  out[!(out %in% c("0", "30", "120"))] <- NA_character_
  out
}

essential_extract_barcode <- function(cell, sample) {
  cell <- essential_clean_character(cell)
  sample <- essential_clean_character(sample)
  barcode <- cell
  has_sample <- !is.na(cell) & !is.na(sample) & nzchar(sample) & startsWith(cell, paste0(sample, "_"))
  barcode[has_sample] <- substring(cell[has_sample], nchar(sample[has_sample]) + 2L)
  fallback <- !has_sample & !is.na(cell) & grepl("_", cell)
  barcode[fallback] <- sub("^[^_]+_", "", cell[fallback])
  barcode
}

essential_require_input_generation_packages <- function() {
  required <- c("readr", "readxl")
  missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Automatic input generation requires R package(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

essential_read_sample_info_map <- function(path) {
  if (!file.exists(path)) stop("Missing sample information workbook: ", path, call. = FALSE)
  raw <- as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE)
  sample_column <- essential_first_existing_column(
    raw,
    c("IDs", "sample", "sample_id", "sampleID"),
    "sample ID column"
  )
  harvest_column <- essential_first_existing_column(raw, c("harvest"), "harvest column")
  dose_matches <- c("Dose", "dose")
  dose_matches <- dose_matches[dose_matches %in% names(raw)]
  dose <- if (length(dose_matches) > 0L) raw[[dose_matches[[1L]]]] else NA_character_
  out <- data.frame(
    sample = essential_clean_character(raw[[sample_column]]),
    harvest = essential_clean_character(raw[[harvest_column]]),
    sample_info_dose = essential_standardize_dose(dose),
    stringsAsFactors = FALSE
  )
  keep <- !is.na(out$sample) & nzchar(out$sample) & !is.na(out$harvest) & nzchar(out$harvest)
  out <- out[keep & !duplicated(out$sample), , drop = FALSE]
  rownames(out) <- NULL
  out
}

essential_read_cell_ploidy_map <- function(path) {
  if (!file.exists(path)) stop("Missing cell ploidy table: ", path, call. = FALSE)
  raw <- as.data.frame(
    readr::read_tsv(
      path,
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE,
      progress = FALSE
    ),
    stringsAsFactors = FALSE
  )
  required <- c("file", "cell_id", "ploidy")
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0L) {
    stop("Cell ploidy table is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  data <- data.frame(
    harvest = sub("\\.sps\\.cbs$", "", essential_clean_character(raw$file)),
    barcode_raw = essential_clean_character(raw$cell_id),
    cell_ploidy = essential_safe_numeric(raw$ploidy),
    ploidy_frac_covered = if ("frac_covered" %in% names(raw)) {
      essential_safe_numeric(raw$frac_covered)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
  keep <- !is.na(data$harvest) & nzchar(data$harvest) &
    !is.na(data$barcode_raw) & nzchar(data$barcode_raw)
  data <- data[keep, , drop = FALSE]
  key <- paste(data$harvest, data$barcode_raw, sep = "\r")
  key_levels <- unique(key)
  group <- match(key, key_levels)
  first_index <- match(seq_along(key_levels), group)
  ploidy_finite <- is.finite(data$cell_ploidy)
  coverage_finite <- is.finite(data$ploidy_frac_covered)
  ploidy_sum <- as.numeric(rowsum(ifelse(ploidy_finite, data$cell_ploidy, 0), group, reorder = FALSE))
  ploidy_n <- as.numeric(rowsum(as.integer(ploidy_finite), group, reorder = FALSE))
  coverage_sum <- as.numeric(rowsum(ifelse(coverage_finite, data$ploidy_frac_covered, 0), group, reorder = FALSE))
  coverage_n <- as.numeric(rowsum(as.integer(coverage_finite), group, reorder = FALSE))
  data.frame(
    harvest = data$harvest[first_index],
    barcode_raw = data$barcode_raw[first_index],
    cell_ploidy = ifelse(ploidy_n > 0, ploidy_sum / ploidy_n, NA_real_),
    ploidy_frac_covered = ifelse(coverage_n > 0, coverage_sum / coverage_n, NA_real_),
    n_ploidy_rows = tabulate(group, nbins = length(key_levels)),
    stringsAsFactors = FALSE
  )
}

essential_build_cell_metrics <- function(scvelo_path, cell_ploidy_path, sample_info_path) {
  required_files <- c(scvelo_path, cell_ploidy_path, sample_info_path)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    stop("Missing automatic-generation input(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
  }
  raw <- as.data.frame(
    readr::read_csv(
      scvelo_path,
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE,
      progress = FALSE
    ),
    stringsAsFactors = FALSE
  )
  required <- c("cell", "TN", "clusters", "sample", "Ploidy", "Dose", "velocity_pseudotime")
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0L) {
    stop("scVelo metrics are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  metrics <- data.frame(
    cell = essential_clean_character(raw$cell),
    TN = essential_clean_character(raw$TN),
    clusters = essential_clean_character(raw$clusters),
    sample = essential_clean_character(raw$sample),
    Ploidy = essential_standardize_ploidy(raw$Ploidy),
    Dose = essential_standardize_dose(raw$Dose),
    pseudotime = essential_safe_numeric(raw$velocity_pseudotime),
    stringsAsFactors = FALSE
  )
  metrics$barcode_raw <- essential_extract_barcode(metrics$cell, metrics$sample)

  sample_map <- essential_read_sample_info_map(sample_info_path)
  sample_index <- match(metrics$sample, sample_map$sample)
  metrics$harvest <- sample_map$harvest[sample_index]
  sample_dose <- sample_map$sample_info_dose[sample_index]
  missing_dose <- is.na(metrics$Dose) | !nzchar(metrics$Dose)
  metrics$Dose[missing_dose] <- sample_dose[missing_dose]

  ploidy_map <- essential_read_cell_ploidy_map(cell_ploidy_path)
  metric_key <- paste(metrics$harvest, metrics$barcode_raw, sep = "\r")
  ploidy_key <- paste(ploidy_map$harvest, ploidy_map$barcode_raw, sep = "\r")
  ploidy_index <- match(metric_key, ploidy_key)
  metrics$cell_ploidy <- ploidy_map$cell_ploidy[ploidy_index]
  metrics
}

essential_normalize_sample_id <- function(x) {
  out <- essential_clean_character(x)
  out <- sub("-Count-HM$", "", out)
  out <- sub("-Count$", "", out)
  out <- sub("-HM$", "", out)
  out
}

essential_infer_initial_ploidy <- function(sample, harvest = NULL) {
  sample <- essential_clean_character(sample)
  harvest <- if (is.null(harvest)) rep(NA_character_, length(sample)) else essential_clean_character(harvest)
  signal <- paste(sample, harvest)
  out <- rep(NA_character_, length(sample))
  out[grepl("(^|-)2N($|-)", signal, ignore.case = TRUE) | grepl("^2N($|-)", sample, ignore.case = TRUE)] <- "2N"
  out[
    grepl("(^|-)4N($|-)", signal, ignore.case = TRUE) |
      grepl("^4N($|-)", sample, ignore.case = TRUE) |
      grepl("^A[0-9]+-4N($|-)", sample, ignore.case = TRUE)
  ] <- "4N"
  out
}

essential_infer_dose_from_harvest <- function(harvest) {
  harvest <- essential_clean_character(harvest)
  dose <- rep(NA_character_, length(harvest))
  matched <- grepl("^SUM159-(2N|4N)-[0-9]+-", harvest, ignore.case = TRUE)
  dose[matched] <- sub(
    "^SUM159-(2N|4N)-([0-9]+)-.*$",
    "\\2",
    harvest[matched],
    ignore.case = TRUE
  )
  essential_standardize_dose(dose)
}

essential_calculate_tgi <- function(
  growth_curve_path,
  control_dose = "0mg/kg",
  baseline_day = "Day_0",
  final_day = NULL
) {
  if (!file.exists(growth_curve_path)) stop("Missing growth curve workbook: ", growth_curve_path, call. = FALSE)
  growth <- as.data.frame(readxl::read_excel(growth_curve_path, sheet = 1), stringsAsFactors = FALSE)
  harvest_column <- essential_resolve_column_case_insensitive(growth, c("harvest"))
  sample_column <- essential_resolve_column_case_insensitive(
    growth,
    c("Sequencing IDs", "Sequencing.IDs", "Sequencing ID", "SequencingIDs")
  )
  if (is.na(harvest_column)) stop("Growth curve workbook is missing a harvest column", call. = FALSE)
  if (is.na(sample_column)) stop("Growth curve workbook is missing a Sequencing IDs column", call. = FALSE)

  day_columns <- grep("^Day_[0-9]+$", names(growth), value = TRUE)
  day_numbers <- essential_safe_numeric(sub("^Day_", "", day_columns))
  day_columns <- day_columns[order(day_numbers)]
  day_numbers <- essential_safe_numeric(sub("^Day_", "", day_columns))
  if (!(baseline_day %in% day_columns)) stop("Missing baseline day: ", baseline_day, call. = FALSE)
  if (!is.null(final_day) && !(final_day %in% day_columns)) stop("Missing final day: ", final_day, call. = FALSE)
  baseline_time <- essential_safe_numeric(sub("^Day_", "", baseline_day))
  endpoint_columns <- day_columns[day_numbers > baseline_time]

  rows <- lapply(seq_len(nrow(growth)), function(index) {
    values <- essential_safe_numeric(unlist(growth[index, day_columns, drop = FALSE], use.names = FALSE))
    names(values) <- day_columns
    baseline_value <- values[[baseline_day]]
    local_final_day <- final_day
    if (is.null(local_final_day)) {
      finite_index <- which(is.finite(values))
      local_final_day <- if (length(finite_index) == 0L) NA_character_ else day_columns[finite_index[[length(finite_index)]]]
    }
    final_value <- if (!is.na(local_final_day)) values[[local_final_day]] else NA_real_
    finite_auc <- is.finite(day_numbers) & is.finite(values) & is.finite(baseline_value)
    auc_raw <- NA_real_
    auc_delta <- NA_real_
    auc_start_day <- NA_character_
    auc_end_day <- NA_character_
    if (sum(finite_auc) >= 2L) {
      time <- day_numbers[finite_auc]
      volume <- values[finite_auc]
      order_index <- order(time)
      time <- time[order_index]
      volume <- volume[order_index]
      auc_raw <- sum(diff(time) * (head(volume, -1L) + tail(volume, -1L)) / 2)
      volume_delta <- volume - baseline_value
      auc_delta <- sum(diff(time) * (head(volume_delta, -1L) + tail(volume_delta, -1L)) / 2)
      auc_start_day <- paste0("Day_", time[[1L]])
      auc_end_day <- paste0("Day_", time[[length(time)]])
    }
    row <- data.frame(
      growth_curve_row = index,
      sample_id_raw = as.character(growth[[sample_column]][index]),
      sample_id = essential_normalize_sample_id(growth[[sample_column]][index]),
      harvest = as.character(growth[[harvest_column]][index]),
      initial_ploidy = essential_infer_initial_ploidy(
        growth[[sample_column]][index],
        growth[[harvest_column]][index]
      ),
      dose = essential_infer_dose_from_harvest(growth[[harvest_column]][index]),
      baseline_day = baseline_day,
      baseline_volume = baseline_value,
      final_day = local_final_day,
      final_volume = final_value,
      tumor_volume_delta = final_value - baseline_value,
      tumor_volume_auc = auc_raw,
      tumor_volume_auc_delta = auc_delta,
      auc_start_day = auc_start_day,
      auc_end_day = auc_end_day,
      auc_n_days = sum(finite_auc),
      stringsAsFactors = FALSE
    )
    for (day_column in endpoint_columns) {
      row[[paste0("tumor_volume_", day_column)]] <- values[[day_column]]
      row[[paste0("tumor_volume_delta_", day_column)]] <- values[[day_column]] - baseline_value
    }
    row
  })
  output <- do.call(rbind, rows)
  keep <- !is.na(output$sample_id) & nzchar(output$sample_id) & is.finite(output$tumor_volume_delta)
  output <- output[keep, , drop = FALSE]
  output$dose <- essential_standardize_dose(output$dose)
  output$control_match_group <- output$initial_ploidy
  output$matched_control_mean_delta <- NA_real_
  output$matched_control_n <- NA_integer_
  output$matched_control_mean_auc_delta <- NA_real_
  output$matched_control_n_auc <- NA_integer_
  for (day_column in endpoint_columns) {
    output[[paste0("matched_control_mean_delta_", day_column)]] <- NA_real_
    output[[paste0("matched_control_n_", day_column)]] <- NA_integer_
  }

  control_dose <- essential_standardize_dose(control_dose)[[1L]]
  groups <- unique(output$control_match_group[!is.na(output$control_match_group)])
  for (group in groups) {
    target <- output$control_match_group == group
    control <- target & output$dose == control_dose & is.finite(output$tumor_volume_delta)
    output$matched_control_mean_delta[target] <- mean(output$tumor_volume_delta[control], na.rm = TRUE)
    output$matched_control_n[target] <- sum(control)
    control_auc <- target & output$dose == control_dose & is.finite(output$tumor_volume_auc_delta)
    output$matched_control_mean_auc_delta[target] <- mean(output$tumor_volume_auc_delta[control_auc], na.rm = TRUE)
    output$matched_control_n_auc[target] <- sum(control_auc)
    for (day_column in endpoint_columns) {
      delta_column <- paste0("tumor_volume_delta_", day_column)
      mean_column <- paste0("matched_control_mean_delta_", day_column)
      count_column <- paste0("matched_control_n_", day_column)
      control_day <- target & output$dose == control_dose & is.finite(output[[delta_column]])
      output[[mean_column]][target] <- mean(output[[delta_column]][control_day], na.rm = TRUE)
      output[[count_column]][target] <- sum(control_day)
    }
  }

  output$TGI_percent <- 100 * (1 - output$tumor_volume_delta / output$matched_control_mean_delta)
  output$TGI_percent[!is.finite(output$TGI_percent)] <- NA_real_
  output$TGI_percent_auc <- 100 * (1 - output$tumor_volume_auc_delta / output$matched_control_mean_auc_delta)
  output$TGI_percent_auc[!is.finite(output$TGI_percent_auc)] <- NA_real_
  for (day_column in endpoint_columns) {
    delta_column <- paste0("tumor_volume_delta_", day_column)
    mean_column <- paste0("matched_control_mean_delta_", day_column)
    tgi_column <- paste0("TGI_percent_", day_column)
    output[[tgi_column]] <- 100 * (1 - output[[delta_column]] / output[[mean_column]])
    output[[tgi_column]][!is.finite(output[[tgi_column]])] <- NA_real_
  }
  output[order(output$initial_ploidy, output$dose, output$sample_id), , drop = FALSE]
}

essential_build_sample_summary <- function(target, tumor) {
  sample_ids <- unique(target$sampleID)
  tumor_counts <- table(tumor$sampleID)
  rows <- lapply(sample_ids, function(sample_id) {
    local <- target[target$sampleID == sample_id, , drop = FALSE]
    mean_ploidy <- mean(local$cell_ploidy, na.rm = TRUE)
    median_ploidy <- stats::median(local$cell_ploidy, na.rm = TRUE)
    mean_pseudotime <- mean(local$pseudotime, na.rm = TRUE)
    median_pseudotime <- stats::median(local$pseudotime, na.rm = TRUE)
    n_tumor <- unname(tumor_counts[[sample_id]])
    data.frame(
      sampleID = sample_id,
      Ploidy = local$Ploidy[[1L]],
      Dose = local$Dose[[1L]],
      DosePanel = local$DosePanel[[1L]],
      n_target_cells = nrow(local),
      mean_cell_ploidy = if (is.nan(mean_ploidy)) NA_real_ else mean_ploidy,
      median_cell_ploidy = if (is.nan(median_ploidy)) NA_real_ else median_ploidy,
      mean_pseudotime = if (is.nan(mean_pseudotime)) NA_real_ else mean_pseudotime,
      median_pseudotime = if (is.nan(median_pseudotime)) NA_real_ else median_pseudotime,
      n_tumor_cells = n_tumor,
      target_cluster_fraction = nrow(local) / n_tumor,
      target_cluster_percent = 100 * nrow(local) / n_tumor,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

essential_prepare_tgi_join <- function(tgi) {
  base <- data.frame(
    sampleID = tgi$sample_id,
    growth_curve_sample_id_raw = tgi$sample_id_raw,
    growth_curve_harvest = tgi$harvest,
    tgi_initial_ploidy = tgi$initial_ploidy,
    tgi_dose = tgi$dose,
    tumor_volume_baseline_day = tgi$baseline_day,
    tumor_volume_baseline = tgi$baseline_volume,
    tumor_volume_final_day = tgi$final_day,
    tumor_volume_final = tgi$final_volume,
    tumor_volume_delta = tgi$tumor_volume_delta,
    matched_control_mean_delta = tgi$matched_control_mean_delta,
    matched_control_n = tgi$matched_control_n,
    TGI_percent = tgi$TGI_percent,
    stringsAsFactors = FALSE
  )
  extra_columns <- grep(
    "^(tumor_volume_Day_|tumor_volume_delta_Day_|matched_control_mean_delta_Day_|matched_control_n_Day_|TGI_percent_Day_|tumor_volume_auc|tumor_volume_auc_delta|matched_control_mean_auc_delta|matched_control_n_auc|TGI_percent_auc|auc_)",
    names(tgi),
    value = TRUE
  )
  cbind(base, tgi[, extra_columns, drop = FALSE])
}

essential_build_compartment_table <- function(tumor, target_clusters, tgi) {
  target <- tumor[tumor$clusters %in% target_clusters, , drop = FALSE]
  target <- target[
    !is.na(target$sampleID) & nzchar(target$sampleID) & is.finite(target$pseudotime),
    ,
    drop = FALSE
  ]
  if (nrow(target) == 0L) {
    stop("No target Tumor cells were found for clusters: ", paste(target_clusters, collapse = ", "), call. = FALSE)
  }
  sample_summary <- essential_build_sample_summary(target, tumor)
  summary_index <- match(target$sampleID, sample_summary$sampleID)
  tgi_join <- essential_prepare_tgi_join(tgi)
  tgi_index <- match(target$sampleID, tgi_join$sampleID)
  if (anyNA(tgi_index)) {
    missing_samples <- unique(target$sampleID[is.na(tgi_index)])
    warning("TGI was not matched for samples: ", paste(missing_samples, collapse = ", "), call. = FALSE)
  }

  output <- data.frame(
    cell_id = target$cell,
    sample_id = target$sampleID,
    cluster = target$clusters,
    initial_ploidy = target$Ploidy,
    gemcitabine_dose = target$Dose,
    gemcitabine_dose_mg_per_kg = essential_safe_numeric(target$DosePanel),
    pseudotime = target$pseudotime,
    cell_ploidy = target$cell_ploidy,
    average_ploidy = sample_summary$mean_cell_ploidy[summary_index],
    median_ploidy = sample_summary$median_cell_ploidy[summary_index],
    sample_mean_pseudotime = sample_summary$mean_pseudotime[summary_index],
    sample_median_pseudotime = sample_summary$median_pseudotime[summary_index],
    n_target_cells = sample_summary$n_target_cells[summary_index],
    n_tumor_cells = sample_summary$n_tumor_cells[summary_index],
    target_cluster_fraction = sample_summary$target_cluster_fraction[summary_index],
    target_cluster_percent = sample_summary$target_cluster_percent[summary_index],
    stringsAsFactors = FALSE
  )
  output <- cbind(output, tgi_join[tgi_index, setdiff(names(tgi_join), "sampleID"), drop = FALSE])
  order_index <- order(
    output$initial_ploidy,
    output$gemcitabine_dose_mg_per_kg,
    output$sample_id,
    output$pseudotime,
    output$cell_id,
    na.last = TRUE
  )
  output[order_index, , drop = FALSE]
}

essential_generate_analysis_inputs <- function(
  scvelo_path,
  cell_ploidy_path,
  sample_info_path,
  growth_curve_path,
  cellcycle_output = NULL,
  noncellcycle_output = NULL,
  cellcycle_clusters = c("4c", "6", "10")
) {
  essential_require_input_generation_packages()
  metrics <- essential_build_cell_metrics(scvelo_path, cell_ploidy_path, sample_info_path)
  keep <- !is.na(metrics$cell) & nzchar(metrics$cell) & is.finite(metrics$pseudotime)
  metrics <- metrics[keep, , drop = FALSE]
  metrics$DosePanel <- essential_dose_panel(metrics$Dose)
  metrics$sampleID <- metrics$sample
  tumor <- metrics[
    !is.na(metrics$TN) & metrics$TN == "Tumor" &
      metrics$Ploidy %in% c("2N", "4N") &
      metrics$DosePanel %in% c("0", "30", "120"),
    ,
    drop = FALSE
  ]
  if (nrow(tumor) == 0L) stop("No eligible Tumor cells were found in scVelo metrics", call. = FALSE)
  tgi <- essential_calculate_tgi(growth_curve_path)
  requested <- c(
    if (!is.null(cellcycle_output)) "CellCycle" else character(),
    if (!is.null(noncellcycle_output)) "NonCellCycle" else character()
  )
  result <- list()
  if (!is.null(cellcycle_output)) {
    result$CellCycle <- essential_build_compartment_table(tumor, cellcycle_clusters, tgi)
    essential_ensure_dir(dirname(cellcycle_output))
    readr::write_csv(result$CellCycle, cellcycle_output, na = "NA")
  }
  if (!is.null(noncellcycle_output)) {
    noncellcycle_clusters <- setdiff(unique(tumor$clusters), cellcycle_clusters)
    result$NonCellCycle <- essential_build_compartment_table(tumor, noncellcycle_clusters, tgi)
    essential_ensure_dir(dirname(noncellcycle_output))
    readr::write_csv(result$NonCellCycle, noncellcycle_output, na = "NA")
  }
  if (length(requested) == 0L) stop("No generated compartment output was requested", call. = FALSE)
  result
}

essential_save_pdf <- function(plot, path, width, height) {
  essential_ensure_dir(dirname(path))
  suppressMessages(ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    device = grDevices::cairo_pdf
  ))
  invisible(path)
}

essential_format_number <- function(x, digits = 3L) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "fg"), "NA")
}

essential_format_p <- function(x) {
  ifelse(
    is.finite(x),
    ifelse(x < 0.001, formatC(x, digits = 2, format = "e"), formatC(x, digits = 3, format = "fg")),
    "NA"
  )
}

essential_endpoint_tgi_columns <- function(names_vector) {
  columns <- grep("^TGI_percent_Day_[0-9]+$", names_vector, value = TRUE)
  days <- essential_safe_numeric(sub("^TGI_percent_Day_", "", columns))
  columns[order(days)]
}

essential_endpoint_day <- function(column) {
  essential_safe_numeric(sub("^TGI_percent_Day_([0-9]+).*$", "\\1", column))
}

essential_tgi_spec <- function(outcome = "day", control_summary = "mean", tgi_day = 17L) {
  outcome <- tolower(trimws(as.character(outcome)))
  control_summary <- tolower(trimws(as.character(control_summary)))
  if (!outcome %in% c("auc", "day")) {
    stop("tgi-outcome must be one of: auc, day", call. = FALSE)
  }
  if (!control_summary %in% c("mean", "median", "max")) {
    stop("control-summary must be one of: mean, median, max", call. = FALSE)
  }
  day_value <- if (is.null(tgi_day) || !nzchar(as.character(tgi_day))) {
    NA_real_
  } else {
    essential_safe_numeric(tgi_day)
  }
  day <- if (
    is.finite(day_value) && day_value >= 0 && day_value == floor(day_value)
  ) as.integer(day_value) else NA_integer_
  if (identical(outcome, "day") && is.na(day)) {
    stop("A non-negative integer --tgi-day is required when --tgi-outcome=day", call. = FALSE)
  }
  if (identical(outcome, "auc")) day <- NA_integer_

  base_measure <- if (identical(outcome, "auc")) {
    "TGI_percent_auc"
  } else {
    paste0("TGI_percent_Day_", day)
  }
  measure <- if (identical(control_summary, "mean")) {
    base_measure
  } else {
    paste0(base_measure, "_control_", control_summary)
  }
  delta_column <- if (identical(outcome, "auc")) {
    "tumor_volume_auc_delta"
  } else {
    paste0("tumor_volume_delta_Day_", day)
  }
  short_label <- if (identical(outcome, "auc")) "AUC TGI" else paste0("Day ", day, " TGI")
  reference_label <- paste0(
    control_summary,
    " of initial-ploidy-matched untreated controls"
  )
  axis_label <- if (identical(control_summary, "mean")) {
    paste0(short_label, " (%)")
  } else {
    paste0(short_label, " (%; ", reference_label, ")")
  }
  title_label <- if (identical(control_summary, "mean")) {
    short_label
  } else {
    paste0(short_label, " [", reference_label, "]")
  }
  list(
    outcome = outcome,
    control_summary = control_summary,
    day = day,
    base_measure = base_measure,
    measure = measure,
    delta_column = delta_column,
    short_label = short_label,
    title_label = title_label,
    axis_label = axis_label,
    reference_label = reference_label
  )
}

essential_tgi_spec_from_measure <- function(measure) {
  measure <- as.character(measure)[[1L]]
  auc_match <- regexec("^TGI_percent_auc(?:_control_(mean|median|max))?$", measure, perl = TRUE)
  auc_parts <- regmatches(measure, auc_match)[[1L]]
  if (length(auc_parts) > 0L) {
    summary <- if (length(auc_parts) >= 2L && nzchar(auc_parts[[2L]])) auc_parts[[2L]] else "mean"
    return(essential_tgi_spec("auc", summary, NULL))
  }
  day_match <- regexec(
    "^TGI_percent_Day_([0-9]+)(?:_control_(mean|median|max))?$",
    measure,
    perl = TRUE
  )
  day_parts <- regmatches(measure, day_match)[[1L]]
  if (length(day_parts) > 0L) {
    summary <- if (length(day_parts) >= 3L && nzchar(day_parts[[3L]])) day_parts[[3L]] else "mean"
    return(essential_tgi_spec("day", summary, day_parts[[2L]]))
  }
  stop("Unsupported TGI measure in saved tables: ", measure, call. = FALSE)
}

essential_tgi_columns <- function(names_vector) {
  columns <- grep(
    "^TGI_percent_(auc(?:_control_(?:mean|median|max))?|Day_[0-9]+(?:_control_(?:mean|median|max))?)$",
    names_vector,
    value = TRUE,
    perl = TRUE
  )
  auc <- columns[grepl("^TGI_percent_auc", columns)]
  endpoint <- columns[grepl("^TGI_percent_Day_", columns)]
  endpoint_days <- essential_endpoint_day(endpoint)
  unique(c(sort(auc), endpoint[order(endpoint_days, endpoint)]))
}

essential_apply_tgi_spec <- function(cellcycle, noncellcycle, tgi_spec) {
  if (identical(tgi_spec$control_summary, "mean")) {
    missing <- Filter(
      function(data) !tgi_spec$measure %in% names(data),
      list(cellcycle, noncellcycle)
    )
    if (length(missing) > 0L) {
      stop("Input is missing selected TGI column: ", tgi_spec$measure, call. = FALSE)
    }
    return(list(cellcycle = cellcycle, noncellcycle = noncellcycle))
  }

  required <- c(
    "sample_id", "initial_ploidy", "gemcitabine_dose_mg_per_kg", tgi_spec$delta_column
  )
  missing <- setdiff(required, intersect(names(cellcycle), names(noncellcycle)))
  if (length(missing) > 0L) {
    stop(
      "Input is missing columns required for ", tgi_spec$measure, ": ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  union <- rbind(
    cellcycle[, required, drop = FALSE],
    noncellcycle[, required, drop = FALSE]
  )
  sample_ids <- sort(unique(as.character(union$sample_id)))
  sample_rows <- lapply(sample_ids, function(sample_id) {
    data.frame(
      sample_id = sample_id,
      initial_ploidy = essential_unique_sample_value(union, "initial_ploidy", sample_id),
      dose_mg = essential_unique_sample_value(
        union, "gemcitabine_dose_mg_per_kg", sample_id, numeric = TRUE
      ),
      delta = essential_unique_sample_value(union, tgi_spec$delta_column, sample_id, numeric = TRUE),
      stringsAsFactors = FALSE
    )
  })
  samples <- do.call(rbind, sample_rows)
  samples$reference_delta <- NA_real_
  for (ploidy in unique(samples$initial_ploidy)) {
    control_values <- samples$delta[
      samples$initial_ploidy == ploidy & samples$dose_mg == 0 & is.finite(samples$delta)
    ]
    if (length(control_values) == 0L) {
      stop("No finite 0 mg/kg AUC/day controls for initial ploidy ", ploidy, call. = FALSE)
    }
    reference <- switch(
      tgi_spec$control_summary,
      median = stats::median(control_values),
      max = max(control_values)
    )
    samples$reference_delta[samples$initial_ploidy == ploidy] <- reference
  }
  samples$selected_tgi <- 100 * (1 - samples$delta / samples$reference_delta)
  samples$selected_tgi[!is.finite(samples$selected_tgi)] <- NA_real_
  selected <- setNames(samples$selected_tgi, samples$sample_id)
  cellcycle[[tgi_spec$measure]] <- unname(selected[as.character(cellcycle$sample_id)])
  noncellcycle[[tgi_spec$measure]] <- unname(selected[as.character(noncellcycle$sample_id)])
  list(cellcycle = cellcycle, noncellcycle = noncellcycle)
}

essential_add_tgi_plot_metadata <- function(data, tgi_spec, value_column = NULL) {
  if (!is.null(value_column)) {
    if (!value_column %in% names(data)) stop("Missing TGI plot value column: ", value_column, call. = FALSE)
    data$tgi_value <- essential_safe_numeric(data[[value_column]])
  } else if (!"tgi_value" %in% names(data)) {
    if (!tgi_spec$measure %in% names(data)) stop("Missing selected TGI measure: ", tgi_spec$measure, call. = FALSE)
    data$tgi_value <- essential_safe_numeric(data[[tgi_spec$measure]])
  }
  data$tgi_measure <- tgi_spec$measure
  data$tgi_outcome <- tgi_spec$outcome
  data$tgi_control_summary <- tgi_spec$control_summary
  data$tgi_day <- tgi_spec$day
  data$tgi_axis_label <- tgi_spec$axis_label
  data$tgi_title_label <- tgi_spec$title_label
  data
}

essential_method_specs <- function() {
  list(
    initial_ploidy = list(
      method = "initial_ploidy",
      type = "initial_ploidy",
      threshold = NA_real_,
      group_levels = c("2N", "4N"),
      group_label = "Initial ploidy",
      group_column = "initial_ploidy",
      group_numeric_column = "initial_ploidy_numeric"
    ),
    ETP_fixed_threshold_2_25 = list(
      method = "ETP_fixed_threshold_2_25",
      type = "ETP",
      threshold = 2.25,
      group_levels = c("ETP-lower", "ETP-higher"),
      group_label = "EndTimePoint ploidy",
      group_column = "end_timepoint_ploidy_group",
      group_numeric_column = "end_timepoint_ploidy_group_numeric"
    ),
    ETP_boundary_stress_threshold_2_375 = list(
      method = "ETP_boundary_stress_threshold_2_375",
      type = "ETP",
      threshold = 2.375,
      group_levels = c("ETP-lower", "ETP-higher"),
      group_label = "EndTimePoint ploidy",
      group_column = "end_timepoint_ploidy_group",
      group_numeric_column = "end_timepoint_ploidy_group_numeric"
    ),
    ETP_reference_balanced_threshold_2_24 = list(
      method = "ETP_reference_balanced_threshold_2_24",
      type = "ETP",
      threshold = 2.24,
      group_levels = c("ETP-lower", "ETP-higher"),
      group_label = "EndTimePoint ploidy",
      group_column = "end_timepoint_ploidy_group",
      group_numeric_column = "end_timepoint_ploidy_group_numeric"
    )
  )
}

essential_required_input_columns <- function() {
  c(
    "cell_id", "sample_id", "initial_ploidy", "gemcitabine_dose",
    "gemcitabine_dose_mg_per_kg", "pseudotime", "cell_ploidy", "TGI_percent_auc"
  )
}

essential_read_input <- function(path, compartment) {
  if (!file.exists(path)) stop("Missing input CSV: ", path, call. = FALSE)
  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing <- setdiff(essential_required_input_columns(), names(data))
  if (length(missing) > 0L) {
    stop("Input is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  numeric_columns <- unique(c(
    "gemcitabine_dose_mg_per_kg", "pseudotime", "cell_ploidy", "TGI_percent_auc",
    grep("^tumor_volume_(auc_delta|delta_Day_[0-9]+)$", names(data), value = TRUE),
    essential_endpoint_tgi_columns(names(data))
  ))
  for (column in intersect(numeric_columns, names(data))) {
    data[[column]] <- essential_safe_numeric(data[[column]])
  }
  data$compartment <- compartment
  data
}

essential_unique_sample_value <- function(data, column, sample_id, numeric = FALSE) {
  values <- data[data$sample_id == sample_id, column]
  if (numeric) {
    values <- essential_safe_numeric(values)
    values <- unique(values[is.finite(values)])
  } else {
    values <- unique(as.character(values[!is.na(values)]))
  }
  if (length(values) > 1L) {
    stop("Inconsistent sample-level column ", column, " for sample ", sample_id, call. = FALSE)
  }
  if (length(values) == 0L) {
    if (numeric) NA_real_ else NA_character_
  } else {
    values[[1L]]
  }
}

essential_derive_assignments <- function(cellcycle, noncellcycle) {
  needed <- c(
    "cell_id", "sample_id", "initial_ploidy", "gemcitabine_dose",
    "gemcitabine_dose_mg_per_kg", "cell_ploidy"
  )
  union_cells <- rbind(
    cellcycle[, needed, drop = FALSE],
    noncellcycle[, needed, drop = FALSE]
  )
  union_cells$cell_ploidy <- essential_safe_numeric(union_cells$cell_ploidy)
  union_cells$gemcitabine_dose_mg_per_kg <- essential_safe_numeric(union_cells$gemcitabine_dose_mg_per_kg)
  union_cells <- union_cells[is.finite(union_cells$cell_ploidy), , drop = FALSE]
  if (nrow(union_cells) == 0L) stop("No finite cell ploidy values were found", call. = FALSE)

  duplicated_ids <- unique(union_cells$cell_id[duplicated(union_cells$cell_id)])
  for (cell_id in duplicated_ids) {
    local <- union_cells[union_cells$cell_id == cell_id, , drop = FALSE]
    if (length(unique(local$sample_id)) != 1L || length(unique(local$cell_ploidy)) != 1L) {
      stop("Conflicting duplicated cell_id: ", cell_id, call. = FALSE)
    }
  }
  union_cells <- union_cells[!duplicated(union_cells$cell_id), , drop = FALSE]

  sample_ids <- sort(unique(union_cells$sample_id))
  rows <- lapply(sample_ids, function(sample_id) {
    local <- union_cells[union_cells$sample_id == sample_id, , drop = FALSE]
    data.frame(
      sample_id = sample_id,
      initial_ploidy = essential_unique_sample_value(local, "initial_ploidy", sample_id),
      dose = essential_unique_sample_value(local, "gemcitabine_dose", sample_id),
      dose_mg = essential_unique_sample_value(local, "gemcitabine_dose_mg_per_kg", sample_id, numeric = TRUE),
      n_endpoint_ploidy_cells = nrow(local),
      sample_mean_endpoint_ploidy = mean(local$cell_ploidy),
      sample_median_endpoint_ploidy = median(local$cell_ploidy),
      sample_max_endpoint_ploidy = max(local$cell_ploidy),
      stringsAsFactors = FALSE
    )
  })
  assignments <- do.call(rbind, rows)
  unexpected <- setdiff(unique(assignments$initial_ploidy), c("2N", "4N"))
  if (length(unexpected) > 0L) {
    stop("Unexpected initial ploidy labels: ", paste(unexpected, collapse = ", "), call. = FALSE)
  }
  assignments
}

essential_prepare_cellcycle <- function(cellcycle, assignments, spec) {
  data <- cellcycle[is.finite(cellcycle$pseudotime), , drop = FALSE]
  assignment_index <- match(data$sample_id, assignments$sample_id)
  if (anyNA(assignment_index)) {
    stop("Missing assignment for CellCycle samples", call. = FALSE)
  }
  joined <- assignments[assignment_index, , drop = FALSE]
  if (any(as.character(data$initial_ploidy) != as.character(joined$initial_ploidy))) {
    stop("Initial ploidy metadata disagree between inputs", call. = FALSE)
  }
  data$dose <- as.character(data$gemcitabine_dose)
  data$dose_mg <- essential_safe_numeric(data$gemcitabine_dose_mg_per_kg)
  data$sample_mean_endpoint_ploidy <- joined$sample_mean_endpoint_ploidy
  data$sample_median_endpoint_ploidy <- joined$sample_median_endpoint_ploidy
  data$sample_max_endpoint_ploidy <- joined$sample_max_endpoint_ploidy
  data$initial_ploidy <- factor(as.character(data$initial_ploidy), levels = c("2N", "4N"))
  data$initial_ploidy_numeric <- ifelse(as.character(data$initial_ploidy) == "4N", 1, 0)
  if (identical(spec$type, "ETP")) {
    etp_group <- ifelse(
      data$sample_mean_endpoint_ploidy > spec$threshold,
      "ETP-higher",
      "ETP-lower"
    )
    data$end_timepoint_ploidy_group <- factor(etp_group, levels = spec$group_levels)
    data$end_timepoint_ploidy_group_numeric <- ifelse(etp_group == "ETP-higher", 1, 0)
    data$analysis_group <- data$end_timepoint_ploidy_group
    data$analysis_group_numeric <- data$end_timepoint_ploidy_group_numeric
  } else {
    data$analysis_group <- factor(as.character(data$initial_ploidy), levels = spec$group_levels)
    data$analysis_group_numeric <- data$initial_ploidy_numeric
  }
  data
}

essential_build_sample_meta <- function(data, tgi_column = "TGI_percent_auc") {
  sample_ids <- sort(unique(data$sample_id))
  tgi_columns <- unique(c(essential_tgi_columns(names(data)), tgi_column))
  rows <- lapply(sample_ids, function(sample_id) {
    local <- data[data$sample_id == sample_id, , drop = FALSE]
    base <- data.frame(
      sample_id = sample_id,
      initial_ploidy = as.character(local$initial_ploidy[[1L]]),
      initial_ploidy_numeric = local$initial_ploidy_numeric[[1L]],
      analysis_group = as.character(local$analysis_group[[1L]]),
      analysis_group_numeric = local$analysis_group_numeric[[1L]],
      dose = as.character(local$dose[[1L]]),
      dose_mg = local$dose_mg[[1L]],
      sample_mean_endpoint_ploidy = local$sample_mean_endpoint_ploidy[[1L]],
      n_cells = nrow(local),
      mean_pseudotime = mean(local$pseudotime),
      median_pseudotime = median(local$pseudotime),
      mean_cell_ploidy = mean(local$cell_ploidy, na.rm = TRUE),
      median_cell_ploidy = median(local$cell_ploidy, na.rm = TRUE),
      p90_cell_ploidy = as.numeric(stats::quantile(local$cell_ploidy, 0.90, na.rm = TRUE, names = FALSE)),
      stringsAsFactors = FALSE
    )
    for (column in tgi_columns) {
      base[[column]] <- essential_unique_sample_value(local, column, sample_id, numeric = TRUE)
    }
    base
  })
  meta <- do.call(rbind, rows)
  meta <- meta[order(meta$dose_mg, meta$analysis_group_numeric, meta$sample_id), , drop = FALSE]
  rownames(meta) <- meta$sample_id
  meta
}

essential_make_grid <- function(data, n = 501L) {
  range_values <- range(data$pseudotime, na.rm = TRUE)
  if (!all(is.finite(range_values)) || diff(range_values) == 0) {
    return(sort(unique(data$pseudotime)))
  }
  seq(range_values[[1L]], range_values[[2L]], length.out = n)
}

essential_sample_ecdf_matrix <- function(data, sample_ids, grid) {
  matrix_values <- t(vapply(sample_ids, function(sample_id) {
    values <- data$pseudotime[data$sample_id == sample_id]
    if (length(values) == 0L) rep(NA_real_, length(grid)) else stats::ecdf(values)(grid)
  }, numeric(length(grid))))
  rownames(matrix_values) <- sample_ids
  matrix_values
}

essential_calculate_shift_metrics <- function(
  data,
  sample_meta,
  pooled_reference = FALSE,
  control_scope = c("analysis_group", "all_untreated")
) {
  control_scope <- match.arg(control_scope)
  grid <- essential_make_grid(data)
  sample_ids <- rownames(sample_meta)
  ecdf_matrix <- essential_sample_ecdf_matrix(data, sample_ids, grid)
  rows <- lapply(sample_ids, function(sample_id) {
    group <- sample_meta[sample_id, "analysis_group"]
    dose_mg <- sample_meta[sample_id, "dose_mg"]
    controls <- if (identical(control_scope, "all_untreated")) {
      sample_meta$sample_id[sample_meta$dose_mg == 0]
    } else {
      sample_meta$sample_id[
        sample_meta$analysis_group == group & sample_meta$dose_mg == 0
      ]
    }
    controls <- sort(as.character(controls))
    reference_samples <- if (dose_mg == 0) setdiff(controls, sample_id) else controls
    if (length(reference_samples) == 0L) {
      delta <- rep(NA_real_, length(grid))
      reference_mean <- NA_real_
      reference_cells <- 0L
    } else if (!pooled_reference) {
      reference_ecdf <- colMeans(ecdf_matrix[reference_samples, , drop = FALSE], na.rm = TRUE)
      delta <- ecdf_matrix[sample_id, ] - reference_ecdf
      reference_mean <- mean(sample_meta[reference_samples, "mean_pseudotime"])
      reference_cells <- sum(sample_meta[reference_samples, "n_cells"])
    } else {
      sample_values <- data$pseudotime[data$sample_id == sample_id]
      reference_values <- data$pseudotime[data$sample_id %in% reference_samples]
      delta <- stats::ecdf(sample_values)(grid) - stats::ecdf(reference_values)(grid)
      reference_mean <- mean(reference_values)
      reference_cells <- length(reference_values)
    }
    finite_delta <- delta[is.finite(delta)]
    row <- data.frame(
      sample_id = sample_id,
      reference_type = if (pooled_reference) {
        "sensitivity_pooled_cell_reference"
      } else if (identical(control_scope, "all_untreated")) {
        "all_untreated_equal_sample_reference"
      } else {
        "primary_equal_sample_reference"
      },
      control_scope = control_scope,
      initial_ploidy = sample_meta[sample_id, "initial_ploidy"],
      initial_ploidy_numeric = sample_meta[sample_id, "initial_ploidy_numeric"],
      analysis_group = group,
      analysis_group_numeric = sample_meta[sample_id, "analysis_group_numeric"],
      dose = sample_meta[sample_id, "dose"],
      dose_mg = dose_mg,
      cell_count = sample_meta[sample_id, "n_cells"],
      sample_mean_endpoint_ploidy = sample_meta[sample_id, "sample_mean_endpoint_ploidy"],
      mean_cell_ploidy = sample_meta[sample_id, "mean_cell_ploidy"],
      median_cell_ploidy = sample_meta[sample_id, "median_cell_ploidy"],
      p90_cell_ploidy = sample_meta[sample_id, "p90_cell_ploidy"],
      mean_pseudotime = sample_meta[sample_id, "mean_pseudotime"],
      median_pseudotime = sample_meta[sample_id, "median_pseudotime"],
      reference_n_samples = length(reference_samples),
      reference_cell_count = reference_cells,
      ecdf_rmse = if (length(finite_delta) > 0L) sqrt(mean(finite_delta^2)) else NA_real_,
      ecdf_ks = if (length(finite_delta) > 0L) max(abs(finite_delta)) else NA_real_,
      ecdf_mean_abs = if (length(finite_delta) > 0L) mean(abs(finite_delta)) else NA_real_,
      signed_mean_shift = sample_meta[sample_id, "mean_pseudotime"] - reference_mean,
      stringsAsFactors = FALSE
    )
    for (column in essential_tgi_columns(names(sample_meta))) {
      row[[column]] <- sample_meta[sample_id, column]
    }
    row
  })
  output <- do.call(rbind, rows)
  output <- output[order(output$dose_mg, output$analysis_group_numeric, output$sample_id), , drop = FALSE]
  rownames(output) <- NULL
  output
}

essential_safe_cor <- function(x, y, method = "pearson") {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(data.frame(n = length(x), estimate = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(stats::cor.test(x, y, method = method, exact = FALSE))
  data.frame(n = length(x), estimate = unname(test$estimate), p_value = test$p.value)
}

essential_permutation_indices <- local({
  cache <- new.env(parent = emptyenv())
  function(n) {
    key <- as.character(n)
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    total <- as.integer(round(factorial(n)))
    matrix_values <- matrix(NA_integer_, nrow = total, ncol = n)
    row_index <- 0L
    recurse <- function(prefix, remaining) {
      if (length(remaining) == 0L) {
        row_index <<- row_index + 1L
        matrix_values[row_index, ] <<- prefix
      } else {
        for (i in seq_along(remaining)) {
          recurse(c(prefix, remaining[[i]]), remaining[-i])
        }
      }
    }
    recurse(integer(0), seq_len(n))
    assign(key, matrix_values, envir = cache)
    matrix_values
  }
})

essential_stratified_permutation_indices <- function(strata) {
  split_indices <- split(seq_along(strata), as.character(strata))
  local_permutations <- lapply(split_indices, function(indices) {
    local <- essential_permutation_indices(length(indices))
    t(apply(local, 1L, function(order_index) indices[order_index]))
  })
  output <- list()
  recurse <- function(i, current) {
    if (i > length(local_permutations)) {
      output[[length(output) + 1L]] <<- current
      return(invisible(NULL))
    }
    for (row in seq_len(nrow(local_permutations[[i]]))) {
      next_index <- current
      indices <- split_indices[[i]]
      next_index[indices] <- local_permutations[[i]][row, ]
      recurse(i + 1L, next_index)
    }
    invisible(NULL)
  }
  recurse(1L, seq_along(strata))
  do.call(rbind, output)
}

essential_permutation_cor <- function(
  x,
  y,
  method = "pearson",
  n_perm = 10000L,
  exact = FALSE,
  strata = NULL
) {
  keep <- is.finite(x) & is.finite(y)
  if (!is.null(strata)) keep <- keep & !is.na(strata)
  x <- x[keep]
  y <- y[keep]
  if (!is.null(strata)) strata <- as.character(strata[keep])
  n <- length(x)
  if (n < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(data.frame(
      n = n, estimate = NA_real_, asymptotic_p = NA_real_,
      permutation_p_two_sided = NA_real_, permutation_p_positive = NA_real_,
      n_permutations = NA_integer_, permutation_mode = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  observed <- suppressWarnings(stats::cor(x, y, method = method))
  asymptotic <- essential_safe_cor(x, y, method)
  exact_count <- if (is.null(strata)) factorial(n) else prod(vapply(table(strata), factorial, numeric(1L)))
  if (isTRUE(exact) && is.finite(exact_count) && exact_count <= 400000) {
    index_matrix <- if (is.null(strata)) {
      essential_permutation_indices(n)
    } else {
      essential_stratified_permutation_indices(strata)
    }
    permuted <- vapply(seq_len(nrow(index_matrix)), function(i) {
      suppressWarnings(stats::cor(x, y[index_matrix[i, ]], method = method))
    }, numeric(1L))
    mode <- if (is.null(strata)) "exact_TGI_label_enumeration" else "exact_within_stratum_TGI_label_enumeration"
    p_two_sided <- mean(abs(permuted) >= abs(observed) - 1e-15, na.rm = TRUE)
    p_positive <- mean(permuted >= observed - 1e-15, na.rm = TRUE)
  } else {
    permuted <- replicate(n_perm, {
      permuted_y <- y
      if (is.null(strata)) {
        permuted_y <- sample(y)
      } else {
        for (level in unique(strata)) {
          indices <- which(strata == level)
          permuted_y[indices] <- sample(y[indices])
        }
      }
      suppressWarnings(stats::cor(x, permuted_y, method = method))
    })
    mode <- if (is.null(strata)) "monte_carlo_TGI_label_permutation" else "monte_carlo_within_stratum_TGI_label_permutation"
    exceed_two <- sum(abs(permuted) >= abs(observed) - 1e-15, na.rm = TRUE)
    exceed_positive <- sum(permuted >= observed - 1e-15, na.rm = TRUE)
    p_two_sided <- (exceed_two + 1) / (sum(is.finite(permuted)) + 1)
    p_positive <- (exceed_positive + 1) / (sum(is.finite(permuted)) + 1)
  }
  data.frame(
    n = n,
    estimate = observed,
    asymptotic_p = asymptotic$p_value,
    permutation_p_two_sided = p_two_sided,
    permutation_p_positive = p_positive,
    n_permutations = length(permuted),
    permutation_mode = mode,
    stringsAsFactors = FALSE
  )
}

essential_choice_lists <- function(groups, group_a, strata) {
  split_indices <- split(seq_along(groups), strata)
  lapply(split_indices, function(indices) {
    k <- sum(groups[indices] == group_a)
    if (k == 0L) return(list(integer(0)))
    if (k == length(indices)) return(list(indices))
    utils::combn(indices, k, simplify = FALSE)
  })
}

essential_exact_group_vectors <- function(groups, group_a, group_b, strata) {
  choices <- essential_choice_lists(groups, group_a, strata)
  output <- list()
  recurse <- function(i, selected) {
    if (i > length(choices)) {
      vector <- rep(group_b, length(groups))
      vector[unlist(selected, use.names = FALSE)] <- group_a
      output[[length(output) + 1L]] <<- vector
    } else {
      for (choice in choices[[i]]) recurse(i + 1L, c(selected, list(choice)))
    }
  }
  recurse(1L, list())
  output
}

essential_sample_group_vector <- function(groups, group_a, group_b, strata) {
  output <- rep(group_b, length(groups))
  for (indices in split(seq_along(groups), strata)) {
    k <- sum(groups[indices] == group_a)
    if (k > 0L) output[sample(indices, k)] <- group_a
  }
  output
}

essential_ecdf_distance <- function(ecdf_matrix, groups, group_a, group_b) {
  if (!any(groups == group_a) || !any(groups == group_b)) {
    return(c(ecdf_rmse = NA_real_, ecdf_ks = NA_real_, ecdf_mean_abs = NA_real_))
  }
  delta <- colMeans(ecdf_matrix[groups == group_a, , drop = FALSE], na.rm = TRUE) -
    colMeans(ecdf_matrix[groups == group_b, , drop = FALSE], na.rm = TRUE)
  c(
    ecdf_rmse = sqrt(mean(delta^2, na.rm = TRUE)),
    ecdf_ks = max(abs(delta), na.rm = TRUE),
    ecdf_mean_abs = mean(abs(delta), na.rm = TRUE)
  )
}

essential_ecdf_group_test <- function(data, group_column, group_a, group_b, stratify, n_perm) {
  local <- data[data[[group_column]] %in% c(group_a, group_b), , drop = FALSE]
  sample_ids <- sort(unique(local$sample_id))
  meta_columns <- unique(c("sample_id", "analysis_group", group_column))
  meta <- local[!duplicated(local$sample_id), meta_columns, drop = FALSE]
  meta <- meta[match(sample_ids, meta$sample_id), , drop = FALSE]
  groups <- as.character(meta[[group_column]])
  not_estimable <- length(sample_ids) < 3L || !all(c(group_a, group_b) %in% groups)
  if (not_estimable) {
    return(data.frame(
      group_a = group_a, group_b = group_b, stratified_by_analysis_group = stratify,
      n_samples = length(sample_ids), n_group_a = sum(groups == group_a), n_group_b = sum(groups == group_b),
      observed_ecdf_rmse = NA_real_, observed_ecdf_ks = NA_real_, observed_ecdf_mean_abs = NA_real_,
      p_ecdf_rmse = NA_real_, p_ecdf_ks = NA_real_, p_ecdf_mean_abs = NA_real_,
      n_permutations = NA_integer_, permutation_mode = "not_estimable", stringsAsFactors = FALSE
    ))
  }
  strata <- if (stratify) as.character(meta$analysis_group) else rep("all", length(groups))
  grid <- essential_make_grid(local)
  ecdf_matrix <- essential_sample_ecdf_matrix(local, sample_ids, grid)
  observed <- essential_ecdf_distance(ecdf_matrix, groups, group_a, group_b)
  choices <- essential_choice_lists(groups, group_a, strata)
  n_exact <- prod(vapply(choices, length, integer(1L)))
  if (is.finite(n_exact) && n_exact <= n_perm) {
    permuted_groups <- essential_exact_group_vectors(groups, group_a, group_b, strata)
    mode <- "exact_label_enumeration"
  } else {
    permuted_groups <- replicate(
      n_perm,
      essential_sample_group_vector(groups, group_a, group_b, strata),
      simplify = FALSE
    )
    mode <- "monte_carlo_label_permutation"
  }
  permuted_stats <- t(vapply(permuted_groups, function(vector) {
    essential_ecdf_distance(ecdf_matrix, vector, group_a, group_b)
  }, numeric(3L)))
  data.frame(
    group_a = group_a,
    group_b = group_b,
    stratified_by_analysis_group = stratify,
    n_samples = length(sample_ids),
    n_group_a = sum(groups == group_a),
    n_group_b = sum(groups == group_b),
    observed_ecdf_rmse = observed[["ecdf_rmse"]],
    observed_ecdf_ks = observed[["ecdf_ks"]],
    observed_ecdf_mean_abs = observed[["ecdf_mean_abs"]],
    p_ecdf_rmse = mean(permuted_stats[, "ecdf_rmse"] >= observed[["ecdf_rmse"]] - 1e-15, na.rm = TRUE),
    p_ecdf_ks = mean(permuted_stats[, "ecdf_ks"] >= observed[["ecdf_ks"]] - 1e-15, na.rm = TRUE),
    p_ecdf_mean_abs = mean(permuted_stats[, "ecdf_mean_abs"] >= observed[["ecdf_mean_abs"]] - 1e-15, na.rm = TRUE),
    n_permutations = nrow(permuted_stats),
    permutation_mode = mode,
    stringsAsFactors = FALSE
  )
}

essential_run_dose_tests <- function(data, n_perm) {
  local <- data
  local$combined_treatment_group <- ifelse(local$dose_mg == 0, "0mg/kg", "treated")
  output <- list(
    cbind(
      comparison = "0_vs_30plus120",
      essential_ecdf_group_test(local, "combined_treatment_group", "0mg/kg", "treated", FALSE, n_perm)
    ),
    cbind(
      comparison = "0_vs_30plus120",
      essential_ecdf_group_test(local, "combined_treatment_group", "0mg/kg", "treated", TRUE, n_perm)
    )
  )
  dose_values <- sort(unique(local$dose_mg))
  for (pair in split(utils::combn(dose_values, 2L), col(utils::combn(dose_values, 2L)))) {
    dose_a <- pair[[1L]]
    dose_b <- pair[[2L]]
    subset <- local[local$dose_mg %in% c(dose_a, dose_b), , drop = FALSE]
    group_a <- unique(subset$dose[subset$dose_mg == dose_a])[[1L]]
    group_b <- unique(subset$dose[subset$dose_mg == dose_b])[[1L]]
    output[[length(output) + 1L]] <- cbind(
      comparison = paste0(dose_a, "_vs_", dose_b),
      essential_ecdf_group_test(subset, "dose", group_a, group_b, FALSE, n_perm)
    )
  }
  do.call(rbind, output)
}

essential_run_tgi_associations <- function(
  equal_shift,
  pooled_shift,
  method,
  n_perm,
  primary_tgi_column = "TGI_percent_auc"
) {
  references <- list(
    primary_equal_sample_reference = equal_shift,
    sensitivity_pooled_cell_reference = pooled_shift
  )
  shift_columns <- c("ecdf_rmse", "ecdf_ks", "ecdf_mean_abs", "signed_mean_shift")
  tgi_columns <- essential_tgi_columns(names(equal_shift))
  rows <- list()
  for (reference_name in names(references)) {
    treated <- references[[reference_name]][references[[reference_name]]$dose_mg > 0, , drop = FALSE]
    for (shift_column in shift_columns) {
      for (tgi_column in tgi_columns) {
        exact <- identical(tgi_column, primary_tgi_column)
        pearson <- essential_permutation_cor(treated[[shift_column]], treated[[tgi_column]], "pearson", n_perm, exact)
        spearman <- essential_permutation_cor(treated[[shift_column]], treated[[tgi_column]], "spearman", n_perm, exact)
        rows[[length(rows) + 1L]] <- data.frame(
          method = method,
          compartment = "CellCycle",
          sample_set = "treated",
          reference_type = reference_name,
          shift_metric = shift_column,
          tgi_measure = tgi_column,
          tgi_day = if (grepl("^TGI_percent_Day_", tgi_column)) essential_endpoint_day(tgi_column) else NA_real_,
          n = pearson$n,
          pearson_r = pearson$estimate,
          pearson_p_asymptotic = pearson$asymptotic_p,
          pearson_p_permutation_two_sided = pearson$permutation_p_two_sided,
          pearson_p_permutation_positive = pearson$permutation_p_positive,
          pearson_n_permutations = pearson$n_permutations,
          pearson_permutation_mode = pearson$permutation_mode,
          spearman_rho = spearman$estimate,
          spearman_p_asymptotic = spearman$asymptotic_p,
          spearman_p_permutation_two_sided = spearman$permutation_p_two_sided,
          spearman_p_permutation_positive = spearman$permutation_p_positive,
          spearman_n_permutations = spearman$n_permutations,
          spearman_permutation_mode = spearman$permutation_mode,
          pre_specified_primary = reference_name == "primary_equal_sample_reference" &&
            shift_column == "ecdf_rmse" && tgi_column == primary_tgi_column,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  output <- do.call(rbind, rows)
  output$pearson_q_exploratory <- stats::p.adjust(output$pearson_p_permutation_two_sided, method = "BH")
  output$spearman_q_exploratory <- stats::p.adjust(output$spearman_p_permutation_two_sided, method = "BH")
  output
}

essential_run_mean_etp_associations <- function(
  equal_shift,
  method,
  n_perm,
  primary_tgi_column = "TGI_percent_auc"
) {
  treated <- equal_shift[equal_shift$dose_mg > 0, , drop = FALSE]
  tgi_columns <- essential_tgi_columns(names(treated))
  rows <- lapply(tgi_columns, function(tgi_column) {
    exact <- identical(tgi_column, primary_tgi_column)
    pearson <- essential_permutation_cor(
      treated$sample_mean_endpoint_ploidy,
      treated[[tgi_column]],
      "pearson",
      n_perm,
      exact
    )
    spearman <- essential_permutation_cor(
      treated$sample_mean_endpoint_ploidy,
      treated[[tgi_column]],
      "spearman",
      n_perm,
      exact
    )
    data.frame(
      method = method,
      compartment = "CellCycle",
      sample_set = "treated",
      predictor = "sample_mean_endpoint_ploidy",
      predictor_label = "sample_mean_ETP",
      tgi_measure = tgi_column,
      tgi_day = if (grepl("^TGI_percent_Day_", tgi_column)) essential_endpoint_day(tgi_column) else NA_real_,
      n = pearson$n,
      pearson_r = pearson$estimate,
      pearson_p_asymptotic = pearson$asymptotic_p,
      pearson_p_permutation_two_sided = pearson$permutation_p_two_sided,
      pearson_p_permutation_positive = pearson$permutation_p_positive,
      pearson_n_permutations = pearson$n_permutations,
      pearson_permutation_mode = pearson$permutation_mode,
      spearman_rho = spearman$estimate,
      spearman_p_asymptotic = spearman$asymptotic_p,
      spearman_p_permutation_two_sided = spearman$permutation_p_two_sided,
      spearman_p_permutation_positive = spearman$permutation_p_positive,
      spearman_n_permutations = spearman$n_permutations,
      spearman_permutation_mode = spearman$permutation_mode,
      multiplicity_family = "CellCycle_mean_ETP_TGI",
      analysis_role = if (exact) "selected_TGI_outcome" else "exploratory_mean_ETP",
      pre_specified_primary = exact,
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, rows)
  output$pearson_q_exploratory <- stats::p.adjust(output$pearson_p_permutation_two_sided, method = "BH")
  output$spearman_q_exploratory <- stats::p.adjust(output$spearman_p_permutation_two_sided, method = "BH")
  output
}

essential_leave_one_out <- function(data, predictor, response = "TGI_percent_auc") {
  full_pearson <- essential_safe_cor(data[[predictor]], data[[response]], "pearson")
  rows <- lapply(data$sample_id, function(sample_id) {
    local <- data[data$sample_id != sample_id, , drop = FALSE]
    pearson <- essential_safe_cor(local[[predictor]], local[[response]], "pearson")
    spearman <- essential_safe_cor(local[[predictor]], local[[response]], "spearman")
    data.frame(
      tgi_measure = response,
      omitted_sample_id = sample_id,
      n = pearson$n,
      pearson_r = pearson$estimate,
      pearson_p = pearson$p_value,
      spearman_rho = spearman$estimate,
      spearman_p = spearman$p_value,
      sign_matches_full_pearson = sign(pearson$estimate) == sign(full_pearson$estimate),
      delta_pearson_r_from_full = pearson$estimate - full_pearson$estimate,
      influential_flag = is.finite(pearson$estimate) && is.finite(full_pearson$estimate) &&
        (sign(pearson$estimate) != sign(full_pearson$estimate) || abs(pearson$estimate - full_pearson$estimate) > 0.20),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

essential_bootstrap <- function(data, predictor, n_boot, response = "TGI_percent_auc") {
  rows <- lapply(seq_len(n_boot), function(iteration) {
    indices <- sample(seq_len(nrow(data)), replace = TRUE)
    local <- data[indices, , drop = FALSE]
    pearson <- essential_safe_cor(local[[predictor]], local[[response]], "pearson")
    spearman <- essential_safe_cor(local[[predictor]], local[[response]], "spearman")
    data.frame(
      tgi_measure = response,
      iteration = iteration,
      n = pearson$n,
      pearson_r = pearson$estimate,
      pearson_p = pearson$p_value,
      spearman_rho = spearman$estimate,
      spearman_p = spearman$p_value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

essential_quantile <- function(x, probability) {
  finite <- x[is.finite(x)]
  if (length(finite) == 0L) return(NA_real_)
  as.numeric(stats::quantile(finite, probability, names = FALSE))
}

essential_primary_robustness_summary <- function(association, leave_one_out, bootstrap) {
  primary <- association[association$pre_specified_primary, , drop = FALSE]
  data.frame(
    tgi_measure = primary$tgi_measure[[1L]],
    full_n = primary$n[[1L]],
    full_pearson_r = primary$pearson_r[[1L]],
    full_pearson_perm_p = primary$pearson_p_permutation_two_sided[[1L]],
    full_spearman_rho = primary$spearman_rho[[1L]],
    full_spearman_perm_p = primary$spearman_p_permutation_two_sided[[1L]],
    leave_one_out_all_positive = all(leave_one_out$pearson_r > 0, na.rm = TRUE),
    leave_one_out_min_pearson_r = min(leave_one_out$pearson_r, na.rm = TRUE),
    leave_one_out_max_pearson_p = max(leave_one_out$pearson_p, na.rm = TRUE),
    n_influential_flags = sum(leave_one_out$influential_flag, na.rm = TRUE),
    bootstrap_pearson_r_ci025 = essential_quantile(bootstrap$pearson_r, 0.025),
    bootstrap_pearson_r_ci975 = essential_quantile(bootstrap$pearson_r, 0.975),
    bootstrap_spearman_rho_ci025 = essential_quantile(bootstrap$spearman_rho, 0.025),
    bootstrap_spearman_rho_ci975 = essential_quantile(bootstrap$spearman_rho, 0.975),
    stringsAsFactors = FALSE
  )
}

essential_mean_etp_robustness_summary <- function(
  method,
  association,
  leave_one_out,
  bootstrap,
  n_boot,
  primary_tgi_column = "TGI_percent_auc"
) {
  primary <- association[association$tgi_measure == primary_tgi_column, , drop = FALSE]
  data.frame(
    method = method,
    predictor = "sample_mean_endpoint_ploidy",
    tgi_measure = primary_tgi_column,
    full_n = primary$n[[1L]],
    full_pearson_r = primary$pearson_r[[1L]],
    full_pearson_asymptotic_p = primary$pearson_p_asymptotic[[1L]],
    full_pearson_perm_p = primary$pearson_p_permutation_two_sided[[1L]],
    full_pearson_q = primary$pearson_q_exploratory[[1L]],
    full_spearman_rho = primary$spearman_rho[[1L]],
    full_spearman_asymptotic_p = primary$spearman_p_asymptotic[[1L]],
    full_spearman_perm_p = primary$spearman_p_permutation_two_sided[[1L]],
    full_spearman_q = primary$spearman_q_exploratory[[1L]],
    leave_one_out_all_same_sign = all(leave_one_out$sign_matches_full_pearson, na.rm = TRUE),
    leave_one_out_min_pearson_r = min(leave_one_out$pearson_r, na.rm = TRUE),
    leave_one_out_max_pearson_r = max(leave_one_out$pearson_r, na.rm = TRUE),
    leave_one_out_max_pearson_p = max(leave_one_out$pearson_p, na.rm = TRUE),
    leave_one_out_max_abs_delta_r = max(abs(leave_one_out$delta_pearson_r_from_full), na.rm = TRUE),
    n_influential_flags = sum(leave_one_out$influential_flag, na.rm = TRUE),
    n_bootstrap_requested = n_boot,
    n_bootstrap_finite_pearson = sum(is.finite(bootstrap$pearson_r)),
    n_bootstrap_finite_spearman = sum(is.finite(bootstrap$spearman_rho)),
    bootstrap_pearson_r_ci025 = essential_quantile(bootstrap$pearson_r, 0.025),
    bootstrap_pearson_r_ci975 = essential_quantile(bootstrap$pearson_r, 0.975),
    bootstrap_spearman_rho_ci025 = essential_quantile(bootstrap$spearman_rho, 0.025),
    bootstrap_spearman_rho_ci975 = essential_quantile(bootstrap$spearman_rho, 0.975),
    stringsAsFactors = FALSE
  )
}

essential_lm_terms <- function(fit, analysis, model_label, n) {
  coefficient_matrix <- summary(fit)$coefficients
  data.frame(
    compartment = "CellCycle",
    analysis = analysis,
    model = model_label,
    term = rownames(coefficient_matrix),
    estimate = coefficient_matrix[, 1L],
    std_error = coefficient_matrix[, 2L],
    statistic = coefficient_matrix[, 3L],
    p_value = coefficient_matrix[, 4L],
    n = n,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

essential_fit_shift_models <- function(treated, spec, tgi_column = "TGI_percent_auc") {
  local <- treated[
    is.finite(treated$ecdf_rmse) & is.finite(essential_safe_numeric(treated[[tgi_column]])),
    ,
    drop = FALSE
  ]
  group_column <- spec$group_column
  local[[group_column]] <- factor(as.character(local$analysis_group), levels = spec$group_levels)
  local$dose_mg_factor <- factor(local$dose_mg, levels = sort(unique(local$dose_mg)))
  outcome <- paste0("`", tgi_column, "`")
  formulas <- c(
    shift_only = paste0(outcome, " ~ scale(ecdf_rmse)"),
    shift_ploidy_dose_adjusted = paste0(
      outcome, " ~ scale(ecdf_rmse) + ", group_column, " + dose_mg_factor"
    ),
    shift_by_ploidy_dose_adjusted = paste0(
      outcome, " ~ scale(ecdf_rmse) * ", group_column, " + dose_mg_factor"
    )
  )
  models <- list()
  summaries <- list()
  for (analysis in names(formulas)) {
    fit <- stats::lm(stats::as.formula(formulas[[analysis]]), data = local)
    models[[analysis]] <- essential_lm_terms(fit, analysis, formulas[[analysis]], nrow(local))
    summary_fit <- summary(fit)
    summaries[[analysis]] <- data.frame(
      compartment = "CellCycle",
      analysis = analysis,
      model = formulas[[analysis]],
      n = stats::nobs(fit),
      residual_df = stats::df.residual(fit),
      r_squared = unname(summary_fit$r.squared),
      adjusted_r_squared = unname(summary_fit$adj.r.squared),
      sigma = unname(summary_fit$sigma),
      aic = stats::AIC(fit),
      stringsAsFactors = FALSE
    )
  }
  model_table <- do.call(rbind, models)
  model_table$tgi_measure <- tgi_column
  model_table$pre_specified_primary <- model_table$analysis == "shift_ploidy_dose_adjusted" &
    model_table$term == "scale(ecdf_rmse)"
  summary_table <- do.call(rbind, summaries)
  summary_table$tgi_measure <- tgi_column
  list(models = model_table, summaries = summary_table)
}

essential_confounding <- function(
  equal_shift,
  spec,
  n_perm,
  tgi_column = "TGI_percent_auc",
  ploidy_shift = equal_shift
) {
  data <- ploidy_shift
  data <- data[order(as.character(data$sample_id)), , drop = FALSE]
  data[[spec$group_numeric_column]] <- data$analysis_group_numeric
  ploidy_columns <- c(
    "mean_cell_ploidy", "median_cell_ploidy", "p90_cell_ploidy", spec$group_numeric_column
  )
  rows <- list()
  for (sample_set in c("all", "treated")) {
    local <- if (sample_set == "treated") data[data$dose_mg > 0, , drop = FALSE] else data
    for (shift_metric in c("ecdf_rmse", "signed_mean_shift")) {
      for (ploidy_column in ploidy_columns) {
        for (method in c("pearson", "spearman")) {
          asymptotic <- essential_safe_cor(local[[shift_metric]], local[[ploidy_column]], method)
          permutation <- essential_permutation_cor(
            local[[shift_metric]], local[[ploidy_column]], method, n_perm, exact = FALSE
          )
          rows[[length(rows) + 1L]] <- data.frame(
            compartment = "CellCycle",
            sample_set = sample_set,
            reference_type = as.character(local$reference_type[[1L]]),
            control_scope = as.character(local$control_scope[[1L]]),
            shift_metric = shift_metric,
            ploidy_measure = ploidy_column,
            method = method,
            n = asymptotic$n,
            estimate = asymptotic$estimate,
            p_value = asymptotic$p_value,
            p_permutation_two_sided = permutation$permutation_p_two_sided,
            n_permutations = permutation$n_permutations,
            permutation_mode = permutation$permutation_mode,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  ploidy_tests <- do.call(rbind, rows)

  treated <- equal_shift[equal_shift$dose_mg > 0, , drop = FALSE]
  centered <- treated
  centered$shift_centered <- centered$ecdf_rmse - ave(centered$ecdf_rmse, centered$dose_mg, FUN = mean)
  centered$tgi_centered <- centered[[tgi_column]] - ave(centered[[tgi_column]], centered$dose_mg, FUN = mean)
  centered_cor <- essential_safe_cor(centered$shift_centered, centered$tgi_centered, "pearson")
  centered_perm <- essential_permutation_cor(
    centered$shift_centered,
    centered$tgi_centered,
    "pearson",
    n_perm,
    exact = TRUE,
    strata = centered$dose_mg
  )

  adjusted <- treated[
    is.finite(treated$ecdf_rmse) & is.finite(essential_safe_numeric(treated[[tgi_column]])),
    ,
    drop = FALSE
  ]
  shift_residual <- stats::resid(stats::lm(
    ecdf_rmse ~ factor(dose_mg) + analysis_group_numeric,
    data = adjusted
  ))
  tgi_residual <- stats::resid(stats::lm(
    stats::as.formula(paste0("`", tgi_column, "` ~ factor(dose_mg) + analysis_group_numeric")),
    data = adjusted
  ))
  adjusted_cor <- essential_safe_cor(shift_residual, tgi_residual, "pearson")
  adjusted_perm <- essential_permutation_cor(shift_residual, tgi_residual, "pearson", n_perm, exact = TRUE)
  residualized <- rbind(
    data.frame(
      compartment = "CellCycle",
      analysis = paste0("residualized_against_dose_and_", spec$group_column),
      method = "pearson", n = adjusted_cor$n, estimate = adjusted_cor$estimate,
      tgi_measure = tgi_column,
      p_value = adjusted_cor$p_value,
      p_permutation_two_sided = adjusted_perm$permutation_p_two_sided,
      n_permutations = adjusted_perm$n_permutations,
      permutation_mode = adjusted_perm$permutation_mode,
      stringsAsFactors = FALSE
    ),
    data.frame(
      compartment = "CellCycle", analysis = "within_dose_centered", method = "pearson",
      tgi_measure = tgi_column,
      n = centered_cor$n, estimate = centered_cor$estimate, p_value = centered_cor$p_value,
      p_permutation_two_sided = centered_perm$permutation_p_two_sided,
      n_permutations = centered_perm$n_permutations,
      permutation_mode = centered_perm$permutation_mode,
      stringsAsFactors = FALSE
    )
  )
  list(ploidy_tests = ploidy_tests, residualized = residualized, centered = centered)
}

essential_plot_theme <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey70"),
      legend.key.height = grid::unit(0.45, "cm"),
      aspect.ratio = 1
    )
}

essential_dose_colors <- function() {
  c(
    "0mg/kg" = "#666666", "30mg/kg" = "#d95f02", "120mg/kg" = "#1b9e77",
    "treated" = "#377eb8", "2N" = "#4C78A8", "4N" = "#E45756"
  )
}

essential_group_colors <- function(spec) {
  if (identical(spec$type, "ETP")) {
    c("ETP-lower" = "#4C78A8", "ETP-higher" = "#E45756")
  } else {
    c("2N" = "#4C78A8", "4N" = "#E45756")
  }
}

essential_association_annotation <- function(association) {
  paste0(
    "Pearson r = ", essential_format_number(association$estimate[[1L]], 3),
    "\nPermutation P = ", essential_format_p(association$permutation_p_two_sided[[1L]]),
    "\nn = ", association$n[[1L]]
  )
}

essential_attach_association <- function(data, analysis_id, pearson, spearman) {
  data$analysis_id <- analysis_id
  data$pearson_r <- pearson$estimate[[1L]]
  data$pearson_p_asymptotic <- pearson$asymptotic_p[[1L]]
  data$pearson_p_permutation_two_sided <- pearson$permutation_p_two_sided[[1L]]
  data$pearson_p_permutation_positive <- pearson$permutation_p_positive[[1L]]
  data$pearson_n_permutations <- pearson$n_permutations[[1L]]
  data$pearson_permutation_mode <- pearson$permutation_mode[[1L]]
  data$spearman_rho <- spearman$estimate[[1L]]
  data$spearman_p_asymptotic <- spearman$asymptotic_p[[1L]]
  data$spearman_p_permutation_two_sided <- spearman$permutation_p_two_sided[[1L]]
  data$spearman_p_permutation_positive <- spearman$permutation_p_positive[[1L]]
  data$spearman_n_permutations <- spearman$n_permutations[[1L]]
  data$spearman_permutation_mode <- spearman$permutation_mode[[1L]]
  data$plot_annotation <- essential_association_annotation(pearson)
  data
}

essential_mean_ecdf_curves <- function(data, panel, grid) {
  rows <- lapply(unique(as.character(data$curve_label)), function(curve_label) {
    local <- data[as.character(data$curve_label) == curve_label, , drop = FALSE]
    sample_ids <- unique(local$sample_id)
    if (length(sample_ids) == 0L) return(NULL)
    matrix_values <- t(vapply(sample_ids, function(sample_id) {
      stats::ecdf(local$pseudotime[local$sample_id == sample_id])(grid)
    }, numeric(length(grid))))
    data.frame(
      panel = panel,
      pseudotime = grid,
      mean_ecdf = colMeans(matrix_values, na.rm = TRUE),
      curve_label = curve_label,
      color_group = as.character(local$color_group[[1L]]),
      line_group = as.character(local$line_group[[1L]]),
      n_samples = length(sample_ids),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

essential_direct_ecdf <- function(data, dose_tests, spec, n_perm) {
  data$treatment_group <- ifelse(data$dose_mg == 0, "0mg/kg", "treated")
  grid <- essential_make_grid(data)
  curve_rows <- list()
  add_curves <- function(local, panel, curve_label, color_group, line_group) {
    if (nrow(local) == 0L) return(invisible(NULL))
    local$curve_label <- curve_label
    local$color_group <- color_group
    local$line_group <- line_group
    curve_rows[[length(curve_rows) + 1L]] <<- essential_mean_ecdf_curves(local, panel, grid)
    invisible(NULL)
  }
  group_low <- spec$group_levels[[1L]]
  group_high <- spec$group_levels[[2L]]
  stratified_panel <- paste0("2. 0 vs treated, ", if (identical(spec$type, "ETP")) "ETP" else "initial-ploidy", "-stratified")
  panel_six <- paste0("6. Treated: ", group_low, " vs ", group_high)
  panel_seven <- paste0("7. Control: ", group_low, " vs ", group_high)
  panel_eight <- paste0("8. ", group_high, ": 0 vs treated")
  panel_nine <- paste0("9. ", group_low, ": 0 vs treated")
  panel_ten <- "10. 30 mg/kg: 2N vs 4N"
  panel_eleven <- "11. 120 mg/kg: 2N vs 4N"

  add_curves(data, "1. 0 vs treated", data$treatment_group, data$treatment_group, "All")
  add_curves(
    data,
    stratified_panel,
    paste(data$treatment_group, data$analysis_group, sep = " | "),
    data$treatment_group,
    as.character(data$analysis_group)
  )
  dose_panels <- list(
    list(panel = "3. 0 vs 30 mg/kg", doses = c(0, 30)),
    list(panel = "4. 0 vs 120 mg/kg", doses = c(0, 120)),
    list(panel = "5. 30 vs 120 mg/kg", doses = c(30, 120))
  )
  for (definition in dose_panels) {
    local <- data[data$dose_mg %in% definition$doses, , drop = FALSE]
    add_curves(local, definition$panel, local$dose, local$dose, "All")
  }
  treated <- data[data$dose_mg > 0, , drop = FALSE]
  control <- data[data$dose_mg == 0, , drop = FALSE]
  add_curves(treated, panel_six, treated$analysis_group, "treated", treated$analysis_group)
  add_curves(control, panel_seven, control$analysis_group, "0mg/kg", control$analysis_group)
  high <- data[data$analysis_group == group_high, , drop = FALSE]
  low <- data[data$analysis_group == group_low, , drop = FALSE]
  add_curves(high, panel_eight, high$treatment_group, high$treatment_group, group_high)
  add_curves(low, panel_nine, low$treatment_group, low$treatment_group, group_low)
  for (dose_definition in list(
    list(panel = panel_ten, dose_mg = 30),
    list(panel = panel_eleven, dose_mg = 120)
  )) {
    local <- data[data$dose_mg == dose_definition$dose_mg, , drop = FALSE]
    add_curves(
      local,
      dose_definition$panel,
      as.character(local$initial_ploidy),
      as.character(local$initial_ploidy),
      "All"
    )
  }
  curves <- do.call(rbind, curve_rows)

  get_test <- function(comparison, stratified) {
    hit <- dose_tests[
      dose_tests$comparison == comparison & dose_tests$stratified_by_analysis_group == stratified,
      ,
      drop = FALSE
    ]
    if (nrow(hit) == 0L) data.frame() else hit[1L, , drop = FALSE]
  }
  between_group <- function(local) {
    essential_ecdf_group_test(local, "analysis_group", group_low, group_high, FALSE, n_perm)
  }
  treatment_within_group <- function(group) {
    local <- data[data$analysis_group == group, , drop = FALSE]
    essential_ecdf_group_test(local, "treatment_group", "0mg/kg", "treated", FALSE, n_perm)
  }
  initial_ploidy_within_dose <- function(dose_mg) {
    local <- data[data$dose_mg == dose_mg, , drop = FALSE]
    essential_ecdf_group_test(local, "initial_ploidy", "2N", "4N", FALSE, n_perm)
  }
  test_rows <- list(
    cbind(panel = "1. 0 vs treated", test_design = "sample_label_permutation", get_test("0_vs_30plus120", FALSE)),
    cbind(panel = stratified_panel, test_design = "treatment_label_permutation_within_group", get_test("0_vs_30plus120", TRUE)),
    cbind(panel = "3. 0 vs 30 mg/kg", test_design = "sample_label_permutation", get_test("0_vs_30", FALSE)),
    cbind(panel = "4. 0 vs 120 mg/kg", test_design = "sample_label_permutation", get_test("0_vs_120", FALSE)),
    cbind(panel = "5. 30 vs 120 mg/kg", test_design = "sample_label_permutation", get_test("30_vs_120", FALSE)),
    cbind(panel = panel_six, test_design = "unpaired_group_label_permutation", between_group(treated)),
    cbind(panel = panel_seven, test_design = "unpaired_group_label_permutation", between_group(control)),
    cbind(panel = panel_eight, test_design = "treatment_label_permutation_within_group", treatment_within_group(group_high)),
    cbind(panel = panel_nine, test_design = "treatment_label_permutation_within_group", treatment_within_group(group_low)),
    cbind(panel = panel_ten, test_design = "initial_ploidy_label_permutation_within_dose", initial_ploidy_within_dose(30)),
    cbind(panel = panel_eleven, test_design = "initial_ploidy_label_permutation_within_dose", initial_ploidy_within_dose(120))
  )
  all_columns <- unique(unlist(lapply(test_rows, names), use.names = FALSE))
  test_rows <- lapply(test_rows, function(row) {
    for (column in setdiff(all_columns, names(row))) row[[column]] <- NA
    row[, all_columns, drop = FALSE]
  })
  tests <- do.call(rbind, test_rows)
  tests$method <- spec$method
  tests$q_ecdf_rmse_exploratory <- stats::p.adjust(tests$p_ecdf_rmse, method = "BH")
  tests$q_ecdf_ks_exploratory <- stats::p.adjust(tests$p_ecdf_ks, method = "BH")
  tests$q_ecdf_mean_abs_exploratory <- stats::p.adjust(tests$p_ecdf_mean_abs, method = "BH")
  tests$threshold <- spec$threshold
  tests$estimable <- is.finite(tests$p_ecdf_rmse)
  tests$annotation <- paste0(
    tests$group_a, " n=", tests$n_group_a, "; ", tests$group_b, " n=", tests$n_group_b,
    "\nRMSE=", essential_format_number(tests$observed_ecdf_rmse, 4), ", P=", essential_format_p(tests$p_ecdf_rmse),
    "\nKS=", essential_format_number(tests$observed_ecdf_ks, 4), ", P=", essential_format_p(tests$p_ecdf_ks)
  )
  panel_levels <- c(
    "1. 0 vs treated", stratified_panel, "3. 0 vs 30 mg/kg", "4. 0 vs 120 mg/kg",
    "5. 30 vs 120 mg/kg", panel_six, panel_seven, panel_eight, panel_nine,
    panel_ten, panel_eleven
  )
  curves$panel <- factor(curves$panel, levels = panel_levels)
  curves$color_group <- factor(curves$color_group, levels = names(essential_dose_colors()))
  curves$line_group <- factor(curves$line_group, levels = c("All", spec$group_levels))
  tests$panel <- factor(tests$panel, levels = panel_levels)
  plot_data <- merge(curves, tests, by = "panel", all.x = TRUE, sort = FALSE)
  plot_data$method <- spec$method
  line_values <- c("All" = "solid")
  line_values[spec$group_levels] <- c("solid", "22")
  plot <- ggplot2::ggplot(
    curves,
    ggplot2::aes(x = pseudotime, y = mean_ecdf, color = color_group, linetype = line_group, group = curve_label)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_text(
      data = tests,
      ggplot2::aes(x = 0.98, y = 0.05, label = annotation),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 0,
      size = 2.35,
      lineheight = 0.92
    ) +
    ggplot2::facet_wrap(~panel, ncol = 3L, drop = FALSE) +
    ggplot2::scale_color_manual(values = essential_dose_colors(), name = "Group") +
    ggplot2::scale_linetype_manual(values = line_values, name = spec$group_label) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.05), clip = "off") +
    ggplot2::labs(
      title = paste0("CellCycle: direct group mean ECDF comparisons (", spec$method, ")"),
      subtitle = "Equal-sample mean ECDFs",
      x = "Pseudotime",
      y = "Mean ECDF"
    ) +
    essential_plot_theme() +
    ggplot2::theme(legend.position = "bottom", aspect.ratio = NULL)
  list(plot = plot, tests = tests, plot_data = plot_data)
}

essential_scatter_plot <- function(
  plot_data,
  x_column,
  y_column,
  color_column,
  shape_column,
  color_values,
  title,
  x_label,
  y_label,
  color_label,
  shape_label,
  annotation,
  label_nudge_y = 2.5
) {
  label_layer <- ggrepel::geom_text_repel(
    data = plot_data,
    ggplot2::aes(
      x = .data[[x_column]],
      y = .data[[y_column]],
      label = sample_id,
      color = .data[[color_column]]
    ),
    inherit.aes = FALSE,
    nudge_y = label_nudge_y,
    size = 2.4,
    seed = 1,
    box.padding = 0.22,
    point.padding = 0.12,
    min.segment.length = 0,
    segment.color = "grey65",
    segment.size = 0.25,
    max.overlaps = Inf,
    show.legend = FALSE
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data[[x_column]], y = .data[[y_column]], color = .data[[color_column]], shape = .data[[shape_column]])
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::geom_smooth(
      data = plot_data,
      ggplot2::aes(x = .data[[x_column]], y = .data[[y_column]], group = 1),
      inherit.aes = FALSE,
      method = "lm",
      se = TRUE,
      color = "black",
      linewidth = 0.55
    ) +
    ggplot2::geom_point(size = 2.9) +
    label_layer +
    ggplot2::annotate(
      "label",
      x = Inf,
      y = Inf,
      label = annotation,
      hjust = 1.05,
      vjust = 1.1,
      size = 3.0,
      linewidth = 0.2
    ) +
    ggplot2::scale_color_manual(values = color_values, name = color_label) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.15, 0.12))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.15))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = title, x = x_label, y = y_label, shape = shape_label) +
    essential_plot_theme() +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 14, 8, 14))
}

essential_dose_adjusted_group_effect <- function(tgi_value, groups, dose_mg, group_high) {
  keep <- is.finite(tgi_value) & is.finite(dose_mg) & !is.na(groups)
  tgi_value <- tgi_value[keep]
  groups <- as.character(groups[keep])
  dose_mg <- dose_mg[keep]
  if (length(tgi_value) < 3L || !group_high %in% groups || length(unique(groups)) < 2L) {
    return(NA_real_)
  }
  group_high_indicator <- as.integer(groups == group_high)
  design <- stats::model.matrix(~factor(dose_mg) + group_high_indicator)
  fit <- stats::lm.fit(design, tgi_value)
  coefficient <- unname(fit$coefficients[["group_high_indicator"]])
  if (length(coefficient) != 1L || !is.finite(coefficient)) NA_real_ else coefficient
}

essential_tgi_group_comparison <- function(data, spec, tgi_spec, n_perm = 10000L) {
  required <- c("sample_id", "dose_mg", "analysis_group", "tgi_value")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("TGI group comparison is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  local <- data[
    essential_safe_numeric(data$dose_mg) > 0 &
      is.finite(essential_safe_numeric(data$tgi_value)) &
      as.character(data$analysis_group) %in% spec$group_levels,
    ,
    drop = FALSE
  ]
  local$dose_mg <- essential_safe_numeric(local$dose_mg)
  local$tgi_value <- essential_safe_numeric(local$tgi_value)
  local$analysis_group <- factor(as.character(local$analysis_group), levels = spec$group_levels)
  local <- local[order(local$dose_mg, local$analysis_group, local$sample_id), , drop = FALSE]

  group_low <- spec$group_levels[[1L]]
  group_high <- spec$group_levels[[2L]]
  groups <- as.character(local$analysis_group)
  if (!all(c(group_low, group_high) %in% groups)) {
    stop(
      "TGI group comparison requires treated samples in both groups: ",
      paste(spec$group_levels, collapse = " vs "),
      call. = FALSE
    )
  }
  low_values <- local$tgi_value[groups == group_low]
  high_values <- local$tgi_value[groups == group_high]
  observed <- essential_dose_adjusted_group_effect(
    local$tgi_value,
    groups,
    local$dose_mg,
    group_high
  )

  welch_p <- tryCatch(
    stats::t.test(high_values, low_values, paired = FALSE)$p.value,
    error = function(e) NA_real_
  )
  wilcoxon_p <- tryCatch(
    suppressWarnings(stats::wilcox.test(high_values, low_values, paired = FALSE, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )

  strata <- as.character(local$dose_mg)
  split_indices <- split(seq_along(groups), strata)
  exact_count <- prod(vapply(split_indices, function(indices) {
    choose(length(indices), sum(groups[indices] == group_low))
  }, numeric(1L)))
  if (is.finite(exact_count) && exact_count <= 400000) {
    permuted_groups <- essential_exact_group_vectors(groups, group_low, group_high, strata)
    permuted_effects <- vapply(permuted_groups, function(permuted_group) {
      essential_dose_adjusted_group_effect(
        local$tgi_value,
        permuted_group,
        local$dose_mg,
        group_high
      )
    }, numeric(1L))
    finite_permutations <- is.finite(permuted_effects)
    permutation_p <- if (is.finite(observed) && any(finite_permutations)) {
      mean(abs(permuted_effects[finite_permutations]) >= abs(observed) - 1e-15)
    } else {
      NA_real_
    }
    permutation_mode <- "exact_group_label_enumeration_within_dose"
    n_permutations <- sum(finite_permutations)
  } else {
    permuted_effects <- replicate(n_perm, {
      permuted_group <- groups
      for (indices in split_indices) permuted_group[indices] <- sample(groups[indices])
      essential_dose_adjusted_group_effect(
        local$tgi_value,
        permuted_group,
        local$dose_mg,
        group_high
      )
    })
    finite_permutations <- is.finite(permuted_effects)
    exceed <- if (is.finite(observed)) {
      sum(abs(permuted_effects[finite_permutations]) >= abs(observed) - 1e-15)
    } else {
      NA_real_
    }
    permutation_p <- if (is.finite(exceed)) {
      (exceed + 1) / (sum(finite_permutations) + 1)
    } else {
      NA_real_
    }
    permutation_mode <- "monte_carlo_group_label_permutation_within_dose"
    n_permutations <- sum(finite_permutations)
  }

  data.frame(
    compartment = "CellCycle",
    sample_set = "treated",
    comparison = paste0(group_high, "_minus_", group_low),
    comparison_design = "independent_tumors_with_group_labels_permuted_within_dose",
    pairing_status = "not_paired_no_one_to_one_mouse_key",
    group_label = spec$group_label,
    group_low = group_low,
    group_high = group_high,
    tgi_measure = tgi_spec$measure,
    tgi_outcome = tgi_spec$outcome,
    tgi_control_summary = tgi_spec$control_summary,
    tgi_day = tgi_spec$day,
    n = nrow(local),
    n_group_low = length(low_values),
    n_group_high = length(high_values),
    mean_group_low = mean(low_values),
    mean_group_high = mean(high_values),
    mean_difference_high_minus_low = mean(high_values) - mean(low_values),
    median_group_low = stats::median(low_values),
    median_group_high = stats::median(high_values),
    median_difference_high_minus_low = stats::median(high_values) - stats::median(low_values),
    dose_adjusted_difference_high_minus_low = observed,
    welch_t_p_unpaired = welch_p,
    wilcoxon_rank_sum_p_unpaired = wilcoxon_p,
    permutation_p_two_sided = permutation_p,
    n_permutations = n_permutations,
    permutation_mode = permutation_mode,
    permutation_strata = "dose_mg",
    stringsAsFactors = FALSE
  )
}

essential_tgi_group_boxplot <- function(data, comparison, spec, tgi_spec) {
  required_data <- c("sample_id", "dose", "dose_mg", "analysis_group", "tgi_value")
  required_stats <- c(
    "group_low", "group_high", "n_group_low", "n_group_high",
    "dose_adjusted_difference_high_minus_low", "permutation_p_two_sided",
    "permutation_mode", "pairing_status"
  )
  missing_data <- setdiff(required_data, names(data))
  missing_stats <- setdiff(required_stats, names(comparison))
  if (length(missing_data) > 0L) {
    stop("TGI group boxplot data are missing: ", paste(missing_data, collapse = ", "), call. = FALSE)
  }
  if (length(missing_stats) > 0L) {
    stop("TGI group boxplot statistics are missing: ", paste(missing_stats, collapse = ", "), call. = FALSE)
  }
  data <- essential_restore_scatter_factors(data, spec)
  effect <- essential_safe_numeric(comparison$dose_adjusted_difference_high_minus_low[[1L]])
  p_value <- essential_safe_numeric(comparison$permutation_p_two_sided[[1L]])
  group_low <- as.character(comparison$group_low[[1L]])
  group_high <- as.character(comparison$group_high[[1L]])
  test_label <- if (grepl("^exact_", comparison$permutation_mode[[1L]])) {
    "Exact dose-stratified permutation"
  } else {
    "Monte Carlo dose-stratified permutation"
  }
  annotation <- paste0(
    "Independent tumors; no one-to-one mouse pairing\n",
    test_label, "\n",
    "Adjusted delta (", group_high, " - ", group_low, ") = ",
    essential_format_number(effect, 3L), " percentage points; P = ", essential_format_p(p_value), "\n",
    "n = ", comparison$n_group_low[[1L]], " (", group_low, ") and ",
    comparison$n_group_high[[1L]], " (", group_high, ")"
  )
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = analysis_group, y = tgi_value, fill = analysis_group)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::geom_boxplot(width = 0.58, alpha = 0.35, outlier.shape = NA, linewidth = 0.55) +
    ggplot2::geom_point(
      ggplot2::aes(color = dose, shape = dose),
      position = ggplot2::position_jitter(width = 0.075, height = 0, seed = 1),
      size = 3.1,
      stroke = 0.7
    ) +
    ggplot2::scale_fill_manual(values = essential_group_colors(spec), guide = "none") +
    ggplot2::scale_color_manual(values = essential_dose_colors(), name = "Dose") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.12))) +
    ggplot2::labs(
      title = paste0("Treated CellCycle tumors: ", tgi_spec$title_label, " by ", spec$group_label),
      subtitle = annotation,
      x = spec$group_label,
      y = tgi_spec$axis_label,
      shape = "Dose"
    ) +
    essential_plot_theme() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey25", lineheight = 1.05)
    )
}

essential_figure_filenames <- function() {
  c(
    "CellCycle_direct_group_ecdf_comparisons.pdf",
    "CellCycle_direct_group_ecdf_comparisons_selected_3panel.pdf",
    "CellCycle_TGI_group_boxplot.pdf",
    "CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref.pdf",
    "CellCycle_TGI_AUC_vs_mean_ETP.pdf",
    "CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose.pdf",
    "CellCycle_TGI_association_within_dose_centered.pdf",
    "CellCycle_ecdf_rmse_vs_ploidy.pdf"
  )
}

essential_read_figure_table <- function(path, required_columns, label) {
  if (!file.exists(path)) stop("Missing figure-only table: ", path, call. = FALSE)
  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing <- setdiff(required_columns, names(data))
  if (length(missing) > 0L) {
    stop(
      label, " is missing required columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(data) == 0L) stop(label, " contains no rows: ", path, call. = FALSE)
  data
}

essential_require_single_row <- function(data, keep, label) {
  keep[is.na(keep)] <- FALSE
  selected <- data[keep, , drop = FALSE]
  if (nrow(selected) != 1L) {
    stop(label, " must select exactly one statistics row; found ", nrow(selected), call. = FALSE)
  }
  selected
}

essential_prefixed_association <- function(row, prefix) {
  estimate_column <- if (identical(prefix, "pearson")) "pearson_r" else "spearman_rho"
  required <- c(
    "n", estimate_column, paste0(prefix, "_p_asymptotic"),
    paste0(prefix, "_p_permutation_two_sided"),
    paste0(prefix, "_p_permutation_positive"),
    paste0(prefix, "_n_permutations"), paste0(prefix, "_permutation_mode")
  )
  missing <- setdiff(required, names(row))
  if (length(missing) > 0L) {
    stop("Statistics row is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  data.frame(
    n = essential_safe_numeric(row$n[[1L]]),
    estimate = essential_safe_numeric(row[[estimate_column]][[1L]]),
    asymptotic_p = essential_safe_numeric(row[[paste0(prefix, "_p_asymptotic")]][[1L]]),
    permutation_p_two_sided = essential_safe_numeric(
      row[[paste0(prefix, "_p_permutation_two_sided")]][[1L]]
    ),
    permutation_p_positive = essential_safe_numeric(
      row[[paste0(prefix, "_p_permutation_positive")]][[1L]]
    ),
    n_permutations = essential_safe_numeric(row[[paste0(prefix, "_n_permutations")]][[1L]]),
    permutation_mode = as.character(row[[paste0(prefix, "_permutation_mode")]][[1L]]),
    stringsAsFactors = FALSE
  )
}

essential_generic_association <- function(row) {
  required <- c(
    "n", "estimate", "p_value", "p_permutation_two_sided",
    "n_permutations", "permutation_mode"
  )
  missing <- setdiff(required, names(row))
  if (length(missing) > 0L) {
    stop("Statistics row is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  data.frame(
    n = essential_safe_numeric(row$n[[1L]]),
    estimate = essential_safe_numeric(row$estimate[[1L]]),
    asymptotic_p = essential_safe_numeric(row$p_value[[1L]]),
    permutation_p_two_sided = essential_safe_numeric(row$p_permutation_two_sided[[1L]]),
    permutation_p_positive = NA_real_,
    n_permutations = essential_safe_numeric(row$n_permutations[[1L]]),
    permutation_mode = as.character(row$permutation_mode[[1L]]),
    stringsAsFactors = FALSE
  )
}

essential_plot_data_association <- function(data, prefix) {
  row <- data[1L, , drop = FALSE]
  row$n <- nrow(data)
  essential_prefixed_association(row, prefix)
}

essential_restore_scatter_factors <- function(data, spec) {
  required <- c("sample_id", "dose", "dose_mg", "analysis_group")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Scatter plot data are missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  dose_order <- unique(as.character(data$dose[order(essential_safe_numeric(data$dose_mg))]))
  data$dose <- factor(as.character(data$dose), levels = dose_order)
  data$analysis_group <- factor(as.character(data$analysis_group), levels = spec$group_levels)
  data
}

essential_select_direct_panels <- function(curves, tests) {
  coordinates <- data.frame(
    row = c(1L, 3L, 3L),
    column = c(1L, 2L, 3L)
  )
  panel_indices <- (coordinates$row - 1L) * 3L + coordinates$column
  panel_levels <- unique(as.character(tests$panel))
  if (length(panel_levels) < max(panel_indices)) {
    stop(
      "Direct ECDF tables do not contain the requested panels at positions: ",
      paste(paste0(coordinates$row, ",", coordinates$column), collapse = "; "),
      call. = FALSE
    )
  }
  selected_levels <- panel_levels[panel_indices]
  selected_tests <- tests[match(selected_levels, as.character(tests$panel)), , drop = FALSE]
  selected_curves <- curves[
    as.character(curves$panel) %in% selected_levels,
    ,
    drop = FALSE
  ]
  selected_curves$panel <- factor(as.character(selected_curves$panel), levels = selected_levels)
  selected_tests$panel <- factor(as.character(selected_tests$panel), levels = selected_levels)
  list(
    curves = selected_curves,
    tests = selected_tests,
    panel_levels = selected_levels,
    coordinates = coordinates
  )
}

essential_direct_plot_from_tables <- function(
  curves,
  tests,
  spec,
  title = paste0("CellCycle: direct group mean ECDF comparisons (", spec$method, ")")
) {
  curve_columns <- c("panel", "pseudotime", "mean_ecdf", "curve_label", "color_group", "line_group")
  test_columns <- c("panel", "annotation")
  curve_missing <- setdiff(curve_columns, names(curves))
  test_missing <- setdiff(test_columns, names(tests))
  if (length(curve_missing) > 0L) {
    stop("Direct ECDF plot data are missing columns: ", paste(curve_missing, collapse = ", "), call. = FALSE)
  }
  if (length(test_missing) > 0L) {
    stop("Direct ECDF statistics are missing columns: ", paste(test_missing, collapse = ", "), call. = FALSE)
  }
  panel_levels <- unique(as.character(tests$panel))
  curves$panel <- factor(as.character(curves$panel), levels = panel_levels)
  curves$color_group <- factor(as.character(curves$color_group), levels = names(essential_dose_colors()))
  curves$line_group <- factor(as.character(curves$line_group), levels = c("All", spec$group_levels))
  tests$panel <- factor(as.character(tests$panel), levels = panel_levels)
  line_values <- c("All" = "solid")
  line_values[spec$group_levels] <- c("solid", "22")
  ggplot2::ggplot(
    curves,
    ggplot2::aes(
      x = pseudotime, y = mean_ecdf, color = color_group,
      linetype = line_group, group = curve_label
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_text(
      data = tests,
      ggplot2::aes(x = 0.98, y = 0.05, label = annotation),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 0,
      size = 2.35,
      lineheight = 0.92
    ) +
    ggplot2::facet_wrap(~panel, ncol = 3L, drop = FALSE) +
    ggplot2::scale_color_manual(values = essential_dose_colors(), name = "Group") +
    ggplot2::scale_linetype_manual(values = line_values, name = spec$group_label) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.05), clip = "off") +
    ggplot2::labs(
      title = title,
      subtitle = "Equal-sample mean ECDFs",
      x = "Pseudotime",
      y = "Mean ECDF"
    ) +
    essential_plot_theme() +
    ggplot2::theme(legend.position = "bottom", aspect.ratio = NULL)
}

essential_regenerate_figures_from_tables <- function(tables_root, output_root, spec) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required", call. = FALSE)
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("ggrepel is required", call. = FALSE)

  plot_data_dir <- file.path(tables_root, "plot_data", spec$method)
  stats_dir <- file.path(tables_root, "stats", spec$method)
  figures_dir <- essential_ensure_dir(file.path(output_root, "Figures", spec$method))
  read_plot_data <- function(stem, required) {
    essential_read_figure_table(
      file.path(plot_data_dir, paste0(stem, "_plot_data.csv")),
      required,
      paste0(stem, " plot data")
    )
  }
  read_stats <- function(filename, required) {
    essential_read_figure_table(
      file.path(stats_dir, filename), required, paste0(filename, " statistics")
    )
  }

  association_stats <- read_stats(
    "CellCycle_TGI_associations_ecdf_rmse.csv",
    c(
      "reference_type", "shift_metric", "tgi_measure", "n", "pearson_r",
      "pearson_p_asymptotic", "pearson_p_permutation_two_sided",
      "pearson_p_permutation_positive", "pearson_n_permutations",
      "pearson_permutation_mode", "spearman_rho", "spearman_p_asymptotic",
      "spearman_p_permutation_two_sided", "spearman_p_permutation_positive",
      "spearman_n_permutations", "spearman_permutation_mode", "pre_specified_primary"
    )
  )
  primary_stats <- essential_require_single_row(
    association_stats,
    as.character(association_stats$pre_specified_primary) %in% c("TRUE", "T", "1"),
    "Primary ECDF-TGI association"
  )
  tgi_spec <- essential_tgi_spec_from_measure(primary_stats$tgi_measure[[1L]])
  primary_pearson <- essential_prefixed_association(primary_stats, "pearson")
  primary_spearman <- essential_prefixed_association(primary_stats, "spearman")

  primary_data <- read_plot_data(
    "CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref",
    c("sample_id", "ecdf_rmse", "dose", "dose_mg", "analysis_group")
  )
  primary_data <- essential_add_tgi_plot_metadata(primary_data, tgi_spec)
  primary_data <- essential_restore_scatter_factors(primary_data, spec)
  primary_data <- essential_attach_association(
    primary_data,
    "CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref",
    primary_pearson,
    primary_spearman
  )
  primary_plot <- essential_scatter_plot(
    primary_data,
    "ecdf_rmse",
    "tgi_value",
    "dose",
    "analysis_group",
    essential_dose_colors(),
    paste0(
      "Cell-cycle-associated tumor cells: ", tgi_spec$title_label,
      " vs sample-equal ECDF shift"
    ),
    "ECDF RMSE from group-matched equal-sample 0 mg/kg reference",
    tgi_spec$axis_label,
    "Dose",
    spec$group_label,
    primary_data$plot_annotation[[1L]]
  )
  essential_save_pdf(
    primary_plot,
    file.path(figures_dir, "CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref.pdf"),
    6.8,
    6.8
  )

  group_comparison_all <- read_stats(
    "CellCycle_TGI_group_comparison.csv",
    c(
      "sample_set", "tgi_measure", "group_low", "group_high", "n_group_low",
      "n_group_high", "dose_adjusted_difference_high_minus_low",
      "permutation_p_two_sided", "permutation_mode", "pairing_status"
    )
  )
  group_comparison <- essential_require_single_row(
    group_comparison_all,
    group_comparison_all$sample_set == "treated" &
      group_comparison_all$tgi_measure == tgi_spec$measure,
    "Treated TGI group comparison"
  )
  group_plot_data <- read_plot_data(
    "CellCycle_TGI_group_boxplot",
    c(
      "sample_id", "dose", "dose_mg", "analysis_group", "tgi_value",
      "tgi_measure"
    )
  )
  if (length(unique(as.character(group_plot_data$tgi_measure))) != 1L ||
      !identical(as.character(group_plot_data$tgi_measure[[1L]]), tgi_spec$measure)) {
    stop("TGI group boxplot data do not match the selected TGI measure", call. = FALSE)
  }
  group_plot <- essential_tgi_group_boxplot(
    group_plot_data,
    group_comparison,
    spec,
    tgi_spec
  )
  essential_save_pdf(
    group_plot,
    file.path(figures_dir, "CellCycle_TGI_group_boxplot.pdf"),
    6.8,
    6.4
  )

  direct_data <- read_plot_data(
    "CellCycle_direct_group_ecdf_comparisons",
    c("panel", "pseudotime", "mean_ecdf", "curve_label", "color_group", "line_group")
  )
  direct_stats <- read_stats(
    "CellCycle_direct_group_ecdf_comparisons_11panel_tests.csv",
    c("panel", "annotation", "p_ecdf_rmse", "p_ecdf_ks")
  )
  direct_plot <- essential_direct_plot_from_tables(direct_data, direct_stats, spec)
  essential_save_pdf(
    direct_plot,
    file.path(figures_dir, "CellCycle_direct_group_ecdf_comparisons.pdf"),
    15,
    15
  )
  selected_direct <- essential_select_direct_panels(direct_data, direct_stats)
  selected_direct_plot <- essential_direct_plot_from_tables(
    selected_direct$curves,
    selected_direct$tests,
    spec,
    title = paste0(
      "CellCycle: selected direct group mean ECDF comparisons (", spec$method, ")"
    )
  )
  essential_save_pdf(
    selected_direct_plot,
    file.path(
      figures_dir,
      "CellCycle_direct_group_ecdf_comparisons_selected_3panel.pdf"
    ),
    15,
    5.5
  )

  mean_stats_all <- read_stats(
    "CellCycle_TGI_associations_mean_ETP.csv",
    c(
      "tgi_measure", "n", "pearson_r", "pearson_p_asymptotic",
      "pearson_p_permutation_two_sided", "pearson_p_permutation_positive",
      "pearson_n_permutations", "pearson_permutation_mode", "spearman_rho",
      "spearman_p_asymptotic", "spearman_p_permutation_two_sided",
      "spearman_p_permutation_positive", "spearman_n_permutations",
      "spearman_permutation_mode"
    )
  )
  mean_stats <- essential_require_single_row(
    mean_stats_all,
    mean_stats_all$tgi_measure == tgi_spec$measure,
    "Mean ETP-TGI association"
  )
  mean_data <- read_plot_data(
    "CellCycle_TGI_AUC_vs_mean_ETP",
    c(
      "sample_id", "sample_mean_endpoint_ploidy", "dose", "dose_mg", "analysis_group"
    )
  )
  mean_data <- essential_add_tgi_plot_metadata(mean_data, tgi_spec)
  mean_data <- essential_restore_scatter_factors(mean_data, spec)
  mean_data <- essential_attach_association(
    mean_data,
    "CellCycle_TGI_AUC_vs_mean_ETP",
    essential_prefixed_association(mean_stats, "pearson"),
    essential_prefixed_association(mean_stats, "spearman")
  )
  mean_plot <- essential_scatter_plot(
    mean_data,
    "sample_mean_endpoint_ploidy",
    "tgi_value",
    "dose",
    "analysis_group",
    essential_dose_colors(),
    paste0(
      "Cell-cycle-associated tumor cells: ", tgi_spec$title_label,
      " vs sample mean ETP"
    ),
    "Sample mean ETP",
    tgi_spec$axis_label,
    "Dose",
    spec$group_label,
    mean_data$plot_annotation[[1L]]
  )
  essential_save_pdf(
    mean_plot,
    file.path(figures_dir, "CellCycle_TGI_AUC_vs_mean_ETP.pdf"),
    6.8,
    6.8
  )

  model_data <- read_plot_data(
    "CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose",
    c("sample_id", "ecdf_rmse", "dose", "dose_mg", "analysis_group")
  )
  model_data <- essential_add_tgi_plot_metadata(model_data, tgi_spec)
  model_data <- essential_restore_scatter_factors(model_data, spec)
  model_data <- essential_attach_association(
    model_data,
    "CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose",
    primary_pearson,
    primary_spearman
  )
  model_plot <- essential_scatter_plot(
    model_data,
    "ecdf_rmse",
    "tgi_value",
    "analysis_group",
    "dose",
    essential_group_colors(spec),
    paste0(
      "CellCycle: ", tgi_spec$title_label,
      " association with pseudotime shift and ploidy"
    ),
    "ECDF RMSE from group-matched 0 mg/kg reference",
    tgi_spec$axis_label,
    spec$group_label,
    "Dose",
    model_data$plot_annotation[[1L]]
  )
  essential_save_pdf(
    model_plot,
    file.path(figures_dir, "CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose.pdf"),
    6.8,
    6.8
  )

  centered_stats_all <- read_stats(
    "residualized_TGI_associations.csv",
    c(
      "analysis", "method", "n", "estimate", "p_value",
      "p_permutation_two_sided", "n_permutations", "permutation_mode"
    )
  )
  centered_stats <- essential_require_single_row(
    centered_stats_all,
    centered_stats_all$analysis == "within_dose_centered" & centered_stats_all$method == "pearson",
    "Within-dose-centered TGI association"
  )
  centered_data <- read_plot_data(
    "CellCycle_TGI_association_within_dose_centered",
    c(
      "sample_id", "shift_centered", "tgi_centered", "dose", "dose_mg",
      "analysis_group", "spearman_rho", "spearman_p_asymptotic",
      "spearman_p_permutation_two_sided", "spearman_p_permutation_positive",
      "spearman_n_permutations", "spearman_permutation_mode"
    )
  )
  centered_data <- essential_restore_scatter_factors(centered_data, spec)
  centered_data <- essential_attach_association(
    centered_data,
    "CellCycle_TGI_association_within_dose_centered",
    essential_generic_association(centered_stats),
    essential_plot_data_association(centered_data, "spearman")
  )
  centered_plot <- essential_scatter_plot(
    centered_data,
    "shift_centered",
    "tgi_centered",
    "dose",
    "analysis_group",
    essential_dose_colors(),
    "CellCycle TGI association after within-dose centering",
    "Dose-centered ECDF RMSE",
    paste0("Dose-centered ", tgi_spec$axis_label),
    "Dose",
    spec$group_label,
    centered_data$plot_annotation[[1L]],
    label_nudge_y = 1.0
  ) +
    ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35)
  essential_save_pdf(
    centered_plot,
    file.path(figures_dir, "CellCycle_TGI_association_within_dose_centered.pdf"),
    6.6,
    6.6
  )

  ploidy_stats_all <- read_stats(
    "ploidy_confounding_tests.csv",
    c(
      "sample_set", "shift_metric", "ploidy_measure", "method", "n",
      "estimate", "p_value", "p_permutation_two_sided", "n_permutations",
      "permutation_mode"
    )
  )
  ploidy_primary <- ploidy_stats_all[
    ploidy_stats_all$sample_set == "treated" &
      ploidy_stats_all$shift_metric == "ecdf_rmse" &
      ploidy_stats_all$ploidy_measure == "mean_cell_ploidy",
    ,
    drop = FALSE
  ]
  ploidy_pearson <- essential_require_single_row(
    ploidy_primary, ploidy_primary$method == "pearson", "Ploidy Pearson association"
  )
  ploidy_spearman <- essential_require_single_row(
    ploidy_primary, ploidy_primary$method == "spearman", "Ploidy Spearman association"
  )
  ploidy_data <- read_plot_data(
    "CellCycle_ecdf_rmse_vs_ploidy",
    c(
      "sample_id", "mean_cell_ploidy", "ecdf_rmse", "dose", "dose_mg",
      "analysis_group", "initial_ploidy"
    )
  )
  ploidy_data <- ploidy_data[
    essential_safe_numeric(ploidy_data$dose_mg) > 0 &
      is.finite(essential_safe_numeric(ploidy_data$mean_cell_ploidy)) &
      is.finite(essential_safe_numeric(ploidy_data$ecdf_rmse)),
    ,
    drop = FALSE
  ]
  ploidy_data <- ploidy_data[
    order(
      essential_safe_numeric(ploidy_data$dose_mg),
      as.character(ploidy_data$initial_ploidy),
      as.character(ploidy_data$sample_id)
    ),
    ,
    drop = FALSE
  ]
  ploidy_data <- essential_restore_scatter_factors(ploidy_data, spec)
  ploidy_data$initial_ploidy <- factor(
    as.character(ploidy_data$initial_ploidy),
    levels = c("2N", "4N")
  )
  ploidy_data <- essential_attach_association(
    ploidy_data,
    "CellCycle_ecdf_rmse_vs_ploidy",
    essential_generic_association(ploidy_pearson),
    essential_generic_association(ploidy_spearman)
  )
  ploidy_plot <- essential_scatter_plot(
    ploidy_data,
    "mean_cell_ploidy",
    "ecdf_rmse",
    "dose",
    "initial_ploidy",
    essential_dose_colors(),
    "Treated CellCycle ECDF RMSE vs mean cell ploidy",
    "Mean cell ploidy",
    "ECDF RMSE",
    "Dose",
    "Initial ploidy",
    ploidy_data$plot_annotation[[1L]],
    label_nudge_y = 0.004
  )
  essential_save_pdf(
    ploidy_plot,
    file.path(figures_dir, "CellCycle_ecdf_rmse_vs_ploidy.pdf"),
    6.4,
    6.4
  )

  invisible(list(method = spec$method, figures = length(essential_figure_filenames())))
}

essential_validate_figures_only_inventory <- function(output_root, methods) {
  expected <- sort(essential_figure_filenames())
  for (method in methods) {
    method_dir <- file.path(output_root, "Figures", method)
    actual <- sort(list.files(method_dir, pattern = "[.]pdf$", full.names = FALSE))
    other <- list.files(method_dir, all.files = TRUE, no.. = TRUE, full.names = FALSE)
    other <- setdiff(other, actual)
    if (!identical(actual, expected)) {
      stop("Unexpected figure-only PDF inventory for method: ", method, call. = FALSE)
    }
    if (length(other) > 0L) {
      stop("Unexpected non-PDF figure-only output for method: ", method, call. = FALSE)
    }
  }
  calculation_expected <- sort(essential_tgi_calculation_figure_filenames())
  calculation_dir <- file.path(output_root, "Figures", "TGI_calculation")
  calculation_actual <- sort(list.files(
    calculation_dir,
    pattern = "[.]pdf$",
    full.names = FALSE
  ))
  calculation_other <- list.files(
    calculation_dir,
    all.files = TRUE,
    no.. = TRUE,
    full.names = FALSE
  )
  calculation_other <- setdiff(calculation_other, calculation_actual)
  if (!identical(calculation_actual, calculation_expected)) {
    stop("Unexpected figure-only TGI-calculation PDF inventory", call. = FALSE)
  }
  if (length(calculation_other) > 0L) {
    stop("Unexpected non-PDF TGI-calculation figure-only output", call. = FALSE)
  }
  invisible(list(
    figures = length(expected) * length(methods) + length(calculation_expected)
  ))
}

essential_write_readme <- function(output_root, tgi_spec = essential_tgi_spec()) {
  lines <- c(
    "# Essential pseudotime-TGI analysis",
    "",
    "This directory was generated from the CellCycle and NonCellCycle 04h cell-level input tables.",
    "",
    "When either cell-level input path is omitted, the standalone entrypoint reconstructs the missing table from `scvelo_cell_metrics.csv`, `all_ploidy.tsv`, `sample_info.xlsx`, and `dt_Gem_VT_20241223_v4.xlsx` under `Data/in-vivo`.",
    "",
    "The standalone workflow uses base R plus the `ggplot2`, `ggrepel`, `readr`, and `readxl` packages.",
    "",
    "## Selected TGI outcome",
    "",
    "The CLI default is Day-17 TGI using the mean of initial-ploidy-matched untreated controls. AUC and the other control summaries remain available through explicit CLI options.",
    "",
    paste0("- Measure: `", tgi_spec$measure, "`."),
    paste0("- Outcome: `", tgi_spec$outcome, "`."),
    paste0("- Matched untreated-control summary: `", tgi_spec$control_summary, "`."),
    if (identical(tgi_spec$outcome, "day")) paste0("- TGI day: `", tgi_spec$day, "`.") else NULL,
    "",
    "TGI is calculated as `100 * (1 - mouse tumor-growth delta / matched-control reference delta)`. Controls are matched by initial ploidy. `--control-summary` chooses whether the matched untreated-control deltas are summarized by their mean, median, or maximum; it does not summarize the treated mice.",
    "",
    "CLI options: `--tgi-outcome=auc|day`, `--control-summary=mean|median|max`, and (for the day outcome) `--tgi-day=<available day>`. With no TGI options, these resolve to `--tgi-outcome=day --tgi-day=17 --control-summary=mean`.",
    "",
    "The established output filenames retain `AUC` where present for backward compatibility. The selected outcome is recorded in the statistical `tgi_measure` fields, in scatter-plot metadata, and in plot titles and axes.",
    "",
    "## Contents",
    "",
    "- `Figures/<method>/`: eight CellCycle PDF figures per method.",
    "- `stats/<method>/`: fourteen statistical CSV files per method.",
    "- `plot_data/<method>/`: eight plotting-data CSV files per method, one for each figure.",
    "- `Figures/TGI_calculation/`: two shared PDFs explaining the selected TGI calculation.",
    "- `stats/TGI_calculation/`: the per-mouse calculation components and run metadata.",
    "- `plot_data/TGI_calculation/`: the exact data behind both shared calculation figures.",
    "",
    "## Methods",
    "",
    "- `initial_ploidy`: original 2N/4N sample grouping.",
    "- `ETP_fixed_threshold_2_25`: sample mean ETP threshold 2.25.",
    "- `ETP_boundary_stress_threshold_2_375`: sample mean ETP threshold 2.375.",
    "- `ETP_reference_balanced_threshold_2_24`: sample mean ETP threshold 2.24.",
    "",
    "## Figure data",
    "",
    "Each plotting-data CSV contains the exact rows used by its PDF. Scatter-plot tables also contain the Pearson and Spearman statistics, permutation P values, permutation mode, and the annotation text printed in the figure.",
    "",
    "`CellCycle_direct_group_ecdf_comparisons_selected_3panel.pdf` retains original grid positions `(1,1)`, `(3,2)`, and `(3,3)` (panels 1, 8, and 9). Its plotting-data CSV is a subset of the full direct-ECDF plotting table, and both figures reuse `CellCycle_direct_group_ecdf_comparisons_11panel_tests.csv`; no duplicate statistics file is produced.",
    "",
    "`CellCycle_TGI_group_boxplot.pdf` compares the selected CLI TGI measure between 2N and 4N treated tumors for the initial-ploidy method, or between ETP-lower and ETP-higher treated tumors for an ETP method. The boxes are adjacent independent groups, not one-to-one mouse pairs. The displayed primary P value comes from a group-label permutation stratified by dose, and the reported effect is adjusted for dose.",
    "",
    "`TGI_calculation_growth_trajectories.pdf` shows every mouse's baseline-adjusted tumor-growth trajectory, facets mice by initial ploidy, overlays the selected summary of matched untreated controls, and marks Day 17. `TGI_calculation_components_by_treated_mouse.pdf` shows the treated-mouse growth delta and its matched-control reference for every treated mouse, with TGI recalculated from the plotted values.",
    "",
    "`CellCycle_ecdf_rmse_vs_ploidy.pdf` uses treated mice only for the correlation and one threshold-independent reference formed by averaging the ECDF of each of the eight 0 mg/kg mice with equal sample weight. Initial ploidy, rather than the threshold-defined ETP group, supplies the point shape.",
    "",
    "Figures can be regenerated without the cell-level inputs and without rerunning statistical tests by using `--figures_only=TRUE --tables_root=<existing-output-root>`. In this mode the workflow reads only `plot_data/` and the required files in `stats/`, and writes only `Figures/`.",
    "",
    "The NonCellCycle input is used only when deriving sample mean end-timepoint ploidy. All figures and reported associations use CellCycle cells."
  )
  writeLines(lines, file.path(output_root, "README.md"))
  invisible(file.path(output_root, "README.md"))
}

essential_run_method <- function(
  cellcycle_path,
  noncellcycle_path,
  output_root,
  spec,
  seed,
  n_perm,
  n_boot,
  tgi_spec = essential_tgi_spec()
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required", call. = FALSE)
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("ggrepel is required", call. = FALSE)
  set.seed(seed)
  figures_dir <- essential_ensure_dir(file.path(output_root, "Figures", spec$method))
  stats_dir <- essential_ensure_dir(file.path(output_root, "stats", spec$method))
  plot_data_dir <- essential_ensure_dir(file.path(output_root, "plot_data", spec$method))

  cellcycle <- essential_read_input(cellcycle_path, "CellCycle")
  noncellcycle <- essential_read_input(noncellcycle_path, "NonCellCycle")
  selected_tgi <- essential_apply_tgi_spec(cellcycle, noncellcycle, tgi_spec)
  cellcycle <- selected_tgi$cellcycle
  noncellcycle <- selected_tgi$noncellcycle
  assignments <- essential_derive_assignments(cellcycle, noncellcycle)
  data <- essential_prepare_cellcycle(cellcycle, assignments, spec)
  sample_meta <- essential_build_sample_meta(data, tgi_spec$measure)
  equal_shift <- essential_calculate_shift_metrics(data, sample_meta, pooled_reference = FALSE)
  pooled_shift <- essential_calculate_shift_metrics(data, sample_meta, pooled_reference = TRUE)
  ploidy_shift <- essential_calculate_shift_metrics(
    data,
    sample_meta,
    pooled_reference = FALSE,
    control_scope = "all_untreated"
  )
  dose_tests <- essential_run_dose_tests(data, n_perm)

  association <- essential_run_tgi_associations(
    equal_shift,
    pooled_shift,
    spec$method,
    n_perm,
    primary_tgi_column = tgi_spec$measure
  )
  essential_write_csv(association, file.path(stats_dir, "CellCycle_TGI_associations_ecdf_rmse.csv"))
  primary <- association[association$pre_specified_primary, , drop = FALSE]
  primary_pearson <- data.frame(
    n = primary$n, estimate = primary$pearson_r, asymptotic_p = primary$pearson_p_asymptotic,
    permutation_p_two_sided = primary$pearson_p_permutation_two_sided,
    permutation_p_positive = primary$pearson_p_permutation_positive,
    n_permutations = primary$pearson_n_permutations,
    permutation_mode = primary$pearson_permutation_mode
  )
  primary_spearman <- data.frame(
    n = primary$n, estimate = primary$spearman_rho, asymptotic_p = primary$spearman_p_asymptotic,
    permutation_p_two_sided = primary$spearman_p_permutation_two_sided,
    permutation_p_positive = primary$spearman_p_permutation_positive,
    n_permutations = primary$spearman_n_permutations,
    permutation_mode = primary$spearman_permutation_mode
  )
  treated <- equal_shift[equal_shift$dose_mg > 0, , drop = FALSE]
  treated <- essential_add_tgi_plot_metadata(treated, tgi_spec)
  treated$analysis_group <- factor(treated$analysis_group, levels = spec$group_levels)
  treated$dose <- factor(treated$dose, levels = unique(treated$dose[order(treated$dose_mg)]))
  primary_plot_data <- essential_attach_association(
    treated,
    "CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref",
    primary_pearson,
    primary_spearman
  )
  essential_write_csv(
    primary_plot_data,
    file.path(plot_data_dir, "CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref_plot_data.csv")
  )
  primary_plot <- essential_scatter_plot(
    primary_plot_data,
    "ecdf_rmse",
    "tgi_value",
    "dose",
    "analysis_group",
    essential_dose_colors(),
    paste0(
      "Cell-cycle-associated tumor cells: ", tgi_spec$title_label,
      " vs sample-equal ECDF shift"
    ),
    "ECDF RMSE from group-matched equal-sample 0 mg/kg reference",
    tgi_spec$axis_label,
    "Dose",
    spec$group_label,
    primary_plot_data$plot_annotation[[1L]]
  )
  essential_save_pdf(
    primary_plot,
    file.path(figures_dir, "CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref.pdf"),
    6.8,
    6.8
  )

  set.seed(seed + 51L)
  group_comparison <- essential_tgi_group_comparison(treated, spec, tgi_spec, n_perm)
  essential_write_csv(
    group_comparison,
    file.path(stats_dir, "CellCycle_TGI_group_comparison.csv")
  )
  group_plot_data <- treated[
    order(treated$dose_mg, treated$analysis_group, treated$sample_id),
    ,
    drop = FALSE
  ]
  group_plot_data$analysis_id <- "CellCycle_TGI_group_boxplot"
  essential_write_csv(
    group_plot_data,
    file.path(plot_data_dir, "CellCycle_TGI_group_boxplot_plot_data.csv")
  )
  group_plot <- essential_tgi_group_boxplot(
    group_plot_data,
    group_comparison,
    spec,
    tgi_spec
  )
  essential_save_pdf(
    group_plot,
    file.path(figures_dir, "CellCycle_TGI_group_boxplot.pdf"),
    6.8,
    6.4
  )

  direct <- essential_direct_ecdf(data, dose_tests, spec, n_perm)
  essential_write_csv(
    direct$tests,
    file.path(stats_dir, "CellCycle_direct_group_ecdf_comparisons_11panel_tests.csv")
  )
  essential_write_csv(
    direct$plot_data,
    file.path(plot_data_dir, "CellCycle_direct_group_ecdf_comparisons_plot_data.csv")
  )
  essential_save_pdf(
    direct$plot,
    file.path(figures_dir, "CellCycle_direct_group_ecdf_comparisons.pdf"),
    15,
    15
  )
  selected_direct <- essential_select_direct_panels(direct$plot_data, direct$tests)
  essential_write_csv(
    selected_direct$curves,
    file.path(
      plot_data_dir,
      "CellCycle_direct_group_ecdf_comparisons_selected_3panel_plot_data.csv"
    )
  )
  selected_direct_plot <- essential_direct_plot_from_tables(
    selected_direct$curves,
    selected_direct$tests,
    spec,
    title = paste0(
      "CellCycle: selected direct group mean ECDF comparisons (", spec$method, ")"
    )
  )
  essential_save_pdf(
    selected_direct_plot,
    file.path(
      figures_dir,
      "CellCycle_direct_group_ecdf_comparisons_selected_3panel.pdf"
    ),
    15,
    5.5
  )

  set.seed(seed + 101L)
  mean_association <- essential_run_mean_etp_associations(
    equal_shift,
    spec$method,
    n_perm,
    primary_tgi_column = tgi_spec$measure
  )
  essential_write_csv(mean_association, file.path(stats_dir, "CellCycle_TGI_associations_mean_ETP.csv"))
  mean_primary <- mean_association[mean_association$tgi_measure == tgi_spec$measure, , drop = FALSE]
  mean_pearson <- data.frame(
    n = mean_primary$n, estimate = mean_primary$pearson_r, asymptotic_p = mean_primary$pearson_p_asymptotic,
    permutation_p_two_sided = mean_primary$pearson_p_permutation_two_sided,
    permutation_p_positive = mean_primary$pearson_p_permutation_positive,
    n_permutations = mean_primary$pearson_n_permutations,
    permutation_mode = mean_primary$pearson_permutation_mode
  )
  mean_spearman <- data.frame(
    n = mean_primary$n, estimate = mean_primary$spearman_rho, asymptotic_p = mean_primary$spearman_p_asymptotic,
    permutation_p_two_sided = mean_primary$spearman_p_permutation_two_sided,
    permutation_p_positive = mean_primary$spearman_p_permutation_positive,
    n_permutations = mean_primary$spearman_n_permutations,
    permutation_mode = mean_primary$spearman_permutation_mode
  )
  mean_plot_data <- essential_attach_association(
    treated,
    "CellCycle_TGI_AUC_vs_mean_ETP",
    mean_pearson,
    mean_spearman
  )
  essential_write_csv(
    mean_plot_data,
    file.path(plot_data_dir, "CellCycle_TGI_AUC_vs_mean_ETP_plot_data.csv")
  )
  mean_plot <- essential_scatter_plot(
    mean_plot_data,
    "sample_mean_endpoint_ploidy",
    "tgi_value",
    "dose",
    "analysis_group",
    essential_dose_colors(),
    paste0(
      "Cell-cycle-associated tumor cells: ", tgi_spec$title_label,
      " vs sample mean ETP"
    ),
    "Sample mean ETP",
    tgi_spec$axis_label,
    "Dose",
    spec$group_label,
    mean_plot_data$plot_annotation[[1L]]
  )
  essential_save_pdf(
    mean_plot,
    file.path(figures_dir, "CellCycle_TGI_AUC_vs_mean_ETP.pdf"),
    6.8,
    6.8
  )

  set.seed(seed + 201L)
  primary_loo <- essential_leave_one_out(treated, "ecdf_rmse", tgi_spec$measure)
  primary_boot <- essential_bootstrap(treated, "ecdf_rmse", n_boot, tgi_spec$measure)
  primary_summary <- essential_primary_robustness_summary(association, primary_loo, primary_boot)
  essential_write_csv(primary_loo, file.path(stats_dir, "CellCycle_primary_TGI_leave_one_out.csv"))
  essential_write_csv(primary_boot, file.path(stats_dir, "CellCycle_primary_TGI_bootstrap.csv"))
  essential_write_csv(primary_summary, file.path(stats_dir, "CellCycle_primary_TGI_robustness_summary.csv"))

  set.seed(seed + 301L)
  mean_loo <- essential_leave_one_out(
    treated, "sample_mean_endpoint_ploidy", tgi_spec$measure
  )
  mean_boot <- essential_bootstrap(
    treated, "sample_mean_endpoint_ploidy", n_boot, tgi_spec$measure
  )
  mean_summary <- essential_mean_etp_robustness_summary(
    spec$method,
    mean_association,
    mean_loo,
    mean_boot,
    n_boot,
    primary_tgi_column = tgi_spec$measure
  )
  essential_write_csv(mean_loo, file.path(stats_dir, "CellCycle_mean_ETP_TGI_leave_one_out.csv"))
  essential_write_csv(mean_boot, file.path(stats_dir, "CellCycle_mean_ETP_TGI_bootstrap.csv"))
  essential_write_csv(mean_summary, file.path(stats_dir, "CellCycle_mean_ETP_TGI_robustness_summary.csv"))

  model_results <- essential_fit_shift_models(treated, spec, tgi_spec$measure)
  essential_write_csv(
    model_results$models,
    file.path(stats_dir, "CellCycle_AUC_TGI_shift_ploidy_dose_models.csv")
  )
  essential_write_csv(
    model_results$summaries,
    file.path(stats_dir, "CellCycle_AUC_TGI_shift_ploidy_dose_model_summaries.csv")
  )
  model_plot_data <- primary_plot_data
  model_plot_data$analysis_id <- "CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose"
  essential_write_csv(
    model_plot_data,
    file.path(plot_data_dir, "CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose_plot_data.csv")
  )
  model_plot <- essential_scatter_plot(
    model_plot_data,
    "ecdf_rmse",
    "tgi_value",
    "analysis_group",
    "dose",
    essential_group_colors(spec),
    paste0(
      "CellCycle: ", tgi_spec$title_label,
      " association with pseudotime shift and ploidy"
    ),
    "ECDF RMSE from group-matched 0 mg/kg reference",
    tgi_spec$axis_label,
    spec$group_label,
    "Dose",
    model_plot_data$plot_annotation[[1L]]
  )
  essential_save_pdf(
    model_plot,
    file.path(figures_dir, "CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose.pdf"),
    6.8,
    6.8
  )

  set.seed(seed + 401L)
  confounding <- essential_confounding(
    equal_shift,
    spec,
    n_perm,
    tgi_spec$measure,
    ploidy_shift = ploidy_shift
  )
  essential_write_csv(confounding$ploidy_tests, file.path(stats_dir, "ploidy_confounding_tests.csv"))
  essential_write_csv(confounding$residualized, file.path(stats_dir, "residualized_TGI_associations.csv"))

  centered_result <- confounding$residualized[confounding$residualized$analysis == "within_dose_centered", , drop = FALSE]
  centered_pearson <- data.frame(
    n = centered_result$n,
    estimate = centered_result$estimate,
    asymptotic_p = centered_result$p_value,
    permutation_p_two_sided = centered_result$p_permutation_two_sided,
    permutation_p_positive = NA_real_,
    n_permutations = centered_result$n_permutations,
    permutation_mode = centered_result$permutation_mode
  )
  centered_spearman <- essential_permutation_cor(
    confounding$centered$shift_centered,
    confounding$centered$tgi_centered,
    "spearman",
    n_perm,
    exact = TRUE,
    strata = confounding$centered$dose_mg
  )
  centered_plot_data <- essential_attach_association(
    confounding$centered,
    "CellCycle_TGI_association_within_dose_centered",
    centered_pearson,
    centered_spearman
  )
  centered_plot_data <- essential_add_tgi_plot_metadata(
    centered_plot_data, tgi_spec, value_column = "tgi_centered"
  )
  centered_plot_data$dose <- factor(
    centered_plot_data$dose,
    levels = unique(centered_plot_data$dose[order(centered_plot_data$dose_mg)])
  )
  centered_plot_data$analysis_group <- factor(
    centered_plot_data$analysis_group,
    levels = spec$group_levels
  )
  essential_write_csv(
    centered_plot_data,
    file.path(plot_data_dir, "CellCycle_TGI_association_within_dose_centered_plot_data.csv")
  )
  centered_plot <- essential_scatter_plot(
    centered_plot_data,
    "shift_centered",
    "tgi_centered",
    "dose",
    "analysis_group",
    essential_dose_colors(),
    "CellCycle TGI association after within-dose centering",
    "Dose-centered ECDF RMSE",
    paste0("Dose-centered ", tgi_spec$axis_label),
    "Dose",
    spec$group_label,
    centered_plot_data$plot_annotation[[1L]],
    label_nudge_y = 1.0
  ) +
    ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35)
  essential_save_pdf(
    centered_plot,
    file.path(figures_dir, "CellCycle_TGI_association_within_dose_centered.pdf"),
    6.6,
    6.6
  )

  ploidy_primary <- confounding$ploidy_tests[
    confounding$ploidy_tests$sample_set == "treated" &
      confounding$ploidy_tests$shift_metric == "ecdf_rmse" &
      confounding$ploidy_tests$ploidy_measure == "mean_cell_ploidy",
    ,
    drop = FALSE
  ]
  pearson_row <- ploidy_primary[ploidy_primary$method == "pearson", , drop = FALSE]
  spearman_row <- ploidy_primary[ploidy_primary$method == "spearman", , drop = FALSE]
  ploidy_pearson <- data.frame(
    n = pearson_row$n,
    estimate = pearson_row$estimate,
    asymptotic_p = pearson_row$p_value,
    permutation_p_two_sided = pearson_row$p_permutation_two_sided,
    permutation_p_positive = NA_real_,
    n_permutations = pearson_row$n_permutations,
    permutation_mode = pearson_row$permutation_mode
  )
  ploidy_spearman <- data.frame(
    n = spearman_row$n,
    estimate = spearman_row$estimate,
    asymptotic_p = spearman_row$p_value,
    permutation_p_two_sided = spearman_row$p_permutation_two_sided,
    permutation_p_positive = NA_real_,
    n_permutations = spearman_row$n_permutations,
    permutation_mode = spearman_row$permutation_mode
  )
  ploidy_source <- ploidy_shift[
    ploidy_shift$dose_mg > 0 &
      is.finite(ploidy_shift$mean_cell_ploidy) &
      is.finite(ploidy_shift$ecdf_rmse),
    ,
    drop = FALSE
  ]
  ploidy_source <- ploidy_source[
    order(ploidy_source$dose_mg, ploidy_source$initial_ploidy, ploidy_source$sample_id),
    ,
    drop = FALSE
  ]
  ploidy_source$dose <- factor(
    ploidy_source$dose,
    levels = unique(ploidy_source$dose[order(ploidy_source$dose_mg)])
  )
  ploidy_source$analysis_group <- factor(ploidy_source$analysis_group, levels = spec$group_levels)
  ploidy_source$initial_ploidy <- factor(
    as.character(ploidy_source$initial_ploidy),
    levels = c("2N", "4N")
  )
  ploidy_plot_data <- essential_attach_association(
    ploidy_source,
    "CellCycle_ecdf_rmse_vs_ploidy",
    ploidy_pearson,
    ploidy_spearman
  )
  essential_write_csv(
    ploidy_plot_data,
    file.path(plot_data_dir, "CellCycle_ecdf_rmse_vs_ploidy_plot_data.csv")
  )
  ploidy_plot <- essential_scatter_plot(
    ploidy_plot_data,
    "mean_cell_ploidy",
    "ecdf_rmse",
    "dose",
    "initial_ploidy",
    essential_dose_colors(),
    "Treated CellCycle ECDF RMSE vs mean cell ploidy",
    "Mean cell ploidy",
    "ECDF RMSE",
    "Dose",
    "Initial ploidy",
    ploidy_plot_data$plot_annotation[[1L]],
    label_nudge_y = 0.004
  )
  essential_save_pdf(
    ploidy_plot,
    file.path(figures_dir, "CellCycle_ecdf_rmse_vs_ploidy.pdf"),
    6.4,
    6.4
  )

  invisible(list(method = spec$method, figures = 8L, stats = 14L, plot_data = 8L))
}

essential_validate_inventory <- function(output_root, methods) {
  os_metadata <- list.files(
    output_root,
    pattern = "^[.]DS_Store$",
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE
  )
  if (length(os_metadata) > 0L) unlink(os_metadata, force = TRUE)
  expected_figures <- 8L * length(methods) + length(essential_tgi_calculation_figure_filenames())
  expected_stats <- 14L * length(methods) + length(essential_tgi_calculation_stats_filenames())
  expected_plot_data <- 8L * length(methods) + length(essential_tgi_calculation_plot_data_filenames())
  figures <- list.files(file.path(output_root, "Figures"), pattern = "[.]pdf$", recursive = TRUE, full.names = TRUE)
  stats <- list.files(file.path(output_root, "stats"), pattern = "[.]csv$", recursive = TRUE, full.names = TRUE)
  plot_data <- list.files(file.path(output_root, "plot_data"), pattern = "[.]csv$", recursive = TRUE, full.names = TRUE)
  if (length(figures) != expected_figures) {
    stop("Unexpected PDF count: ", length(figures), "; expected ", expected_figures, call. = FALSE)
  }
  if (length(stats) != expected_stats) {
    stop("Unexpected statistics CSV count: ", length(stats), "; expected ", expected_stats, call. = FALSE)
  }
  if (length(plot_data) != expected_plot_data) {
    stop("Unexpected plot-data CSV count: ", length(plot_data), "; expected ", expected_plot_data, call. = FALSE)
  }
  unexpected <- list.files(
    output_root,
    pattern = "[.](png|ppt|pptx)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(unexpected) > 0L) {
    stop("Unexpected output files: ", paste(unexpected, collapse = ", "), call. = FALSE)
  }
  all_files <- list.files(
    output_root,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = FALSE
  )
  allowed <- grepl("[.](pdf|csv)$", all_files, ignore.case = TRUE) |
    basename(all_files) == "README.md"
  if (any(!allowed)) {
    stop("Unexpected essential output: ", paste(all_files[!allowed], collapse = ", "), call. = FALSE)
  }
  invisible(list(figures = length(figures), stats = length(stats), plot_data = length(plot_data)))
}
