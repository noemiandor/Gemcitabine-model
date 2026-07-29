#!/usr/bin/env Rscript

parse_cli_args <- function(args) {
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

arg_value <- function(args, name, default = NULL) {
  value <- args[[name]]
  if (!is.null(value) && length(value) == 1L && nzchar(value)) value else default
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript Code/in-vivo/figure7/generate_pseudotime_distribution_with_ploidy_dose_tgi.R \\",
      "    --input_root Data/in-vivo",
      "",
      "Default inputs under --input_root:",
      "  scvelo_cell_metrics.csv",
      "  all_ploidy.tsv",
      "  sample_info.xlsx",
      "  dt_Gem_VT_20241223_v4.xlsx",
      "",
      "Default outputs:",
      "  Data/in-vivo/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv",
      "  Data/in-vivo/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv",
      "",
      "Common options:",
      "  --scvelo_metrics_input PATH",
      "  --cell_ploidy_input PATH",
      "  --sample_info_input PATH",
      "  --growth_curve_input PATH",
      "  --cellcycle_output PATH",
      "  --noncellcycle_output PATH",
      "  --cellcycle_clusters 4c,6,10",
      "  --tgi_day 17",
      sep = "\n"
    ),
    "\n"
  )
}

resolve_script_path <- function() {
  file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_argument) > 0L) {
    return(normalizePath(sub("^--file=", "", file_argument[[1L]]), mustWork = TRUE))
  }
  NA_character_
}

is_absolute_path <- function(path) {
  grepl("^/", path)
}

resolve_path <- function(path, repo_root) {
  if (is_absolute_path(path)) {
    normalizePath(path, mustWork = FALSE)
  } else {
    normalizePath(file.path(repo_root, path), mustWork = FALSE)
  }
}

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing required R package: ", pkg, call. = FALSE)
  }
}

split_values <- function(value) {
  value <- trimws(as.character(value))
  if (!nzchar(value)) return(character())
  out <- unlist(strsplit(value, "[,;]", perl = TRUE), use.names = FALSE)
  out <- trimws(out)
  out[!is.na(out) & nzchar(out)]
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

clean_character <- function(x) {
  out <- trimws(as.character(x))
  out[out %in% c("", "NA", "NaN", "NULL", "None")] <- NA_character_
  out
}

first_existing_column <- function(data, candidates, label) {
  matches <- candidates[candidates %in% names(data)]
  if (length(matches) == 0L) {
    stop("Cannot find ", label, ". Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  matches[[1L]]
}

resolve_column_case_insensitive <- function(data, candidates) {
  lower_names <- tolower(names(data))
  for (candidate in candidates) {
    index <- match(tolower(candidate), lower_names)
    if (!is.na(index)) return(names(data)[[index]])
  }
  NA_character_
}

standardize_ploidy <- function(x) {
  out <- toupper(clean_character(x))
  out <- gsub("\\s+", "", out)
  out[out %in% c("2", "2.0", "2N")] <- "2N"
  out[out %in% c("4", "4.0", "4N")] <- "4N"
  out[!(out %in% c("2N", "4N"))] <- NA_character_
  out
}

standardize_dose <- function(x) {
  out <- clean_character(x)
  out[out %in% c("0", "0mg", "0 mg/kg", "0mg/kg")] <- "0mg/kg"
  out[out %in% c("30", "30mg", "30 mg/kg", "30mg/kg")] <- "30mg/kg"
  out[out %in% c("120", "120mg", "120 mg/kg", "120mg/kg")] <- "120mg/kg"
  out
}

dose_panel <- function(x) {
  out <- sub("mg/kg$", "", standardize_dose(x))
  out[!(out %in% c("0", "30", "120"))] <- NA_character_
  out
}

extract_barcode <- function(cell, sample) {
  cell <- clean_character(cell)
  sample <- clean_character(sample)
  barcode <- cell
  has_sample <- !is.na(cell) & !is.na(sample) & nzchar(sample) & startsWith(cell, paste0(sample, "_"))
  barcode[has_sample] <- substring(cell[has_sample], nchar(sample[has_sample]) + 2L)
  fallback <- !has_sample & !is.na(cell) & grepl("_", cell)
  barcode[fallback] <- sub("^[^_]+_", "", cell[fallback])
  barcode
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

read_sample_info_map <- function(path) {
  if (!file.exists(path)) stop("Missing sample information workbook: ", path, call. = FALSE)
  raw <- as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE)
  sample_column <- first_existing_column(raw, c("IDs", "sample", "sample_id", "sampleID"), "sample ID column")
  harvest_column <- first_existing_column(raw, c("harvest"), "harvest column")
  dose_matches <- c("Dose", "dose")
  dose_matches <- dose_matches[dose_matches %in% names(raw)]
  dose <- if (length(dose_matches) > 0L) raw[[dose_matches[[1L]]]] else NA_character_
  out <- data.frame(
    sample = clean_character(raw[[sample_column]]),
    harvest = clean_character(raw[[harvest_column]]),
    sample_info_dose = standardize_dose(dose),
    stringsAsFactors = FALSE
  )
  keep <- !is.na(out$sample) & nzchar(out$sample) & !is.na(out$harvest) & nzchar(out$harvest)
  out <- out[keep & !duplicated(out$sample), , drop = FALSE]
  rownames(out) <- NULL
  out
}

read_cell_ploidy_map <- function(path) {
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
    harvest = sub("\\.sps\\.cbs$", "", clean_character(raw$file)),
    barcode_raw = clean_character(raw$cell_id),
    cell_ploidy = safe_numeric(raw$ploidy),
    ploidy_frac_covered = if ("frac_covered" %in% names(raw)) {
      safe_numeric(raw$frac_covered)
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

build_cell_metrics <- function(scvelo_path, cell_ploidy_path, sample_info_path) {
  required_files <- c(scvelo_path, cell_ploidy_path, sample_info_path)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    stop("Missing input file(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
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
  pseudotime_column <- first_existing_column(raw, c("velocity_pseudotime", "pseudotime"), "pseudotime column")
  required <- c("cell", "TN", "clusters", "sample", "Ploidy", "Dose", pseudotime_column)
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0L) {
    stop("scVelo metrics are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  metrics <- data.frame(
    cell = clean_character(raw$cell),
    TN = clean_character(raw$TN),
    clusters = clean_character(raw$clusters),
    sample = clean_character(raw$sample),
    Ploidy = standardize_ploidy(raw$Ploidy),
    Dose = standardize_dose(raw$Dose),
    pseudotime = safe_numeric(raw[[pseudotime_column]]),
    stringsAsFactors = FALSE
  )
  metrics$barcode_raw <- extract_barcode(metrics$cell, metrics$sample)

  sample_map <- read_sample_info_map(sample_info_path)
  sample_index <- match(metrics$sample, sample_map$sample)
  metrics$harvest <- sample_map$harvest[sample_index]
  sample_dose <- sample_map$sample_info_dose[sample_index]
  missing_dose <- is.na(metrics$Dose) | !nzchar(metrics$Dose)
  metrics$Dose[missing_dose] <- sample_dose[missing_dose]

  ploidy_map <- read_cell_ploidy_map(cell_ploidy_path)
  metric_key <- paste(metrics$harvest, metrics$barcode_raw, sep = "\r")
  ploidy_key <- paste(ploidy_map$harvest, ploidy_map$barcode_raw, sep = "\r")
  ploidy_index <- match(metric_key, ploidy_key)
  metrics$cell_ploidy <- ploidy_map$cell_ploidy[ploidy_index]
  metrics
}

normalize_sample_id <- function(x) {
  out <- clean_character(x)
  out <- sub("-Count-HM$", "", out)
  out <- sub("-Count$", "", out)
  out <- sub("-HM$", "", out)
  out
}

infer_initial_ploidy <- function(sample, harvest = NULL) {
  sample <- clean_character(sample)
  harvest <- if (is.null(harvest)) rep(NA_character_, length(sample)) else clean_character(harvest)
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

infer_dose_from_harvest <- function(harvest) {
  harvest <- clean_character(harvest)
  dose <- rep(NA_character_, length(harvest))
  matched <- grepl("^SUM159-(2N|4N)-[0-9]+-", harvest, ignore.case = TRUE)
  dose[matched] <- sub(
    "^SUM159-(2N|4N)-([0-9]+)-.*$",
    "\\2",
    harvest[matched],
    ignore.case = TRUE
  )
  standardize_dose(dose)
}

calculate_tgi <- function(
  growth_curve_path,
  tgi_day = 17L,
  control_dose = "0mg/kg",
  baseline_day = "Day_0"
) {
  if (!file.exists(growth_curve_path)) stop("Missing growth curve workbook: ", growth_curve_path, call. = FALSE)
  growth <- as.data.frame(readxl::read_excel(growth_curve_path, sheet = 1), stringsAsFactors = FALSE)
  harvest_column <- resolve_column_case_insensitive(growth, c("harvest"))
  sample_column <- resolve_column_case_insensitive(
    growth,
    c("Sequencing IDs", "Sequencing.IDs", "Sequencing ID", "SequencingIDs")
  )
  if (is.na(harvest_column)) stop("Growth curve workbook is missing a harvest column", call. = FALSE)
  if (is.na(sample_column)) stop("Growth curve workbook is missing a Sequencing IDs column", call. = FALSE)

  day_columns <- grep("^Day_[0-9]+$", names(growth), value = TRUE)
  day_numbers <- safe_numeric(sub("^Day_", "", day_columns))
  day_columns <- day_columns[order(day_numbers)]
  day_numbers <- safe_numeric(sub("^Day_", "", day_columns))
  if (!(baseline_day %in% day_columns)) stop("Missing baseline day: ", baseline_day, call. = FALSE)
  tgi_day_column <- paste0("Day_", as.integer(tgi_day))
  if (!(tgi_day_column %in% day_columns)) {
    stop("Missing configured TGI day: ", tgi_day_column, call. = FALSE)
  }
  baseline_time <- safe_numeric(sub("^Day_", "", baseline_day))
  endpoint_columns <- day_columns[day_numbers > baseline_time]

  rows <- lapply(seq_len(nrow(growth)), function(index) {
    values <- safe_numeric(unlist(growth[index, day_columns, drop = FALSE], use.names = FALSE))
    names(values) <- day_columns
    baseline_value <- values[[baseline_day]]
    row <- data.frame(
      growth_curve_row = index,
      sample_id_raw = as.character(growth[[sample_column]][index]),
      sample_id = normalize_sample_id(growth[[sample_column]][index]),
      harvest = as.character(growth[[harvest_column]][index]),
      initial_ploidy = infer_initial_ploidy(
        growth[[sample_column]][index],
        growth[[harvest_column]][index]
      ),
      dose = infer_dose_from_harvest(growth[[harvest_column]][index]),
      baseline_day = baseline_day,
      baseline_volume = baseline_value,
      stringsAsFactors = FALSE
    )
    for (day_column in endpoint_columns) {
      row[[paste0("tumor_volume_", day_column)]] <- values[[day_column]]
      if (identical(day_column, tgi_day_column)) {
        row[[paste0("tumor_volume_delta_", day_column)]] <-
          values[[day_column]] - baseline_value
      }
    }
    row
  })
  output <- do.call(rbind, rows)
  keep <- !is.na(output$sample_id) & nzchar(output$sample_id) &
    is.finite(output$baseline_volume)
  output <- output[keep, , drop = FALSE]
  output$dose <- standardize_dose(output$dose)
  output$control_match_group <- output$initial_ploidy
  delta_column <- paste0("tumor_volume_delta_", tgi_day_column)
  mean_column <- paste0("matched_control_mean_delta_", tgi_day_column)
  count_column <- paste0("matched_control_n_", tgi_day_column)
  tgi_column <- paste0("TGI_percent_", tgi_day_column)
  output[[mean_column]] <- NA_real_
  output[[count_column]] <- NA_integer_

  control_dose <- standardize_dose(control_dose)[[1L]]
  groups <- unique(output$control_match_group[!is.na(output$control_match_group)])
  for (group in groups) {
    target <- output$control_match_group == group
    control_day <- target & output$dose == control_dose &
      is.finite(output[[delta_column]])
    output[[mean_column]][target] <- mean(
      output[[delta_column]][control_day],
      na.rm = TRUE
    )
    output[[count_column]][target] <- sum(control_day)
  }

  output[[tgi_column]] <- 100 * (
    1 - output[[delta_column]] / output[[mean_column]]
  )
  output[[tgi_column]][!is.finite(output[[tgi_column]])] <- NA_real_
  output[order(output$initial_ploidy, output$dose, output$sample_id), , drop = FALSE]
}

prepare_tgi_join <- function(tgi, tgi_day) {
  day_label <- paste0("Day_", as.integer(tgi_day))
  base <- data.frame(
    sampleID = tgi$sample_id,
    tumor_volume_baseline_day = tgi$baseline_day,
    tumor_volume_baseline = tgi$baseline_volume,
    stringsAsFactors = FALSE
  )
  extra_columns <- grep(
    paste0(
      "^(tumor_volume_Day_[0-9]+|",
      "tumor_volume_delta_", day_label, "|",
      "TGI_percent_", day_label, ")$"
    ),
    names(tgi),
    value = TRUE
  )
  cbind(base, tgi[, extra_columns, drop = FALSE])
}

build_compartment_table <- function(
  tumor,
  target_clusters,
  tgi,
  label,
  tgi_day
) {
  target <- tumor[tumor$clusters %in% target_clusters, , drop = FALSE]
  target <- target[
    !is.na(target$sampleID) & nzchar(target$sampleID) & is.finite(target$pseudotime),
    ,
    drop = FALSE
  ]
  if (nrow(target) == 0L) {
    stop("No ", label, " Tumor cells were found for clusters: ", paste(target_clusters, collapse = ", "), call. = FALSE)
  }
  tgi_join <- prepare_tgi_join(tgi, tgi_day)
  tgi_index <- match(target$sampleID, tgi_join$sampleID)
  if (anyNA(tgi_index)) {
    missing_samples <- unique(target$sampleID[is.na(tgi_index)])
    warning("TGI was not matched for ", label, " samples: ", paste(missing_samples, collapse = ", "), call. = FALSE)
  }

  output <- data.frame(
    cell_id = target$cell,
    sample_id = target$sampleID,
    cluster = target$clusters,
    initial_ploidy = target$Ploidy,
    gemcitabine_dose = target$Dose,
    gemcitabine_dose_mg_per_kg = safe_numeric(target$DosePanel),
    pseudotime = target$pseudotime,
    cell_ploidy = target$cell_ploidy,
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

generate_outputs <- function(
  scvelo_path,
  cell_ploidy_path,
  sample_info_path,
  growth_curve_path,
  cellcycle_output,
  noncellcycle_output,
  tgi_day = 17L,
  cellcycle_clusters = c("4c", "6", "10")
) {
  metrics <- build_cell_metrics(scvelo_path, cell_ploidy_path, sample_info_path)
  if (anyNA(metrics$cell) ||
      any(!nzchar(metrics$cell)) ||
      anyDuplicated(metrics$cell) ||
      any(!is.finite(metrics$pseudotime)) ||
      any(metrics$pseudotime < 0 | metrics$pseudotime > 1)) {
    stop(
      "scVelo metrics require unique nonempty cells and one finite ",
      "pseudotime within [0,1] per cell",
      call. = FALSE
    )
  }
  metrics$DosePanel <- dose_panel(metrics$Dose)
  metrics$sampleID <- metrics$sample
  tumor <- metrics[
    !is.na(metrics$TN) & metrics$TN == "Tumor" &
      metrics$Ploidy %in% c("2N", "4N") &
      metrics$DosePanel %in% c("0", "30", "120"),
    ,
    drop = FALSE
  ]
  if (nrow(tumor) == 0L) stop("No eligible Tumor cells were found in scVelo metrics", call. = FALSE)
  tgi <- calculate_tgi(growth_curve_path, tgi_day = tgi_day)

  tumor_clusters <- unique(as.character(tumor$clusters))
  tumor_clusters <- tumor_clusters[!is.na(tumor_clusters) & nzchar(tumor_clusters)]
  cellcycle_clusters <- cellcycle_clusters[cellcycle_clusters %in% tumor_clusters]
  if (length(cellcycle_clusters) == 0L) {
    stop("No CellCycle clusters are present in eligible Tumor cells.", call. = FALSE)
  }
  noncellcycle_clusters <- setdiff(tumor_clusters, cellcycle_clusters)
  if (length(noncellcycle_clusters) == 0L) {
    stop("No NonCellCycle clusters remain after excluding CellCycle clusters.", call. = FALSE)
  }

  cellcycle <- build_compartment_table(
    tumor,
    cellcycle_clusters,
    tgi,
    "CellCycle",
    tgi_day
  )
  noncellcycle <- build_compartment_table(
    tumor,
    noncellcycle_clusters,
    tgi,
    "NonCellCycle",
    tgi_day
  )

  ensure_dir(dirname(cellcycle_output))
  ensure_dir(dirname(noncellcycle_output))
  readr::write_csv(cellcycle, cellcycle_output, na = "NA")
  readr::write_csv(noncellcycle, noncellcycle_output, na = "NA")

  list(
    cellcycle = cellcycle,
    noncellcycle = noncellcycle,
    tumor = tumor,
    cellcycle_clusters = cellcycle_clusters,
    noncellcycle_clusters = noncellcycle_clusters
  )
}

main <- function() {
  args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  if (!is.null(args$help) || !is.null(args$h)) {
    usage()
    quit(save = "no", status = 0L)
  }

  script_path <- resolve_script_path()
  script_dir <- if (!is.na(script_path)) dirname(script_path) else getwd()
  repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE)

  input_root <- resolve_path(arg_value(args, "input_root", file.path("Data", "in-vivo")), repo_root)
  scvelo_metrics_input <- resolve_path(
    arg_value(args, "scvelo_metrics_input", file.path(input_root, "scvelo_cell_metrics.csv")),
    repo_root
  )
  cell_ploidy_input <- resolve_path(
    arg_value(args, "cell_ploidy_input", file.path(input_root, "all_ploidy.tsv")),
    repo_root
  )
  sample_info_input <- resolve_path(
    arg_value(args, "sample_info_input", file.path(input_root, "sample_info.xlsx")),
    repo_root
  )
  growth_curve_input <- resolve_path(
    arg_value(args, "growth_curve_input", file.path(input_root, "dt_Gem_VT_20241223_v4.xlsx")),
    repo_root
  )
  cellcycle_output <- resolve_path(
    arg_value(
      args,
      "cellcycle_output",
      file.path("Data", "in-vivo", "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv")
    ),
    repo_root
  )
  noncellcycle_output <- resolve_path(
    arg_value(
      args,
      "noncellcycle_output",
      file.path("Data", "in-vivo", "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv")
    ),
    repo_root
  )
  cellcycle_clusters <- split_values(arg_value(args, "cellcycle_clusters", "4c,6,10"))
  tgi_day <- suppressWarnings(as.integer(arg_value(args, "tgi_day", "17")))

  if (length(cellcycle_clusters) == 0L) {
    stop("--cellcycle_clusters must contain at least one cluster", call. = FALSE)
  }
  if (length(tgi_day) != 1L || !is.finite(tgi_day) || tgi_day < 0L) {
    stop("--tgi_day must be one non-negative integer", call. = FALSE)
  }
  require_package("readr")
  require_package("readxl")

  message("Generating pseudotime cell-level TGI inputs")
  message("  input_root: ", input_root)
  message("  scVelo metrics: ", scvelo_metrics_input)
  message("  cell ploidy: ", cell_ploidy_input)
  message("  sample info: ", sample_info_input)
  message("  growth curve: ", growth_curve_input)
  message("  CellCycle output: ", cellcycle_output)
  message("  NonCellCycle output: ", noncellcycle_output)
  message("  CellCycle clusters: ", paste(cellcycle_clusters, collapse = ", "))
  message("  TGI day: ", tgi_day)

  result <- generate_outputs(
    scvelo_path = scvelo_metrics_input,
    cell_ploidy_path = cell_ploidy_input,
    sample_info_path = sample_info_input,
    growth_curve_path = growth_curve_input,
    cellcycle_output = cellcycle_output,
    noncellcycle_output = noncellcycle_output,
    tgi_day = tgi_day,
    cellcycle_clusters = cellcycle_clusters
  )

  message("Completed")
  message("  CellCycle rows: ", nrow(result$cellcycle))
  message("  NonCellCycle rows: ", nrow(result$noncellcycle))
  message("  eligible Tumor rows: ", nrow(result$tumor))
  message("  NonCellCycle clusters: ", paste(result$noncellcycle_clusters, collapse = ", "))
}

if (identical(environment(), globalenv())) {
  main()
}
