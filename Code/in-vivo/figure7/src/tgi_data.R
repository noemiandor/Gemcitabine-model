# Input and TGI preparation adapted from the reviewed 04h workflow at 77caec93.

figure7_numeric <- function(x) suppressWarnings(as.numeric(x))

figure7_unique_sample_value <- function(data, column, sample_id, numeric = FALSE) {
  values <- data[data$sample_id == sample_id, column]
  if (numeric) {
    values <- unique(figure7_numeric(values))
    values <- values[is.finite(values)]
  } else {
    values <- unique(trimws(as.character(values)))
    values <- values[!is.na(values) & nzchar(values)]
  }
  if (length(values) != 1L) {
    figure7_stop("Expected one sample-level value for ", column, " in ", sample_id,
                 "; found ", length(values))
  }
  values[[1L]]
}

figure7_read_cell_table <- function(path, compartment, config) {
  delta_measure <- figure7_tgi_delta_measure(config)
  tgi_measure <- figure7_tgi_measure(config)
  required <- c(
    "cell_id", "sample_id", "initial_ploidy", "gemcitabine_dose",
    "gemcitabine_dose_mg_per_kg", "pseudotime", "cell_ploidy",
    "tumor_volume_baseline_day", "tumor_volume_baseline",
    delta_measure, tgi_measure
  )
  if (!file.exists(path)) figure7_stop("Missing ", compartment, " input: ", path)
  data <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing <- setdiff(required, names(data))
  if (length(missing)) figure7_stop(compartment, " input is missing: ", paste(missing, collapse = ", "))
  numeric_columns <- unique(c(
    "gemcitabine_dose_mg_per_kg", "pseudotime", "cell_ploidy",
    "tumor_volume_baseline", delta_measure, tgi_measure,
    grep("^tumor_volume_Day_[0-9]+$", names(data), value = TRUE)
  ))
  for (column in intersect(numeric_columns, names(data))) data[[column]] <- figure7_numeric(data[[column]])
  if (any(!data$initial_ploidy %in% c("2N", "4N"))) figure7_stop("Unexpected initial-ploidy label in ", compartment)
  data$compartment <- compartment
  data
}

figure7_sample_table <- function(cellcycle, noncellcycle, config) {
  delta_measure <- figure7_tgi_delta_measure(config)
  tgi_measure <- figure7_tgi_measure(config)
  embedded_tgi_measure <- paste0("embedded_", tgi_measure)
  columns <- c("cell_id", "sample_id", "initial_ploidy", "gemcitabine_dose",
               "gemcitabine_dose_mg_per_kg", "cell_ploidy")
  union <- rbind(cellcycle[, columns], noncellcycle[, columns])
  union <- union[is.finite(union$cell_ploidy), , drop = FALSE]
  duplicate_ids <- unique(union$cell_id[duplicated(union$cell_id)])
  for (cell_id in duplicate_ids) {
    local <- union[union$cell_id == cell_id, , drop = FALSE]
    if (length(unique(local$sample_id)) != 1L || length(unique(local$cell_ploidy)) != 1L) {
      figure7_stop("Conflicting duplicated cell_id: ", cell_id)
    }
  }
  union <- union[!duplicated(union$cell_id), , drop = FALSE]
  rows <- lapply(sort(unique(union$sample_id)), function(id) {
    local <- union[union$sample_id == id, , drop = FALSE]
    cc <- cellcycle[cellcycle$sample_id == id, , drop = FALSE]
    if (!nrow(cc)) figure7_stop("Sample has no CellCycle cells: ", id)
    row <- data.frame(
      sample_id = id,
      initial_ploidy = figure7_unique_sample_value(local, "initial_ploidy", id),
      dose = figure7_unique_sample_value(local, "gemcitabine_dose", id),
      dose_mg = figure7_unique_sample_value(local, "gemcitabine_dose_mg_per_kg", id, TRUE),
      sample_mean_endpoint_ploidy = mean(local$cell_ploidy),
      sample_median_endpoint_ploidy = stats::median(local$cell_ploidy),
      n_endpoint_ploidy_cells = nrow(local), n_cellcycle_cells = nrow(cc),
      mean_cellcycle_pseudotime = mean(cc$pseudotime, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    row[[delta_measure]] <- figure7_unique_sample_value(cc, delta_measure, id, TRUE)
    row[[embedded_tgi_measure]] <- figure7_unique_sample_value(cc, tgi_measure, id, TRUE)
    row
  })
  samples <- do.call(rbind, rows)
  samples$matched_control_reference_delta <- NA_real_
  samples$matched_control_n <- NA_integer_
  samples$matched_control_sample_ids <- NA_character_
  for (ploidy in c("2N", "4N")) {
    controls <- samples[samples$initial_ploidy == ploidy & samples$dose_mg == 0, , drop = FALSE]
    if (!nrow(controls)) figure7_stop("No untreated controls for initial ploidy ", ploidy)
    hit <- samples$initial_ploidy == ploidy
    samples$matched_control_reference_delta[hit] <- mean(controls[[delta_measure]])
    samples$matched_control_n[hit] <- nrow(controls)
    samples$matched_control_sample_ids[hit] <- paste(sort(controls$sample_id), collapse = ";")
  }
  samples[[tgi_measure]] <- 100 * (1 - samples[[delta_measure]] /
                                      samples$matched_control_reference_delta)
  treated <- samples$dose_mg %in% as.numeric(unlist(config$tgi$treated_doses_mg_per_kg))
  difference <- abs(samples[[tgi_measure]] - samples[[embedded_tgi_measure]])
  tolerance <- as.numeric(config$tgi$numerical_tolerance)
  if (any(!is.finite(samples[[embedded_tgi_measure]][treated]))) {
    figure7_stop("All treated mice require a finite embedded Day-", figure7_tgi_day(config),
                 " TGI for verification")
  }
  if (any(difference[treated] > tolerance)) {
    figure7_stop("Recomputed Day-", figure7_tgi_day(config), " TGI disagrees with embedded values for: ",
                 paste(samples$sample_id[treated & difference > tolerance], collapse = ", "))
  }
  samples$etp_group <- ifelse(samples$sample_mean_endpoint_ploidy > as.numeric(config$etp$threshold),
                              "ETP-higher", "ETP-lower")
  samples <- figure7_add_tgi_metadata(samples, config)
  samples <- samples[order(samples$dose_mg, samples$initial_ploidy, samples$sample_id), , drop = FALSE]
  treated_doses <- as.numeric(unlist(config$tgi$treated_doses_mg_per_kg))
  treated_rows <- samples$dose_mg %in% treated_doses
  if (sum(treated_rows) != 8L || any(samples$dose_mg > 0 & !treated_rows)) {
    figure7_stop("Panels 7C-7E require exactly eight treated mice at 30 or 120 mg/kg")
  }
  if (any(!is.finite(samples[[tgi_measure]][treated_rows]))) {
    figure7_stop("All eight treated mice require finite recomputed Day-", figure7_tgi_day(config), " TGI")
  }
  rownames(samples) <- NULL
  samples
}

figure7_prepare_cellcycle <- function(cellcycle, samples, config) {
  tgi_measure <- figure7_tgi_measure(config)
  data <- cellcycle[is.finite(cellcycle$pseudotime), , drop = FALSE]
  index <- match(data$sample_id, samples$sample_id)
  if (anyNA(index)) figure7_stop("Missing sample metadata for CellCycle rows")
  data$dose <- samples$dose[index]; data$dose_mg <- samples$dose_mg[index]
  data$sample_mean_endpoint_ploidy <- samples$sample_mean_endpoint_ploidy[index]
  data$etp_group <- samples$etp_group[index]
  data[[tgi_measure]] <- samples[[tgi_measure]][index]
  data
}
