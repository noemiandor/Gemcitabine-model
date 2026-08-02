#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(file_arg) != 1L) stop("Cannot resolve test path", call. = FALSE)
test_path <- normalizePath(file_arg)
repo_root <- normalizePath(file.path(dirname(test_path), "../../../.."))
source(file.path(repo_root, "Code/in-vivo/flow_cytometry/reconstruct_endpoint_flow.R"))

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

run_test <- function() {
  assert_true(agreement_tolerance(18000L, 18200L) == 19L,
    "Agreement tolerance changed unexpectedly")
  assert_true(agreement_tolerance(0L, 18000L) == 5L,
    "Sparse-gate tolerance must retain its five-event floor")
  assert_true(identical(deterministic_indices(5L, 10L), 1:5),
    "Small deterministic sample changed")
  assert_true(length(deterministic_indices(1000L, 100L)) == 100L,
    "Deterministic downsampling size changed")

  output_dir <- tempfile("endpoint_flow_integration_")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  dataset <- file.path(repo_root, "Data/in-vivo/flow_cytometry/endpoint_tumors_20250128")
  command <- c(
    file.path(repo_root, "Code/in-vivo/flow_cytometry/run_endpoint_flow.sh"),
    "--crosswalk", file.path(dataset, "crosswalk.tsv"),
    "--paired-sensitivity", file.path(dataset, "paired_acquisition_sensitivity.tsv"),
    "--workspace", file.path(dataset, "workspace/20250129_TumorSamples.wsp"),
    "--expected-input-manifest", file.path(dataset, "content_manifest.tsv"),
    "--fcs-dir", file.path(dataset, "fcs"),
    "--expected-samples", "16",
    "--min-human-cells", "1000",
    "--output-dir", output_dir
  )
  output <- system2("bash", command, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status") %||% 0L
  if (status != 0L) {
    cat(paste(output, collapse = "\n"), "\n", file = stderr())
    stop("Endpoint-flow integration command failed", call. = FALSE)
  }

  expected_files <- c(
    "endpoint_flow_per_sample.tsv",
    "endpoint_flow_group_summary.tsv",
    "input_hashes.tsv",
    "paired_workspace_sensitivity.tsv",
    "endpoint_flow_count_agreement.tsv",
    "endpoint_flow_gate_geometry.tsv",
    "endpoint_flow_dna_histograms.tsv",
    "endpoint_flow_reconstruction_per_sample.tsv",
    "endpoint_flow_representative_selection.tsv",
    "endpoint_flow_named_2n_4n_summary.tsv",
    "metadata/run_config.tsv",
    "figures/panel_SuppFig9_endpoint_flow_cytometry.pdf",
    "figures/panel_SuppFig9_endpoint_flow_cytometry.png"
  )
  assert_true(all(file.exists(file.path(output_dir, expected_files))),
    "One or more endpoint-flow reconstruction outputs are missing")
  assert_true(all(file.info(file.path(output_dir, expected_files))$size > 0),
    "One or more endpoint-flow outputs are empty")

  read_png_dimensions <- function(path) {
    connection <- file(path, "rb")
    on.exit(close(connection))
    prefix <- readBin(connection, what = "raw", n = 16L)
    assert_true(identical(rawToChar(prefix[13:16]), "IHDR"), "Figure PNG lacks an IHDR header")
    c(width = readBin(connection, what = "integer", n = 1L, size = 4L, endian = "big"),
      height = readBin(connection, what = "integer", n = 1L, size = 4L, endian = "big"))
  }
  png_dimensions <- read_png_dimensions(
    file.path(output_dir, "figures/panel_SuppFig9_endpoint_flow_cytometry.png"))
  assert_true(png_dimensions[["width"]] == 2130L && png_dimensions[["height"]] == 2700L,
    "Figure PNG is not the required 7.1-by-9-inch 300-dpi canvas")

  agreement <- read_tsv(file.path(output_dir, "endpoint_flow_count_agreement.tsv"), "agreement output")
  assert_true(nrow(agreement) == 80L, "Agreement output must contain five gates for each of 16 mice")
  assert_true(length(unique(agreement$mouse_id)) == 16L, "Agreement output lost a mouse")
  assert_true(all(agreement$agreement_status == "close_match"), "A gate replay is outside tolerance")
  observed_max <- aggregate(absolute_delta_events ~ population_role, agreement, max)
  expected_max <- c(human_cell_enrichment = 21L, human_cells = 15L,
    reviewed_peak = 15L, workspace_named_2N = 15L, workspace_named_4N = 2L)
  assert_true(all(observed_max$absolute_delta_events == expected_max[observed_max$population_role]),
    "Pinned full-cohort replay deltas changed")

  named <- read_tsv(file.path(output_dir, "endpoint_flow_named_2n_4n_summary.tsv"), "named-gate output")
  assert_true(nrow(named) == 32L, "Named-gate table must contain 2N and 4N rows for all mice")
  low_two <- named[named$mouse_id == "2N-A1-0" & named$population_name == "2N", , drop = FALSE]
  assert_true(nrow(low_two) == 1L && low_two$workspace_parent_count == 174L &&
      low_two$reconstructed_parent_count == 173L && isTRUE(low_two$low_human_cells_flag),
    "The 174-event sample was not retained and flagged correctly")

  histograms <- read_tsv(file.path(output_dir, "endpoint_flow_dna_histograms.tsv"), "histogram output")
  assert_true(nrow(histograms) == 16L * 115L, "Histogram output does not use 115 common bins per mouse")
  assert_true(length(unique(histograms$bin_left)) == 115L &&
      min(histograms$bin_left) == 25000 && max(histograms$bin_right) == 140000,
    "Histogram raw-axis contract changed")
  mass <- rowsum(histograms$probability_mass, histograms$mouse_id)
  assert_true(all(abs(mass - 1) < 1e-12), "Within-mouse histogram mass is not one")
  low_hist <- histograms[histograms$mouse_id == "2N-A1-0", , drop = FALSE]
  assert_true(sum(low_hist$event_count) == 173L && grepl("replay n=173", low_hist$facet_label[[1L]], fixed = TRUE),
    "The low-count distribution is not labelled with its replayed event count")

  geometry <- read_tsv(file.path(output_dir, "endpoint_flow_gate_geometry.tsv"), "gate geometry output")
  assert_true(length(unique(geometry$mouse_id)) == 16L, "Gate geometry does not cover all mice")
  required_roles <- c("human_cell_enrichment", "human_cells", "workspace_named_2N",
    "workspace_named_4N", "reviewed_peak")
  unique_roles <- unique(geometry[, c("mouse_id", "population_role")])
  role_counts <- table(unique_roles$mouse_id)
  assert_true(length(role_counts) == 16L && all(role_counts == 5L) &&
      all(required_roles %in% unique_roles$population_role),
    "Gate geometry does not contain exactly five required roles per mouse")
  hce_geometry <- geometry[geometry$population_role == "human_cell_enrichment", , drop = FALSE]
  human_geometry <- geometry[geometry$population_role == "human_cells", , drop = FALSE]
  dna_geometry <- geometry[geometry$population_role %in% c("workspace_named_2N",
    "workspace_named_4N", "reviewed_peak"), , drop = FALSE]
  assert_true(all(hce_geometry$x_parameter == "FSC-A" & hce_geometry$y_parameter == "SSC-A") &&
      all(human_geometry$x_parameter == "450_50 Violet B-A" & human_geometry$y_parameter == "FSC-A") &&
      all(dna_geometry$x_parameter == "450_50 Violet B-A" & is.na(dna_geometry$y_parameter)),
    "Gate geometry dimensions changed")

  selection <- read_tsv(file.path(output_dir, "endpoint_flow_representative_selection.tsv"), "representative output")
  selected <- selection$mouse_id[selection$selected]
  assert_true(identical(selected, "4N-A8-RR"), "Representative selection is no longer deterministic")
  selected_row <- selection[selection$selected, , drop = FALSE]
  assert_true(selected_row$displayed_all_events == 10000L &&
      selected_row$displayed_human_cell_enrichment_events == 10000L &&
      grepl("all events used for gate counts", selected_row$display_sampling_policy, fixed = TRUE),
    "Representative scatter downsampling is not disclosed correctly")

  bad_crosswalk <- read_tsv(file.path(dataset, "crosswalk.tsv"), "crosswalk")
  bad_crosswalk$fcs_file[[1L]] <- bad_crosswalk$fcs_file[[2L]]
  bad_crosswalk_path <- file.path(output_dir, "bad_crosswalk.tsv")
  write_tsv_atomic(bad_crosswalk, bad_crosswalk_path)
  bad_config <- list(
    crosswalk = bad_crosswalk_path,
    workspace = file.path(dataset, "workspace/20250129_TumorSamples.wsp"),
    fcs_dir = file.path(dataset, "fcs"),
    frozen_table = file.path(output_dir, "endpoint_flow_per_sample.tsv"),
    output_dir = file.path(output_dir, "should_not_exist"),
    expected_samples = 16L,
    min_human_cells = 1000L,
    dna_axis_min = 25000,
    dna_axis_max = 140000,
    dna_bin_width = 1000,
    representative_origin = "4N",
    max_representative_points = 10000L
  )
  mismatch_error <- tryCatch({
    reconstruct_endpoint_flow(bad_config)
    ""
  }, error = function(e) conditionMessage(e))
  assert_true(grepl("Crosswalk and frozen table disagree", mismatch_error, fixed = TRUE),
    "Renderer did not reject a frozen-count/FCS acquisition mismatch")
  cat("Endpoint-flow reconstruction integration test: PASS\n")
}

run_test()
