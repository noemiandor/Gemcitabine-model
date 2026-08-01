#!/usr/bin/env Rscript

test_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(test_file_arg)) stop("Cannot resolve test path", call. = FALSE)
test_path <- normalizePath(
  sub("^--file=", "", test_file_arg[[1L]]),
  mustWork = TRUE
)
repo_root <- normalizePath(
  file.path(dirname(test_path), "..", "..", "..", ".."),
  mustWork = TRUE
)

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop("Copy-number heatmap test requires pheatmap", call. = FALSE)
}

expect_error <- function(expression, pattern) {
  observed <- tryCatch(
    {
      force(expression)
      ""
    },
    error = function(error) conditionMessage(error)
  )
  if (!nzchar(observed) || !grepl(pattern, observed, fixed = TRUE)) {
    stop(
      "Expected error containing `", pattern, "`; observed: ",
      if (nzchar(observed)) observed else "<no error>",
      call. = FALSE
    )
  }
  invisible(observed)
}

helper <- new.env(parent = globalenv())
sys.source(
  file.path(
    repo_root, "Code", "in-vivo", "SI_figures", "copy_number_heatmap.R"
  ),
  envir = helper
)

cache_root <- file.path(repo_root, "Data", "in-vivo", "SIfigures")
endpoint_audit <- utils::read.csv(
  file.path(cache_root, "si_figure6_endpoint_ploidy_join_audit.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
cells <- utils::read.csv(
  file.path(cache_root, "si_figures_cell_metadata.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
included <- tolower(as.character(cells$included_in_si_figures)) %in%
  c("true", "t", "1")
sample_metadata <- unique(data.frame(
  mouse = cells$sample_id[included],
  initial_ploidy = cells$initial_ploidy[included],
  dose = cells$dose[included],
  stringsAsFactors = FALSE
))

cbs_root <- file.path(repo_root, "Data", "in-vivo", "scRNAseq_Numbat")
reference_root <- file.path(cbs_root, "injected_reference")
ploidy_path <- file.path(repo_root, "Data", "in-vivo", "all_ploidy.tsv")
collection <- helper$si_copy_number_read_collection(
  cbs_root,
  ploidy_path,
  endpoint_audit,
  sample_metadata
)

# Figure 7's two frozen transcriptional tables and the SI endpoint audit must
# describe one identical final-QC tumor universe. This prevents downstream CBS
# consumers from silently reintroducing clusters 3, 4, 9, or 9c.
figure7_cellcycle <- utils::read.csv(
  file.path(
    repo_root, "Data", "in-vivo", "figure7", "processed",
    "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
figure7_noncellcycle <- utils::read.csv(
  file.path(
    repo_root,
    "Data", "in-vivo", "figure7", "processed",
    "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
figure7_tumor <- rbind(figure7_cellcycle, figure7_noncellcycle)
figure7_treated <- figure7_tumor$gemcitabine_dose_mg_per_kg > 0
collection_tumor_ids <- paste(
  collection$cell_annotations$sample_id,
  collection$cell_annotations$cell_id,
  sep = "_"
)
collection_treated <- collection$cell_annotations$dose_mg_per_kg > 0
stopifnot(
  nrow(figure7_cellcycle) == 2881L,
  nrow(figure7_noncellcycle) == 6951L,
  nrow(figure7_tumor) == 9832L,
  !anyDuplicated(figure7_tumor$cell_id),
  !anyDuplicated(collection_tumor_ids),
  setequal(figure7_tumor$cell_id, collection$qc_selection$cell),
  setequal(figure7_tumor$cell_id, collection_tumor_ids),
  sum(figure7_treated) == 5335L,
  sum(collection_treated) == 5335L,
  setequal(
    figure7_tumor$cell_id[figure7_treated],
    collection_tumor_ids[collection_treated]
  ),
  identical(
    as.integer(table(factor(
      figure7_tumor$initial_ploidy[figure7_treated],
      levels = c("2N", "4N")
    ))),
    c(2407L, 2928L)
  ),
  !any(as.character(figure7_tumor$cluster) %in% c("3", "4", "9", "9c"))
)

# The reviewed 16-matrix inventory is checksum-pinned independently of Git's
# current working-tree state.
manifest_path <- file.path(cbs_root, "cbs_manifest.tsv")
reviewed_manifest <- utils::read.delim(
  manifest_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
cbs_fixture <- tempfile("reviewed-cbs-")
dir.create(cbs_fixture)
fixture_sources <- c(
  manifest_path,
  file.path(cbs_root, reviewed_manifest$filename)
)
stopifnot(all(file.copy(fixture_sources, cbs_fixture)))
missing_matrix <- file.path(cbs_fixture, reviewed_manifest$filename[[1L]])
unlink(missing_matrix)
expect_error(
  helper$si_copy_number_validate_manifest(cbs_fixture),
  "CBS manifest inventory mismatch"
)
stopifnot(file.copy(
  file.path(cbs_root, reviewed_manifest$filename[[1L]]),
  missing_matrix
))
write("# tampered after review", file = missing_matrix, append = TRUE)
expect_error(
  helper$si_copy_number_validate_manifest(cbs_fixture),
  "CBS matrix differs from its reviewed checksum"
)
stopifnot(file.copy(
  file.path(cbs_root, reviewed_manifest$filename[[1L]]),
  missing_matrix,
  overwrite = TRUE
))
fixture_manifest <- file.path(cbs_fixture, "cbs_manifest.tsv")
utils::write.table(
  reviewed_manifest[-1L, , drop = FALSE],
  fixture_manifest,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
expect_error(
  helper$si_copy_number_validate_manifest(cbs_fixture),
  "exactly 16 unique pinned matrices"
)
utils::write.table(
  rbind(reviewed_manifest, reviewed_manifest[1L, , drop = FALSE]),
  fixture_manifest,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
expect_error(
  helper$si_copy_number_validate_manifest(cbs_fixture),
  "exactly 16 unique pinned matrices"
)
unlink(cbs_fixture, recursive = TRUE)

references <- helper$si_copy_number_read_injected_references(reference_root)
stopifnot(
  nrow(references$manifest) == 2L,
  identical(references$manifest$injected_origin, c("2N", "4N")),
  nrow(references$cells) == 36L,
  sum(references$cells$injected_origin == "2N") == 20L,
  sum(references$cells$injected_origin == "4N") == 16L,
  all(references$cells$source_repository == "miningcloneid"),
  all(references$cells$source_commit ==
    "c505cd9159fa2a8c0974c7379f6aacd09fe19abc"),
  all(references$cells$policy_source_commit ==
    "c0051b17e375703b20e32fd3c9258263138b16dd"),
  all(references$cells$policy_source_locator == paste0(
    "code/beam_search_flip_rate_wgd.py:",
    "load_initial_ploidy_from_cbs"
  )),
  all(grepl("proxy; not the same-passage", references$cells$designation_basis,
            fixed = TRUE)),
  all(references$cells$ploidy ==
    references$cells$assigned_autosomal_ploidy *
      (1 + references$cells$extra_dna_fraction)),
  abs(mean(references$cells$assigned_autosomal_ploidy[
    references$cells$injected_origin == "2N"
  ]) - 2.0099748244977409) < 1e-13,
  abs(mean(references$cells$assigned_autosomal_ploidy[
    references$cells$injected_origin == "4N"
  ]) - 3.5156098144651913) < 1e-13,
  abs(mean(references$cells$ploidy[
    references$cells$injected_origin == "2N"
  ]) - 2.2933485709302355) < 1e-13,
  abs(mean(references$cells$ploidy[
    references$cells$injected_origin == "4N"
  ]) - 4.9862311678488558) < 1e-13,
  abs(min(references$cells$ploidy[
    references$cells$injected_origin == "4N"
  ]) - 4.0171010312680098) < 1e-13,
  abs(max(references$cells$ploidy[
    references$cells$injected_origin == "4N"
  ]) - 6.1486064402620189) < 1e-13
)

reference_fixture <- tempfile("injected-reference-")
dir.create(reference_fixture)
reference_manifest_path <- file.path(reference_root, "reference_manifest.tsv")
reference_sources <- c(
  reference_manifest_path,
  file.path(reference_root, references$manifest$filename)
)
stopifnot(all(file.copy(reference_sources, reference_fixture)))
tampered_reference <- file.path(
  reference_fixture, references$manifest$filename[[2L]]
)
write("# tampered after review", file = tampered_reference, append = TRUE)
expect_error(
  helper$si_copy_number_validate_reference_manifest(reference_fixture),
  "Injected-cell reference matrix differs from its checksum"
)
stopifnot(file.copy(
  file.path(reference_root, references$manifest$filename[[2L]]),
  tampered_reference,
  overwrite = TRUE
))
missing_reference <- file.path(
  reference_fixture, references$manifest$filename[[1L]]
)
unlink(missing_reference)
expect_error(
  helper$si_copy_number_validate_reference_manifest(reference_fixture),
  "inventory differs from its manifest"
)
stopifnot(file.copy(
  file.path(reference_root, references$manifest$filename[[1L]]),
  missing_reference
))
reference_manifest <- utils::read.delim(
  file.path(reference_fixture, "reference_manifest.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
reference_manifest$extra_dna_policy[[1L]] <- "exclude chr999"
utils::write.table(
  reference_manifest,
  file.path(reference_fixture, "reference_manifest.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
expect_error(
  helper$si_copy_number_validate_reference_manifest(reference_fixture),
  "must pin one valid 2N and one valid 4N matrix"
)
unlink(reference_fixture, recursive = TRUE)

matrix_widths <- vapply(
  collection$matrices,
  function(item) ncol(item$values),
  integer(1L)
)
cell_ids <- collection$cell_annotations$cell_id
repeated_barcodes <- unique(cell_ids[duplicated(cell_ids)])
stopifnot(
  length(collection$matrices) == 16L,
  identical(names(collection$matrices), reviewed_manifest$filename),
  nrow(collection$ploidy) == 14125L,
  nrow(collection$qc_selection) == 9832L,
  nrow(collection$cell_annotations) == 9832L,
  sum(collection$cell_annotations$dose_mg_per_kg > 0) == 5335L,
  sum(collection$cell_annotations$dose_mg_per_kg == 0) == 4497L,
  identical(
    as.integer(table(factor(
      collection$cell_annotations$initial_ploidy,
      levels = c("2N", "4N")
    ))),
    c(4880L, 4952L)
  ),
  identical(sort(unique(matrix_widths)), c(37L, 45L)),
  sum(matrix_widths == 37L) == 8L,
  sum(matrix_widths == 45L) == 8L,
  length(repeated_barcodes) == 8L,
  !anyDuplicated(collection$cell_annotations$heatmap_row_id),
  setequal(unique(collection$cell_annotations$initial_ploidy), c("2N", "4N"))
)

harmonized <- helper$si_copy_number_harmonize(collection)
stopifnot(
  identical(dim(harmonized$matrix), c(9832L, 22L)),
  identical(colnames(harmonized$matrix), paste0("chr", seq_len(22L))),
  nrow(harmonized$chromosome_schema_audit) == 44L,
  nrow(harmonized$chromosome_file_audit) == 352L,
  setequal(harmonized$chromosome_schema_audit$schema_id, c("schema_1", "schema_2")),
  all(harmonized$chromosome_schema_audit$n_available_segments >= 1L),
  all(harmonized$chromosome_schema_audit$represented_bp > 0),
  all(harmonized$chromosome_schema_audit$
    represented_bp_fraction_of_exported > 0),
  all(harmonized$chromosome_schema_audit$
    represented_bp_fraction_of_exported <= 1),
  all(harmonized$chromosome_file_audit$minimum_cell_available_bp_fraction > 0),
  all(harmonized$chromosome_file_audit$maximum_cell_available_bp_fraction <= 1),
  all(harmonized$chromosomes$fraction_all_cells_available == 1),
  identical(
    unique(harmonized$cell_annotations$initial_ploidy),
    c("2N", "4N")
  )
)

# The 2N and 4N schemas have materially different represented spans. The
# production reduction is therefore chromosome-level and never a false
# coordinate alignment or direct segment-column bind.
schema_bp <- stats::aggregate(
  represented_bp ~ schema_id + initial_ploidy,
  data = harmonized$chromosome_schema_audit,
  FUN = sum
)
stopifnot(
  nrow(schema_bp) == 2L,
  length(unique(schema_bp$represented_bp)) == 2L
)
# The centromeric exports on chr13-15 are all missing and therefore must be
# reported as exported but unavailable, rather than silently counted as
# represented genomic span.
centromeric_gap_audit <- harmonized$chromosome_schema_audit[
  harmonized$chromosome_schema_audit$chromosome %in% 13:15,
  ,
  drop = FALSE
]
stopifnot(
  nrow(centromeric_gap_audit) == 6L,
  all(centromeric_gap_audit$n_all_missing_exported_segments == 1L),
  all(centromeric_gap_audit$n_available_segments ==
    centromeric_gap_audit$n_exported_segments - 1L),
  all(centromeric_gap_audit$represented_bp_fraction_of_exported < 1)
)

heatmap <- helper$si_copy_number_heatmap(harmonized)
stopifnot(
  inherits(heatmap$gtable, "gtable"),
  identical(names(heatmap$row_annotation), c("Injected origin", "Mouse")),
  identical(heatmap$labels_col, paste0("chr", seq_len(22L))),
  identical(heatmap$gaps_col, seq_len(21L)),
  identical(heatmap$cluster_rows, FALSE),
  identical(heatmap$cluster_cols, FALSE),
  length(heatmap$gaps_row) == 15L
)

sample_summary <- helper$si_copy_number_sample_summary(
  harmonized$cell_annotations
)
endpoint_summaries <- helper$si_copy_number_endpoint_summaries(sample_summary)
reduction <- helper$si_copy_number_reduction_summary(
  references$cells,
  sample_summary,
  harmonized$cell_annotations
)
separation <- helper$si_copy_number_separation_summary(reduction)
stopifnot(
  nrow(sample_summary) == 16L,
  "mean_postprocessed_copy_number_score" %in% names(sample_summary),
  !any(grepl("endpoint_ploidy", names(sample_summary), fixed = TRUE)),
  identical(endpoint_summaries$origin_summary$analysis_type,
            c("descriptive_only", "descriptive_only")),
  all(!endpoint_summaries$origin_summary$formal_test_performed),
  all(!endpoint_summaries$origin_dose_summary$formal_test_performed),
  !any(grepl(
    "p_value|permutation|coefficient|statistic",
    c(
      names(endpoint_summaries$samples),
      names(endpoint_summaries$origin_summary),
      names(endpoint_summaries$origin_dose_summary)
    ),
    ignore.case = TRUE
  )),
  isTRUE(all.equal(
    endpoint_summaries$origin_summary$mean_of_mouse_means,
    c(2.1355132970092279, 2.3185529494402664),
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    endpoint_summaries$origin_summary$min_mouse_mean,
    c(1.9414587666013141, 2.2136453794185855),
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    endpoint_summaries$origin_summary$max_mouse_mean,
    c(2.3157000073747946, 2.5855181343704836),
    tolerance = 1e-14
  )),
  identical(reduction$injected_origin, c("2N", "4N")),
  identical(reduction$n_reference_cells, c(20L, 16L)),
  identical(reduction$n_endpoint_mice, c(8L, 8L)),
  identical(reduction$n_endpoint_cells, c(4880L, 4952L)),
  all(reduction$analysis_type == "descriptive_only"),
  all(!reduction$formal_test_performed),
  isTRUE(all.equal(
    reduction$reference_mean_assigned_autosomal_ploidy,
    c(2.0099748244977409, 3.5156098144651913),
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    reduction$reference_mean_chr999_extra_dna_fraction,
    c(0.14126812865550156, 0.43090318524437865),
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    reduction$reference_mean_ploidy,
    c(2.2933485709302355, 4.9862311678488558),
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    reduction$endpoint_mouse_balanced_mean_ploidy,
    c(2.1355132970092279, 2.3185529494402664),
    tolerance = 1e-14
  )),
  abs(reduction$absolute_change[reduction$injected_origin == "2N"] -
        (-0.1578352739210076)) < 1e-13,
  abs(reduction$relative_change_percent[
    reduction$injected_origin == "2N"
  ] - (-6.8823063323943838)) < 1e-11,
  abs(reduction$absolute_change[reduction$injected_origin == "4N"] -
        (-2.6676782184085894)) < 1e-13,
  abs(reduction$relative_change_percent[
    reduction$injected_origin == "4N"
  ] - (-53.500893332217302)) < 1e-11,
  isTRUE(reduction$all_endpoint_mouse_means_below_reference_min[
    reduction$injected_origin == "4N"
  ]),
  isTRUE(reduction$all_endpoint_cells_below_reference_min[
    reduction$injected_origin == "4N"
  ]),
  identical(separation$analysis_type, "descriptive_only"),
  identical(separation$formal_test_performed, FALSE),
  abs(separation$reference_4n_minus_2n_mean_ploidy -
        2.6928825969186203) < 1e-13,
  abs(separation$endpoint_4n_minus_2n_mouse_balanced_mean_ploidy -
        0.18303965243103848) < 1e-13,
  abs(separation$separation_contraction_percent -
        93.202835777523873) < 1e-11
)

# The endpoint audit is the frozen QC mask. It must identify exactly the final
# 9,832-cell tumor universe and cannot silently admit one of the discarded
# clusters (or omit a retained cell).
bad_endpoint_audit <- endpoint_audit
matched_position <- which(helper$si_copy_number_truthy(
  bad_endpoint_audit$matched
))[[1L]]
bad_endpoint_audit$matched[[matched_position]] <- FALSE
expect_error(
  helper$si_copy_number_read_collection(
    cbs_root,
    ploidy_path,
    bad_endpoint_audit,
    sample_metadata
  ),
  "exact 9,832-cell QC-passed tumor universe"
)

# A same-count swap to a discarded CBS barcode must also fail. Counts alone
# are not an inclusion contract: the full metadata cell ID must bind the exact
# sample-specific endpoint barcode.
bad_endpoint_audit <- endpoint_audit
matched_positions <- which(helper$si_copy_number_truthy(
  bad_endpoint_audit$matched
))
matched_position <- matched_positions[[1L]]
matched_file <- bad_endpoint_audit$endpoint_file[[matched_position]]
selected_barcodes <- bad_endpoint_audit$endpoint_cell_id[
  helper$si_copy_number_truthy(bad_endpoint_audit$matched) &
    bad_endpoint_audit$endpoint_file == matched_file
]
unused_barcode <- setdiff(
  collection$ploidy$cell_id[collection$ploidy$file == matched_file],
  selected_barcodes
)[[1L]]
unused_ploidy <- collection$ploidy$ploidy[
  collection$ploidy$file == matched_file &
    collection$ploidy$cell_id == unused_barcode
][[1L]]
bad_endpoint_audit$endpoint_cell_id[[matched_position]] <- unused_barcode
bad_endpoint_audit$endpoint_ploidy[[matched_position]] <- unused_ploidy
expect_error(
  helper$si_copy_number_read_collection(
    cbs_root,
    ploidy_path,
    bad_endpoint_audit,
    sample_metadata
  ),
  "exact 9,832-cell QC-passed tumor universe"
)

# Filename-encoded dose must agree with the independently reviewed metadata.
bad_sample_metadata <- sample_metadata
bad_sample_metadata$dose[[1L]] <- "120mg/kg"
expect_error(
  helper$si_copy_number_read_collection(
    cbs_root,
    ploidy_path,
    endpoint_audit,
    bad_sample_metadata
  ),
  "disagrees with reviewed sample metadata"
)

# The matrices recompute and validate the tracked ploidy/coverage values.
tampered_ploidy <- helper$si_copy_number_read_ploidy(ploidy_path)
tampered_ploidy$ploidy[[1L]] <- tampered_ploidy$ploidy[[1L]] + 0.01
tampered_path <- tempfile(fileext = ".tsv")
utils::write.table(
  tampered_ploidy,
  tampered_path,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
expect_error(
  helper$si_copy_number_read_collection(
    cbs_root,
    tampered_path,
    endpoint_audit,
    sample_metadata
  ),
  "CBS-derived ploidy disagrees"
)

generator_text <- paste(readLines(file.path(
  repo_root,
  "Code", "in-vivo", "SI_figures", "generate_supplementary_figures.R"
), warn = FALSE), collapse = "\n")
stopifnot(
  grepl("si_copy_number_harmonize(cbs_collection)", generator_text, fixed = TRUE),
  !grepl("bin_size", generator_text, fixed = TRUE),
  grepl("numbat_cbs_matrix_", generator_text, fixed = TRUE),
  grepl("endpoint_ploidy_derivation_helper", generator_text, fixed = TRUE),
  grepl("mean_postprocessed_copy_number_score", generator_text, fixed = TRUE),
  !grepl("Endpoint ploidy in", generator_text, fixed = TRUE),
  !grepl('name = "Endpoint ploidy"', generator_text, fixed = TRUE),
  grepl(
    "Supplementary Figure 6 | Endpoint tumor ploidy and",
    generator_text,
    fixed = TRUE
  ),
  !grepl("exact within-dose P", generator_text, fixed = TRUE),
  !grepl(
    "si6_postprocessed_copy_number_score_p_value",
    generator_text,
    fixed = TRUE
  ),
  !grepl("si_copy_number_endpoint_statistics", generator_text, fixed = TRUE),
  grepl("si_copy_number_read_injected_references", generator_text, fixed = TRUE),
  grepl("si_copy_number_reduction_summary", generator_text, fixed = TRUE),
  grepl("si_copy_number_separation_summary", generator_text, fixed = TRUE),
  grepl("Injected-reference to endpoint ploidy", generator_text, fixed = TRUE),
  grepl("project-designated proxy cells (chr999 included)",
        generator_text, fixed = TRUE),
  grepl("descriptive_only", generator_text, fixed = TRUE),
  grepl("relative_change_percent", generator_text, fixed = TRUE),
  grepl("separation_contraction_percent", generator_text, fixed = TRUE),
  grepl("seed = plot_seed + 61L", generator_text, fixed = TRUE),
  grepl("seed = plot_seed + 62L", generator_text, fixed = TRUE),
  grepl(
    'shared_context_add_tag(s6e_mouse, "F")',
    generator_text,
    fixed = TRUE
  )
)

message("SI6E copy-number heatmap and ploidy-reduction tests passed.")
