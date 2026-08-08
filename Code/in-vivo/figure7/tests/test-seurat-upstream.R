source(file.path(module_dir, "src", "seurat_upstream.R"), local = FALSE)
if (!exists("figure7_select_seurat_source", mode = "function")) {
  source(
    file.path(module_dir, "src", "seurat_upstream_selection.R"),
    local = FALSE
  )
}

testthat::test_that("manual merge and final metadata preserve reviewed labels", {
  refined <- c(
    "0", "1", "7", "2", "3", "4", "4c", "5", "6", "8",
    "9", "9c", "10", "11", "12", "13", "14"
  )
  merged <- figure7_upstream_manual_merge_values(refined, refined)
  testthat::expect_identical(
    as.character(merged),
    c(
      "0", "0", "0", "2", "3", "4", "4c", "5", "6", "8",
      "9", "9c", "10", "10", "10", "13", "14"
    )
  )
  testthat::expect_identical(
    levels(merged),
    c("0", "2", "3", "4", "4c", "5", "6", "8", "9", "9c",
      "10", "13", "14")
  )

  metadata <- data.frame(
    IDs = c(
      "2N-Cell-Culture", "4N-Cell-Culture", "2N-A1-0", "4N-A5-0"
    ),
    stringsAsFactors = FALSE
  )
  observed <- figure7_upstream_derive_tn_ploidy(metadata)
  testthat::expect_identical(
    as.character(observed$Ploidy),
    c("2N", "4N", "2N", "4N")
  )
  testthat::expect_identical(
    as.character(observed$TN),
    c("CellLine", "CellLine", "Tumor", "Tumor")
  )
  testthat::expect_error(
    figure7_upstream_derive_tn_ploidy(
      data.frame(IDs = "unknown", stringsAsFactors = FALSE)
    ),
    "Cannot derive"
  )
})

testthat::test_that("final object mutation filters only the reviewed clusters", {
  testthat::skip_if_not_installed("Seurat")
  labels <- c(
    "0", "2", "3", "4", "4c", "5", "6", "8", "9", "9c",
    "10", "13", "14"
  )
  cells <- paste0("cell_", seq_along(labels))
  counts <- Matrix::Matrix(
    matrix(
      seq_len(3L * length(cells)),
      nrow = 3L,
      dimnames = list(paste0("gene", 1:3), cells)
    ),
    sparse = TRUE
  )
  object <- Seurat::CreateSeuratObject(counts)
  object$sample <- "sample"
  object$IDs <- ifelse(
    seq_along(labels) %% 2L,
    "2N-A1-0",
    "4N-Cell-Culture"
  )
  object$manual_merge_test <- factor(labels, levels = labels)
  final <- figure7_upstream_finalize_clusters(
    object,
    rerun_reductions = FALSE
  )
  testthat::expect_identical(
    as.character(final$clusters),
    setdiff(labels, c("3", "4", "9", "9c"))
  )
  testthat::expect_false("manual_merge_test" %in% colnames(final@meta.data))
  testthat::expect_identical(
    levels(final$clusters),
    setdiff(labels, c("3", "4", "9", "9c"))
  )
  testthat::expect_setequal(as.character(final$Ploidy), c("2N", "4N"))
  testthat::expect_setequal(as.character(final$TN), c("Tumor", "CellLine"))
})

testthat::test_that("cell-cycle candidate rule retains all three reviewed methods", {
  clusters <- rep(c("0", "1", "2", "6"), each = 8L)
  s_score <- ifelse(clusters == "6", 10, 0)
  metadata <- data.frame(
    seurat_clusters = clusters,
    S.Score = s_score,
    G2M.Score = rep(0, length(clusters)),
    Phase = ifelse(clusters == "6", "S", "G1"),
    stringsAsFactors = FALSE
  )
  pca <- matrix(0, nrow = nrow(metadata), ncol = 20L)
  pca[, 1L] <- s_score
  colnames(pca) <- paste0("PC_", seq_len(ncol(pca)))

  candidates <- suppressWarnings(
    figure7_upstream_cell_cycle_candidate_table(metadata, pca)
  )
  testthat::expect_identical(candidates$cluster, c("0", "1", "2", "6"))
  testthat::expect_identical(
    candidates$final_candidate_flag,
    c(FALSE, FALSE, FALSE, TRUE)
  )
  testthat::expect_identical(candidates$n_methods_flagged, c(0L, 0L, 0L, 3L))
  testthat::expect_true(candidates$flag_score_outlier[[4L]])
  testthat::expect_true(candidates$flag_phase_enrichment[[4L]])
  testthat::expect_true(candidates$flag_pc_separation_support[[4L]])
})

testthat::test_that("UMAP refinement applies the reviewed hull and component rule", {
  clusters <- c("0", "3", "6", "10", "11", "12", "4", "9")
  coordinates <- rbind(
    c(-1, -1),
    c(1, -1),
    c(-1, 1),
    c(1, 1),
    c(0, 0.8),
    c(0, -0.8),
    c(0, 0),
    c(10, 10)
  )
  refined <- figure7_upstream_refine_cluster_values(
    clusters,
    coordinates,
    knn_k = 3L
  )
  testthat::expect_identical(
    as.character(refined$refined),
    c("0", "3", "6", "10", "11", "12", "4c", "9")
  )
  testthat::expect_identical(refined$seed_index, 7L)
  testthat::expect_true(refined$inside_hull[[7L]])
  testthat::expect_false(refined$inside_hull[[8L]])
})

testthat::test_that("published stage acceptance counts are internally closed", {
  testthat::expect_identical(
    sum(figure7_upstream_expected_counts("integrated")),
    42884L
  )
  testthat::expect_identical(
    sum(figure7_upstream_expected_counts("refined")),
    42884L
  )
  testthat::expect_identical(
    sum(figure7_upstream_expected_counts("merged")),
    42884L
  )
  testthat::expect_identical(
    sum(figure7_upstream_expected_counts("final")),
    35513L
  )
  expanded <- rep(
    names(figure7_upstream_expected_counts("final")),
    figure7_upstream_expected_counts("final")
  )
  testthat::expect_silent(
    figure7_upstream_assert_counts(expanded, "final")
  )
  expanded[[1L]] <- "unexpected"
  testthat::expect_error(
    figure7_upstream_assert_counts(expanded, "final"),
    "differ from the reviewed contract"
  )
})

testthat::test_that("final metadata validator fixes the reviewed cell universe", {
  final_counts <- figure7_upstream_expected_counts("final")
  clusters <- rep(names(final_counts), final_counts)
  context <- c(rep("CellLine", 25681L), rep("Tumor", 9832L))
  ploidy <- c(
    rep("2N", 14836L),
    rep("4N", 10845L),
    rep("2N", 4880L),
    rep("4N", 4952L)
  )
  dose <- c(
    rep(NA_real_, 25681L),
    rep(0, 2473L), rep(30, 822L), rep(120, 1585L),
    rep(0, 2024L), rep(30, 849L), rep(120, 2079L)
  )
  expected_samples <- figure7_upstream_expected_tumor_samples()
  testthat::expect_identical(
    stats::setNames(expected_samples$n_cells, expected_samples$sample),
    figure7_curated_endpoint_counts()
  )
  testthat::expect_identical(
    expected_samples$Ploidy,
    c(rep("2N", 8L), rep("4N", 8L))
  )
  testthat::expect_identical(
    expected_samples$Dose,
    c(
      rep(0, 4L), rep(30, 2L), rep(120, 2L),
      rep(0, 4L), rep(30, 2L), rep(120, 2L)
    )
  )
  sample <- c(
    rep("2N-Cell-Culture", 14836L),
    rep("4N-Cell-Culture", 10845L),
    rep(expected_samples$sample, expected_samples$n_cells)
  )
  metadata <- data.frame(
    clusters = clusters,
    TN = context,
    Ploidy = ploidy,
    Dose = dose,
    sample = sample,
    stringsAsFactors = FALSE
  )

  testthat::expect_silent(
    figure7_upstream_validate_final_metadata(metadata)
  )

  too_few_treated <- metadata
  treated_2n <- which(
    too_few_treated$TN == "Tumor" &
      too_few_treated$Ploidy == "2N" &
      too_few_treated$Dose == 30
  )[[1L]]
  too_few_treated$Dose[[treated_2n]] <- 0
  testthat::expect_error(
    figure7_upstream_validate_final_metadata(too_few_treated),
    "treated Tumor count differs"
  )

  wrong_treated_split <- metadata
  control_4n <- which(
    wrong_treated_split$TN == "Tumor" &
      wrong_treated_split$Ploidy == "4N" &
      wrong_treated_split$Dose == 0
  )[[1L]]
  wrong_treated_split$Dose[[treated_2n]] <- 0
  wrong_treated_split$Dose[[control_4n]] <- 30
  testthat::expect_error(
    figure7_upstream_validate_final_metadata(wrong_treated_split),
    "treated Tumor Ploidy counts differ"
  )

  wrong_sample_count <- metadata
  sample_row <- which(wrong_sample_count$sample == "2N-A2-0")[[1L]]
  wrong_sample_count$sample[[sample_row]] <- "2N-A2-L"
  testthat::expect_error(
    figure7_upstream_validate_final_metadata(wrong_sample_count),
    "Tumor sample counts differ"
  )

  wrong_sample_dose <- metadata
  wrong_sample_dose$Dose[[sample_row]] <- 120
  testthat::expect_error(
    figure7_upstream_validate_final_metadata(wrong_sample_dose),
    "sample-to-Ploidy/Dose mapping differs"
  )

  leaked_qc_cluster <- metadata
  leaked_qc_cluster$clusters[[1L]] <- "3"
  testthat::expect_error(
    figure7_upstream_validate_final_metadata(leaked_qc_cluster),
    "retains discarded QC cluster"
  )

  dosed_cell_line <- metadata
  dosed_cell_line$Dose[[1L]] <- 30
  testthat::expect_error(
    figure7_upstream_validate_final_metadata(dosed_cell_line),
    "CellLine Dose values must be missing"
  )
})

testthat::test_that("Cell Ranger inventory requires exactly 18 reviewed samples", {
  ids <- paste0("sample_", seq_len(16L))
  sample_info <- data.frame(
    harvest = paste0("harvest_", ids),
    sequencing = ids,
    IDs = ids,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  testthat::expect_identical(
    figure7_upstream_expected_samples(sample_info),
    sort(c(ids, "2N-Cell-Culture", "4N-Cell-Culture"))
  )
  duplicated <- rbind(sample_info, sample_info[1L, , drop = FALSE])
  testthat::expect_error(
    figure7_upstream_expected_samples(duplicated),
    "duplicate"
  )
  empty_root <- tempfile("figure7_empty_cellranger_")
  dir.create(empty_root)
  testthat::expect_error(
    figure7_upstream_h5_inventory(empty_root, sample_info),
    "18-sample contract"
  )
})

testthat::test_that("reviewed Seurat execution remains deterministic", {
  testthat::skip_if_not_installed("future")
  previous_plan <- future::plan()
  on.exit(future::plan(previous_plan), add = TRUE)
  observed <- figure7_upstream_configure_future(jobs = 8L)
  testthat::expect_identical(observed$strategy, "sequential")
  testthat::expect_identical(observed$workers, 1L)
  testthat::expect_identical(observed$requested_workers, 8L)

  figure7_upstream_set_single_thread()
  thread_variables <- c(
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
  )
  testthat::expect_true(all(Sys.getenv(thread_variables) == "1"))
  testthat::expect_identical(Sys.getenv("KMP_DUPLICATE_LIB_OK"), "TRUE")
  testthat::expect_identical(Sys.getenv("KMP_INIT_AT_FORK"), "FALSE")

  seeded_functions <- c(
    annotate = "figure7_upstream_annotate_cell_cycle",
    refine = "figure7_upstream_refine_clusters",
    merge = "figure7_upstream_merge_clusters",
    final = "figure7_upstream_finalize_clusters"
  )
  expected_seeds <- c(annotate = "12345", refine = "12345",
                      merge = "12345", final = "1234")
  for (label in names(seeded_functions)) {
    function_text <- paste(
      deparse(get(seeded_functions[[label]], mode = "function")),
      collapse = "\n"
    )
    testthat::expect_match(
      function_text,
      paste0("set.seed\\(", expected_seeds[[label]], "\\)")
    )
  }
})

testthat::test_that("upstream cache contract is isolated to scientific inputs", {
  generator_environment <- new.env(parent = globalenv())
  sys.source(
    file.path(module_dir, "generate_final_seurat_from_cellranger.R"),
    envir = generator_environment
  )
  lock_path <- file.path(module_dir, "environment_lock.tsv")
  config_path <- file.path(module_dir, "figure7_config.yaml")
  config <- figure7_read_config(config_path)
  baseline <- generator_environment$figure7_upstream_common_contract(
    lock_path,
    config,
    jobs = 1L
  )
  testthat::expect_false(any(c(
    "common_io_sha256",
    "environment_lock_sha256",
    "figure7_config_sha256",
    "helper_sha256",
    "generator_sha256",
    "legacy_source_contract_sha256"
  ) %in% names(baseline)))

  lock <- figure7_read_environment_lock(lock_path)
  python_only_lock <- rbind(
    lock,
    data.frame(
      ecosystem = "python",
      stage = "scvelo",
      package = "test-renderer-only-pin",
      version = "9.9.9",
      stringsAsFactors = FALSE
    )
  )
  python_only_path <- tempfile("figure7_python_only_lock_", fileext = ".tsv")
  utils::write.table(
    python_only_lock,
    python_only_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  python_only <- generator_environment$figure7_upstream_common_contract(
    python_only_path,
    config,
    jobs = 16L
  )
  testthat::expect_identical(python_only, baseline)

  unrelated_config <- config
  unrelated_config$si_figures$plot_shuffle_seed <-
    as.integer(config$si_figures$plot_shuffle_seed) + 1L
  unrelated_config$si_figures$si7_feature_species_policy <-
    "unrelated renderer/species policy test"
  testthat::expect_identical(
    generator_environment$figure7_upstream_common_contract(
      lock_path,
      unrelated_config,
      jobs = 4L
    ),
    baseline
  )

  upstream_lock <- rbind(
    lock,
    data.frame(
      ecosystem = "r",
      stage = "seurat_upstream",
      package = "test-upstream-science-pin",
      version = "1.0.0",
      stringsAsFactors = FALSE
    )
  )
  upstream_lock_path <- tempfile(
    "figure7_upstream_lock_",
    fileext = ".tsv"
  )
  utils::write.table(
    upstream_lock,
    upstream_lock_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  testthat::expect_false(identical(
    generator_environment$figure7_upstream_common_contract(
      upstream_lock_path,
      config,
      jobs = 1L
    ),
    baseline
  ))

  changed_source <- config
  changed_source$versioned_source_artifacts$sample_info$sha256 <-
    paste(rep("9", 64L), collapse = "")
  testthat::expect_false(identical(
    generator_environment$figure7_upstream_common_contract(
      lock_path,
      changed_source,
      jobs = 1L
    ),
    baseline
  ))
})

testthat::test_that("stage code contracts invalidate only affected descendants", {
  generator <- new.env(parent = globalenv())
  sys.source(
    file.path(module_dir, "generate_final_seurat_from_cellranger.R"),
    envir = generator
  )
  stages <- generator$figure7_upstream_stage_order()
  baseline <- stats::setNames(
    vapply(
      stages,
      generator$figure7_upstream_stage_code_contract,
      character(1L)
    ),
    stages
  )
  config <- figure7_read_config(
    file.path(module_dir, "figure7_config.yaml")
  )
  common <- generator$figure7_upstream_common_contract(
    file.path(module_dir, "environment_lock.tsv"),
    config,
    jobs = 1L
  )
  output_root <- tempfile("figure7_stage_code_boundaries_")
  dir.create(output_root)
  definitions <- generator$figure7_upstream_stage_definitions(output_root)
  dependencies <- list()
  h5_samples <- sprintf("sample_%02d", seq_len(18L))
  h5_hashes <- vapply(seq_len(18L), function(index) {
    paste(rep(as.character(index %% 10L), 64L), collapse = "")
  }, character(1L))
  h5_inventory <- generator$figure7_upstream_h5_inventory_digest(
    h5_samples,
    h5_hashes
  )
  for (stage in stages) {
    definition <- definitions[[stage]]
    dir.create(definition$directory, recursive = TRUE)
    writeLines(
      paste("synthetic", stage),
      generator$figure7_upstream_stage_output(definition)
    )
    roles <- generator$figure7_upstream_stage_dependency_roles(stage)
    dependencies[[stage]] <- stats::setNames(
      vapply(seq_along(roles), function(index) {
        paste(rep(as.character(index %% 10L), 64L), collapse = "")
      }, character(1L)),
      roles
    )
    extra <- character()
    if (identical(stage, "integrated")) {
      dependencies[[stage]][["h5_inventory"]] <- h5_inventory
      extra <- c(
        h5_inventory_sha256 = h5_inventory,
        h5_file_count = "18",
        stats::setNames(
          h5_hashes,
          paste0("h5_sha256:", h5_samples)
        )
      )
    }
    generator$figure7_upstream_write_manifest(
      stage,
      definition,
      common,
      dependencies[[stage]],
      extra
    )
  }
  final_manifest_path <- generator$figure7_upstream_stage_manifest(
    definitions$final
  )
  final_manifest <- generator$figure7_upstream_manifest_values(
    final_manifest_path
  )
  changed_audit <- final_manifest
  changed_audit[["audit_helper_sha256"]] <-
    paste(rep("7", 64L), collapse = "")
  changed_audit[["audit_generator_sha256"]] <-
    paste(rep("8", 64L), collapse = "")
  figure7_write_tsv(
    data.frame(
      key = names(changed_audit),
      value = unname(changed_audit),
      stringsAsFactors = FALSE
    ),
    final_manifest_path
  )
  testthat::expect_true(generator$figure7_upstream_manifest_matches(
    "final",
    definitions$final,
    common,
    dependencies$final
  ))
  missing_audit <- changed_audit[
    names(changed_audit) != "audit_helper_sha256"
  ]
  figure7_write_tsv(
    data.frame(
      key = names(missing_audit),
      value = unname(missing_audit),
      stringsAsFactors = FALSE
    ),
    final_manifest_path
  )
  testthat::expect_false(generator$figure7_upstream_manifest_matches(
    "final",
    definitions$final,
    common,
    dependencies$final
  ))
  figure7_write_tsv(
    data.frame(
      key = names(final_manifest),
      value = unname(final_manifest),
      stringsAsFactors = FALSE
    ),
    final_manifest_path
  )
  original_final <- figure7_upstream_finalize_clusters
  original_merge <- figure7_upstream_merge_clusters
  original_integrated <- figure7_upstream_build_integrated
  original_parameter_contract <-
    generator$figure7_upstream_stage_parameter_contract
  on.exit({
    assign(
      "figure7_upstream_finalize_clusters",
      original_final,
      envir = .GlobalEnv
    )
    assign(
      "figure7_upstream_merge_clusters",
      original_merge,
      envir = .GlobalEnv
    )
    assign(
      "figure7_upstream_build_integrated",
      original_integrated,
      envir = .GlobalEnv
    )
    generator$figure7_upstream_stage_parameter_contract <-
      original_parameter_contract
  }, add = TRUE)

  changed_final <- original_final
  body(changed_final) <- substitute(
    {
      stage_contract_test_marker <- TRUE
      BODY
    },
    list(BODY = body(original_final))
  )
  assign(
    "figure7_upstream_finalize_clusters",
    changed_final,
    envir = .GlobalEnv
  )
  final_changed <- stats::setNames(
    vapply(
      stages,
      generator$figure7_upstream_stage_code_contract,
      character(1L)
    ),
    stages
  )
  testthat::expect_identical(
    final_changed[stages != "final"],
    baseline[stages != "final"]
  )
  testthat::expect_false(identical(
    unname(final_changed[["final"]]),
    unname(baseline[["final"]])
  ))
  testthat::expect_true(generator$figure7_upstream_manifest_matches(
    "merged",
    definitions$merged,
    common,
    dependencies$merged
  ))
  testthat::expect_false(generator$figure7_upstream_manifest_matches(
    "final",
    definitions$final,
    common,
    dependencies$final
  ))

  assign(
    "figure7_upstream_finalize_clusters",
    original_final,
    envir = .GlobalEnv
  )
  changed_merge <- original_merge
  formals(changed_merge)$jobs <- 2L
  assign(
    "figure7_upstream_merge_clusters",
    changed_merge,
    envir = .GlobalEnv
  )
  merge_changed <- stats::setNames(
    vapply(
      stages,
      generator$figure7_upstream_stage_code_contract,
      character(1L)
    ),
    stages
  )
  testthat::expect_identical(
    merge_changed[c("integrated", "cell_cycle", "refined", "final")],
    baseline[c("integrated", "cell_cycle", "refined", "final")]
  )
  testthat::expect_false(identical(
    unname(merge_changed[["merged"]]),
    unname(baseline[["merged"]])
  ))
  testthat::expect_false(generator$figure7_upstream_manifest_matches(
    "merged",
    definitions$merged,
    common,
    dependencies$merged
  ))
  testthat::expect_false(generator$figure7_upstream_manifest_matches(
    "final",
    definitions$final,
    common,
    dependencies$final
  ))

  cumulative <- generator$figure7_upstream_cumulative_code_contracts(
    "final"
  )
  testthat::expect_identical(
    unname(cumulative[["scientific_code_contract:merged"]]),
    unname(merge_changed[["merged"]])
  )
  testthat::expect_identical(
    unname(cumulative[["scientific_code_contract:final"]]),
    unname(baseline[["final"]])
  )

  assign(
    "figure7_upstream_merge_clusters",
    original_merge,
    envir = .GlobalEnv
  )
  changed_integrated <- original_integrated
  body(changed_integrated) <- substitute(
    {
      integrated_contract_test_marker <- TRUE
      BODY
    },
    list(BODY = body(original_integrated))
  )
  assign(
    "figure7_upstream_build_integrated",
    changed_integrated,
    envir = .GlobalEnv
  )
  integrated_changed <- stats::setNames(
    vapply(
      stages,
      generator$figure7_upstream_stage_code_contract,
      character(1L)
    ),
    stages
  )
  testthat::expect_false(identical(
    unname(integrated_changed[["integrated"]]),
    unname(baseline[["integrated"]])
  ))
  testthat::expect_identical(
    integrated_changed[stages != "integrated"],
    baseline[stages != "integrated"]
  )
  for (stage in stages) {
    testthat::expect_false(
      generator$figure7_upstream_manifest_matches(
        stage,
        definitions[[stage]],
        common,
        dependencies[[stage]]
      )
    )
  }

  assign(
    "figure7_upstream_build_integrated",
    original_integrated,
    envir = .GlobalEnv
  )
  generator$figure7_upstream_stage_parameter_contract <- function(stage) {
    value <- original_parameter_contract(stage)
    if (identical(stage, "integrated")) {
      paste0(value, ";parameter_contract_test_marker=true")
    } else {
      value
    }
  }
  for (stage in stages) {
    testthat::expect_false(
      generator$figure7_upstream_manifest_matches(
        stage,
        definitions[[stage]],
        common,
        dependencies[[stage]]
      )
    )
  }

  generator$figure7_upstream_stage_parameter_contract <-
    original_parameter_contract
  generator$figure7_upstream_stage_parameter_contract <- function(stage) {
    value <- original_parameter_contract(stage)
    if (identical(stage, "final")) {
      paste0(value, ";parameter_contract_test_marker=true")
    } else {
      value
    }
  }
  testthat::expect_true(generator$figure7_upstream_manifest_matches(
    "merged",
    definitions$merged,
    common,
    dependencies$merged
  ))
  testthat::expect_false(generator$figure7_upstream_manifest_matches(
    "final",
    definitions$final,
    common,
    dependencies$final
  ))
})

testthat::test_that("stage manifests bind the complete dependency chain", {
  generator_environment <- new.env(parent = globalenv())
  sys.source(
    file.path(module_dir, "generate_final_seurat_from_cellranger.R"),
    envir = generator_environment
  )
  output_root <- tempfile("figure7_upstream_manifest_")
  dir.create(output_root)
  definitions <- generator_environment$figure7_upstream_stage_definitions(
    output_root
  )
  definition <- definitions$integrated
  dir.create(definition$directory, recursive = TRUE)
  writeLines(
    "synthetic-cache",
    generator_environment$figure7_upstream_stage_output(definition),
    useBytes = TRUE
  )
  common <- c(
    seurat_upstream_environment_contract_sha256 =
      paste(rep("d", 64L), collapse = ""),
    source_config_contract_sha256 =
      paste(rep("e", 64L), collapse = ""),
    deterministic_workers = "1"
  )
  h5_samples <- sprintf("sample_%02d", seq_len(18L))
  h5_hashes <- vapply(seq_len(18L), function(index) {
    paste(rep(as.character(index %% 10L), 64L), collapse = "")
  }, character(1L))
  inventory_digest <-
    generator_environment$figure7_upstream_h5_inventory_digest(
      h5_samples,
      h5_hashes
    )
  raw_dependencies <- c(
    all_ploidy = paste(rep("1", 64L), collapse = ""),
    sample_info = paste(rep("2", 64L), collapse = ""),
    h5_inventory = inventory_digest
  )
  h5_extra <- c(
    h5_inventory_sha256 = inventory_digest,
    h5_file_count = "18",
    stats::setNames(h5_hashes, paste0("h5_sha256:", h5_samples))
  )
  generator_environment$figure7_upstream_write_manifest(
    "integrated",
    definition,
    common,
    raw_dependencies,
    h5_extra
  )
  manifest_path <- generator_environment$figure7_upstream_stage_manifest(
    definition
  )
  manifest_text <- readLines(manifest_path, warn = FALSE)
  testthat::expect_false(any(grepl(output_root, manifest_text, fixed = TRUE)))
  testthat::expect_true(
    generator_environment$figure7_upstream_manifest_matches(
      "integrated",
      definition,
      common,
      raw_dependencies
    )
  )
  changed_raw <- raw_dependencies
  changed_raw[["sample_info"]] <- paste(rep("9", 64L), collapse = "")
  testthat::expect_false(
    generator_environment$figure7_upstream_manifest_matches(
      "integrated",
      definition,
      common,
      changed_raw
    )
  )
  changed_h5 <- raw_dependencies
  changed_h5[["h5_inventory"]] <- paste(rep("8", 64L), collapse = "")
  testthat::expect_false(
    generator_environment$figure7_upstream_manifest_matches(
      "integrated",
      definition,
      common,
      changed_h5
    )
  )
  manifest <- generator_environment$figure7_upstream_manifest_values(
    manifest_path
  )
  testthat::expect_true(all(c(
    "creator_r_runtime_version",
    "creator_r_platform",
    "creator_r_blas"
  ) %in% names(manifest)))

  stages <- c("cell_cycle", "refined", "merged", "final")
  parents <- c("integrated", "cell_cycle", "refined", "merged")
  stage_dependencies <- list()
  for (index in seq_along(stages)) {
    stage <- stages[[index]]
    stage_definition <- definitions[[stage]]
    dir.create(stage_definition$directory, recursive = TRUE)
    writeLines(
      paste0("synthetic-", stage, "-cache"),
      generator_environment$figure7_upstream_stage_output(
        stage_definition
      ),
      useBytes = TRUE
    )
    child_dependencies <-
      generator_environment$figure7_upstream_descendant_dependencies(
        parents[[index]],
        definitions
      )
    stage_dependencies[[stage]] <- child_dependencies
    testthat::expect_setequal(
      names(child_dependencies),
      generator_environment$figure7_upstream_stage_dependency_roles(stage)
    )
    generator_environment$figure7_upstream_write_manifest(
      stage,
      stage_definition,
      common,
      child_dependencies
    )
    testthat::expect_true(
      generator_environment$figure7_upstream_manifest_matches(
        stage,
        stage_definition,
        common,
        child_dependencies
      )
    )
  }
  writeLines(
    "tampered-cache",
    generator_environment$figure7_upstream_stage_output(definition),
    useBytes = TRUE
  )
  testthat::expect_false(
    generator_environment$figure7_upstream_manifest_matches(
      "integrated",
      definition,
      common,
      raw_dependencies
    )
  )
  changed_integrated_hash <- figure7_sha256(
    generator_environment$figure7_upstream_stage_output(definition)
  )
  for (stage in stages) {
    changed_parent <- stage_dependencies[[stage]]
    changed_parent[["integrated_rds"]] <- changed_integrated_hash
    testthat::expect_false(
      generator_environment$figure7_upstream_manifest_matches(
        stage,
        definitions[[stage]],
        common,
        changed_parent
      )
    )
  }
  testthat::expect_identical(
    unname(generator_environment$figure7_upstream_legacy_sources[
      c(
        "01_data.R", "01a_cell_cycle.R", "02b_cluster_refine.R",
        "02d_manual_cluster_merge.R", "03_final_cluster.R", "Utils.R"
      )
    ]),
    c(
      "1148912de8feee2321e89440570a4b51b628f5421e12cf4ba8d47ba0d757c6d9",
      "558d40178fbba21fac2eac3f09972b464bdb5b90a716eed9e071cfa8a72a5c87",
      "c37da49ed1f5bf8740652cf3cd51d1b36e6d635b45104ef99192baa36e476551",
      "ef033f9984e3937b90c12886a2f2b1682c411825d01bc4194049c9cbec3d56a6",
      "47e255458ad6054a502b15db13b60ee21c25ae48c960507cfcfd4090fc3e8b6e",
      "0145103a7e2a1c16bac2531037e182e0010a4ddf9cc3f944f2076f77108a2ba0"
    )
  )
})

testthat::test_that("reuse provenance separates creator and validator runtime", {
  generator_environment <- new.env(parent = globalenv())
  sys.source(
    file.path(module_dir, "generate_final_seurat_from_cellranger.R"),
    envir = generator_environment
  )
  output_root <- tempfile("figure7_upstream_runtime_")
  dir.create(output_root)
  definitions <- generator_environment$figure7_upstream_stage_definitions(
    output_root
  )
  definition <- definitions$final
  dir.create(definition$directory, recursive = TRUE)
  final_path <- generator_environment$figure7_upstream_stage_output(
    definition
  )
  writeLines("synthetic-final-cache", final_path, useBytes = TRUE)
  common <- c(
    seurat_upstream_environment_contract_sha256 =
      paste(rep("d", 64L), collapse = ""),
    source_config_contract_sha256 =
      paste(rep("e", 64L), collapse = ""),
    deterministic_workers = "1"
  )
  roles <- generator_environment$figure7_upstream_stage_dependency_roles(
    "final"
  )
  dependencies <- stats::setNames(
    vapply(seq_along(roles), function(index) {
      paste(rep(as.character(index), 64L), collapse = "")
    }, character(1L)),
    roles
  )
  generator_environment$figure7_upstream_write_manifest(
    "final",
    definition,
    common,
    dependencies
  )
  stage_manifest <- generator_environment$figure7_upstream_manifest_values(
    generator_environment$figure7_upstream_stage_manifest(definition)
  )
  creator_version <- stage_manifest[["creator_r_runtime_version"]]
  generator_environment$figure7_r_runtime_provenance <- function() {
    c(
      audit_r_runtime_version = "validator-runtime",
      audit_r_platform = "validator-platform",
      audit_r_blas = "validator-blas"
    )
  }
  generator_environment$figure7_upstream_write_run_manifest(
    output_root,
    final_path,
    common,
    requested_jobs = 8L,
    dependency_validation =
      generator_environment$figure7_upstream_dependency_validation_summary(
        "final",
        dependencies
      )
  )
  run_manifest <- generator_environment$figure7_upstream_manifest_values(
    file.path(output_root, "reconstruction_manifest.tsv")
  )
  testthat::expect_identical(
    run_manifest[["creator_r_runtime_version"]],
    creator_version
  )
  testthat::expect_identical(
    run_manifest[["validation_r_runtime_version"]],
    "validator-runtime"
  )
  testthat::expect_identical(run_manifest[["scientific_workers"]], "1")
  testthat::expect_identical(run_manifest[["requested_jobs"]], "8")
  testthat::expect_false(
    "audit_r_runtime_version" %in% names(run_manifest)
  )
})

testthat::test_that("selector reuses final or merged cache without H5 access", {
  generator <- new.env(parent = globalenv())
  sys.source(
    file.path(module_dir, "generate_final_seurat_from_cellranger.R"),
    envir = generator
  )
  config_path <- file.path(module_dir, "figure7_config.yaml")
  lock_path <- file.path(module_dir, "environment_lock.tsv")
  config <- figure7_read_config(config_path)
  all_ploidy <- file.path(
    repo_root,
    config$versioned_source_artifacts$endpoint_ploidy$default_path
  )
  sample_info <- file.path(
    repo_root,
    config$versioned_source_artifacts$sample_info$default_path
  )
  common <- generator$figure7_upstream_common_contract(
    lock_path,
    config,
    jobs = 1L
  )
  deposited <- tempfile("figure7_absent_deposited_", fileext = ".rds")
  stale_h5 <- tempfile("figure7_absent_cellranger_")

  final_root <- tempfile("figure7_final_only_")
  dir.create(final_root)
  final_definitions <- generator$figure7_upstream_stage_definitions(
    final_root
  )
  final_definition <- final_definitions$final
  dir.create(final_definition$directory, recursive = TRUE)
  final_path <- generator$figure7_upstream_stage_output(
    final_definition
  )
  writeLines("synthetic manifested final", final_path)
  final_roles <- generator$figure7_upstream_stage_dependency_roles(
    "final"
  )
  final_dependencies <- stats::setNames(
    rep(paste(rep("c", 64L), collapse = ""), length(final_roles)),
    final_roles
  )
  final_dependencies[["all_ploidy"]] <- figure7_sha256(all_ploidy)
  final_dependencies[["sample_info"]] <- figure7_sha256(sample_info)
  generator$figure7_upstream_write_manifest(
    "final",
    final_definition,
    common,
    final_dependencies
  )
  generator$figure7_upstream_write_run_manifest(
    final_root,
    final_path,
    common,
    requested_jobs = 1L,
    dependency_validation =
      generator$figure7_upstream_dependency_validation_summary(
        "final",
        c(
          all_ploidy = figure7_sha256(all_ploidy),
          sample_info = figure7_sha256(sample_info)
        )
      )
  )

  final_selection <- figure7_select_seurat_source(
    needs_seurat = TRUE,
    explicit_rds = "",
    deposited_rds = deposited,
    deposited_sha256 = as.character(config$raw_data$seurat_rds_sha256),
    upstream_dir = final_root,
    cellranger_root = stale_h5,
    module_dir = module_dir,
    environment_lock = lock_path,
    config = config,
    all_ploidy = all_ploidy,
    sample_info = sample_info
  )
  testthat::expect_identical(
    final_selection$source,
    "reused_reconstructed_rds"
  )
  testthat::expect_identical(
    final_selection$rds,
    normalizePath(final_path, mustWork = TRUE)
  )
  testthat::expect_false(file.exists(deposited))
  testthat::expect_false(dir.exists(stale_h5))

  explicit_selection <- figure7_select_seurat_source(
    needs_seurat = TRUE,
    explicit_rds = final_path,
    deposited_rds = deposited,
    deposited_sha256 = as.character(config$raw_data$seurat_rds_sha256),
    upstream_dir = final_root,
    cellranger_root = stale_h5,
    module_dir = module_dir,
    environment_lock = lock_path,
    config = config,
    all_ploidy = all_ploidy,
    sample_info = sample_info
  )
  testthat::expect_identical(
    explicit_selection$source,
    "explicit_reconstructed_rds"
  )
  equal_byte_hash <- figure7_sha256(final_path)
  equal_byte_reconstructed <- figure7_select_seurat_source(
    needs_seurat = TRUE,
    explicit_rds = final_path,
    deposited_rds = deposited,
    deposited_sha256 = equal_byte_hash,
    upstream_dir = final_root,
    cellranger_root = stale_h5,
    module_dir = module_dir,
    environment_lock = lock_path,
    config = config,
    all_ploidy = all_ploidy,
    sample_info = sample_info,
    required_source_kind = "manifested_reconstruction"
  )
  testthat::expect_identical(
    equal_byte_reconstructed$source,
    "explicit_reconstructed_rds"
  )
  equal_byte_auto <- figure7_select_seurat_source(
    needs_seurat = TRUE,
    explicit_rds = final_path,
    deposited_rds = deposited,
    deposited_sha256 = equal_byte_hash,
    upstream_dir = final_root,
    cellranger_root = stale_h5,
    module_dir = module_dir,
    environment_lock = lock_path,
    config = config,
    all_ploidy = all_ploidy,
    sample_info = sample_info
  )
  testthat::expect_identical(
    equal_byte_auto$source,
    "explicit_reconstructed_rds"
  )
  equal_byte_deposited <- figure7_select_seurat_source(
    needs_seurat = TRUE,
    explicit_rds = final_path,
    deposited_rds = deposited,
    deposited_sha256 = equal_byte_hash,
    upstream_dir = final_root,
    cellranger_root = stale_h5,
    module_dir = module_dir,
    environment_lock = lock_path,
    config = config,
    all_ploidy = all_ploidy,
    sample_info = sample_info,
    required_source_kind = "deposited"
  )
  testthat::expect_identical(
    equal_byte_deposited$source,
    "explicit_deposited_rds"
  )
  testthat::expect_length(
    equal_byte_deposited$transitive_dependencies,
    0L
  )

  merged_root <- tempfile("figure7_merged_only_")
  dir.create(merged_root)
  merged_definitions <- generator$figure7_upstream_stage_definitions(
    merged_root
  )
  merged_definition <- merged_definitions$merged
  dir.create(merged_definition$directory, recursive = TRUE)
  writeLines(
    "synthetic manifested merged",
    generator$figure7_upstream_stage_output(merged_definition)
  )
  merged_roles <- generator$figure7_upstream_stage_dependency_roles(
    "merged"
  )
  merged_dependencies <- stats::setNames(
    rep(paste(rep("d", 64L), collapse = ""), length(merged_roles)),
    merged_roles
  )
  merged_dependencies[["all_ploidy"]] <- figure7_sha256(all_ploidy)
  merged_dependencies[["sample_info"]] <- figure7_sha256(sample_info)
  generator$figure7_upstream_write_manifest(
    "merged",
    merged_definition,
    common,
    merged_dependencies
  )

  merged_selection <- figure7_select_seurat_source(
    needs_seurat = TRUE,
    explicit_rds = "",
    deposited_rds = deposited,
    deposited_sha256 = as.character(config$raw_data$seurat_rds_sha256),
    upstream_dir = merged_root,
    cellranger_root = stale_h5,
    module_dir = module_dir,
    environment_lock = lock_path,
    config = config,
    all_ploidy = all_ploidy,
    sample_info = sample_info
  )
  testthat::expect_identical(merged_selection$action, "build_upstream")
  testthat::expect_identical(merged_selection$resumable_stage, "merged")
  testthat::expect_identical(
    merged_selection$source,
    "planned_upstream_resume_merged"
  )
  testthat::expect_false(file.exists(deposited))
})

testthat::test_that("trimmed upstream source excludes non-figure diagnostics", {
  source_text <- c(
    readLines(
      file.path(module_dir, "src", "seurat_upstream.R"),
      warn = FALSE
    ),
    readLines(
      file.path(
        module_dir,
        "generate_final_seurat_from_cellranger.R"
      ),
      warn = FALSE
    )
  )
  prohibited_calls <- c(
    "FindMarkers\\(",
    "FindAllMarkers\\(",
    "ggsave\\(",
    "DimPlot\\(",
    "FeaturePlot\\(",
    "VlnPlot\\(",
    "DotPlot\\(",
    "pdf\\("
  )
  for (pattern in prohibited_calls) {
    testthat::expect_false(any(grepl(pattern, source_text)))
  }
  testthat::expect_false(any(grepl("/Volumes/", source_text, fixed = TRUE)))
})

testthat::test_that("standalone cluster completion chain is reusable downstream", {
  root <- tempfile("standalone_cluster_")
  provenance <- file.path(root, "00_provenance")
  object_dir <- file.path(root, "03_final_cluster", "03_objects")
  h5_root <- file.path(root, "raw", "SUM-159", "A02_cellRanger")
  dir.create(provenance, recursive = TRUE)
  dir.create(object_dir, recursive = TRUE)
  h5_paths <- vapply(sprintf("sample%02d", seq_len(18L)), function(sample) {
    path <- file.path(
      h5_root,
      paste0(sample, "-Count-HM"),
      "outs",
      paste0(sample, "-Count-HM_filtered_feature_bc_matrix.h5")
    )
    dir.create(dirname(path), recursive = TRUE)
    writeBin(charToRaw(paste0("h5-", sample)), path)
    path
  }, character(1))
  all_ploidy <- file.path(root, "all_ploidy.tsv")
  sample_info <- file.path(root, "sample_info.xlsx")
  writeLines("cell\tploidy", all_ploidy)
  writeBin(charToRaw("xlsx fixture"), sample_info)
  inputs <- c(h5_paths, all_ploidy, sample_info)
  input_types <- c(
    rep("filtered_feature_bc_matrix.h5", 18L),
    "all_ploidy.tsv",
    "sample_info.xlsx"
  )
  run_manifest <- data.frame(
    input_type = input_types,
    path = normalizePath(inputs, mustWork = TRUE),
    size_bytes = as.character(file.info(inputs)$size),
    md5 = unname(tools::md5sum(inputs)),
    modified_time = rep("2026-08-05T00:00:00-0400", 20L),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    run_manifest,
    file.path(provenance, "run_manifest.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  final_rds <- file.path(
    object_dir,
    "integrated_sct_cca_seurat_final_reclustered.rds"
  )
  writeBin(charToRaw("standalone rds fixture"), final_rds)
  runtime <- c(
    status = "PASS",
    phase = "post_filter",
    active_sif_md5 = paste(rep("a", 32L), collapse = ""),
    R = "4.5.0",
    Seurat = "5.3.0",
    private_r_library = "",
    final_object = normalizePath(final_rds, mustWork = TRUE),
    final_object_size_bytes = as.character(file.info(final_rds)$size),
    final_object_md5 = unname(tools::md5sum(final_rds))
  )
  utils::write.table(
    data.frame(field = names(runtime), value = unname(runtime)),
    file.path(provenance, "final_artifact_runtime.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  writeLines(
    c(
      "status=PASS",
      "final_object=03_final_cluster/03_objects/integrated_sct_cca_seurat_final_reclustered.rds"
    ),
    file.path(provenance, "PIPELINE_COMPLETE.txt")
  )
  audit_dir <- file.path(
    root,
    "00_validation",
    "rds_semantic_audit"
  )
  dir.create(audit_dir, recursive = TRUE)
  reference_rds <- file.path(root, "zenodo_reference.rds")
  writeBin(charToRaw("different Zenodo reference serialization"), reference_rds)
  audit_summary <- data.frame(
    group = "object",
    check = "synthetic_semantic_contract",
    status = "PASS",
    observed = "equivalent",
    reference = "equivalent",
    threshold = "fixture",
    details = "standalone selection fixture",
    stringsAsFactors = FALSE
  )
  audit_summary_path <- file.path(audit_dir, "audit_summary.tsv")
  utils::write.table(
    audit_summary,
    audit_summary_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  identities <- data.frame(
    artifact = c("generated", "zenodo_reference"),
    path = normalizePath(c(final_rds, reference_rds), mustWork = TRUE),
    size_bytes = as.character(file.info(c(final_rds, reference_rds))$size),
    md5 = unname(tools::md5sum(c(final_rds, reference_rds))),
    sha256 = vapply(
      c(final_rds, reference_rds),
      figure7_sha256,
      character(1L)
    ),
    role = c("candidate_downstream_input", "semantic_reference"),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    identities,
    file.path(audit_dir, "rds_file_identity.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  audit_report_names <- c(
    "audit_summary.tsv", "metadata_comparison.tsv",
    "cluster_comparison.tsv", "graph_comparison.tsv",
    "assay_numeric_comparison.tsv", "pca_comparison.tsv",
    "umap_displacement_summary.tsv", "umap_largest_displacements.tsv",
    "command_comparison.tsv", "rds_file_identity.tsv"
  )
  missing_report_names <- setdiff(
    audit_report_names,
    c("audit_summary.tsv", "rds_file_identity.tsv")
  )
  for (name in missing_report_names) {
    writeLines("check\tstatus\nsynthetic\tPASS", file.path(audit_dir, name))
  }
  audit_report_hash_lines <- vapply(audit_report_names, function(name) {
    paste0(
      "report_sha256.", name, "=",
      figure7_sha256(file.path(audit_dir, name))
    )
  }, character(1L))
  writeLines(
    c(
      "schema_version=semantic_rds_audit_v1",
      "status=PASS",
      "checks_total=1",
      "checks_passed=1",
      "checks_failed=0",
      paste0("generated_rds=", normalizePath(final_rds)),
      paste0("generated_rds_size_bytes=", file.info(final_rds)$size),
      paste0("generated_rds_md5=", unname(tools::md5sum(final_rds))),
      paste0("generated_rds_sha256=", figure7_sha256(final_rds)),
      paste0("zenodo_reference_rds=", normalizePath(reference_rds)),
      paste0("zenodo_reference_rds_size_bytes=", file.info(reference_rds)$size),
      paste0("zenodo_reference_rds_md5=", unname(tools::md5sum(reference_rds))),
      paste0("zenodo_reference_rds_sha256=", figure7_sha256(reference_rds)),
      paste0("downstream_rds=", normalizePath(final_rds)),
      paste0("audit_summary=", normalizePath(audit_summary_path)),
      paste0("audit_summary_sha256=", figure7_sha256(audit_summary_path)),
      audit_report_hash_lines,
      "validated_at=2026-08-05T00:00:00-0400"
    ),
    file.path(audit_dir, "AUDIT_COMPLETE.txt")
  )

  validation <- figure7_validate_standalone_seurat_artifact(
    output_root = root,
    module_dir = module_dir,
    all_ploidy = all_ploidy,
    sample_info = sample_info,
    expected_rds = final_rds,
    cellranger_root = h5_root
  )
  testthat::expect_true(validation$valid)
  testthat::expect_identical(validation$rds, normalizePath(final_rds))
  testthat::expect_length(validation$dependencies, 21L)
  testthat::expect_identical(
    unname(validation$dependencies[["semantic_audit_reference_rds"]]),
    figure7_sha256(reference_rds)
  )
  testthat::expect_identical(
    figure7_infer_seurat_upstream_root(final_rds),
    normalizePath(root)
  )

  invalid_runtime <- runtime
  invalid_runtime[["Seurat"]] <- ""
  utils::write.table(
    data.frame(field = names(invalid_runtime), value = unname(invalid_runtime)),
    file.path(provenance, "final_artifact_runtime.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  testthat::expect_error(
    figure7_validate_standalone_seurat_artifact(
      output_root = root,
      module_dir = module_dir,
      all_ploidy = all_ploidy,
      sample_info = sample_info,
      expected_rds = final_rds,
      cellranger_root = h5_root
    ),
    "invalid field/value schema"
  )
  utils::write.table(
    data.frame(field = names(runtime), value = unname(runtime)),
    file.path(provenance, "final_artifact_runtime.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  writeLines("tampered", audit_summary_path)
  testthat::expect_error(
    figure7_validate_standalone_seurat_artifact(
      output_root = root,
      module_dir = module_dir,
      all_ploidy = all_ploidy,
      sample_info = sample_info,
      expected_rds = final_rds,
      cellranger_root = h5_root
    ),
    "audit summary is inconsistent"
  )
  utils::write.table(
    audit_summary,
    audit_summary_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  writeLines("tampered", h5_paths[[1L]])
  testthat::expect_error(
    figure7_validate_standalone_seurat_artifact(
      output_root = root,
      module_dir = module_dir,
      all_ploidy = all_ploidy,
      sample_info = sample_info,
      expected_rds = final_rds,
      cellranger_root = h5_root
    ),
    "size/MD5 verification failed"
  )
})

testthat::test_that("standalone completion contract ignores only the run timestamp", {
  first <- c(
    status = "PASS",
    completed_at = "2026-08-05T15:00:00-0400",
    final_object = paste0(
      "03_final_cluster/03_objects/",
      "integrated_sct_cca_seurat_final_reclustered.rds"
    ),
    removed_clusters_selected_dynamically = "3,4,9,9c"
  )
  second <- first
  second[["completed_at"]] <- "2026-08-05T16:00:00-0400"
  changed <- second
  changed[["removed_clusters_selected_dynamically"]] <- "3,4,9"

  testthat::expect_identical(
    figure7_standalone_completion_contract_sha256(first),
    figure7_standalone_completion_contract_sha256(second)
  )
  testthat::expect_false(identical(
    figure7_standalone_completion_contract_sha256(first),
    figure7_standalone_completion_contract_sha256(changed)
  ))
})
