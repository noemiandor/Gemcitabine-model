#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
)

resolve_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(hit) != 1L) stop("Cannot resolve script directory.", call. = FALSE)
  dirname(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = TRUE))
}

script_dir <- resolve_script_dir()
source(file.path(script_dir, "Utils.R"))
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")
required_packages <- c("yaml", "Seurat", "future")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
input_rds <- file.path(results_root, "05b_qc", "mouse_seurat_list_post_qc.rds")
out_dir <- file.path(results_root, "05c_cca_integration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds)) stop("Missing 05b input: ", input_rds, call. = FALSE)

icfg <- cfg$integration
nfeatures <- as.integer(icfg$nfeatures)
reference_samples <- unlist(icfg$reference_samples, use.names = FALSE)
future_strategy <- as.character(icfg$future_strategy)
future_workers <- as.integer(icfg$future_workers)
future_globals_maxsize_gb <- as.numeric(icfg$future_globals_maxsize_gb)
seed <- as.integer(icfg$seed)

seurat_list <- readRDS(input_rds)
if (!is.list(seurat_list) || length(seurat_list) < 2L) stop("CCA requires at least two non-empty sample objects.", call. = FALSE)
if (any(vapply(seurat_list, function(x) any(startsWith(rownames(x), cfg$species$human_prefix)), logical(1)))) {
  stop("Human-prefixed features are present before CCA.", call. = FALSE)
}
if (length(reference_samples) < 1L || anyDuplicated(reference_samples)) {
  stop("integration.reference_samples must contain one or more unique sample names.", call. = FALSE)
}
reference_indices <- match(reference_samples, names(seurat_list))
if (anyNA(reference_indices)) {
  stop(
    "Configured CCA reference samples are absent after QC: ",
    paste(reference_samples[is.na(reference_indices)], collapse = ", "),
    call. = FALSE
  )
}

future_plan <- configure_future_for_seurat(
  max_size_gb = future_globals_maxsize_gb,
  strategy = future_strategy,
  workers = future_workers
)
set.seed(seed)
prepared_checkpoint <- file.path(out_dir, "checkpoint_sct_prepared.rds")
if (file.exists(prepared_checkpoint)) {
  message("[05c] Reusing prepared-SCT checkpoint: ", prepared_checkpoint)
  checkpoint <- readRDS(prepared_checkpoint)
  expected_cell_counts <- vapply(seurat_list, ncol, integer(1))
  if (
    !is.list(checkpoint) ||
    !identical(checkpoint$sample_names, names(seurat_list)) ||
    !identical(checkpoint$cell_counts, expected_cell_counts) ||
    !identical(as.integer(checkpoint$nfeatures), nfeatures)
  ) {
    stop("Prepared-SCT checkpoint does not match the current 05b input/configuration.", call. = FALSE)
  }
  seurat_list <- checkpoint$seurat_list
  integration_features <- checkpoint$integration_features
  rm(checkpoint)
  invisible(gc())
} else {
  message("[05c] Running SCTransform for ", length(seurat_list), " samples.")
  seurat_list <- lapply(seurat_list, function(obj) {
    Seurat::SCTransform(
      object = obj,
      assay = "RNA",
      new.assay.name = "SCT",
      vars.to.regress = "percent.mt",
      return.only.var.genes = FALSE,
      verbose = FALSE
    )
  })

  message("[05c] Selecting ", nfeatures, " features and preparing SCT integration.")
  integration_features <- Seurat::SelectIntegrationFeatures(object.list = seurat_list, nfeatures = nfeatures)
  seurat_list <- Seurat::PrepSCTIntegration(
    object.list = seurat_list,
    anchor.features = integration_features,
    verbose = FALSE
  )
  checkpoint <- list(
    sample_names = names(seurat_list),
    cell_counts = vapply(seurat_list, ncol, integer(1)),
    nfeatures = nfeatures,
    integration_features = integration_features,
    seurat_list = seurat_list
  )
  saveRDS(checkpoint, prepared_checkpoint, compress = FALSE)
  rm(checkpoint)
  invisible(gc())
}

set.seed(seed)
message(
  "[05c] Finding reference-based SCT-CCA anchors using: ",
  paste(reference_samples, collapse = ", ")
)
minimum_sample_cells <- min(vapply(seurat_list, ncol, integer(1)))
if (minimum_sample_cells < 6L) {
  stop("Every sample requires at least six post-QC mouse cells for SCT-CCA integration.", call. = FALSE)
}
integration_k_weight <- min(100L, minimum_sample_cells - 1L)
anchors_checkpoint <- file.path(out_dir, "checkpoint_reference_cca_anchors.rds")
if (file.exists(anchors_checkpoint)) {
  message("[05c] Reusing reference-CCA anchors checkpoint: ", anchors_checkpoint)
  anchor_checkpoint <- readRDS(anchors_checkpoint)
  if (
    !is.list(anchor_checkpoint) ||
    !identical(anchor_checkpoint$sample_names, names(seurat_list)) ||
    !identical(anchor_checkpoint$reference_samples, reference_samples) ||
    !identical(as.integer(anchor_checkpoint$nfeatures), nfeatures)
  ) {
    stop("Reference-CCA anchors checkpoint does not match the current input/configuration.", call. = FALSE)
  }
  anchors <- anchor_checkpoint$anchors
  rm(anchor_checkpoint)
  invisible(gc())
} else {
  anchors <- Seurat::FindIntegrationAnchors(
    object.list = seurat_list,
    normalization.method = "SCT",
    anchor.features = integration_features,
    reduction = "cca",
    reference = reference_indices,
    verbose = FALSE
  )
  anchor_checkpoint <- list(
    sample_names = names(seurat_list),
    reference_samples = reference_samples,
    nfeatures = nfeatures,
    anchors = anchors
  )
  saveRDS(anchor_checkpoint, anchors_checkpoint, compress = FALSE)
  rm(anchor_checkpoint)
  invisible(gc())
}
integrated <- Seurat::IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  k.weight = integration_k_weight,
  verbose = FALSE
)
Seurat::DefaultAssay(integrated) <- "integrated"
integrated@misc$mouse05_integration_contract <- list(
  normalization_method = "SCT",
  reduction = "cca",
  integration_mode = "reference_based",
  reference_samples = reference_samples,
  integration_features = nfeatures,
  k_weight = integration_k_weight,
  seed = seed,
  samples = names(seurat_list),
  human_features_retained_in_RNA = sum(startsWith(rownames(integrated[["RNA"]]), cfg$species$human_prefix))
)
if (integrated@misc$mouse05_integration_contract$human_features_retained_in_RNA != 0L) {
  stop("Human features were detected after CCA integration.", call. = FALSE)
}

saveRDS(seurat_list, file.path(out_dir, "mouse_seurat_list_post_sct.rds"), compress = FALSE)
saveRDS(integrated, file.path(out_dir, "mouse_integrated_sct_cca.rds"), compress = FALSE)
utils::write.csv(
  data.frame(
    parameter = c("samples", "cells", "nfeatures", "normalization_method", "reduction", "integration_mode", "reference_samples", "minimum_sample_cells", "k_weight", "future_strategy", "future_workers", "future_globals_maxsize_gb", "seed"),
    value = c(length(seurat_list), ncol(integrated), nfeatures, "SCT", "cca", "reference_based", paste(reference_samples, collapse = ";"), minimum_sample_cells, integration_k_weight, future_plan$strategy, future_plan$workers, future_globals_maxsize_gb, seed)
  ),
  file.path(out_dir, "integration_parameters.csv"), row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05c.txt"))
writeLines(
  c("05c SCT-CCA integration completed.", paste0("Samples: ", length(seurat_list)), paste0("Cells: ", ncol(integrated))),
  file.path(out_dir, "completion.txt")
)
message("[05c] Completed. Output: ", out_dir)
