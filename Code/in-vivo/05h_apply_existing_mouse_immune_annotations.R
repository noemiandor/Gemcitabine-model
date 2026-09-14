#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
)

resolve_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(hit) != 1L) stop("Cannot resolve script directory.", call. = FALSE)
  dirname(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = TRUE))
}

cluster_order <- function(x) {
  x <- unique(as.character(x))
  numeric_x <- suppressWarnings(as.numeric(x))
  if (all(is.finite(numeric_x))) x[order(numeric_x)] else sort(x)
}

parse_marker_list <- function(x) {
  x <- trimws(as.character(x))
  if (length(x) != 1L || is.na(x) || !nzchar(x)) return(character())
  genes <- trimws(unlist(strsplit(x, "[;,|]", perl = TRUE), use.names = FALSE))
  unique(genes[nzchar(genes)])
}

validate_support_markers <- function(marker_string, candidate, resource, min_n, max_n, evidence_class, cluster_id) {
  genes <- parse_marker_list(marker_string)
  if (length(genes) < min_n || length(genes) > max_n) {
    stop(
      "Cluster ", cluster_id, " ", evidence_class, " support requires ", min_n, "-", max_n,
      " markers; observed ", length(genes), ".", call. = FALSE
    )
  }
  allowed <- resource$marker_gene[resource$candidate_cell_type == candidate]
  bad <- genes[!toupper(genes) %in% toupper(allowed)]
  if (length(bad) > 0L) {
    stop(
      "Cluster ", cluster_id, " ", evidence_class, " support markers are absent from the approved ",
      candidate, " marker set: ", paste(bad, collapse = ", "), call. = FALSE
    )
  }
  genes
}

collapse_unique <- function(x) {
  x <- trimws(as.character(x))
  x <- unique(x[!is.na(x) & nzchar(x) & x != "NA"])
  paste(x, collapse = "; ")
}

candidate_reference_evidence <- function(candidate, support_markers, resource) {
  if (!nzchar(candidate) || length(support_markers) == 0L) {
    return(list(titles = "", pmids = "", dois = "", urls = ""))
  }
  z <- resource[
    as.character(resource$candidate_cell_type) == candidate &
      toupper(resource$marker_gene) %in% toupper(support_markers),
    , drop = FALSE
  ]
  z <- z[match(toupper(support_markers), toupper(z$marker_gene)), , drop = FALSE]
  pmids <- trimws(as.character(z$pmid))
  pmids <- sub("[.]0$", "", pmids)
  pmids <- unique(pmids[!is.na(pmids) & nzchar(pmids) & pmids != "NA"])
  dois <- trimws(as.character(z$doi))
  dois <- unique(dois[!is.na(dois) & nzchar(dois) & dois != "NA"])
  urls <- c(
    if (length(pmids) > 0L) paste0("https://pubmed.ncbi.nlm.nih.gov/", pmids, "/") else character(),
    if (length(dois) > 0L) paste0("https://doi.org/", dois) else character()
  )
  list(
    titles = collapse_unique(z$reference),
    pmids = paste(pmids, collapse = ";"),
    dois = paste(dois, collapse = ";"),
    urls = paste(unique(urls), collapse = ";")
  )
}

cluster_marker_evidence <- function(cluster_id, candidate, support_markers, review_expression) {
  if (!nzchar(candidate) || length(support_markers) == 0L) {
    return(list(expression = "", positive_degs = "", all_markers_observed = TRUE))
  }
  z <- review_expression[
    review_expression$cluster == cluster_id &
      as.character(review_expression$candidate_cell_type) == candidate &
      toupper(review_expression$marker_gene) %in% toupper(support_markers),
    , drop = FALSE
  ]
  z <- z[match(toupper(support_markers), toupper(z$marker_gene)), , drop = FALSE]
  is_positive <- as.logical(z$is_positive_cluster_deg)
  is_positive[is.na(is_positive)] <- FALSE
  list(
    expression = paste0(
      z$marker_gene, "=", sprintf("%.1f%%", z$percent_cells_expressing),
      "[cluster_DEG=", ifelse(is_positive, "yes", "no"), "]",
      collapse = ";"
    ),
    positive_degs = paste(z$marker_gene[is_positive], collapse = ";"),
    all_markers_observed = nrow(z) == length(support_markers) && all(z$percent_cells_expressing > 0)
  )
}

write_plot_pair <- function(plot, stem, out_dir, width, height) {
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = width, height = height, limitsize = FALSE)
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300, limitsize = FALSE)
}

script_dir <- resolve_script_dir()
script_path <- file.path(script_dir, "05h_apply_existing_mouse_immune_annotations.R")
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")
required_packages <- c("yaml", "Seurat", "ggplot2", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
input_rds <- file.path(results_root, "05f_broad_celltype_annotation", "mouse_clustered_broad_celltype_annotated.rds")
review_dir <- file.path(results_root, as.character(cfg$immune_fine_annotation$existing_review_output_directory))
review_completion <- file.path(review_dir, "completion.txt")
review_manifest <- file.path(review_dir, "output_manifest_sha256.csv")
review_expression_file <- file.path(review_dir, "existing_immune_cluster_marker_expression.csv")
mapping_file <- file.path(script_dir, as.character(cfg$immune_fine_annotation$approved_existing_cluster_mapping))
resource_file <- file.path(script_dir, as.character(cfg$annotation$immune_marker_resource))
out_dir <- file.path(results_root, as.character(cfg$immune_fine_annotation$existing_annotation_output_directory))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  input_rds, review_completion, review_manifest, review_expression_file,
  mapping_file, resource_file
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required input(s): ", paste(missing_inputs, collapse = ", "),
    ". 05h must not run until the 05g manual-review template has been completed, approved, and installed as the configured mapping file.",
    call. = FALSE
  )
}

completion_lines <- readLines(review_completion, warn = FALSE)
required_review_contract <- c(
  "Annotation level under review: existing 05d cluster only.",
  "No automated immune lineage or state label was calculated.",
  "No fine-annotation label was added to Seurat metadata.",
  "Manual review is required before 05h can apply any cluster label.",
  "05h and later steps were not invoked."
)
if (!all(required_review_contract %in% completion_lines)) {
  stop("05g completion file does not satisfy the review-only/manual-gate contract.", call. = FALSE)
}

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("05f input is not a Seurat object.", call. = FALSE)
required_meta <- c(
  "sample", "barcode_raw", "seurat_clusters", "mouse_fraction", "scDblFinder.class",
  "mouse_broad_cell_type", "mouse_broad_annotation_level", "mouse_broad_review_status"
)
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0L) stop("05f object lacks metadata: ", paste(missing_meta, collapse = ", "), call. = FALSE)
if (!"RNA" %in% names(obj@assays)) stop("05f object lacks RNA.", call. = FALSE)
if (!"umap" %in% names(obj@reductions)) stop("05f object lacks UMAP.", call. = FALSE)
human_feature_n <- sum(startsWith(rownames(obj[["RNA"]]), as.character(cfg$species$human_prefix)))
if (human_feature_n != 0L) stop("Human features detected before existing-cluster immune annotation: ", human_feature_n, call. = FALSE)

cluster_vector <- as.character(obj$seurat_clusters)
broad_vector <- as.character(obj$mouse_broad_cell_type)
confirmed_types <- as.character(unlist(cfg$immune_fine_annotation$confirmed_broad_cell_types, use.names = FALSE))
candidate_types <- as.character(unlist(cfg$immune_fine_annotation$candidate_broad_cell_types, use.names = FALSE))
review_mask <- broad_vector %in% c(confirmed_types, candidate_types)
review_clusters <- cluster_order(cluster_vector[review_mask])
review_cells <- colnames(obj)[review_mask]
cluster_to_broad <- tapply(broad_vector[review_mask], cluster_vector[review_mask], unique)
cluster_counts <- table(factor(cluster_vector[review_mask], levels = review_clusters))
if (any(as.character(obj$scDblFinder.class[review_mask]) != "singlet")) stop("The reviewed immune universe contains non-singlets.", call. = FALSE)
if (any(as.numeric(obj$mouse_fraction[review_mask]) < as.numeric(cfg$species$mouse_fraction_min), na.rm = TRUE)) {
  stop("The reviewed immune universe violates the mouse-fraction threshold.", call. = FALSE)
}

resource <- utils::read.csv(resource_file, stringsAsFactors = FALSE, check.names = FALSE, na.strings = character())
required_resource <- c("annotation_level", "parent_cell_type", "candidate_cell_type", "marker_gene", "marker_role")
missing_resource <- setdiff(required_resource, colnames(resource))
if (length(missing_resource) > 0L) stop("Immune marker resource lacks: ", paste(missing_resource, collapse = ", "), call. = FALSE)
if (any(tolower(resource$marker_role) != "positive")) stop("Immune marker resource contains non-positive markers.", call. = FALSE)
lineage_candidates <- unique(resource$candidate_cell_type[resource$annotation_level == "lineage"])
state_resource <- resource[resource$annotation_level == "state", , drop = FALSE]
state_candidates <- unique(state_resource$candidate_cell_type)
allowed_lineages <- c(lineage_candidates, "Ambiguous", "Unresolved")
allowed_states <- c(state_candidates, "Mixed", "Unresolved")

mapping <- utils::read.csv(mapping_file, stringsAsFactors = FALSE, check.names = FALSE, na.strings = character())
required_mapping <- c(
  "cluster", "n_cells", "source_broad_cell_type", "review_scope", "top_10_positive_degs",
  "reviewer_immune_lineage", "reviewer_lineage_confidence", "reviewer_immune_state",
  "reviewer_state_confidence", "lineage_support_markers", "state_support_markers",
  "reviewer_notes", "review_status"
)
missing_mapping <- setdiff(required_mapping, colnames(mapping))
if (length(missing_mapping) > 0L) stop("Approved mapping lacks columns: ", paste(missing_mapping, collapse = ", "), call. = FALSE)
mapping_text_columns <- c(
  "source_broad_cell_type", "review_scope", "top_10_positive_degs",
  "reviewer_immune_lineage", "reviewer_lineage_confidence", "reviewer_immune_state",
  "reviewer_state_confidence", "lineage_support_markers", "state_support_markers",
  "reviewer_notes", "review_status"
)
for (column in mapping_text_columns) {
  mapping[[column]] <- as.character(mapping[[column]])
  mapping[[column]][is.na(mapping[[column]])] <- ""
}
mapping$cluster <- as.character(mapping$cluster)
if (anyDuplicated(mapping$cluster)) stop("Approved mapping contains duplicate clusters.", call. = FALSE)
if (!identical(cluster_order(mapping$cluster), review_clusters)) {
  stop(
    "Approved mapping and review cluster sets differ. Missing mapping: ", paste(setdiff(review_clusters, mapping$cluster), collapse = ","),
    "; extra mapping: ", paste(setdiff(mapping$cluster, review_clusters), collapse = ","), call. = FALSE
  )
}
mapping <- mapping[match(review_clusters, mapping$cluster), , drop = FALSE]
if (any(mapping$review_status != "manually_approved")) stop("Every mapping row must be manually_approved.", call. = FALSE)
if (any(!mapping$reviewer_lineage_confidence %in% c("high", "medium", "low"))) stop("Invalid lineage confidence value.", call. = FALSE)
if (any(!mapping$reviewer_immune_lineage %in% allowed_lineages)) {
  stop("Invalid immune lineage(s): ", paste(setdiff(unique(mapping$reviewer_immune_lineage), allowed_lineages), collapse = ", "), call. = FALSE)
}
if (any(as.integer(mapping$n_cells) != as.integer(cluster_counts))) stop("Approved mapping n_cells differs from the 05f object.", call. = FALSE)
expected_broad <- unname(cluster_to_broad[mapping$cluster])
if (!identical(as.character(mapping$source_broad_cell_type), as.character(expected_broad))) {
  stop("Approved mapping source_broad_cell_type differs from the 05f object.", call. = FALSE)
}

min_support <- as.integer(cfg$immune_fine_annotation$min_support_markers)
max_support <- as.integer(cfg$immune_fine_annotation$max_support_markers)
if (!is.finite(min_support) || !is.finite(max_support) || min_support < 1L || max_support < min_support) {
  stop("Invalid support-marker limits in immune_fine_annotation.", call. = FALSE)
}

mapping$lineage_support_markers_normalized <- ""
mapping$state_support_markers_normalized <- ""
for (i in seq_len(nrow(mapping))) {
  cluster_id <- mapping$cluster[[i]]
  lineage <- mapping$reviewer_immune_lineage[[i]]
  state <- trimws(mapping$reviewer_immune_state[[i]])
  if (!lineage %in% c("Ambiguous", "Unresolved")) {
    lineage_genes <- validate_support_markers(
      mapping$lineage_support_markers[[i]], lineage, resource,
      min_support, max_support, "lineage", cluster_id
    )
    mapping$lineage_support_markers_normalized[[i]] <- paste(lineage_genes, collapse = ";")
  } else {
    optional_lineage <- parse_marker_list(mapping$lineage_support_markers[[i]])
    if (length(optional_lineage) > max_support) stop("Cluster ", cluster_id, " has too many optional lineage support markers.", call. = FALSE)
    mapping$lineage_support_markers_normalized[[i]] <- paste(optional_lineage, collapse = ";")
  }

  if (!nzchar(state)) {
    if (nzchar(trimws(mapping$reviewer_state_confidence[[i]])) || length(parse_marker_list(mapping$state_support_markers[[i]])) > 0L) {
      stop("Cluster ", cluster_id, " has state confidence/support markers but no immune state.", call. = FALSE)
    }
    next
  }
  if (!state %in% allowed_states) stop("Cluster ", cluster_id, " has invalid immune state: ", state, call. = FALSE)
  if (!mapping$reviewer_state_confidence[[i]] %in% c("high", "medium", "low", "unresolved")) {
    stop("Cluster ", cluster_id, " has invalid state confidence.", call. = FALSE)
  }
  if (state %in% state_candidates) {
    expected_parent <- unique(state_resource$parent_cell_type[state_resource$candidate_cell_type == state])
    if (length(expected_parent) != 1L || lineage != expected_parent) {
      stop(
        "Cluster ", cluster_id, " state ", state, " requires parent lineage ", paste(expected_parent, collapse = "/"),
        ", observed ", lineage, ".", call. = FALSE
      )
    }
    state_genes <- validate_support_markers(
      mapping$state_support_markers[[i]], state, resource,
      min_support, max_support, "state", cluster_id
    )
    mapping$state_support_markers_normalized[[i]] <- paste(state_genes, collapse = ";")
  } else {
    optional_state <- parse_marker_list(mapping$state_support_markers[[i]])
    if (length(optional_state) > max_support) stop("Cluster ", cluster_id, " has too many optional state support markers.", call. = FALSE)
    mapping$state_support_markers_normalized[[i]] <- paste(optional_state, collapse = ";")
  }
}

review_expression <- utils::read.csv(review_expression_file, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", ""))
required_expression <- c("cluster", "candidate_cell_type", "marker_gene", "percent_cells_expressing", "is_positive_cluster_deg")
missing_expression <- setdiff(required_expression, colnames(review_expression))
if (length(missing_expression) > 0L) stop("05g marker-expression table lacks: ", paste(missing_expression, collapse = ", "), call. = FALSE)
review_expression$cluster <- as.character(review_expression$cluster)
for (i in seq_len(nrow(mapping))) {
  cluster_id <- mapping$cluster[[i]]
  lineage <- mapping$reviewer_immune_lineage[[i]]
  if (!lineage %in% c("Ambiguous", "Unresolved")) {
    support <- parse_marker_list(mapping$lineage_support_markers_normalized[[i]])
    z <- review_expression[
      review_expression$cluster == cluster_id &
        as.character(review_expression$candidate_cell_type) == lineage &
        toupper(review_expression$marker_gene) %in% toupper(support),
      , drop = FALSE
    ]
    if (nrow(z) != length(support) || any(z$percent_cells_expressing <= 0, na.rm = TRUE)) {
      stop("Cluster ", cluster_id, " lineage support markers lack observed RNA expression.", call. = FALSE)
    }
  }
  state <- trimws(mapping$reviewer_immune_state[[i]])
  if (state %in% state_candidates) {
    support <- parse_marker_list(mapping$state_support_markers_normalized[[i]])
    z <- review_expression[
      review_expression$cluster == cluster_id &
        as.character(review_expression$candidate_cell_type) == state &
        toupper(review_expression$marker_gene) %in% toupper(support),
      , drop = FALSE
    ]
    if (nrow(z) != length(support) || any(z$percent_cells_expressing <= 0, na.rm = TRUE)) {
      stop("Cluster ", cluster_id, " state support markers lack observed RNA expression.", call. = FALSE)
    }
  }
}

annotation_evidence <- do.call(rbind, lapply(seq_len(nrow(mapping)), function(i) {
  cluster_id <- mapping$cluster[[i]]
  lineage <- mapping$reviewer_immune_lineage[[i]]
  state <- trimws(mapping$reviewer_immune_state[[i]])
  lineage_markers <- parse_marker_list(mapping$lineage_support_markers_normalized[[i]])
  state_markers <- parse_marker_list(mapping$state_support_markers_normalized[[i]])
  lineage_expression <- cluster_marker_evidence(cluster_id, lineage, lineage_markers, review_expression)
  state_expression <- cluster_marker_evidence(cluster_id, state, state_markers, review_expression)
  lineage_references <- candidate_reference_evidence(lineage, lineage_markers, resource)
  state_references <- candidate_reference_evidence(state, state_markers, resource)

  traceability_status <- if (lineage %in% c("Ambiguous", "Unresolved")) {
    "data_traceable_no_specific_lineage_claim"
  } else if (state == "Unresolved") {
    "lineage_literature_traceable_state_data_unresolved"
  } else if (!nzchar(state)) {
    "lineage_literature_traceable_no_state_claim"
  } else {
    "lineage_and_state_literature_traceable"
  }
  annotation_basis <- if (lineage %in% c("Ambiguous", "Unresolved")) {
    "Cluster DEGs, marker-expression review, UMAP context, and manual decision; no specific lineage claim"
  } else if (state == "Unresolved") {
    "Literature-backed lineage marker panel plus cluster DEGs; state intentionally unresolved"
  } else if (!nzchar(state)) {
    "Literature-backed lineage marker panel plus cluster DEGs; no separate state claim"
  } else {
    "Literature-backed lineage and state marker panels, marker expression, cluster DEGs, and manual review"
  }

  data.frame(
    annotation_level = "existing_05d_cluster",
    cluster = cluster_id,
    n_cells = as.integer(mapping$n_cells[[i]]),
    source_broad_cell_type = mapping$source_broad_cell_type[[i]],
    approved_immune_lineage = lineage,
    lineage_confidence = mapping$reviewer_lineage_confidence[[i]],
    approved_immune_state = state,
    state_confidence = mapping$reviewer_state_confidence[[i]],
    top_10_positive_degs = mapping$top_10_positive_degs[[i]],
    lineage_support_markers = mapping$lineage_support_markers_normalized[[i]],
    lineage_marker_expression_and_deg_status = lineage_expression$expression,
    lineage_support_positive_cluster_degs = lineage_expression$positive_degs,
    state_support_markers = mapping$state_support_markers_normalized[[i]],
    state_marker_expression_and_deg_status = state_expression$expression,
    state_support_positive_cluster_degs = state_expression$positive_degs,
    annotation_reason = mapping$reviewer_notes[[i]],
    annotation_basis = annotation_basis,
    traceability_status = traceability_status,
    lineage_reference_titles = lineage_references$titles,
    lineage_reference_pmids = lineage_references$pmids,
    lineage_reference_dois = lineage_references$dois,
    lineage_reference_urls = lineage_references$urls,
    state_reference_titles = state_references$titles,
    state_reference_pmids = state_references$pmids,
    state_reference_dois = state_references$dois,
    state_reference_urls = state_references$urls,
    all_support_markers_observed = lineage_expression$all_markers_observed && state_expression$all_markers_observed,
    review_status = mapping$review_status[[i]],
    stringsAsFactors = FALSE
  )
}))
rownames(annotation_evidence) <- NULL
if (nrow(annotation_evidence) != nrow(mapping) || any(!annotation_evidence$all_support_markers_observed)) {
  stop("Per-cluster annotation evidence export failed its completeness checks.", call. = FALSE)
}
specific_lineage <- !annotation_evidence$approved_immune_lineage %in% c("Ambiguous", "Unresolved")
specific_state <- nzchar(annotation_evidence$approved_immune_state) & annotation_evidence$approved_immune_state != "Unresolved"
if (any(specific_lineage & !nzchar(annotation_evidence$lineage_reference_titles))) {
  stop("At least one specific lineage annotation lacks a literature reference.", call. = FALSE)
}
if (any(specific_state & !nzchar(annotation_evidence$state_reference_titles))) {
  stop("At least one specific state annotation lacks a literature reference.", call. = FALSE)
}

lookup <- function(column) unname(stats::setNames(mapping[[column]], mapping$cluster)[cluster_vector])
obj$mouse_parent_immune_review_candidate <- review_mask
obj$mouse_parent_immune_lineage <- NA_character_
obj$mouse_parent_immune_state <- NA_character_
obj$mouse_parent_immune_label <- NA_character_
obj$mouse_parent_immune_lineage_confidence <- NA_character_
obj$mouse_parent_immune_state_confidence <- NA_character_
obj$mouse_parent_immune_lineage_support_markers <- NA_character_
obj$mouse_parent_immune_state_support_markers <- NA_character_
obj$mouse_parent_immune_annotation_notes <- NA_character_
obj$mouse_parent_immune_review_status <- NA_character_
obj$mouse_parent_immune_annotation_level <- NA_character_

obj$mouse_parent_immune_lineage[review_mask] <- lookup("reviewer_immune_lineage")[review_mask]
obj$mouse_parent_immune_state[review_mask] <- lookup("reviewer_immune_state")[review_mask]
obj$mouse_parent_immune_lineage_confidence[review_mask] <- lookup("reviewer_lineage_confidence")[review_mask]
obj$mouse_parent_immune_state_confidence[review_mask] <- lookup("reviewer_state_confidence")[review_mask]
obj$mouse_parent_immune_lineage_support_markers[review_mask] <- lookup("lineage_support_markers_normalized")[review_mask]
obj$mouse_parent_immune_state_support_markers[review_mask] <- lookup("state_support_markers_normalized")[review_mask]
obj$mouse_parent_immune_annotation_notes[review_mask] <- lookup("reviewer_notes")[review_mask]
obj$mouse_parent_immune_review_status[review_mask] <- lookup("review_status")[review_mask]
obj$mouse_parent_immune_annotation_level[review_mask] <- "existing_05d_cluster"
obj$mouse_parent_immune_label[review_mask] <- ifelse(
  nzchar(obj$mouse_parent_immune_state[review_mask]),
  obj$mouse_parent_immune_state[review_mask],
  obj$mouse_parent_immune_lineage[review_mask]
)

lineage_uniformity <- tapply(obj$mouse_parent_immune_lineage[review_mask], cluster_vector[review_mask], function(x) length(unique(x)))
state_uniformity <- tapply(obj$mouse_parent_immune_state[review_mask], cluster_vector[review_mask], function(x) length(unique(x)))
if (any(lineage_uniformity != 1L) || any(state_uniformity != 1L)) stop("Applied immune annotations are not cluster-uniform.", call. = FALSE)
if (any(is.na(obj$mouse_parent_immune_lineage[review_mask]))) stop("At least one reviewed cell lacks an immune lineage label.", call. = FALSE)

cluster_annotation <- mapping
cluster_annotation$n_cells_verified <- as.integer(cluster_counts)
cluster_annotation$percent_full_mouse_object <- 100 * cluster_annotation$n_cells_verified / ncol(obj)
lineage_counts <- as.data.frame(table(
  immune_lineage = factor(obj$mouse_parent_immune_lineage[review_mask], levels = as.character(cfg$immune_fine_annotation$lineage_order))
), stringsAsFactors = FALSE)
colnames(lineage_counts)[2L] <- "n_cells"
lineage_counts <- lineage_counts[lineage_counts$n_cells > 0L, , drop = FALSE]
lineage_counts$percent_review_cells <- 100 * lineage_counts$n_cells / sum(review_mask)
state_counts <- as.data.frame(table(
  immune_state = obj$mouse_parent_immune_state[review_mask], useNA = "no"
), stringsAsFactors = FALSE)
colnames(state_counts)[2L] <- "n_cells"
state_counts <- state_counts[nzchar(state_counts$immune_state) & state_counts$n_cells > 0L, , drop = FALSE]

utils::write.csv(cluster_annotation, file.path(out_dir, "existing_immune_cluster_annotation.csv"), row.names = FALSE, na = "")
utils::write.csv(annotation_evidence, file.path(out_dir, "mouse_existing_immune_cluster_annotation_evidence.csv"), row.names = FALSE, na = "")
utils::write.csv(lineage_counts, file.path(out_dir, "existing_immune_lineage_counts.csv"), row.names = FALSE, na = "")
utils::write.csv(state_counts, file.path(out_dir, "existing_immune_state_counts.csv"), row.names = FALSE, na = "")
utils::write.csv(mapping, file.path(out_dir, basename(mapping_file)), row.names = FALSE, na = "")

review_obj <- subset(obj, cells = review_cells)
lineage_order <- as.character(unlist(cfg$immune_fine_annotation$lineage_order, use.names = FALSE))
lineage_colors <- unlist(cfg$immune_fine_annotation$lineage_colors, use.names = TRUE)
observed_lineages <- lineage_order[lineage_order %in% unique(as.character(review_obj$mouse_parent_immune_lineage))]
if (anyNA(lineage_colors[observed_lineages])) stop("Fine-immune lineage palette is incomplete.", call. = FALSE)
review_obj$mouse_parent_immune_lineage <- factor(as.character(review_obj$mouse_parent_immune_lineage), levels = observed_lineages)
p_lineage <- Seurat::DimPlot(
  review_obj, reduction = "umap", group.by = "mouse_parent_immune_lineage",
  cols = unname(lineage_colors[observed_lineages]), label = TRUE, repel = TRUE, raster = FALSE
) +
  ggplot2::labs(
    title = "Existing 05d clusters: manually approved fine immune lineages",
    subtitle = "One lineage per existing cluster; no cell-level classifier was used",
    color = "Immune lineage"
  ) +
  ggplot2::theme_bw(base_size = 11)

state_order <- as.character(unlist(cfg$immune_fine_annotation$state_order, use.names = FALSE))
state_colors <- unlist(cfg$immune_fine_annotation$state_colors, use.names = TRUE)
lineage_or_state <- ifelse(
  nzchar(as.character(review_obj$mouse_parent_immune_state)),
  as.character(review_obj$mouse_parent_immune_state),
  as.character(review_obj$mouse_parent_immune_lineage)
)
combined_order <- unique(c(
  state_order[state_order %in% lineage_or_state],
  lineage_order[lineage_order %in% lineage_or_state]
))
combined_colors <- c(state_colors, lineage_colors)
combined_colors <- combined_colors[!duplicated(names(combined_colors))]
if (anyNA(combined_colors[combined_order])) stop("Fine-immune lineage/state palette is incomplete.", call. = FALSE)
review_obj$mouse_parent_immune_lineage_or_state <- factor(lineage_or_state, levels = combined_order)
p_state <- Seurat::DimPlot(
  review_obj, reduction = "umap", group.by = "mouse_parent_immune_lineage_or_state",
  cols = unname(combined_colors[combined_order]), label = TRUE, repel = TRUE, raster = FALSE
) +
  ggplot2::labs(
    title = "Existing 05d clusters: manually approved immune lineage and state",
    subtitle = "State labels are displayed only after parent-lineage validation",
    color = "Lineage or state"
  ) +
  ggplot2::theme_bw(base_size = 11)
write_plot_pair(p_lineage, "umap_existing_immune_lineage_approved", out_dir, width = 11, height = 8)
write_plot_pair(p_state, "umap_existing_immune_lineage_state_approved", out_dir, width = 13, height = 9)

obj@misc$mouse05_existing_immune_annotation_contract <- list(
  annotation_level = "existing_05d_cluster",
  annotation_method = "manually approved mapping after 05g literature-marker and existing-cluster DEG review",
  review_candidate_cells = length(review_cells),
  review_clusters = review_clusters,
  state_requires_validated_parent_lineage = TRUE,
  nsg_lymphoid_qc_used_as_expected_lineage = FALSE,
  cell_level_classifier_used = FALSE,
  human_features_retained_in_RNA = human_feature_n,
  mapping_source = normalizePath(mapping_file, mustWork = TRUE),
  review_completion = normalizePath(review_completion, mustWork = TRUE)
)

output_rds <- file.path(out_dir, "mouse_existing_immune_clusters_annotated.rds")
saveRDS(obj, output_rds, compress = FALSE)

input_files <- c(
  input_rds, review_completion, review_manifest, review_expression_file,
  normalizePath(config_path, mustWork = TRUE), script_path, mapping_file, resource_file
)
input_manifest <- data.frame(
  role = c(
    "05f_manually_approved_broad_object", "05g_review_completion", "05g_output_manifest",
    "05g_cluster_marker_expression", "analysis_config", "analysis_script",
    "approved_existing_cluster_mapping", "immune_marker_resource"
  ),
  path = normalizePath(input_files, mustWork = TRUE),
  size_bytes = file.info(input_files)$size,
  sha256 = vapply(input_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(input_manifest, file.path(out_dir, "input_manifest_sha256.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05h_existing_immune_annotation.txt"))
writeLines(
  c(
    "05h manually approved existing-cluster immune annotation completed.",
    paste0("Full mouse cells: ", ncol(obj)),
    paste0("Reviewed immune-candidate cells: ", length(review_cells)),
    paste0("Reviewed existing clusters: ", length(review_clusters)),
    paste0("Human features in RNA assay: ", human_feature_n),
    "Annotation level: existing 05d cluster only.",
    "Exactly one immune lineage label was applied per reviewed existing cluster.",
    "Lineage-specific states were accepted only after parent-lineage validation.",
    "Per-cluster DEGs, marker-expression evidence, annotation reasons, PMID/DOI references, and URLs were exported to mouse_existing_immune_cluster_annotation_evidence.csv.",
    "No cell-level classifier was used.",
    "Labels stored in cell metadata are direct copies of the approved cluster mapping."
  ),
  file.path(out_dir, "completion.txt")
)

output_files <- setdiff(list.files(out_dir, full.names = TRUE), file.path(out_dir, "output_manifest_sha256.csv"))
output_manifest <- data.frame(
  file = basename(output_files),
  size_bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(output_manifest, file.path(out_dir, "output_manifest_sha256.csv"), row.names = FALSE)

message("[05h] Manually approved existing immune-cluster annotation completed. Output: ", out_dir)
