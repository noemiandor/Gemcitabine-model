#!/usr/bin/env Rscript

# Raw table builder for Supplementary Figures 4-7.
#
# Descriptive UMAP/composition views plus standalone Hallmark cluster analysis.
# The deposited Seurat RDS supplies the complete cell universe, reviewed
# metadata, UMAP, and expression. all_ploidy.tsv supplies endpoint tumor
# ploidy. scVelo metrics are optional audit-only input and never supply values
# used by a Supplementary Figure.
#
# SI7 has an explicit inferential boundary: it extracts the RNA counts, retains
# only exact GRCh38-prefixed features, constructs a new SI7-only Seurat object,
# and normalizes those human counts from scratch before differential expression
# and enrichment. The generated cache remains noncanonical until a genuine raw
# rerun has been reviewed and must never overwrite the frozen Data cache.

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
      stop("Unexpected positional argument: ", token, call. = FALSE)
    }
  }
  out
}

arg_value <- function(args, name, default = NULL) {
  value <- args[[name]]
  if (!is.null(value) && length(value) == 1L && nzchar(value)) value else default
}

arg_flag <- function(args, name, default = FALSE) {
  value <- arg_value(args, name, NULL)
  if (is.null(value)) return(default)
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

resolve_work_dependency_fingerprint <- function(
  supplied,
  authoritative
) {
  authoritative <- tolower(as.character(authoritative))
  if (length(authoritative) != 1L ||
      !grepl("^[0-9a-f]{64}$", authoritative)) {
    stop("Authoritative work dependency fingerprint is invalid", call. = FALSE)
  }
  if (is.null(supplied) || !nzchar(supplied)) return(authoritative)
  supplied <- tolower(as.character(supplied))
  if (length(supplied) != 1L ||
      !grepl("^[0-9a-f]{64}$", supplied)) {
    stop("--work-dependency-fingerprint must be one SHA-256 value", call. = FALSE)
  }
  if (!identical(supplied, authoritative)) {
    stop(
      "Supplied --work-dependency-fingerprint does not match the ",
      "authoritative fingerprint derived from actual inputs, the scoped ",
      "SI config, and the exact environment lock",
      call. = FALSE
    )
  }
  authoritative
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript Code/in-vivo/SI_figures/build_raw_supplementary_tables.R \\",
      "    --all-ploidy Data/in-vivo/all_ploidy.tsv \\",
      "    --seurat-rds /absolute/path/to/integrated_sct_cca_seurat_final_reclustered.rds \\",
      "    [--scvelo-metrics /path/to/scvelo_cell_metrics.csv] \\",
      "    --config Code/in-vivo/figure7/figure7_config.yaml \\",
      "    --output-dir Results/in-vivo/SI_figures/raw_table_build",
      "",
      "Options:",
      "  --audit-only      Validate optional scVelo metrics against the RDS and exit.",
      "  --overwrite       Replace a previous raw-table build at output-dir.",
      "  --help             Show this help message.",
      sep = "\n"
    ),
    "\n"
  )
}

resolve_script_path <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command, value = TRUE)
  if (!length(file_arg)) return(NA_character_)
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
}

resolve_builder_source_path <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      if (is.null(frame$ofile)) NA_character_ else as.character(frame$ofile)
    },
    character(1L)
  )
  command_file <- sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)
  )
  candidates <- c(
    rev(frame_files[!is.na(frame_files) & nzchar(frame_files)]),
    command_file,
    file.path(
      getwd(),
      "Code", "in-vivo", "SI_figures",
      "build_raw_supplementary_tables.R"
    )
  )
  candidates <- candidates[
    basename(candidates) == "build_raw_supplementary_tables.R"
  ]
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) {
    stop("Cannot resolve SI raw-table builder path", call. = FALSE)
  }
  normalizePath(candidates[[1L]], mustWork = TRUE)
}

si7_builder_source_path <- resolve_builder_source_path()
si7_feature_species_helper_path <- file.path(
  dirname(dirname(si7_builder_source_path)),
  "figure7", "src", "feature_species_policy.R"
)
if (!file.exists(si7_feature_species_helper_path)) {
  stop(
    "Missing shared feature-species policy helper: ",
    si7_feature_species_helper_path,
    call. = FALSE
  )
}
sys.source(si7_feature_species_helper_path, envir = environment())

is_absolute_path <- function(path) grepl("^/", path)

resolve_path <- function(path, repo_root, must_work = FALSE) {
  resolved <- if (is_absolute_path(path)) path else file.path(repo_root, path)
  normalizePath(resolved, mustWork = must_work)
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Missing required R package: ", package, call. = FALSE)
  }
}

read_csv_checked <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
  if (file.info(path)$size <= 0) stop(label, " is empty: ", path, call. = FALSE)
  out <- utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN")
  )
  if (!nrow(out)) stop(label, " has no data rows: ", path, call. = FALSE)
  out
}

read_tsv_checked <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
  if (file.info(path)$size <= 0) stop(label, " is empty: ", path, call. = FALSE)
  out <- utils::read.delim(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    quote = "", comment.char = "", na.strings = c("", "NA", "NaN")
  )
  if (!nrow(out)) stop(label, " has no data rows: ", path, call. = FALSE)
  out
}

assert_single_column <- function(df, column, label) {
  n <- sum(names(df) == column)
  if (n != 1L) {
    stop(label, " must contain exactly one `", column, "` column; found ", n, call. = FALSE)
  }
  invisible(column)
}

resolve_column <- function(df, candidates, label, required = TRUE) {
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  exact <- candidates[candidates %in% names(df)]
  if (length(exact)) {
    assert_single_column(df, exact[[1L]], label)
    return(exact[[1L]])
  }
  lower_names <- tolower(names(df))
  for (candidate in candidates) {
    index <- which(lower_names == tolower(candidate))
    if (length(index) == 1L) return(names(df)[index])
    if (length(index) > 1L) {
      stop(label, " has ambiguous case-insensitive matches for `", candidate, "`", call. = FALSE)
    }
  }
  if (required) {
    stop(label, " is missing a required column; tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  NA_character_
}

clean_character <- function(values) {
  out <- trimws(enc2utf8(as.character(values)))
  out[is.na(values) | out == "" | out %in% c("NA", "NaN", "None")] <- NA_character_
  out
}

numeric_strict <- function(values, label) {
  raw <- clean_character(values)
  out <- suppressWarnings(as.numeric(raw))
  bad <- !is.na(raw) & !is.finite(out)
  if (any(bad)) {
    stop(label, " contains nonnumeric/nonfinite values, including: ",
         paste(head(unique(raw[bad]), 10L), collapse = ", "), call. = FALSE)
  }
  out
}

standardize_dose <- function(values, label = "Dose") {
  raw <- clean_character(values)
  normalized <- tolower(gsub("[[:space:]]+", "", raw))
  out <- rep(NA_character_, length(raw))
  out[normalized %in% c("0", "0mg", "0mg/kg")] <- "0mg/kg"
  out[normalized %in% c("30", "30mg", "30mg/kg")] <- "30mg/kg"
  out[normalized %in% c("120", "120mg", "120mg/kg")] <- "120mg/kg"
  unknown <- !is.na(raw) & is.na(out)
  if (any(unknown)) {
    stop(label, " contains unexpected values: ",
         paste(sort(unique(raw[unknown])), collapse = ", "), call. = FALSE)
  }
  out
}

standardize_ploidy <- function(values) {
  raw <- clean_character(values)
  upper <- toupper(raw)
  out <- rep(NA_character_, length(raw))
  out[!is.na(upper) & grepl("2N", upper, fixed = TRUE)] <- "2N"
  out[is.na(out) & !is.na(upper) & grepl("4N", upper, fixed = TRUE)] <- "4N"
  out
}

standardize_context <- function(values) {
  raw <- tolower(clean_character(values))
  out <- rep(NA_character_, length(raw))
  out[!is.na(raw) & grepl("cell[-_ ]?(line|culture)|cellline", raw)] <- "CellLine"
  out[is.na(out) & !is.na(raw) & grepl("tumou?r", raw)] <- "Tumor"
  out
}

sort_cluster_levels <- function(values) {
  values <- sort(unique(clean_character(values)))
  numeric_prefix <- suppressWarnings(as.numeric(sub("^([0-9]+).*$", "\\1", values)))
  has_prefix <- grepl("^[0-9]+", values)
  order_key <- ifelse(has_prefix, numeric_prefix, Inf)
  values[order(order_key, values, na.last = TRUE, method = "radix")]
}

si_assert_reviewed_final_universe <- function(
  data,
  included_in_si_figures = NULL
) {
  required <- c(
    "cell", "mouse", "cluster", "context", "initial_ploidy", "dose"
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "Final-QC SI universe lacks required field(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  cell <- clean_character(data$cell)
  mouse <- clean_character(data$mouse)
  cluster <- clean_character(data$cluster)
  context <- clean_character(data$context)
  initial_ploidy <- clean_character(data$initial_ploidy)
  dose <- clean_character(data$dose)
  if (nrow(data) != 35513L || anyNA(cell) || anyDuplicated(cell)) {
    stop(
      "Raw SI analyses require the exact 35,513-cell reviewed final-QC ",
      "Seurat universe with unique cell IDs",
      call. = FALSE
    )
  }

  expected_cluster_counts <- c(
    "0" = 22670L,
    "2" = 3644L,
    "4c" = 103L,
    "5" = 1823L,
    "6" = 1685L,
    "8" = 1579L,
    "10" = 3095L,
    "13" = 704L,
    "14" = 210L
  )
  discarded_clusters <- c("3", "4", "9", "9c")
  observed_discarded <- intersect(unique(cluster), discarded_clusters)
  if (length(observed_discarded)) {
    stop(
      "Raw SI universe still contains discarded QC cluster(s): ",
      paste(sort_cluster_levels(observed_discarded), collapse = ", "),
      call. = FALSE
    )
  }
  observed_clusters <- sort_cluster_levels(cluster)
  observed_cluster_counts <- table(factor(
    cluster,
    levels = names(expected_cluster_counts)
  ))
  if (anyNA(cluster) ||
      !identical(observed_clusters, names(expected_cluster_counts)) ||
      !identical(
        as.integer(observed_cluster_counts),
        as.integer(expected_cluster_counts)
      )) {
    stop(
      "Raw SI universe differs from the reviewed final cluster counts ",
      "after excluding clusters 3, 4, 9, and 9c",
      call. = FALSE
    )
  }

  context_counts <- table(factor(
    context,
    levels = c("CellLine", "Tumor")
  ))
  if (anyNA(context) ||
      !setequal(unique(context), c("CellLine", "Tumor")) ||
      !identical(as.integer(context_counts), c(25681L, 9832L))) {
    stop(
      "Raw SI universe must have the reviewed context counts: ",
      "25,681 CellLine and 9,832 Tumor cells",
      call. = FALSE
    )
  }

  tumor <- context == "Tumor"
  expected_sample_counts <- c(
    "2N-A1-0" = 413L,
    "2N-A1-LR" = 369L,
    "2N-A1-R" = 196L,
    "2N-A1-RR" = 1495L,
    "2N-A2-0" = 317L,
    "2N-A2-L" = 505L,
    "2N-A4-R" = 305L,
    "2N-A4-RL" = 1280L,
    "4N-A5-0" = 888L,
    "4N-A5-RR" = 393L,
    "A5-4N-L" = 385L,
    "A5-4N-R" = 358L,
    "A6-4N-O" = 189L,
    "A6-4N-RR" = 660L,
    "4N-A8-RL" = 1832L,
    "4N-A8-RR" = 247L
  )
  expected_sample_origins <- stats::setNames(
    c(rep("2N", 8L), rep("4N", 8L)),
    names(expected_sample_counts)
  )
  expected_sample_doses <- stats::setNames(
    c(
      rep("0mg/kg", 4L), rep("30mg/kg", 2L), rep("120mg/kg", 2L),
      rep("0mg/kg", 4L), rep("30mg/kg", 2L), rep("120mg/kg", 2L)
    ),
    names(expected_sample_counts)
  )
  observed_sample_counts <- table(factor(
    mouse[tumor],
    levels = names(expected_sample_counts)
  ))
  mapped_sample_origins <- unname(expected_sample_origins[mouse[tumor]])
  mapped_sample_doses <- unname(expected_sample_doses[mouse[tumor]])
  if (anyNA(mouse[tumor]) ||
      !setequal(unique(mouse[tumor]), names(expected_sample_counts)) ||
      !identical(
        as.integer(observed_sample_counts),
        as.integer(expected_sample_counts)
      ) ||
      anyNA(mapped_sample_origins) || anyNA(mapped_sample_doses) ||
      any(initial_ploidy[tumor] != mapped_sample_origins) ||
      any(dose[tumor] != mapped_sample_doses)) {
    stop(
      "Raw SI universe must retain the exact reviewed 16-sample cell counts ",
      "and one injected-origin/dose assignment per tumor sample",
      call. = FALSE
    )
  }

  allowed_tumor_doses <- c("0mg/kg", "30mg/kg", "120mg/kg")
  treated <- tumor & dose %in% c("30mg/kg", "120mg/kg")
  tumor_ploidy_counts <- table(factor(
    initial_ploidy[tumor],
    levels = c("2N", "4N")
  ))
  treated_ploidy_counts <- table(factor(
    initial_ploidy[treated],
    levels = c("2N", "4N")
  ))
  if (anyNA(dose[tumor]) ||
      any(!dose[tumor] %in% allowed_tumor_doses) ||
      anyNA(initial_ploidy[tumor]) ||
      any(!initial_ploidy[tumor] %in% c("2N", "4N")) ||
      sum(treated) != 5335L ||
      sum(tumor & dose == "0mg/kg") != 4497L ||
      !identical(as.integer(tumor_ploidy_counts), c(4880L, 4952L)) ||
      !identical(as.integer(treated_ploidy_counts), c(2407L, 2928L))) {
    stop(
      "Raw SI universe must retain the exact reviewed treated-tumor ",
      "contract: 5,335 treated cells (2N=2,407; 4N=2,928) and ",
      "4,497 untreated cells",
      call. = FALSE
    )
  }

  if (!is.null(included_in_si_figures)) {
    if (!is.logical(included_in_si_figures) ||
        length(included_in_si_figures) != nrow(data) ||
        anyNA(included_in_si_figures) ||
        !identical(included_in_si_figures, tumor) ||
        sum(included_in_si_figures & treated) != 5335L) {
      stop(
        "`included_in_si_figures` must select exactly all 9,832 final-QC ",
        "Tumor cells, including the same 5,335 treated cells, and no ",
        "CellLine cells",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

mapping_by_group <- function(group, value, value_label) {
  group <- clean_character(group)
  value <- clean_character(value)
  keep <- !is.na(group) & !is.na(value)
  pieces <- split(value[keep], group[keep])
  conflicts <- names(pieces)[vapply(pieces, function(x) length(unique(x)) > 1L, logical(1))]
  if (length(conflicts)) {
    stop(
      "Each mouse must map to one ", value_label, "; conflicts include: ",
      paste(head(conflicts, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  if (!length(pieces)) return(setNames(character(0), character(0)))
  setNames(vapply(pieces, function(x) unique(x)[[1L]], character(1)), names(pieces))
}

assert_consistent <- function(left, right, label) {
  conflict <- !is.na(left) & !is.na(right) & left != right
  if (any(conflict)) {
    stop(label, " disagrees for ", sum(conflict), " matched cells", call. = FALSE)
  }
  invisible(TRUE)
}

infer_sample_ploidy <- function(values) {
  out <- standardize_ploidy(values)
  raw <- clean_character(values)
  out[is.na(out) & !is.na(raw) & grepl("^A5", raw, ignore.case = TRUE)] <- "4N"
  out[is.na(out) & !is.na(raw) & grepl("^A6", raw, ignore.case = TRUE)] <- "4N"
  out
}

infer_id_ploidy_context <- function(
  metadata,
  id_col,
  sample_folder_col = NA_character_,
  sample_col = "mouse"
) {
  for (column in c(id_col, sample_col)) {
    assert_single_column(metadata, column, "Seurat metadata")
  }
  if (!is.na(sample_folder_col)) {
    assert_single_column(metadata, sample_folder_col, "Seurat metadata")
  }

  initial_ploidy <- standardize_ploidy(metadata[[id_col]])
  if (!is.na(sample_folder_col)) {
    fallback <- infer_sample_ploidy(metadata[[sample_folder_col]])
    initial_ploidy[is.na(initial_ploidy)] <- fallback[is.na(initial_ploidy)]
  }
  sample_fallback <- infer_sample_ploidy(metadata[[sample_col]])
  initial_ploidy[is.na(initial_ploidy)] <-
    sample_fallback[is.na(initial_ploidy)]

  id_text <- clean_character(metadata[[id_col]])
  context <- rep(NA_character_, nrow(metadata))
  known_id <- !is.na(id_text)
  context[known_id] <- ifelse(
    grepl("cell-culture", tolower(id_text[known_id]), fixed = TRUE),
    "CellLine",
    "Tumor"
  )
  data.frame(
    initial_ploidy = initial_ploidy,
    context = context,
    stringsAsFactors = FALSE
  )
}

resolve_reviewed_ploidy_context <- function(
  metadata,
  ploidy_col,
  context_col,
  sample_col,
  id_col,
  sample_folder_col = NA_character_
) {
  for (column in c(ploidy_col, context_col, sample_col, id_col)) {
    assert_single_column(metadata, column, "Seurat metadata")
  }
  sample <- clean_character(metadata[[sample_col]])
  initial_ploidy <- standardize_ploidy(metadata[[ploidy_col]])
  context <- standardize_context(metadata[[context_col]])
  if (anyNA(sample)) {
    stop("Seurat sample identities cannot be missing", call. = FALSE)
  }
  if (anyNA(initial_ploidy) || anyNA(context)) {
    stop(
      "Configured RDS Ploidy/TN fields contain unresolved values",
      call. = FALSE
    )
  }

  sample_ploidy <- mapping_by_group(
    sample,
    initial_ploidy,
    "reviewed initial ploidy"
  )
  sample_context <- mapping_by_group(
    sample,
    context,
    "reviewed context"
  )
  initial_ploidy <- unname(sample_ploidy[sample])
  context <- unname(sample_context[sample])

  inferred <- infer_id_ploidy_context(
    metadata,
    id_col = id_col,
    sample_folder_col = sample_folder_col,
    sample_col = sample_col
  )
  assert_consistent(
    initial_ploidy,
    inferred$initial_ploidy,
    "Configured RDS ploidy and ID-derived ploidy"
  )
  assert_consistent(
    context,
    inferred$context,
    "Configured RDS context and ID-derived context"
  )
  data.frame(
    initial_ploidy = initial_ploidy,
    context = context,
    stringsAsFactors = FALSE
  )
}

audit_optional_scvelo <- function(
  seurat,
  metrics = NULL,
  canonical_cluster_col
) {
  if (is.null(metrics)) return(seurat)
  for (column in c(
    "cell", "initial_ploidy", "context", "dose", "cluster"
  )) {
    assert_single_column(seurat, column, "Canonical Seurat metadata")
  }

  assert_single_column(metrics, "cell", "scVelo cell metrics")
  metrics$cell <- clean_character(metrics$cell)
  if (anyNA(metrics$cell) || anyDuplicated(metrics$cell) ||
      any(!metrics$cell %in% seurat$cell)) {
    stop(
      "Optional scVelo metrics have invalid or foreign cell IDs",
      call. = FALSE
    )
  }
  metric_index <- match(metrics$cell, seurat$cell)
  metric_ploidy_col <- resolve_column(
    metrics, c("Ploidy", "ploidy"), "scVelo ploidy"
  )
  metric_context_col <- resolve_column(
    metrics,
    c("TN", "trajectory_context", "trajectory_group", "sample_type"),
    "scVelo context"
  )
  metric_dose_col <- resolve_column(
    metrics, c("Dose", "Dose_DEG"), "scVelo dose"
  )
  metric_cluster_col <- resolve_column(
    metrics, canonical_cluster_col, "scVelo canonical cluster"
  )
  metric_ploidy <- standardize_ploidy(metrics[[metric_ploidy_col]])
  metric_context <- standardize_context(metrics[[metric_context_col]])
  metric_dose <- standardize_dose(
    metrics[[metric_dose_col]],
    "scVelo Dose"
  )
  metric_cluster <- clean_character(metrics[[metric_cluster_col]])
  if (anyNA(metric_ploidy) || anyNA(metric_context) ||
      anyNA(metric_cluster)) {
    stop(
      "Optional scVelo metrics contain unresolved ploidy, context, or cluster",
      call. = FALSE
    )
  }
  if (any(metric_context == "Tumor" & is.na(metric_dose))) {
    stop(
      "Optional scVelo metrics contain missing tumor dose",
      call. = FALSE
    )
  }
  assert_consistent(
    seurat$initial_ploidy[metric_index],
    metric_ploidy,
    "Seurat/scVelo ploidy audit"
  )
  assert_consistent(
    seurat$context[metric_index],
    metric_context,
    "Seurat/scVelo context audit"
  )
  assert_consistent(
    seurat$dose[metric_index],
    metric_dose,
    "Seurat/scVelo dose audit"
  )
  assert_consistent(
    seurat$cluster[metric_index],
    metric_cluster,
    "Seurat/scVelo cluster audit"
  )
  seurat
}

file_sha256 <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(file = path, algo = "sha256", serialize = FALSE))
  }
  sha256sum <- Sys.which("sha256sum")
  if (nzchar(sha256sum)) {
    result <- system2(sha256sum, shQuote(path), stdout = TRUE, stderr = TRUE)
  } else {
    shasum <- Sys.which("shasum")
    if (!nzchar(shasum)) return(NA_character_)
    result <- system2(shasum, c("-a", "256", shQuote(path)), stdout = TRUE, stderr = TRUE)
  }
  if (!length(result)) return(NA_character_)
  strsplit(trimws(result[[1L]]), "[[:space:]]+")[[1L]][[1L]]
}

write_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, na = "", quote = TRUE)
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Failed to write output table: ", path, call. = FALSE)
  }
  invisible(path)
}

write_tsv <- function(df, path) {
  utils::write.table(
    df,
    path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = ""
  )
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Failed to write output table: ", path, call. = FALSE)
  }
  invisible(path)
}

si7_join_rna_layers <- function(object) {
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    object <- tryCatch(
      SeuratObject::JoinLayers(object, assay = "RNA"),
      error = function(error) {
        stop(
          "SI7 cannot join RNA layers before counts extraction: ",
          conditionMessage(error),
          call. = FALSE
        )
      }
    )
  }
  object
}

si7_rna_counts <- function(object) {
  if (!inherits(object, "Seurat") || !("RNA" %in% names(object@assays))) {
    stop("SI7 requires a Seurat object with an RNA assay", call. = FALSE)
  }
  object <- si7_join_rna_layers(object)
  if ("Layers" %in% getNamespaceExports("SeuratObject")) {
    count_layers <- SeuratObject::Layers(
      object[["RNA"]],
      search = "^counts"
    )
    if (length(count_layers) != 1L) {
      stop(
        "SI7 requires exactly one joined RNA counts layer; observed ",
        length(count_layers),
        call. = FALSE
      )
    }
  }
  counts <- tryCatch(
    SeuratObject::LayerData(object, assay = "RNA", layer = "counts"),
    error = function(error) {
      tryCatch(
        Seurat::GetAssayData(object, assay = "RNA", slot = "counts"),
        error = function(second_error) NULL
      )
    }
  )
  if (is.null(counts) || !nrow(counts) || !ncol(counts)) {
    stop("SI7 cannot read nonempty RNA counts", call. = FALSE)
  }
  counts
}

si7_build_human_only_object <- function(
  source_object,
  human_prefix = "GRCh38-",
  mouse_prefix = "GRCm39-",
  scale_factor = 10000
) {
  if (!requireNamespace("Seurat", quietly = TRUE) ||
      !requireNamespace("SeuratObject", quietly = TRUE)) {
    stop(
      "SI7 human-only reconstruction requires Seurat/SeuratObject",
      call. = FALSE
    )
  }
  counts <- si7_rna_counts(source_object)
  filtered <- figure7_filter_human_feature_matrix(
    counts,
    analysis = "SI7 RNA counts",
    human_prefix = human_prefix,
    mouse_prefix = mouse_prefix
  )
  cells <- colnames(filtered$matrix)
  metadata <- source_object@meta.data
  if (is.null(cells) || !identical(cells, rownames(metadata))) {
    stop(
      "SI7 RNA counts must contain the exact ordered source-cell inventory",
      call. = FALSE
    )
  }
  metadata <- metadata[cells, , drop = FALSE]
  metadata <- metadata[
    ,
    setdiff(
      colnames(metadata),
      c("nCount_RNA", "nFeature_RNA", "percent.mt")
    ),
    drop = FALSE
  ]
  object <- Seurat::CreateSeuratObject(
    counts = filtered$matrix,
    meta.data = metadata,
    project = "SI7_human_only",
    min.cells = 0,
    min.features = 0
  )
  object <- Seurat::NormalizeData(
    object,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = scale_factor,
    verbose = FALSE
  )
  if (!identical(colnames(object), cells)) {
    stop("SI7 human-only reconstruction changed the cell inventory",
         call. = FALSE)
  }
  figure7_assert_human_feature_names(
    rownames(object[["RNA"]]),
    analysis = "SI7 normalized RNA assay",
    human_prefix = human_prefix,
    mouse_prefix = mouse_prefix
  )
  audit <- filtered$audit
  audit$n_cells <- ncol(object)
  audit$normalization_method <- "LogNormalize"
  audit$normalization_scale_factor <- scale_factor
  audit$source_data_layer_reused <- FALSE
  list(object = object, audit = audit)
}

si7_case_sensitive_intersection <- function(left, right) {
  intersect(as.character(left), as.character(right))
}

si7_deg_cache_paths <- function(work_dir, cluster_id) {
  cluster_stub <- gsub(
    "[^A-Za-z0-9._-]+",
    "_",
    as.character(cluster_id)
  )
  cluster_stub <- sub(
    "^_|_$",
    "",
    gsub("_+", "_", cluster_stub)
  )
  stem <- paste0("si_figure7_cluster_", cluster_stub, "_vs_rest_DEG")
  list(
    table = file.path(work_dir, paste0(stem, ".csv")),
    contract = file.path(work_dir, paste0(stem, ".cache.tsv"))
  )
}

si7_read_reusable_cluster_deg <- function(
  work_dir,
  cluster_id,
  dependency_fingerprint
) {
  paths <- si7_deg_cache_paths(work_dir, cluster_id)
  if (!file.exists(paths$table) || !file.exists(paths$contract)) {
    return(NULL)
  }
  contract <- tryCatch(
    utils::read.delim(
      paths$contract,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(error) NULL
  )
  if (is.null(contract) ||
      !identical(names(contract), c("key", "value")) ||
      anyDuplicated(contract$key)) {
    return(NULL)
  }
  values <- stats::setNames(as.character(contract$value), contract$key)
  if (!identical(values[["schema_version"]], "1") ||
      !identical(values[["cluster"]], as.character(cluster_id)) ||
      !identical(
        values[["dependency_fingerprint"]],
        dependency_fingerprint
      ) ||
      !identical(values[["table_sha256"]], file_sha256(paths$table))) {
    return(NULL)
  }
  de <- tryCatch(
    read_csv_checked(paths$table, paste0("cached DEG cluster ", cluster_id)),
    error = function(error) NULL
  )
  required <- c("gene", "gene_symbol", "cluster", "p_val_adj")
  lfc_columns <- c("avg_log2FC", "avg_logFC", "log2FC", "logFC")
  if (is.null(de) ||
      !all(required %in% names(de)) ||
      !any(lfc_columns %in% names(de)) ||
      !nrow(de) ||
      any(as.character(de$cluster) != as.character(cluster_id)) ||
      anyNA(de$gene) ||
      any(!nzchar(as.character(de$gene))) ||
      anyDuplicated(de$gene)) {
    return(NULL)
  }
  valid_human <- tryCatch(
    {
      figure7_assert_human_feature_names(
        as.character(de$gene),
        analysis = "reusable SI7 DEG cache"
      )
      parsed <- figure7_classify_feature_species(as.character(de$gene))
      identical(as.character(de$gene_symbol), parsed$symbol)
    },
    error = function(error) FALSE
  )
  if (!isTRUE(valid_human)) return(NULL)
  list(cluster = as.character(cluster_id), de = de, reused = TRUE)
}

si7_write_cluster_deg_cache <- function(
  de,
  work_dir,
  cluster_id,
  dependency_fingerprint
) {
  if (!all(c("gene", "gene_symbol") %in% names(de))) {
    stop("SI7 DEG cache requires gene and gene_symbol columns",
         call. = FALSE)
  }
  figure7_assert_human_feature_names(
    as.character(de$gene),
    analysis = "SI7 DEG cache write"
  )
  parsed <- figure7_classify_feature_species(as.character(de$gene))
  if (!identical(as.character(de$gene_symbol), parsed$symbol)) {
    stop("SI7 DEG cache symbols disagree with exact GRCh feature IDs",
         call. = FALSE)
  }
  paths <- si7_deg_cache_paths(work_dir, cluster_id)
  table_tmp <- paste0(paths$table, ".tmp.", Sys.getpid())
  contract_tmp <- paste0(paths$contract, ".tmp.", Sys.getpid())
  on.exit(unlink(c(table_tmp, contract_tmp), force = TRUE), add = TRUE)
  write_csv(de, table_tmp)
  contract <- data.frame(
    key = c(
      "schema_version",
      "cluster",
      "dependency_fingerprint",
      "table_sha256"
    ),
    value = c(
      "1",
      as.character(cluster_id),
      dependency_fingerprint,
      file_sha256(table_tmp)
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(contract, contract_tmp)
  if (file.exists(paths$table) && unlink(paths$table, force = TRUE) != 0L) {
    stop("Cannot replace stale cluster DEG cache: ", paths$table, call. = FALSE)
  }
  if (!file.rename(table_tmp, paths$table)) {
    stop("Cannot atomically promote cluster DEG cache: ", paths$table, call. = FALSE)
  }
  if (file.exists(paths$contract) &&
      unlink(paths$contract, force = TRUE) != 0L) {
    stop(
      "Cannot replace stale cluster DEG contract: ",
      paths$contract,
      call. = FALSE
    )
  }
  if (!file.rename(contract_tmp, paths$contract)) {
    stop(
      "Cannot atomically promote cluster DEG contract: ",
      paths$contract,
      call. = FALSE
    )
  }
  invisible(paths)
}

resolve_si_seurat_source_lineage <- function(
  seurat_rds_path,
  deposited_sha256,
  upstream_dir = "",
  module_dir,
  environment_lock,
  config,
  all_ploidy,
  sample_info,
  cellranger_root = ""
) {
  seurat_rds_sha256 <- file_sha256(seurat_rds_path)
  if (!nzchar(upstream_dir) &&
      identical(seurat_rds_sha256, as.character(deposited_sha256))) {
    return(list(
      kind = "deposited",
      dependencies = c(source_seurat_rds = seurat_rds_sha256),
      validation = NULL
    ))
  }
  if (!nzchar(upstream_dir)) {
    stop(
      "A non-deposited Seurat RDS requires --seurat-upstream-dir and ",
      "a strict final/reconstruction manifest chain",
      call. = FALSE
    )
  }
  validation <- figure7_validate_any_seurat_upstream_artifact(
    output_root = upstream_dir,
    module_dir = module_dir,
    environment_lock = environment_lock,
    config = config,
    all_ploidy = all_ploidy,
    sample_info = sample_info,
    expected_rds = seurat_rds_path,
    cellranger_root = if (
      nzchar(cellranger_root) && dir.exists(cellranger_root)
    ) {
      cellranger_root
    } else {
      ""
    }
  )
  list(
    kind = "manifested_reconstruction",
    dependencies = figure7_seurat_source_dependency_values(
      validation$rds_sha256,
      validation$dependencies,
      validation$scientific_code_contracts,
      validation$upstream_common_contract
    ),
    validation = validation
  )
}

load_si_seurat_validation_runtime <- function(
  environment_validator_path,
  seurat_selection_path
) {
  runtime_environment <- environment(resolve_si_seurat_source_lineage)
  sys.source(environment_validator_path, envir = runtime_environment)
  sys.source(seurat_selection_path, envir = runtime_environment)
  required_functions <- c(
    "figure7_read_config",
    "figure7_validate_any_seurat_upstream_artifact"
  )
  loaded <- vapply(
    required_functions,
    exists,
    logical(1L),
    envir = runtime_environment,
    mode = "function",
    inherits = FALSE
  )
  if (any(!loaded)) {
    stop(
      "Failed to load SI Seurat validation helpers into the builder runtime",
      call. = FALSE
    )
  }
  invisible(runtime_environment)
}

si_science_marker <- function(name) {
  invisible(name)
}

si_raw_scientific_code_contract <- function(path = NULL) {
  if (is.null(path)) path <- resolve_script_path()
  if (is.na(path) || !file.exists(path)) {
    stop("Cannot resolve SI builder for its scientific contract", call. = FALSE)
  }
  helper_names <- c(
    "file_sha256",
    "read_csv_checked",
    "read_tsv_checked",
    "write_csv",
    "write_tsv",
    "assert_single_column",
    "resolve_column",
    "clean_character",
    "numeric_strict",
    "standardize_dose",
    "standardize_ploidy",
    "standardize_context",
    "sort_cluster_levels",
    "si_assert_reviewed_final_universe",
    "mapping_by_group",
    "assert_consistent",
    "infer_sample_ploidy",
    "infer_id_ploidy_context",
    "resolve_reviewed_ploidy_context",
    "figure7_human_feature_policy_id",
    "figure7_human_feature_policy_description",
    "figure7_classify_feature_species",
    "figure7_select_human_features",
    "figure7_filter_human_feature_matrix",
    "figure7_assert_human_feature_names",
    "si7_join_rna_layers",
    "si7_rna_counts",
    "si7_build_human_only_object",
    "si7_case_sensitive_intersection",
    "si7_deg_cache_paths",
    "si7_read_reusable_cluster_deg",
    "si7_write_cluster_deg_cache"
  )
  helper_contracts <- vapply(
    helper_names,
    function(name) {
      fn <- get(name, mode = "function", inherits = TRUE)
      paste(deparse(fn, width.cutoff = 500L), collapse = "\n")
    },
    character(1L)
  )
  names(helper_contracts) <- paste0("function:", helper_names)

  expressions <- parse(path, keep.source = FALSE)
  main_assignment <- Filter(
    function(expression) {
      is.call(expression) &&
        identical(as.character(expression[[1L]]), "<-") &&
        identical(as.character(expression[[2L]]), "main")
    },
    as.list(expressions)
  )
  if (length(main_assignment) != 1L) {
    stop("SI scientific contract cannot locate exactly one main()", call. = FALSE)
  }
  main_function <- main_assignment[[1L]][[3L]]
  main_body <- as.list(main_function[[3L]])[-1L]
  marker_name <- vapply(
    main_body,
    function(expression) {
      if (!is.call(expression) ||
          !identical(
            as.character(expression[[1L]]),
            "si_science_marker"
          ) ||
          length(expression) != 2L) {
        return("")
      }
      as.character(expression[[2L]])
    },
    character(1L)
  )
  regions <- c(
    "metadata",
    "endpoint_composition",
    "si7",
    "final_tables",
    "plot_table_outputs"
  )
  region_contracts <- vapply(
    regions,
    function(region) {
      begin <- which(marker_name == paste0("BEGIN:", region))
      end <- which(marker_name == paste0("END:", region))
      if (length(begin) != 1L || length(end) != 1L ||
          begin >= end) {
        stop(
          "SI scientific contract has invalid markers for ",
          region,
          call. = FALSE
        )
      }
      paste(
        vapply(
          main_body[seq.int(begin + 1L, end - 1L)],
          function(expression) {
            paste(
              deparse(expression, width.cutoff = 500L),
              collapse = "\n"
            )
          },
          character(1L)
        ),
        collapse = "\n"
      )
    },
    character(1L)
  )
  names(region_contracts) <- paste0("region:", regions)
  values <- c(helper_contracts, region_contracts)
  values <- values[order(names(values))]
  digest::digest(
    paste(names(values), values, sep = "=", collapse = "\n"),
    algo = "sha256",
    serialize = FALSE
  )
}

main <- function(args = parse_cli_args(commandArgs(trailingOnly = TRUE))) {
if (arg_flag(args, "help", FALSE)) {
  usage()
  quit(save = "no", status = 0L)
}

require_package("yaml")
require_package("digest")

script_path <- resolve_script_path()
if (is.na(script_path)) stop("Cannot resolve the script path", call. = FALSE)
script_dir <- dirname(script_path)
repo_root <- normalizePath(dirname(dirname(dirname(script_dir))), mustWork = TRUE)
environment_validator_path <- file.path(
  repo_root,
  "Code", "in-vivo", "figure7", "src", "common_io.R"
)
seurat_selection_path <- file.path(
  repo_root,
  "Code", "in-vivo", "figure7", "src",
  "seurat_upstream_selection.R"
)
load_si_seurat_validation_runtime(
  environment_validator_path,
  seurat_selection_path
)
seurat_generator_path <- file.path(
  repo_root,
  "Code", "in-vivo", "figure7",
  "generate_final_seurat_from_cellranger.R"
)
seurat_upstream_helper_path <- file.path(
  repo_root,
  "Code", "in-vivo", "figure7", "src",
  "seurat_upstream.R"
)
config_path <- resolve_path(
  arg_value(args, "config", file.path("Code", "in-vivo", "figure7", "figure7_config.yaml")),
  repo_root,
  must_work = TRUE
)
config <- figure7_read_config(config_path)
si_config <- config$si_figures
if (is.null(si_config)) stop("Figure 7 config is missing si_figures", call. = FALSE)
feature_species_contract <- figure7_feature_species_contract_values(config)
si7_human_prefix <- unname(feature_species_contract[["human_prefix"]])
si7_mouse_prefix <- unname(feature_species_contract[["mouse_prefix"]])
analysis_seed <- as.integer(si_config$raw_rebuild_analysis_seed)
if (length(analysis_seed) != 1L || !is.finite(analysis_seed)) {
  stop("si_figures.raw_rebuild_analysis_seed must be one finite integer", call. = FALSE)
}
fgsea_nperm_simple <- figure7_si_raw_fgsea_nperm_simple(
  config,
  arg_value(args, "fgsea-nperm-simple", NULL)
)

all_ploidy_path <- resolve_path(
  arg_value(
    args,
    "all-ploidy",
    as.character(
      config$versioned_source_artifacts$endpoint_ploidy$default_path
    )
  ),
  repo_root,
  must_work = FALSE
)
all_ploidy_sha256 <- if (file.exists(all_ploidy_path)) {
  file_sha256(all_ploidy_path)
} else {
  NA_character_
}
if (!identical(
  all_ploidy_sha256,
  as.character(config$versioned_source_artifacts$endpoint_ploidy$sha256)
)) {
  stop(
    "Endpoint-ploidy input is missing or differs from the frozen checksum: ",
    all_ploidy_path,
    call. = FALSE
  )
}
seurat_rds_arg <- arg_value(args, "seurat-rds", NULL)
if (is.null(seurat_rds_arg)) {
  stop("SI Figures 4-7 require --seurat-rds", call. = FALSE)
}
seurat_rds_path <- resolve_path(
  seurat_rds_arg,
  repo_root,
  must_work = TRUE
)
environment_lock_path <- resolve_path(
  as.character(config$raw_data$environment_lock),
  repo_root,
  must_work = TRUE
)
sample_info_path <- resolve_path(
  arg_value(
    args,
    "sample-info",
    as.character(
      config$versioned_source_artifacts$sample_info$default_path
    )
  ),
  repo_root,
  must_work = TRUE
)
upstream_arg <- arg_value(args, "seurat-upstream-dir", "")
upstream_dir <- if (nzchar(upstream_arg)) {
  resolve_path(
    upstream_arg,
    repo_root,
    must_work = TRUE
  )
} else {
  ""
}
cellranger_arg <- arg_value(args, "cellranger-root", "")
cellranger_root <- if (nzchar(cellranger_arg)) {
  resolve_path(cellranger_arg, repo_root, must_work = FALSE)
} else {
  ""
}
seurat_lineage <- resolve_si_seurat_source_lineage(
  seurat_rds_path = seurat_rds_path,
  deposited_sha256 = as.character(config$raw_data$seurat_rds_sha256),
  upstream_dir = upstream_dir,
  module_dir = file.path(repo_root, "Code", "in-vivo", "figure7"),
  environment_lock = environment_lock_path,
  config = config,
  all_ploidy = all_ploidy_path,
  sample_info = sample_info_path,
  cellranger_root = cellranger_root
)
seurat_rds_sha256 <- unname(
  seurat_lineage$dependencies[["source_seurat_rds"]]
)
seurat_source_kind <- seurat_lineage$kind
seurat_source_dependencies <- seurat_lineage$dependencies
seurat_validation <- seurat_lineage$validation
builder_sha256 <- file_sha256(script_path)
environment_validator_sha256 <- file_sha256(environment_validator_path)
environment_lock_sha256 <- file_sha256(environment_lock_path)
config_sha256 <- file_sha256(config_path)
si_raw_scientific_code_contract_sha256 <-
  si_raw_scientific_code_contract(script_path)
si_raw_config_contract_sha256 <-
  figure7_si_raw_config_contract_sha256(config)
si_raw_runtime_contract_sha256 <-
  figure7_environment_stage_contract_sha256(
    environment_lock_path,
    "r",
    "si_raw"
  )
worker_budget <- suppressWarnings(as.numeric(arg_value(args, "workers", "1")))
if (length(worker_budget) != 1L || !is.finite(worker_budget) ||
    worker_budget < 1L || worker_budget > 64L ||
    worker_budget != floor(worker_budget)) {
  stop("--workers must be one integer from 1 to 64", call. = FALSE)
}
worker_budget <- as.integer(worker_budget)
work_dependency_values <- figure7_si_raw_dependency_values(
  si_raw_scientific_code_contract_sha256,
  si_raw_config_contract_sha256,
  si_raw_runtime_contract_sha256,
  all_ploidy_sha256,
  seurat_source_dependencies
)
authoritative_work_dependency_fingerprint <-
  figure7_si_raw_work_fingerprint(work_dependency_values)
work_dependency_fingerprint <- resolve_work_dependency_fingerprint(
  arg_value(args, "work-dependency-fingerprint", NULL),
  authoritative_work_dependency_fingerprint
)
scvelo_arg <- arg_value(args, "scvelo-metrics", NULL)
scvelo_metrics_path <- if (is.null(scvelo_arg)) {
  NA_character_
} else {
  resolve_path(scvelo_arg, repo_root, must_work = TRUE)
}
scvelo_metrics_sha256 <- if (is.na(scvelo_metrics_path)) {
  NA_character_
} else {
  file_sha256(scvelo_metrics_path)
}
audit_only <- arg_flag(args, "audit-only", FALSE)
if (audit_only && is.na(scvelo_metrics_path)) {
  stop("--audit-only requires --scvelo-metrics", call. = FALSE)
}
locked_r_versions <- figure7_validate_r_environment(
  environment_lock_path,
  "si_raw",
  if (audit_only) {
    c("R", "digest", "yaml", "Seurat", "SeuratObject")
  } else {
    NULL
  }
)
r_runtime_provenance <- figure7_r_runtime_provenance()
output_dir <- resolve_path(
  arg_value(args, "output-dir", file.path("Results", "in-vivo", "SI_figures")),
  repo_root,
  must_work = FALSE
)
overwrite <- arg_flag(args, "overwrite", FALSE)

require_package("Seurat")
require_package("SeuratObject")
message("Reading selected Seurat RDS: ", seurat_rds_path)
object <- readRDS(seurat_rds_path)
if (!inherits(object, "Seurat")) {
  stop("--seurat-rds is not a Seurat object", call. = FALSE)
}
si_science_marker("BEGIN:metadata")
seurat_metadata_raw <- object@meta.data
seurat_cells <- rownames(seurat_metadata_raw)
if (anyNA(seurat_cells) || any(!nzchar(seurat_cells)) ||
    anyDuplicated(seurat_cells)) {
  stop("Seurat RDS has invalid cell identifiers", call. = FALSE)
}
umap_name <- as.character(si_config$umap_reduction)
if (!(umap_name %in% names(object@reductions))) {
  stop("Seurat RDS is missing UMAP reduction: ", umap_name, call. = FALSE)
}
umap <- as.data.frame(
  Seurat::Embeddings(object, reduction = umap_name),
  check.names = FALSE
)
if (ncol(umap) < 2L || any(!seurat_cells %in% rownames(umap))) {
  stop("Seurat RDS UMAP cannot be aligned to all cells", call. = FALSE)
}
umap <- umap[seurat_cells, seq_len(2L), drop = FALSE]
names(umap) <- c("UMAP_1", "UMAP_2")
seurat <- data.frame(
  cell = seurat_cells,
  seurat_metadata_raw,
  UMAP_1 = as.numeric(umap$UMAP_1),
  UMAP_2 = as.numeric(umap$UMAP_2),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

for (column in c("cell", "UMAP_1", "UMAP_2", "S.Score", "harvest", "barcode_raw")) {
  assert_single_column(seurat, column, "Seurat metadata")
}

canonical_cluster_col <- resolve_column(
  seurat,
  as.character(si_config$cluster_id_field),
  "Seurat canonical cluster"
)
annotation_col <- resolve_column(
  seurat,
  as.character(si_config$cluster_annotation_field),
  "Seurat canonical cluster annotation"
)
reviewed_ploidy_col <- resolve_column(
  seurat,
  as.character(si_config$ploidy_field),
  "Seurat reviewed ploidy"
)
reviewed_context_col <- resolve_column(
  seurat,
  as.character(si_config$context_field),
  "Seurat reviewed context"
)
dose_col <- resolve_column(
  seurat,
  as.character(si_config$dose_field),
  "Seurat reviewed dose"
)

seurat$cell <- clean_character(seurat$cell)
if (anyNA(seurat$cell)) stop("Input cell identifiers cannot be missing", call. = FALSE)
if (anyDuplicated(seurat$cell)) stop("Seurat metadata contains duplicated cell identifiers", call. = FALSE)

sample_col <- resolve_column(
  seurat,
  as.character(si_config$sample_field),
  "Seurat metadata sample identity"
)

seurat$UMAP_1 <- numeric_strict(seurat$UMAP_1, "Seurat UMAP_1")
seurat$UMAP_2 <- numeric_strict(seurat$UMAP_2, "Seurat UMAP_2")
seurat$s_phase_score <- numeric_strict(seurat$S.Score, "Seurat S.Score")
if (any(!is.finite(seurat$UMAP_1)) || any(!is.finite(seurat$UMAP_2))) {
  stop("Seurat metadata contains missing or nonfinite UMAP coordinates", call. = FALSE)
}

seurat$mouse <- clean_character(seurat[[sample_col]])
seurat$cluster <- clean_character(seurat[[canonical_cluster_col]])
seurat$cluster_annotation <- clean_character(seurat[[annotation_col]])
seurat$dose <- standardize_dose(
  seurat[[dose_col]],
  paste0("Seurat ", dose_col)
)
if (anyNA(seurat$mouse)) stop("Seurat metadata contains missing mouse/sample identities", call. = FALSE)
if (anyNA(seurat$cluster)) stop("Seurat metadata contains missing canonical cluster labels", call. = FALSE)
if (anyNA(seurat$cluster_annotation)) {
  stop("Seurat metadata contains missing canonical cluster annotations", call. = FALSE)
}
cellcycle_mapping <- unlist(si_config$cellcycle_mapping, use.names = TRUE)
unknown_annotations <- setdiff(unique(seurat$cluster_annotation), names(cellcycle_mapping))
if (length(unknown_annotations)) {
  stop(
    "Cluster annotations are absent from the reviewed CellCycle mapping: ",
    paste(unknown_annotations, collapse = ", "),
    call. = FALSE
  )
}
seurat$cellcycle_classification <- unname(cellcycle_mapping[seurat$cluster_annotation])

id_col <- resolve_column(
  seurat,
  c("IDs", "ID"),
  "Seurat ID used for ploidy/context inference"
)
sample_folder_col <- resolve_column(
  seurat,
  c("sample_folder", "Sequencing.IDs", "sample", "orig.ident", "IDs"),
  "Seurat sample folder",
  required = FALSE
)
ploidy_context <- resolve_reviewed_ploidy_context(
  seurat,
  ploidy_col = reviewed_ploidy_col,
  context_col = reviewed_context_col,
  sample_col = "mouse",
  id_col = id_col,
  sample_folder_col = sample_folder_col
)
seurat$initial_ploidy <- ploidy_context$initial_ploidy
seurat$context <- ploidy_context$context
si_assert_reviewed_final_universe(seurat)
si_science_marker("END:metadata")

if (!is.na(scvelo_metrics_path)) {
  message("Auditing optional scVelo metrics: ", scvelo_metrics_path)
  metrics <- read_csv_checked(scvelo_metrics_path, "scVelo cell metrics")
  seurat <- audit_optional_scvelo(
    seurat,
    metrics,
    canonical_cluster_col = canonical_cluster_col
  )
}

if (audit_only) {
  message("Optional scVelo metrics are consistent with deposited Seurat metadata")
  return(invisible(TRUE))
}

if (dir.exists(output_dir)) {
  existing <- list.files(output_dir, all.files = TRUE, no.. = TRUE)
  if (length(existing) && !overwrite) {
    stop("Output directory is not empty; pass --overwrite to replace known outputs: ", output_dir, call. = FALSE)
  }
  if (length(existing) && overwrite) unlink(output_dir, recursive = TRUE)
}

table_dir <- file.path(output_dir, "tables")
metadata_dir <- file.path(output_dir, "metadata")
work_cache_arg <- arg_value(args, "work-cache-dir", NULL)
work_dir <- if (is.null(work_cache_arg)) {
  file.path(output_dir, "work")
} else {
  resolve_path(work_cache_arg, repo_root, must_work = FALSE)
}
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

si_science_marker("BEGIN:endpoint_composition")
message("Reading endpoint ploidy: ", all_ploidy_path)
endpoint_ploidy <- read_tsv_checked(all_ploidy_path, "endpoint ploidy")
for (column in c("file", "cell_id", "ploidy")) {
  assert_single_column(endpoint_ploidy, column, "endpoint ploidy")
}

endpoint_ploidy$file <- clean_character(endpoint_ploidy$file)
endpoint_ploidy$cell_id <- clean_character(endpoint_ploidy$cell_id)
endpoint_ploidy$ploidy <- numeric_strict(endpoint_ploidy$ploidy, "endpoint ploidy value")
endpoint_ploidy$join_key <- paste(endpoint_ploidy$file, endpoint_ploidy$cell_id, sep = "\r")
if (anyNA(endpoint_ploidy$file) || anyNA(endpoint_ploidy$cell_id) ||
    any(!is.finite(endpoint_ploidy$ploidy)) || anyDuplicated(endpoint_ploidy$join_key)) {
  stop("Endpoint ploidy requires unique, nonmissing file + cell_id keys and finite ploidy", call. = FALSE)
}
seurat$endpoint_file <- paste0(clean_character(seurat$harvest), ".sps.cbs")
seurat$endpoint_cell_id <- clean_character(seurat$barcode_raw)
seurat$endpoint_join_key <- paste(seurat$endpoint_file, seurat$endpoint_cell_id, sep = "\r")
endpoint_index <- match(seurat$endpoint_join_key, endpoint_ploidy$join_key)
seurat$endpoint_ploidy <- endpoint_ploidy$ploidy[endpoint_index]
tumor_endpoint_missing <- seurat$context == "Tumor" & !is.finite(seurat$endpoint_ploidy)
cellline_endpoint_present <- seurat$context == "CellLine" & is.finite(seurat$endpoint_ploidy)
if (any(tumor_endpoint_missing)) {
  stop(
    "Endpoint ploidy is missing for ", sum(tumor_endpoint_missing),
    " tumor cells; first cells: ",
    paste(head(seurat$cell[tumor_endpoint_missing], 10L), collapse = ", "),
    call. = FALSE
  )
}
if (any(cellline_endpoint_present)) {
  stop("Endpoint ploidy unexpectedly matched CellLine cells", call. = FALSE)
}

ploidy_levels <- as.character(unlist(si_config$inclusion$ploidy_levels))
dose_levels <- as.character(unlist(si_config$inclusion$dose_levels))
cluster_levels <- as.character(unlist(si_config$cluster_order))
observed_cluster_levels <- sort_cluster_levels(seurat$cluster)
if (!length(cluster_levels) || !identical(observed_cluster_levels, cluster_levels)) {
  stop(
    "RDS canonical clusters differ from configured cluster_order: observed=",
    paste(observed_cluster_levels, collapse = ","),
    "; configured=",
    paste(cluster_levels, collapse = ","),
    call. = FALSE
  )
}
cluster_order <- stats::setNames(seq_along(cluster_levels), cluster_levels)
cluster_colors <- setNames(
  grDevices::hcl.colors(max(length(cluster_levels), 3L), palette = "Dark 3")[seq_along(cluster_levels)],
  cluster_levels
)
cluster_annotation_map <- mapping_by_group(
  seurat$cluster, seurat$cluster_annotation, "cluster annotation"
)[cluster_levels]
cluster_cellcycle_map <- mapping_by_group(
  seurat$cluster, seurat$cellcycle_classification, "CellCycle classification"
)[cluster_levels]
if (anyNA(cluster_annotation_map) || anyNA(cluster_cellcycle_map)) {
  stop("Every cluster must have one nonmissing annotation and CellCycle classification", call. = FALSE)
}

seurat$cluster_order <- unname(cluster_order[seurat$cluster])
seurat$cluster_color <- unname(cluster_colors[seurat$cluster])
seurat$dose_mg_per_kg <- as.numeric(sub("mg/kg$", "", seurat$dose))
seurat$treatment <- ifelse(
  is.na(seurat$dose), NA_character_,
  ifelse(seurat$dose == "0mg/kg", "Control", "Gemcitabine")
)
required_nonmissing <- as.character(unlist(si_config$inclusion$required_nonmissing))
canonical_required <- list(
  cell_id = seurat$cell,
  UMAP_1 = seurat$UMAP_1,
  UMAP_2 = seurat$UMAP_2,
  sample_id = seurat$mouse,
  cluster_id = seurat$cluster,
  cluster_annotation = seurat$cluster_annotation,
  initial_ploidy = seurat$initial_ploidy,
  dose = seurat$dose,
  cellcycle_classification = seurat$cellcycle_classification
)
unknown_required <- setdiff(required_nonmissing, names(canonical_required))
if (length(unknown_required)) {
  stop("Unknown configured SI Figure 5 required field(s): ", paste(unknown_required, collapse = ", "), call. = FALSE)
}
nonmissing_required <- rep(TRUE, nrow(seurat))
for (field in required_nonmissing) {
  values <- canonical_required[[field]]
  nonmissing_required <- nonmissing_required & !is.na(values)
  if (is.character(values)) nonmissing_required <- nonmissing_required & nzchar(values)
}
included_context <- as.character(si_config$inclusion$context)
seurat$included_in_si_figures <- (
  seurat$context == included_context &
    seurat$initial_ploidy %in% ploidy_levels &
    seurat$dose %in% dose_levels &
    nonmissing_required
)
seurat$exclusion_reason <- ifelse(
  seurat$included_in_si_figures,
  "",
  ifelse(
    is.na(seurat$context) | seurat$context != included_context,
    "context_not_Tumor",
    ifelse(
      is.na(seurat$initial_ploidy) | !seurat$initial_ploidy %in% ploidy_levels,
      "invalid_initial_ploidy",
      ifelse(
        is.na(seurat$dose) | !seurat$dose %in% dose_levels,
        "invalid_dose",
        "missing_required_metadata"
      )
    )
  )
)
si_assert_reviewed_final_universe(
  seurat,
  included_in_si_figures = seurat$included_in_si_figures
)

canonical_cells <- data.frame(
  cell_id = seurat$cell,
  UMAP_1 = seurat$UMAP_1,
  UMAP_2 = seurat$UMAP_2,
  sample_id = seurat$mouse,
  cluster_id = seurat$cluster,
  cluster_annotation = seurat$cluster_annotation,
  cluster_order = seurat$cluster_order,
  cluster_color = seurat$cluster_color,
  initial_ploidy = seurat$initial_ploidy,
  treatment = seurat$treatment,
  dose = seurat$dose,
  dose_mg_per_kg = seurat$dose_mg_per_kg,
  s_phase_score = seurat$s_phase_score,
  endpoint_ploidy = seurat$endpoint_ploidy,
  endpoint_file = seurat$endpoint_file,
  endpoint_cell_id = seurat$endpoint_cell_id,
  cellcycle_classification = seurat$cellcycle_classification,
  context = seurat$context,
  included_in_si_figures = seurat$included_in_si_figures,
  exclusion_reason = seurat$exclusion_reason,
  stringsAsFactors = FALSE
)
canonical_cells <- canonical_cells[
  order(canonical_cells$cell_id, method = "radix"),
  ,
  drop = FALSE
]
if (anyDuplicated(canonical_cells$cell_id)) stop("Canonical cell table contains duplicated cell IDs", call. = FALSE)

tumor <- seurat[seurat$included_in_si_figures, , drop = FALSE]
if (!nrow(tumor)) stop("No tumor cells remain after context filtering", call. = FALSE)
if (anyNA(tumor$initial_ploidy)) stop("Tumor cells contain unresolved initial ploidy", call. = FALSE)
if (any(!tumor$initial_ploidy %in% c("2N", "4N"))) stop("Unexpected tumor ploidy value", call. = FALSE)
if (anyNA(tumor$dose)) stop("Tumor cells contain missing Dose values", call. = FALSE)

mouse_dose <- mapping_by_group(tumor$mouse, tumor$dose, "dose")
mouse_ploidy <- mapping_by_group(tumor$mouse, tumor$initial_ploidy, "initial ploidy")
mouse_context <- mapping_by_group(tumor$mouse, tumor$context, "context")
if (any(mouse_context != "Tumor")) stop("Tumor subset contains a non-tumor mouse", call. = FALSE)

sample_metadata <- unique(tumor[, c("mouse", "initial_ploidy", "dose"), drop = FALSE])
if (nrow(sample_metadata) != length(unique(tumor$mouse))) {
  stop("Mouse metadata are not unique after ploidy/dose resolution", call. = FALSE)
}
sample_metadata$initial_ploidy <- factor(sample_metadata$initial_ploidy, levels = ploidy_levels)
sample_metadata$dose <- factor(sample_metadata$dose, levels = dose_levels)
sample_metadata <- sample_metadata[
  order(sample_metadata$initial_ploidy, sample_metadata$dose, sample_metadata$mouse, method = "radix"),
  ,
  drop = FALSE
]
mouse_levels <- sample_metadata$mouse

tumor$mouse <- factor(tumor$mouse, levels = mouse_levels)
tumor$initial_ploidy <- factor(tumor$initial_ploidy, levels = ploidy_levels)
tumor$dose <- factor(tumor$dose, levels = dose_levels)
tumor$cluster <- factor(tumor$cluster, levels = cluster_levels)

ploidy_colors <- c("2N" = "#4C78A8", "4N" = "#E45756")
dose_colors <- c("0mg/kg" = "#666666", "30mg/kg" = "#d95f02", "120mg/kg" = "#1b9e77")

cluster_key <- data.frame(
  cluster_id = cluster_levels,
  cluster_annotation = unname(cluster_annotation_map[cluster_levels]),
  cellcycle_classification = unname(cluster_cellcycle_map[cluster_levels]),
  cluster_order = seq_along(cluster_levels),
  color = unname(cluster_colors[cluster_levels]),
  n_all_cells = as.integer(table(factor(seurat$cluster, levels = cluster_levels))),
  n_included_tumor_cells = as.integer(table(factor(tumor$cluster, levels = cluster_levels))),
  stringsAsFactors = FALSE
)
if (anyNA(cluster_key) || any(!nzchar(cluster_key$cluster_annotation)) ||
    anyDuplicated(cluster_key$cluster_order) || anyDuplicated(cluster_key$color) ||
    any(!grepl("^#[0-9A-Fa-f]{6}$", cluster_key$color))) {
  stop("Cluster key violates annotation, order, or color requirements", call. = FALSE)
}

composition_table <- as.data.frame(
  table(
    mouse = factor(tumor$mouse, levels = mouse_levels),
    cluster_final = factor(tumor$cluster, levels = cluster_levels)
  ),
  stringsAsFactors = FALSE
)
names(composition_table)[names(composition_table) == "Freq"] <- "n_cells"
composition_table$mouse <- as.character(composition_table$mouse)
composition_table$cluster_final <- as.character(composition_table$cluster_final)
composition_table$total_cells <- ave(composition_table$n_cells, composition_table$mouse, FUN = sum)
composition_table$proportion <- composition_table$n_cells / composition_table$total_cells
composition_table$initial_ploidy <- as.character(mouse_ploidy[composition_table$mouse])
composition_table$dose <- as.character(mouse_dose[composition_table$mouse])
composition_table$cluster_annotation <- unname(cluster_annotation_map[composition_table$cluster_final])
composition_table$cluster_order <- unname(cluster_order[composition_table$cluster_final])
composition_table$denominator_definition <- "all canonical cells with included_in_si_figures=true within mouse"
composition_table <- composition_table[, c(
  "mouse", "initial_ploidy", "dose", "cluster_final", "cluster_annotation",
  "cluster_order", "n_cells", "total_cells", "proportion", "denominator_definition"
)]

mouse_sums <- tapply(composition_table$proportion, composition_table$mouse, sum)
if (any(abs(mouse_sums - 1) > 1e-12)) {
  stop("Cluster proportions do not sum to one within every mouse", call. = FALSE)
}
canonical_included_counts <- table(canonical_cells$sample_id[canonical_cells$included_in_si_figures])
composition_denominators <- stats::setNames(
  composition_table$total_cells[!duplicated(composition_table$mouse)],
  composition_table$mouse[!duplicated(composition_table$mouse)]
)
if (!identical(
  as.integer(canonical_included_counts[names(composition_denominators)]),
  as.integer(composition_denominators)
)) {
  stop("Composition denominators do not reconcile to the canonical included-cell table", call. = FALSE)
}

group_key <- interaction(
  composition_table$initial_ploidy,
  composition_table$dose,
  composition_table$cluster_final,
  drop = TRUE,
  lex.order = TRUE
)
group_pieces <- split(composition_table, group_key)
group_summary <- do.call(rbind, lapply(group_pieces, function(piece) {
  data.frame(
    initial_ploidy = piece$initial_ploidy[[1L]],
    dose = piece$dose[[1L]],
    cluster_final = piece$cluster_final[[1L]],
    cluster_annotation = piece$cluster_annotation[[1L]],
    cluster_order = piece$cluster_order[[1L]],
    n_mice = length(unique(piece$mouse)),
    sum_n_cells = sum(piece$n_cells),
    sum_total_cells = sum(piece$total_cells),
    mean_proportion = mean(piece$proportion),
    sd_proportion = if (nrow(piece) > 1L) stats::sd(piece$proportion) else NA_real_,
    min_proportion = min(piece$proportion),
    max_proportion = max(piece$proportion),
    stringsAsFactors = FALSE
  )
}))
rownames(group_summary) <- NULL
group_summary$initial_ploidy <- factor(group_summary$initial_ploidy, levels = ploidy_levels)
group_summary$dose <- factor(group_summary$dose, levels = dose_levels)
group_summary$cluster_final <- factor(group_summary$cluster_final, levels = cluster_levels)
group_summary <- group_summary[order(
  group_summary$initial_ploidy,
  group_summary$dose,
  group_summary$cluster_order
), , drop = FALSE]

group_sums <- aggregate(mean_proportion ~ initial_ploidy + dose, data = group_summary, sum)
if (any(abs(group_sums$mean_proportion - 1) > 1e-12)) {
  stop("Equal-mouse group mean proportions do not sum to one", call. = FALSE)
}
si_science_marker("END:endpoint_composition")

si_science_marker("BEGIN:si7")
composition_by <- function(data, group_field, fill_field, group_levels, fill_levels) {
  output <- as.data.frame(
    table(
      group_value = factor(data[[group_field]], levels = group_levels),
      fill_value = factor(data[[fill_field]], levels = fill_levels)
    ),
    stringsAsFactors = FALSE
  )
  names(output)[3L] <- "n_cells"
  output$group_value <- as.character(output$group_value)
  output$fill_value <- as.character(output$fill_value)
  output$total_cells <- ave(output$n_cells, output$group_value, FUN = sum)
  output$proportion <- ifelse(
    output$total_cells > 0,
    output$n_cells / output$total_cells,
    0
  )
  output$group_field <- group_field
  output$fill_field <- fill_field
  output$denominator_definition <- paste0(
    "all plotted cells within ",
    group_field
  )
  output
}

s4_context <- composition_by(
  seurat, "cluster", "context", cluster_levels, c("Tumor", "CellLine")
)
s4_ploidy <- composition_by(
  seurat, "cluster", "initial_ploidy", cluster_levels, ploidy_levels
)
s5_dose_composition <- composition_by(
  tumor, "cluster", "dose", cluster_levels, dose_levels
)
s5_ploidy_composition <- composition_by(
  tumor, "cluster", "initial_ploidy", cluster_levels, ploidy_levels
)

safe_cluster_stub <- function(x) {
  out <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  out <- gsub("_+", "_", out)
  sub("^_|_$", "", out)
}

clean_gene_symbol <- function(x) {
  feature <- as.character(x)
  figure7_assert_human_feature_names(
    feature,
    analysis = "SI7 differential-expression genes",
    human_prefix = si7_human_prefix,
    mouse_prefix = si7_mouse_prefix
  )
  figure7_classify_feature_species(
    feature,
    human_prefix = si7_human_prefix,
    mouse_prefix = si7_mouse_prefix
  )$symbol
}

lfc_column <- function(df) {
  hit <- c("avg_log2FC", "avg_logFC", "log2FC", "logFC")
  hit <- hit[hit %in% names(df)]
  if (!length(hit)) stop("Cannot find a log-fold-change column in FindMarkers output", call. = FALSE)
  hit[[1L]]
}

hallmark_sets <- function() {
  hallmark <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  )
  if (!("db_version" %in% names(hallmark))) {
    stop("msigdbr Hallmark result does not report db_version", call. = FALSE)
  }
  releases <- unique(clean_character(hallmark$db_version))
  expected_release <- as.character(config$gene_sets$database_release)
  if (!identical(releases, expected_release)) {
    stop(
      "Legacy SI7 raw rebuild requires MSigDB ",
      expected_release,
      "; observed ",
      paste(releases, collapse = ","),
      call. = FALSE
    )
  }
  sets <- lapply(
    split(hallmark$gene_symbol, hallmark$gs_name),
    function(genes) sort(unique(clean_character(genes)))
  )
  membership <- unlist(
    lapply(
      sort(names(sets)),
      function(pathway) paste(pathway, sets[[pathway]], sep = "\t")
    ),
    use.names = FALSE
  )
  attr(sets, "database_release") <- expected_release
  attr(sets, "membership_sha256") <- digest::digest(
    paste(membership, collapse = "\n"),
    algo = "sha256",
    serialize = FALSE
  )
  sets
}

hallmark_label <- function(x) {
  tools::toTitleCase(tolower(gsub("_", " ", sub("^HALLMARK_", "", x))))
}

prepare_ora_markers <- function(de, lfc_col) {
  out <- de
  out$gene_symbol <- clean_gene_symbol(out$gene)
  out$gene_key <- out$gene_symbol
  out$lfc_value <- as.numeric(out[[lfc_col]])
  out$p_val_adj_num <- as.numeric(out$p_val_adj)
  out$delta_pct <- as.numeric(out$pct.1) - as.numeric(out$pct.2)
  out$abs_logfc <- abs(out$lfc_value)
  out$abs_delta_pct <- abs(out$delta_pct)
  out <- out[
    !is.na(out$gene_key) & is.finite(out$lfc_value) &
      is.finite(out$p_val_adj_num) & is.finite(out$delta_pct) &
      out$p_val_adj_num < 0.05 & out$abs_logfc >= 0.25 &
      out$abs_delta_pct >= 0.05 & out$lfc_value > 0,
    ,
    drop = FALSE
  ]
  out <- out[order(-out$lfc_value, out$p_val_adj_num), , drop = FALSE]
  out <- out[!duplicated(out$gene_key), , drop = FALSE]
  head(out, 100L)
}

run_ora <- function(query_df, universe, pathways) {
  query <- unique(query_df$gene_symbol)
  query <- query[!is.na(query)]
  results <- lapply(names(pathways), function(pathway) {
    genes <- unique(si7_case_sensitive_intersection(
      pathways[[pathway]],
      universe
    ))
    if (length(genes) < 15L || length(genes) > 500L) return(NULL)
    overlap <- si7_case_sensitive_intersection(query, genes)
    if (length(overlap) < 3L) return(NULL)
    p_value <- stats::phyper(
      length(overlap) - 1L,
      length(genes),
      length(universe) - length(genes),
      length(query),
      lower.tail = FALSE
    )
    data.frame(
      pathway = pathway,
      hallmark_label = hallmark_label(pathway),
      set_size = length(genes),
      query_size = length(query),
      overlap = length(overlap),
      p_value = p_value,
      overlap_genes = paste(sort(overlap), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, results)
  if (is.null(out) || !nrow(out)) return(data.frame())
  out$p_adj <- stats::p.adjust(out$p_value, method = "BH")
  out <- out[out$p_adj < 0.05, , drop = FALSE]
  if (!nrow(out)) return(out)
  overlap_stats <- lapply(strsplit(out$overlap_genes, ";", fixed = TRUE), function(genes) {
    match_rows <- query_df[query_df$gene_symbol %in% genes, , drop = FALSE]
    data.frame(
      expression_score_raw = mean(match_rows$abs_logfc, na.rm = TRUE),
      detection_score_raw = mean(as.numeric(match_rows$pct.1), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  score <- do.call(rbind, overlap_stats)
  scale01 <- function(x) {
    if (!length(x) || !any(is.finite(x))) return(rep(0, length(x)))
    limits <- range(x[is.finite(x)])
    if (diff(limits) == 0) return(rep(1, length(x)))
    (x - limits[[1L]]) / diff(limits)
  }
  out <- cbind(out, score)
  out$overlap_score_raw <- rowMeans(cbind(
    out$overlap / pmax(1, out$query_size),
    out$overlap / pmax(1, out$set_size)
  ))
  out$expression_score_scaled <- scale01(out$expression_score_raw)
  out$detection_score_scaled <- scale01(out$detection_score_raw)
  out$overlap_score_scaled <- scale01(out$overlap_score_raw)
  out$annotation_score <- rowMeans(out[, c(
    "expression_score_scaled", "detection_score_scaled", "overlap_score_scaled"
  )])
  out <- out[order(-out$annotation_score, out$p_adj, out$pathway), , drop = FALSE]
  out$annotation_rank <- seq_len(nrow(out))
  out
}

prepare_rank_stats <- function(de, lfc_col) {
  gene <- clean_gene_symbol(de$gene)
  lfc <- as.numeric(de[[lfc_col]])
  p_value <- as.numeric(de$p_val_adj)
  finite_positive <- p_value[is.finite(p_value) & p_value > 0]
  p_floor <- if (length(finite_positive)) max(min(finite_positive) * 0.1, 1e-300) else 1e-300
  p_value[!is.finite(p_value) | p_value <= 0] <- p_floor
  rank_metric <- sign(lfc) * (-log10(p_value) + abs(lfc) * 1e-6)
  ranked <- data.frame(gene = gene, rank_metric = rank_metric, lfc = lfc)
  ranked <- ranked[
    !is.na(ranked$gene) & is.finite(ranked$rank_metric) & ranked$rank_metric != 0,
    ,
    drop = FALSE
  ]
  ranked <- ranked[order(-abs(ranked$rank_metric), -abs(ranked$lfc), ranked$gene), , drop = FALSE]
  ranked <- ranked[!duplicated(ranked$gene), , drop = FALSE]
  ranked <- ranked[order(-ranked$rank_metric, ranked$gene), , drop = FALSE]
  stats <- ranked$rank_metric
  names(stats) <- ranked$gene
  stats
}

parse_memory_bytes <- function(value, default_unit = 1) {
  if (is.null(value) || !length(value)) return(NA_real_)
  value <- trimws(as.character(value[[1L]]))
  if (!nzchar(value) || tolower(value) == "max") return(NA_real_)
  match <- regexec("^([0-9]+[.]?[0-9]*)([KMGT]?)$", toupper(value))
  pieces <- regmatches(toupper(value), match)[[1L]]
  if (!length(pieces)) return(NA_real_)
  multiplier <- switch(
    pieces[[3L]],
    K = 1024,
    M = 1024^2,
    G = 1024^3,
    T = 1024^4,
    default_unit
  )
  as.numeric(pieces[[2L]]) * multiplier
}

detect_memory_limit_bytes <- function(detected_cpus) {
  candidates <- numeric(0)
  slurm_node <- parse_memory_bytes(Sys.getenv("SLURM_MEM_PER_NODE"), 1024^2)
  if (is.finite(slurm_node) && slurm_node > 0) {
    candidates <- c(candidates, slurm_node)
  }
  slurm_cpu <- parse_memory_bytes(Sys.getenv("SLURM_MEM_PER_CPU"), 1024^2)
  if (is.finite(slurm_cpu) && slurm_cpu > 0) {
    candidates <- c(candidates, slurm_cpu * detected_cpus)
  }
  for (path in c("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes")) {
    if (!file.exists(path)) next
    value <- tryCatch(readLines(path, n = 1L, warn = FALSE), error = function(e) NA_character_)
    parsed <- parse_memory_bytes(value, 1)
    if (is.finite(parsed) && parsed > 0 && parsed < 2^60) {
      candidates <- c(candidates, parsed)
    }
  }
  meminfo <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) character(0))
  memtotal <- grep("^MemTotal:", meminfo, value = TRUE)
  if (length(memtotal)) {
    parsed <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", memtotal[[1L]]))) * 1024
    if (is.finite(parsed) && parsed > 0) candidates <- c(candidates, parsed)
  }
  if (!length(candidates)) return(NA_real_)
  min(candidates)
}

deg_cache_paths <- function(work_dir, cluster_id) {
  si7_deg_cache_paths(work_dir, cluster_id)
}

read_reusable_cluster_deg <- function(
  work_dir,
  cluster_id,
  dependency_fingerprint
) {
  si7_read_reusable_cluster_deg(
    work_dir,
    cluster_id,
    dependency_fingerprint
  )
}

write_cluster_deg_cache <- function(
  de,
  work_dir,
  cluster_id,
  dependency_fingerprint
) {
  si7_write_cluster_deg_cache(
    de,
    work_dir,
    cluster_id,
    dependency_fingerprint
  )
}

find_cluster_markers_parallel <- function(
  cluster_index,
  cluster_levels,
  inner_workers,
  multicore_enabled,
  seurat_object,
  table_dir,
  dependency_fingerprint
) {
  cluster_id <- cluster_levels[[cluster_index]]
  assigned_inner_workers <- inner_workers[[cluster_index]]
  worker_previous_plan <- future::plan()
  on.exit(future::plan(worker_previous_plan), add = TRUE)
  if (multicore_enabled && assigned_inner_workers > 1L) {
    future::plan(future::multicore, workers = assigned_inner_workers)
  } else {
    future::plan(future::sequential)
  }
  message("FindMarkers: cluster ", cluster_id, " versus rest")
  de <- Seurat::FindMarkers(
    object = seurat_object,
    ident.1 = cluster_id,
    ident.2 = NULL,
    assay = "RNA",
    slot = "data",
    verbose = FALSE
  )
  de$gene <- rownames(de)
  de$gene_symbol <- clean_gene_symbol(de$gene)
  de$cluster <- cluster_id
  lfc_col <- lfc_column(de)
  de <- de[order(de$p_val_adj, -abs(de[[lfc_col]])), , drop = FALSE]
  rownames(de) <- NULL
  write_cluster_deg_cache(
    de,
    table_dir,
    cluster_id,
    dependency_fingerprint
  )
  list(cluster = cluster_id, de = de, reused = FALSE)
}

run_si7 <- function() {
  require_package("Seurat")
  require_package("future")
  require_package("future.apply")
  require_package("fgsea")
  require_package("msigdbr")
  observed_msigdbr <- as.character(
    utils::packageVersion("msigdbr")
  )
  expected_msigdbr <- as.character(config$gene_sets$package_version)
  if (!identical(observed_msigdbr, expected_msigdbr)) {
    stop(
      "Generated human-only SI7 rebuild requires msigdbr ",
      expected_msigdbr,
      "; observed ",
      observed_msigdbr,
      call. = FALSE
    )
  }
  previous_future_plan <- future::plan()
  previous_future_max_size <- getOption("future.globals.maxSize")
  on.exit(
    {
      future::plan(previous_future_plan)
      options(future.globals.maxSize = previous_future_max_size)
    },
    add = TRUE
  )
  options(future.globals.maxSize = 8 * 1024^3)
  set.seed(analysis_seed)
  message("Generating Supplementary Figure 7 from raw Seurat RDS: ", seurat_rds_path)
  future_global_limit_bytes <- max(
    8 * 1024^3,
    as.numeric(utils::object.size(object)) * 3 + 2 * 1024^3
  )
  options(future.globals.maxSize = future_global_limit_bytes)
  if (!(canonical_cluster_col %in% names(object@meta.data))) {
    stop("Raw Seurat object is missing cluster field: ", canonical_cluster_col, call. = FALSE)
  }
  if (ncol(object) != nrow(seurat)) {
    stop(
      "Raw Seurat RDS cell count differs from extracted metadata: ",
      ncol(object), " versus ", nrow(seurat),
      call. = FALSE
    )
  }
  if (!("RNA" %in% names(object@assays))) {
    stop("Raw Seurat object is missing RNA assay", call. = FALSE)
  }
  n_source_features <- nrow(object[["RNA"]])
  human_only <- si7_build_human_only_object(
    object,
    human_prefix = si7_human_prefix,
    mouse_prefix = si7_mouse_prefix
  )
  object <- human_only$object
  species_audit <- human_only$audit
  # Release the large mixed-species source and reconstruction temporaries
  # before FindMarkers allocates its cluster-vs-rest working vectors.
  source_environment <- environment(run_si7)
  if (!exists("object", envir = source_environment, inherits = FALSE) ||
      !inherits(
        get("object", envir = source_environment, inherits = FALSE),
        "Seurat"
      )) {
    stop("SI7 cannot release the enclosing source Seurat object",
         call. = FALSE)
  }
  rm(list = "object", envir = source_environment)
  rm(human_only)
  invisible(gc(full = TRUE))
  if (!identical(
        as.integer(species_audit$n_input_features),
        as.integer(n_source_features)
      ) ||
      !identical(
        as.integer(species_audit$n_human_features_retained),
        as.integer(nrow(object))
      ) ||
      !identical(
        as.integer(species_audit$n_input_features),
        as.integer(
          species_audit$n_human_features_retained +
            species_audit$n_mouse_features_excluded +
            species_audit$n_ambiguous_features
        )
      ) ||
      !identical(as.integer(species_audit$n_ambiguous_features), 0L)) {
    stop("SI7 feature-species audit is internally inconsistent",
         call. = FALSE)
  }
  write_tsv(
    species_audit,
    file.path(work_dir, "si_figure7_feature_species_audit.tsv")
  )
  Seurat::DefaultAssay(object) <- "RNA"
  rds_clusters <- sort_cluster_levels(object@meta.data[[canonical_cluster_col]])
  if (!identical(rds_clusters, cluster_levels)) {
    stop(
      "Raw RDS and Seurat metadata cluster levels differ: RDS=",
      paste(rds_clusters, collapse = ","), "; metadata=", paste(cluster_levels, collapse = ","),
      call. = FALSE
    )
  }
  Seurat::Idents(object) <- factor(
    as.character(object@meta.data[[canonical_cluster_col]]),
    levels = cluster_levels
  )
  available_cpus <- max(1L, as.integer(future::availableCores()))
  detected_cpus <- min(available_cpus, worker_budget)
  reserved_cpus <- min(2L, max(0L, detected_cpus - 1L))
  usable_cpus <- max(1L, detected_cpus - reserved_cpus)
  object_size_bytes <- as.numeric(utils::object.size(object))
  memory_limit_bytes <- detect_memory_limit_bytes(detected_cpus)
  memory_per_worker_bytes <- object_size_bytes * 1.5
  memory_worker_budget <- if (is.finite(memory_limit_bytes)) {
    max(
      1L,
      as.integer(floor(
        max(0, memory_limit_bytes * 0.90 - object_size_bytes) /
          max(memory_per_worker_bytes, 1)
      ))
    )
  } else {
    usable_cpus
  }
  execution_worker_budget <- max(
    1L,
    min(usable_cpus, memory_worker_budget)
  )
  multicore_enabled <- future::supportsMulticore() && execution_worker_budget > 1L
  cluster_workers <- if (multicore_enabled) {
    min(length(cluster_levels), execution_worker_budget)
  } else {
    1L
  }
  active_inner_workers <- rep(
    execution_worker_budget %/% cluster_workers,
    cluster_workers
  )
  remainder <- execution_worker_budget %% cluster_workers
  if (remainder > 0L) {
    active_inner_workers[seq_len(remainder)] <-
      active_inner_workers[seq_len(remainder)] + 1L
  }
  inner_workers <- rep(active_inner_workers, length.out = length(cluster_levels))
  names(inner_workers) <- cluster_levels
  message(
    "SI Figure 7 differential expression layout: detected_cpus=", detected_cpus,
    "; reserved_cpus=", reserved_cpus,
    "; usable_cpus=", usable_cpus,
    "; memory_limit_bytes=", format(memory_limit_bytes, scientific = FALSE),
    "; object_size_bytes=", format(object_size_bytes, scientific = FALSE),
    "; memory_per_worker_bytes=", format(memory_per_worker_bytes, scientific = FALSE),
    "; memory_worker_budget=", memory_worker_budget,
    "; execution_worker_budget=", execution_worker_budget,
    "; cluster_workers=", cluster_workers,
    "; inner_workers=", paste(inner_workers, collapse = ",")
  )

  deg_results <- lapply(
    cluster_levels,
    function(cluster_id) {
      read_reusable_cluster_deg(
        work_dir,
        cluster_id,
        work_dependency_fingerprint
      )
    }
  )
  missing_cluster_indices <- which(vapply(
    deg_results,
    is.null,
    logical(1L)
  ))
  if (length(missing_cluster_indices)) {
    message(
      "SI Figure 7 DEG cache: reusing ",
      length(cluster_levels) - length(missing_cluster_indices),
      "; computing ",
      length(missing_cluster_indices),
      " cluster(s)"
    )
  } else {
    message("SI Figure 7 DEG cache: reusing all clusters")
  }

  computed_results <- list()
  if (cluster_workers > 1L && length(missing_cluster_indices)) {
    future::plan(future::multicore, workers = cluster_workers)
    computed_results <- future.apply::future_lapply(
      missing_cluster_indices,
      find_cluster_markers_parallel,
      cluster_levels = cluster_levels,
      inner_workers = inner_workers,
      multicore_enabled = multicore_enabled,
      seurat_object = object,
      table_dir = work_dir,
      dependency_fingerprint = work_dependency_fingerprint,
      future.seed = TRUE,
      future.scheduling = 1
    )
  } else if (length(missing_cluster_indices)) {
    future::plan(future::sequential)
    computed_results <- lapply(
      missing_cluster_indices,
      find_cluster_markers_parallel,
      cluster_levels = cluster_levels,
      inner_workers = inner_workers,
      multicore_enabled = multicore_enabled,
      seurat_object = object,
      table_dir = work_dir,
      dependency_fingerprint = work_dependency_fingerprint
    )
  }
  if (length(missing_cluster_indices)) {
    deg_results[missing_cluster_indices] <- computed_results
  }
  deg_list <- stats::setNames(
    lapply(deg_results, function(result) result$de),
    vapply(deg_results, function(result) result$cluster, character(1))
  )
  universe <- sort(unique(unlist(lapply(deg_list, function(de) de$gene_symbol))))
  universe <- universe[!is.na(universe) & nzchar(universe)]
  write_csv(
    data.frame(gene_symbol = universe, stringsAsFactors = FALSE),
    file.path(work_dir, "si_figure7_hallmark_ORA_universe.csv")
  )
  pathways <- hallmark_sets()
  hallmark_release <- attr(pathways, "database_release")
  hallmark_membership_sha256 <- attr(pathways, "membership_sha256")
  ora_list <- list()
  gsea_list <- list()
  nperm_simple <- fgsea_nperm_simple
  for (cluster_id in cluster_levels) {
    de <- deg_list[[cluster_id]]
    lfc_col <- lfc_column(de)
    ora_markers <- prepare_ora_markers(de, lfc_col)
    write_csv(
      ora_markers,
      file.path(
        work_dir,
        paste0("si_figure7_cluster_", safe_cluster_stub(cluster_id), "_top100_up_ORA_input.csv")
      )
    )
    ora <- run_ora(ora_markers, universe, pathways)
    if (nrow(ora)) ora$cluster <- cluster_id
    ora_list[[cluster_id]] <- ora
    stats <- prepare_rank_stats(de, lfc_col)
    set.seed(analysis_seed + match(cluster_id, cluster_levels))
    gsea <- fgsea::fgseaMultilevel(
      pathways = pathways,
      stats = stats,
      minSize = 15L,
      maxSize = 500L,
      nPermSimple = nperm_simple,
      eps = 0
    )
    gsea <- as.data.frame(gsea, stringsAsFactors = FALSE)
    if ("leadingEdge" %in% names(gsea)) {
      gsea$leadingEdge <- vapply(
        gsea$leadingEdge,
        function(x) paste(as.character(x), collapse = ";"),
        character(1)
      )
    }
    gsea$hallmark_label <- hallmark_label(gsea$pathway)
    gsea$cluster <- cluster_id
    gsea <- gsea[order(gsea$padj, -abs(gsea$NES), gsea$pathway), , drop = FALSE]
    gsea_list[[cluster_id]] <- gsea
  }
  ora_nonempty <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, ora_list)
  ora_all <- if (length(ora_nonempty)) do.call(rbind, ora_nonempty) else data.frame()
  gsea_all <- do.call(rbind, gsea_list)
  if (is.null(ora_all) || !nrow(ora_all)) stop("No significant Hallmark ORA results", call. = FALSE)
  if (is.null(gsea_all) || !nrow(gsea_all)) stop("No Hallmark GSEA results", call. = FALSE)
  rownames(ora_all) <- NULL
  rownames(gsea_all) <- NULL
  write_csv(ora_all, file.path(work_dir, "si_figure7_cluster_Hallmark_ORA_all.csv"))
  write_csv(gsea_all, file.path(work_dir, "si_figure7_cluster_Hallmark_GSEA_all.csv"))

  ora_top <- aggregate(annotation_score ~ hallmark_label, data = ora_all, FUN = max)
  ora_top <- ora_top[order(-ora_top$annotation_score, ora_top$hallmark_label), , drop = FALSE]
  ora_pathways <- head(ora_top$hallmark_label, 20L)
  ora_matrix <- matrix(
    0,
    nrow = length(ora_pathways),
    ncol = length(cluster_levels),
    dimnames = list(ora_pathways, cluster_levels)
  )
  for (i in seq_len(nrow(ora_all))) {
    row_id <- match(ora_all$hallmark_label[[i]], rownames(ora_matrix))
    col_id <- match(ora_all$cluster[[i]], colnames(ora_matrix))
    if (!is.na(row_id) && !is.na(col_id)) {
      ora_matrix[row_id, col_id] <- max(
        ora_matrix[row_id, col_id], ora_all$annotation_score[[i]], na.rm = TRUE
      )
    }
  }
  gsea_heatmap <- gsea_all[is.finite(gsea_all$NES), , drop = FALSE]
  if (!nrow(gsea_heatmap)) {
    stop("No finite Hallmark GSEA NES values are available for the heatmap", call. = FALSE)
  }
  gsea_heatmap$abs_NES <- abs(gsea_heatmap$NES)
  gsea_top <- aggregate(abs_NES ~ hallmark_label, data = gsea_heatmap, FUN = max)
  names(gsea_top)[2L] <- "best_abs_nes"
  gsea_top <- gsea_top[order(-gsea_top$best_abs_nes, gsea_top$hallmark_label), , drop = FALSE]
  gsea_pathways <- head(gsea_top$hallmark_label, 20L)
  gsea_matrix <- matrix(
    0,
    nrow = length(gsea_pathways),
    ncol = length(cluster_levels),
    dimnames = list(gsea_pathways, cluster_levels)
  )
  for (i in seq_len(nrow(gsea_heatmap))) {
    row_id <- match(gsea_heatmap$hallmark_label[[i]], rownames(gsea_matrix))
    col_id <- match(gsea_heatmap$cluster[[i]], colnames(gsea_matrix))
    if (!is.na(row_id) && !is.na(col_id)) {
      current <- gsea_matrix[row_id, col_id]
      candidate <- gsea_heatmap$NES[[i]]
      if (abs(candidate) >= abs(current)) gsea_matrix[row_id, col_id] <- candidate
    }
  }
  list(
    ora_matrix = ora_matrix,
    gsea_matrix = gsea_matrix,
    ora_output = data.frame(
      pathway = rownames(ora_matrix),
      ora_matrix,
      check.names = FALSE
    ),
    gsea_output = data.frame(
      pathway = rownames(gsea_matrix),
      gsea_matrix,
      check.names = FALSE
    ),
    nperm_simple = nperm_simple,
    requested_worker_budget = worker_budget,
    available_cpus = available_cpus,
    detected_cpus = detected_cpus,
    reserved_cpus = reserved_cpus,
    usable_cpus = usable_cpus,
    cluster_workers = cluster_workers,
    inner_workers = paste(inner_workers, collapse = ","),
    memory_limit_bytes = memory_limit_bytes,
    object_size_bytes = object_size_bytes,
    memory_per_worker_bytes = memory_per_worker_bytes,
    memory_worker_budget = memory_worker_budget,
    execution_worker_budget = execution_worker_budget,
    deg_clusters_reused = sum(vapply(
      deg_results,
      function(result) isTRUE(result$reused),
      logical(1L)
    )),
    deg_clusters_computed = sum(vapply(
      deg_results,
      function(result) !isTRUE(result$reused),
      logical(1L)
    )),
    future_global_limit_bytes = future_global_limit_bytes,
    hallmark_release = hallmark_release,
    hallmark_membership_sha256 = hallmark_membership_sha256,
    n_object_cells = ncol(object),
    n_source_features = n_source_features,
    n_genes = nrow(object),
    feature_species_audit = species_audit
  )
}

si7_result <- run_si7()
si_science_marker("END:si7")

si_science_marker("BEGIN:final_tables")
group_summary_output <- group_summary
group_summary_output$initial_ploidy <- as.character(
  group_summary_output$initial_ploidy
)
group_summary_output$dose <- as.character(group_summary_output$dose)
group_summary_output$cluster_final <- as.character(
  group_summary_output$cluster_final
)
endpoint_audit <- data.frame(
  cell = seurat$cell,
  context = seurat$context,
  initial_ploidy = seurat$initial_ploidy,
  endpoint_file = seurat$endpoint_file,
  endpoint_cell_id = seurat$endpoint_cell_id,
  endpoint_ploidy = seurat$endpoint_ploidy,
  matched = is.finite(seurat$endpoint_ploidy),
  stringsAsFactors = FALSE
)
si_science_marker("END:final_tables")

si_science_marker("BEGIN:plot_table_outputs")
write_csv(canonical_cells, file.path(table_dir, "si_figures_cell_metadata.csv"))
write_tsv(cluster_key, file.path(table_dir, "si_figures_cluster_key.tsv"))
write_csv(
  s4_context,
  file.path(table_dir, "si_figure4_cluster_context_composition.csv")
)
write_csv(
  s4_ploidy,
  file.path(table_dir, "si_figure4_cluster_initial_ploidy_composition.csv")
)
write_csv(
  composition_table,
  file.path(table_dir, "si_figure5_cluster_composition_by_mouse.csv")
)
write_csv(
  group_summary_output,
  file.path(
    table_dir,
    "si_figure5_mouse_weighted_composition_by_initial_ploidy_dose.csv"
  )
)
write_csv(
  s5_dose_composition,
  file.path(table_dir, "si_figure5_cluster_dose_composition.csv")
)
write_csv(
  s5_ploidy_composition,
  file.path(table_dir, "si_figure5_cluster_initial_ploidy_composition.csv")
)
write_csv(
  endpoint_audit,
  file.path(table_dir, "si_figure6_endpoint_ploidy_join_audit.csv")
)
write_tsv(
  si7_result$ora_output,
  file.path(
    table_dir,
    paste0(
      "si_figure7_cluster_Hallmark_ORA_annotation_score_",
      "heatmap_top20_matrix.tsv"
    )
  )
)
write_tsv(
  si7_result$gsea_output,
  file.path(
    table_dir,
    "si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv"
  )
)

table_files <- c(
  "si_figure4_cluster_context_composition.csv",
  "si_figure4_cluster_initial_ploidy_composition.csv",
  "si_figure5_cluster_composition_by_mouse.csv",
  "si_figure5_cluster_dose_composition.csv",
  "si_figure5_cluster_initial_ploidy_composition.csv",
  "si_figure5_mouse_weighted_composition_by_initial_ploidy_dose.csv",
  "si_figure6_endpoint_ploidy_join_audit.csv",
  "si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv",
  "si_figure7_cluster_Hallmark_ORA_annotation_score_heatmap_top20_matrix.tsv",
  "si_figures_cell_metadata.csv",
  "si_figures_cluster_key.tsv"
)
si_science_marker("END:plot_table_outputs")
table_paths <- file.path(table_dir, table_files)
missing_tables <- table_files[!file.exists(table_paths)]
if (length(missing_tables)) {
  stop(
    "Raw SI table build is incomplete: ",
    paste(missing_tables, collapse = ", "),
    call. = FALSE
  )
}
unexpected_tables <- setdiff(
  list.files(table_dir),
  c(table_files, "manifest.tsv")
)
if (length(unexpected_tables)) {
  stop(
    "Raw SI table directory contains non-plot-facing files: ",
    paste(unexpected_tables, collapse = ", "),
    call. = FALSE
  )
}

source_revision <- tryCatch(
  system2(
    "git",
    c("-C", shQuote(repo_root), "rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
  )[[1L]],
  error = function(error) "unavailable"
)
table_manifest <- data.frame(
  filename = table_files,
  bytes = as.numeric(file.info(table_paths)$size),
  sha256 = vapply(table_paths, file_sha256, character(1)),
  source_revision = paste0("raw-generated-human-only@", source_revision),
  notes = ifelse(
    grepl("^si_figure7_cluster_Hallmark", table_files),
    paste(
      "generated human-only GRCh policy: exact GRCh38 features retained before",
      "fresh normalization, DE, ORA, and GSEA; not approved for canonical publication"
    ),
    "raw-rebuilt plot-facing table; not approved for canonical publication"
  ),
  stringsAsFactors = FALSE
)
write_tsv(table_manifest, file.path(table_dir, "manifest.tsv"))

  portable_locator <- function(path) {
  normalized <- normalizePath(path, mustWork = TRUE)
  root_prefix <- paste0(
    normalizePath(repo_root, mustWork = TRUE),
    .Platform$file.sep
  )
  if (startsWith(normalized, root_prefix)) {
    substring(normalized, nchar(root_prefix) + 1L)
  } else {
      paste0("external:", basename(normalized))
    }
  }
  dependency_locators <- stats::setNames(
    paste0("contract:", names(work_dependency_values)),
    names(work_dependency_values)
  )
  dependency_bytes <- stats::setNames(
    rep(0, length(work_dependency_values)),
    names(work_dependency_values)
  )
  dependency_locators[["versioned_endpoint_ploidy"]] <-
    portable_locator(all_ploidy_path)
  dependency_bytes[["versioned_endpoint_ploidy"]] <-
    as.numeric(file.info(all_ploidy_path)$size)
  dependency_locators[["source_seurat_rds"]] <-
    portable_locator(seurat_rds_path)
  dependency_bytes[["source_seurat_rds"]] <-
    as.numeric(file.info(seurat_rds_path)$size)
  computational_input_manifest <- data.frame(
    role = names(work_dependency_values),
    locator = unname(dependency_locators[names(work_dependency_values)]),
    sha256 = unname(work_dependency_values),
    bytes = unname(dependency_bytes[names(work_dependency_values)]),
    stringsAsFactors = FALSE
  )
audit_paths <- c(
  script_path,
  si7_feature_species_helper_path,
  environment_validator_path,
    environment_lock_path,
    config_path,
    seurat_selection_path
  )
audit_roles <- c(
  "audit_si_table_builder",
  "audit_feature_species_helper",
  "audit_environment_validator",
    "audit_environment_lock",
    "audit_figure7_config",
    "audit_seurat_selection"
  )
audit_sha256 <- c(
  builder_sha256,
  file_sha256(si7_feature_species_helper_path),
  environment_validator_sha256,
    environment_lock_sha256,
    config_sha256,
    file_sha256(seurat_selection_path)
  )
  if (!is.null(seurat_validation)) {
    audit_paths <- c(
      audit_paths,
      seurat_generator_path,
      seurat_upstream_helper_path,
      seurat_validation$final_stage_manifest,
      seurat_validation$reconstruction_manifest
    )
    audit_roles <- c(
      audit_roles,
      "audit_seurat_upstream_generator",
      "audit_seurat_upstream_helper",
      "audit_seurat_final_stage_manifest",
      "audit_seurat_reconstruction_manifest"
    )
    audit_sha256 <- c(
      audit_sha256,
      file_sha256(seurat_generator_path),
      file_sha256(seurat_upstream_helper_path),
      seurat_validation$final_stage_manifest_sha256,
      seurat_validation$reconstruction_manifest_sha256
    )
  }
  if (!is.na(scvelo_metrics_path)) {
  audit_paths <- c(audit_paths, scvelo_metrics_path)
  audit_roles <- c(audit_roles, "audit_optional_scvelo")
  audit_sha256 <- c(audit_sha256, scvelo_metrics_sha256)
}
audit_input_manifest <- data.frame(
  role = audit_roles,
  locator = vapply(audit_paths, portable_locator, character(1)),
  sha256 = audit_sha256,
  bytes = as.numeric(file.info(audit_paths)$size),
  stringsAsFactors = FALSE
  )
  input_manifest <- rbind(
    computational_input_manifest,
    audit_input_manifest
  )
write_tsv(input_manifest, file.path(metadata_dir, "input_manifest.tsv"))

work_paths <- list.files(work_dir, full.names = TRUE)
work_manifest <- data.frame(
  filename = basename(work_paths),
  bytes = as.numeric(file.info(work_paths)$size),
  sha256 = vapply(work_paths, file_sha256, character(1)),
  dependency_sha256 = paste(
    sort(input_manifest$sha256[!startsWith(input_manifest$role, "audit_")]),
    collapse = ";"
  ),
  stringsAsFactors = FALSE
)
write_tsv(work_manifest, file.path(metadata_dir, "work_cache_manifest.tsv"))

run_config <- data.frame(
  key = c(
    "schema_version",
    "module",
    "analysis_mode",
    "reviewed_ploidy_field",
    "reviewed_context_field",
    "ploidy_context_policy",
    "output_table_count",
    "figures_supported",
    "si7_species_policy",
    "si7_species_policy_id",
    "si7_human_feature_prefix",
    "si7_mouse_feature_prefix",
    "si7_canonical_publication_allowed",
    "si7_input_features",
    "si7_human_features_retained",
    "si7_mouse_features_excluded",
    "si7_ambiguous_features",
    "si7_normalization_method",
    "si7_normalization_scale_factor",
    "si7_source_data_layer_reused",
    "si7_feature_species_audit_sha256",
    "feature_species_helper_sha256",
    "si7_gene_set_source",
    "si7_gene_set_release",
    "si7_gene_set_membership_sha256",
    "fgsea_nperm_simple",
    "analysis_seed",
    "seurat_source_kind",
    "si_raw_scientific_code_contract_sha256",
    "environment_lock_sha256",
    "environment_validator_sha256",
    "si_raw_config_contract_sha256",
    "si_raw_runtime_contract_sha256",
    "locked_r_versions",
    "r_runtime_version",
    "r_platform",
    "r_blas",
    "work_dependency_fingerprint",
    "requested_worker_budget",
    "available_cpus",
    "deg_clusters_reused",
    "deg_clusters_computed",
    "cluster_order",
    "seurat_cells",
    "tumor_cells",
    "seurat_rds_cells",
    "seurat_rds_features",
    "si7_analyzed_features",
    "audit_optional_scvelo_sha256"
  ),
  value = c(
    "1",
    "si_figures_raw_table_build",
    paste0(seurat_source_kind, " source reconstruction"),
    reviewed_ploidy_col,
    reviewed_context_col,
    paste(
      "configured RDS fields are authoritative and sample-invariant;",
      "ID-derived values are contradiction checks only"
    ),
    as.character(length(table_files)),
    "4,5,6,7",
    figure7_human_feature_policy_description(),
    as.character(si7_result$feature_species_audit$policy_id),
    as.character(si7_result$feature_species_audit$human_prefix),
    as.character(si7_result$feature_species_audit$mouse_prefix),
    "false",
    as.character(si7_result$feature_species_audit$n_input_features),
    as.character(
      si7_result$feature_species_audit$n_human_features_retained
    ),
    as.character(
      si7_result$feature_species_audit$n_mouse_features_excluded
    ),
    as.character(si7_result$feature_species_audit$n_ambiguous_features),
    as.character(si7_result$feature_species_audit$normalization_method),
    as.character(
      si7_result$feature_species_audit$normalization_scale_factor
    ),
    tolower(as.character(
      si7_result$feature_species_audit$source_data_layer_reused
    )),
    file_sha256(file.path(
      work_dir,
      "si_figure7_feature_species_audit.tsv"
    )),
    file_sha256(si7_feature_species_helper_path),
    "live msigdbr Homo sapiens Hallmark; case-sensitive symbol matching",
    as.character(si7_result$hallmark_release),
    as.character(si7_result$hallmark_membership_sha256),
    as.character(si7_result$nperm_simple),
    as.character(analysis_seed),
    seurat_source_kind,
    si_raw_scientific_code_contract_sha256,
    environment_lock_sha256,
    environment_validator_sha256,
    si_raw_config_contract_sha256,
    si_raw_runtime_contract_sha256,
    paste(
      paste(names(locked_r_versions), locked_r_versions, sep = "="),
      collapse = ";"
    ),
    unname(r_runtime_provenance[["audit_r_runtime_version"]]),
    unname(r_runtime_provenance[["audit_r_platform"]]),
    unname(r_runtime_provenance[["audit_r_blas"]]),
    work_dependency_fingerprint,
    as.character(si7_result$requested_worker_budget),
    as.character(si7_result$available_cpus),
    as.character(si7_result$deg_clusters_reused),
    as.character(si7_result$deg_clusters_computed),
    paste(cluster_levels, collapse = ","),
    as.character(nrow(seurat)),
    as.character(nrow(tumor)),
    as.character(si7_result$n_object_cells),
    as.character(si7_result$n_source_features),
    as.character(si7_result$n_genes),
    if (is.na(scvelo_metrics_path)) {
      "not_supplied"
    } else {
      scvelo_metrics_sha256
    }
  ),
  stringsAsFactors = FALSE
)
write_tsv(run_config, file.path(metadata_dir, "run_config.tsv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(metadata_dir, "sessionInfo.txt")
)

validator <- file.path(
  repo_root,
  "Code",
  "tools",
  "validate_si_figures_table_cache.py"
)
validation_status <- system2(
  "python3",
  c(
    shQuote(validator),
    "--cache-dir",
    shQuote(table_dir),
    "--si7-policy",
    "generated-human-only"
  )
)
if (!identical(validation_status, 0L)) {
  stop("Raw SI plot-facing tables failed validation", call. = FALSE)
}

message(
  "Built exactly 11 raw SI plot-facing tables with generated human-only ",
  "SI7 policy: ",
  table_dir
)
}

if (sys.nframe() == 0L) main()
