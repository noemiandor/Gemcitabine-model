# Shared feature-species boundary for Figure 7F and Supplementary Figure 7.
#
# The xenograft RNA assay is expected to carry an explicit reference prefix on
# every feature. Only the exact configured GRCh prefix identifies human-tumor
# measurements; the exact configured GRCm prefix identifies host-mouse
# measurements. An unprefixed or otherwise unsupported feature is not guessed
# from symbol capitalization and therefore fails closed.

figure7_human_feature_policy_id <- function() {
  "grch_human_tumor_only_v2"
}

figure7_human_feature_policy_description <- function() {
  paste(
    "retain only explicitly GRCh-prefixed human-tumor features;",
    "exclude GRCm-prefixed mouse features before expression filtering,",
    "normalization, symbol cleanup, duplicate resolution, differential",
    "expression, model fitting, ORA, and GSEA; reject ambiguous or",
    "unprefixed feature names"
  )
}

figure7_classify_feature_species <- function(
  features,
  human_prefix = "GRCh38-",
  mouse_prefix = "GRCm39-"
) {
  if (length(human_prefix) != 1L || is.na(human_prefix) ||
      !nzchar(human_prefix) ||
      length(mouse_prefix) != 1L || is.na(mouse_prefix) ||
      !nzchar(mouse_prefix) ||
      identical(human_prefix, mouse_prefix)) {
    stop(
      "Feature-species policy requires distinct, nonempty human/mouse prefixes",
      call. = FALSE
    )
  }
  feature <- as.character(features)
  if (!length(feature)) {
    stop("Feature-species policy received no features", call. = FALSE)
  }
  if (anyNA(feature) || any(!nzchar(feature))) {
    stop(
      "Feature-species policy requires nonempty, nonmissing feature names",
      call. = FALSE
    )
  }
  human <- startsWith(feature, human_prefix)
  mouse <- startsWith(feature, mouse_prefix)
  species <- rep("ambiguous", length(feature))
  species[human & !mouse] <- "human"
  species[mouse & !human] <- "mouse"

  symbol <- rep(NA_character_, length(feature))
  recognized <- species != "ambiguous"
  symbol[human] <- substring(feature[human], nchar(human_prefix) + 1L)
  symbol[mouse] <- substring(feature[mouse], nchar(mouse_prefix) + 1L)
  symbol[recognized] <- sub("\\.[0-9]+$", "", symbol[recognized])
  empty_symbol <- recognized & (!nzchar(symbol) | is.na(symbol))
  species[empty_symbol] <- "ambiguous"
  symbol[species == "ambiguous"] <- NA_character_

  data.frame(
    feature = feature,
    species = species,
    symbol = symbol,
    stringsAsFactors = FALSE
  )
}

figure7_select_human_features <- function(
  features,
  analysis = "Figure 7F/SI7",
  human_prefix = "GRCh38-",
  mouse_prefix = "GRCm39-"
) {
  annotation <- figure7_classify_feature_species(
    features,
    human_prefix = human_prefix,
    mouse_prefix = mouse_prefix
  )
  ambiguous <- annotation$feature[annotation$species == "ambiguous"]
  if (length(ambiguous)) {
    examples <- paste(head(ambiguous, 8L), collapse = ", ")
    stop(
      analysis,
      " requires an explicit GRCh/GRCm prefix on every RNA feature; ",
      length(ambiguous),
      " ambiguous or unprefixed feature(s) include: ",
      examples,
      call. = FALSE
    )
  }
  keep <- annotation$species == "human"
  if (!any(keep)) {
    stop(
      analysis,
      " has no explicitly GRCh-prefixed human features",
      call. = FALSE
    )
  }
  audit <- data.frame(
    policy_id = figure7_human_feature_policy_id(),
    policy = figure7_human_feature_policy_description(),
    human_prefix = human_prefix,
    mouse_prefix = mouse_prefix,
    n_input_features = nrow(annotation),
    n_human_features_retained = sum(keep),
    n_mouse_features_excluded = sum(annotation$species == "mouse"),
    n_ambiguous_features = sum(annotation$species == "ambiguous"),
    stringsAsFactors = FALSE
  )
  list(
    keep = keep,
    features = annotation$feature[keep],
    symbols = annotation$symbol[keep],
    annotation = annotation,
    audit = audit
  )
}

figure7_filter_human_feature_matrix <- function(
  matrix,
  analysis = "Figure 7F/SI7",
  human_prefix = "GRCh38-",
  mouse_prefix = "GRCm39-"
) {
  features <- rownames(matrix)
  if (is.null(features)) {
    stop(analysis, " expression matrix has no feature names", call. = FALSE)
  }
  selection <- figure7_select_human_features(
    features,
    analysis,
    human_prefix = human_prefix,
    mouse_prefix = mouse_prefix
  )
  filtered <- matrix[selection$keep, , drop = FALSE]
  if (!identical(rownames(filtered), selection$features)) {
    stop(analysis, " changed feature order at the human-only boundary",
         call. = FALSE)
  }
  figure7_assert_human_feature_names(
    rownames(filtered),
    analysis,
    human_prefix = human_prefix,
    mouse_prefix = mouse_prefix
  )
  list(matrix = filtered, selection = selection, audit = selection$audit)
}

figure7_assert_human_feature_names <- function(
  features,
  analysis = "Figure 7F/SI7",
  human_prefix = "GRCh38-",
  mouse_prefix = "GRCm39-"
) {
  selection <- figure7_select_human_features(
    features,
    analysis,
    human_prefix = human_prefix,
    mouse_prefix = mouse_prefix
  )
  if (!all(selection$keep)) {
    stop(
      analysis,
      " still contains nonhuman features after the GRCh-only boundary",
      call. = FALSE
    )
  }
  invisible(selection)
}
