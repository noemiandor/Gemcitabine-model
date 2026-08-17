#!/usr/bin/env Rscript

# Reconstruct the reviewed endpoint-tumor FlowJo hierarchy from raw FCS events.
#
# The relevant FlowJo channels are linear and the embedded spillover matrix is
# the identity.  This script therefore keeps recorded event intensities on the
# common raw detector scale and applies each sample's workspace-defined geometry
# in parent-to-child order.  Frozen FlowJo counts remain the reference; replay
# counts are reported alongside them rather than substituted for them.

suppressPackageStartupMessages({
  library(flowCore)
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(xml2)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) return(y)
  if (length(x) == 1L && (is.na(x) || (is.character(x) && !nzchar(x)))) return(y)
  x
}

abort <- function(...) {
  stop(paste0(...), call. = FALSE)
}

parse_cli <- function(args) {
  defaults <- list(
    expected_samples = 16L,
    min_human_cells = 1000L,
    dna_axis_min = 25000,
    dna_axis_max = 140000,
    dna_bin_width = 1000,
    representative_origin = "4N",
    max_representative_points = 10000L
  )
  value <- defaults
  required <- c("crosswalk", "workspace", "fcs_dir", "frozen_table", "output_dir")
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (token %in% c("-h", "--help")) {
      cat(paste(
        "Usage:",
        "  reconstruct_endpoint_flow.R --crosswalk FILE --workspace FILE",
        "    --fcs-dir DIR --frozen-table FILE --output-dir DIR [options]",
        "",
        "Options:",
        "  --expected-samples N             default: 16",
        "  --min-human-cells N              default: 1000",
        "  --dna-axis-min X                 default: 25000 raw a.u.",
        "  --dna-axis-max X                 default: 140000 raw a.u.",
        "  --dna-bin-width X                default: 1000 raw a.u.",
        "  --representative-origin 2N|4N    default: 4N",
        "  --max-representative-points N    default: 10000",
        sep = "\n"
      ))
      return(structure(value, help = TRUE))
    }
    if (!startsWith(token, "--")) abort("Unexpected argument: ", token)
    if (i == length(args)) abort("Missing value for ", token)
    key <- gsub("-", "_", substring(token, 3L), fixed = TRUE)
    if (!key %in% c(names(defaults), required)) abort("Unknown option: ", token)
    value[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  missing <- required[!vapply(required, function(k) nzchar(value[[k]] %||% ""), logical(1))]
  if (length(missing)) abort("Missing required option(s): --", paste(gsub("_", "-", missing), collapse = ", --"))
  integer_keys <- c("expected_samples", "min_human_cells", "max_representative_points")
  numeric_keys <- c("dna_axis_min", "dna_axis_max", "dna_bin_width")
  for (key in integer_keys) {
    parsed <- suppressWarnings(as.integer(value[[key]]))
    if (is.na(parsed) || parsed < 0L) abort("--", gsub("_", "-", key), " must be a nonnegative integer")
    value[[key]] <- parsed
  }
  for (key in numeric_keys) {
    parsed <- suppressWarnings(as.numeric(value[[key]]))
    if (!is.finite(parsed)) abort("--", gsub("_", "-", key), " must be finite")
    value[[key]] <- parsed
  }
  if (value$dna_axis_max <= value$dna_axis_min || value$dna_bin_width <= 0) {
    abort("DNA axis bounds and bin width are invalid")
  }
  span_bins <- (value$dna_axis_max - value$dna_axis_min) / value$dna_bin_width
  if (abs(span_bins - round(span_bins)) > 1e-9) abort("DNA axis span must be an integer number of bins")
  if (!value$representative_origin %in% c("2N", "4N")) abort("--representative-origin must be 2N or 4N")
  if (value$max_representative_points < 100L) abort("--max-representative-points must be at least 100")
  value
}

read_tsv <- function(path, label) {
  if (!file.exists(path)) abort("Missing ", label, ": ", path)
  tryCatch(
    read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
      check.names = FALSE, quote = "", comment.char = "", na.strings = "NA"),
    error = function(e) abort("Cannot read ", label, " ", path, ": ", conditionMessage(e))
  )
}

assert_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing)) abort(label, " is missing column(s): ", paste(missing, collapse = ", "))
}

write_tsv_atomic <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  old_options <- options(digits = 17)
  on.exit(options(old_options), add = TRUE)
  write.table(data, temporary, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "NA", eol = "\n")
  if (!file.rename(temporary, path)) abort("Cannot atomically replace output: ", path)
}

direct_children <- function(node, local_name) {
  xml_find_all(node, sprintf("./*[local-name()='%s']", local_name))
}

one_node <- function(nodes, description) {
  if (length(nodes) != 1L) abort("Expected one ", description, "; found ", length(nodes))
  nodes[[1L]]
}

direct_population <- function(parent, population_name, context) {
  subpops <- one_node(direct_children(parent, "Subpopulations"), paste0(context, "/Subpopulations"))
  populations <- direct_children(subpops, "Population")
  matches <- populations[xml_attr(populations, "name") == population_name]
  one_node(matches, paste0(context, "/", population_name))
}

population_count <- function(population, context) {
  value <- suppressWarnings(as.integer(xml_attr(population, "count")))
  if (length(value) != 1L || is.na(value) || value < 0L) abort("Invalid workspace count for ", context)
  value
}

sample_node_for_name <- function(sample_nodes, sample_name) {
  matches <- sample_nodes[xml_attr(sample_nodes, "name") == sample_name]
  one_node(matches, paste0("workspace SampleNode named ", sample_name))
}

sample_container_for_node <- function(sample_node) {
  parent <- xml_parent(sample_node)
  if (xml_name(parent) != "Sample") abort("SampleNode is not directly contained by Sample")
  parent
}

parse_gate <- function(population, population_name, population_role, mouse_id) {
  wrapper <- one_node(direct_children(population, "Gate"), paste0(mouse_id, "/", population_name, "/Gate"))
  gate_nodes <- xml_children(wrapper)
  if (length(gate_nodes) != 1L) abort("Expected one gate definition for ", mouse_id, "/", population_name)
  gate_node <- gate_nodes[[1L]]
  gate_type <- xml_name(gate_node)
  if (!gate_type %in% c("PolygonGate", "RectangleGate")) {
    abort("Unsupported gate type ", gate_type, " for ", mouse_id, "/", population_name)
  }
  events_inside <- xml_attr(gate_node, "eventsInside") %||% "1"
  if (events_inside != "1") abort("Only eventsInside=1 gates are supported: ", mouse_id, "/", population_name)
  dimensions <- direct_children(gate_node, "dimension")
  if (!length(dimensions)) abort("Gate has no dimensions: ", mouse_id, "/", population_name)
  parameter_names <- vapply(dimensions, function(dimension) {
    parameter <- one_node(direct_children(dimension, "fcs-dimension"), "gate fcs-dimension")
    xml_attr(parameter, "name") %||% abort("Gate dimension lacks a parameter name")
  }, character(1))
  if (anyDuplicated(parameter_names)) abort("Gate repeats a parameter: ", mouse_id, "/", population_name)

  parsed <- list(
    mouse_id = mouse_id,
    population_name = population_name,
    population_role = population_role,
    gate_type = gate_type,
    events_inside = TRUE,
    parameters = parameter_names,
    gate_resolution = suppressWarnings(as.integer(xml_attr(gate_node, "gateResolution") %||% NA_character_))
  )
  if (gate_type == "PolygonGate") {
    vertices <- direct_children(gate_node, "vertex")
    if (length(parameter_names) != 2L || length(vertices) < 3L) {
      abort("Polygon gate must have two dimensions and at least three vertices: ", mouse_id, "/", population_name)
    }
    coordinates <- t(vapply(vertices, function(vertex) {
      coordinate_nodes <- direct_children(vertex, "coordinate")
      values <- suppressWarnings(as.numeric(xml_attr(coordinate_nodes, "value")))
      if (length(values) != 2L || any(!is.finite(values))) abort("Invalid polygon vertex for ", mouse_id, "/", population_name)
      values
    }, numeric(2)))
    colnames(coordinates) <- parameter_names
    parsed$coordinates <- coordinates
  } else {
    bounds <- lapply(seq_along(dimensions), function(index) {
      minimum <- suppressWarnings(as.numeric(xml_attr(dimensions[[index]], "min") %||% "-Inf"))
      maximum <- suppressWarnings(as.numeric(xml_attr(dimensions[[index]], "max") %||% "Inf"))
      if (is.na(minimum) || is.na(maximum) || minimum >= maximum) abort("Invalid rectangle bounds for ", mouse_id, "/", population_name)
      c(min = minimum, max = maximum)
    })
    names(bounds) <- parameter_names
    parsed$bounds <- bounds
  }
  parsed
}

gate_geometry_rows <- function(gate) {
  common <- data.frame(
    mouse_id = gate$mouse_id,
    population_role = gate$population_role,
    population_name = gate$population_name,
    gate_type = gate$gate_type,
    x_parameter = gate$parameters[[1L]],
    y_parameter = if (length(gate$parameters) > 1L) gate$parameters[[2L]] else NA_character_,
    events_inside = TRUE,
    gate_resolution = gate$gate_resolution,
    stringsAsFactors = FALSE
  )
  if (gate$gate_type == "PolygonGate") {
    data.frame(common[rep(1L, nrow(gate$coordinates)), , drop = FALSE],
      vertex_index = seq_len(nrow(gate$coordinates)),
      x = gate$coordinates[, 1L], y = gate$coordinates[, 2L],
      x_min = NA_real_, x_max = NA_real_, y_min = NA_real_, y_max = NA_real_,
      stringsAsFactors = FALSE)
  } else {
    x_bounds <- gate$bounds[[1L]]
    y_bounds <- if (length(gate$bounds) > 1L) gate$bounds[[2L]] else c(min = NA_real_, max = NA_real_)
    data.frame(common, vertex_index = NA_integer_, x = NA_real_, y = NA_real_,
      x_min = unname(x_bounds[["min"]]), x_max = unname(x_bounds[["max"]]),
      y_min = unname(y_bounds[["min"]]), y_max = unname(y_bounds[["max"]]),
      stringsAsFactors = FALSE)
  }
}

resolve_parameter <- function(workspace_name, observed_names, context) {
  slash_name <- gsub("([0-9]+)_([0-9]+)", "\\1/\\2", workspace_name, perl = TRUE)
  candidates <- unique(c(workspace_name, slash_name))
  matches <- observed_names[observed_names %in% candidates]
  if (length(matches) != 1L) {
    abort("Cannot map workspace parameter ", workspace_name, " to one FCS channel in ", context)
  }
  matches[[1L]]
}

validate_linear_transforms <- function(sample_container, parameters, context) {
  transformations <- one_node(direct_children(sample_container, "Transformations"), paste0(context, "/Transformations"))
  transform_nodes <- xml_children(transformations)
  for (parameter_name in unique(parameters)) {
    matches <- transform_nodes[vapply(transform_nodes, function(node) {
      parameter <- direct_children(node, "parameter")
      length(parameter) == 1L && identical(xml_attr(parameter[[1L]], "name"), parameter_name)
    }, logical(1))]
    transform <- one_node(matches, paste0(context, " transform for ", parameter_name))
    if (xml_name(transform) != "linear") abort("Relevant channel is not linear in ", context, ": ", parameter_name)
    minimum <- as.numeric(xml_attr(transform, "minRange"))
    maximum <- as.numeric(xml_attr(transform, "maxRange"))
    gain <- as.numeric(xml_attr(transform, "gain"))
    if (!isTRUE(all.equal(minimum, 0)) || !isTRUE(all.equal(maximum, 262144)) || !isTRUE(all.equal(gain, 1))) {
      abort("Unexpected linear transform for ", context, "/", parameter_name)
    }
  }
  invisible(TRUE)
}

validate_identity_spill <- function(flow_frame, context) {
  keys <- keyword(flow_frame)
  spill <- keys[["$SPILL"]] %||% keys[["SPILL"]] %||% keys[["$SPILLOVER"]]
  if (is.null(spill) || !is.matrix(spill) || nrow(spill) != ncol(spill)) {
    abort("Missing or malformed FCS spillover matrix in ", context)
  }
  identity <- diag(nrow(spill))
  if (!isTRUE(all.equal(unname(spill), identity, tolerance = 1e-12))) {
    abort("Non-identity compensation is unsupported for the reviewed endpoint-flow module: ", context)
  }
  invisible(TRUE)
}

apply_polygon_gate <- function(flow_frame, gate, context) {
  observed <- colnames(exprs(flow_frame))
  coordinates <- gate$coordinates
  colnames(coordinates) <- vapply(gate$parameters, resolve_parameter,
    observed_names = observed, context = context, character(1))
  as.logical(flowCore::filter(flow_frame, polygonGate(.gate = coordinates))@subSet)
}

apply_rectangle_gate <- function(expression_matrix, gate, context) {
  mask <- rep(TRUE, nrow(expression_matrix))
  for (parameter in gate$parameters) {
    observed <- resolve_parameter(parameter, colnames(expression_matrix), context)
    bounds <- gate$bounds[[parameter]]
    mask <- mask & expression_matrix[, observed] >= bounds[["min"]] & expression_matrix[, observed] <= bounds[["max"]]
  }
  mask
}

agreement_tolerance <- function(workspace_count, workspace_parent_count) {
  # A small, fixed numerical tolerance; every discrepancy remains visible.
  max(5L, min(ceiling(0.001 * workspace_parent_count), ceiling(0.0025 * max(workspace_count, 1L))))
}

agreement_row <- function(mouse, role, population_name, gate_type,
    workspace_count, reconstructed_count, workspace_parent_count,
    reconstructed_parent_count, parent_name) {
  tolerance <- agreement_tolerance(workspace_count, workspace_parent_count)
  workspace_pct <- 100 * workspace_count / workspace_parent_count
  reconstructed_pct <- 100 * reconstructed_count / reconstructed_parent_count
  delta <- reconstructed_count - workspace_count
  data.frame(
    mouse_id = mouse$mouse_id,
    injected_origin = mouse$injected_origin,
    dose_mg_kg = mouse$dose_mg_kg,
    population_role = role,
    population_name = population_name,
    parent_population = parent_name,
    gate_type = gate_type,
    workspace_count = workspace_count,
    reconstructed_count = reconstructed_count,
    delta_events = delta,
    absolute_delta_events = abs(delta),
    workspace_parent_count = workspace_parent_count,
    reconstructed_parent_count = reconstructed_parent_count,
    workspace_pct_parent = workspace_pct,
    reconstructed_pct_parent = reconstructed_pct,
    delta_percentage_points = reconstructed_pct - workspace_pct,
    tolerance_events = tolerance,
    agreement_status = if (abs(delta) <= tolerance) "close_match" else "review_required",
    low_human_cells_flag = mouse$low_human_cells_flag,
    stringsAsFactors = FALSE
  )
}

deterministic_indices <- function(n, maximum) {
  if (n <= maximum) return(seq_len(n))
  unique(as.integer(round(seq.int(1, n, length.out = maximum))))
}

select_representative <- function(frozen, origin, min_human_cells) {
  eligible <- frozen[frozen$injected_origin == origin & frozen$human_cells_count >= min_human_cells, , drop = FALSE]
  if (!nrow(eligible)) abort("No QC-pass representative candidates for origin ", origin)
  target <- median(eligible$human_cells_count)
  eligible$distance_from_origin_median_human_cells <- abs(eligible$human_cells_count - target)
  eligible <- eligible[order(eligible$distance_from_origin_median_human_cells, eligible$mouse_id), , drop = FALSE]
  eligible$selected <- FALSE
  eligible$selected[[1L]] <- TRUE
  data.frame(
    selection_policy = "QC-pass requested-origin sample closest to the frozen HumanCells median; ties by mouse_id",
    representative_origin = origin,
    origin_median_frozen_human_cells = target,
    mouse_id = eligible$mouse_id,
    frozen_human_cells_count = eligible$human_cells_count,
    distance_from_origin_median_human_cells = eligible$distance_from_origin_median_human_cells,
    selected = eligible$selected,
    stringsAsFactors = FALSE
  )
}

format_dose <- function(x) {
  ifelse(as.numeric(x) == 0, "vehicle", paste0(format(as.numeric(x), trim = TRUE), " mg/kg"))
}

manuscript_theme <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text = element_text(color = "#333333", size = base_size - 1),
      axis.title = element_text(color = "#222222", size = base_size),
      axis.line = element_line(color = "#333333", linewidth = 0.3),
      axis.ticks = element_line(color = "#333333", linewidth = 0.3),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.3, color = "#333333", hjust = 0),
      strip.background = element_rect(fill = "#F2F2F2", color = "#D0D0D0", linewidth = 0.25),
      strip.text = element_text(size = base_size - 1, face = "bold", color = "#222222"),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size - 0.5),
      plot.margin = margin(4, 4, 4, 4)
    )
}

comma_labels <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

close_polygon <- function(coordinates) {
  rbind(coordinates, coordinates[1L, , drop = FALSE])
}

build_figure <- function(histograms, count_agreement, sample_summary,
    representative, representative_events, representative_gates,
    axis_min, axis_max, bin_width, max_display_points) {
  origin_colors <- c("2N" = "#0072B2", "4N" = "#D55E00")
  gate_colors <- c("workspace_named_2N" = "#0072B2", "workspace_named_4N" = "#D55E00")
  representative_mouse <- representative$mouse_id[representative$selected][[1L]]
  rep_summary <- sample_summary[sample_summary$mouse_id == representative_mouse, , drop = FALSE]
  rep_count <- count_agreement[count_agreement$mouse_id == representative_mouse, , drop = FALSE]
  hce_gate <- representative_gates[["human_cell_enrichment"]]
  human_gate <- representative_gates[["human_cells"]]
  two_gate <- representative_gates[["workspace_named_2N"]]
  four_gate <- representative_gates[["workspace_named_4N"]]
  peak_gate <- representative_gates[["reviewed_peak"]]

  hce_polygon <- as.data.frame(close_polygon(hce_gate$coordinates))
  names(hce_polygon) <- c("x", "y")
  human_polygon <- as.data.frame(close_polygon(human_gate$coordinates))
  names(human_polygon) <- c("x", "y")

  a1_counts <- rep_count[rep_count$population_role == "human_cell_enrichment", ]
  a2_counts <- rep_count[rep_count$population_role == "human_cells", ]
  p_a1 <- ggplot(representative_events$all, aes(FSC_A, SSC_A)) +
    geom_point(aes(color = inside_hce), size = 0.18, alpha = 0.42, stroke = 0) +
    geom_path(data = hce_polygon, aes(x, y), inherit.aes = FALSE,
      color = "#009E73", linewidth = 0.65, linejoin = "round") +
    annotate("label", x = 8000, y = 252000, hjust = 0, vjust = 1,
      label = paste0("FlowJo ", comma_labels(a1_counts$workspace_count),
        " | replay ", comma_labels(a1_counts$reconstructed_count)),
      size = 2.25, linewidth = 0.15, fill = "white") +
    scale_color_manual(values = c("FALSE" = "#BDBDBD", "TRUE" = "#009E73"), guide = "none") +
    scale_x_continuous(breaks = c(0, 100000, 200000), labels = comma_labels, expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 100000, 200000), labels = comma_labels, expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 262144), ylim = c(0, 262144), expand = FALSE) +
    labs(title = "1. Human Cell Enrichment", x = "FSC-A (raw a.u.)", y = "SSC-A (raw a.u.)") +
    manuscript_theme(7.6)

  p_a2 <- ggplot(representative_events$enrichment, aes(DNA, FSC_A)) +
    geom_point(aes(color = inside_human_cells), size = 0.18, alpha = 0.42, stroke = 0) +
    geom_path(data = human_polygon, aes(x, y), inherit.aes = FALSE,
      color = "#CC79A7", linewidth = 0.65, linejoin = "round") +
    annotate("label", x = axis_min + 2000, y = 252000, hjust = 0, vjust = 1,
      label = paste0("FlowJo ", comma_labels(a2_counts$workspace_count),
        " | replay ", comma_labels(a2_counts$reconstructed_count)),
      size = 2.25, linewidth = 0.15, fill = "white") +
    scale_color_manual(values = c("FALSE" = "#BDBDBD", "TRUE" = "#CC79A7"), guide = "none") +
    scale_x_continuous(breaks = c(25000, 75000, 125000), labels = comma_labels, expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 100000, 200000), labels = comma_labels, expand = c(0, 0)) +
    coord_cartesian(xlim = c(axis_min, axis_max), ylim = c(0, 262144), expand = FALSE) +
    labs(title = "2. HumanCells", x = "450/50 Violet B-A (raw a.u.)", y = "FSC-A (raw a.u.)") +
    manuscript_theme(7.6)

  rep_hist <- histograms[histograms$mouse_id == representative_mouse, , drop = FALSE]
  gate_band <- function(gate, fill, label) data.frame(
    xmin = gate$bounds[[1L]][["min"]], xmax = gate$bounds[[1L]][["max"]],
    fill = fill, label = label, stringsAsFactors = FALSE)
  generic_bands <- rbind(
    gate_band(two_gate, gate_colors[["workspace_named_2N"]], "2N"),
    gate_band(four_gate, gate_colors[["workspace_named_4N"]], "4N")
  )
  peak_band <- gate_band(peak_gate, "#009E73", paste0("reviewed peak: ", peak_gate$population_name))
  p_a3 <- ggplot(rep_hist, aes(bin_midpoint, percent_human_cells)) +
    geom_rect(data = generic_bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
      inherit.aes = FALSE, alpha = 0.12, color = NA) +
    geom_col(width = bin_width, fill = "#777777", color = "white", linewidth = 0.08) +
    geom_vline(xintercept = unlist(two_gate$bounds[[1L]]),
      color = gate_colors[["workspace_named_2N"]], linewidth = 0.45) +
    geom_vline(xintercept = unlist(four_gate$bounds[[1L]]),
      color = gate_colors[["workspace_named_4N"]], linewidth = 0.45) +
    geom_rect(data = peak_band, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = NA, color = "#009E73", linewidth = 0.6, linetype = 2) +
    scale_fill_identity() +
    annotate("text", x = c(mean(c(generic_bands$xmin[[1L]], generic_bands$xmax[[1L]])),
        mean(c(generic_bands$xmin[[2L]], generic_bands$xmax[[2L]]))),
      y = Inf, label = c("2N", "4N"), hjust = 0.5, vjust = 1.2,
      color = generic_bands$fill, size = 2.3, fontface = "bold") +
    annotate("text", x = peak_band$xmin + 500, y = Inf,
      label = paste0("peak ", peak_gate$population_name), hjust = 0, vjust = 2.7,
      color = peak_band$fill, size = 2.15) +
    scale_x_continuous(limits = c(axis_min, axis_max), breaks = c(25000, 75000, 125000), labels = comma_labels, expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = "3. Workspace DNA gates", x = "450/50 Violet B-A (raw a.u.)",
      y = sprintf("%% HumanCells / %s a.u.", comma_labels(bin_width))) +
    manuscript_theme(7.6)

  panel_a <- p_a1 + p_a2 + p_a3 + plot_layout(widths = c(1, 1, 1.05)) +
    plot_annotation(
      title = paste0("Representative raw-event hierarchy: ", representative_mouse),
      subtitle = paste0(rep_summary$injected_origin, " origin; ", format_dose(rep_summary$dose_mg_kg),
        "; workspace-defined per-sample coordinates\n",
        "Scatter display: up to ", comma_labels(max_display_points),
        " evenly spaced events; counts use every event"),
      theme = theme(
        plot.title = element_text(family = "sans", face = "bold", size = 9, hjust = 0),
        plot.subtitle = element_text(family = "sans", size = 7.5, hjust = 0)
      )
    )

  histograms$facet_label <- factor(histograms$facet_label,
    levels = unique(sample_summary$facet_label))
  panel_b <- ggplot(histograms, aes(bin_midpoint, percent_human_cells, group = mouse_id)) +
    geom_col(aes(fill = injected_origin), width = bin_width, color = "white", linewidth = 0.04) +
    facet_wrap(~facet_label, ncol = 4, scales = "fixed") +
    scale_fill_manual(values = origin_colors, guide = "none") +
    scale_x_continuous(limits = c(axis_min, axis_max), breaks = c(25000, 75000, 125000),
      labels = function(x) format(x / 1000, trim = TRUE), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Within-mouse-normalized DNA-content distributions",
      subtitle = paste0(
        "Replayed HumanCells use identical 1,000-a.u. bins and one raw fluorescence window;",
        " only y is normalized within mouse.\n",
        "Each facet reports its reviewed FlowJo peak annotation."
      ),
      x = "450/50 Violet B-A fluorescence (raw a.u., ×1,000)",
      y = "% HumanCells per bin") +
    manuscript_theme(7.2) +
    theme(panel.spacing = unit(2.5, "pt"), strip.text = element_text(size = 6.5),
      axis.text = element_text(size = 6.2), axis.title = element_text(size = 7.2),
      plot.title = element_text(size = 8.5), plot.subtitle = element_text(size = 6.8))

  wrap_elements(panel_a) / wrap_elements(panel_b) +
    plot_layout(heights = c(2.15, 3.25)) +
    plot_annotation(
      tag_levels = "A",
      caption = paste0(
        "† 2N-A1-0 was retained but flagged: 174 frozen FlowJo HumanCells (173 by replay).\n",
        "Human Cell Enrichment and HumanCells are workspace population names; no explicit singlet or viability gate is present.\n",
        "The eight 4N-origin reviewed peak annotations span 1.88N–2.20N and are printed in their Panel B facets.\n",
        "Peak-gate names are analyst-supplied FlowJo annotations, not newly calibrated or NUMBAT-equivalent ploidy estimates."
      ),
      theme = theme(
        plot.tag = element_text(family = "sans", face = "bold", size = 11, color = "#111111"),
        plot.margin = margin(t = 10, r = 4, b = 4, l = 4),
        plot.caption = element_text(family = "sans", size = 7, hjust = 0, color = "#333333",
          margin = margin(t = 3))
      )
    )
}

save_figure_atomic <- function(plot, pdf_path, png_path) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  pdf_temp <- file.path(dirname(pdf_path), paste0(".", tools::file_path_sans_ext(basename(pdf_path)), ".tmp.pdf"))
  png_temp <- file.path(dirname(png_path), paste0(".", tools::file_path_sans_ext(basename(png_path)), ".tmp.png"))
  on.exit(unlink(c(pdf_temp, png_temp)), add = TRUE)
  ggsave(pdf_temp, plot = plot, width = 7.1, height = 6.6, units = "in",
    device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
  ggsave(png_temp, plot = plot, width = 7.1, height = 6.6, units = "in",
    dpi = 300, device = ragg::agg_png, background = "white", limitsize = FALSE)
  if (!file.rename(pdf_temp, pdf_path) || !file.rename(png_temp, png_path)) abort("Cannot atomically publish figure outputs")
}

reconstruct_endpoint_flow <- function(config) {
  required_crosswalk <- c("mouse_id", "injected_origin", "dose_mg_kg", "fcs_file", "wsp_sample_name", "peak_gate_name")
  required_frozen <- c("mouse_id", "injected_origin", "dose_mg_kg", "fcs_total_count",
    "human_cell_enrichment_count", "human_cells_count", "human_cells_qc_pass",
    "fcs_file", "wsp_sample_name", "peak_gate_name", "peak_count", "two_n_count", "four_n_count")
  crosswalk <- read_tsv(config$crosswalk, "crosswalk")
  frozen <- read_tsv(config$frozen_table, "frozen per-sample table")
  assert_columns(crosswalk, required_crosswalk, "Crosswalk")
  assert_columns(frozen, required_frozen, "Frozen per-sample table")
  if (nrow(crosswalk) != config$expected_samples || nrow(frozen) != config$expected_samples) {
    abort("Crosswalk and frozen table must each contain exactly ", config$expected_samples, " rows")
  }
  if (anyDuplicated(crosswalk$mouse_id) || anyDuplicated(crosswalk$wsp_sample_name) || anyDuplicated(frozen$mouse_id)) {
    abort("Mouse and workspace sample identifiers must be unique")
  }
  if (!setequal(crosswalk$mouse_id, frozen$mouse_id)) abort("Crosswalk/frozen mouse universes differ")
  frozen_match <- frozen[match(crosswalk$mouse_id, frozen$mouse_id), , drop = FALSE]
  if (any(crosswalk$injected_origin != frozen_match$injected_origin) ||
      any(as.numeric(crosswalk$dose_mg_kg) != as.numeric(frozen_match$dose_mg_kg)) ||
      any(crosswalk$fcs_file != frozen_match$fcs_file) ||
      any(crosswalk$wsp_sample_name != frozen_match$wsp_sample_name) ||
      any(crosswalk$peak_gate_name != frozen_match$peak_gate_name)) {
    abort("Crosswalk and frozen table disagree on origin, dose, FCS file, WSP sample, or reviewed peak gate")
  }
  crosswalk <- merge(crosswalk, frozen[, c("mouse_id", "human_cells_count", "human_cells_qc_pass")],
    by = "mouse_id", all.x = TRUE, sort = FALSE)
  crosswalk <- crosswalk[match(frozen$mouse_id, crosswalk$mouse_id), , drop = FALSE]
  crosswalk$dose_mg_kg <- as.numeric(crosswalk$dose_mg_kg)
  frozen$dose_mg_kg <- as.numeric(frozen$dose_mg_kg)
  for (field in c("fcs_total_count", "human_cell_enrichment_count", "human_cells_count", "peak_count", "two_n_count", "four_n_count")) {
    frozen[[field]] <- as.integer(frozen[[field]])
  }
  selection <- select_representative(frozen, config$representative_origin, config$min_human_cells)
  representative_mouse <- selection$mouse_id[selection$selected][[1L]]

  workspace <- read_xml(config$workspace)
  sample_nodes <- xml_find_all(workspace,
    "/*[local-name()='Workspace']/*[local-name()='SampleList']/*[local-name()='Sample']/*[local-name()='SampleNode']")
  if (anyDuplicated(xml_attr(sample_nodes, "name"))) abort("Workspace SampleNode names are not unique")
  breaks <- seq(config$dna_axis_min, config$dna_axis_max, by = config$dna_bin_width)

  agreement_parts <- list()
  geometry_parts <- list()
  histogram_parts <- list()
  summary_parts <- list()
  representative_events <- NULL
  representative_gates <- NULL

  for (row_index in seq_len(nrow(crosswalk))) {
    mouse <- as.list(crosswalk[row_index, , drop = FALSE])
    mouse$mouse_id <- as.character(mouse$mouse_id)
    mouse$injected_origin <- as.character(mouse$injected_origin)
    mouse$dose_mg_kg <- as.numeric(mouse$dose_mg_kg)
    frozen_row <- frozen[frozen$mouse_id == mouse$mouse_id, , drop = FALSE]
    sample_node <- sample_node_for_name(sample_nodes, mouse$wsp_sample_name)
    sample_container <- sample_container_for_node(sample_node)
    workspace_total <- suppressWarnings(as.integer(xml_attr(sample_node, "count")))
    hce_population <- direct_population(sample_node, "Human Cell Enrichment", mouse$mouse_id)
    human_population <- direct_population(hce_population, "HumanCells", paste0(mouse$mouse_id, "/Human Cell Enrichment"))
    two_population <- direct_population(human_population, "2N", paste0(mouse$mouse_id, "/HumanCells"))
    four_population <- direct_population(human_population, "4N", paste0(mouse$mouse_id, "/HumanCells"))
    peak_population <- direct_population(human_population, mouse$peak_gate_name, paste0(mouse$mouse_id, "/HumanCells"))
    gates <- list(
      human_cell_enrichment = parse_gate(hce_population, "Human Cell Enrichment", "human_cell_enrichment", mouse$mouse_id),
      human_cells = parse_gate(human_population, "HumanCells", "human_cells", mouse$mouse_id),
      workspace_named_2N = parse_gate(two_population, "2N", "workspace_named_2N", mouse$mouse_id),
      workspace_named_4N = parse_gate(four_population, "4N", "workspace_named_4N", mouse$mouse_id),
      reviewed_peak = parse_gate(peak_population, mouse$peak_gate_name, "reviewed_peak", mouse$mouse_id)
    )
    if (!identical(gates$human_cell_enrichment$parameters, c("FSC-A", "SSC-A")) ||
        !identical(gates$human_cells$parameters, c("450_50 Violet B-A", "FSC-A")) ||
        !identical(gates$workspace_named_2N$parameters, "450_50 Violet B-A") ||
        !identical(gates$workspace_named_4N$parameters, "450_50 Violet B-A") ||
        !identical(gates$reviewed_peak$parameters, "450_50 Violet B-A")) {
      abort("Unexpected gate dimensions or ordering for ", mouse$mouse_id)
    }
    relevant_parameters <- unlist(lapply(gates, `[[`, "parameters"), use.names = FALSE)
    validate_linear_transforms(sample_container, relevant_parameters, mouse$mouse_id)

    fcs_path <- file.path(config$fcs_dir, mouse$fcs_file)
    if (!file.exists(fcs_path)) abort("Missing FCS for ", mouse$mouse_id, ": ", fcs_path)
    frame <- read.FCS(fcs_path, transformation = FALSE, alter.names = FALSE,
      truncate_max_range = FALSE, emptyValue = FALSE)
    validate_identity_spill(frame, mouse$mouse_id)
    if (nrow(exprs(frame)) != workspace_total || workspace_total != frozen_row$fcs_total_count) {
      abort("FCS/FlowJo/frozen total-event mismatch for ", mouse$mouse_id)
    }
    hce_mask <- apply_polygon_gate(frame, gates$human_cell_enrichment, paste0(mouse$mouse_id, "/HCE"))
    hce_frame <- frame[hce_mask, ]
    human_within_hce <- apply_polygon_gate(hce_frame, gates$human_cells, paste0(mouse$mouse_id, "/HumanCells"))
    human_indices <- which(hce_mask)[which(human_within_hce)]
    human_expression <- exprs(frame)[human_indices, , drop = FALSE]
    dna_name <- resolve_parameter("450_50 Violet B-A", colnames(human_expression), mouse$mouse_id)
    dna <- human_expression[, dna_name]
    two_mask <- apply_rectangle_gate(human_expression, gates$workspace_named_2N, paste0(mouse$mouse_id, "/2N"))
    four_mask <- apply_rectangle_gate(human_expression, gates$workspace_named_4N, paste0(mouse$mouse_id, "/4N"))
    peak_mask <- apply_rectangle_gate(human_expression, gates$reviewed_peak, paste0(mouse$mouse_id, "/peak"))
    reconstructed <- c(
      human_cell_enrichment = sum(hce_mask),
      human_cells = length(human_indices),
      workspace_named_2N = sum(two_mask),
      workspace_named_4N = sum(four_mask),
      reviewed_peak = sum(peak_mask)
    )
    workspace_counts <- c(
      human_cell_enrichment = population_count(hce_population, paste0(mouse$mouse_id, "/HCE")),
      human_cells = population_count(human_population, paste0(mouse$mouse_id, "/HumanCells")),
      workspace_named_2N = population_count(two_population, paste0(mouse$mouse_id, "/2N")),
      workspace_named_4N = population_count(four_population, paste0(mouse$mouse_id, "/4N")),
      reviewed_peak = population_count(peak_population, paste0(mouse$mouse_id, "/", mouse$peak_gate_name))
    )
    frozen_expected <- c(
      human_cell_enrichment = frozen_row$human_cell_enrichment_count,
      human_cells = frozen_row$human_cells_count,
      workspace_named_2N = frozen_row$two_n_count,
      workspace_named_4N = frozen_row$four_n_count,
      reviewed_peak = frozen_row$peak_count
    )
    if (!identical(as.integer(workspace_counts), as.integer(frozen_expected))) {
      abort("Workspace gate counts disagree with frozen extractor for ", mouse$mouse_id)
    }
    mouse$low_human_cells_flag <- workspace_counts[["human_cells"]] < config$min_human_cells
    agreement_parts[[mouse$mouse_id]] <- rbind(
      agreement_row(mouse, "human_cell_enrichment", "Human Cell Enrichment", "PolygonGate",
        workspace_counts[["human_cell_enrichment"]], reconstructed[["human_cell_enrichment"]],
        workspace_total, workspace_total, "All events"),
      agreement_row(mouse, "human_cells", "HumanCells", "PolygonGate",
        workspace_counts[["human_cells"]], reconstructed[["human_cells"]],
        workspace_counts[["human_cell_enrichment"]], reconstructed[["human_cell_enrichment"]], "Human Cell Enrichment"),
      agreement_row(mouse, "workspace_named_2N", "2N", "RectangleGate",
        workspace_counts[["workspace_named_2N"]], reconstructed[["workspace_named_2N"]],
        workspace_counts[["human_cells"]], reconstructed[["human_cells"]], "HumanCells"),
      agreement_row(mouse, "workspace_named_4N", "4N", "RectangleGate",
        workspace_counts[["workspace_named_4N"]], reconstructed[["workspace_named_4N"]],
        workspace_counts[["human_cells"]], reconstructed[["human_cells"]], "HumanCells"),
      agreement_row(mouse, "reviewed_peak", mouse$peak_gate_name, "RectangleGate",
        workspace_counts[["reviewed_peak"]], reconstructed[["reviewed_peak"]],
        workspace_counts[["human_cells"]], reconstructed[["human_cells"]], "HumanCells")
    )
    geometry_parts[[mouse$mouse_id]] <- do.call(rbind, lapply(gates, gate_geometry_rows))
    underflow <- sum(dna < config$dna_axis_min)
    overflow <- sum(dna > config$dna_axis_max)
    if (underflow || overflow) abort("Common DNA display window excludes reconstructed HumanCells events for ", mouse$mouse_id)
    bins <- cut(dna, breaks = breaks, include.lowest = TRUE, right = TRUE)
    counts <- as.integer(table(factor(bins, levels = levels(bins))))
    if (sum(counts) != length(dna)) abort("Histogram binning lost events for ", mouse$mouse_id)
    dose_label <- format_dose(mouse$dose_mg_kg)
    low_mark <- if (mouse$low_human_cells_flag) "† " else ""
    facet_label <- paste0(
      low_mark, mouse$mouse_id,
      "\n", mouse$injected_origin, " · ", dose_label,
      "\nn=", format(reconstructed[["human_cells"]], big.mark = ","),
      " | peak=", mouse$peak_gate_name
    )
    histogram_parts[[mouse$mouse_id]] <- data.frame(
      mouse_id = mouse$mouse_id,
      injected_origin = mouse$injected_origin,
      dose_mg_kg = mouse$dose_mg_kg,
      bin_left = head(breaks, -1L),
      bin_right = tail(breaks, -1L),
      bin_midpoint = (head(breaks, -1L) + tail(breaks, -1L)) / 2,
      event_count = counts,
      percent_human_cells = 100 * counts / length(dna),
      probability_mass = counts / length(dna),
      facet_label = facet_label,
      stringsAsFactors = FALSE
    )
    summary_parts[[mouse$mouse_id]] <- data.frame(
      mouse_id = mouse$mouse_id,
      injected_origin = mouse$injected_origin,
      dose_mg_kg = mouse$dose_mg_kg,
      dose_label = dose_label,
      fcs_total_events = nrow(exprs(frame)),
      frozen_human_cell_enrichment_count = workspace_counts[["human_cell_enrichment"]],
      reconstructed_human_cell_enrichment_count = reconstructed[["human_cell_enrichment"]],
      frozen_human_cells_count = workspace_counts[["human_cells"]],
      reconstructed_human_cells_count = reconstructed[["human_cells"]],
      low_human_cells_flag = mouse$low_human_cells_flag,
      reviewed_peak_annotation = mouse$peak_gate_name,
      minimum_reconstructed_human_cells_dna = min(dna),
      maximum_reconstructed_human_cells_dna = max(dna),
      display_underflow_events = underflow,
      display_overflow_events = overflow,
      facet_label = facet_label,
      comparison_label = paste0(low_mark, mouse$mouse_id, "  [", mouse$injected_origin, ", ", dose_label, "]"),
      stringsAsFactors = FALSE
    )

    if (mouse$mouse_id == representative_mouse) {
      all_expression <- exprs(frame)
      all_indices <- deterministic_indices(nrow(all_expression), config$max_representative_points)
      hce_expression <- exprs(hce_frame)
      hce_indices <- deterministic_indices(nrow(hce_expression), config$max_representative_points)
      fsc_name <- resolve_parameter("FSC-A", colnames(all_expression), mouse$mouse_id)
      ssc_name <- resolve_parameter("SSC-A", colnames(all_expression), mouse$mouse_id)
      dna_all_name <- resolve_parameter("450_50 Violet B-A", colnames(all_expression), mouse$mouse_id)
      representative_events <- list(
        all = data.frame(FSC_A = all_expression[all_indices, fsc_name],
          SSC_A = all_expression[all_indices, ssc_name], inside_hce = hce_mask[all_indices]),
        enrichment = data.frame(DNA = hce_expression[hce_indices, dna_all_name],
          FSC_A = hce_expression[hce_indices, fsc_name], inside_human_cells = human_within_hce[hce_indices])
      )
      representative_gates <- gates
    }
  }

  count_agreement <- do.call(rbind, agreement_parts)
  rownames(count_agreement) <- NULL
  gate_geometry <- do.call(rbind, geometry_parts)
  rownames(gate_geometry) <- NULL
  histograms <- do.call(rbind, histogram_parts)
  rownames(histograms) <- NULL
  sample_summary <- do.call(rbind, summary_parts)
  rownames(sample_summary) <- NULL
  sample_summary$origin_order <- match(sample_summary$injected_origin, c("2N", "4N"))
  sample_summary <- sample_summary[order(sample_summary$origin_order, sample_summary$dose_mg_kg, sample_summary$mouse_id), , drop = FALSE]
  sample_summary$origin_order <- NULL
  # Reorder explicitly without relying on list or row-name order.
  histograms$order_key <- match(histograms$mouse_id, sample_summary$mouse_id)
  histograms <- histograms[order(histograms$order_key, histograms$bin_left), , drop = FALSE]
  histograms$order_key <- NULL
  if (!all(abs(rowsum(histograms$probability_mass, histograms$mouse_id) - 1) < 1e-12)) {
    abort("Within-mouse histogram probabilities do not sum to one")
  }
  if (is.null(representative_events) || is.null(representative_gates)) abort("Representative sample was not reconstructed")

  if (any(count_agreement$agreement_status != "close_match")) {
    abort("One or more raw-event gate replays require review under the fixed agreement tolerance")
  }

  output_dir <- config$output_dir
  figure_dir <- file.path(output_dir, "figures")
  metadata_dir <- file.path(output_dir, "metadata")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv_atomic(count_agreement, file.path(output_dir, "endpoint_flow_count_agreement.tsv"))
  write_tsv_atomic(gate_geometry, file.path(output_dir, "endpoint_flow_gate_geometry.tsv"))
  histograms_output <- histograms
  histograms_output$facet_label <- gsub("\n", " | ", histograms_output$facet_label, fixed = TRUE)
  sample_summary_output <- sample_summary
  sample_summary_output$facet_label <- gsub("\n", " | ", sample_summary_output$facet_label, fixed = TRUE)
  write_tsv_atomic(histograms_output, file.path(output_dir, "endpoint_flow_dna_histograms.tsv"))
  write_tsv_atomic(sample_summary_output, file.path(output_dir, "endpoint_flow_reconstruction_per_sample.tsv"))
  selected_index <- which(selection$selected)
  selection$display_sampling_policy <- "evenly spaced event indices for plotting only; all events used for gate counts"
  selection$maximum_display_events_per_scatter <- config$max_representative_points
  selection$displayed_all_events <- NA_integer_
  selection$displayed_human_cell_enrichment_events <- NA_integer_
  selection$displayed_all_events[selected_index] <- nrow(representative_events$all)
  selection$displayed_human_cell_enrichment_events[selected_index] <- nrow(representative_events$enrichment)
  write_tsv_atomic(selection, file.path(output_dir, "endpoint_flow_representative_selection.tsv"))

  named_gate_summary <- count_agreement[count_agreement$population_role %in% c("workspace_named_2N", "workspace_named_4N"),
    c("mouse_id", "injected_origin", "dose_mg_kg", "population_name", "parent_population",
      "workspace_count", "reconstructed_count", "delta_events", "workspace_parent_count",
      "reconstructed_parent_count", "workspace_pct_parent", "reconstructed_pct_parent",
      "delta_percentage_points", "tolerance_events", "agreement_status", "low_human_cells_flag"), drop = FALSE]
  write_tsv_atomic(named_gate_summary, file.path(output_dir, "endpoint_flow_named_2n_4n_summary.tsv"))

  max_delta <- aggregate(absolute_delta_events ~ population_role, count_agreement, max)
  config_rows <- data.frame(
    key = c(
      "analysis_universe", "event_intensity_policy", "relevant_wsp_transform",
      "spillover_policy", "dna_parameter", "dna_axis_units", "dna_axis_min",
      "dna_axis_max", "dna_bin_width", "histogram_normalization",
      "representative_policy", "representative_mouse", "representative_scatter_display",
      "low_count_policy",
      "agreement_tolerance", "all_gate_replays_within_tolerance",
      "R_version", "flowCore_version", "xml2_version", "ggplot2_version",
      "patchwork_version", "ragg_version", "png_device", "pdf_device", "cairo_version"
    ),
    value = c(
      paste0(config$expected_samples, " explicitly crosswalked endpoint tumors"),
      "recorded FCS values; no centering, peak alignment, P7BS offset, or x rescaling",
      "linear 0-262144 with gain 1 for FSC-A, SSC-A, and 450_50 Violet B-A",
      "identity matrix verified independently in every selected FCS file",
      "450/50 Violet B-A",
      "unscaled arbitrary detector units",
      format(config$dna_axis_min, scientific = FALSE, trim = TRUE),
      format(config$dna_axis_max, scientific = FALSE, trim = TRUE),
      format(config$dna_bin_width, scientific = FALSE, trim = TRUE),
      "probability mass sums to 1 independently within each reconstructed HumanCells sample",
      selection$selection_policy[[1L]],
      representative_mouse,
      paste0("up to ", config$max_representative_points,
        " evenly spaced event indices per scatter; all events used for counts"),
      paste0("retain every mouse; flag frozen HumanCells count <", config$min_human_cells),
      "max(5 events, min(0.1% of frozen parent, 0.25% of frozen population)); numerical tolerance without attributing a specific software cause",
      if (all(count_agreement$agreement_status == "close_match")) "TRUE" else "FALSE",
      R.version.string,
      as.character(packageVersion("flowCore")), as.character(packageVersion("xml2")),
      as.character(packageVersion("ggplot2")), as.character(packageVersion("patchwork")),
      as.character(packageVersion("ragg")), "ragg::agg_png at 300 dpi",
      "grDevices::cairo_pdf", as.character(grSoftVersion()[["cairo"]] %||% "unavailable")
    ),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(max_delta))) {
    config_rows <- rbind(config_rows, data.frame(
      key = paste0("max_absolute_delta_events_", max_delta$population_role[[i]]),
      value = as.character(max_delta$absolute_delta_events[[i]]), stringsAsFactors = FALSE))
  }
  write_tsv_atomic(config_rows, file.path(metadata_dir, "run_config.tsv"))

  plot <- build_figure(histograms, count_agreement, sample_summary, selection,
    representative_events, representative_gates, config$dna_axis_min,
    config$dna_axis_max, config$dna_bin_width, config$max_representative_points)
  pdf_path <- file.path(figure_dir, "panel_SuppFig9_endpoint_flow_cytometry.pdf")
  png_path <- file.path(figure_dir, "panel_SuppFig9_endpoint_flow_cytometry.png")
  save_figure_atomic(plot, pdf_path, png_path)
  cat("Endpoint-flow raw-event reconstruction complete\n")
  cat("Representative mouse: ", representative_mouse, "\n", sep = "")
  cat("Gate replay close matches: ", sum(count_agreement$agreement_status == "close_match"),
    "/", nrow(count_agreement), "\n", sep = "")
  cat("Maximum absolute count differences by population:\n")
  for (i in seq_len(nrow(max_delta))) cat("- ", max_delta$population_role[[i]], ": ", max_delta$absolute_delta_events[[i]], "\n", sep = "")
  cat("Figure PDF: ", pdf_path, "\n", sep = "")
  cat("Figure PNG: ", png_path, "\n", sep = "")
  invisible(list(count_agreement = count_agreement, histograms = histograms,
    sample_summary = sample_summary, selection = selection, figure = plot))
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  config <- parse_cli(args)
  if (isTRUE(attr(config, "help"))) return(invisible(NULL))
  reconstruct_endpoint_flow(config)
}

if (sys.nframe() == 0L) main()
