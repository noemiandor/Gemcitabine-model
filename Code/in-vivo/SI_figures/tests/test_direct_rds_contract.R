#!/usr/bin/env Rscript

test_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(test_file_arg)) {
  stop("Cannot resolve test path", call. = FALSE)
}
test_path <- normalizePath(
  sub("^--file=", "", test_file_arg[[1L]]),
  mustWork = TRUE
)
repo_root <- normalizePath(
  file.path(dirname(test_path), "..", "..", "..", ".."),
  mustWork = TRUE
)

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

builder_env <- new.env(parent = globalenv())
sys.source(
  file.path(
    repo_root,
    "Code", "in-vivo", "SI_figures",
    "build_raw_supplementary_tables.R"
  ),
  envir = builder_env
)

authoritative_fingerprint <- strrep("a", 64L)
stopifnot(
  identical(
    builder_env$resolve_work_dependency_fingerprint(
      NULL,
      authoritative_fingerprint
    ),
    authoritative_fingerprint
  ),
  identical(
    builder_env$resolve_work_dependency_fingerprint(
      authoritative_fingerprint,
      authoritative_fingerprint
    ),
    authoritative_fingerprint
  )
)
expect_error(
  builder_env$resolve_work_dependency_fingerprint(
    strrep("b", 64L),
    authoritative_fingerprint
  ),
  "does not match the authoritative fingerprint"
)

canonical <- data.frame(
  cell = c("cell-a", "cell-b", "cell-c"),
  initial_ploidy = c("2N", "4N", "4N"),
  context = c("Tumor", "Tumor", "CellLine"),
  dose = c("0mg/kg", "30mg/kg", "120mg/kg"),
  cluster = c("0", "2", "4c"),
  stringsAsFactors = FALSE
)
consistent_metrics <- data.frame(
  cell = c("cell-c", "cell-a", "cell-b"),
  Ploidy = c("4N", "2N", "4N"),
  TN = c("cell-culture", "tumor", "tumor"),
  Dose = c("120", "0", "30"),
  clusters = c("4c", "0", "2"),
  stringsAsFactors = FALSE
)

without_metrics <- builder_env$audit_optional_scvelo(
  canonical,
  metrics = NULL,
  canonical_cluster_col = "clusters"
)
with_consistent_metrics <- builder_env$audit_optional_scvelo(
  canonical,
  metrics = consistent_metrics,
  canonical_cluster_col = "clusters"
)
stopifnot(identical(without_metrics, canonical))
stopifnot(identical(with_consistent_metrics, canonical))

contradictory_metrics <- consistent_metrics
contradictory_metrics$Ploidy[contradictory_metrics$cell == "cell-a"] <- "4N"
expect_error(
  builder_env$audit_optional_scvelo(
    canonical,
    metrics = contradictory_metrics,
    canonical_cluster_col = "clusters"
  ),
  "Seurat/scVelo ploidy audit disagrees"
)

unresolved_metrics <- consistent_metrics
unresolved_metrics$Ploidy[unresolved_metrics$cell == "cell-a"] <- "unknown"
expect_error(
  builder_env$audit_optional_scvelo(
    canonical,
    metrics = unresolved_metrics,
    canonical_cluster_col = "clusters"
  ),
  "Optional scVelo metrics contain unresolved"
)

missing_tumor_dose <- consistent_metrics
missing_tumor_dose$Dose[missing_tumor_dose$cell == "cell-a"] <- NA_character_
expect_error(
  builder_env$audit_optional_scvelo(
    canonical,
    metrics = missing_tumor_dose,
    canonical_cluster_col = "clusters"
  ),
  "Optional scVelo metrics contain missing tumor dose"
)

configured <- builder_env$resolve_reviewed_ploidy_context(
  data.frame(
    Ploidy = c("2N", "4N", "4N"),
    TN = c("Tumor", "CellLine", "Tumor"),
    IDs = c("2N-tumor", "sample-cell-culture", "tumor"),
    sample_folder = c("A1_2N", "A5_cellline", "A6_tumor"),
    mouse = c("A1", "A5", "A6"),
    stringsAsFactors = FALSE
  ),
  ploidy_col = "Ploidy",
  context_col = "TN",
  sample_col = "mouse",
  id_col = "IDs",
  sample_folder_col = "sample_folder"
)
stopifnot(
  identical(configured$initial_ploidy, c("2N", "4N", "4N")),
  identical(configured$context, c("Tumor", "CellLine", "Tumor"))
)

# Reviewed configured fields remain authoritative when ID text has no usable
# inference; inference can only detect contradictions, never supply values.
configured_without_inference <- builder_env$resolve_reviewed_ploidy_context(
  data.frame(
    Ploidy = "2N",
    TN = "Tumor",
    IDs = NA_character_,
    mouse = "reviewed-sample",
    stringsAsFactors = FALSE
  ),
  ploidy_col = "Ploidy",
  context_col = "TN",
  sample_col = "mouse",
  id_col = "IDs"
)
stopifnot(
  identical(configured_without_inference$initial_ploidy, "2N"),
  identical(configured_without_inference$context, "Tumor")
)

expect_error(
  builder_env$resolve_reviewed_ploidy_context(
    data.frame(
      Ploidy = c("2N", "4N"),
      TN = c("Tumor", "Tumor"),
      IDs = c("2N-tumor", "4N-tumor"),
      sample_folder = c("A1_2N", "A1_4N"),
      mouse = c("A1", "A1"),
      stringsAsFactors = FALSE
    ),
    ploidy_col = "Ploidy",
    context_col = "TN",
    sample_col = "mouse",
    id_col = "IDs",
    sample_folder_col = "sample_folder"
  ),
  "Each mouse must map to one reviewed initial ploidy"
)

expect_error(
  builder_env$resolve_reviewed_ploidy_context(
    data.frame(
      Ploidy = "2N",
      TN = "Tumor",
      IDs = "4N-tumor",
      mouse = "A1",
      stringsAsFactors = FALSE
    ),
    ploidy_col = "Ploidy",
    context_col = "TN",
    sample_col = "mouse",
    id_col = "IDs"
  ),
  "Configured RDS ploidy and ID-derived ploidy disagrees"
)

expect_error(
  builder_env$resolve_reviewed_ploidy_context(
    data.frame(
      Ploidy = "2N",
      TN = "Tumor",
      IDs = "2N-cell-culture",
      mouse = "A1",
      stringsAsFactors = FALSE
    ),
    ploidy_col = "Ploidy",
    context_col = "TN",
    sample_col = "mouse",
    id_col = "IDs"
  ),
  "Configured RDS context and ID-derived context disagrees"
)

# A raw SI7 retry may reuse only a byte-verified, dependency-bound cluster DEG
# table.  This contract is deliberately testable without loading Seurat.
deg_cache_dir <- tempfile("si7-deg-resume-")
dir.create(deg_cache_dir, recursive = TRUE, showWarnings = FALSE)
deg_fingerprint <- strrep("a", 64L)
valid_de <- data.frame(
  p_val = c(1e-8, 2e-6),
  avg_log2FC = c(1.5, -0.8),
  pct.1 = c(0.8, 0.2),
  pct.2 = c(0.3, 0.5),
  p_val_adj = c(2e-7, 3e-5),
  gene = c("GRCh38_GENE1", "GRCm39_Gene2"),
  gene_symbol = c("GENE1", "Gene2"),
  cluster = c("0", "0"),
  stringsAsFactors = FALSE
)
builder_env$si7_write_cluster_deg_cache(
  valid_de,
  deg_cache_dir,
  "0",
  deg_fingerprint
)
reused_de <- builder_env$si7_read_reusable_cluster_deg(
  deg_cache_dir,
  "0",
  deg_fingerprint
)
reused_de$de$cluster <- as.character(reused_de$de$cluster)
stopifnot(
  !is.null(reused_de),
  identical(reused_de$cluster, "0"),
  isTRUE(all.equal(
    reused_de$de,
    valid_de,
    check.attributes = FALSE,
    tolerance = 0
  )),
  identical(reused_de$reused, TRUE)
)

deg_paths <- builder_env$si7_deg_cache_paths(deg_cache_dir, "0")
write("tampered", file = deg_paths$table, append = TRUE)
stopifnot(is.null(builder_env$si7_read_reusable_cluster_deg(
  deg_cache_dir,
  "0",
  deg_fingerprint
)))

builder_env$si7_write_cluster_deg_cache(
  valid_de,
  deg_cache_dir,
  "0",
  deg_fingerprint
)
stopifnot(is.null(builder_env$si7_read_reusable_cluster_deg(
  deg_cache_dir,
  "0",
  strrep("b", 64L)
)))

duplicate_de <- rbind(valid_de, valid_de[1L, , drop = FALSE])
builder_env$si7_write_cluster_deg_cache(
  duplicate_de,
  deg_cache_dir,
  "0",
  deg_fingerprint
)
stopifnot(is.null(builder_env$si7_read_reusable_cluster_deg(
  deg_cache_dir,
  "0",
  deg_fingerprint
)))

malformed_de <- valid_de[, setdiff(names(valid_de), "p_val_adj"), drop = FALSE]
builder_env$si7_write_cluster_deg_cache(
  malformed_de,
  deg_cache_dir,
  "0",
  deg_fingerprint
)
stopifnot(is.null(builder_env$si7_read_reusable_cluster_deg(
  deg_cache_dir,
  "0",
  deg_fingerprint
)))

# The scoped SI scientific contract changes for table-affecting helpers and
# value-to-canonical-filename mappings, but not for audit-only prose.
builder_path <- file.path(
  repo_root,
  "Code", "in-vivo", "SI_figures",
  "build_raw_supplementary_tables.R"
)
contract_for_script <- function(path) {
  environment <- new.env(parent = globalenv())
  sys.source(path, envir = environment)
  environment$si_raw_scientific_code_contract(path)
}
mutated_builder <- function(from, to) {
  lines <- readLines(builder_path, warn = FALSE)
  hits <- grepl(from, lines, fixed = TRUE)
  stopifnot(sum(hits) == 1L)
  lines[hits] <- sub(from, to, lines[hits], fixed = TRUE)
  path <- tempfile("si-builder-contract-", fileext = ".R")
  writeLines(lines, path)
  path
}
baseline_science_contract <- contract_for_script(builder_path)
audit_only_mutation <- mutated_builder(
  "Built exactly 11 raw SI plot-facing tables",
  "Audit: built exactly 11 raw SI plot-facing tables"
)
mapping_mutation <- mutated_builder(
  "write_csv(canonical_cells, file.path(table_dir, \"si_figures_cell_metadata.csv\"))",
  "write_csv(cluster_key, file.path(table_dir, \"si_figures_cell_metadata.csv\"))"
)
helper_mutation <- mutated_builder(
  "standardize_dose <- function(values, label = \"Dose\") {",
  "standardize_dose <- function(values, label = \"Reviewed dose\") {"
)
stopifnot(
  identical(
    baseline_science_contract,
    contract_for_script(audit_only_mutation)
  ),
  !identical(
    baseline_science_contract,
    contract_for_script(mapping_mutation)
  ),
  !identical(
    baseline_science_contract,
    contract_for_script(helper_mutation)
  )
)

# A non-deposited RDS is accepted only through the strict upstream validator.
nondeposited_rds <- tempfile("si-nondeposited-", fileext = ".rds")
saveRDS(list(test = TRUE), nondeposited_rds)
expect_error(
  builder_env$resolve_si_seurat_source_lineage(
    seurat_rds_path = nondeposited_rds,
    deposited_sha256 = strrep("0", 64L),
    upstream_dir = "",
    module_dir = "/module",
    environment_lock = "/lock",
    config = list(),
    all_ploidy = "/ploidy",
    sample_info = "/samples"
  ),
  "requires --seurat-upstream-dir"
)
validated_hash <- builder_env$file_sha256(nondeposited_rds)
builder_env$figure7_validate_seurat_upstream_artifact <- function(...) {
  list(
    rds_sha256 = validated_hash,
    dependencies = c(all_ploidy = strrep("1", 64L)),
    scientific_code_contracts = c(
      "scientific_code_contract:final" = strrep("2", 64L)
    ),
    upstream_common_contract = c(runtime = strrep("3", 64L))
  )
}
builder_env$figure7_seurat_source_dependency_values <- function(
  rds_sha256,
  ...
) {
  c(
    source_seurat_rds = rds_sha256,
    "seurat_transitive_dependency:all_ploidy" = strrep("1", 64L),
    "seurat_scientific_code_contract:final" = strrep("2", 64L),
    seurat_upstream_scientific_common_contract = strrep("3", 64L)
  )
}
manifested_lineage <- builder_env$resolve_si_seurat_source_lineage(
  seurat_rds_path = nondeposited_rds,
  deposited_sha256 = strrep("0", 64L),
  upstream_dir = "/manifested",
  module_dir = "/module",
  environment_lock = "/lock",
  config = list(),
  all_ploidy = "/ploidy",
  sample_info = "/samples"
)
stopifnot(
  identical(manifested_lineage$kind, "manifested_reconstruction"),
  identical(
    unname(manifested_lineage$dependencies[["source_seurat_rds"]]),
    validated_hash
  )
)
equal_byte_manifested_lineage <-
  builder_env$resolve_si_seurat_source_lineage(
    seurat_rds_path = nondeposited_rds,
    deposited_sha256 = validated_hash,
    upstream_dir = "/manifested",
    module_dir = "/module",
    environment_lock = "/lock",
    config = list(),
    all_ploidy = "/ploidy",
    sample_info = "/samples"
  )
stopifnot(identical(
  equal_byte_manifested_lineage$kind,
  "manifested_reconstruction"
))

wrapper_env <- new.env(parent = globalenv())
wrapper_path <- file.path(
  repo_root,
  "Code", "in-vivo", "SI_figures",
  "run_supplementary_figures.R"
)
sys.source(wrapper_path, envir = wrapper_env)

hashes <- list(
  science_contract = strrep("1", 64L),
  config_contract = strrep("2", 64L),
  runtime_contract = strrep("3", 64L),
  ploidy = strrep("4", 64L),
  rds = strrep("5", 64L),
  metrics = strrep("6", 64L),
  audit_builder = strrep("7", 64L),
  audit_common = strrep("8", 64L),
  audit_lock = strrep("9", 64L)
)
absent_dependencies <- wrapper_env$computational_dependencies(
  hashes$science_contract,
  hashes$config_contract,
  hashes$runtime_contract,
  hashes$ploidy,
  c(source_seurat_rds = hashes$rds)
)
audited_dependencies <- wrapper_env$computational_dependencies(
  hashes$science_contract,
  hashes$config_contract,
  hashes$runtime_contract,
  hashes$ploidy,
  c(source_seurat_rds = hashes$rds),
  audit_optional_scvelo = hashes$metrics
)
stopifnot(identical(absent_dependencies, audited_dependencies))
stopifnot(identical(
  wrapper_env$computational_fingerprint(absent_dependencies),
  wrapper_env$computational_fingerprint(audited_dependencies)
))
changed_science_dependencies <- absent_dependencies
changed_science_dependencies[["si_raw_scientific_code_contract"]] <-
  strrep("a", 64L)
stopifnot(!identical(
  wrapper_env$computational_fingerprint(absent_dependencies),
  wrapper_env$computational_fingerprint(changed_science_dependencies)
))
valid_manifest <- data.frame(
  role = c(names(absent_dependencies), "audit_optional_scvelo"),
  locator = c(
    paste0("contract:", names(absent_dependencies)),
    "external:metrics.csv"
  ),
  sha256 = c(unname(absent_dependencies), hashes$metrics),
  bytes = c(rep(0, length(absent_dependencies)), 10),
  stringsAsFactors = FALSE
)
duplicate_manifest <- rbind(valid_manifest, valid_manifest[1L, , drop = FALSE])
rogue_manifest <- rbind(
  valid_manifest,
  data.frame(
    role = "unscoped_extra",
    locator = "contract:unscoped_extra",
    sha256 = strrep("f", 64L),
    bytes = 0,
    stringsAsFactors = FALSE
  )
)
stopifnot(
  wrapper_env$manifest_matches(valid_manifest, absent_dependencies),
  identical(
    wrapper_env$manifest_computational_dependencies(valid_manifest),
    absent_dependencies
  ),
  !wrapper_env$manifest_matches(duplicate_manifest, absent_dependencies),
  !wrapper_env$manifest_matches(rogue_manifest, absent_dependencies),
  is.null(wrapper_env$manifest_computational_dependencies(
    duplicate_manifest
  ))
)

contract_env <- new.env(parent = globalenv())
sys.source(
  file.path(
    repo_root,
    "Code", "in-vivo", "figure7", "src", "common_io.R"
  ),
  envir = contract_env
)
config_path <- file.path(
  repo_root,
  "Code", "in-vivo", "figure7", "figure7_config.yaml"
)
config <- contract_env$figure7_read_config(config_path)
baseline_contract <- contract_env$figure7_si_raw_config_contract_sha256(config)
stopifnot(
  identical(
    contract_env$figure7_si_raw_fgsea_nperm_simple(config),
    10000L
  ),
  identical(
    unname(
      contract_env$figure7_si_raw_config_contract_values(config)[[
        "si_figures.raw_rebuild_fgsea_nperm_simple"
      ]]
    ),
    "10000"
  )
)
expect_error(
  contract_env$figure7_si_raw_fgsea_nperm_simple(config, "9999"),
  "--fgsea-nperm-simple must equal the reviewed value 10000"
)
changed_nperm <- config
changed_nperm$si_figures$raw_rebuild_fgsea_nperm_simple <- 9999L
expect_error(
  contract_env$figure7_si_raw_config_contract_sha256(changed_nperm),
  "raw_rebuild_fgsea_nperm_simple must equal the reviewed 10000"
)
changed_consumed <- config
changed_consumed$si_figures$ploidy_field <- "changed_ploidy"
stopifnot(!identical(
  baseline_contract,
  contract_env$figure7_si_raw_config_contract_sha256(changed_consumed)
))
changed_renderer_only <- config
changed_renderer_only$si_figures$plot_shuffle_seed <- 999L
stopifnot(identical(
  baseline_contract,
  contract_env$figure7_si_raw_config_contract_sha256(changed_renderer_only)
))
changed_unrelated <- config
changed_unrelated$tgi$day <- 99L
stopifnot(identical(
  baseline_contract,
  contract_env$figure7_si_raw_config_contract_sha256(changed_unrelated)
))
stopifnot(identical(
  wrapper_env$computational_fingerprint(absent_dependencies),
  contract_env$figure7_si_raw_work_fingerprint(absent_dependencies)
))

# A generated raw-table cache is its own deepest checkpoint. Missing source
# artifacts/manifests are attested by its strict input manifest; every present
# explicit/final/partial/H5 source must still match.
sys.source(
  file.path(
    repo_root,
    "Code", "in-vivo", "figure7", "src", "common_io.R"
  ),
  envir = wrapper_env
)
sys.source(
  file.path(
    repo_root,
    "Code", "in-vivo", "figure7", "src",
    "seurat_upstream_selection.R"
  ),
  envir = wrapper_env
)
candidate_source <- tempfile("si-source-", fileext = ".rds")
writeLines("selected reconstructed source", candidate_source)
candidate_source_hash <- wrapper_env$sha256(candidate_source)
base_candidate_dependencies <- absent_dependencies[
  setdiff(names(absent_dependencies), "source_seurat_rds")
]
reconstruction_roles <- c(
  "integrated_rds",
  "cell_cycle_rds",
  "refined_rds",
  "merged_rds",
  "all_ploidy",
  "sample_info",
  "h5_inventory"
)
science_roles <- paste0(
  "seurat_scientific_code_contract:",
  c("integrated", "cell_cycle", "refined", "merged", "final")
)
reconstruction_contract <- list(
  transitive_roles = paste0(
    "seurat_transitive_dependency:",
    reconstruction_roles
  ),
  scientific_roles = science_roles,
  fixed = c(
    "seurat_transitive_dependency:all_ploidy" = strrep("a", 64L),
    "seurat_transitive_dependency:sample_info" = strrep("b", 64L),
    stats::setNames(
      vapply(seq_along(science_roles), function(i) {
        strrep(sprintf("%x", i + 1L), 64L)
      }, character(1L)),
      science_roles
    ),
    seurat_upstream_scientific_common_contract = strrep("c", 64L)
  )
)
h5_inventory_hash <- strrep("d", 64L)
source_dependencies <- c(
  source_seurat_rds = candidate_source_hash,
  stats::setNames(
    c(
      strrep("e", 64L),
      strrep("f", 64L),
      strrep("1", 64L),
      strrep("2", 64L),
      unname(reconstruction_contract$fixed[[
        "seurat_transitive_dependency:all_ploidy"
      ]]),
      unname(reconstruction_contract$fixed[[
        "seurat_transitive_dependency:sample_info"
      ]]),
      h5_inventory_hash
    ),
    reconstruction_contract$transitive_roles
  ),
  reconstruction_contract$fixed[science_roles],
  reconstruction_contract$fixed[[
    "seurat_upstream_scientific_common_contract"
  ]]
)
names(source_dependencies)[length(source_dependencies)] <-
  "seurat_upstream_scientific_common_contract"
candidate_dependencies <- c(
  base_candidate_dependencies,
  source_dependencies
)
candidate_fingerprint <- wrapper_env$computational_fingerprint(
  candidate_dependencies
)
candidate_root <- tempfile("si-generated-cache-root-")
candidate_dir <- file.path(
  candidate_root,
  paste0("legacy_mixed_", substr(candidate_fingerprint, 1L, 20L))
)
dir.create(file.path(candidate_dir, "metadata"), recursive = TRUE)
candidate_audit_roles <- c(
  "audit_si_table_builder",
  "audit_environment_validator",
  "audit_environment_lock",
  "audit_figure7_config",
  "audit_seurat_selection",
  "audit_seurat_upstream_generator",
  "audit_seurat_upstream_helper",
  "audit_seurat_final_stage_manifest",
  "audit_seurat_reconstruction_manifest"
)
candidate_manifest <- data.frame(
  role = c(names(candidate_dependencies), candidate_audit_roles),
  locator = c(
    paste0("contract:", names(candidate_dependencies)),
    paste0("audit:", candidate_audit_roles)
  ),
  sha256 = c(
    unname(candidate_dependencies),
    vapply(seq_along(candidate_audit_roles), function(index) {
      strrep(sprintf("%x", index), 64L)
    }, character(1L))
  ),
  bytes = c(
    rep(0, length(candidate_dependencies)),
    rep(123, length(candidate_audit_roles))
  ),
  stringsAsFactors = FALSE
)
utils::write.table(
  candidate_manifest,
  file.path(candidate_dir, "metadata", "input_manifest.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
candidate_run_config <- c(
  schema_version = "1",
  module = "si_figures_raw_table_build",
  output_table_count = "11",
  figures_supported = "4,5,6,7",
  si7_canonical_publication_allowed = "false",
  seurat_source_kind = "manifested_reconstruction",
  si_raw_scientific_code_contract_sha256 =
    unname(base_candidate_dependencies[[
      "si_raw_scientific_code_contract"
    ]]),
  si_raw_config_contract_sha256 =
    unname(base_candidate_dependencies[["si_raw_config_contract"]]),
  si_raw_runtime_contract_sha256 =
    unname(base_candidate_dependencies[["si_raw_runtime_contract"]]),
  work_dependency_fingerprint = candidate_fingerprint
)
write_candidate_run_config <- function(values = candidate_run_config) {
  utils::write.table(
    data.frame(
      key = names(values),
      value = unname(values),
      stringsAsFactors = FALSE
    ),
    file.path(candidate_dir, "metadata", "run_config.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}
write_candidate_run_config()
fake_generator <- new.env(parent = emptyenv())
fake_generator$figure7_upstream_stage_definitions <- function(root) {
  stats::setNames(
    lapply(
      c("integrated", "cell_cycle", "refined", "merged"),
      function(stage) list(path = file.path(root, paste0(stage, ".rds")))
    ),
    c("integrated", "cell_cycle", "refined", "merged")
  )
}
fake_generator$figure7_upstream_stage_output <- function(definition) {
  definition$path
}
fake_generator$figure7_upstream_read_sample_info <- function(path) {
  data.frame(id = "ignored")
}
fake_generator$figure7_upstream_h5_inventory <- function(...) {
  data.frame(sample = "ignored", sha256 = strrep("0", 64L))
}
fake_generator$figure7_upstream_h5_inventory_digest <- function(...) {
  h5_inventory_hash
}
absent_upstream <- file.path(candidate_root, "archived-upstream")
stale_h5 <- file.path(candidate_root, "archived-cellranger")
candidate_result <- wrapper_env$si_cache_candidate_dependencies(
  candidate_dir,
  base_candidate_dependencies,
  strrep("0", 64L),
  reconstruction_contract,
  fake_generator,
  explicit_rds = "",
  deposited_rds = "",
  upstream_dir = absent_upstream,
  cellranger_root = stale_h5,
  sample_info = "/archived/sample_info.xlsx"
)
stopifnot(
  !is.null(candidate_result),
  identical(candidate_result$source_hash, candidate_source_hash)
)
stopifnot(isTRUE(wrapper_env$si_cache_candidate_dependencies(
  candidate_dir,
  base_candidate_dependencies,
  candidate_source_hash,
  reconstruction_contract,
  fake_generator,
  upstream_dir = absent_upstream,
  cellranger_root = stale_h5,
  sample_info = "/archived/sample_info.xlsx"
)$reconstructed))

# Audit hashes do not affect reuse.
candidate_manifest$sha256[
  candidate_manifest$role == "audit_si_table_builder"
] <- strrep("9", 64L)
utils::write.table(
  candidate_manifest,
  file.path(candidate_dir, "metadata", "input_manifest.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
stopifnot(!is.null(wrapper_env$si_cache_candidate_dependencies(
  candidate_dir,
  base_candidate_dependencies,
  strrep("0", 64L),
  reconstruction_contract,
  fake_generator,
  upstream_dir = absent_upstream,
  cellranger_root = stale_h5,
  sample_info = "/archived/sample_info.xlsx"
)))

# A run-scoped optional scVelo audit replaces any historical metrics row
# without mutating the reusable cache manifest.
old_metrics <- tempfile("si-old-metrics-", fileext = ".csv")
new_metrics <- tempfile("si-new-metrics-", fileext = ".csv")
writeLines("cell,value\nold,1", old_metrics)
writeLines("cell,value\nnew,2", new_metrics)
historical_manifest <- rbind(
  candidate_manifest,
  data.frame(
    role = "audit_optional_scvelo",
    locator = "external:old.csv",
    sha256 = wrapper_env$sha256(old_metrics),
    bytes = as.numeric(file.info(old_metrics)$size),
    stringsAsFactors = FALSE
  )
)
historical_manifest_path <- tempfile(
  "si-historical-manifest-",
  fileext = ".tsv"
)
utils::write.table(
  historical_manifest,
  historical_manifest_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
augmented_manifest_path <- tempfile(
  "si-augmented-manifest-",
  fileext = ".tsv"
)
wrapper_env$augment_manifest_with_scvelo_audit(
  historical_manifest_path,
  new_metrics,
  augmented_manifest_path,
  repo_root
)
augmented_manifest <- wrapper_env$read_manifest(
  augmented_manifest_path
)
metric_rows <- augmented_manifest$role == "audit_optional_scvelo"
stopifnot(
  sum(metric_rows) == 1L,
  identical(
    augmented_manifest$sha256[metric_rows],
    wrapper_env$sha256(new_metrics)
  ),
  identical(
    wrapper_env$read_manifest(historical_manifest_path)$sha256[
      historical_manifest$role == "audit_optional_scvelo"
    ],
    wrapper_env$sha256(old_metrics)
  )
)

# Present explicit source, H5 inventory, and partial ancestor checks.
stopifnot(!is.null(wrapper_env$si_cache_candidate_dependencies(
  candidate_dir,
  base_candidate_dependencies,
  strrep("0", 64L),
  reconstruction_contract,
  fake_generator,
  explicit_rds = candidate_source,
  upstream_dir = absent_upstream,
  cellranger_root = stale_h5,
  sample_info = "/archived/sample_info.xlsx"
)))
changed_source <- tempfile("si-changed-source-", fileext = ".rds")
writeLines("changed source", changed_source)
stopifnot(is.null(wrapper_env$si_cache_candidate_dependencies(
  candidate_dir,
  base_candidate_dependencies,
  strrep("0", 64L),
  reconstruction_contract,
  fake_generator,
  explicit_rds = changed_source,
  upstream_dir = absent_upstream,
  cellranger_root = stale_h5,
  sample_info = "/archived/sample_info.xlsx"
)))
dir.create(stale_h5)
stopifnot(!is.null(wrapper_env$si_cache_candidate_dependencies(
  candidate_dir,
  base_candidate_dependencies,
  strrep("0", 64L),
  reconstruction_contract,
  fake_generator,
  upstream_dir = absent_upstream,
  cellranger_root = stale_h5,
  sample_info = "/archived/sample_info.xlsx"
)))
h5_inventory_hash <- strrep("9", 64L)
stopifnot(is.null(wrapper_env$si_cache_candidate_dependencies(
  candidate_dir,
  base_candidate_dependencies,
  strrep("0", 64L),
  reconstruction_contract,
  fake_generator,
  upstream_dir = absent_upstream,
  cellranger_root = stale_h5,
  sample_info = "/archived/sample_info.xlsx"
)))
h5_inventory_hash <- strrep("d", 64L)
dir.create(absent_upstream)
writeLines("changed ancestor", file.path(absent_upstream, "integrated.rds"))
stopifnot(is.null(wrapper_env$si_cache_candidate_dependencies(
  candidate_dir,
  base_candidate_dependencies,
  strrep("0", 64L),
  reconstruction_contract,
  fake_generator,
  upstream_dir = absent_upstream,
  cellranger_root = "",
  sample_info = "/archived/sample_info.xlsx"
)))

# Scientific/config/runtime and run-config provenance changes invalidate.
changed_base <- base_candidate_dependencies
changed_base[["si_raw_config_contract"]] <- strrep("9", 64L)
stopifnot(is.null(wrapper_env$si_cache_candidate_dependencies(
  candidate_dir,
  changed_base,
  strrep("0", 64L),
  reconstruction_contract,
  fake_generator,
  upstream_dir = file.path(candidate_root, "missing-again"),
  cellranger_root = "",
  sample_info = "/archived/sample_info.xlsx"
)))
contradictory_run_config <- candidate_run_config
contradictory_run_config[["seurat_source_kind"]] <- "deposited"
write_candidate_run_config(contradictory_run_config)
stopifnot(is.null(wrapper_env$si_cache_candidate_dependencies(
  candidate_dir,
  base_candidate_dependencies,
  strrep("0", 64L),
  reconstruction_contract,
  fake_generator,
  upstream_dir = file.path(candidate_root, "missing-again"),
  cellranger_root = "",
  sample_info = "/archived/sample_info.xlsx"
)))
write_candidate_run_config()

fallback_arguments <- wrapper_env$raw_builder_args(
  "/input/all_ploidy.tsv",
  "/input/deposited.rds",
  "/input/config.yaml",
  "/output/staging"
)
stopifnot(
  identical(fallback_arguments[["seurat-rds"]], "/input/deposited.rds"),
  is.null(fallback_arguments[["scvelo-metrics"]])
)
audit_arguments <- wrapper_env$raw_audit_args(
  "/input/all_ploidy.tsv",
  "/input/deposited.rds",
  "/input/config.yaml",
  "/audit/metrics.csv"
)
stopifnot(
  identical(audit_arguments[["seurat-rds"]], "/input/deposited.rds"),
  identical(audit_arguments[["scvelo-metrics"]], "/audit/metrics.csv"),
  identical(audit_arguments[["audit-only"]], "true"),
  is.null(audit_arguments[["output-dir"]])
)

default_offline <- wrapper_env$raw_rds_acquisition(
  "",
  "/cache",
  "deposited.rds",
  allow_download = FALSE
)
explicit_rds <- wrapper_env$raw_rds_acquisition(
  "/read-only/deposited.rds",
  "/cache",
  "deposited.rds",
  allow_download = TRUE
)
stopifnot(
  identical(default_offline$path, "/cache/deposited.rds"),
  identical(default_offline$verify_with_downloader, TRUE),
  identical(default_offline$allow_download, FALSE),
  identical(explicit_rds$path, "/read-only/deposited.rds"),
  identical(explicit_rds$verify_with_downloader, FALSE)
)
stopifnot(
  !wrapper_env$raw_rds_needed(TRUE),
  wrapper_env$raw_rds_needed(FALSE),
  wrapper_env$raw_rds_needed(TRUE, "/audit/metrics.csv")
)

# Prove the valid canonical-cache branch returns before resolving, checking, or
# downloading the deliberately missing RDS.
wrapper_env$script_path <- function() wrapper_path
wrapper_env$validate_cache <- function(validator, cache_dir, legacy = FALSE) {
  !legacy
}
wrapper_env$rendered <- NULL
wrapper_env$render_cache <- function(
  renderer,
  cache,
  config,
  output,
  legacy,
  log_name = "00_render.log"
) {
  wrapper_env$rendered <- list(
    cache = cache,
    output = output,
    legacy = legacy
  )
  invisible(output)
}
output_dir <- file.path(tempdir(), "si-direct-rds-contract-output")
wrapper_env$main(list(
  mode = "plot-only",
  "output-dir" = output_dir,
  "table-cache-dir" = file.path(tempdir(), "reviewed-cache"),
  "seurat-rds" = "/definitely/missing/deposited.rds",
  "download-missing-raw" = "false"
))
stopifnot(
  !is.null(wrapper_env$rendered),
  identical(wrapper_env$rendered$legacy, FALSE),
  identical(wrapper_env$rendered$output, normalizePath(
    output_dir,
    mustWork = FALSE
  ))
)

blocked_intermediate <- file.path(
  tempdir(),
  paste0("si-plot-only-must-not-fallback-", Sys.getpid())
)
blocked_output <- paste0(blocked_intermediate, "-output")
wrapper_env$validate_cache <- function(validator, cache_dir, legacy = FALSE) {
  FALSE
}
wrapper_env$run_stage <- function(...) {
  stop("Raw stage must not run in plot-only mode", call. = FALSE)
}
expect_error(
  wrapper_env$main(list(
    mode = "plot-only",
    "output-dir" = blocked_output,
    "table-cache-dir" = file.path(blocked_intermediate, "missing-cache"),
    "intermediate-dir" = blocked_intermediate,
    "seurat-rds" = "/definitely/missing/deposited.rds"
  )),
  "plot-only mode will not acquire raw data"
)
stopifnot(
  !dir.exists(blocked_intermediate),
  !dir.exists(blocked_output)
)

expect_error(
  wrapper_env$main(list(
    mode = "unknown",
    "output-dir" = output_dir
  )),
  "--mode must be auto, plot-only, or full-workflow"
)

message("SI direct-RDS and optional scVelo audit contract tests passed")
