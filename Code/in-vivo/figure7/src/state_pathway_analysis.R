# Retired optional full-analysis guard.
# Corrected panel-7F recomputation is implemented by the cache-aware
# full-workflow support/export path, not this older artifact-driven entrypoint.

figure7_full_analysis <- function(seurat_rds, gene_set_artifact, reference, output_dir, config) {
  if (is.null(seurat_rds) || !nzchar(seurat_rds) || !file.exists(seurat_rds)) {
    figure7_stop("Full analysis requires an existing explicit --seurat-rds")
  }
  if (is.null(gene_set_artifact) || !nzchar(gene_set_artifact) || !file.exists(gene_set_artifact)) {
    figure7_stop("Full analysis requires an existing pinned --gene-set-artifact; live msigdbr queries are forbidden")
  }
  invisible(vapply(c(seurat_rds, gene_set_artifact), figure7_sha256, character(1L)))
  figure7_stop(
    "This older --mode=full-analysis entrypoint is retired. Use full-workflow ",
    "to recompute panel 7F with exact GRCh-only counts and a run-scoped, ",
    "noncanonical generated reference."
  )
}
