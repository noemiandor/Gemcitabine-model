# Optional full pathway recomputation guard.
# The full state-pathway model is deliberately not invoked until the canonical compact
# reference and its mixed-species feature policy have been scientifically approved.

figure7_full_analysis <- function(seurat_rds, gene_set_artifact, reference, output_dir, config) {
  if (is.null(seurat_rds) || !nzchar(seurat_rds) || !file.exists(seurat_rds)) {
    figure7_stop("Full analysis requires an existing explicit --seurat-rds")
  }
  if (is.null(gene_set_artifact) || !nzchar(gene_set_artifact) || !file.exists(gene_set_artifact)) {
    figure7_stop("Full analysis requires an existing pinned --gene-set-artifact; live msigdbr queries are forbidden")
  }
  invisible(vapply(c(seurat_rds, gene_set_artifact), figure7_sha256, character(1L)))
  figure7_stop(
    "Full panel-7F recomputation is guarded: the canonical state-pathway compact reference is not yet available ",
    "and the mixed human/mouse feature policy has not been approved. No live gene-set fallback or ",
    "reference replacement will be performed."
  )
}
