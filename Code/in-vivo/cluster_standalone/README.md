# Standalone in-vivo cluster workflow

This directory contains the shortest standalone path from the three accepted
input sources to the refined and final Seurat objects. It intentionally omits
annotation, ORA/GSEA, Dose/Ploidy/TN differential expression, trajectory, and
other downstream analyses.

## Public entry point

```bash
bash run_cluster_standalone.sh INPUT_DIR OUTPUT_DIR
```

`INPUT_DIR` must contain exactly the required source locations:

```text
INPUT_DIR/
├── .../
│   └── A02_cellRanger/
│       └── *-Count-HM/outs/filtered_feature_bc_matrix.h5
├── all_ploidy.tsv
└── sample_info.xlsx
```

`A02_cellRanger` may be nested below `INPUT_DIR` (for example under a model
name), but exactly one such directory must be present. The two auxiliary files
must be directly under `INPUT_DIR`.

The workflow does not read any external pre-existing Seurat object, DEG table,
QC flag table, annotation result, or gene-set result. A retry with the same
`OUTPUT_DIR` may resume from its own initial, cell-cycle, refined, manual-merge,
and final RDS only after the corresponding cell-order, metadata, label,
cluster-count, assay, reduction, and command contracts pass. A verified final
RDS resume runs only missing final cluster-vs-rest DEG outputs and never
overwrites the object. The workflow may also reuse its own validated DEG chunk
caches and complete DEG CSV files.
Paths are supplied only through the two positional arguments above; no
machine-specific input or output path is embedded in the source.

Before analysis, the entry point creates a runtime-specific private library
under `OUTPUT_DIR/cluster_standalone_r_library/R-<version>-<platform>` and
installs seven vendored source dependencies: RcppAnnoy 0.0.22, irlba 2.3.5.1,
sctransform 0.4.2, uwot 0.2.3, xgboost 1.7.11.1,
BiocNeighbors 2.2.0, and assorthead 1.2.0.
This requires no network access. The assorthead source contains an explicit
Annoy floating-point accumulation order. The sctransform source contains
macOS arm64 `libsystem_m` compatibility routines for the sparse geometric mean
and `log10`, plus the R 4.5.1 Apple Silicon fused-multiply-add path used by
the density-weighted step-1 gene sampler. This preserves the original
SCTransform step-1 genes and variable features without storing a gene
whitelist. The two cell-culture PCA matrices
are oriented by their maximum-absolute gene loading and rounded to 10 decimal
places. These platform-stability rules were verified to retain exactly the
original 16,050 2N and 12,709 4N singlet barcode sets on both the source Mac
and the target Linux container; they contain no barcode whitelist or cluster
deletion list.

The vendored uwot source implements the Apple `libsystem_m` single-precision
`powf` path used by the reference macOS run. UMAP is executed in a small
subprocess with the container's reference BLAS/LAPACK preloaded and uwot's
irlba spectral initializer forced, matching the R 4.5.1 macOS reference rather
than the container's default OpenBLAS/RSpectra path. The private package and
subprocess do not replace or modify any package inside the SIF. Given the exact
reference PCA matrix, this helper reproduces the 42,884-cell reference initial
UMAP with zero numeric difference. The Linux integration/PCA matrix itself is
not bitwise identical to the Apple Accelerate result, so this helper contract
must not be interpreted as an end-to-end claim that every stored floating-point
matrix is byte-identical.

Because irlba component signs are mathematically arbitrary and differed between
Apple Accelerate and the Linux reference BLAS, the stored 50-component initial
and final PCA reductions are oriented to the reference sign convention using
the sign of each component's maximum-absolute gene loading. This contract
contains only 50 signs per PCA stage, not gene names, barcodes, or cluster
assignments. After sign alignment, the observed platform residual was about
`3.22e-7` for refined PCA coordinates and `2.41e-7` for final PCA coordinates;
those remaining differences reflect BLAS arithmetic and are not hidden by the
validator.

For the platform-sensitive 4c/9c boundary only, the pipeline creates three
auxiliary compatibility UMAPs. It orients each PCA component by requiring the
cell with the largest absolute coordinate to be positive, then quantizes the
oriented PCA at fixed steps. Cluster 4 uses the intersection of the standard
refinement calls at steps `0.0125` and `0.03`; cluster 9 uses the union of the
calls at steps `0.03` and `0.09`. All three quantized PCA matrices were
verified to be element-for-element identical between the reference macOS PCA
and the standalone Linux PCA. The resulting sets reproduce the reference
4c=103 and 9c=409 cell IDs exactly, without a barcode/cell whitelist. These
auxiliary embeddings are used only for refinement membership; the original
single main UMAP remains the UMAP stored in the Seurat object.

## Analysis order

1. Read the `filtered_feature_bc_matrix.h5` files and the two auxiliary files.
2. Apply the original cell-culture QC/doublet rules and the `all_ploidy.tsv`
   barcode whitelist for tumor samples.
3. Run per-sample SCTransform, SCT-CCA integration, PCA, graph clustering, and
   the initial UMAP.
4. Score the cell cycle and add the original cluster-level cell-cycle flag.
5. On all initial clusters, run cluster-vs-rest DEG and calculate signed-LFC
   cosine similarity under eight feature definitions. A filesystem gate must
   pass before any fixed merge label can be created.
6. Reproduce the UMAP-neighborhood split of clusters 4 and 9 into 4c and 9c.
   Use the three platform-stable auxiliary UMAP calls described above to resolve
   only the numerical boundary cells, then save
   `integrated_sct_cca_seurat_cluster_refine.rds`.
7. Compute refined-cluster QC flags without adding QC columns to the saved
   Seurat objects.
8. After the expression-similarity gate, apply the fixed merges 0/1/7 -> 0 and
   10/11/12 -> 10.
9. On all pre-filter cells, run merged/refined cluster-vs-rest DEG. Select
   clusters dynamically when either `concern_level == "High"` or the cluster
   has zero qualifying upregulated genes.
10. Assert that the selected set is the reference set `3,4,9,9c`; this is a
    reproduction check only and is never used as a fallback deletion rule.
11. Subset, retain the merge labels as `clusters`, derive the reference
    `Ploidy` (`2N`/`4N`, with the reference's unused `Unknown` level) and `TN`
    (`CellLine`/`Tumor`) factor metadata from `IDs`, recompute PCA/UMAP, save
    the final object, and run final
    cluster-vs-rest DEG.

Independent cluster-vs-rest Wilcoxon comparisons are decomposed into about 62
`(cluster, feature-chunk)` tasks and use up to 8 Unix fork workers on HPC. The
allocation provides 64 CPUs but has a 512 GiB memory cgroup limit. A 56-worker
run recorded 39 OOM-killed calls, and a subsequent 32-worker run reached
`MaxRSS=1,116,308,464 KiB` before Slurm marked the step `OUT_OF_MEMORY`; active
workers used approximately 39--42 GiB each. The 8-worker cap therefore keeps
substantial memory headroom while preserving the same deterministic 62-chunk
decomposition; all 62 allocated CPUs remain available to the Slurm step.
Chunk-level Bonferroni correction still uses the full RNA-assay feature count;
the reassembled CSV was verified byte-for-byte against the sequential result.
Every worker explicitly receives the registered stage seed; fork RNG stream
generation is disabled. Completed feature chunks are cached until their cluster
CSV is assembled, and completed DEG CSV files are validated and reused on a
retry. Before each large RDS write, the complete cluster-count table is printed
to the live log; the initial, refined, and final count contracts must pass
before their corresponding object is saved.

Final cluster-vs-rest DEG files use the original `03a_DEGs.R` 17-column schema
(`scope` through `p_val_adj`) and `readr::write_csv` serialization. All nine
files were verified byte-for-byte against the formal reference outputs after
the R 4.5.1 HPC run.

The qualifying upregulated-gene rule is:

```text
p_val_adj < 0.05
avg_log2FC > 0
abs(avg_log2FC) >= 0.25
abs(pct.1 - pct.2) >= 0.05
```

The DEG check is performed before any cluster deletion and therefore covers all
42,884 cells expected in the refined object.

## Main outputs

```text
OUTPUT_DIR/
├── 02b_cluster_refine/objects/
│   └── integrated_sct_cca_seurat_cluster_refine.rds
├── 02c_cluster_quality/summaries/
│   └── cluster_qc_outlier_flags.csv
├── 02e_prefilter_DEGs/summaries/
│   └── data_driven_cluster_removal_decision.csv
├── 03_final_cluster/03_objects/
│   └── integrated_sct_cca_seurat_final_reclustered.rds
└── 03a_DEGs/01_clusters_vs_rest/
    └── cluster_*_vs_rest_markers.csv
```

`run_manifest.tsv`, `seed_registry.tsv`, stage completion markers, and
`sessionInfo.txt` are written under `OUTPUT_DIR/00_provenance`. The same folder
also records source checksums and the private dependency build contract.

## Reproduction validator

After a run, compare the two generated objects against reference objects with:

```bash
Rscript validate_reproduction.R \
  OUTPUT_DIR \
  REFERENCE_REFINED_RDS \
  REFERENCE_FINAL_RDS
```

The validator writes a machine-readable report under
`OUTPUT_DIR/00_validation` and exits nonzero if a required comparison fails.
The analysis entry point never overwrites an existing final RDS; it resumes
only after the final-object contract passes and then fills missing final DEG
outputs.

On the validated R 4.5.1 Linux HPC run, 128 of 144 strict object checks passed.
Cell IDs and order, cluster values and levels, active identities, metadata
schema and values (including `Ploidy` and `TN`), assay dimensions/features, and
Seurat command order all matched the formal macOS reference. The 16 remaining
strict failures are the exact hashes of the SCT/integrated floating-point
matrices, a `1.192092895507812e-7` `scDblFinder.score` residual, PCA coordinate
residuals of `3.214471973045363e-7` (refined) and
`2.408718611790484e-7` (final), non-identical PCA-loading hashes, and the UMAP
coordinates amplified from those platform-specific inputs. The strict report
therefore intentionally remains `FAIL`; the workflow does not weaken the
tolerances or label a non-bitwise object as identical.
