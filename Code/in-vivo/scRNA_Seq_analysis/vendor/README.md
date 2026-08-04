# Vendored R source dependencies

## Purpose and runtime behavior

This directory makes the two-argument standalone workflow independent of
network access and of package versions or builds already present in the SIF.
It does **not** replace the whole R environment in the image and it does not
modify the SIF.

`run_cluster_standalone.sh` creates a runtime-specific private library at:

```text
$OUTPUT_DIR/scRNA_Seq_analysis_r_library/R-<version>-<platform>
```

It prepends that library to `.libPaths()`, then calls
`bootstrap_dependencies.R`. The bootstrap performs the following operations
before any analysis begins:

1. verifies the MD5 of every archive in this directory;
2. installs or rebuilds the required version in the private library;
3. applies the package-specific build contracts described below;
4. verifies all seven installed versions and loads all seven namespaces; and
5. records the source and build contracts under
   `$OUTPUT_DIR/00_provenance`.

Installing these packages manually into the SIF is therefore unnecessary.
With the current bootstrap, doing so is also redundant: the bootstrap checks
the output-specific private library, not merely the SIF site library.

## Source audit summary

The vendored archives were compared directly with the official CRAN archives
and the Bioconductor 3.21 source repository. Four archives are byte-for-byte
upstream sources. Three archives contain deliberate compatibility changes.

| Archive | Source status | Vendor MD5 | Vendor SHA-256 |
|---|---|---|---|
| `RcppAnnoy_0.0.22.tar.gz` | Unmodified CRAN archive | `8ae334c4634f6fdbc066f6a0e93dcb95` | `9f2121d787c4d3e7beccdd65f5d1de81f31c99d57d5d61ca3cc5af7169dd8f65` |
| `irlba_2.3.5.1.tar.gz` | Unmodified CRAN archive | `f738200d5272c7258ee59f7074dd4a6a` | `2cfe6384fef91c223a9920895ce89496f990d1450d731e44309fdbec2bb5c5cf` |
| `uwot_0.2.3.tar.gz` | Modified numerical-compatibility build | `a66c56f91d39b984e0bf02399441004b` | `7bd5ec304f3c9ee31bfd91f9cb90d021306eb970b69621d58a2fa91c38a7ddfd` |
| `sctransform_0.4.2.tar.gz` | Modified compatibility build | `02e80f35d5be2993cbc1740c628cc05a` | `989e8d3ecd51741d7ed5330abed564cf5440a8d2cfdd380f4e5db977b20bb3c3` |
| `assorthead_1.2.0.tar.gz` | Modified compatibility header | `07ba54164ad32cde38585327e94cafc0` | `70bb57a760a3c4ace0b9b1e9ae0b0bb19175512dc3bdf4f80922082ca162d6bd` |
| `xgboost_1.7.11.1.tar.gz` | Unmodified CRAN archive | `fcb32ea43b53faf9625515a4f5197750` | `c15b631be77ae17fef46166abffb3336738a812c219d82c601bc70664a88d7c5` |
| `BiocNeighbors_2.2.0.tar.gz` | Unmodified Bioconductor archive; custom build contract | `a23b27aadce180b89bd71792608006f4` | `3c5a6b91c7318702e51d3defc47591eee0a3f6babfc4c91a4b28f953e0a2d605` |

The source URLs used for the comparison are:

- CRAN archive directories for
  [RcppAnnoy](https://cran.r-project.org/src/contrib/Archive/RcppAnnoy/),
  [irlba](https://cran.r-project.org/src/contrib/Archive/irlba/),
  [uwot](https://cran.r-project.org/src/contrib/Archive/uwot/),
  [sctransform](https://cran.r-project.org/src/contrib/Archive/sctransform/),
  and [xgboost](https://cran.r-project.org/src/contrib/Archive/xgboost/).
- Bioconductor 3.21 source packages for
  [BiocNeighbors 2.2.0](https://bioconductor.org/packages/3.21/bioc/src/contrib/BiocNeighbors_2.2.0.tar.gz)
  and
  [assorthead 1.2.0](https://bioconductor.org/packages/3.21/bioc/src/contrib/assorthead_1.2.0.tar.gz).

## Package-by-package details

### `RcppAnnoy 0.0.22`

Source changes: **none**. The vendored archive has the same MD5 and SHA-256 as
the official CRAN archive.

Why it is pinned:

- Seurat 4.4.0 uses RcppAnnoy for Annoy nearest-neighbor searches in the
  initial graph construction.
- A fixed version prevents a newer Annoy wrapper or bundled Annoy
  implementation from changing neighbor ordering at floating-point ties.

Build behavior:

- It is installed before any compatibility `R_MAKEVARS_USER` is enabled.
- It therefore uses the normal compiler settings from the R 4.5.1 SIF.
- The bootstrap verifies version `0.0.22` in the private library.

### `irlba 2.3.5.1`

Source changes: **none**. The vendored archive has the same MD5 and SHA-256 as
the official CRAN archive.

Why it is pinned:

- Seurat's approximate `RunPCA()` path uses `irlba` for truncated singular
  value decomposition.
- PCA coordinates feed the neighbor graph, clustering, UMAP, and downstream
  cluster refinement, so the implementation version is part of the numerical
  reproduction contract.

Build behavior:

- It is installed with the normal compiler settings from the SIF, before the
  BiocNeighbors-specific compiler contract is enabled.
- The bootstrap verifies version `2.3.5.1` in the private library.

### `uwot 0.2.3`

This is a modified build of the official CRAN 0.2.3 source.

- Official CRAN archive MD5:
  `62b3c53b64bc841664db019f59362ef3`
- Official CRAN archive SHA-256:
  `eeb9e49e1765db61d462577b0c8426de4c52761180e8104ae70d3631b3672ecc`
- Vendored archive MD5:
  `a66c56f91d39b984e0bf02399441004b`
- Vendored archive SHA-256:
  `7bd5ec304f3c9ee31bfd91f9cb90d021306eb970b69621d58a2fa91c38a7ddfd`

The graph construction, epoch schedule, PCG random-number generator, negative
sampling, gradient equations, and all R-level UMAP defaults remain unchanged.
The functional change is limited to the power calculation used by the standard
UMAP attractive and repulsive gradients.

#### Functional source changes

`inst/include/uwot/gradient.h`

- Includes the new `apple_powf.h` compatibility implementation.
- Changes the standard `umap_gradient` template from host `std::pow` to
  `apple_powf_compat`.
- The approximate-power, t-UMAP, LargeVis, and unrelated gradient classes are
  unchanged.

`inst/include/uwot/apple_powf.h` (new)

- Reconstructs the positive, finite-input path of macOS arm64
  `libsystem_m::powf`, which is the only path reached by UMAP squared
  distances.
- Uses the same 128-entry log2 and exp2 lookup tables, range reduction,
  polynomial order, clamp, round-to-nearest-away conversion, and float return.
- Uses explicit `std::fma()` exactly where the Apple arm64 routine contracts a
  multiply-add. Inputs outside the UMAP path fall back to `std::pow`.

`inst/include/uwot/apple_powf_log_data.inc` and
`inst/include/uwot/apple_powf_exp_data.inc` (new)

- Store the lookup tables and polynomial constants as exact bytes, so compiler
  decimal parsing cannot change a bit.
- The installed `apple_powf.h` MD5 is
  `69f469487fb4d13ada4caa9c5f9473f1`.

`inst/cluster_standalone_contract.txt` (new)

- Records the package version, gradient implementation, spectral backend, and
  the standalone numerical probe result.

The upstream generated `MD5` file was removed because the source files above
are deliberately changed; the bootstrap verifies the outer archive and
installed-header MD5 values instead.

#### Numerical validation and runtime contract

- The reconstructed `powf` was compared with 200,000 macOS reference inputs:
  200,000 bitwise matches, maximum ULP difference 0.
- With the exact formal fuzzy graph and spectral initialization, the patched
  Linux optimizer reproduced all 85,768 final UMAP coordinates bit-for-bit.
- Linux OpenBLAS still changed the last bits of the irlba spectral
  initialization. `run_compatible_umap.R` therefore launches only UMAP in a
  subprocess preloaded with the reference BLAS/LAPACK already shipped in the
  SIF and forces uwot's irlba branch. No SIF file or installed SIF package is
  changed.
- When supplied with the exact formal 42,884-cell PCA matrix, the full initial
  UMAP reproduced the formal macOS embedding bit-for-bit: maximum absolute
  difference 0 and RMSE 0. This validates the UMAP helper; it does not make the
  independently recomputed Linux integration/PCA matrix bit-for-bit identical
  to Apple Accelerate output.

The bootstrap rejects a version-only `uwot 0.2.3` installation unless the
installed contract and compatibility-header MD5 also match.

### `sctransform 0.4.2`

This is a modified build of the official CRAN 0.4.2 source.

- Official CRAN archive MD5:
  `5d84b733622a79d6f6fc53b355a48a6b`
- Official CRAN archive SHA-256:
  `56180aebc7f91345e63210df76541681d3673a59d1a388523dc7ba69886ba0d0`
- Vendored archive MD5:
  `02e80f35d5be2993cbc1740c628cc05a`
- Vendored archive SHA-256:
  `989e8d3ecd51741d7ed5330abed564cf5440a8d2cfdd380f4e5db977b20bb3c3`

The model, residual formula, default arguments, and SCTransform workflow are
unchanged. The modifications make a small set of platform-sensitive numerical
operations reproduce the R 4.5.1 Apple Silicon/macOS path that generated the
reference object.

#### Functional source changes

`R/vst.R`

- Replaces the step-1 gene-sampling sequence
  `density(..., bw = "nrd")` plus `approx(...)` with
  `apple_r451_density_sampling_probability()`.
- The mathematical sampling rule remains
  `1 / (interpolated_density + .Machine$double.eps)`.
- The change controls only the floating-point evaluation of the bandwidth,
  density grid, FFT convolution, interpolation, and therefore the seeded gene
  sample.

`R/apple_compat_math.R` (new)

- Defines a package-internal `log10()` wrapper. It calls the compatibility C++
  implementation and restores the input attributes. It does not replace
  `base::log10()` outside the `sctransform` namespace.
- Implements `apple_r451_var_double()`, a two-pass double-precision variance
  calculation. This avoids Linux R's long-double accumulation changing the
  `bw.nrd` bandwidth by a few ULPs.
- Implements `apple_r451_bw_nrd()` with the R rule
  `1.06 * min(sd, IQR / 1.34) * n^(-1/5)`.
- Implements the R 4.5.1 density algorithm with the same `n >= 512`, power-of-2
  expansion, `cut = 3`, and `ext = 4` structure, while routing the
  platform-sensitive steps through compatibility helpers:
  equal-weight binning, grid construction, Singleton FFT, and linear
  interpolation.
- `apple_r451_density_sampling_probability()` returns the probabilities used
  by the modified `vst()` call.

`src/utils.cpp`

- Adds seven Rcpp entry points:
  `apple_r451_fft_numeric_cpp()`, `apple_current_exp_numeric_cpp()`,
  `apple_current_dnorm_zero_numeric_cpp()`, `apple_r451_seq_numeric_cpp()`,
  `apple_r451_bindist_equal_numeric_cpp()`,
  `apple_r451_approx_linear_numeric_cpp()`, and
  `apple_current_log10_numeric_cpp()`.
- `apple_r451_seq_numeric_cpp()` explicitly uses `std::fma(index, by, from)`
  for interior sequence values.
- `apple_r451_bindist_equal_numeric_cpp()` uses explicit FMA for both weighted
  contributions to adjacent bins.
- `apple_r451_approx_linear_numeric_cpp()` performs the same binary interval
  search as linear interpolation and uses explicit FMA for the interpolated
  value.
- Changes both `row_gmean_dgcmatrix()` and
  `row_gmean_grouped_dgcmatrix()` to use the compatibility `log()` and `exp()`
  implementations instead of the host C library functions.

`src/apple_current_log_compat.h`,
`src/apple_current_log10_compat.h`, and
`src/apple_current_exp_compat.h` (new)

- Store the captured lookup-table constants as exact `uint64_t` bit patterns.
- Reconstruct the relevant macOS arm64 `libsystem_m` range reduction,
  polynomial evaluation, split high/low terms, and rounding order.
- Use explicit `std::fma()` where the Apple implementation contracts a
  multiply-add.
- Fall back to the platform standard-library function for special values or
  inputs outside the implemented normal range.

`src/apple_r451_fft_compat.c` (new)

- Contains R 4.5.1's Singleton mixed-radix FFT core, with package-local symbol
  names `apple_r451_fft_factor()` and `apple_r451_fft_work()`.
- It is compiled inside the package so the target contraction behavior is
  controlled instead of depending on the host R binary's precompiled FFT.

`src/Makevars`

- Adds `PKG_CFLAGS = -mfma -ffp-contract=fast` for the C compatibility code.
- The explicit C++ `std::fma()` calls remain explicit regardless of compiler
  reassociation.

`R/RcppExports.R` and `src/RcppExports.cpp`

- Are regenerated registration wrappers for the seven new C++ entry points.
- They do not add a separate statistical model or user-facing workflow.

Packaging-only differences:

- `DESCRIPTION` changes only the `Packaged:` timestamp/builder.
- `build/partial.rdb` was regenerated during package assembly.
- The upstream archive's generated `MD5` file is not retained in the repacked
  source archive; the outer archive MD5 is enforced by the bootstrap instead.

#### Active versus available compatibility helpers

The current standalone SCTransform path actively uses the package-internal
`log10()`, compatibility `log()`/`exp()` in sparse geometric means, and the
custom density helpers for sequence generation, binning, FFT, and
interpolation. `apple_current_exp_numeric_cpp()` and
`apple_current_dnorm_zero_numeric_cpp()` are registered compatibility helpers,
but the current R density wrapper does not call them directly; it continues to
use `stats::dnorm()` for the kernel values.

#### Installation contract

The bootstrap requires both version `0.4.2` and a marker containing:

```text
sctransform=0.4.2
source_md5=02e80f35d5be2993cbc1740c628cc05a
row_gmean_math=macos_arm64_libsystem_m_log_exp_compat
log10=macos_arm64_libsystem_m_compat
density_sampling=R_4.5.1_apple_silicon_fma_compat
```

A matching version number without this marker is rejected and rebuilt.

### `assorthead 1.2.0`

This is a modified build of the Bioconductor 3.21 source.

- Official archive MD5: `aa08d774cd62ced76fe5cce39c2d7d9e`
- Official archive SHA-256:
  `6d00461796b54ee1da063bc23045a7d745c5b6318ed9c34dbe656645a9407295`
- Vendored archive MD5: `07ba54164ad32cde38585327e94cafc0`
- Vendored archive SHA-256:
  `70bb57a760a3c4ace0b9b1e9ae0b0bb19175512dc3bdf4f80922082ca162d6bd`

Functional changes are confined to
`inst/include/annoy/annoylib.h`:

- `dot()` now processes complete blocks of 16 dimensions by first calculating
  16 products into a `volatile` array and then accumulating those products in
  index order. Remaining dimensions use `std::fma(x[z], y[z], sum)`.
- `euclidean_distance()` now processes complete blocks of 16 dimensions by
  first calculating 16 squared differences into a `volatile` array and then
  accumulating them in index order. Remaining dimensions use
  `std::fma(delta, delta, distance)`.
- The `volatile` block temporaries prevent the x86 compiler from silently
  reassociating or vectorizing the full reduction into a different summation
  tree. The explicit tail FMA reproduces the contraction used by the reference
  ARM/macOS build.

The public Annoy API, tree construction logic, distance definitions, random
number generator, and package R code are unchanged. `DESCRIPTION` differs only
in the `Packaged:` timestamp/builder.

`assorthead` is a header provider (`NeedsCompilation: no`). BiocNeighbors 2.2.0
declares `LinkingTo: Rcpp, assorthead`, so this patched header is compiled into
the private BiocNeighbors installation. The bootstrap verifies the installed
header MD5:

```text
3706ec09e55d95fac8306d48f2bd0ba6
```

### `xgboost 1.7.11.1`

Source changes: **none**. The vendored archive has the same MD5 and SHA-256 as
the official CRAN archive.

Why it is pinned:

- `scDblFinder 1.22.0` imports `xgboost` and uses it for doublet
  classification.
- Doublet calls determine which cells enter SCTransform and all later cluster
  analysis, so the classifier implementation is part of the input-cell
  reproduction contract.

Build behavior:

- `R_MAKEVARS_USER` is explicitly unset before xgboost is installed.
- xgboost is compiled with the SIF's normal compiler defaults, not with the
  Annoy/BiocNeighbors compatibility flags.
- The bootstrap verifies version `1.7.11.1` in the private library.

### `BiocNeighbors 2.2.0`

Source changes in the archive: **none**. The vendored archive has the same MD5
and SHA-256 as the official Bioconductor 3.21 archive.

The installed build is nevertheless intentionally different from a default
installation:

- It is compiled only after the patched `assorthead` header has been installed.
- On x86_64/amd64, the bootstrap writes the following
  `R_MAKEVARS_USER` setting before compilation:

  ```text
  CXX17FLAGS = -g -O2 -fno-tree-vectorize -mfma -ffp-contract=fast
  ```

- `-fno-tree-vectorize` prevents a new reduction tree from replacing the
  accumulation order encoded in `annoylib.h`.
- `-mfma` and `-ffp-contract=fast` make the intended FMA path available for the
  compatibility tail operations.
- On non-x86 targets, native compiler defaults are retained because the
  compatibility arithmetic is explicit in the patched header.

Why it is pinned:

- `scDblFinder 1.22.0` imports BiocNeighbors for nearest-neighbor operations.
- Those neighbor results affect doublet filtering and therefore the exact cell
  population passed into integration and clustering.

The bootstrap requires version `2.2.0` plus a build marker containing the
BiocNeighbors source MD5, installed patched-header MD5, machine identifier, and
the exact Makevars line. A version-only match is rejected if that marker is
missing or different.

## Reproducibility scope

The modified archives contain no cell barcodes, cluster labels, sample names,
or analysis results. They reproduce platform-sensitive numerical evaluation;
they do not hard-code the expected 15 clusters.

With R 4.5.1, Seurat 4.4.0, and these private dependencies, the standalone
initial object was validated against the formal reference object as follows:

- 42,884 cell IDs: identical set and identical order;
- all 42,884 cell-to-cluster labels: identical;
- all 15 per-cluster member-ID sets: identical; and
- mismatched cells: 0.

The full runtime records are written to:

```text
$OUTPUT_DIR/00_provenance/dependency_bootstrap.tsv
$OUTPUT_DIR/00_provenance/uwot_build_contract.txt
$OUTPUT_DIR/00_provenance/sctransform_build_contract.txt
$OUTPUT_DIR/00_provenance/biocneighbors_build_contract.txt
```

The upstream source archives retain their original license files. The added
R 4.5.1 FFT compatibility source retains the R Core copyright and GPL notice
from the original FFT implementation.
