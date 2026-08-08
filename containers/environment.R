#!/usr/bin/env Rscript

read_versions <- function(path) {
  versions <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("ecosystem", "package", "version")
  if (!identical(names(versions), required)) {
    stop("packages.tsv must contain exactly: ecosystem, package, version", call. = FALSE)
  }
  if (any(!nzchar(versions$ecosystem)) || any(!nzchar(versions$package)) ||
      any(!nzchar(versions$version)) || anyDuplicated(versions[c("ecosystem", "package")])) {
    stop("packages.tsv contains blank or duplicate package records", call. = FALSE)
  }
  versions
}

verify_environment <- function(mode, versions_path) {
  if (!mode %in% c("base", "builder", "full")) {
    stop("verification mode must be base, builder, or full", call. = FALSE)
  }
  architecture <- system2("uname", "-m", stdout = TRUE)
  stopifnot(identical(architecture, "x86_64"))
  stopifnot(getRversion() == package_version("4.5.0"))
  versions <- read_versions(versions_path)
  r_versions <- setNames(
    versions$version[versions$ecosystem == "r"],
    versions$package[versions$ecosystem == "r"]
  )
  stopifnot(requireNamespace("BiocManager", quietly = TRUE))
  stopifnot(requireNamespace("BiocVersion", quietly = TRUE))
  stopifnot(packageVersion("BiocManager") == package_version(r_versions[["BiocManager"]]))
  stopifnot(packageVersion("BiocVersion") == package_version(r_versions[["BiocVersion"]]))
  expected_bioconductor <- sub("\\.0$", "", r_versions[["BiocVersion"]])
  stopifnot(package_version(as.character(BiocManager::version())) == package_version(expected_bioconductor))
  cat("R=", as.character(getRversion()), "\n", sep = "")
  cat("BiocManager=", as.character(packageVersion("BiocManager")), "\n", sep = "")
  cat("Bioconductor=", as.character(BiocManager::version()), "\n", sep = "")

  if (mode != "base") {
    lock <- versions[versions$ecosystem == "r", c("package", "version")]
    if (nrow(lock) != 700L || anyDuplicated(lock$package)) {
      stop("R version lock must contain 700 unique packages", call. = FALSE)
    }
    inventories <- lapply(.libPaths(), function(library) {
      result <- installed.packages(lib.loc = library, noCache = TRUE)
      if (!nrow(result)) return(NULL)
      data.frame(package = result[, "Package"], version = result[, "Version"], stringsAsFactors = FALSE)
    })
    inventory <- do.call(rbind, inventories[!vapply(inventories, is.null, logical(1L))])
    visible <- inventory[!duplicated(inventory$package), , drop = FALSE]
    actual <- visible$version[match(lock$package, visible$package)]
    mismatch <- is.na(actual) | actual != lock$version
    if (any(mismatch)) {
      print(data.frame(package = lock$package[mismatch], expected = lock$version[mismatch], actual = actual[mismatch]), row.names = FALSE)
      stop(sum(mismatch), " locked R package versions do not match", call. = FALSE)
    }
    cat("locked_R_packages=", nrow(lock), "\n", sep = "")

    suppressPackageStartupMessages(library(Signac))
    signac_namespace <- asNamespace("Signac")
    layer_data <- getS3method("LayerData", "default", envir = signac_namespace)
    set_layer_data <- getS3method("LayerData<-", "default", envir = signac_namespace)
    layers <- getS3method("Layers", "default", envir = signac_namespace)
    counts <- Matrix::sparseMatrix(i = c(1L, 2L), j = c(1L, 2L), x = c(1, 2), dims = c(2L, 2L), dimnames = list(c("feature1", "feature2"), c("cell1", "cell2")))
    assay <- SeuratObject::CreateAssayObject(counts = counts)
    stopifnot(identical(dim(layer_data(assay, layer = "counts")), c(2L, 2L)))
    assay <- set_layer_data(assay, layer = "data", value = counts)
    stopifnot("counts" %in% layers(assay), "data" %in% layers(assay))
    cat("signac_v4_compatibility_smoke=PASS\n")
    two_bit_path <- system.file("extdata", "single_sequences.2bit", package = "BSgenome.Hsapiens.UCSC.hg38")
    stopifnot(nzchar(two_bit_path))
    two_bit <- rtracklayer::TwoBitFile(two_bit_path)
    stopifnot(length(GenomeInfoDb::seqinfo(two_bit)) > 0L)
    cat("hg38_twobit_smoke=PASS\n")
  }

  python_code <- paste(c(
    "import csv, importlib, importlib.metadata, platform, sys",
    "with open(sys.argv[1], newline=\"\") as handle:",
    "    rows = list(csv.DictReader(handle, delimiter=\"\t\"))",
    "expected = {row[\"package\"]: row[\"version\"] for row in rows if row[\"ecosystem\"] == \"python\"}",
    "if sys.version_info[:3] != (3, 10, 13): raise SystemExit(f\"Python version mismatch: {platform.python_version()}\")",
    "for package, version in expected.items():",
    "    actual = importlib.metadata.version(package)",
    "    if actual != version: raise SystemExit(f\"{package}: expected {version}, found {actual}\")",
    "imports = {\"anndata\": \"anndata\", \"h5py\": \"h5py\", \"igraph\": \"igraph\", \"leidenalg\": \"leidenalg\", \"loompy\": \"loompy\", \"matplotlib\": \"matplotlib\", \"numpy\": \"numpy\", \"pandas\": \"pandas\", \"scanpy\": \"scanpy\", \"scikit-learn\": \"sklearn\", \"scipy\": \"scipy\", \"scvelo\": \"scvelo\", \"velocyto\": \"velocyto\"}",
    "for module in imports.values(): importlib.import_module(module)",
    "print(f\"Python={platform.python_version()}\")",
    "print(f\"Python_distributions_verified={len(expected)}\")"
  ), collapse = "\n")
  python_output <- system2("python", c("-c", shQuote(python_code), shQuote(versions_path)), stdout = TRUE, stderr = TRUE)
  python_status <- attr(python_output, "status")
  cat(python_output, sep = "\n")
  cat("\n")
  if (!is.null(python_status) && python_status != 0L) stop("Python verification failed", call. = FALSE)

  if (mode == "full") {
    tools <- Sys.which(c("gcc", "g++", "gfortran", "make"))
    if (any(!nzchar(tools))) stop("compiler toolchain is incomplete", call. = FALSE)
    cat("rocker_build_toolchain=PASS\n")
  }
  samtools_output <- system2("samtools", "--version", stdout = TRUE, stderr = TRUE)
  stopifnot(sub("^samtools ", "", samtools_output[[1L]]) == "1.23.1")
  stopifnot(sub("^Using htslib ", "", samtools_output[[2L]]) == "1.23.1")
  cat("samtools=1.23.1\nhtslib=1.23.1\n")
  cat("architecture=", architecture, "\n", sep = "")
  cat("environment_verification=", mode, "_PASS\n", sep = "")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && identical(args[[1L]], "verify")) {
  if (length(args) != 3L) stop("usage: environment.R verify MODE PACKAGES_TSV", call. = FALSE)
  verify_environment(args[[2L]], args[[3L]])
  quit(save = "no", status = 0L)
}
if (length(args) != 3L || !identical(args[[1L]], "install")) {
  stop("usage: environment.R install PACKAGES_TSV INSTALL_LIBRARY", call. = FALSE)
}

lock_path <- args[[2L]]
install_library <- args[[3L]]
stopifnot(getRversion() == package_version("4.5.0"))
options(timeout = 3600, pak.no_extra_messages = TRUE, Ncpus = 4L)
versions <- read_versions(lock_path)
r_versions <- versions[versions$ecosystem == "r", c("package", "version")]
if (nrow(r_versions) != 700L || anyDuplicated(r_versions$package)) stop("R version lock must contain 700 unique packages", call. = FALSE)
stopifnot(installed.packages()["pak", "Version"] == r_versions$version[r_versions$package == "pak"])
source_urls <-
c(abind = "https://cran.r-project.org/src/contrib/abind_1.4-8.tar.gz",
actuar = "https://cran.r-project.org/src/contrib/actuar_3.3-7.tar.gz",
adbcdrivermanager = "https://cran.r-project.org/src/contrib/adbcdrivermanager_0.23.0-2.tar.gz",
AER = "https://cran.r-project.org/src/contrib/AER_1.2-17.tar.gz",
alabaster.base = "https://bioconductor.org/packages/3.20/bioc/src/contrib/alabaster.base_1.6.1.tar.gz",
alabaster.matrix = "https://bioconductor.org/packages/3.20/bioc/src/contrib/alabaster.matrix_1.6.1.tar.gz",
alabaster.ranges = "https://bioconductor.org/packages/3.20/bioc/src/contrib/alabaster.ranges_1.6.0.tar.gz",
alabaster.schemas = "https://bioconductor.org/packages/3.20/bioc/src/contrib/alabaster.schemas_1.6.0.tar.gz",
alabaster.se = "https://bioconductor.org/packages/3.20/bioc/src/contrib/alabaster.se_1.6.0.tar.gz",
alfakR = "https://github.com/TaoLee0510/alfakR/archive/604ab0c5a71a078b72a24fa598a7e9f992f1138b.tar.gz",
alluvial = "https://cran.r-project.org/src/contrib/alluvial_0.1-2.tar.gz",
annotate = "https://bioconductor.org/packages/3.18/bioc/src/contrib/annotate_1.80.0.tar.gz",
AnnotationDbi = "https://bioconductor.org/packages/3.20/bioc/src/contrib/AnnotationDbi_1.68.0.tar.gz",
AnnotationFilter = "https://bioconductor.org/packages/3.20/bioc/src/contrib/AnnotationFilter_1.30.0.tar.gz",
AnnotationHub = "https://bioconductor.org/packages/3.20/bioc/src/contrib/AnnotationHub_3.14.0.tar.gz",
ape = "https://cran.r-project.org/src/contrib/ape_5.8-1.tar.gz",
aplot = "https://cran.r-project.org/src/contrib/Archive/aplot/aplot_0.2.9.tar.gz",
AsioHeaders = "https://cran.r-project.org/src/contrib/AsioHeaders_1.30.2-1.tar.gz",
askpass = "https://cran.r-project.org/src/contrib/askpass_1.2.1.tar.gz",
assertthat = "https://cran.r-project.org/src/contrib/assertthat_0.2.1.tar.gz",
assorthead = "https://bioconductor.org/packages/3.22/bioc/src/contrib/assorthead_1.4.0.tar.gz",
AUC = "https://cran.r-project.org/src/contrib/AUC_0.3.2.tar.gz",
audio = "https://cran.r-project.org/src/contrib/audio_0.1-12.tar.gz",
babelgene = "https://cran.r-project.org/src/contrib/babelgene_22.9.tar.gz",
backports = "https://cran.r-project.org/src/contrib/Archive/backports/backports_1.5.0.tar.gz",
base64enc = "https://cran.r-project.org/src/contrib/base64enc_0.1-6.tar.gz",
batchelor = "https://bioconductor.org/packages/3.22/bioc/src/contrib/batchelor_1.26.0.tar.gz",
bbmle = "https://cran.r-project.org/src/contrib/bbmle_1.0.25.1.tar.gz",
bdsmatrix = "https://cran.r-project.org/src/contrib/bdsmatrix_1.3-7.tar.gz",
beachmat = "https://bioconductor.org/packages/3.22/bioc/src/contrib/beachmat_2.26.0.tar.gz",
beepr = "https://cran.r-project.org/src/contrib/beepr_2.0.tar.gz",
beeswarm = "https://cran.r-project.org/src/contrib/beeswarm_0.4.0.tar.gz",
bench = "https://cran.r-project.org/src/contrib/bench_1.1.4.tar.gz",
betareg = "https://cran.r-project.org/src/contrib/betareg_3.2-5.tar.gz",
BH = "https://cran.r-project.org/src/contrib/BH_1.90.0-1.tar.gz",
BiasedUrn = "https://cran.r-project.org/src/contrib/BiasedUrn_2.0.12.tar.gz",
biglm = "https://cran.r-project.org/src/contrib/biglm_0.9-3.tar.gz",
Biobase = "https://bioconductor.org/packages/3.22/bioc/src/contrib/Biobase_2.70.0.tar.gz",
BiocFileCache = "https://bioconductor.org/packages/3.20/bioc/src/contrib/BiocFileCache_2.14.0.tar.gz",
BiocGenerics = "https://bioconductor.org/packages/3.22/bioc/src/contrib/BiocGenerics_0.56.0.tar.gz",
BiocIO = "https://bioconductor.org/packages/3.20/bioc/src/contrib/BiocIO_1.16.0.tar.gz",
BiocManager = "https://cran.r-project.org/src/contrib/BiocManager_1.30.27.tar.gz",
BiocNeighbors = "https://bioconductor.org/packages/3.22/bioc/src/contrib/BiocNeighbors_2.4.0.tar.gz",
BiocParallel = "https://bioconductor.org/packages/3.22/bioc/src/contrib/BiocParallel_1.44.0.tar.gz",
BiocSingular = "https://bioconductor.org/packages/3.22/bioc/src/contrib/BiocSingular_1.26.1.tar.gz",
BiocVersion = "https://bioconductor.org/packages/3.22/bioc/src/contrib/BiocVersion_3.22.0.tar.gz",
biomaRt = "https://bioconductor.org/packages/3.20/bioc/src/contrib/biomaRt_2.62.1.tar.gz",
Biostrings = "https://bioconductor.org/packages/3.20/bioc/src/contrib/Biostrings_2.74.1.tar.gz",
bit = "https://cran.r-project.org/src/contrib/bit_4.6.0.tar.gz",
bit64 = "https://cran.r-project.org/src/contrib/bit64_4.8.2.tar.gz",
bitops = "https://cran.r-project.org/src/contrib/bitops_1.0-9.tar.gz",
blob = "https://cran.r-project.org/src/contrib/Archive/blob/blob_1.2.4.tar.gz",
bluster = "https://bioconductor.org/packages/3.20/bioc/src/contrib/bluster_1.16.0.tar.gz",
bookdown = "https://cran.r-project.org/src/contrib/bookdown_0.47.tar.gz",
BPCells = "https://github.com/bnprks/BPCells/archive/28759cdd512578b6cbe549e226e1cd52a2d2308c.tar.gz",
brew = "https://cran.r-project.org/src/contrib/brew_1.0-10.tar.gz",
brio = "https://cran.r-project.org/src/contrib/brio_1.1.5.tar.gz",
broom = "https://cran.r-project.org/src/contrib/Archive/broom/broom_1.0.8.tar.gz",
broom.mixed = "https://cran.r-project.org/src/contrib/broom.mixed_0.2.9.7.tar.gz",
BSgenome = "https://bioconductor.org/packages/3.20/bioc/src/contrib/BSgenome_1.74.0.tar.gz",
BSgenome.Hsapiens.UCSC.hg38 = "https://bioconductor.org/packages/3.22/data/annotation/src/contrib/BSgenome.Hsapiens.UCSC.hg38_1.4.5.tar.gz",
bslib = "https://cran.r-project.org/src/contrib/bslib_0.11.0.tar.gz",
bspm = "https://cran.r-project.org/src/contrib/Archive/bspm/bspm_0.5.7.tar.gz",
btergm = "https://cran.r-project.org/src/contrib/btergm_1.11.1.tar.gz",
ca = "https://cran.r-project.org/src/contrib/ca_0.71.1.tar.gz",
cachem = "https://cran.r-project.org/src/contrib/cachem_1.1.0.tar.gz",
Cairo = "https://cran.r-project.org/src/contrib/Cairo_1.7-0.tar.gz",
callr = "https://cran.r-project.org/src/contrib/Archive/callr/callr_3.7.6.tar.gz",
car = "https://cran.r-project.org/src/contrib/Archive/car/car_3.1-3.tar.gz",
carData = "https://cran.r-project.org/src/contrib/Archive/carData/carData_3.0-5.tar.gz",
caret = "https://cran.r-project.org/src/contrib/caret_7.0-1.tar.gz",
caTools = "https://cran.r-project.org/src/contrib/Archive/caTools/caTools_1.18.3.tar.gz",
celldex = "https://bioconductor.org/packages/3.20/data/experiment/src/contrib/celldex_1.16.0.tar.gz",
cellranger = "https://cran.r-project.org/src/contrib/cellranger_1.1.0.tar.gz",
checkmate = "https://cran.r-project.org/src/contrib/checkmate_2.3.4.tar.gz",
circlize = "https://cran.r-project.org/src/contrib/Archive/circlize/circlize_0.4.16.tar.gz",
classInt = "https://cran.r-project.org/src/contrib/classInt_0.4-11.tar.gz",
cli = "https://cran.r-project.org/src/contrib/cli_3.6.6.tar.gz",
clipr = "https://cran.r-project.org/src/contrib/clipr_0.8.1.tar.gz",
clock = "https://cran.r-project.org/src/contrib/Archive/clock/clock_0.7.3.tar.gz",
clue = "https://cran.r-project.org/src/contrib/Archive/clue/clue_0.3-66.tar.gz",
cmprsk = "https://cran.r-project.org/src/contrib/cmprsk_2.2-12.tar.gz",
CNEr = "https://bioconductor.org/packages/3.20/bioc/src/contrib/CNEr_1.42.0.tar.gz",
cobs = "https://cran.r-project.org/src/contrib/cobs_1.3-9-1.tar.gz",
coda = "https://cran.r-project.org/src/contrib/coda_0.19-4.1.tar.gz",
collapse = "https://cran.r-project.org/src/contrib/collapse_2.1.7.tar.gz",
collections = "https://cran.r-project.org/src/contrib/Archive/collections/collections_0.3.8.tar.gz",
colorspace = "https://cran.r-project.org/src/contrib/Archive/colorspace/colorspace_2.1-1.tar.gz",
commonmark = "https://cran.r-project.org/src/contrib/commonmark_2.0.0.tar.gz",
ComplexHeatmap = "https://bioconductor.org/packages/3.20/bioc/src/contrib/ComplexHeatmap_2.22.0.tar.gz",
CompQuadForm = "https://cran.r-project.org/src/contrib/CompQuadForm_1.4.4.tar.gz",
conflicted = "https://cran.r-project.org/src/contrib/conflicted_1.2.0.tar.gz",
copykat = "https://github.com/navinlabcode/copykat/archive/d7d6569ae9e30bf774908301af312f626de4cbd5.tar.gz",
corrplot = "https://cran.r-project.org/src/contrib/corrplot_0.95.tar.gz",
covr = "https://cran.r-project.org/src/contrib/covr_3.6.5.tar.gz",
cowplot = "https://cran.r-project.org/src/contrib/cowplot_1.2.0.tar.gz",
cpp11 = "https://cran.r-project.org/src/contrib/cpp11_0.5.5.tar.gz",
crayon = "https://cran.r-project.org/src/contrib/crayon_1.5.3.tar.gz",
credentials = "https://cran.r-project.org/src/contrib/Archive/credentials/credentials_2.0.2.tar.gz",
crosstalk = "https://cran.r-project.org/src/contrib/crosstalk_1.2.2.tar.gz",
curl = "https://cran.r-project.org/src/contrib/curl_7.1.0.tar.gz",
cyclocomp = "https://cran.r-project.org/src/contrib/cyclocomp_1.1.2.tar.gz",
data.table = "https://cran.r-project.org/src/contrib/data.table_1.18.4.tar.gz",
data.tree = "https://cran.r-project.org/src/contrib/data.tree_1.2.0.tar.gz",
DBI = "https://cran.r-project.org/src/contrib/DBI_1.3.0.tar.gz",
DBItest = "https://cran.r-project.org/src/contrib/DBItest_1.8.2.tar.gz",
dbplyr = "https://cran.r-project.org/src/contrib/Archive/dbplyr/dbplyr_2.5.0.tar.gz",
dbscan = "https://cran.r-project.org/src/contrib/Archive/dbscan/dbscan_1.2.2.tar.gz",
decor = "https://cran.r-project.org/src/contrib/decor_1.0.2.tar.gz",
DelayedArray = "https://bioconductor.org/packages/3.22/bioc/src/contrib/DelayedArray_0.36.1.tar.gz",
DelayedMatrixStats = "https://bioconductor.org/packages/3.22/bioc/src/contrib/DelayedMatrixStats_1.32.0.tar.gz",
deldir = "https://cran.r-project.org/src/contrib/deldir_2.0-4.tar.gz",
dendextend = "https://cran.r-project.org/src/contrib/dendextend_1.19.1.tar.gz",
DEoptim = "https://cran.r-project.org/src/contrib/DEoptim_2.2-8.tar.gz",
DEoptimR = "https://cran.r-project.org/src/contrib/DEoptimR_1.2-0.tar.gz",
Deriv = "https://cran.r-project.org/src/contrib/Archive/Deriv/Deriv_4.1.6.tar.gz",
desc = "https://cran.r-project.org/src/contrib/desc_1.4.3.tar.gz",
deSolve = "https://cran.r-project.org/src/contrib/Archive/deSolve/deSolve_1.40.tar.gz",
devtools = "https://cran.r-project.org/src/contrib/devtools_2.5.2.tar.gz",
dfidx = "https://cran.r-project.org/src/contrib/dfidx_0.2-0.tar.gz",
diagram = "https://cran.r-project.org/src/contrib/diagram_1.6.5.tar.gz",
diffobj = "https://cran.r-project.org/src/contrib/Archive/diffobj/diffobj_0.3.6.tar.gz",
digest = "https://cran.r-project.org/src/contrib/digest_0.6.39.tar.gz",
distributional = "https://cran.r-project.org/src/contrib/Archive/distributional/distributional_0.7.0.tar.gz",
DistributionUtils = "https://cran.r-project.org/src/contrib/DistributionUtils_0.6-2.tar.gz",
distro = "https://cran.r-project.org/src/contrib/distro_0.1.1.tar.gz",
dlm = "https://cran.r-project.org/src/contrib/dlm_1.1-6.1.tar.gz",
doBy = "https://cran.r-project.org/src/contrib/Archive/doBy/doBy_4.6.27.tar.gz",
docopt = "https://cran.r-project.org/src/contrib/docopt_0.7.2.tar.gz",
doFuture = "https://cran.r-project.org/src/contrib/doFuture_1.2.2.tar.gz",
doParallel = "https://cran.r-project.org/src/contrib/doParallel_1.0.17.tar.gz",
doSNOW = "https://cran.r-project.org/src/contrib/doSNOW_1.0.20.tar.gz",
dotCall64 = "https://cran.r-project.org/src/contrib/dotCall64_1.2.tar.gz",
downlit = "https://cran.r-project.org/src/contrib/Archive/downlit/downlit_0.4.4.tar.gz",
dplyr = "https://cran.r-project.org/src/contrib/dplyr_1.2.1.tar.gz",
dqrng = "https://cran.r-project.org/src/contrib/dqrng_0.4.1.tar.gz",
drc = "https://cran.r-project.org/src/contrib/drc_3.0-1.tar.gz",
dreamerr = "https://cran.r-project.org/src/contrib/dreamerr_1.5.0.tar.gz",
DT = "https://cran.r-project.org/src/contrib/DT_0.34.0.tar.gz",
dtplyr = "https://cran.r-project.org/src/contrib/Archive/dtplyr/dtplyr_1.3.1.tar.gz",
duckdb = "https://cran.r-project.org/src/contrib/Archive/duckdb/duckdb_1.2.2.tar.gz",
duckdbfs = "https://cran.r-project.org/src/contrib/Archive/duckdbfs/duckdbfs_0.1.0.tar.gz",
e1071 = "https://cran.r-project.org/src/contrib/e1071_1.7-17.tar.gz",
edgeR = "https://bioconductor.org/packages/3.22/bioc/src/contrib/edgeR_4.8.2.tar.gz",
egg = "https://cran.r-project.org/src/contrib/egg_0.4.5.tar.gz",
ellipsis = "https://cran.r-project.org/src/contrib/ellipsis_0.3.3.tar.gz",
emmeans = "https://cran.r-project.org/src/contrib/emmeans_2.0.4.tar.gz",
emojifont = "https://cran.r-project.org/src/contrib/Archive/emojifont/emojifont_0.5.5.tar.gz",
EnhancedVolcano = "https://bioconductor.org/packages/3.20/bioc/src/contrib/EnhancedVolcano_1.24.0.tar.gz",
EnsDb.Hsapiens.v86 = "https://bioconductor.org/packages/3.22/data/annotation/src/contrib/EnsDb.Hsapiens.v86_2.99.0.tar.gz",
ensembldb = "https://bioconductor.org/packages/3.20/bioc/src/contrib/ensembldb_2.30.0.tar.gz",
ergm = "https://cran.r-project.org/src/contrib/ergm_4.12.0.tar.gz",
estimability = "https://cran.r-project.org/src/contrib/estimability_2.0.0.tar.gz",
evaluate = "https://cran.r-project.org/src/contrib/evaluate_1.0.5.tar.gz",
evd = "https://cran.r-project.org/src/contrib/evd_2.3-7.1.tar.gz",
ExperimentHub = "https://bioconductor.org/packages/3.20/bioc/src/contrib/ExperimentHub_2.14.0.tar.gz",
expint = "https://cran.r-project.org/src/contrib/expint_0.2-1.tar.gz",
ExPosition = "https://cran.r-project.org/src/contrib/ExPosition_2.11.0.tar.gz",
fansi = "https://cran.r-project.org/src/contrib/Archive/fansi/fansi_1.0.6.tar.gz",
farver = "https://cran.r-project.org/src/contrib/farver_2.1.2.tar.gz",
fastcluster = "https://cran.r-project.org/src/contrib/fastcluster_1.3.0.tar.gz",
fastDummies = "https://cran.r-project.org/src/contrib/Archive/fastDummies/fastDummies_1.7.5.tar.gz",
fastmap = "https://cran.r-project.org/src/contrib/fastmap_1.2.0.tar.gz",
fastmatch = "https://cran.r-project.org/src/contrib/fastmatch_1.1-8.tar.gz",
fastshap = "https://cran.r-project.org/src/contrib/Archive/fastshap/fastshap_0.1.1.tar.gz",
fftwtools = "https://cran.r-project.org/src/contrib/fftwtools_0.9-11.tar.gz",
fgsea = "https://bioconductor.org/packages/3.22/bioc/src/contrib/fgsea_1.36.2.tar.gz",
fields = "https://cran.r-project.org/src/contrib/Archive/fields/fields_16.3.tar.gz",
filelock = "https://cran.r-project.org/src/contrib/filelock_1.0.3.tar.gz",
fit.models = "https://cran.r-project.org/src/contrib/fit.models_0.64.tar.gz",
fitdistrplus = "https://cran.r-project.org/src/contrib/fitdistrplus_1.2-6.tar.gz",
fixest = "https://cran.r-project.org/src/contrib/fixest_0.14.2.tar.gz",
flexmix = "https://cran.r-project.org/src/contrib/flexmix_2.3-20.tar.gz",
flextable = "https://cran.r-project.org/src/contrib/flextable_0.10.0.tar.gz",
FNN = "https://cran.r-project.org/src/contrib/FNN_1.1.4.1.tar.gz",
fontawesome = "https://cran.r-project.org/src/contrib/fontawesome_0.5.3.tar.gz",
fontBitstreamVera = "https://cran.r-project.org/src/contrib/fontBitstreamVera_0.1.1.tar.gz",
fontLiberation = "https://cran.r-project.org/src/contrib/fontLiberation_0.1.0.tar.gz",
fontquiver = "https://cran.r-project.org/src/contrib/fontquiver_0.2.1.tar.gz",
forcats = "https://cran.r-project.org/src/contrib/Archive/forcats/forcats_1.0.0.tar.gz",
foreach = "https://cran.r-project.org/src/contrib/foreach_1.5.2.tar.gz",
forecast = "https://cran.r-project.org/src/contrib/forecast_9.0.2.tar.gz",
formatR = "https://cran.r-project.org/src/contrib/formatR_1.14.tar.gz",
Formula = "https://cran.r-project.org/src/contrib/Formula_1.2-5.tar.gz",
fracdiff = "https://cran.r-project.org/src/contrib/fracdiff_1.5-4.tar.gz",
fs = "https://cran.r-project.org/src/contrib/fs_2.1.0.tar.gz",
furrr = "https://cran.r-project.org/src/contrib/furrr_0.4.0.tar.gz",
futile.logger = "https://cran.r-project.org/src/contrib/futile.logger_1.4.9.tar.gz",
futile.options = "https://cran.r-project.org/src/contrib/futile.options_1.0.1.tar.gz",
future = "https://cran.r-project.org/src/contrib/future_1.70.0.tar.gz",
future.apply = "https://cran.r-project.org/src/contrib/future.apply_1.20.2.tar.gz",
gam = "https://cran.r-project.org/src/contrib/gam_1.22-7.tar.gz",
gamlss.dist = "https://cran.r-project.org/src/contrib/gamlss.dist_6.1-1.tar.gz",
gargle = "https://cran.r-project.org/src/contrib/Archive/gargle/gargle_1.5.2.tar.gz",
gclus = "https://cran.r-project.org/src/contrib/gclus_1.3.3.tar.gz",
gdtools = "https://cran.r-project.org/src/contrib/gdtools_0.5.1.tar.gz",
gee = "https://cran.r-project.org/src/contrib/gee_4.13-29.tar.gz",
geepack = "https://cran.r-project.org/src/contrib/geepack_1.3.13.tar.gz",
GeneralizedHyperbolic = "https://cran.r-project.org/src/contrib/GeneralizedHyperbolic_0.8-7.tar.gz",
generics = "https://cran.r-project.org/src/contrib/generics_0.1.4.tar.gz",
GenomeInfoDb = "https://bioconductor.org/packages/3.20/bioc/src/contrib/GenomeInfoDb_1.42.3.tar.gz",
GenomeInfoDbData = "https://bioconductor.org/packages/3.20/data/annotation/src/contrib/GenomeInfoDbData_1.2.13.tar.gz",
GenomicAlignments = "https://bioconductor.org/packages/3.20/bioc/src/contrib/GenomicAlignments_1.42.0.tar.gz",
GenomicFeatures = "https://bioconductor.org/packages/3.20/bioc/src/contrib/GenomicFeatures_1.58.0.tar.gz",
GenomicRanges = "https://bioconductor.org/packages/3.22/bioc/src/contrib/GenomicRanges_1.62.1.tar.gz",
gert = "https://cran.r-project.org/src/contrib/Archive/gert/gert_2.1.5.tar.gz",
getopt = "https://cran.r-project.org/src/contrib/Archive/getopt/getopt_1.20.4.tar.gz",
GetoptLong = "https://cran.r-project.org/src/contrib/Archive/GetoptLong/GetoptLong_1.0.5.tar.gz",
ggalluvial = "https://cran.r-project.org/src/contrib/ggalluvial_0.12.6.tar.gz",
ggbeeswarm = "https://cran.r-project.org/src/contrib/ggbeeswarm_0.7.3.tar.gz",
ggdist = "https://cran.r-project.org/src/contrib/ggdist_3.3.3.tar.gz",
ggforce = "https://cran.r-project.org/src/contrib/ggforce_0.5.0.tar.gz",
ggfun = "https://cran.r-project.org/src/contrib/Archive/ggfun/ggfun_0.2.0.tar.gz",
ggh4x = "https://cran.r-project.org/src/contrib/ggh4x_0.3.1.tar.gz",
ggiraph = "https://cran.r-project.org/src/contrib/Archive/ggiraph/ggiraph_0.9.2.tar.gz",
ggnewscale = "https://cran.r-project.org/src/contrib/ggnewscale_0.5.2.tar.gz",
ggplot2 = "https://cran.r-project.org/src/contrib/ggplot2_4.0.3.tar.gz",
ggplot2movies = "https://cran.r-project.org/src/contrib/ggplot2movies_0.0.1.tar.gz",
ggplotify = "https://cran.r-project.org/src/contrib/ggplotify_0.1.3.tar.gz",
ggpubr = "https://cran.r-project.org/src/contrib/Archive/ggpubr/ggpubr_0.6.2.tar.gz",
ggrastr = "https://cran.r-project.org/src/contrib/ggrastr_1.0.2.tar.gz",
ggrepel = "https://cran.r-project.org/src/contrib/ggrepel_0.9.8.tar.gz",
ggridges = "https://cran.r-project.org/src/contrib/ggridges_0.5.7.tar.gz",
ggsci = "https://cran.r-project.org/src/contrib/Archive/ggsci/ggsci_4.2.0.tar.gz",
ggsignif = "https://cran.r-project.org/src/contrib/ggsignif_0.6.4.tar.gz",
ggtree = "https://github.com/YuLab-SMU/ggtree/archive/5d8e3a43481dbb6f49fcf6806b69993b6aacb0b4.tar.gz",
gh = "https://cran.r-project.org/src/contrib/Archive/gh/gh_1.4.1.tar.gz",
gitcreds = "https://cran.r-project.org/src/contrib/gitcreds_0.1.2.tar.gz",
glmGamPoi = "https://bioconductor.org/packages/3.20/bioc/src/contrib/glmGamPoi_1.18.0.tar.gz",
glmnet = "https://cran.r-project.org/src/contrib/Archive/glmnet/glmnet_4.1-10.tar.gz",
glmnetUtils = "https://cran.r-project.org/src/contrib/glmnetUtils_1.1.9.tar.gz",
GlobalOptions = "https://cran.r-project.org/src/contrib/Archive/GlobalOptions/GlobalOptions_0.1.2.tar.gz",
globals = "https://cran.r-project.org/src/contrib/globals_0.19.1.tar.gz",
glue = "https://cran.r-project.org/src/contrib/glue_1.8.1.tar.gz",
gmm = "https://cran.r-project.org/src/contrib/gmm_1.9-1.tar.gz",
GO.db = "https://bioconductor.org/packages/3.18/data/annotation/src/contrib/GO.db_3.18.0.tar.gz",
goftest = "https://cran.r-project.org/src/contrib/goftest_1.2-3.tar.gz",
googledrive = "https://cran.r-project.org/src/contrib/Archive/googledrive/googledrive_2.1.1.tar.gz",
googlesheets4 = "https://cran.r-project.org/src/contrib/Archive/googlesheets4/googlesheets4_1.1.1.tar.gz",
gower = "https://cran.r-project.org/src/contrib/gower_1.0.2.tar.gz",
GPArotation = "https://cran.r-project.org/src/contrib/GPArotation_2026.6-1.tar.gz",
gplots = "https://cran.r-project.org/src/contrib/gplots_3.3.0.tar.gz",
graph = "https://bioconductor.org/packages/3.18/bioc/src/contrib/graph_1.80.0.tar.gz",
gridExtra = "https://cran.r-project.org/src/contrib/Archive/gridExtra/gridExtra_2.3.tar.gz",
gridGraphics = "https://cran.r-project.org/src/contrib/gridGraphics_0.5-1.tar.gz",
grr = "https://cran.r-project.org/src/contrib/Archive/grr/grr_0.9.5.tar.gz",
GSEABase = "https://bioconductor.org/packages/3.18/bioc/src/contrib/GSEABase_1.64.0.tar.gz",
GSVA = "https://bioconductor.org/packages/3.18/bioc/src/contrib/GSVA_1.50.5.tar.gz",
gtable = "https://cran.r-project.org/src/contrib/gtable_0.3.6.tar.gz",
gtools = "https://cran.r-project.org/src/contrib/gtools_3.9.5.tar.gz",
gypsum = "https://bioconductor.org/packages/3.20/bioc/src/contrib/gypsum_1.2.0.tar.gz",
h5mread = "https://bioconductor.org/packages/3.22/bioc/src/contrib/h5mread_1.2.1.tar.gz",
hardhat = "https://cran.r-project.org/src/contrib/Archive/hardhat/hardhat_1.4.2.tar.gz",
harmony = "https://cran.r-project.org/src/contrib/Archive/harmony/harmony_1.2.3.tar.gz",
haven = "https://cran.r-project.org/src/contrib/Archive/haven/haven_2.5.4.tar.gz",
HDF5Array = "https://bioconductor.org/packages/3.22/bioc/src/contrib/HDF5Array_1.38.0.tar.gz",
hdf5r = "https://cran.r-project.org/src/contrib/hdf5r_1.3.12.tar.gz",
heatmaply = "https://cran.r-project.org/src/contrib/heatmaply_1.6.0.tar.gz",
here = "https://cran.r-project.org/src/contrib/Archive/here/here_1.0.1.tar.gz",
hexbin = "https://cran.r-project.org/src/contrib/hexbin_1.28.5.tar.gz",
highr = "https://cran.r-project.org/src/contrib/highr_0.12.tar.gz",
Hmisc = "https://cran.r-project.org/src/contrib/Hmisc_5.2-6.tar.gz",
hms = "https://cran.r-project.org/src/contrib/hms_1.1.4.tar.gz",
htmlTable = "https://cran.r-project.org/src/contrib/htmlTable_2.5.0.tar.gz",
htmltools = "https://cran.r-project.org/src/contrib/htmltools_0.5.9.tar.gz",
htmlwidgets = "https://cran.r-project.org/src/contrib/htmlwidgets_1.6.4.tar.gz",
httpgd = "https://cran.r-project.org/src/contrib/Archive/httpgd/httpgd_2.0.4.tar.gz",
httpuv = "https://cran.r-project.org/src/contrib/httpuv_1.6.17.tar.gz",
httr = "https://cran.r-project.org/src/contrib/httr_1.4.8.tar.gz",
httr2 = "https://cran.r-project.org/src/contrib/Archive/httr2/httr2_1.1.2.tar.gz",
hunspell = "https://cran.r-project.org/src/contrib/hunspell_3.0.6.tar.gz",
ica = "https://cran.r-project.org/src/contrib/ica_1.0-3.tar.gz",
ids = "https://cran.r-project.org/src/contrib/ids_1.0.1.tar.gz",
igraph = "https://cran.r-project.org/src/contrib/igraph_2.3.3.tar.gz",
ini = "https://cran.r-project.org/src/contrib/ini_0.3.1.tar.gz",
insight = "https://cran.r-project.org/src/contrib/insight_1.5.2.tar.gz",
interp = "https://cran.r-project.org/src/contrib/interp_1.1-6.tar.gz",
ipred = "https://cran.r-project.org/src/contrib/ipred_0.9-15.tar.gz",
IRanges = "https://bioconductor.org/packages/3.22/bioc/src/contrib/IRanges_2.44.0.tar.gz",
irlba = "https://cran.r-project.org/src/contrib/irlba_2.3.7.tar.gz",
isoband = "https://cran.r-project.org/src/contrib/isoband_0.3.0.tar.gz",
iterators = "https://cran.r-project.org/src/contrib/iterators_1.0.14.tar.gz",
JASPAR2020 = "https://bioconductor.org/packages/3.22/data/annotation/src/contrib/JASPAR2020_0.99.10.tar.gz",
joineRML = "https://cran.r-project.org/src/contrib/joineRML_0.4.8.tar.gz",
jose = "https://cran.r-project.org/src/contrib/jose_2.0.0.tar.gz",
jpeg = "https://cran.r-project.org/src/contrib/jpeg_0.1-11.tar.gz",
jquerylib = "https://cran.r-project.org/src/contrib/jquerylib_0.1.4.tar.gz",
jsonlite = "https://cran.r-project.org/src/contrib/jsonlite_2.0.0.tar.gz",
jsonvalidate = "https://cran.r-project.org/src/contrib/Archive/jsonvalidate/jsonvalidate_1.3.2.tar.gz",
KEGGREST = "https://bioconductor.org/packages/3.20/bioc/src/contrib/KEGGREST_1.46.0.tar.gz",
Kendall = "https://cran.r-project.org/src/contrib/Kendall_2.2.2.tar.gz",
kernlab = "https://cran.r-project.org/src/contrib/kernlab_0.9-33.tar.gz",
knitr = "https://cran.r-project.org/src/contrib/knitr_1.51.tar.gz",
ks = "https://cran.r-project.org/src/contrib/ks_1.15.2.tar.gz",
labeling = "https://cran.r-project.org/src/contrib/labeling_0.4.3.tar.gz",
Lahman = "https://cran.r-project.org/src/contrib/Lahman_14.0-0.tar.gz",
lambda.r = "https://cran.r-project.org/src/contrib/lambda.r_1.2.4.tar.gz",
languageserver = "https://cran.r-project.org/src/contrib/Archive/languageserver/languageserver_0.3.16.tar.gz",
later = "https://cran.r-project.org/src/contrib/later_1.4.8.tar.gz",
lattice = "https://cran.r-project.org/src/contrib/Archive/lattice/lattice_0.22-7.tar.gz",
latticeExtra = "https://cran.r-project.org/src/contrib/latticeExtra_0.6-31.tar.gz",
lava = "https://cran.r-project.org/src/contrib/Archive/lava/lava_1.8.1.tar.gz",
lavaan = "https://cran.r-project.org/src/contrib/lavaan_0.7-2.tar.gz",
lazyeval = "https://cran.r-project.org/src/contrib/lazyeval_0.2.3.tar.gz",
leaps = "https://cran.r-project.org/src/contrib/leaps_3.2.tar.gz",
LearnBayes = "https://cran.r-project.org/src/contrib/LearnBayes_2.15.2.tar.gz",
leiden = "https://cran.r-project.org/src/contrib/leiden_0.4.3.1.tar.gz",
leidenbase = "https://cran.r-project.org/src/contrib/leidenbase_0.1.37.tar.gz",
lfe = "https://cran.r-project.org/src/contrib/lfe_3.1.1.tar.gz",
lhs = "https://cran.r-project.org/src/contrib/Archive/lhs/lhs_1.2.0.tar.gz",
lifecycle = "https://cran.r-project.org/src/contrib/lifecycle_1.0.5.tar.gz",
limma = "https://bioconductor.org/packages/3.22/bioc/src/contrib/limma_3.66.0.tar.gz",
lintr = "https://cran.r-project.org/src/contrib/Archive/lintr/lintr_3.2.0.tar.gz",
listenv = "https://cran.r-project.org/src/contrib/Archive/listenv/listenv_0.10.1.tar.gz",
litedown = "https://cran.r-project.org/src/contrib/litedown_0.10.tar.gz",
littler = "https://cran.r-project.org/src/contrib/Archive/littler/littler_0.3.21.tar.gz",
lm.beta = "https://cran.r-project.org/src/contrib/lm.beta_1.7-3.tar.gz",
lme4 = "https://cran.r-project.org/src/contrib/Archive/lme4/lme4_2.0-1.tar.gz",
lmerTest = "https://cran.r-project.org/src/contrib/Archive/lmerTest/lmerTest_3.2-0.tar.gz",
lmodel2 = "https://cran.r-project.org/src/contrib/lmodel2_1.7-4.tar.gz",
lmtest = "https://cran.r-project.org/src/contrib/lmtest_0.9-40.tar.gz",
locfit = "https://cran.r-project.org/src/contrib/locfit_1.5-9.12.tar.gz",
lpSolve = "https://cran.r-project.org/src/contrib/lpSolve_5.6.23.tar.gz",
lpSolveAPI = "https://cran.r-project.org/src/contrib/lpSolveAPI_5.5.2.0-17.15.tar.gz",
lsmeans = "https://cran.r-project.org/src/contrib/lsmeans_2.30-2.tar.gz",
lubridate = "https://cran.r-project.org/src/contrib/Archive/lubridate/lubridate_1.9.4.tar.gz",
magrittr = "https://cran.r-project.org/src/contrib/magrittr_2.0.5.tar.gz",
maps = "https://cran.r-project.org/src/contrib/Archive/maps/maps_3.4.2.1.tar.gz",
marginaleffects = "https://cran.r-project.org/src/contrib/marginaleffects_0.32.0.tar.gz",
margins = "https://cran.r-project.org/src/contrib/margins_0.3.28.tar.gz",
markdown = "https://cran.r-project.org/src/contrib/markdown_2.0.tar.gz",
mathjaxr = "https://cran.r-project.org/src/contrib/mathjaxr_2.0-0.tar.gz",
MatrixGenerics = "https://bioconductor.org/packages/3.22/bioc/src/contrib/MatrixGenerics_1.22.0.tar.gz",
MatrixModels = "https://cran.r-project.org/src/contrib/MatrixModels_0.5-4.tar.gz",
matrixStats = "https://cran.r-project.org/src/contrib/matrixStats_1.5.0.tar.gz",
maxLik = "https://cran.r-project.org/src/contrib/maxLik_1.5-2.2.tar.gz",
mc2d = "https://cran.r-project.org/src/contrib/mc2d_0.2.2.tar.gz",
mclust = "https://cran.r-project.org/src/contrib/mclust_6.1.3.tar.gz",
mcmc = "https://cran.r-project.org/src/contrib/mcmc_0.9-8.tar.gz",
MCMCpack = "https://cran.r-project.org/src/contrib/MCMCpack_1.7-1.tar.gz",
mediation = "https://cran.r-project.org/src/contrib/mediation_4.5.1.tar.gz",
memoise = "https://cran.r-project.org/src/contrib/memoise_2.0.1.tar.gz",
metadat = "https://cran.r-project.org/src/contrib/metadat_1.6-0.tar.gz",
metafor = "https://cran.r-project.org/src/contrib/metafor_5.0-1.tar.gz",
metapod = "https://bioconductor.org/packages/3.20/bioc/src/contrib/metapod_1.14.0.tar.gz",
mfx = "https://cran.r-project.org/src/contrib/mfx_1.2-4.tar.gz",
mgcv = "https://cran.r-project.org/src/contrib/Archive/mgcv/mgcv_1.9-3.tar.gz",
microbenchmark = "https://cran.r-project.org/src/contrib/microbenchmark_1.5.0.tar.gz",
micsr = "https://cran.r-project.org/src/contrib/micsr_0.1-5.tar.gz",
mime = "https://cran.r-project.org/src/contrib/mime_0.13.tar.gz",
minioclient = "https://cran.r-project.org/src/contrib/minioclient_0.0.6.tar.gz",
miniUI = "https://cran.r-project.org/src/contrib/miniUI_0.1.2.tar.gz",
minqa = "https://cran.r-project.org/src/contrib/minqa_1.2.8.tar.gz",
miscTools = "https://cran.r-project.org/src/contrib/miscTools_0.6-30.tar.gz",
mitools = "https://cran.r-project.org/src/contrib/mitools_2.4.tar.gz",
mixtools = "https://cran.r-project.org/src/contrib/mixtools_2.0.0.1.tar.gz",
mlogit = "https://cran.r-project.org/src/contrib/mlogit_2.0-0.tar.gz",
mnormt = "https://cran.r-project.org/src/contrib/mnormt_2.1.2.tar.gz",
mockery = "https://cran.r-project.org/src/contrib/mockery_0.4.5.tar.gz",
modeldata = "https://cran.r-project.org/src/contrib/modeldata_1.5.1.tar.gz",
ModelMetrics = "https://cran.r-project.org/src/contrib/ModelMetrics_1.2.2.2.tar.gz",
modelr = "https://cran.r-project.org/src/contrib/modelr_0.1.11.tar.gz",
modeltests = "https://cran.r-project.org/src/contrib/modeltests_0.1.8.tar.gz",
modeltools = "https://cran.r-project.org/src/contrib/modeltools_0.2-24.tar.gz",
monocle3 = "https://github.com/cole-trapnell-lab/monocle3/archive/536f1033d6de7c957f26a1f403f81efbd825e0db.tar.gz",
msigdbr = "https://cran.r-project.org/src/contrib/msigdbr_26.1.0.tar.gz",
muhaz = "https://cran.r-project.org/src/contrib/muhaz_1.2.6.4.tar.gz",
multcomp = "https://cran.r-project.org/src/contrib/multcomp_1.4-31.tar.gz",
multicool = "https://cran.r-project.org/src/contrib/multicool_1.0.1.tar.gz",
munsell = "https://cran.r-project.org/src/contrib/munsell_0.5.1.tar.gz",
mvtnorm = "https://cran.r-project.org/src/contrib/mvtnorm_1.4-2.tar.gz",
nanoarrow = "https://cran.r-project.org/src/contrib/nanoarrow_0.8.0-1.tar.gz",
nanonext = "https://cran.r-project.org/src/contrib/nanonext_1.10.1.tar.gz",
network = "https://cran.r-project.org/src/contrib/network_1.20.0.tar.gz",
nloptr = "https://cran.r-project.org/src/contrib/nloptr_2.2.1.tar.gz",
ntfy = "https://cran.r-project.org/src/contrib/ntfy_0.1.0.tar.gz",
numDeriv = "https://cran.r-project.org/src/contrib/numDeriv_2016.8-1.1.tar.gz",
officer = "https://cran.r-project.org/src/contrib/officer_0.7.6.tar.gz",
openssl = "https://cran.r-project.org/src/contrib/openssl_2.4.2.tar.gz",
optparse = "https://cran.r-project.org/src/contrib/Archive/optparse/optparse_1.7.5.tar.gz",
ordinal = "https://cran.r-project.org/src/contrib/ordinal_2025.12-29.tar.gz",
org.Hs.eg.db = "https://bioconductor.org/packages/3.20/data/annotation/src/contrib/org.Hs.eg.db_3.20.0.tar.gz",
otel = "https://cran.r-project.org/src/contrib/otel_0.2.0.tar.gz",
packrat = "https://cran.r-project.org/src/contrib/packrat_0.9.3.tar.gz",
pacman = "https://cran.r-project.org/src/contrib/pacman_0.5.1.tar.gz",
pak = "https://cran.r-project.org/src/contrib/Archive/pak/pak_0.9.5.tar.gz",
palmerpenguins = "https://cran.r-project.org/src/contrib/palmerpenguins_0.1.1.tar.gz",
pander = "https://cran.r-project.org/src/contrib/pander_0.6.6.tar.gz",
parallelDist = "https://cran.r-project.org/src/contrib/Archive/parallelDist/parallelDist_0.2.6.tar.gz",
parallelly = "https://cran.r-project.org/src/contrib/Archive/parallelly/parallelly_1.47.0.tar.gz",
PASWR = "https://cran.r-project.org/src/contrib/PASWR_1.3.tar.gz",
patchwork = "https://cran.r-project.org/src/contrib/patchwork_1.3.2.tar.gz",
patrick = "https://cran.r-project.org/src/contrib/patrick_0.3.1.tar.gz",
paws.common = "https://cran.r-project.org/src/contrib/paws.common_0.8.10.tar.gz",
pbapply = "https://cran.r-project.org/src/contrib/pbapply_1.7-4.tar.gz",
pbivnorm = "https://cran.r-project.org/src/contrib/pbivnorm_0.6.0.tar.gz",
pbkrtest = "https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.5.4.tar.gz",
pbmcapply = "https://cran.r-project.org/src/contrib/pbmcapply_1.5.1.tar.gz",
pcaPP = "https://cran.r-project.org/src/contrib/pcaPP_2.0-5.tar.gz",
permute = "https://cran.r-project.org/src/contrib/permute_0.9-10.tar.gz",
pheatmap = "https://cran.r-project.org/src/contrib/pheatmap_1.0.13.tar.gz",
pillar = "https://cran.r-project.org/src/contrib/pillar_1.11.1.tar.gz",
pkgbuild = "https://cran.r-project.org/src/contrib/pkgbuild_1.4.8.tar.gz",
pkgconfig = "https://cran.r-project.org/src/contrib/pkgconfig_2.0.3.tar.gz",
pkgdown = "https://cran.r-project.org/src/contrib/Archive/pkgdown/pkgdown_2.2.0.tar.gz",
pkgKitten = "https://cran.r-project.org/src/contrib/pkgKitten_0.2.4.tar.gz",
pkgload = "https://cran.r-project.org/src/contrib/Archive/pkgload/pkgload_1.5.2.tar.gz",
plm = "https://cran.r-project.org/src/contrib/plm_2.6-7.tar.gz",
plogr = "https://cran.r-project.org/src/contrib/Archive/plogr/plogr_0.2.0.tar.gz",
plotly = "https://cran.r-project.org/src/contrib/plotly_4.12.0.tar.gz",
plotrix = "https://cran.r-project.org/src/contrib/plotrix_3.8-14.tar.gz",
plyr = "https://cran.r-project.org/src/contrib/plyr_1.8.9.tar.gz",
png = "https://cran.r-project.org/src/contrib/png_0.1-9.tar.gz",
poLCA = "https://cran.r-project.org/src/contrib/poLCA_1.6.0.2.tar.gz",
polyclip = "https://cran.r-project.org/src/contrib/polyclip_1.10-7.tar.gz",
polynom = "https://cran.r-project.org/src/contrib/polynom_1.4-1.tar.gz",
poweRlaw = "https://cran.r-project.org/src/contrib/poweRlaw_1.0.0.tar.gz",
pracma = "https://cran.r-project.org/src/contrib/Archive/pracma/pracma_2.4.4.tar.gz",
praise = "https://cran.r-project.org/src/contrib/praise_1.0.0.tar.gz",
prediction = "https://cran.r-project.org/src/contrib/prediction_0.3.18.tar.gz",
presto = "https://satijalab.r-universe.dev/src/contrib/presto_1.0.0.tar.gz",
prettycode = "https://cran.r-project.org/src/contrib/prettycode_1.1.0.tar.gz",
prettydoc = "https://cran.r-project.org/src/contrib/prettydoc_0.4.1.tar.gz",
prettyGraphs = "https://cran.r-project.org/src/contrib/prettyGraphs_2.2.0.tar.gz",
prettyunits = "https://cran.r-project.org/src/contrib/prettyunits_1.2.0.tar.gz",
pROC = "https://cran.r-project.org/src/contrib/pROC_1.19.0.1.tar.gz",
processx = "https://cran.r-project.org/src/contrib/processx_3.9.0.tar.gz",
prodlim = "https://cran.r-project.org/src/contrib/Archive/prodlim/prodlim_2025.04.28.tar.gz",
profmem = "https://cran.r-project.org/src/contrib/profmem_0.7.0.tar.gz",
profvis = "https://cran.r-project.org/src/contrib/profvis_0.4.0.tar.gz",
progress = "https://cran.r-project.org/src/contrib/progress_1.2.3.tar.gz",
progressr = "https://cran.r-project.org/src/contrib/Archive/progressr/progressr_0.19.0.tar.gz",
promises = "https://cran.r-project.org/src/contrib/promises_1.5.0.tar.gz",
ProtGenerics = "https://bioconductor.org/packages/3.20/bioc/src/contrib/ProtGenerics_1.38.0.tar.gz",
proto = "https://cran.r-project.org/src/contrib/proto_1.0.0.tar.gz",
proxy = "https://cran.r-project.org/src/contrib/proxy_0.4-29.tar.gz",
ps = "https://cran.r-project.org/src/contrib/ps_1.9.3.tar.gz",
pscl = "https://cran.r-project.org/src/contrib/pscl_1.5.9.tar.gz",
psych = "https://cran.r-project.org/src/contrib/psych_2.6.5.tar.gz",
purrr = "https://cran.r-project.org/src/contrib/purrr_1.2.2.tar.gz",
pwalign = "https://bioconductor.org/packages/3.20/bioc/src/contrib/pwalign_1.2.0.tar.gz",
qap = "https://cran.r-project.org/src/contrib/qap_0.1-2.tar.gz",
quadprog = "https://cran.r-project.org/src/contrib/quadprog_1.5-8.tar.gz",
quantmod = "https://cran.r-project.org/src/contrib/quantmod_0.4.29.tar.gz",
quantreg = "https://cran.r-project.org/src/contrib/quantreg_6.1.tar.gz",
quarto = "https://cran.r-project.org/src/contrib/Archive/quarto/quarto_1.4.4.tar.gz",
R.cache = "https://cran.r-project.org/src/contrib/R.cache_0.17.0.tar.gz",
R.methodsS3 = "https://cran.r-project.org/src/contrib/R.methodsS3_1.8.2.tar.gz",
R.oo = "https://cran.r-project.org/src/contrib/R.oo_1.27.1.tar.gz",
R.rsp = "https://cran.r-project.org/src/contrib/R.rsp_0.46.0.tar.gz",
R.utils = "https://cran.r-project.org/src/contrib/R.utils_2.13.0.tar.gz",
r2d2 = "https://cran.r-project.org/src/contrib/r2d2_1.0.2.tar.gz",
R6 = "https://cran.r-project.org/src/contrib/R6_2.6.1.tar.gz",
ragg = "https://cran.r-project.org/src/contrib/ragg_1.5.2.tar.gz",
randomForest = "https://cran.r-project.org/src/contrib/randomForest_4.7-1.2.tar.gz",
randtoolbox = "https://cran.r-project.org/src/contrib/randtoolbox_2.0.5.tar.gz",
RANN = "https://cran.r-project.org/src/contrib/RANN_2.6.2.tar.gz",
rappdirs = "https://cran.r-project.org/src/contrib/rappdirs_0.3.4.tar.gz",
rbibutils = "https://cran.r-project.org/src/contrib/rbibutils_2.4.1.tar.gz",
rcmdcheck = "https://cran.r-project.org/src/contrib/rcmdcheck_1.4.0.tar.gz",
RColorBrewer = "https://cran.r-project.org/src/contrib/RColorBrewer_1.1-3.tar.gz",
Rcpp = "https://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_1.1.1-1.1.tar.gz",
RcppAnnoy = "https://cran.r-project.org/src/contrib/RcppAnnoy_0.0.23.tar.gz",
RcppArmadillo = "https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_15.2.6-1.tar.gz",
RcppEigen = "https://cran.r-project.org/src/contrib/RcppEigen_0.3.4.0.2.tar.gz",
RcppHNSW = "https://cran.r-project.org/src/contrib/RcppHNSW_0.7.0.tar.gz",
RcppML = "https://cran.r-project.org/src/contrib/Archive/RcppML/RcppML_0.3.7.tar.gz",
RcppParallel = "https://cran.r-project.org/src/contrib/Archive/RcppParallel/RcppParallel_5.1.10.tar.gz",
RcppProgress = "https://cran.r-project.org/src/contrib/RcppProgress_0.4.2.tar.gz",
RcppRoll = "https://cran.r-project.org/src/contrib/Archive/RcppRoll/RcppRoll_0.3.1.tar.gz",
RcppTOML = "https://cran.r-project.org/src/contrib/RcppTOML_0.2.3.tar.gz",
RCurl = "https://cran.r-project.org/src/contrib/RCurl_1.98-1.19.tar.gz",
Rdpack = "https://cran.r-project.org/src/contrib/Rdpack_2.6.6.tar.gz",
readr = "https://cran.r-project.org/src/contrib/readr_2.2.0.tar.gz",
readxl = "https://cran.r-project.org/src/contrib/Archive/readxl/readxl_1.4.5.tar.gz",
recipes = "https://cran.r-project.org/src/contrib/Archive/recipes/recipes_1.3.1.tar.gz",
reformulas = "https://cran.r-project.org/src/contrib/reformulas_0.4.4.tar.gz",
registry = "https://cran.r-project.org/src/contrib/registry_0.5-1.tar.gz",
rematch = "https://cran.r-project.org/src/contrib/rematch_2.0.0.tar.gz",
rematch2 = "https://cran.r-project.org/src/contrib/rematch2_2.1.2.tar.gz",
remotes = "https://cran.r-project.org/src/contrib/remotes_2.5.0.tar.gz",
renv = "https://cran.r-project.org/src/contrib/renv_1.2.3.tar.gz",
reprex = "https://cran.r-project.org/src/contrib/reprex_2.1.1.tar.gz",
reshape2 = "https://cran.r-project.org/src/contrib/reshape2_1.4.5.tar.gz",
ResidualMatrix = "https://bioconductor.org/packages/3.22/bioc/src/contrib/ResidualMatrix_1.20.0.tar.gz",
restfulr = "https://cran.r-project.org/src/contrib/Archive/restfulr/restfulr_0.0.16.tar.gz",
reticulate = "https://cran.r-project.org/src/contrib/Archive/reticulate/reticulate_1.42.0.tar.gz",
rex = "https://cran.r-project.org/src/contrib/Archive/rex/rex_1.2.1.tar.gz",
rgenoud = "https://cran.r-project.org/src/contrib/rgenoud_5.9-0.11.tar.gz",
rhdf5 = "https://bioconductor.org/packages/3.22/bioc/src/contrib/rhdf5_2.54.1.tar.gz",
rhdf5filters = "https://bioconductor.org/packages/3.22/bioc/src/contrib/rhdf5filters_1.22.0.tar.gz",
Rhdf5lib = "https://bioconductor.org/packages/3.22/bioc/src/contrib/Rhdf5lib_1.32.0.tar.gz",
RhpcBLASctl = "https://cran.r-project.org/src/contrib/RhpcBLASctl_0.23-42.tar.gz",
Rhtslib = "https://bioconductor.org/packages/3.20/bioc/src/contrib/Rhtslib_3.2.0.tar.gz",
rjson = "https://cran.r-project.org/src/contrib/rjson_0.2.23.tar.gz",
rlang = "https://cran.r-project.org/src/contrib/Archive/rlang/rlang_1.2.0.tar.gz",
rle = "https://cran.r-project.org/src/contrib/rle_0.10.0.tar.gz",
rmarkdown = "https://cran.r-project.org/src/contrib/rmarkdown_2.31.tar.gz",
rngWELL = "https://cran.r-project.org/src/contrib/rngWELL_0.10-10.tar.gz",
rnndescent = "https://cran.r-project.org/src/contrib/Archive/rnndescent/rnndescent_0.1.8.tar.gz",
robust = "https://cran.r-project.org/src/contrib/robust_0.7-5.tar.gz",
robustbase = "https://cran.r-project.org/src/contrib/robustbase_0.99-7.tar.gz",
ROCR = "https://cran.r-project.org/src/contrib/ROCR_1.0-12.tar.gz",
roxygen2 = "https://cran.r-project.org/src/contrib/roxygen2_8.0.0.tar.gz",
rprojroot = "https://cran.r-project.org/src/contrib/rprojroot_2.1.1.tar.gz",
RPushbullet = "https://cran.r-project.org/src/contrib/RPushbullet_0.3.5.tar.gz",
rrcov = "https://cran.r-project.org/src/contrib/rrcov_1.7-7.tar.gz",
rsample = "https://cran.r-project.org/src/contrib/rsample_1.3.2.tar.gz",
Rsamtools = "https://bioconductor.org/packages/3.20/bioc/src/contrib/Rsamtools_2.22.0.tar.gz",
rsconnect = "https://cran.r-project.org/src/contrib/rsconnect_1.10.1.tar.gz",
RSpectra = "https://cran.r-project.org/src/contrib/RSpectra_0.16-2.tar.gz",
RSQLite = "https://cran.r-project.org/src/contrib/Archive/RSQLite/RSQLite_2.3.9.tar.gz",
rstatix = "https://cran.r-project.org/src/contrib/Archive/rstatix/rstatix_0.7.2.tar.gz",
rstudioapi = "https://cran.r-project.org/src/contrib/rstudioapi_0.19.0.tar.gz",
rsvd = "https://cran.r-project.org/src/contrib/rsvd_1.0.5.tar.gz",
rtracklayer = "https://bioconductor.org/packages/3.20/bioc/src/contrib/rtracklayer_1.66.0.tar.gz",
Rtsne = "https://cran.r-project.org/src/contrib/Rtsne_0.17.tar.gz",
rversions = "https://cran.r-project.org/src/contrib/rversions_3.0.0.tar.gz",
rvest = "https://cran.r-project.org/src/contrib/Archive/rvest/rvest_1.0.4.tar.gz",
s2 = "https://cran.r-project.org/src/contrib/s2_1.1.11.tar.gz",
S4Arrays = "https://bioconductor.org/packages/3.22/bioc/src/contrib/S4Arrays_1.10.1.tar.gz",
S4Vectors = "https://bioconductor.org/packages/3.22/bioc/src/contrib/S4Vectors_0.48.1.tar.gz",
S7 = "https://cran.r-project.org/src/contrib/S7_0.2.2.tar.gz",
sandwich = "https://cran.r-project.org/src/contrib/sandwich_3.1-2.tar.gz",
SASmixed = "https://cran.r-project.org/src/contrib/SASmixed_1.0-5.tar.gz",
sass = "https://cran.r-project.org/src/contrib/sass_0.4.10.tar.gz",
ScaledMatrix = "https://bioconductor.org/packages/3.22/bioc/src/contrib/ScaledMatrix_1.18.0.tar.gz",
scales = "https://cran.r-project.org/src/contrib/scales_1.4.0.tar.gz",
scater = "https://bioconductor.org/packages/3.20/bioc/src/contrib/scater_1.34.1.tar.gz",
scattermore = "https://cran.r-project.org/src/contrib/scattermore_1.2.tar.gz",
scatterplot3d = "https://cran.r-project.org/src/contrib/scatterplot3d_0.3-45.tar.gz",
scDblFinder = "https://bioconductor.org/packages/3.21/bioc/src/contrib/scDblFinder_1.22.0.tar.gz",
scGSVA = "https://github.com/guokai8/scGSVA/archive/1bd086f1ef20a3376a5513189192a53e3ae4c9cb.tar.gz",
scran = "https://bioconductor.org/packages/3.20/bioc/src/contrib/scran_1.34.0.tar.gz",
sctransform = "https://cran.r-project.org/src/contrib/sctransform_0.4.3.tar.gz",
scuttle = "https://bioconductor.org/packages/3.22/bioc/src/contrib/scuttle_1.20.0.tar.gz",
segmented = "https://cran.r-project.org/src/contrib/Archive/segmented/segmented_2.1-4.tar.gz",
selectr = "https://cran.r-project.org/src/contrib/Archive/selectr/selectr_0.4-2.tar.gz",
Seqinfo = "https://bioconductor.org/packages/3.22/bioc/src/contrib/Seqinfo_1.0.0.tar.gz",
seqLogo = "https://bioconductor.org/packages/3.20/bioc/src/contrib/seqLogo_1.72.0.tar.gz",
seriation = "https://cran.r-project.org/src/contrib/seriation_1.5.8.tar.gz",
sessioninfo = "https://cran.r-project.org/src/contrib/Archive/sessioninfo/sessioninfo_1.2.3.tar.gz",
sets = "https://cran.r-project.org/src/contrib/sets_1.0-25.tar.gz",
Seurat = "https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_4.4.0.tar.gz",
SeuratData = "https://github.com/satijalab/seurat-data/archive/3e51f44303069b64f5dc4d68e6a3d4a343f55c39.tar.gz",
SeuratObject = "https://cran.r-project.org/src/contrib/Archive/SeuratObject/SeuratObject_4.1.4.tar.gz",
SeuratWrappers = "https://github.com/satijalab/seurat-wrappers/archive/a1eb0d8b039ad6d5ef0ff1332fd8eb1c0c223553.tar.gz",
sf = "https://cran.r-project.org/src/contrib/sf_1.1-1.tar.gz",
shadowtext = "https://cran.r-project.org/src/contrib/shadowtext_0.1.6.tar.gz",
shape = "https://cran.r-project.org/src/contrib/shape_1.4.6.1.tar.gz",
shiny = "https://cran.r-project.org/src/contrib/Archive/shiny/shiny_1.13.0.tar.gz",
shinyBS = "https://cran.r-project.org/src/contrib/Archive/shinyBS/shinyBS_0.61.1.tar.gz",
shinydashboard = "https://cran.r-project.org/src/contrib/shinydashboard_0.7.3.tar.gz",
shinyjs = "https://cran.r-project.org/src/contrib/Archive/shinyjs/shinyjs_2.1.0.tar.gz",
showtext = "https://cran.r-project.org/src/contrib/Archive/showtext/showtext_0.9-7.tar.gz",
showtextdb = "https://cran.r-project.org/src/contrib/showtextdb_3.0.tar.gz",
Signac = "https://cran.r-project.org/src/contrib/Archive/Signac/Signac_1.15.0.tar.gz",
simplermarkdown = "https://cran.r-project.org/src/contrib/simplermarkdown_0.0.6.tar.gz",
SingleCellExperiment = "https://bioconductor.org/packages/3.22/bioc/src/contrib/SingleCellExperiment_1.32.0.tar.gz",
SingleR = "https://bioconductor.org/packages/3.20/bioc/src/contrib/SingleR_2.8.0.tar.gz",
sitmo = "https://cran.r-project.org/src/contrib/sitmo_2.0.2.tar.gz",
slam = "https://cran.r-project.org/src/contrib/Archive/slam/slam_0.1-55.tar.gz",
slider = "https://cran.r-project.org/src/contrib/slider_0.3.3.tar.gz",
sm = "https://cran.r-project.org/src/contrib/sm_2.2-6.0.tar.gz",
sna = "https://cran.r-project.org/src/contrib/sna_2.8.tar.gz",
SNFtool = "https://cran.r-project.org/src/contrib/SNFtool_2.3.1.tar.gz",
snow = "https://cran.r-project.org/src/contrib/snow_0.4-4.tar.gz",
snowflakeauth = "https://cran.r-project.org/src/contrib/snowflakeauth_0.2.2.tar.gz",
sourcetools = "https://cran.r-project.org/src/contrib/sourcetools_0.1.7-2.tar.gz",
sp = "https://cran.r-project.org/src/contrib/Archive/sp/sp_2.2-1.tar.gz",
spam = "https://cran.r-project.org/src/contrib/Archive/spam/spam_2.11-3.tar.gz",
spam64 = "https://cran.r-project.org/src/contrib/spam64_2.11-4.tar.gz",
SparseArray = "https://bioconductor.org/packages/3.22/bioc/src/contrib/SparseArray_1.10.10.tar.gz",
SparseM = "https://cran.r-project.org/src/contrib/SparseM_1.84-2.tar.gz",
sparseMatrixStats = "https://bioconductor.org/packages/3.22/bioc/src/contrib/sparseMatrixStats_1.22.0.tar.gz",
sparsevctrs = "https://cran.r-project.org/src/contrib/Archive/sparsevctrs/sparsevctrs_0.3.4.tar.gz",
spatstat = "https://cran.r-project.org/src/contrib/spatstat_3.6-1.tar.gz",
spatstat.data = "https://cran.r-project.org/src/contrib/spatstat.data_3.1-9.tar.gz",
spatstat.explore = "https://cran.r-project.org/src/contrib/Archive/spatstat.explore/spatstat.explore_3.8-0.tar.gz",
spatstat.geom = "https://cran.r-project.org/src/contrib/Archive/spatstat.geom/spatstat.geom_3.7-3.tar.gz",
spatstat.linnet = "https://cran.r-project.org/src/contrib/spatstat.linnet_3.5-1.tar.gz",
spatstat.model = "https://cran.r-project.org/src/contrib/spatstat.model_3.7-1.tar.gz",
spatstat.random = "https://cran.r-project.org/src/contrib/Archive/spatstat.random/spatstat.random_3.4-5.tar.gz",
spatstat.sparse = "https://cran.r-project.org/src/contrib/spatstat.sparse_3.2-0.tar.gz",
spatstat.univar = "https://cran.r-project.org/src/contrib/Archive/spatstat.univar/spatstat.univar_3.1-7.tar.gz",
spatstat.utils = "https://cran.r-project.org/src/contrib/Archive/spatstat.utils/spatstat.utils_3.2-2.tar.gz",
spData = "https://cran.r-project.org/src/contrib/spData_2.3.5.tar.gz",
spdep = "https://cran.r-project.org/src/contrib/spdep_1.4-2.tar.gz",
speedglm = "https://github.com/cole-trapnell-lab/speedglm/archive/ca34b4e53319424b60c442bb550adf3574b4bfec.tar.gz",
spelling = "https://cran.r-project.org/src/contrib/spelling_2.3.2.tar.gz",
SQUAREM = "https://cran.r-project.org/src/contrib/Archive/SQUAREM/SQUAREM_2021.1.tar.gz",
statmod = "https://cran.r-project.org/src/contrib/statmod_1.5.2.tar.gz",
statnet.common = "https://cran.r-project.org/src/contrib/statnet.common_4.13.0.tar.gz",
stringi = "https://cran.r-project.org/src/contrib/stringi_1.8.7.tar.gz",
stringmagic = "https://cran.r-project.org/src/contrib/stringmagic_1.2.0.tar.gz",
stringr = "https://cran.r-project.org/src/contrib/stringr_1.6.0.tar.gz",
strucchange = "https://cran.r-project.org/src/contrib/strucchange_1.5-4.tar.gz",
styler = "https://cran.r-project.org/src/contrib/Archive/styler/styler_1.10.3.tar.gz",
SummarizedExperiment = "https://bioconductor.org/packages/3.22/bioc/src/contrib/SummarizedExperiment_1.40.0.tar.gz",
survey = "https://cran.r-project.org/src/contrib/survey_4.5.tar.gz",
sys = "https://cran.r-project.org/src/contrib/sys_3.4.3.tar.gz",
sysfonts = "https://cran.r-project.org/src/contrib/sysfonts_0.8.9.tar.gz",
systemfit = "https://cran.r-project.org/src/contrib/systemfit_1.1-30.tar.gz",
systemfonts = "https://cran.r-project.org/src/contrib/systemfonts_1.3.2.tar.gz",
tensor = "https://cran.r-project.org/src/contrib/tensor_1.5.1.tar.gz",
testthat = "https://cran.r-project.org/src/contrib/testthat_3.3.2.tar.gz",
textshaping = "https://cran.r-project.org/src/contrib/textshaping_1.0.5.tar.gz",
TFMPvalue = "https://cran.r-project.org/src/contrib/Archive/TFMPvalue/TFMPvalue_0.0.9.tar.gz",
TH.data = "https://cran.r-project.org/src/contrib/TH.data_1.1-5.tar.gz",
tibble = "https://cran.r-project.org/src/contrib/tibble_3.3.1.tar.gz",
tidyr = "https://cran.r-project.org/src/contrib/tidyr_1.3.2.tar.gz",
tidyselect = "https://cran.r-project.org/src/contrib/tidyselect_1.2.1.tar.gz",
tidytree = "https://cran.r-project.org/src/contrib/Archive/tidytree/tidytree_0.4.7.tar.gz",
tidyverse = "https://cran.r-project.org/src/contrib/tidyverse_2.0.0.tar.gz",
timechange = "https://cran.r-project.org/src/contrib/Archive/timechange/timechange_0.3.0.tar.gz",
timeDate = "https://cran.r-project.org/src/contrib/Archive/timeDate/timeDate_4041.110.tar.gz",
tinytest = "https://cran.r-project.org/src/contrib/tinytest_1.4.3.tar.gz",
tinytex = "https://cran.r-project.org/src/contrib/Archive/tinytex/tinytex_0.59.tar.gz",
TMB = "https://cran.r-project.org/src/contrib/Archive/TMB/TMB_1.9.19.tar.gz",
transport = "https://cran.r-project.org/src/contrib/transport_0.15-4.tar.gz",
treeio = "https://bioconductor.org/packages/3.20/bioc/src/contrib/treeio_1.30.0.tar.gz",
truncdist = "https://cran.r-project.org/src/contrib/truncdist_1.0-2.tar.gz",
trust = "https://cran.r-project.org/src/contrib/trust_0.1-9.tar.gz",
tseries = "https://cran.r-project.org/src/contrib/tseries_0.10-62.tar.gz",
TSP = "https://cran.r-project.org/src/contrib/TSP_1.2.7.tar.gz",
TTR = "https://cran.r-project.org/src/contrib/TTR_0.24.4.tar.gz",
tufte = "https://cran.r-project.org/src/contrib/tufte_0.15.0.tar.gz",
tweenr = "https://cran.r-project.org/src/contrib/tweenr_2.0.3.tar.gz",
tzdb = "https://cran.r-project.org/src/contrib/tzdb_0.5.0.tar.gz",
ucminf = "https://cran.r-project.org/src/contrib/ucminf_1.2.3.tar.gz",
UCSC.utils = "https://bioconductor.org/packages/3.20/bioc/src/contrib/UCSC.utils_1.2.0.tar.gz",
UCSCXenaTools = "https://cran.r-project.org/src/contrib/Archive/UCSCXenaTools/UCSCXenaTools_1.6.1.tar.gz",
umap = "https://cran.r-project.org/src/contrib/umap_0.2.10.0.tar.gz",
unigd = "https://cran.r-project.org/src/contrib/Archive/unigd/unigd_0.1.3.tar.gz",
units = "https://cran.r-project.org/src/contrib/units_1.0-1.tar.gz",
urca = "https://cran.r-project.org/src/contrib/urca_1.3-4.tar.gz",
urlchecker = "https://cran.r-project.org/src/contrib/Archive/urlchecker/urlchecker_1.0.1.tar.gz",
usethis = "https://cran.r-project.org/src/contrib/usethis_3.2.1.tar.gz",
utf8 = "https://cran.r-project.org/src/contrib/utf8_1.2.6.tar.gz",
uuid = "https://cran.r-project.org/src/contrib/Archive/uuid/uuid_1.2-1.tar.gz",
uwot = "https://cran.r-project.org/src/contrib/uwot_0.2.4.tar.gz",
V8 = "https://cran.r-project.org/src/contrib/Archive/V8/V8_6.0.3.tar.gz",
vars = "https://cran.r-project.org/src/contrib/vars_1.6-1.tar.gz",
vctrs = "https://cran.r-project.org/src/contrib/vctrs_0.7.3.tar.gz",
vdiffr = "https://cran.r-project.org/src/contrib/vdiffr_1.0.9.tar.gz",
vegan = "https://cran.r-project.org/src/contrib/vegan_2.7-5.tar.gz",
VennDiagram = "https://cran.r-project.org/src/contrib/VennDiagram_1.8.2.tar.gz",
vipor = "https://cran.r-project.org/src/contrib/vipor_0.4.7.tar.gz",
viridis = "https://cran.r-project.org/src/contrib/viridis_0.6.5.tar.gz",
viridisLite = "https://cran.r-project.org/src/contrib/viridisLite_0.4.3.tar.gz",
vroom = "https://cran.r-project.org/src/contrib/vroom_1.7.1.tar.gz",
waldo = "https://cran.r-project.org/src/contrib/waldo_0.6.2.tar.gz",
warp = "https://cran.r-project.org/src/contrib/warp_0.2.3.tar.gz",
webfakes = "https://cran.r-project.org/src/contrib/webfakes_1.5.0.tar.gz",
webshot = "https://cran.r-project.org/src/contrib/webshot_0.5.5.tar.gz",
whisker = "https://cran.r-project.org/src/contrib/whisker_0.4.1.tar.gz",
whoami = "https://cran.r-project.org/src/contrib/whoami_1.3.0.tar.gz",
withr = "https://cran.r-project.org/src/contrib/Archive/withr/withr_3.0.2.tar.gz",
wk = "https://cran.r-project.org/src/contrib/wk_0.9.5.tar.gz",
xfun = "https://cran.r-project.org/src/contrib/Archive/xfun/xfun_0.58.tar.gz",
xgboost = "https://cran.r-project.org/src/contrib/xgboost_3.2.1.1.tar.gz",
XML = "https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.18.tar.gz",
xml2 = "https://cran.r-project.org/src/contrib/Archive/xml2/xml2_1.3.8.tar.gz",
xmlparsedata = "https://cran.r-project.org/src/contrib/xmlparsedata_1.0.5.tar.gz",
xopen = "https://cran.r-project.org/src/contrib/xopen_1.0.1.tar.gz",
xtable = "https://cran.r-project.org/src/contrib/xtable_1.8-8.tar.gz",
xts = "https://cran.r-project.org/src/contrib/xts_0.14.2.tar.gz",
XVector = "https://bioconductor.org/packages/3.22/bioc/src/contrib/XVector_0.50.0.tar.gz",
yaml = "https://cran.r-project.org/src/contrib/yaml_2.3.12.tar.gz",
yulab.utils = "https://cran.r-project.org/src/contrib/Archive/yulab.utils/yulab.utils_0.2.3.tar.gz",
zip = "https://cran.r-project.org/src/contrib/Archive/zip/zip_2.3.3.tar.gz",
zlibbioc = "https://bioconductor.org/packages/3.20/bioc/src/contrib/zlibbioc_1.52.0.tar.gz",
zoo = "https://cran.r-project.org/src/contrib/zoo_1.8-15.tar.gz"
)
source_url <- unname(source_urls[r_versions$package])
source_url[is.na(source_url)] <- ""
lock <- data.frame(
  Package = r_versions$package,
  SelectedVersion = r_versions$version,
  SourceCategory = ifelse(nzchar(source_url), "source", "base"),
  SourceURL = source_url,
  DockerInstallRef = ifelse(nzchar(source_url), paste0("url::", source_url), ""),
  stringsAsFactors = FALSE
)
dir.create(install_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(install_library, .libPaths()))

installed_in_library <- function(library) {
  result <- installed.packages(lib.loc = library, noCache = TRUE)
  if (!nrow(result)) {
    return(data.frame(Package = character(), Version = character()))
  }
  data.frame(
    Package = result[, "Package"],
    Version = result[, "Version"],
    stringsAsFactors = FALSE
  )
}

visible_installed <- function() {
  inventories <- lapply(.libPaths(), function(library) {
    result <- installed.packages(lib.loc = library, noCache = TRUE)
    if (!nrow(result)) {
      return(NULL)
    }
    data.frame(
      Package = result[, "Package"],
      Version = result[, "Version"],
      Priority = result[, "Priority"],
      LibPath = library,
      stringsAsFactors = FALSE
    )
  })
  inventory <- do.call(rbind, inventories[!vapply(inventories, is.null, logical(1L))])
  inventory[!duplicated(inventory$Package), , drop = FALSE]
}

# BuildKit preserves this library between failed attempts. Remove stale rows so
# a resumed build cannot silently inherit packages outside or older than the lock.
cached <- installed_in_library(install_library)
if (nrow(cached)) {
  target_version <- setNames(lock$SelectedVersion, lock$Package)[cached$Package]
  stale <- is.na(target_version) | cached$Version != target_version
  if (any(stale)) {
    message("removing ", sum(stale), " stale cached package directories")
    for (package in cached$Package[stale]) {
      remove.packages(package, lib = install_library)
    }
  }

  rtracklayer_marker <- file.path(
    install_library,
    "rtracklayer",
    "gemcitabine-compatibility"
  )
  if (
    "rtracklayer" %in% cached$Package &&
      !file.exists(rtracklayer_marker)
  ) {
    message("removing cached rtracklayer built before the Bioconductor generic compatibility fix")
    remove.packages("rtracklayer", lib = install_library)
  }

  signac_marker <- file.path(
    install_library,
    "Signac",
    "gemcitabine-seuratobject-v4-compatibility"
  )
  if ("Signac" %in% cached$Package && !file.exists(signac_marker)) {
    message("removing cached Signac built before the registered S3 compatibility fix")
    remove.packages("Signac", lib = install_library)
  }
}

visible <- visible_installed()
visible_version <- setNames(visible$Version, visible$Package)
already_exact <- !is.na(visible_version[lock$Package]) &
  visible_version[lock$Package] == lock$SelectedVersion

sources <- lock[
  lock$SourceCategory != "provided_by_base_image" & nzchar(lock$SourceURL),
  , drop = FALSE
]
sources <- sources[!already_exact[match(sources$Package, lock$Package)], , drop = FALSE]

version_satisfies <- function(actual, operator, required) {
  comparison <- utils::compareVersion(actual, required)
  switch(
    operator,
    ">=" = comparison >= 0L,
    ">" = comparison > 0L,
    "<=" = comparison <= 0L,
    "<" = comparison < 0L,
    "==" = comparison == 0L,
    "=" = comparison == 0L,
    TRUE
  )
}

apply_source_compatibility_patches <- function(package, package_directory) {
  if (package %in% c("CNEr", "rtracklayer")) {
    # C23 gives f() a no-argument prototype.  The bundled UCSC code predates
    # that rule and calls this callback with one pointer, so state its actual
    # prototype explicitly in both declaration and definition.
    common_paths <- file.path(
      package_directory,
      "src",
      "ucsc",
      c("common.h", "common.c")
    )
    for (path in common_paths) {
      contents <- readLines(path, warn = FALSE)
      contents <- gsub(
        "void (*free)()",
        "void (*free)(void *)",
        contents,
        fixed = TRUE
      )
      writeLines(contents, path)
    }
    if (package == "CNEr") {
      makevars_path <- file.path(package_directory, "src", "Makevars")
      makevars <- if (file.exists(makevars_path)) {
        readLines(makevars_path, warn = FALSE)
      } else {
        character()
      }
      writeLines(
        c(makevars, "PKG_CFLAGS += -std=gnu17"),
        makevars_path
      )
    }
    if (package == "rtracklayer") {
      namespace_path <- file.path(package_directory, "NAMESPACE")
      namespace <- readLines(namespace_path, warn = FALSE)
      genome_info_index <- which(namespace == "import(GenomeInfoDb)")
      genomic_ranges_index <- which(namespace == "import(GenomicRanges)")
      if (length(genome_info_index) != 1L || length(genomic_ranges_index) != 1L) {
        stop("could not locate rtracklayer Bioconductor namespace imports")
      }
      namespace[c(genome_info_index, genomic_ranges_index)] <- c(
        "import(GenomicRanges)",
        "import(GenomeInfoDb)"
      )
      writeLines(namespace, namespace_path)
      marker_directory <- file.path(package_directory, "inst")
      dir.create(marker_directory, showWarnings = FALSE)
      writeLines(
        "rtracklayer 1.66.0 C23 and Bioconductor generic compatibility patch",
        file.path(marker_directory, "gemcitabine-compatibility")
      )
    }
    message("applied ", package, " UCSC callback prototype fix for C23")
    return(invisible(TRUE))
  }

  if (package != "Signac") return(invisible(FALSE))

  # The merged HPC snapshot deliberately keeps SeuratObject 4.1.4 from the
  # R-4.5 user library, but its R-4.4-only Signac 1.15.0 imports the v5-only
  # LayerData/LayerData<-/Layers API.  Preserve both selected versions and
  # provide the narrow v4 equivalents needed to build and load Signac.
  namespace_path <- file.path(package_directory, "NAMESPACE")
  namespace <- readLines(namespace_path, warn = FALSE)
  incompatible_imports <- c(
    'importFrom(SeuratObject,"LayerData<-")',
    "importFrom(SeuratObject,CheckFeaturesNames)",
    "importFrom(SeuratObject,LayerData)",
    "importFrom(SeuratObject,Layers)"
  )
  namespace <- namespace[!namespace %in% incompatible_imports]
  namespace <- c(
    namespace,
    "S3method(LayerData,default)",
    'S3method("LayerData<-",default)',
    "S3method(Layers,default)"
  )
  writeLines(namespace, namespace_path)

  r_paths <- list.files(
    file.path(package_directory, "R"),
    pattern = "[.]R$",
    full.names = TRUE
  )
  for (path in r_paths) {
    contents <- readLines(path, warn = FALSE)
    contents <- gsub("SeuratObject::LayerData", "LayerData", contents, fixed = TRUE)
    contents <- gsub("SeuratObject::Layers", "Layers", contents, fixed = TRUE)
    writeLines(contents, path)
  }

  compatibility_path <- file.path(
    package_directory,
    "R",
    "000-seuratobject-v4-compatibility.R"
  )
  writeLines(c(
    'CheckFeaturesNames <- function(data) {',
    '  if (any(grepl("_", rownames(data), fixed = TRUE))) {',
    '    warning("Feature names cannot have underscores; replacing with dashes", call. = FALSE)',
    '    rownames(data) <- gsub("_", "-", rownames(data), fixed = TRUE)',
    '  }',
    '  if (any(grepl("|", rownames(data), fixed = TRUE))) {',
    '    warning("Feature names cannot have pipe characters; replacing with dashes", call. = FALSE)',
    '    rownames(data) <- gsub("|", "-", rownames(data), fixed = TRUE)',
    '  }',
    '  data',
    '}',
    '',
    'LayerData <- function(object, ...) UseMethod("LayerData")',
    '',
    'LayerData.default <- function(object, layer = "data", assay = NULL, ...) {',
    '  arguments <- c(list(object = object, slot = layer), list(...))',
    '  if (!is.null(assay)) arguments$assay <- assay',
    '  do.call(SeuratObject::GetAssayData, arguments)',
    '}',
    '',
    '`LayerData<-` <- function(object, layer = "data", assay = NULL, ..., value) {',
    '  UseMethod("LayerData<-")',
    '}',
    '',
    '`LayerData<-.default` <- function(',
    '    object, layer = "data", assay = NULL, ..., value',
    ') {',
    '  arguments <- c(list(object = object, slot = layer, new.data = value), list(...))',
    '  if (!is.null(assay)) arguments$assay <- assay',
    '  do.call(SeuratObject::SetAssayData, arguments)',
    '}',
    '',
    'Layers <- function(object, search = NA, ...) UseMethod("Layers")',
    '',
    'Layers.default <- function(object, search = NA, ...) {',
    '  if (inherits(object, "Seurat")) {',
    '    object <- object[[SeuratObject::DefaultAssay(object)]]',
    '  }',
    '  candidates <- intersect(c("counts", "data", "scale.data"), methods::slotNames(object))',
    '  present <- vapply(candidates, function(layer) {',
    '    value <- methods::slot(object, layer)',
    '    !is.null(value) && length(value) > 0L',
    '  }, logical(1L))',
    '  candidates[present]',
    '}'
  ), compatibility_path)
  description_path <- file.path(package_directory, "DESCRIPTION")
  description <- read.dcf(description_path)
  if ("Collate" %in% colnames(description)) {
    description[1L, "Collate"] <- paste(
      "'000-seuratobject-v4-compatibility.R'",
      description[1L, "Collate"]
    )
    write.dcf(description, description_path)
  }
  marker_directory <- file.path(package_directory, "inst")
  dir.create(marker_directory, showWarnings = FALSE)
  writeLines(
    "Signac 1.15.0 registered SeuratObject 4 compatibility methods",
    file.path(
      marker_directory,
      "gemcitabine-seuratobject-v4-compatibility"
    )
  )
  message("applied Signac 1.15.0 compatibility layer for SeuratObject 4.1.4")
  invisible(TRUE)
}

install_relaxed_locked_source <- function(
    package,
    relaxed_dependencies,
    makevars_lines = character()
) {
  row <- sources[sources$Package == package, , drop = FALSE]
  if (nrow(row) != 1L || !nzchar(row$SourceURL)) {
    stop("cannot locate one exact source URL for compatibility build: ", package)
  }

  work <- tempfile(paste0(package, "-compatibility-source-"))
  dir.create(work)
  archive <- file.path(work, "source.tar.gz")
  downloaded <- FALSE
  for (attempt in seq_len(3L)) {
    downloaded <- tryCatch({
      utils::download.file(row$SourceURL, archive, mode = "wb", quiet = FALSE)
      TRUE
    }, error = function(error) FALSE)
    if (downloaded) break
    Sys.sleep(2^attempt)
  }
  if (!downloaded) {
    stop("failed to download exact source for compatibility build: ", package)
  }

  source_root <- file.path(work, "unpacked")
  dir.create(source_root)
  utils::untar(archive, exdir = source_root)
  description_paths <- list.files(
    source_root,
    pattern = "^DESCRIPTION$",
    recursive = TRUE,
    full.names = TRUE
  )
  description_packages <- vapply(description_paths, function(path) {
    value <- tryCatch(read.dcf(path, fields = "Package")[1L, 1L], error = identity)
    if (inherits(value, "error")) "" else value
  }, character(1L))
  description_path <- description_paths[description_packages == package]
  if (length(description_path) != 1L) {
    stop("could not identify package DESCRIPTION in exact source for: ", package)
  }

  apply_source_compatibility_patches(package, dirname(description_path))

  description <- read.dcf(description_path)
  for (field in intersect(c("Depends", "Imports", "LinkingTo"), colnames(description))) {
    entries <- trimws(strsplit(description[1L, field], ",", fixed = TRUE)[[1L]])
    dependency_names <- sub("[[:space:]]*\\(.*$", "", entries)
    relax <- dependency_names %in% relaxed_dependencies
    entries[relax] <- dependency_names[relax]
    description[1L, field] <- paste(entries, collapse = ", ")
  }
  write.dcf(description, description_path)
  source_md5 <- file.path(dirname(description_path), "MD5")
  if (file.exists(source_md5)) unlink(source_md5)

  if (length(relaxed_dependencies)) {
    message(
      "installing ", package,
      " with locked dependency lower-bound relaxation: ",
      paste(relaxed_dependencies, collapse = ", ")
    )
  } else {
    message("installing ", package, " from its exact source package directory")
  }
  original_makevars_user <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  compatibility_makevars <- NA_character_
  if (length(makevars_lines)) {
    compatibility_makevars <- tempfile(paste0(package, "-Makevars-"))
    writeLines(makevars_lines, compatibility_makevars)
    Sys.setenv(R_MAKEVARS_USER = compatibility_makevars)
  }
  status <- system2(
    file.path(R.home("bin"), "R"),
    c(
      "CMD", "INSTALL", "--no-lock",
      paste0("--library=", install_library),
      dirname(description_path)
    )
  )
  if (length(makevars_lines)) {
    if (is.na(original_makevars_user)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = original_makevars_user)
    }
    unlink(compatibility_makevars)
  }
  unlink(work, recursive = TRUE)
  if (!identical(status, 0L)) {
    stop("compatibility source build failed for: ", package)
  }

  installed <- installed_in_library(install_library)
  actual <- installed$Version[match(package, installed$Package)]
  if (is.na(actual) || actual != row$SelectedVersion) {
    stop(
      "compatibility build version mismatch for ", package,
      ": expected ", row$SelectedVersion, ", found ", actual
    )
  }
}

install_cxx14_locked_source <- function(package) {
  message("installing ", package, " with CXX11 mapped to GNU++14")
  install_relaxed_locked_source(
    package,
    relaxed_dependencies = character(),
    makevars_lines = "CXX11STD = -std=gnu++14"
  )
}

# pak 0.9.5 refuses to plan an explicit URL that updates an R recommended
# package already present in .Library. Install those exact sources first into
# the higher-priority site library, then let pak solve the remaining graph.
base_inventory <- installed.packages(lib.loc = .Library, noCache = TRUE)
base_priority <- setNames(base_inventory[, "Priority"], base_inventory[, "Package"])
recommended_override <- sources$Package[
  base_priority[sources$Package] %in% c("base", "recommended")
]
recommended_override <- unique(recommended_override[!is.na(recommended_override)])

for (package in recommended_override) {
  row <- sources[sources$Package == package, , drop = FALSE]
  message("installing exact R recommended override: ", package, " ", row$SelectedVersion)
  install.packages(
    row$SourceURL,
    lib = install_library,
    repos = NULL,
    type = "source",
    Ncpus = 1L
  )
  installed <- installed_in_library(install_library)
  actual <- installed$Version[match(package, installed$Package)]
  if (is.na(actual) || actual != row$SelectedVersion) {
    stop(
      "recommended override version mismatch for ", package,
      ": expected ", row$SelectedVersion, ", found ", actual,
      call. = FALSE
    )
  }
}

visible <- visible_installed()
visible_version <- setNames(visible$Version, visible$Package)
remaining <- sources[
  is.na(visible_version[sources$Package]) |
    visible_version[sources$Package] != sources$SelectedVersion,
  , drop = FALSE
]

# RcppParallel's bundled TBB invokes a nested make.  When pak builds several
# packages at once, that nested process can inherit closed GNU make jobserver
# descriptors.  Install this one source alone with an explicit single-job make;
# leave every other package free to use pak's normal parallel build settings.
if ("RcppParallel" %in% remaining$Package) {
  options(pkg.sysreqs = FALSE, pkg.sysreqs_db_update = FALSE)
  row <- remaining[remaining$Package == "RcppParallel", , drop = FALSE]
  message("installing RcppParallel with isolated single-job make")
  original_makeflags <- Sys.getenv("MAKEFLAGS", unset = NA_character_)
  Sys.setenv(MAKEFLAGS = "-j1")
  pak::pkg_install(
    row$DockerInstallRef,
    lib = install_library,
    upgrade = FALSE,
    ask = FALSE,
    dependencies = FALSE
  )
  if (is.na(original_makeflags)) {
    Sys.unsetenv("MAKEFLAGS")
  } else {
    Sys.setenv(MAKEFLAGS = original_makeflags)
  }

  visible <- visible_installed()
  visible_version <- setNames(visible$Version, visible$Package)
  remaining <- sources[
    is.na(visible_version[sources$Package]) |
      visible_version[sources$Package] != sources$SelectedVersion,
    , drop = FALSE
  ]
}

# Rhdf5lib recursively builds its bundled HDF5 source.  Parallel nested make
# can start compilation before Automake creates subdirectory .deps folders,
# producing intermittent *.Tpo "No such file or directory" failures.
if ("Rhdf5lib" %in% remaining$Package) {
  options(pkg.sysreqs = FALSE, pkg.sysreqs_db_update = FALSE)
  row <- remaining[remaining$Package == "Rhdf5lib", , drop = FALSE]
  message("installing Rhdf5lib with isolated single-job make")
  original_makeflags <- Sys.getenv("MAKEFLAGS", unset = NA_character_)
  Sys.setenv(MAKEFLAGS = "-j1")
  pak::pkg_install(
    row$DockerInstallRef,
    lib = install_library,
    upgrade = FALSE,
    ask = FALSE,
    dependencies = FALSE
  )
  if (is.na(original_makeflags)) {
    Sys.unsetenv("MAKEFLAGS")
  } else {
    Sys.setenv(MAKEFLAGS = original_makeflags)
  }

  visible <- visible_installed()
  visible_version <- setNames(visible$Version, visible$Package)
  remaining <- sources[
    is.na(visible_version[sources$Package]) |
      visible_version[sources$Package] != sources$SelectedVersion,
    , drop = FALSE
  ]
}

message("locked packages already exact: ", nrow(lock) - nrow(remaining))
message("locked source packages remaining for pak: ", nrow(remaining))
if (nrow(remaining)) {
  if (any(!nzchar(remaining$DockerInstallRef))) {
    stop("locked source package is missing DockerInstallRef", call. = FALSE)
  }
  refs <- remaining$DockerInstallRef
  # BPCells is an R package under r/ in a larger repository archive.  A plain
  # fixed tarball URL keeps the exact commit reproducible, but pak cannot infer
  # package metadata from the archive root.  Record the hard dependencies from
  # that commit's r/DESCRIPTION explicitly and install its directory source
  # through install_relaxed_locked_source() below.
  manual_dependencies_by_package <- list(
    BPCells = c(
      "Rcpp", "RcppEigen", "methods", "grDevices", "magrittr", "Matrix",
      "rlang", "tools", "vctrs", "lifecycle", "stringr", "tibble",
      "dplyr", "tidyr", "readr", "ggplot2", "scales", "patchwork",
      "scattermore", "ggrepel", "RColorBrewer", "hexbin"
    )
  )
  manual_dependencies_by_package <- manual_dependencies_by_package[
    intersect(names(manual_dependencies_by_package), remaining$Package)
  ]
  metadata_remaining <- remaining[
    !remaining$Package %in% names(manual_dependencies_by_package),
    , drop = FALSE
  ]
  # System dependencies are generated and audited separately in
  # ubuntu_noble_r_system_packages.txt.  Do not let pak reinterpret free-form
  # SystemRequirements text or mutate apt sources during the package build.
  options(
    pkg.sysreqs = FALSE,
    pkg.sysreqs_db_update = FALSE
  )

  # The captured HPC library is not a solver-consistent repository snapshot:
  # some packages declare newer minimum versions than the versions that were
  # actually installed together on the HPC.  Resolve DESCRIPTION metadata only,
  # then install in package-name dependency order without re-solving versions.
  resolve_metadata <- function(batch_refs) {
    attempts <- if (length(batch_refs) == 1L) 3L else 1L
    last_error <- NULL
    for (attempt in seq_len(attempts)) {
      result <- tryCatch(
        pak::pkg_deps(batch_refs, dependencies = FALSE),
        error = identity
      )
      if (!inherits(result, "error")) {
        return(result)
      }
      last_error <- result
      if (attempt < attempts) {
        message("retrying source metadata: ", batch_refs, " (attempt ", attempt + 1L, ")")
        Sys.sleep(2^attempt)
      }
    }
    if (length(batch_refs) == 1L) {
      stop(
        "failed to resolve source metadata after retries for ", batch_refs,
        ": ", conditionMessage(last_error),
        call. = FALSE
      )
    }
    midpoint <- length(batch_refs) %/% 2L
    message("splitting failed metadata batch of ", length(batch_refs), " sources")
    rbind(
      resolve_metadata(batch_refs[seq_len(midpoint)]),
      resolve_metadata(batch_refs[seq.int(midpoint + 1L, length(batch_refs))])
    )
  }

  metadata_refs <- metadata_remaining$DockerInstallRef
  metadata_batches <- split(metadata_refs, ceiling(seq_along(metadata_refs) / 8L))
  metadata <- do.call(rbind, lapply(seq_along(metadata_batches), function(index) {
    message(
      "resolving source metadata batch ", index, "/", length(metadata_batches),
      " (", length(metadata_batches[[index]]), " packages)"
    )
    message("metadata refs: ", paste(unname(metadata_batches[[index]]), collapse = ", "))
    resolve_metadata(unname(metadata_batches[[index]]))
  }))
  metadata <- metadata[match(metadata_remaining$Package, metadata$package), , drop = FALSE]
  metadata_version <- as.character(metadata$version)
  invalid_metadata <- is.na(metadata$package) |
    metadata$package != metadata_remaining$Package |
    metadata_version != metadata_remaining$SelectedVersion
  if (any(invalid_metadata)) {
    print(data.frame(
      ExpectedPackage = metadata_remaining$Package[invalid_metadata],
      ExpectedVersion = metadata_remaining$SelectedVersion[invalid_metadata],
      MetadataPackage = metadata$package[invalid_metadata],
      MetadataVersion = metadata_version[invalid_metadata],
      stringsAsFactors = FALSE
    ), row.names = FALSE)
    stop("pak source metadata does not match the Docker lock", call. = FALSE)
  }

  locked_version <- setNames(lock$SelectedVersion, lock$Package)
  relaxed_dependencies_by_package <- setNames(
    vector("list", nrow(remaining)),
    remaining$Package
  )
  for (index in seq_len(nrow(metadata_remaining))) {
    dependencies <- metadata$deps[[index]]
    candidates <- dependencies$type %in% c("Depends", "Imports", "LinkingTo") &
      dependencies$package %in% names(locked_version) &
      nzchar(dependencies$op) & nzchar(dependencies$version)
    if (any(candidates)) {
      candidate_rows <- which(candidates)
      incompatible <- !vapply(candidate_rows, function(dependency_index) {
        dependency <- dependencies$package[[dependency_index]]
        version_satisfies(
          locked_version[[dependency]],
          dependencies$op[[dependency_index]],
          dependencies$version[[dependency_index]]
        )
      }, logical(1L))
      relaxed_dependencies_by_package[[metadata_remaining$Package[[index]]]] <- unique(
        dependencies$package[candidate_rows[incompatible]]
      )
    }
  }
  relaxed_packages <- names(relaxed_dependencies_by_package)[vapply(
    relaxed_dependencies_by_package,
    length,
    integer(1L)
  ) > 0L]
  if (length(relaxed_packages)) {
    message(
      "locked snapshot has unsatisfied hard-dependency minimums for: ",
      paste(relaxed_packages, collapse = ", ")
    )
  }
  directory_source_packages <- union(
    relaxed_packages,
    intersect(c("BPCells", "CNEr", "rtracklayer"), remaining$Package)
  )
  # These locked releases need package-scoped C++14 flags.  Keep them inside
  # the dependency graph so LinkingTo dependencies are installed first.
  cxx14_compatibility_packages <- intersect(
    c("glmGamPoi", "parallelDist"),
    remaining$Package
  )
  c17_compatibility_packages <- intersect("CNEr", remaining$Package)
  isolated_make_packages <- intersect("duckdb", remaining$Package)
  compatibility_source_packages <- union(
    directory_source_packages,
    union(
      cxx14_compatibility_packages,
      union(c17_compatibility_packages, isolated_make_packages)
    )
  )

  install_isolated_make_source <- function(package, jobs = 4L) {
    row <- remaining[remaining$Package == package, , drop = FALSE]
    message("installing ", package, " with isolated ", jobs, "-job make")
    original_makeflags <- Sys.getenv("MAKEFLAGS", unset = NA_character_)
    on.exit({
      if (is.na(original_makeflags)) {
        Sys.unsetenv("MAKEFLAGS")
      } else {
        Sys.setenv(MAKEFLAGS = original_makeflags)
      }
    }, add = TRUE)
    Sys.setenv(MAKEFLAGS = paste0("-j", jobs))
    pak::pkg_install(
      row$DockerInstallRef,
      lib = install_library,
      upgrade = FALSE,
      ask = FALSE,
      dependencies = FALSE
    )
  }

  hard_dependency_types <- c("Depends", "Imports", "LinkingTo")
  dependencies_by_package <- setNames(vector("list", nrow(remaining)), remaining$Package)
  for (index in seq_len(nrow(metadata_remaining))) {
    dependencies <- metadata$deps[[index]]
    hard <- dependencies$package[
      dependencies$type %in% hard_dependency_types & dependencies$package != "R"
    ]
    dependencies_by_package[[metadata_remaining$Package[[index]]]] <- unique(hard)
  }
  for (package in names(manual_dependencies_by_package)) {
    dependencies_by_package[[package]] <- manual_dependencies_by_package[[package]]
  }

  visible_names <- visible_installed()$Package
  required_names <- unique(unlist(dependencies_by_package, use.names = FALSE))
  missing_hard_dependencies <- setdiff(
    required_names,
    union(remaining$Package, visible_names)
  )
  if (length(missing_hard_dependencies)) {
    stop(
      "locked sources have hard dependencies absent from the image and lock: ",
      paste(sort(missing_hard_dependencies), collapse = ", "),
      call. = FALSE
    )
  }

  pending <- remaining$Package
  dependency_levels <- list()
  while (length(pending)) {
    ready <- pending[vapply(
      pending,
      function(package) !any(dependencies_by_package[[package]] %in% pending),
      logical(1L)
    )]
    if (!length(ready)) {
      stop(
        "hard-dependency cycle among locked packages: ",
        paste(sort(pending), collapse = ", "),
        call. = FALSE
      )
    }
    dependency_levels[[length(dependency_levels) + 1L]] <- ready
    pending <- setdiff(pending, ready)
  }

  message("installing locked sources in ", length(dependency_levels), " dependency levels")
  batch_size <- 32L
  for (level_index in seq_along(dependency_levels)) {
    level <- dependency_levels[[level_index]]
    batches <- split(level, ceiling(seq_along(level) / batch_size))
    for (batch_index in seq_along(batches)) {
      packages <- unname(batches[[batch_index]])
      rows <- match(packages, remaining$Package)
      message(
        "installing dependency level ", level_index, "/", length(dependency_levels),
        ", batch ", batch_index, "/", length(batches),
        " (", length(packages), " packages)"
      )
      regular <- !packages %in% compatibility_source_packages
      for (package in packages[!regular]) {
        if (package %in% isolated_make_packages) {
          install_isolated_make_source(package)
        } else if (package %in% c17_compatibility_packages) {
          message("installing ", package, " with its legacy C sources compiled as C17")
          install_relaxed_locked_source(
            package,
            relaxed_dependencies = character(),
            makevars_lines = "C_STD = C17"
          )
        } else if (package %in% cxx14_compatibility_packages) {
          install_cxx14_locked_source(package)
        } else {
          install_relaxed_locked_source(
            package,
            relaxed_dependencies_by_package[[package]]
          )
        }
      }
      if (any(regular)) {
        pak::pkg_install(
          refs[rows[regular]],
          lib = install_library,
          upgrade = FALSE,
          ask = FALSE,
          dependencies = FALSE
        )
      }
    }
  }
}

visible <- visible_installed()
actual_version <- visible$Version[match(lock$Package, visible$Package)]
mismatch <- is.na(actual_version) | actual_version != lock$SelectedVersion
if (any(mismatch)) {
  failures <- data.frame(
    Package = lock$Package[mismatch],
    ExpectedVersion = lock$SelectedVersion[mismatch],
    ActualVersion = actual_version[mismatch],
    stringsAsFactors = FALSE
  )
  print(failures, row.names = FALSE)
  stop(sum(mismatch), " locked R package versions do not match", call. = FALSE)
}

cached <- installed_in_library(install_library)
extra_cached <- setdiff(cached$Package, lock$Package)
if (length(extra_cached)) {
  stop(
    "installation cache contains packages outside the lock: ",
    paste(extra_cached, collapse = ", "),
    call. = FALSE
  )
}

message("verified all ", nrow(lock), " locked R package versions under R 4.5.0")
