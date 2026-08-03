# Gemcitabine-model container

This directory defines Linux `amd64` analysis images. The `full` environment
uses R 4.5.0, Bioconductor 3.22, Python 3.10.13, samtools/htslib 1.23.1,
aria2c, curl, and the system CA certificate bundle. The focused
`in-vivo-cluster` environment reproduces the local R 4.5.1 setup needed by the
in-vivo scripts from stages 01 through 03.

## Files

```text
containers/
├── Dockerfile                         # Full multi-stage image definition
├── Dockerfile.in-vivo-cluster         # Focused R 4.5.1 image definition
├── build.sh                           # Shared build and verification entrypoint
├── environment.R                      # Full-image installer and verifier
├── environment-in-vivo-cluster.R      # Focused-image installer and verifier
├── packages.tsv                       # Full-image package lock
├── packages-in-vivo-cluster.tsv       # Focused complete R dependency lock
├── restore-in-vivo-cluster.R           # Exact CRAN source-version restorer
└── README.md
```

Each package manifest contains three columns:
`ecosystem`, `package`, and `version`. It records final package versions only;
it does not record where an environment or package installation originated.

Public source locations required to reproduce non-registry packages are handled
inside `environment.R`. System build and runtime libraries are declared directly
in the Dockerfile.

## Build

From the repository root:

```bash
bash containers/build.sh full
```

The default image is `gemcitabine-model:full`. Override the tag or platform with
`IMAGE_TAG` or `TARGET_PLATFORM`. For the `base` and `full` targets, a custom
certificate can be passed through the optional `CORPORATE_CA_FILE` variable; it
is mounted as a build secret and is not copied into the image. Its SHA-256
fingerprint is used only to invalidate network layers when the certificate
changes; the temporary trust entry is removed before the final image is created.

## In-vivo-cluster image

Build the focused image from the repository root:

```bash
bash containers/build.sh in-vivo-cluster
```

The default tag is `gemcitabine-model:in-vivo-cluster-r4.5.1`. The build and
runtime verification both reject any platform other than `linux/amd64`. The
Rocker base is pinned to the Linux `amd64` digest for R 4.5.1, and the image
uses Bioconductor 3.21. The manifest records the audited direct dependencies and
their complete CRAN `Depends`/`Imports`/`LinkingTo` closure from the local
`/usr/local/bin/R` installation. Runtime verification checks 187 R package
versions rather than only checking the top-level packages.

This target is intended for public distribution. It rejects
`CORPORATE_CA_FILE`, does not add a repository/account URL label, and runs a
post-build privacy audit. The audit checks image metadata, layer history, bundled
configuration files, an empty `/work` directory, common credential locations,
host paths, and the local build username. A failed audit prevents the build
command from reporting success.

The stage 01-03 scripts do not invoke Python or `reticulate`. The image therefore
provides a clean Python 3.12 virtual environment with `pip`, but intentionally
does not install Scanpy, scVelo, velocyto, or other stage-04 Python packages.

Run a script while mounting the repository read-only and a separate results
directory read-write:

```bash
mkdir -p docker-results

docker run --rm --platform linux/amd64 \
  -v "$PWD:/work:ro" \
  -v "$PWD/docker-results:/results" \
  -w /work \
  gemcitabine-model:in-vivo-cluster-r4.5.1 \
  Rscript Code/in-vivo/01_data.R
```

Input and output paths in the current R scripts/configuration are host-specific.
Mount those host directories into the container at matching paths or override
the relevant configuration/environment variables before running an analysis.
The image contains software only; it does not copy the Seurat objects or project
data into an image layer.

## Published Docker Hub image

The verified full image for Linux `amd64` is published on Docker Hub. The tag is
convenient for following future updates, while the digest-pinned reference always
selects this exact published image:

```text
zafiro/gemcitabine-model:full
zafiro/gemcitabine-model@sha256:a9dd54da6bd9e3394f4f604306c6ce8e1904649c59d86353d43ceebcd163d14b
```

Pull the exact image:

```bash
docker pull --platform linux/amd64 \
  zafiro/gemcitabine-model@sha256:a9dd54da6bd9e3394f4f604306c6ce8e1904649c59d86353d43ceebcd163d14b
```

Verify the exact image:

```bash
docker run --rm --platform linux/amd64 \
  zafiro/gemcitabine-model@sha256:a9dd54da6bd9e3394f4f604306c6ce8e1904649c59d86353d43ceebcd163d14b \
  Rscript /opt/gemcitabine-container/environment.R verify full \
  /opt/gemcitabine-container/packages.tsv
```

Run repository code with the exact image:

```bash
mkdir -p docker-results

docker run --rm --platform linux/amd64 \
  -v "$PWD:/work:ro" \
  -v "$PWD/docker-results:/results" \
  -w /work \
  zafiro/gemcitabine-model@sha256:a9dd54da6bd9e3394f4f604306c6ce8e1904649c59d86353d43ceebcd163d14b \
  Rscript path/to/script.R
```

## Verify an existing image

```bash
docker run --rm --platform linux/amd64 gemcitabine-model:full \
  Rscript /opt/gemcitabine-container/environment.R verify full \
  /opt/gemcitabine-container/packages.tsv

docker run --rm --platform linux/amd64 gemcitabine-model:full sh -c \
  'set -eu; aria2c --version; curl --version; ! ldd "$(command -v curl)" | grep -q /opt/samtools; test -s /etc/ssl/certs/ca-certificates.crt; echo ca_certificates=PASS'
```

Verification checks the architecture, exact R and Python package versions,
Bioconductor, samtools/htslib, the compiler toolchain, Signac compatibility, and
the hg38 TwoBitFile path. It also requires aria2c, curl, and a non-empty system
CA certificate bundle, and verifies that curl does not load samtools libraries.

## Run repository code

```bash
mkdir -p docker-results

docker run --rm --platform linux/amd64 \
  -v "$PWD:/work:ro" \
  -v "$PWD/docker-results:/results" \
  -w /work \
  gemcitabine-model:full \
  Rscript path/to/script.R
```

Change dependencies through `packages.tsv` and rebuild the image instead of
installing packages interactively.
