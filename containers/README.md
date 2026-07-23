# Gemcitabine-model container

This directory defines the Linux `amd64` analysis image. The environment uses
R 4.5.0, Bioconductor 3.22, Python 3.10.13, samtools/htslib 1.23.1, aria2c,
curl, and the system CA certificate bundle.

## Files

```text
containers/
├── Dockerfile       # Multi-stage image definition
├── build.sh         # Build and verification entrypoint
├── environment.R    # R package installer and complete environment verifier
├── packages.tsv     # R and Python package names and exact versions
└── README.md
```

`packages.tsv` is the only package version input. It contains three columns:
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
`IMAGE_TAG` or `TARGET_PLATFORM`. A custom certificate can be passed through the
optional `CORPORATE_CA_FILE` variable; it is mounted as a build secret and is not
copied into the image. Its SHA-256 fingerprint is used only to invalidate network
layers when the certificate changes; the temporary trust entry is removed before
the final image is created.

## Published Docker Hub image

The verified full image for Linux `amd64` is published on Docker Hub. The tag is
convenient for following future updates, while the digest-pinned reference always
selects this exact published image:

```text
zafiro/gemcitabine-model:full
zafiro/gemcitabine-model@sha256:ff1edc60b05ac1e0303a25799efe8c2df5eef22ed875109d12bf490b209668c8
```

Pull the exact image:

```bash
docker pull --platform linux/amd64 \
  zafiro/gemcitabine-model@sha256:ff1edc60b05ac1e0303a25799efe8c2df5eef22ed875109d12bf490b209668c8
```

Verify the exact image:

```bash
docker run --rm --platform linux/amd64 \
  zafiro/gemcitabine-model@sha256:ff1edc60b05ac1e0303a25799efe8c2df5eef22ed875109d12bf490b209668c8 \
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
  zafiro/gemcitabine-model@sha256:ff1edc60b05ac1e0303a25799efe8c2df5eef22ed875109d12bf490b209668c8 \
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
