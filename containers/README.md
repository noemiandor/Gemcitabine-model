# Gemcitabine-model container

This directory defines the Linux `amd64` analysis image. The environment uses
R 4.5.0, Bioconductor 3.22, Python 3.10.13, and samtools/htslib 1.23.1.

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
copied into the image.

## Verify an existing image

```bash
docker run --rm --platform linux/amd64 gemcitabine-model:full \
  Rscript /opt/gemcitabine-container/environment.R verify full \
  /opt/gemcitabine-container/packages.tsv
```

Verification checks the architecture, exact R and Python package versions,
Bioconductor, samtools/htslib, the compiler toolchain, Signac compatibility, and
the hg38 TwoBitFile path.

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
