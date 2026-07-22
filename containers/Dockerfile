# syntax=docker/dockerfile:1.7

ARG TARGETPLATFORM=linux/amd64

FROM --platform=$TARGETPLATFORM python@sha256:1326d0fd281d283b077fd249e618339a44c9ca5aae6e05cb4f069a087e827922 AS python-dependencies

ENV PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONNOUSERSITE=1

COPY packages.tsv /tmp/packages.tsv

RUN --mount=type=secret,id=corporate_ca,target=/usr/local/share/ca-certificates/build-ca.crt \
    set -eux; \
    awk -F '\t' 'NR > 1 && $1 == "python" { print $2 "==" $3 }' \
        /tmp/packages.tsv > /tmp/python-requirements.txt; \
    apt-get update; \
    apt-get install --yes --no-install-recommends build-essential; \
    update-ca-certificates; \
    grep -v '^velocyto==' /tmp/python-requirements.txt > /tmp/python-before-velocyto.txt; \
    velocyto_requirement="$(grep '^velocyto==' /tmp/python-requirements.txt)"; \
    python -m pip install --no-cache-dir --requirement /tmp/python-before-velocyto.txt; \
    python -m pip install --no-cache-dir --no-build-isolation "$velocyto_requirement"; \
    python -m pip check; \
    rm -rf /var/lib/apt/lists/* /tmp/packages.tsv /tmp/python-*.txt /root/.cache

FROM --platform=$TARGETPLATFORM quay.io/biocontainers/samtools@sha256:23cda33a3a42125872766df9aaf1d2db67cdb8c85314b793465188435af31ba6 AS samtools-runtime

FROM --platform=$TARGETPLATFORM rocker/r-ver@sha256:f97b4944221bdc309583dfbc37dc0bf2654cfbd1f330cded047aeec4ae5df7e3 AS base-core

ARG TARGETPLATFORM
ARG BUILD_DATE
ARG VCS_REF

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    LD_LIBRARY_PATH=/opt/samtools/lib \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONNOUSERSITE=1 \
    TERMINFO=/opt/samtools/share/terminfo

COPY packages.tsv /opt/gemcitabine-container/packages.tsv

RUN --mount=type=secret,id=corporate_ca,target=/usr/local/share/ca-certificates/build-ca.crt \
    set -eux; \
    test "$TARGETPLATFORM" = "linux/amd64"; \
    apt-get update; \
    apt-get install --yes --no-install-recommends \
        bash ca-certificates libbz2-1.0 libexpat1 libffi8 libfontconfig1 \
        libfreetype6 libgomp1 liblzma5 libncursesw6 libreadline8t64 \
        libsqlite3-0 libssl3t64 zlib1g; \
    update-ca-certificates; \
    bioc_manager_version="$(awk -F '\t' '$1 == "r" && $2 == "BiocManager" { print $3 }' /opt/gemcitabine-container/packages.tsv)"; \
    bioc_release="$(awk -F '\t' '$1 == "r" && $2 == "BiocVersion" { sub(/\.0$/, "", $3); print $3 }' /opt/gemcitabine-container/packages.tsv)"; \
    Rscript -e 'version <- commandArgs(TRUE)[1]; install.packages(sprintf("https://cran.r-project.org/src/contrib/BiocManager_%s.tar.gz", version), repos = NULL, type = "source")' "$bioc_manager_version"; \
    Rscript -e 'args <- commandArgs(TRUE); stopifnot(as.character(packageVersion("BiocManager")) == args[1]); BiocManager::install(version = args[2], ask = FALSE, update = FALSE)' "$bioc_manager_version" "$bioc_release"; \
    rm -rf /var/lib/apt/lists/* /tmp/*

COPY --from=python-dependencies /usr/local/ /usr/local/
COPY --from=samtools-runtime /usr/local/bin/samtools /usr/local/bin/samtools
COPY --from=samtools-runtime /usr/local/lib/ /opt/samtools/lib/
COPY --from=samtools-runtime /usr/local/share/terminfo/ /opt/samtools/share/terminfo/

FROM base-core AS base

COPY environment.R /opt/gemcitabine-container/environment.R

RUN Rscript /opt/gemcitabine-container/environment.R verify base \
        /opt/gemcitabine-container/packages.tsv

LABEL org.opencontainers.image.title="Gemcitabine-model analysis base" \
      org.opencontainers.image.description="R 4.5.0, Bioconductor 3.22, Python 3.10.13, and samtools/htslib 1.23.1" \
      org.opencontainers.image.created="$BUILD_DATE" \
      org.opencontainers.image.revision="$VCS_REF" \
      org.opencontainers.image.source="https://github.com/noemiandor/Gemcitabine-model"

WORKDIR /work
CMD ["R"]

FROM base-core AS r-package-tools

ENV MAKEFLAGS=-j4 \
    PKG_NUM_PROCS=4 \
    R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

RUN --mount=type=secret,id=corporate_ca,target=/usr/local/share/ca-certificates/build-ca.crt \
    set -eux; \
    apt-get update; \
    apt-get install --yes --no-install-recommends \
        build-essential cmake curl gfortran git libabsl-dev libbz2-dev \
        libcairo2-dev libcurl4-openssl-dev libfftw3-dev libfontconfig1-dev \
        libfreetype6-dev libfribidi-dev libgdal-dev libgeos-dev libgit2-dev \
        libglpk-dev libharfbuzz-dev libhdf5-dev libicu-dev libjpeg-dev \
        liblzma-dev libmagick++-dev libnode-dev libpng-dev libproj-dev \
        libsqlite3-dev libssl-dev libtiff-dev libudunits2-dev libuv1-dev \
        libwebp-dev libx11-dev libxml2-dev libzstd-dev make pkg-config \
        zlib1g-dev gdal-bin gsfonts imagemagick-6-common libabsl20220623t64 \
        libbz2-1.0 libcairo-gobject2 libcairo2 libcurl4t64 libfftw3-double3 \
        libfftw3-long3 libfftw3-quad3 libfftw3-single3 libfontconfig1 \
        libfreetype6 libfribidi0 libgdal34t64 libgeos-c1t64 libgit2-1.7 \
        libglpk40 libharfbuzz-cairo0 libharfbuzz-gobject0 libharfbuzz-icu0 \
        libharfbuzz-subset0 libharfbuzz0b libhdf5-103-1t64 \
        libhdf5-cpp-103-1t64 libhdf5-fortran-102t64 libhdf5-hl-100t64 \
        libhdf5-hl-cpp-100t64 libhdf5-hl-fortran-100t64 libicu74 \
        libjpeg-turbo8 liblzma5 libmagick++-6.q16-9t64 \
        libmagickcore-6.q16-7-extra libmagickcore-6.q16-7t64 \
        libmagickwand-6.q16-7t64 libnode109 libpng16-16t64 libproj25 \
        libsqlite3-0 libssl3t64 libtiff6 libtiffxx6 libudunits2-0 \
        libuv1t64 libwebp7 libwebpdecoder3 libwebpdemux2 libwebpmux3 \
        libx11-6 libxml2 libzstd1 pandoc perl xz-utils zlib1g; \
    update-ca-certificates; \
    pak_version="$(awk -F '\t' '$1 == "r" && $2 == "pak" { print $3 }' /opt/gemcitabine-container/packages.tsv)"; \
    Rscript -e 'version <- commandArgs(TRUE)[1]; install.packages(sprintf("https://cran.r-project.org/src/contrib/Archive/pak/pak_%s.tar.gz", version), repos = NULL, type = "source")' "$pak_version"; \
    Rscript -e 'stopifnot(as.character(packageVersion("pak")) == commandArgs(TRUE)[1])' "$pak_version"; \
    rm -rf /var/lib/apt/lists/* /tmp/*

FROM r-package-tools AS r-package-builder

COPY environment.R /opt/gemcitabine-container/environment.R

RUN --mount=type=cache,id=gemcitabine-r-library-v1,target=/opt/r-locked-library \
    --mount=type=cache,id=gemcitabine-pak-cache-v1,target=/root/.cache/R/pkgcache \
    set -eux; \
    Rscript /opt/gemcitabine-container/environment.R install \
        /opt/gemcitabine-container/packages.tsv \
        /opt/r-locked-library; \
    cp -a /opt/r-locked-library/. /usr/local/lib/R/site-library/; \
    Rscript /opt/gemcitabine-container/environment.R verify builder \
        /opt/gemcitabine-container/packages.tsv

FROM base AS full

RUN set -eux; \
    apt-get update; \
    apt-get install --yes --no-install-recommends \
        gdal-bin gsfonts imagemagick-6-common libabsl20220623t64 libbz2-1.0 \
        libcairo-gobject2 libcairo2 libcurl4t64 libfftw3-double3 \
        libfftw3-long3 libfftw3-quad3 libfftw3-single3 libfontconfig1 \
        libfreetype6 libfribidi0 libgdal34t64 libgeos-c1t64 libgit2-1.7 \
        libglpk40 libharfbuzz-cairo0 libharfbuzz-gobject0 libharfbuzz-icu0 \
        libharfbuzz-subset0 libharfbuzz0b libhdf5-103-1t64 \
        libhdf5-cpp-103-1t64 libhdf5-fortran-102t64 libhdf5-hl-100t64 \
        libhdf5-hl-cpp-100t64 libhdf5-hl-fortran-100t64 libicu74 \
        libjpeg-turbo8 liblzma5 libmagick++-6.q16-9t64 \
        libmagickcore-6.q16-7-extra libmagickcore-6.q16-7t64 \
        libmagickwand-6.q16-7t64 libnode109 libpng16-16t64 libproj25 \
        libsqlite3-0 libssl3t64 libtiff6 libtiffxx6 libudunits2-0 \
        libuv1t64 libwebp7 libwebpdecoder3 libwebpdemux2 libwebpmux3 \
        libx11-6 libxml2 libzstd1 pandoc perl xz-utils zlib1g; \
    rm -rf /var/lib/apt/lists/* /tmp/*

COPY --from=r-package-builder /usr/local/lib/R/site-library/ /usr/local/lib/R/site-library/

RUN Rscript /opt/gemcitabine-container/environment.R verify full \
        /opt/gemcitabine-container/packages.tsv

LABEL org.opencontainers.image.title="Gemcitabine-model analysis environment" \
      org.opencontainers.image.description="R 4.5.0 analysis image with exact R and Python package versions"
