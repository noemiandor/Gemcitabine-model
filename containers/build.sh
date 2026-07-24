#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
target="${1:-full}"
platform="${TARGET_PLATFORM:-linux/amd64}"

case "$target" in
  base|full) ;;
  *)
    echo "usage: bash containers/build.sh [base|full]" >&2
    exit 2
    ;;
esac

image_tag="${IMAGE_TAG:-gemcitabine-model:$target}"
vcs_ref="$(git -C "$repo_root" rev-parse HEAD)"
build_date="${BUILD_DATE:-$(git -C "$repo_root" show -s --format=%cI "$vcs_ref")}"
secret_args=()
corporate_ca_cache_bust=none

if [[ -n "${CORPORATE_CA_FILE:-}" ]]; then
  if [[ ! -f "$CORPORATE_CA_FILE" ]]; then
    echo "CORPORATE_CA_FILE does not exist: $CORPORATE_CA_FILE" >&2
    exit 1
  fi
  secret_args=(--secret "id=corporate_ca,src=$CORPORATE_CA_FILE")
  corporate_ca_cache_bust="$(openssl dgst -sha256 -r "$CORPORATE_CA_FILE" | awk '{print $1}')"
fi

docker buildx build \
  --platform "$platform" \
  --target "$target" \
  "${secret_args[@]}" \
  --build-arg "CORPORATE_CA_CACHE_BUST=$corporate_ca_cache_bust" \
  --build-arg "BUILD_DATE=$build_date" \
  --build-arg "VCS_REF=$vcs_ref" \
  --tag "$image_tag" \
  --load \
  --file "$script_dir/Dockerfile" \
  "$script_dir"

docker run --rm --platform "$platform" "$image_tag" \
  Rscript /opt/gemcitabine-container/environment.R verify "$target" \
  /opt/gemcitabine-container/packages.tsv

if [[ "$target" == "full" ]]; then
  docker run --rm --platform "$platform" "$image_tag" sh -c \
    'set -eu; aria2c --version; curl --version; ! ldd "$(command -v curl)" | grep -q /opt/samtools; test -s /etc/ssl/certs/ca-certificates.crt; echo ca_certificates=PASS'
fi

docker image inspect "$image_tag" \
  --format 'image={{.RepoTags}} id={{.Id}} size_bytes={{.Size}} architecture={{.Architecture}}'
