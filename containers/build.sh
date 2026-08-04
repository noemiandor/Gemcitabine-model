#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
target="${1:-full}"
platform="${TARGET_PLATFORM:-linux/amd64}"
dockerfile="$script_dir/Dockerfile"
default_image_tag="gemcitabine-model:$target"

audit_public_image() {
  local image="$1"
  local host_home="${HOME:-}"
  local host_user
  local image_metadata
  local image_history
  host_user="$(id -un)"
  image_metadata="$(docker image inspect "$image" --format '{{json .Config}}')"
  image_history="$(docker history --no-trunc --format '{{.CreatedBy}}' "$image")"

  for value in "$host_home" /Users/ /Volumes/ github.com/; do
    [[ -z "$value" ]] && continue
    if grep -Fqi -- "$value" <<<"$image_metadata"$'\n'"$image_history"; then
      echo "public image audit failed: metadata or history contains forbidden host/account data" >&2
      exit 1
    fi
  done
  if [[ "$host_user" != "root" ]] && grep -Fqi -- "$host_user" <<<"$image_metadata"$'\n'"$image_history"; then
    echo "public image audit failed: metadata or history contains the host username" >&2
    exit 1
  fi
  if grep -Eiq '(password|passwd|token|api[_-]?key|private[_-]?key|secret)[[:space:]]*[:=][[:space:]]*[^[:space:],}"]+' \
      <<<"$image_metadata"$'\n'"$image_history"; then
    echo "public image audit failed: metadata or history resembles an embedded credential" >&2
    exit 1
  fi

  docker run --rm --platform "$platform" "$image" sh -eu -c '
    for path in \
      /Users /Volumes /root/.ssh /root/.aws /root/.docker /root/.netrc \
      /root/.gitconfig /root/.Renviron /root/.Rhistory /root/.config/gh \
      /usr/local/share/ca-certificates/build-ca.crt /tmp/build-ca-bundle.crt; do
      if [ -e "$path" ]; then
        echo "public image audit failed: forbidden path exists: $path" >&2
        exit 1
      fi
    done
    if find /work -mindepth 1 -print -quit | grep -q .; then
      echo "public image audit failed: /work is not empty" >&2
      exit 1
    fi
    for value in "$1" "$2" /Users/ /Volumes/ github.com/; do
      [ -z "$value" ] && continue
      if grep -IRFqi -- "$value" /opt/gemcitabine-container; then
        echo "public image audit failed: bundled files contain host/account data" >&2
        exit 1
      fi
    done
    if grep -IREiq \
        "(password|passwd|token|api[_-]?key|private[_-]?key|secret)[[:space:]]*[:=][[:space:]]*[^[:space:],}\"]+" \
        /opt/gemcitabine-container; then
      echo "public image audit failed: bundled files resemble an embedded credential" >&2
      exit 1
    fi
    echo public_image_audit=PASS
  ' sh "$host_user" "$host_home"
}

case "$target" in
  base|full) ;;
  in-vivo-cluster)
    dockerfile="$script_dir/Dockerfile.in-vivo-cluster"
    default_image_tag="gemcitabine-model:in-vivo-cluster-r4.5.1"
    ;;
  *)
    echo "usage: bash containers/build.sh [base|full|in-vivo-cluster]" >&2
    exit 2
    ;;
esac

if [[ "$platform" != "linux/amd64" ]]; then
  echo "unsupported platform: $platform (only linux/amd64 is allowed)" >&2
  exit 2
fi

image_tag="${IMAGE_TAG:-$default_image_tag}"
vcs_ref="$(git -C "$repo_root" rev-parse HEAD)"
build_date="${BUILD_DATE:-$(git -C "$repo_root" show -s --format=%cI "$vcs_ref")}"
corporate_ca_cache_bust=none
build_args=(
  docker buildx build
  --platform "$platform"
  --target "$target"
)

if [[ "$target" == "in-vivo-cluster" && -n "${CORPORATE_CA_FILE:-}" ]]; then
  echo "CORPORATE_CA_FILE is not allowed for the public in-vivo-cluster image" >&2
  exit 2
fi

if [[ -n "${CORPORATE_CA_FILE:-}" ]]; then
  if [[ ! -f "$CORPORATE_CA_FILE" ]]; then
    echo "CORPORATE_CA_FILE does not exist: $CORPORATE_CA_FILE" >&2
    exit 1
  fi
  build_args+=(--secret "id=corporate_ca,src=$CORPORATE_CA_FILE")
  corporate_ca_cache_bust="$(openssl dgst -sha256 -r "$CORPORATE_CA_FILE" | awk '{print $1}')"
fi

if [[ "$target" != "in-vivo-cluster" ]]; then
  build_args+=(--build-arg "CORPORATE_CA_CACHE_BUST=$corporate_ca_cache_bust")
fi

build_args+=(
  --build-arg "BUILD_DATE=$build_date"
  --build-arg "VCS_REF=$vcs_ref"
  --tag "$image_tag"
  --load
  --file "$dockerfile"
  "$script_dir"
)

"${build_args[@]}"

if [[ "$target" == "in-vivo-cluster" ]]; then
  docker run --rm --platform "$platform" "$image_tag" \
    Rscript /opt/gemcitabine-container/environment-in-vivo-cluster.R verify \
    /opt/gemcitabine-container/packages-in-vivo-cluster.tsv
  audit_public_image "$image_tag"
else
  docker run --rm --platform "$platform" "$image_tag" \
    Rscript /opt/gemcitabine-container/environment.R verify "$target" \
    /opt/gemcitabine-container/packages.tsv
fi

if [[ "$target" == "full" ]]; then
  docker run --rm --platform "$platform" "$image_tag" sh -c \
    'set -eu; aria2c --version; curl --version; ! ldd "$(command -v curl)" | grep -q /opt/samtools; test -s /etc/ssl/certs/ca-certificates.crt; echo ca_certificates=PASS'
fi

docker image inspect "$image_tag" \
  --format 'image={{.RepoTags}} id={{.Id}} size_bytes={{.Size}} architecture={{.Architecture}}'
