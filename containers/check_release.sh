#!/usr/bin/env bash
# Release consistency, run by .github/workflows/container.yml on v* tags only,
# after the image is built and checked and BEFORE it is pushed:
#
#   check_release.sh VERSION IMAGE_VERSION [NEXTFLOW_CONFIG]
#
#   VERSION          the release, ${GITHUB_REF_NAME#v} (v0.3.0 -> 0.3.0)
#   IMAGE_VERSION    `relocaTE3 --version` inside the built image
#   NEXTFLOW_CONFIG  default nextflow.config
#
# Fails unless IMAGE_VERSION, the config's manifest.version and its default
# params.container_tag all equal VERSION (the image tag the release publishes,
# docker/metadata-action type=pep440). Prints PASS/FAIL per item, then
# RELEASE CHECK: PASS|FAIL. tests/release_version_test.py.
set -euo pipefail

VERSION="${1:?usage: check_release.sh VERSION IMAGE_VERSION [NEXTFLOW_CONFIG]}"
IMAGE_VERSION="${2:?usage: check_release.sh VERSION IMAGE_VERSION [NEXTFLOW_CONFIG]}"
CONFIG="${3:-nextflow.config}"
[[ -f "${CONFIG}" ]] || { echo "ERROR: no ${CONFIG}" >&2; exit 1; }

# The single-quoted value of KEY inside the top-level BLOCK { ... }.
config_value() {
  awk -v block="$1" '$0 ~ "^"block"[[:space:]]*[{]" {b=1; next} b && /^[}]/ {b=0} b' "${CONFIG}" \
    | sed -nE "s/^[[:space:]]*$2[[:space:]]*=[[:space:]]*'([^']*)'.*/\1/p" | head -1
}
MANIFEST="$(config_value manifest version)"
CONTAINER_TAG="$(config_value params container_tag)"

echo "check_release: $(date '+%Y-%m-%d %H:%M:%S') release=${VERSION} config=${CONFIG}"
FAIL=0
compare() {
  if [[ "$2" == "${VERSION}" ]]; then
    printf 'PASS  %-24s %s\n' "$1" "$2"
  else
    printf 'FAIL  %-24s %s (release %s)\n' "$1" "${2:-<unset>}" "${VERSION}"; FAIL=1
  fi
}
compare "relocaTE3 --version" "${IMAGE_VERSION}"
compare "manifest.version" "${MANIFEST}"
compare "params.container_tag" "${CONTAINER_TAG}"
if [[ "${FAIL}" -eq 0 ]]; then echo "RELEASE CHECK: PASS"; else echo "RELEASE CHECK: FAIL"; fi
exit "${FAIL}"
