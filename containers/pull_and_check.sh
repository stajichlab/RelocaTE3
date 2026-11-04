#!/usr/bin/bash -l
#SBATCH -p intel
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH -J rt3-image-check
#SBATCH -o logs/image_check.%j.log
#
# Pull the published RelocaTE3 image with Apptainer, without logging in to
# ghcr.io (so this also checks that the package is public), and run
# containers/check_image.sh inside it.
#
# Usage, from the repository root:
#   sbatch containers/pull_and_check.sh [--force-pull] [TAG]   (default TAG: dev)
#
# The image is cached as a SIF under NXF_APPTAINER_CACHEDIR (default
# results/apptainer-cache), the same cache the pipeline's apptainer profile
# uses, so a later pipeline run does not pull it again.
#
# A SIF already in the cache is NOT pulled again, so the anonymous-pull (public
# package) check is skipped and a mutable tag (dev, pr-<n>) is not refreshed.
# --force-pull (or PULL=1) pulls anyway and replaces the cached SIF on success.
# Pulls go to <sif>.tmp, renamed on success; a stale .tmp is removed first and
# a failed pull's .tmp is removed on exit.
set -euo pipefail

usage() { sed -n '/^# Usage/,/^#   sbatch/{s/^# \{0,1\}//;p;}' "${BASH_SOURCE[0]}"; }
FORCE_PULL="${PULL:-0}"
TAG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --force-pull) FORCE_PULL=1 ;;
    -h|--help) usage; exit 0 ;;
    -*) echo "ERROR: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    *) [[ -z "${TAG}" ]] || { echo "ERROR: too many arguments" >&2; usage >&2; exit 2; }
       TAG="$1" ;;
  esac
  shift
done
TAG="${TAG:-dev}"
[[ "${TAG}" =~ ^[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}$ ]] || { echo "ERROR: invalid tag '${TAG}'" >&2; exit 2; }
IMAGE="docker://ghcr.io/stajichlab/relocate3:${TAG}"
CACHE="${NXF_APPTAINER_CACHEDIR:-$(pwd)/results/apptainer-cache}"
# Nextflow's cache file name for this image (ghcr.io/stajichlab/relocate3:<tag>).
SIF="${CACHE}/ghcr.io-stajichlab-relocate3-${TAG}.img"
TMP="${SIF}.tmp"

echo "date   : $(date -Is)"
echo "image  : ${IMAGE}"
echo "sif    : ${SIF}"
echo "pull   : $([[ "${FORCE_PULL}" == 1 ]] && echo forced || echo 'only if not cached')"

command -v apptainer > /dev/null || module load apptainer/1.5.3
apptainer --version
mkdir -p "${CACHE}"
unset APPTAINER_DOCKER_USERNAME APPTAINER_DOCKER_PASSWORD

if [[ -s "${SIF}" && "${FORCE_PULL}" != 1 ]]; then
  echo "SIF already present; not pulled again. The anonymous (public) pull is NOT"
  echo "checked and a mutable tag is not refreshed; rerun with --force-pull for that."
else
  trap 'rm -f "${TMP}"' EXIT
  rm -f "${TMP}"
  apptainer pull "${TMP}" "${IMAGE}"
  mv -f "${TMP}" "${SIF}"
  trap - EXIT
  echo "pulled anonymously: the package is public"
fi
ls -l "${SIF}"
apptainer inspect "${SIF}" | grep -E "org.opencontainers.image.(revision|version|created)" || true

if apptainer exec --cleanenv "${SIF}" bash /opt/relocate3/check_image.sh; then
  echo "RESULT: PASS"
else
  echo "RESULT: FAIL"; exit 1
fi
