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
#   sbatch containers/pull_and_check.sh [TAG]      (default TAG: dev)
#
# The image is cached as a SIF under NXF_APPTAINER_CACHEDIR (default
# results/apptainer-cache), the same cache the pipeline's apptainer profile
# uses, so a later pipeline run does not pull it again.
set -euo pipefail

TAG="${1:-dev}"
IMAGE="docker://ghcr.io/stajichlab/relocate3:${TAG}"
CACHE="${NXF_APPTAINER_CACHEDIR:-$(pwd)/results/apptainer-cache}"
# Nextflow's cache file name for this image (ghcr.io/stajichlab/relocate3:<tag>).
SIF="${CACHE}/ghcr.io-stajichlab-relocate3-${TAG}.img"

echo "date   : $(date -Is)"
echo "image  : ${IMAGE}"
echo "sif    : ${SIF}"

module load apptainer/1.5.3
apptainer --version
mkdir -p "${CACHE}"
unset APPTAINER_DOCKER_USERNAME APPTAINER_DOCKER_PASSWORD

if [[ -s "${SIF}" ]]; then
  echo "SIF already present; not pulling again"
else
  # Build into a temporary name, then rename, so a failed pull leaves nothing.
  apptainer pull "${SIF}.tmp" "${IMAGE}"
  mv "${SIF}.tmp" "${SIF}"
fi
ls -l "${SIF}"
apptainer inspect "${SIF}" | grep -E "org.opencontainers.image.(revision|version|created)" || true

if apptainer exec --cleanenv "${SIF}" bash /opt/relocate3/check_image.sh; then
  echo "RESULT: PASS"
else
  echo "RESULT: FAIL"; exit 1
fi
