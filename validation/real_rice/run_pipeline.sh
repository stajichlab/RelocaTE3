#!/usr/bin/bash -l
#SBATCH -p intel
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH -J rt3-nf-real-rice
#SBATCH -o logs/real_rice_pipeline.%j.log
#
# Real-rice validation through the Nextflow pipeline: the head job. Nextflow
# runs here and submits one SLURM job per task (-profile hpcc,apptainer).
#
# Usage, from the repository root:
#   sbatch validation/real_rice/run_pipeline.sh CONTAINER_TAG [OUTDIR]
#
#   CONTAINER_TAG  image tag of ghcr.io/stajichlab/relocate3: sha-<short> or a
#                  release version. dev and pr-<n> are refused (CI re-pushes
#                  them, so -resume could mix images) unless ALLOW_MUTABLE_TAG=1
#   OUTDIR         default results/nextflow-real-rice/CONTAINER_TAG
#
# Settings come from the config run_relocate3.sh (the run-all validation) reads:
# CONFIG, default validation/real_rice/config.toml, else config.example.toml.
# make_samplesheet.py --print-params maps it to pipeline params (--size from
# insert_size, --te_aligner from aligner, --genotype from characterize, ...).
# Settings run_relocate3.sh does not pass (te_opts, genome_opts, distance,
# require_both_junctions, excision) stay at the pipeline defaults, which equal
# run-all's (tests/nextflow_schema_test.py). Each sample's genotype_bam is the
# config's external CRAM, as --genotype-bam in run_relocate3.sh.
#
# Layout:
#   OUTDIR/                   published results, samplesheet.csv, container_tag.txt,
#                             container_image.txt (one record per run, after it:
#                             tag, cached SIF sha256, `apptainer inspect` OCI
#                             labels incl. org.opencontainers.image.revision)
#   OUTDIR.nextflow/          launch directory: .nextflow/ (run history), .nextflow.log
#   OUTDIR.nextflow/work/     Nextflow work directory
#
# Rerunning: the run uses -resume, so rerunning the same command after a failure
# or a time-out reuses every finished task. Nextflow resumes the last run
# recorded in the launch directory's .nextflow/history, which is why each
# OUTDIR has its own launch directory. A rerun with an existing OUTDIR is the
# intended use; refused are an OUTDIR made with a different container tag
# (container_tag.txt) and a sample sheet that no longer matches the config
# (the inputs changed: use a new OUTDIR).
#
# Environment:
#   NXF_APPTAINER_CACHEDIR  image cache, default results/apptainer-cache
#   CONFIG                  real-rice config (see above)
#   PYTHON                  Python >= 3.11 for the config, default .pixi/envs/blat/bin/python
#   DRY_RUN=1               check the inputs and print the command; run nothing
#   ALLOW_MUTABLE_TAG=1     accept dev / pr-<n> (container_image.txt still
#                           records which image each run got)
set -euo pipefail

usage() {
  sed -n '/^# Usage/,/^#   OUTDIR /{s/^# \{0,1\}//;p;}' "${BASH_SOURCE[0]}"
  echo "  Reruns resume (-resume); see the script header."
}
die_usage() { echo "ERROR: $*" >&2; usage >&2; exit 2; }

case "${1:-}" in
  -h|--help) usage; exit 0 ;;
  "") die_usage "missing CONTAINER_TAG" ;;
  -*) die_usage "unknown option '$1'" ;;
esac
[[ $# -le 2 ]] || die_usage "too many arguments"
TAG="$1"
[[ "${TAG}" =~ ^[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}$ ]] || die_usage "invalid image tag '${TAG}'"
# dev and pr-<n> move with every push (.github/workflows/container.yml).
MUTABLE=0
if [[ "${TAG}" == dev || "${TAG}" =~ ^pr-[0-9]+$ ]]; then
  MUTABLE=1
  if [[ "${ALLOW_MUTABLE_TAG:-0}" != 1 ]]; then
    echo "ERROR: container tag '${TAG}' is mutable (CI re-pushes it), so a rerun or -resume" >&2
    echo "       could use a different image. Use sha-<short> (the commit's image, e.g." >&2
    echo "       sha-$(git rev-parse --short=7 HEAD 2>/dev/null || echo '<commit>')) or a release version," >&2
    echo "       or set ALLOW_MUTABLE_TAG=1." >&2
    exit 2
  fi
fi

DRY_RUN="${DRY_RUN:-0}"
REPO="$(pwd)"
HARNESS="${REPO}/validation/real_rice"
[[ -f "${REPO}/main.nf" && -f "${HARNESS}/make_samplesheet.py" ]] \
  || { echo "ERROR: run from the repository root (no main.nf or ${HARNESS}/make_samplesheet.py)" >&2; exit 1; }
if [[ -z "${CONFIG:-}" ]]; then
  CONFIG="${HARNESS}/config.toml"
  [[ -f "${CONFIG}" ]] || CONFIG="${HARNESS}/config.example.toml"
fi
PYTHON="${PYTHON:-${REPO}/.pixi/envs/blat/bin/python}"
[[ -x "${PYTHON}" ]] || { echo "ERROR: no Python at ${PYTHON}; run 'pixi install -e blat' or set PYTHON" >&2; exit 1; }
OUTDIR="$(realpath -m "${2:-results/nextflow-real-rice/${TAG}}")"
LAUNCH="${OUTDIR}.nextflow"
WORK="${LAUNCH}/work"
SHEET="${OUTDIR}/samplesheet.csv"
export NXF_APPTAINER_CACHEDIR="${NXF_APPTAINER_CACHEDIR:-${REPO}/results/apptainer-cache}"

mapfile -t PARAM_LINES < <("${PYTHON}" "${HARNESS}/make_samplesheet.py" --config "${CONFIG}" --print-params)
[[ ${#PARAM_LINES[@]} -gt 0 ]] || { echo "ERROR: no parameters read from ${CONFIG}" >&2; exit 1; }

echo "date            : $(date -Is)"
echo "repo            : ${REPO} ($(git -C "${REPO}" describe --tags --always --dirty 2>/dev/null || echo unknown))"
echo "job             : ${SLURM_JOB_ID:-none}"
echo "container tag   : ${TAG}$([[ ${MUTABLE} == 1 ]] && echo ' (mutable; ALLOW_MUTABLE_TAG=1)')"
echo "config          : ${CONFIG}"
echo "outdir          : ${OUTDIR}"
echo "launch dir      : ${LAUNCH}"
echo "work dir        : ${WORK}"
echo "apptainer cache : ${NXF_APPTAINER_CACHEDIR}"
NF_PARAMS=()
for line in "${PARAM_LINES[@]}"; do
  name="${line%%$'\t'*}"; value="${line#*$'\t'}"
  printf '  %-14s: %s\n' "${name}" "${value}"
  NF_PARAMS+=("--${name}" "${value}")
done
[[ "${DRY_RUN}" == 1 ]] && echo "dry run         : the command is printed, not run"

# Inputs (all reported together), before anything is created.
if [[ "${DRY_RUN}" == 1 ]]; then
  "${PYTHON}" "${HARNESS}/make_samplesheet.py" --config "${CONFIG}" -o - > /dev/null
fi
TAG_FILE="${OUTDIR}/container_tag.txt"
IMAGE_FILE="${OUTDIR}/container_image.txt"
# Nextflow's file name for ghcr.io/stajichlab/relocate3:<tag> in its cache.
SIF="${NXF_APPTAINER_CACHEDIR}/ghcr.io-stajichlab-relocate3-${TAG}.img"
if [[ -f "${TAG_FILE}" && "$(cat "${TAG_FILE}")" != "${TAG}" ]]; then
  echo "ERROR: ${OUTDIR} holds a run with container tag $(cat "${TAG_FILE}"), not ${TAG}; use another OUTDIR" >&2
  exit 1
fi

NF_CMD=(pixi run --manifest-path "${REPO}/pixi.toml" --frozen -e nf nextflow
  run "${REPO}" -profile "hpcc,apptainer" --container_tag "${TAG}"
  --input "${SHEET}" --outdir "${OUTDIR}" -work-dir "${WORK}"
  "${NF_PARAMS[@]}" -resume)

if [[ "${DRY_RUN}" == 1 ]]; then
  show() {
    local w out=()
    for w in "$@"; do
      if [[ "${w}" =~ ^[A-Za-z0-9_./:=,@%+-]+$ ]]; then out+=("${w}"); else out+=("$(printf '%q' "${w}")"); fi
    done
    echo "${out[*]}"
  }
  show "${NF_CMD[@]}"
  exit 0
fi

if ! command -v pixi >/dev/null; then
  if [[ -x "${HOME}/.pixi/bin/pixi" ]]; then
    export PATH="${HOME}/.pixi/bin:${PATH}"
  else
    echo "ERROR: pixi not on PATH" >&2; exit 1
  fi
fi
if ! command -v apptainer >/dev/null; then
  module load apptainer/1.5.3 \
    || { echo "ERROR: apptainer not on PATH and 'module load apptainer/1.5.3' failed" >&2; exit 1; }
fi
echo "pixi            : $(command -v pixi)"
echo "apptainer       : $(apptainer --version)"

mkdir -p "${OUTDIR}" "${LAUNCH}" "${NXF_APPTAINER_CACHEDIR}"
"${PYTHON}" "${HARNESS}/make_samplesheet.py" --config "${CONFIG}" -o "${SHEET}"
[[ -f "${TAG_FILE}" ]] || echo "${TAG}" > "${TAG_FILE}"

# Which image this run used: Nextflow pulls the tag into the cache (or reuses
# the SIF already there), so after the run the cached SIF is the image. One
# record per run, appended: a rerun shows if the image changed.
record_image() {
  {
    echo "== $(date -Is) job=${SLURM_JOB_ID:-none} nextflow_exit=$1"
    echo "tag: ${TAG}"
    echo "image: ghcr.io/stajichlab/relocate3:${TAG}"
    if [[ -s "${SIF}" ]]; then
      echo "sif: ${SIF}"
      echo "sif_sha256: $(sha256sum "${SIF}" | cut -d' ' -f1)"
      apptainer inspect "${SIF}" 2>&1 | grep -E 'org\.opencontainers\.image\.(revision|version|created|source)' \
        || echo "labels: none (apptainer inspect ${SIF})"
    else
      echo "sif: not found at ${SIF} (NXF_APPTAINER_CACHEDIR changed?)"
    fi
  } >> "${IMAGE_FILE}"
}

echo "== nextflow ($(date -Is))"
cd "${LAUNCH}"
set +e
"${NF_CMD[@]}"
NF_RC=$?
set -e
record_image "${NF_RC}"
echo "image record    : ${IMAGE_FILE}"
tail -n 8 "${IMAGE_FILE}"
if [[ "${NF_RC}" -ne 0 ]]; then
  echo "== nextflow FAILED (exit ${NF_RC}, $(date -Is)); rerun the same command to resume" >&2
  exit "${NF_RC}"
fi
echo "== done ($(date -Is)); compare with:"
echo "   sbatch validation/real_rice/compare_pipeline_to_rc2.sh ${OUTDIR}"
