#!/usr/bin/bash -l
#SBATCH -p intel
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH -J rt3-nf-parity
#SBATCH -o logs/nextflow_parity.%j.log
#
# Parity gate: the Nextflow pipeline must produce exactly the files
# `relocaTE3 run-all` produces on the Chr3 2 Mb fixture, with the same settings
# (tests/acceptance_test.py test_run_all_cli_end_to_end; conf/test.config).
#
# Text files are compared byte for byte, BAMs record by record
# (validation/blat_memory/compare_outputs.sh); the file lists, .bai included,
# must also be equal. Prints RESULT: PASS or RESULT: FAIL; non-zero on FAIL.
#
# Usage, from the repository root:
#   sbatch tests/nextflow/parity_chr3.sh [pixi [OUTDIR]]
#   sbatch tests/nextflow/parity_chr3.sh apptainer TAG [OUTDIR]
#   bash   tests/nextflow/parity_chr3.sh docker TAG [OUTDIR]      (CI)
#
# Modes:
#   pixi       (default) run-all from the blat env (PIXI_BIN, default
#              .pixi/envs/blat/bin); the pipeline with -profile test,pixi,
#              i.e. from the pipeline env (blat env + MultiQC, pixi.toml), so
#              a PASS also shows the two envs' tools behave the same.
#   docker     run-all, the pipeline (-profile test,docker --container_tag TAG)
#   apptainer  and the comparison all run in ghcr.io/stajichlab/relocate3:TAG
#              (CONTAINER_REGISTRY overrides the repository), so both sides
#              use identical tools. apptainer is loaded with
#              `module load apptainer/1.5.3` when it is not on PATH.
# OUTDIR default: results/nextflow-parity/<date>_<jobid or pid>; an existing
# OUTDIR is refused. Nextflow comes from the pixi `nf` env in every mode.
#
# NXF_EXTRA_PARAMS="--param value ..." adds pipeline params to the nextflow
# run only (run-all is unchanged), for params that must not change the
# outputs. Forced BLAT scatter (many chunks and chunk tasks on Chr3):
#   NXF_EXTRA_PARAMS="--blat_chunk_reads 5000 --blat_chunks_per_task 3" \
#     sbatch tests/nextflow/parity_chr3.sh
# Only --name value pairs of [A-Za-z0-9_.,:+-] words; params the script sets
# (outdir, container_tag, container_registry) and Nextflow options are refused.
#
# DRY_RUN=1 checks the arguments and prints the commands without running them
# (tests/nextflow_parity_script_test.py).
set -euo pipefail

usage() {
  sed -n '/^# Usage/,/^# DRY_RUN/{/^# DRY_RUN/d;s/^# \{0,1\}//;p;}' "${BASH_SOURCE[0]}"
}
die_usage() { echo "ERROR: $*" >&2; usage >&2; exit 2; }

MODE="pixi"
TAG=""
case "${1:-}" in
  -h|--help) usage; exit 0 ;;
  "" ) ;;
  pixi) MODE="pixi"; shift ;;
  docker|apptainer)
    MODE="$1"; shift
    [[ $# -ge 1 ]] || die_usage "${MODE} mode needs an image TAG"
    TAG="$1"; shift ;;
  *) die_usage "unknown mode '$1' (pixi, docker or apptainer)" ;;
esac
[[ $# -le 1 ]] || die_usage "too many arguments"
if [[ -n "${TAG}" && ! "${TAG}" =~ ^[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}$ ]]; then
  die_usage "invalid image tag '${TAG}'"
fi

# NXF_EXTRA_PARAMS -> EXTRA_PARAMS: --name value pairs, each word checked, so
# nothing reaches the command line but plain words.
EXTRA_PARAMS=()
read -r -a EXTRA_WORDS <<< "${NXF_EXTRA_PARAMS:-}"
for (( i = 0; i < ${#EXTRA_WORDS[@]}; i += 2 )); do
  name="${EXTRA_WORDS[i]}"
  value="${EXTRA_WORDS[i + 1]:-}"
  [[ "${name}" =~ ^--[a-z][a-z0-9_]*$ ]] \
    || die_usage "NXF_EXTRA_PARAMS: '${name}' is not a --param name"
  [[ -n "${value}" && "${value}" =~ ^[A-Za-z0-9_.,:+][A-Za-z0-9_.,:+-]*$ ]] \
    || die_usage "NXF_EXTRA_PARAMS: ${name} needs a plain value, got '${value}'"
  case "${name}" in
    --outdir|--container_tag|--container_registry)
      die_usage "NXF_EXTRA_PARAMS: ${name} is set by this script" ;;
  esac
  EXTRA_PARAMS+=("${name}" "${value}")
done

DRY_RUN="${DRY_RUN:-0}"
REPO="$(pwd)"
OUT="${1:-results/nextflow-parity/$(date +%Y%m%d-%H%M%S)_${SLURM_JOB_ID:-$$}}"
OUT_ABS="$(realpath -m "${OUT}")"
BLAT_BIN="${PIXI_BIN:-${REPO}/.pixi/envs/blat/bin}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
IMAGE="${CONTAINER_REGISTRY:-ghcr.io/stajichlab/relocate3}:${TAG}"
IMAGE_ENV_BIN="/opt/relocate3/.pixi/envs/pipeline/bin"   # containers/Dockerfile ENV_DIR
DATA="${REPO}/tests/data"
R1="${DATA}/sim_reads/MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"
R2="${DATA}/sim_reads/MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"
TE="${DATA}/TE_lib/RiceTE.fa"
GENOME="${DATA}/sim_genome/MSU7.Chr3_2M.fa"
RM="${DATA}/sim_genome/MSU7.Chr3_2M.fa.RepeatMasker.out"

echo "date      : $(date -Is)"
echo "repo      : ${REPO} ($(git -C "${REPO}" describe --tags --always --dirty 2>/dev/null || echo unknown))"
echo "mode      : ${MODE}"
[[ "${MODE}" == pixi ]] && echo "blat env  : ${BLAT_BIN}" || echo "image     : ${IMAGE}"
echo "outdir    : ${OUT}"
echo "threads   : ${THREADS}"
echo "extra     : ${EXTRA_PARAMS[*]:-(none)}"
[[ "${DRY_RUN}" == 1 ]] && echo "dry run   : commands are printed, not run"

for f in "${R1}" "${R2}" "${TE}" "${GENOME}" "${RM}" main.nf; do
  [[ -f "${f}" ]] || { echo "ERROR: missing ${f} (run from the repository root)" >&2; exit 1; }
done
if [[ -e "${OUT}" ]]; then
  echo "ERROR: ${OUT} exists; refusing to overwrite" >&2; exit 1
fi

# Prints a command on one line, quoting only the words that need it.
show() {
  local w out=()
  for w in "$@"; do
    if [[ "${w}" =~ ^[A-Za-z0-9_./:=,@%+-]+$ ]]; then out+=("${w}"); else out+=("$(printf '%q' "${w}")"); fi
  done
  echo "${out[*]}"
}
# Runs "$@", or prints it under DRY_RUN=1.
run() {
  if [[ "${DRY_RUN}" == 1 ]]; then
    show "$@"
  else
    "$@"
  fi
}
# run_log LOG CMD...: run CMD with stdout and stderr to LOG.
run_log() {
  local log="$1"; shift
  if [[ "${DRY_RUN}" == 1 ]]; then
    echo "$(show "$@") > ${log} 2>&1"
  else
    "$@" > "${log}" 2>&1
  fi
}

# Prefix that runs a command inside the image (empty in pixi mode). The repo
# and OUTDIR are mounted at their host paths, so every path is the same inside.
IN_IMAGE=()
case "${MODE}" in
  docker)
    IN_IMAGE=(docker run --rm -u "$(id -u):$(id -g)"
              -v "${REPO}:${REPO}" -v "${OUT_ABS}:${OUT_ABS}" -w "${REPO}" "${IMAGE}") ;;
  apptainer)
    IN_IMAGE=(apptainer exec -B "${REPO}" -B "${OUT_ABS}" --pwd "${REPO}" "docker://${IMAGE}") ;;
esac

if [[ "${DRY_RUN}" != 1 ]]; then
  command -v pixi >/dev/null || { echo "ERROR: pixi not on PATH" >&2; exit 1; }
  case "${MODE}" in
    pixi)
      [[ -x "${BLAT_BIN}/relocaTE3" ]] \
        || { echo "ERROR: no relocaTE3 in ${BLAT_BIN}; run 'pixi install -e blat'" >&2; exit 1; } ;;
    docker)
      command -v docker >/dev/null || { echo "ERROR: docker not on PATH" >&2; exit 1; } ;;
    apptainer)
      if ! command -v apptainer >/dev/null; then
        module load apptainer/1.5.3 \
          || { echo "ERROR: apptainer not on PATH and 'module load apptainer/1.5.3' failed" >&2; exit 1; }
      fi ;;
  esac
  mkdir -p "${OUT}/ref"
fi

# 1. run-all, on a copy of the genome (run-all indexes in place).
RUN_ALL_ARGS=(
  -1 "${R1}" -2 "${R2}" -T "${TE}" -g "${OUT_ABS}/ref/$(basename "${GENOME}")"
  -n chr3 -o "${OUT_ABS}/run_all/chr3" --threads "${THREADS}"
  --mismatch 2 --tsd UNK --te-name repeat --require-both-junctions
  --repeatmasker "${RM}" --genotype
)
run cp "${GENOME}" "${OUT}/ref/"
echo "== run-all"
if [[ "${MODE}" == pixi ]]; then
  run_log "${OUT}/run_all.log" env PATH="${BLAT_BIN}:${PATH}" relocaTE3 run-all "${RUN_ALL_ARGS[@]}"
else
  run_log "${OUT}/run_all.log" "${IN_IMAGE[@]}" relocaTE3 run-all "${RUN_ALL_ARGS[@]}"
fi

# 2. The pipeline, test profile (same inputs and settings).
echo "== nextflow"
NF_PROFILE=(-profile "test,${MODE}")
[[ "${MODE}" == pixi ]] || NF_PROFILE+=(--container_tag "${TAG}")
[[ -z "${CONTAINER_REGISTRY:-}" ]] || NF_PROFILE+=(--container_registry "${CONTAINER_REGISTRY}")
run_log "${OUT}/nextflow.log" pixi run -e nf nextflow run . "${NF_PROFILE[@]}" \
  --outdir "${OUT}/nf" -work-dir "${OUT}/work" ${EXTRA_PARAMS[@]+"${EXTRA_PARAMS[@]}"}
if [[ "${DRY_RUN}" != 1 ]]; then
  # How the TE search ran: one TE_SEARCH task, or TE_SPLIT + n chunk tasks + merge.
  TRACE="${OUT}/nf/pipeline_info/execution_trace.txt"
  for p in RELOCATE3_TE_SEARCH RELOCATE3_TE_SPLIT RELOCATE3_TE_SEARCH_CHUNK RELOCATE3_TE_MERGE_TRIM; do
    echo "tasks     : ${p} $(awk -F'\t' -v p="${p}" '{ split($4, w, " "); n = split(w[1], q, ":"); if (q[n] == p) c++ } END { print c + 0 }' "${TRACE}")"
  done
fi

# 3. Compare. In container mode samtools (compare_outputs.sh) is the image's.
echo "== compare"
if [[ "${DRY_RUN}" != 1 ]]; then
  (cd "${OUT}/run_all/chr3" && find . -type f | sort) > "${OUT}/files_run_all.txt"
  (cd "${OUT}/nf/chr3" && find . -type f | sort) > "${OUT}/files_nf.txt"
  if ! diff "${OUT}/files_run_all.txt" "${OUT}/files_nf.txt"; then
    echo "RESULT: FAIL (file lists differ)"; exit 1
  fi
  echo "files     : $(wc -l < "${OUT}/files_nf.txt") in both"
fi
COMPARE=(bash validation/blat_memory/compare_outputs.sh "${OUT_ABS}/run_all/chr3" "${OUT_ABS}/nf/chr3")
if [[ "${MODE}" == pixi ]]; then
  COMPARE=(env PIXI_BIN="${BLAT_BIN}" "${COMPARE[@]}")
else
  COMPARE=("${IN_IMAGE[@]}" env PIXI_BIN="${IMAGE_ENV_BIN}" "${COMPARE[@]}")
fi
if run "${COMPARE[@]}"; then
  [[ "${DRY_RUN}" == 1 ]] || echo "RESULT: PASS"
else
  echo "RESULT: FAIL (contents differ)"; exit 1
fi
