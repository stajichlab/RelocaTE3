#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --job-name=lr_small_panel
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH -o logs/lr_small_panel.%j.log

# Build the small Chr3 long-read panel. Submit from the repository root:
#   mkdir -p logs && sbatch validation/longread_small/build_panel.slurm
# Set RELOCATE3_LR_SMALL_PANEL to build somewhere other than panel_v1 here.
# Rerun-safe: finished steps are skipped; incomplete output is refused.

set -euo pipefail

BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${BASE_DIR}"
echo "start: $(date)"

THREADS="${SLURM_CPUS_PER_TASK:-1}"
SIM_PROJECT_DIR="${SIM_PROJECT_DIR:-/rhome/nmath020/bigdata/github/github_tools/data_sim/simulate_data}"
PANEL="${RELOCATE3_LR_SMALL_PANEL:-validation/longread_small/panel_v1}"
SCRIPTS="validation/longread_small"
GFF="tests/data/sim_genome/MSU7.Chr3_2M.ALL.gff"
GENOME="tests/data/sim_genome/MSU7.Chr3_2M.fa"
TE_LIB="tests/data/TE_lib/RiceTE.fa"
SHORT_R1="tests/data/sim_reads/MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"
SHORT_R2="tests/data/sim_reads/MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"

# check_panel.py refuses to run twice, so a checked panel ends the job here.
if [[ -f "${PANEL}/.complete" ]]; then
  echo "panel already built and checked: ${PANEL}"
  echo "end: $(date)"
  exit 0
fi

echo "BASE_DIR=${BASE_DIR} HOST=$(hostname)"
echo "PANEL=${PANEL} THREADS=${THREADS} SIM_PROJECT_DIR=${SIM_PROJECT_DIR}"
echo "GFF=${GFF}"
echo "GENOME=${GENOME}"
echo "TE_LIB=${TE_LIB}"
echo "SHORT_R1=${SHORT_R1}"
echo "SHORT_R2=${SHORT_R2}"

for input in "${GFF}" "${GENOME}" "${TE_LIB}" "${SHORT_R1}" "${SHORT_R2}" \
    "${SIM_PROJECT_DIR}/pyproject.toml" "${SIM_PROJECT_DIR}/src/simulate_data/cli.py"; do
  if [[ ! -f "${input}" ]]; then
    echo "ERROR: Missing input: ${input}" >&2
    exit 1
  fi
done
command -v pixi >/dev/null || { echo "ERROR: pixi not on PATH" >&2; exit 1; }

# Provenance. Both repos must be git checkouts; the counts are captured
# after the guards so a failing git command cannot hide behind the pipe.
for repo in "${BASE_DIR}" "${SIM_PROJECT_DIR}"; do
  if ! git -C "${repo}" rev-parse --git-dir >/dev/null 2>&1; then
    echo "ERROR: ${repo} is not a git checkout" >&2
    exit 1
  fi
done
RELOCATE3_COMMIT="$(git -C "${BASE_DIR}" rev-parse HEAD)"
RELOCATE3_DIRTY="$(git -C "${BASE_DIR}" status --porcelain | wc -l)"
SIM_COMMIT="$(git -C "${SIM_PROJECT_DIR}" rev-parse HEAD)"
SIM_DIRTY="$(git -C "${SIM_PROJECT_DIR}" status --porcelain | wc -l)"
echo "RelocaTE3 commit: ${RELOCATE3_COMMIT} (${RELOCATE3_DIRTY} uncommitted paths)"
echo "toolkit commit: ${SIM_COMMIT} (${SIM_DIRTY} uncommitted paths; working tree may be dirty)"

# Run a command in the toolkit's pixi env. Its editable install of
# simulate_data is missing (no import, no simulate-data launcher) and the
# toolkit is read-only to us, so put its src/ on PYTHONPATH instead; no
# bytecode is written into it.
sim() {
  PYTHONPATH="${SIM_PROJECT_DIR}/src" PYTHONDONTWRITEBYTECODE=1 \
    pixi run --as-is --manifest-path "${SIM_PROJECT_DIR}/pyproject.toml" "$@"
}

mkdir -p "${PANEL}"

# Provenance, written once and never overwritten.
if [[ ! -f "${PANEL}/build_info.tsv" ]]; then
  printf 'key\tvalue\nrelocate3_commit\t%s\nrelocate3_dirty_paths\t%s\ntoolkit_commit\t%s\ntoolkit_dirty_paths\t%s\nslurm_job_id\t%s\ndate\t%s\n' \
    "${RELOCATE3_COMMIT}" "${RELOCATE3_DIRTY}" "${SIM_COMMIT}" "${SIM_DIRTY}" \
    "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%dT%H:%M:%S%z')" \
    > "${PANEL}/build_info.tsv"
fi

echo "== step 1: te-panel config"
if [[ ! -d "${PANEL}/config" ]]; then
  # RelocaTE3's env guarantees Python >= 3.10; the node python3 may not.
  # --locked aborts rather than rewriting pixi.lock if it is out of date.
  pixi run --locked python "${SCRIPTS}/make_panel_config.py" --gff "${GFF}" --genome "${GENOME}" \
    --te-library "${TE_LIB}" --outdir "${PANEL}/config"
else
  echo "skip: ${PANEL}/config exists"
fi

# cp -n exits 0 on coreutils 8.30 even when the target differs, so compare.
if [[ -e "${PANEL}/truth_events.tsv" ]]; then
  if ! cmp -s "${PANEL}/config/truth_events.tsv" "${PANEL}/truth_events.tsv"; then
    echo "ERROR: ${PANEL}/truth_events.tsv differs from config/truth_events.tsv;" \
      "delete ${PANEL}/truth_events.tsv and resubmit" >&2
    exit 1
  fi
else
  cp "${PANEL}/config/truth_events.tsv" "${PANEL}/truth_events.tsv.tmp"
  mv "${PANEL}/truth_events.tsv.tmp" "${PANEL}/truth_events.tsv"
fi

echo "== step 2: te-panel truth genome"
if [[ ! -f "${PANEL}/truth/.complete" ]]; then
  if [[ -e "${PANEL}/truth" ]]; then
    echo "ERROR: incomplete ${PANEL}/truth; move it aside and resubmit" >&2
    exit 1
  fi
  sim python -m simulate_data.cli te-panel --config "${PANEL}/config/panel.json" \
    --output "${PANEL}/truth"
else
  echo "skip: ${PANEL}/truth/.complete exists"
fi

echo "== step 3: long reads"
sim python "${SCRIPTS}/simulate_reads.py" --panel "${PANEL}" \
  --sample-genome "${PANEL}/truth/genomes/sample_h1.fa" \
  --reference "${GENOME}" --threads "${THREADS}"

echo "== step 4: panel check"
pixi run --locked python "${SCRIPTS}/check_panel.py" --panel "${PANEL}" \
  --reference "${GENOME}" --short-r1 "${SHORT_R1}" --short-r2 "${SHORT_R2}" \
  --threads "${THREADS}"

# Record which commit finished the panel (append only; the header rows above
# describe the job that started it).
printf 'completed_relocate3_commit\t%s\ncompleted_slurm_job_id\t%s\n' \
  "${RELOCATE3_COMMIT}" "${SLURM_JOB_ID:-none}" >> "${PANEL}/build_info.tsv"

echo "end: $(date)"
