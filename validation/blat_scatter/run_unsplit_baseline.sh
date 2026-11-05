#!/usr/bin/bash -l
#SBATCH -p intel
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH -J rt3-te-unsplit
#SBATCH -o logs/blat_unsplit_baseline.%j.log
#
# Unsplit reference for the BLAT scatter (M4) gate: `relocaTE3 run` (BLAT
# TE search + trim, all chunks inside one job) on the riceTElib 30x sample.
# Its outputs (TE BAMs, te_containing/, flanking/, te_portions/) are what the
# scatter path (te-split -> te-search-chunk -> te-merge-trim) must reproduce.
#
# Usage, from the repository root:
#   sbatch validation/blat_scatter/run_unsplit_baseline.sh [OUTDIR]
set -euo pipefail

PANEL="/bigdata/wesslerlab/shared/Rice/Nathan/rice/make_simulated_genome/make_simulation_new/results/riceTElib_benchmark"
SAMPLE="${SAMPLE:-cov30x_rep1}"
R1="${PANEL}/reads/${SAMPLE}/R1.fastq.gz"
R2="${PANEL}/reads/${SAMPLE}/R2.fastq.gz"
TE_LIBRARY="${PANEL}/selected_te_library.fa"
OUT="${1:-$(pwd)/results/nextflow-m4/baseline_unsplit/${SAMPLE}}"
ENV_BIN="${PIXI_BIN:-$(pwd)/.pixi/envs/pipeline/bin}"
# Run a frozen copy of the code, not the worktree: the pixi env installs
# RelocaTE3 editable, so edits made while this ~13 h job runs would be loaded.
#   git archive <commit> src | tar -x -C <dir>
SRC="${SRC:?set SRC to a frozen src/ copy (git archive <commit> src)}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

echo "date     : $(date -Is)"
echo "src      : ${SRC}"
echo "sample   : ${SAMPLE}"
echo "outdir   : ${OUT}"
echo "threads  : ${THREADS}"
for f in "${R1}" "${R2}" "${TE_LIBRARY}" "${ENV_BIN}/relocaTE3"; do
  [[ -e "${f}" ]] || { echo "ERROR: missing ${f}" >&2; exit 1; }
done
[[ -e "${OUT}" ]] && { echo "ERROR: ${OUT} exists; refusing to overwrite" >&2; exit 1; }
mkdir -p "${OUT}"
export PATH="${ENV_BIN}:${PATH}" TMPDIR="${OUT}/.tmp" PYTHONNOUSERSITE=1 PYTHONPATH="${SRC}"
mkdir -p "${TMPDIR}"

python -c "import RelocaTE3, sys; print('RelocaTE3 from', RelocaTE3.__file__)"
/usr/bin/time -v python -m RelocaTE3 run --left "${R1}" --right "${R2}" \
  --te-library "${TE_LIBRARY}" --name "${SAMPLE}" --outdir "${OUT}" \
  --threads "${THREADS}" 2> "${OUT}/time.txt"
rmdir "${TMPDIR}" 2>/dev/null || true
grep -E "Elapsed|Maximum resident" "${OUT}/time.txt"
echo "done     : $(date -Is)"
