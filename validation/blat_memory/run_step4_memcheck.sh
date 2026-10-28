#!/usr/bin/bash -l
#SBATCH -p intel
#SBATCH --mem=48G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH -J rt3-step4-mem
#SBATCH -o logs/blat_memory/step4_memcheck.%j.log
#
# Step 4 (`relocaTE3 align-genome`, bwa aln) on cached step-3 outputs of one
# benchmark sample (relocate-benchmark checkpoint-213a0df, riceTElib
# relocate3-blat-bwaaln), measuring peak memory the same two ways as
# run_step3_memcheck.sh. The checkpoint is only read: its flanking FASTQs are
# passed by path and its te_containing/ is symlinked into OUTDIR.
#
# Usage (from the repo root):
#   sbatch validation/blat_memory/run_step4_memcheck.sh <RT3_SRC_REPO> <OUTDIR>
set -euo pipefail

RT3_REPO="$(realpath "${1:?RT3 repo path (its src/ is used)}")"
OUTDIR="$(realpath -m "${2:?output directory}")"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
PIXI_BIN="${PIXI_BIN:-/bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/RelocaTE3_jason/RelocaTE3/results/release-validation/v0.2.0rc1/clone/.pixi/envs/blat/bin}"
PANEL="/bigdata/wesslerlab/shared/Rice/Nathan/rice/make_simulated_genome/make_simulation_new"
SAMPLE="${SAMPLE:-cov30x_rep1}"
CACHE="/bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/relocate_benchmark/relocate-benchmark/runs/checkpoint-213a0df/ricetelib/relocate3-blat-bwaaln/${SAMPLE}/raw"
R1="${PANEL}/results/riceTElib_benchmark/reads/${SAMPLE}/R1.fastq.gz"
R2="${PANEL}/results/riceTElib_benchmark/reads/${SAMPLE}/R2.fastq.gz"
REFERENCE="${PANEL}/input/ref_genome/MSU_r7.fa"
FLANKS=("${CACHE}/flanking/${SAMPLE}.left.flankingReads.fq" "${CACHE}/flanking/${SAMPLE}.right.flankingReads.fq")

echo "date      : $(date -Is)"
echo "host      : $(hostname)  job ${SLURM_JOB_ID:-none}"
echo "rt3 repo  : ${RT3_REPO} ($(git -C "${RT3_REPO}" rev-parse --short HEAD))"
echo "sample    : ${SAMPLE}  threads ${THREADS}"
echo "cache     : ${CACHE}"
echo "outdir    : ${OUTDIR}"

for f in "${R1}" "${R2}" "${REFERENCE}" "${FLANKS[@]}" "${CACHE}/te_containing/${SAMPLE}.te_hit_names.txt" "${PIXI_BIN}/python"; do
  [[ -e "${f}" ]] || { echo "ERROR: missing ${f}" >&2; exit 1; }
done
if [[ -e "${OUTDIR}/.complete" ]]; then
  echo "already complete: ${OUTDIR}"; exit 0
fi
mkdir -p "${OUTDIR}"
[[ -L "${OUTDIR}/te_containing" ]] || ln -s "${CACHE}/te_containing" "${OUTDIR}/te_containing"
export PATH="${PIXI_BIN}:${PATH}"
export PYTHONPATH="${RT3_REPO}/src"

/usr/bin/time -v -o "${OUTDIR}/time.txt" \
  python -m RelocaTE3 align-genome \
    -g "${REFERENCE}" -f "${FLANKS[@]}" -n "${SAMPLE}" -o "${OUTDIR}" \
    --threads "${THREADS}" --genome-aligner bwaaln -1 "${R1}" -2 "${R2}" &
MAIN=$!

tree_rss() {
  ps -eo pid=,ppid=,rss= | awk -v root="$1" '
    { parent[$1]=$2; rss[$1]=$3 }
    END { for (p in parent) { q=p; while (q != "" && q != 0 && q != 1) { if (q == root) { s += rss[p]; break } q = parent[q] } } print s+0 }'
}
PEAK=0
while kill -0 "${MAIN}" 2>/dev/null; do
  NOW="$(tree_rss "${MAIN}")"
  if (( NOW > PEAK )); then PEAK="${NOW}"; fi
  echo -e "$(date +%s)\t${NOW}" >> "${OUTDIR}/tree_rss.tsv"
  sleep 5
done
wait "${MAIN}"
echo "tree_rss_peak_kb=${PEAK}" | tee "${OUTDIR}/tree_rss_peak.txt"
grep -E "Maximum resident|Elapsed|Exit status" "${OUTDIR}/time.txt"
touch "${OUTDIR}/.complete"
echo "done      : $(date -Is)"
