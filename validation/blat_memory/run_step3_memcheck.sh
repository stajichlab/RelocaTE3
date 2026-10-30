#!/usr/bin/bash -l
#SBATCH -p intel
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -J rt3-step3-mem
#SBATCH -o logs/blat_memory/step3_memcheck.%j.log
#
# Step 3 (`relocaTE3 run`: BLAT TE search + trim) on one benchmark sample,
# measuring peak memory two ways:
#   - /usr/bin/time -v "Maximum resident set size": the largest single process
#     (Python plus any in-process pysam/samtools work), NOT the sum of children;
#   - a 5 s sampler summing RSS over the whole process tree (python, blat,
#     samtools, seqtk), which is what the node actually has to provide.
# Slurm MaxRSS on this cluster includes page cache, so it is not used.
#
# Usage (from the repo root):
#   sbatch validation/blat_memory/run_step3_memcheck.sh <RT3_SRC_REPO> <OUTDIR>
# e.g. RT3_SRC_REPO=../RelocaTE3-release (old) or . (this branch).
set -euo pipefail

RT3_REPO="$(realpath "${1:?RT3 repo path (its src/ is used)}")"
OUTDIR="${2:?output directory}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
PIXI_BIN="${PIXI_BIN:-/bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/RelocaTE3_jason/RelocaTE3/results/release-validation/v0.2.0rc1/clone/.pixi/envs/blat/bin}"
PANEL="/bigdata/wesslerlab/shared/Rice/Nathan/rice/make_simulated_genome/make_simulation_new/results/riceTElib_benchmark"
SAMPLE="${SAMPLE:-cov30x_rep1}"
R1="${PANEL}/reads/${SAMPLE}/R1.fastq.gz"
R2="${PANEL}/reads/${SAMPLE}/R2.fastq.gz"
TE_LIBRARY="${PANEL}/selected_te_library.fa"

echo "date      : $(date -Is)"
echo "host      : $(hostname)  job ${SLURM_JOB_ID:-none}"
echo "rt3 repo  : ${RT3_REPO} ($(git -C "${RT3_REPO}" rev-parse --short HEAD))"
echo "sample    : ${SAMPLE}  threads ${THREADS}"
echo "outdir    : ${OUTDIR}"

for f in "${R1}" "${R2}" "${TE_LIBRARY}" "${PIXI_BIN}/python"; do
  [[ -e "${f}" ]] || { echo "ERROR: missing ${f}" >&2; exit 1; }
done
if [[ -e "${OUTDIR}/.complete" ]]; then
  echo "already complete: ${OUTDIR}"; exit 0
fi
mkdir -p "${OUTDIR}"
export PATH="${PIXI_BIN}:${PATH}"
export PYTHONPATH="${RT3_REPO}/src"
for tool in blat seqtk samtools; do command -v "${tool}" >/dev/null || { echo "ERROR: ${tool} not on PATH" >&2; exit 1; }; done

/usr/bin/time -v -o "${OUTDIR}/time.txt" \
  python -m RelocaTE3 run \
    --left "${R1}" --right "${R2}" --te-library "${TE_LIBRARY}" \
    --name "${SAMPLE}" --outdir "${OUTDIR}" --threads "${THREADS}" \
    --te-aligner blat --min-match 10 --min-trimmed 10 --mismatch 2 &
MAIN=$!

# Sum RSS (kB) over MAIN and all its descendants every 5 s.
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
