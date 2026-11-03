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
# must also be equal. Exits non-zero on any difference.
#
# Usage, from the repository root:
#   sbatch tests/nextflow/parity_chr3.sh [OUTDIR]
#   (default OUTDIR: results/nextflow-parity/<date>_<jobid or pid>)
set -euo pipefail

REPO="$(pwd)"
OUT="${1:-results/nextflow-parity/$(date +%Y%m%d-%H%M%S)_${SLURM_JOB_ID:-$$}}"
BLAT_BIN="${PIXI_BIN:-${REPO}/.pixi/envs/blat/bin}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
DATA="${REPO}/tests/data"
R1="${DATA}/sim_reads/MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"
R2="${DATA}/sim_reads/MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"
TE="${DATA}/TE_lib/RiceTE.fa"
GENOME="${DATA}/sim_genome/MSU7.Chr3_2M.fa"
RM="${DATA}/sim_genome/MSU7.Chr3_2M.fa.RepeatMasker.out"

echo "date      : $(date -Is)"
echo "repo      : ${REPO} ($(git -C "${REPO}" describe --tags --always --dirty))"
echo "outdir    : ${OUT}"
echo "blat env  : ${BLAT_BIN}"
echo "threads   : ${THREADS}"

for f in "${R1}" "${R2}" "${TE}" "${GENOME}" "${RM}" main.nf; do
  [[ -f "${f}" ]] || { echo "ERROR: missing ${f} (run from the repository root)" >&2; exit 1; }
done
command -v pixi >/dev/null || { echo "ERROR: pixi not on PATH" >&2; exit 1; }
[[ -x "${BLAT_BIN}/relocaTE3" ]] || { echo "ERROR: no relocaTE3 in ${BLAT_BIN}; run 'pixi install -e blat'" >&2; exit 1; }
if [[ -e "${OUT}" ]]; then
  echo "ERROR: ${OUT} exists; refusing to overwrite" >&2; exit 1
fi
mkdir -p "${OUT}/ref"

# 1. run-all, on a copy of the genome (run-all indexes in place).
cp "${GENOME}" "${OUT}/ref/"
echo "== run-all"
PATH="${BLAT_BIN}:${PATH}" relocaTE3 run-all \
  -1 "${R1}" -2 "${R2}" -T "${TE}" -g "${OUT}/ref/$(basename "${GENOME}")" \
  -n chr3 -o "${OUT}/run_all/chr3" --threads "${THREADS}" \
  --mismatch 2 --tsd UNK --te-name repeat --require-both-junctions \
  --repeatmasker "${RM}" --genotype > "${OUT}/run_all.log" 2>&1

# 2. The pipeline, test profile (same inputs and settings).
echo "== nextflow"
pixi run -e nf nextflow run . -profile test,pixi \
  --outdir "${OUT}/nf" -work-dir "${OUT}/work" > "${OUT}/nextflow.log" 2>&1

# 3. Compare.
echo "== compare"
(cd "${OUT}/run_all/chr3" && find . -type f | sort) > "${OUT}/files_run_all.txt"
(cd "${OUT}/nf/chr3" && find . -type f | sort) > "${OUT}/files_nf.txt"
if ! diff "${OUT}/files_run_all.txt" "${OUT}/files_nf.txt"; then
  echo "RESULT: FAIL (file lists differ)"; exit 1
fi
echo "files     : $(wc -l < "${OUT}/files_nf.txt") in both"
if PIXI_BIN="${BLAT_BIN}" bash validation/blat_memory/compare_outputs.sh \
     "${OUT}/run_all/chr3" "${OUT}/nf/chr3"; then
  echo "RESULT: PASS"
else
  echo "RESULT: FAIL (contents differ)"; exit 1
fi
