#!/usr/bin/bash -l
#SBATCH -p intel
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH -J rt3-cohort-report
#SBATCH -o logs/real_rice_cohort_report.%j.log
#
# Cohort table, genotype matrix and MultiQC report for a finished real-rice
# pipeline output directory, using the same commands as the pipeline's
# RELOCATE3_SAMPLE_STATS, RELOCATE3_COHORT and MULTIQC processes. Then checks
# that the cohort table holds every sample's characterized rows unchanged.
#
# Sample names are the pipeline's ([A-Za-z0-9._-], so never '='); they are
# passed to `relocaTE3 cohort` as NAME=PATH.
#
# Usage, from the repository root:
#   sbatch validation/real_rice/cohort_report.sh <pipeline_outdir> [outdir]
#   (default outdir: results/nextflow-m3/real-rice-report)
set -euo pipefail

M2="$(realpath "${1:?usage: cohort_report.sh <pipeline_outdir> [outdir]}")"
OUT="${2:-$(pwd)/results/nextflow-m3/real-rice-report}"
E="${PIXI_ENV:-$(pwd)/.pixi/envs/pipeline}/bin"
TE_NAME="${TE_NAME:-mping}"
TARGET="${TARGET:-ALL}"
export PYTHONNOUSERSITE=1

echo "date     : $(date -Is)"
echo "input    : ${M2}"
echo "outdir   : ${OUT}"
echo "env bin  : ${E}"
echo "te_name  : ${TE_NAME}   target: ${TARGET}"
[[ -x "${E}/relocaTE3" && -x "${E}/multiqc" ]] || { echo "ERROR: relocaTE3/multiqc missing in ${E}" >&2; exit 1; }
[[ -e "${OUT}" ]] && { echo "ERROR: ${OUT} exists; refusing to overwrite" >&2; exit 1; }
mkdir -p "${OUT}/mqc_in"

INPUTS=()
for d in "${M2}"/*/results/"${TARGET}.${TE_NAME}".all_nonref_insert.characTErized.txt; do
  s=$(basename "$(dirname "$(dirname "$d")")")
  "${E}/relocaTE3" sample-stats -d "${M2}/${s}" -n "${s}" --target "${TARGET}" --te-name "${TE_NAME}" \
    -o "${OUT}/mqc_in/${s}.stats_mqc.tsv"
  INPUTS+=("${s}=${d}")
done
echo "samples  : ${#INPUTS[@]}"
"${E}/relocaTE3" cohort -i "${INPUTS[@]}" -o "${OUT}/cohort"
cp "${OUT}/cohort/cohort.summary_mqc.tsv" "${OUT}/mqc_in/"
"${E}/multiqc" --no-ansi --config "$(pwd)/assets/multiqc_config.yml" \
  --filename relocate3_report.html -o "${OUT}/multiqc" "${OUT}/mqc_in"

# Check: cohort rows == concatenation of per-sample data rows, in sample order
# (byte order of the sample name, as relocaTE3 cohort sorts; G_1 < G_104).
expected=$(mktemp -p "${OUT}" .expected.XXXXXX); trap 'rm -f "${expected}"' EXIT
printf '%s\n' "${INPUTS[@]}" | LC_ALL=C sort -t= -k1,1 | while IFS= read -r pair; do
  tail -n +2 "${pair#*=}"
done > "${expected}"
if cmp -s "${expected}" <(tail -n +2 "${OUT}/cohort/cohort.characTErized.txt"); then
  echo "cohort rows : $(wc -l < "${expected}") identical to the per-sample tables"
else
  echo "RESULT: FAIL (cohort rows differ from the per-sample tables)"; exit 1
fi
echo "matrix      : $(($(wc -l < "${OUT}/cohort/cohort.genotype_matrix.tsv") - 1)) sites x ${#INPUTS[@]} samples"
awk -F'\t' 'NR>1 {n=0; for (i=6;i<=NF;i++) if ($i!=".") n++; if (n>=2) s++} END {print "shared      : " s+0 " sites called in >=2 samples"}' \
  "${OUT}/cohort/cohort.genotype_matrix.tsv"
ls -l "${OUT}/multiqc/relocate3_report.html"
echo "RESULT: PASS"
