#!/usr/bin/bash -l
#SBATCH -p intel
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH -J rt3-nf-rice-compare
#SBATCH -o logs/real_rice_compare.%j.log
#
# Compare the pipeline's real-rice output (run_pipeline.sh) with the v0.2.0rc2
# run-all output, sample by sample.
#
# Usage, from the repository root (sbatch, or bash directly: the rc2 outputs
# are 27-82 MB per sample, about 0.6 GB in all):
#   sbatch validation/real_rice/compare_pipeline_to_rc2.sh PIPELINE_OUTDIR [RC2_RESULTS]
#
#   PIPELINE_OUTDIR  run_pipeline.sh's OUTDIR; its samplesheet.csv lists the samples
#   RC2_RESULTS      default: the v0.2.0rc2 release-validation run (RC2_DEFAULT below)
#   PIXI_BIN         directory holding samtools, default .pixi/envs/blat/bin
#
# For each sample:
#   1. File lists (`find . -type f | sort`, indexes included) must be equal;
#      every file only on one side is printed.
#   2. validation/blat_memory/compare_outputs.sh RC2/<s> PIPELINE/<s>: text files
#      byte for byte, BAMs record by record (`samtools view`; @PG header lines
#      differ legitimately). It does not compare the contents of .bai/.fai
#      indexes (nor of time.txt, tree_rss*, .complete, which rc2 does not have),
#      so step 1 is what checks those: they are counted as presence_only.
# Summary line per sample: <s> files=N identical=M differing=K presence_only=P,
# N = files in rc2. Then RESULT: PASS (every sample equal) or RESULT: FAIL;
# the exit status is non-zero on FAIL. Nothing is excluded silently: a file
# missing, extra or different on either side fails its sample. Samples in
# RC2_RESULTS but not in the sample sheet fail too.
#
# Expected per-sample file set (checked 2026-11-03: all 10 rc2 samples hold
# these 31 files, <s> = sample):
#   <s>.left.bam  <s>.left.bam.bai  <s>.right.bam  <s>.right.bam.bai
#   <s>.repeat.bwaaln.sorted.bam  <s>.repeat.bwaaln.sorted.bam.bai
#   existingTE.bed
#   flanking/<s>.left.flankingReads.fq  flanking/<s>.right.flankingReads.fq
#   genome_aln/<s>.fullreads.genome.bam  genome_aln/<s>.fullreads.genome.bam.bai
#   results/ALL.mping.all_nonref_insert.{all,characTErized,high_conf,raw}.{gff,txt}
#   results/ALL.mping.all_nonref_insert.{gff,txt}
#   results/ALL.mping.all_nonref_supporting.{gff,txt}
#   results/<s>.all_ref_insert.{gff,txt}
#   te_containing/<s>.{left,right}.ContainingReads.fq
#   te_containing/<s>.read_repeat_name.txt  te_containing/<s>.te_hit_names.txt
#   te_portions/<s>.five_prime.fa  te_portions/<s>.three_prime.fa
# No genome_aln/<s>.reads.genome.bam on either side: rc2 ran run-all with
# --genotype-bam <CRAM>, so run-all aligned no reads for genotyping, and the
# pipeline skips ALIGN_READS for samples with a genotype_bam. If one appears,
# the comparison fails (only in pipeline).
# Not compared: PIPELINE_OUTDIR's top level (pipeline_info/, samplesheet.csv,
# container_tag.txt), which run-all has no counterpart for.
set -euo pipefail

RC2_DEFAULT="/bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/RelocaTE3_jason/RelocaTE3/results/release-validation/v0.2.0rc2/clone/validation/real_rice/results"

usage() { sed -n '/^# Usage/,/^#   PIXI_BIN/{s/^# \{0,1\}//;p;}' "${BASH_SOURCE[0]}"; }
case "${1:-}" in
  -h|--help) usage; exit 0 ;;
  "") echo "ERROR: missing PIPELINE_OUTDIR" >&2; usage >&2; exit 2 ;;
esac
[[ $# -le 2 ]] || { echo "ERROR: too many arguments" >&2; usage >&2; exit 2; }
NF="$(realpath -m "$1")"
RC2="$(realpath -m "${2:-${RC2_DEFAULT}}")"
REPO="$(pwd)"
COMPARE_OUTPUTS="${REPO}/validation/blat_memory/compare_outputs.sh"
export PIXI_BIN="${PIXI_BIN:-${REPO}/.pixi/envs/blat/bin}"
SHEET="${NF}/samplesheet.csv"

echo "date      : $(date -Is)"
echo "job       : ${SLURM_JOB_ID:-none}"
echo "pipeline  : ${NF}"
echo "rc2       : ${RC2}"
echo "samtools  : ${PIXI_BIN}/samtools"

[[ -f "${COMPARE_OUTPUTS}" ]] || { echo "ERROR: run from the repository root (no ${COMPARE_OUTPUTS})" >&2; exit 1; }
[[ -x "${PIXI_BIN}/samtools" ]] || { echo "ERROR: no samtools in PIXI_BIN=${PIXI_BIN}" >&2; exit 1; }
[[ -d "${RC2}" ]] || { echo "ERROR: rc2 results not found: ${RC2}" >&2; exit 1; }
[[ -f "${SHEET}" ]] || { echo "ERROR: ${SHEET} not found (written by run_pipeline.sh)" >&2; exit 1; }

# Sample names, first-seen order (lanes of one sample share a name).
mapfile -t SAMPLES < <(tail -n +2 "${SHEET}" | cut -d, -f1 | awk 'NF && !seen[$0]++')
echo "samples   : ${#SAMPLES[@]} (${SAMPLES[*]})"

list_files() { (cd "$1" && find . -type f | sort); }

fail=0
SUMMARY=()
# rc2 samples absent from the sheet.
while read -r d; do
  s="$(basename "${d}")"
  if ! printf '%s\n' "${SAMPLES[@]}" | grep -qxF "${s}"; then
    SUMMARY+=("${s} missing: in rc2 but not in ${SHEET}"); fail=1
  fi
done < <(find "${RC2}" -mindepth 1 -maxdepth 1 -type d | sort)

TMP="$(mktemp -d)"
trap 'rm -rf "${TMP}"' EXIT
for s in "${SAMPLES[@]}"; do
  echo "== ${s}"
  if [[ ! -d "${NF}/${s}" || ! -d "${RC2}/${s}" ]]; then
    [[ -d "${NF}/${s}" ]] || SUMMARY+=("${s} missing: no pipeline output ${NF}/${s}")
    [[ -d "${RC2}/${s}" ]] || SUMMARY+=("${s} missing: no rc2 output ${RC2}/${s}")
    fail=1; continue
  fi
  list_files "${RC2}/${s}" > "${TMP}/rc2.txt"
  list_files "${NF}/${s}" > "${TMP}/nf.txt"
  lists_equal=1
  while read -r f; do echo "  only in rc2: ${f}"; lists_equal=0; done < <(comm -23 "${TMP}/rc2.txt" "${TMP}/nf.txt")
  while read -r f; do echo "  only in pipeline: ${f}"; lists_equal=0; done < <(comm -13 "${TMP}/rc2.txt" "${TMP}/nf.txt")

  set +e
  bash "${COMPARE_OUTPUTS}" "${RC2}/${s}" "${NF}/${s}" > "${TMP}/cmp.txt" 2>&1
  cmp_status=$?
  set -e
  sed 's/^/  /' "${TMP}/cmp.txt"
  counts="$(grep -E '^identical=[0-9]+ differing=[0-9]+$' "${TMP}/cmp.txt" | tail -n 1 || true)"
  if [[ -z "${counts}" ]]; then
    SUMMARY+=("${s} error: compare_outputs.sh exited ${cmp_status} without counts"); fail=1; continue
  fi
  identical="${counts#identical=}"; identical="${identical%% *}"
  differing="${counts##*differing=}"
  files="$(wc -l < "${TMP}/rc2.txt")"
  presence=$(( files - identical - differing ))
  (( presence >= 0 )) || presence=0
  line="${s} files=${files} identical=${identical} differing=${differing} presence_only=${presence}"
  if [[ "${lists_equal}" -ne 1 ]]; then line+=" file_lists=DIFFER"; fail=1; fi
  if [[ "${cmp_status}" -ne 0 || "${differing}" -ne 0 ]]; then fail=1; fi
  SUMMARY+=("${line}")
done

echo "== summary"
printf '%s\n' "${SUMMARY[@]}"
if [[ "${fail}" -eq 0 ]]; then
  echo "RESULT: PASS"
else
  echo "RESULT: FAIL"; exit 1
fi
