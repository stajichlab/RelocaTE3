#!/usr/bin/bash -l
#SBATCH -p intel
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=06:00:00
#SBATCH -J rt3-compare
#SBATCH -o logs/blat_memory/compare_outputs.%j.log
#
# Compare every output file of two RelocaTE3 output directories.
# Text files: byte-for-byte. BAMs: record-for-record (`samtools view`) and their
# @SQ header lines (`samtools view -H`), since @PG lines legitimately differ.
# A BAM samtools cannot read (either side) is UNREADABLE, a failure: its empty
# output must never compare equal. Each stream is checksummed (sha256) with
# pipefail, so nothing is held in memory or on disk.
# Indexes (.bai/.fai), the measurement files written by the memcheck scripts,
# and symlinked inputs are skipped.
#
# Usage: compare_outputs.sh <DIR_A> <DIR_B>   (or sbatch for large outputs)
set -euo pipefail
A="${1:?first output directory}"
B="${2:?second output directory}"
PIXI_BIN="${PIXI_BIN:-/bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/RelocaTE3_jason/RelocaTE3/results/release-validation/v0.2.0rc1/clone/.pixi/envs/blat/bin}"
export PATH="${PIXI_BIN}:${PATH}"
command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH" >&2; exit 1; }
echo "date: $(date -Is)  A=${A}  B=${B}"

list_files() {
  (cd "$1" && find . -type f ! -name '*.bai' ! -name '*.fai' ! -name 'time.txt' \
     ! -name 'tree_rss*' ! -name '.complete' | sort)
}
# sha256 of a BAM's records (`bam_sum <bam>`) or @SQ lines (`bam_sum <bam> -H`);
# non-zero when samtools fails (pipefail). awk, unlike grep, exits 0 on no match.
bam_sum() {
  local out
  if [[ "${2:-}" == -H ]]; then
    out="$(samtools view -H "$1" | awk '/^@SQ/' | sha256sum)" || return 1
  else
    out="$(samtools view "$1" | sha256sum)" || return 1
  fi
  echo "${out%% *}"
}
same=0; differ=0
while read -r f; do
  if [[ ! -f "${B}/${f}" ]]; then echo "MISSING in B: ${f}"; differ=$((differ + 1)); continue; fi
  case "${f}" in
    *.bam)
      unreadable=""
      ra="$(bam_sum "${A}/${f}")" && ha="$(bam_sum "${A}/${f}" -H)" || unreadable+=" A"
      rb="$(bam_sum "${B}/${f}")" && hb="$(bam_sum "${B}/${f}" -H)" || unreadable+=" B"
      if [[ -n "${unreadable}" ]]; then
        echo "UNREADABLE (samtools view failed in${unreadable}) ${f}"; differ=$((differ + 1))
      elif [[ "${ha}" != "${hb}" ]]; then
        echo "DIFF (@SQ header) ${f}"; differ=$((differ + 1))
      elif [[ "${ra}" != "${rb}" ]]; then
        echo "DIFF (records) ${f}"; differ=$((differ + 1))
      else
        same=$((same + 1))
      fi ;;
    *)
      if cmp -s "${A}/${f}" "${B}/${f}"; then
        same=$((same + 1)); else echo "DIFF ${f}"; differ=$((differ + 1)); fi ;;
  esac
done < <(list_files "${A}")
while read -r f; do [[ -f "${A}/${f}" ]] || { echo "EXTRA in B: ${f}"; differ=$((differ + 1)); }; done < <(list_files "${B}")
echo "identical=${same} differing=${differ}"
[[ "${differ}" -eq 0 ]]
