#!/usr/bin/bash -l
#SBATCH -p intel
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=06:00:00
#SBATCH -J rt3-compare
#SBATCH -o logs/blat_memory/compare_outputs.%j.log
#
# Compare every output file of two RelocaTE3 output directories.
# Text files: byte-for-byte. BAMs: record-for-record (`samtools view`), since
# @PG header lines legitimately differ. Indexes (.bai/.fai), the measurement
# files written by the memcheck scripts, and symlinked inputs are skipped.
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
same=0; differ=0
while read -r f; do
  if [[ ! -f "${B}/${f}" ]]; then echo "MISSING in B: ${f}"; differ=$((differ + 1)); continue; fi
  case "${f}" in
    *.bam)
      if cmp -s <(samtools view "${A}/${f}") <(samtools view "${B}/${f}"); then
        same=$((same + 1)); else echo "DIFF (records) ${f}"; differ=$((differ + 1)); fi ;;
    *)
      if cmp -s "${A}/${f}" "${B}/${f}"; then
        same=$((same + 1)); else echo "DIFF ${f}"; differ=$((differ + 1)); fi ;;
  esac
done < <(list_files "${A}")
while read -r f; do [[ -f "${A}/${f}" ]] || { echo "EXTRA in B: ${f}"; differ=$((differ + 1)); }; done < <(list_files "${B}")
echo "identical=${same} differing=${differ}"
[[ "${differ}" -eq 0 ]]
