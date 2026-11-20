#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=04:00:00
#SBATCH -J rt3-lr-l3-gate
#SBATCH -o logs/longread_l3_gate.%j.log
#
# Long-read L3 gate: call recall and precision of lr-align -> lr-signatures ->
# lr-cluster -> lr-annotate -> lr-call on the small simulated Chr3 long-read
# panel.
#
# For each panel_manifest.tsv row, the sample reads and the reference-only
# control reads are each run through the five stages into
# $OUT/<platform>_<sample>[_control]/, then score_calls.py scores each sample's
# calls with its control into $OUT/gate.tsv, plus a per-family table
# $OUT/<platform>_<sample>.per_family.tsv. A match is a PASS spanning call
# within MAX_DIST bp of a truth site naming its family (greedy one-to-one).
# GATE: PASS when recall >= MIN_RECALL and precision >= MIN_PRECISION for the
# 20x row of every platform (ont-hq and hifi); otherwise GATE: FAIL. The script
# exits 0 either way: read the verdict from the log.
#
# Usage, from the repo root (see validation/longread_l3/README.md):
#   mkdir -p logs && sbatch validation/longread_l3/run_gate.slurm
# Overridable: PANEL, GENOME, TE_LIBRARY, OUT, MAX_DIST, MIN_RECALL,
# MIN_PRECISION (environment variables).
#
# Rerun-safe, per stage: each stage is skipped when its output (<set>.lr.bam +
# .bai, <set>.signatures.tsv.gz, <set>.clusters.tsv.gz, <set>.te_hits.tsv.gz,
# <set>.calls.tsv.gz) exists. A BAM without its .bai is refused, and so is a
# skipped stage's output that is older than an input it was made from (an
# earlier stage was rebuilt): remove the set directory and resubmit.

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

THREADS="${SLURM_CPUS_PER_TASK:-1}"
PANEL="${PANEL:-tests/data/longread_small}"
GENOME="${GENOME:-tests/data/sim_genome/MSU7.Chr3_2M.fa}"
TE_LIBRARY="${TE_LIBRARY:-$PANEL/config/te_library.fa}"
OUT="${OUT:-validation/longread_l3/output}"
MAX_DIST="${MAX_DIST:-50}"
MIN_RECALL="${MIN_RECALL:-0.95}"
MIN_PRECISION="${MIN_PRECISION:-0.95}"
GATE_COVERAGE=20
GATE_PLATFORMS=(ont-hq hifi)
MANIFEST="$PANEL/panel_manifest.tsv"
TRUTH="$PANEL/truth_events.tsv"
SCORER="validation/longread_l3/score_calls.py"
GATE_TSV="$OUT/gate.tsv"

echo "[$(date '+%F %T')] RelocaTE3 long-read L3 gate"

# --- checks -----------------------------------------------------------------
command -v pixi >/dev/null 2>&1 || { echo "ERROR: pixi not on PATH" >&2; exit 127; }
for f in "$MANIFEST" "$TRUTH" "$GENOME" "$TE_LIBRARY" "$SCORER"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing input: $f" >&2
    exit 1
  fi
done

COMMIT="$(git rev-parse HEAD)"
DIRTY="$(git status --porcelain | wc -l)"

echo "  host          : $(hostname)"
echo "  job id        : ${SLURM_JOB_ID:-none}"
echo "  repo          : $(pwd)"
echo "  commit        : $COMMIT (dirty files: $DIRTY)"
echo "  panel         : $PANEL -> $(readlink -f "$PANEL")"
echo "  genome        : $GENOME"
echo "  te library    : $TE_LIBRARY"
echo "  out           : $OUT"
echo "  threads       : $THREADS"
echo "  max dist      : $MAX_DIST bp"
echo "  min recall    : $MIN_RECALL at ${GATE_COVERAGE}x for ${GATE_PLATFORMS[*]}"
echo "  min precision : $MIN_PRECISION at ${GATE_COVERAGE}x for ${GATE_PLATFORMS[*]}"
echo "  version       : $(pixi run --frozen relocaTE3 --version 2>/dev/null || echo unknown)"

mkdir -p "$OUT"

# --- per read set -----------------------------------------------------------
# require_fresh <set dir> <output> <input>...: refuse a skipped stage's output
# that is older than an input it was made from (that input was rebuilt since).
require_fresh() {
  local dir="$1" out="$2" input
  shift 2
  for input in "$@"; do
    if [[ "$out" -ot "$input" ]]; then
      echo "ERROR: $out is older than its input $input (an earlier stage was" \
        "rebuilt); remove $dir and resubmit" >&2
      exit 1
    fi
  done
}

# run_set <set> <platform> <reads, relative to PANEL>
run_set() {
  local name="$1" platform="$2" reads="$PANEL/$3"
  local dir="$OUT/$name"
  local bam="$dir/$name.lr.bam"
  local sigs="$dir/$name.signatures.tsv.gz"
  local clusters="$dir/$name.clusters.tsv.gz"
  local hits="$dir/$name.te_hits.tsv.gz"
  local calls="$dir/$name.calls.tsv.gz"

  # Each stage refuses to overwrite its output and writes it atomically, so
  # an existing output means that stage finished: skip it.
  if [[ -f "$bam" && ! -f "$bam.bai" ]]; then
    echo "ERROR: $bam exists without $bam.bai (interrupted lr-align);" \
      "remove $dir and resubmit" >&2
    exit 1
  fi
  if [[ -f "$bam" ]]; then
    echo "[$(date '+%F %T')] $name: lr-align skipped ($bam exists)"
  else
    if [[ ! -f "$reads" ]]; then
      echo "ERROR: missing reads for $name: $reads" >&2
      exit 1
    fi
    echo "[$(date '+%F %T')] $name: lr-align ($platform, $reads)"
    pixi run --frozen relocaTE3 lr-align --sample "$name" -o "$dir" \
      --genome "$GENOME" --read-type "$platform" --reads "$reads" \
      --threads "$THREADS"
  fi
  if [[ -f "$sigs" ]]; then
    require_fresh "$dir" "$sigs" "$bam"
    echo "[$(date '+%F %T')] $name: lr-signatures skipped ($sigs exists)"
  else
    echo "[$(date '+%F %T')] $name: lr-signatures"
    pixi run --frozen relocaTE3 lr-signatures --sample "$name" -o "$dir" \
      --bam "$bam" --read-type "$platform" --genome "$GENOME"
  fi
  if [[ -f "$clusters" ]]; then
    require_fresh "$dir" "$clusters" "$sigs"
    echo "[$(date '+%F %T')] $name: lr-cluster skipped ($clusters exists)"
  else
    echo "[$(date '+%F %T')] $name: lr-cluster"
    pixi run --frozen relocaTE3 lr-cluster --sample "$name" -o "$dir" \
      --signatures "$sigs"
  fi
  if [[ -f "$hits" ]]; then
    require_fresh "$dir" "$hits" "$sigs"
    echo "[$(date '+%F %T')] $name: lr-annotate skipped ($hits exists)"
  else
    echo "[$(date '+%F %T')] $name: lr-annotate"
    pixi run --frozen relocaTE3 lr-annotate --sample "$name" -o "$dir" \
      --signatures "$sigs" --te-library "$TE_LIBRARY" --threads "$THREADS"
  fi
  if [[ -f "$calls" ]]; then
    require_fresh "$dir" "$calls" "$clusters" "$sigs" "$hits"
    echo "[$(date '+%F %T')] $name: lr-call skipped ($calls exists)"
  else
    echo "[$(date '+%F %T')] $name: lr-call"
    pixi run --frozen relocaTE3 lr-call --sample "$name" -o "$dir" \
      --clusters "$clusters" --signatures "$sigs" --te-hits "$hits"
  fi
  echo "[$(date '+%F %T')] $name: done"
}

# Manifest columns: sample platform coverage replicate reads control_reads.
# Read on fd 3 so the commands in the loop cannot consume the manifest from
# stdin; "|| [[ -n $sample ]]" keeps a last line without a trailing newline.
ROWS=()
while IFS=$'\t' read -r -u 3 sample platform coverage replicate reads control_reads \
  || [[ -n "$sample" ]]; do
  [[ "$sample" == "sample" || -z "$sample" ]] && continue
  ROWS+=("$sample"$'\t'"$platform"$'\t'"$coverage")
  run_set "${platform}_${sample}" "$platform" "$reads"
  run_set "${platform}_${sample}_control" "$platform" "$control_reads"
done 3<"$MANIFEST"

if [[ ${#ROWS[@]} -eq 0 ]]; then
  echo "ERROR: no rows in $MANIFEST" >&2
  exit 1
fi

# --- scoring ----------------------------------------------------------------
# gate.tsv and each per-family table are written to .tmp and moved into place.
TMP="$GATE_TSV.tmp"
rm -f "$TMP"
HEADER=(--header)
for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r sample platform coverage <<<"$row"
  name="${platform}_${sample}"
  per_family="$OUT/$name.per_family.tsv"
  rm -f "$per_family.tmp"
  pixi run --frozen python "$SCORER" --truth "$TRUTH" \
    --calls "$OUT/$name/$name.calls.tsv.gz" \
    --control "$OUT/${name}_control/${name}_control.calls.tsv.gz" \
    --max-dist "$MAX_DIST" --label "$name" "${HEADER[@]}" \
    --per-family "$per_family.tmp" >>"$TMP"
  mv "$per_family.tmp" "$per_family"
  HEADER=()
done
mv "$TMP" "$GATE_TSV"

echo
echo "[$(date '+%F %T')] $GATE_TSV"
column -t -s $'\t' "$GATE_TSV"

for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r sample platform coverage <<<"$row"
  name="${platform}_${sample}"
  echo
  echo "Per-family ($OUT/$name.per_family.tsv):"
  column -t -s $'\t' "$OUT/$name.per_family.tsv"
done

# --- timing -----------------------------------------------------------------
echo
echo "Per-stage timing (from each set's timing.tsv; last ok row per stage):"
{
  printf 'set\tstage\twall_s\tmax_rss_mb\trecords_in\trecords_out\tstatus\n'
  for row in "${ROWS[@]}"; do
    IFS=$'\t' read -r sample platform coverage <<<"$row"
    for name in "${platform}_${sample}" "${platform}_${sample}_control"; do
      timing="$OUT/$name/timing.tsv"
      [[ -f "$timing" ]] || continue
      awk -F'\t' -v OFS='\t' -v set="$name" '
        NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        $col["status"] != "ok" { next }  # e.g. FileExistsError rows from reruns
        { last[$col["stage"]] = set OFS $col["stage"] OFS $col["wall_s"] OFS \
            $col["max_rss_mb"] OFS $col["records_in"] OFS $col["records_out"] \
            OFS $col["status"]; order[$col["stage"]] = NR }
        END {
          n = 0
          for (s in order) keys[++n] = s
          for (i = 1; i <= n; i++)
            for (j = i + 1; j <= n; j++)
              if (order[keys[j]] < order[keys[i]]) { t = keys[i]; keys[i] = keys[j]; keys[j] = t }
          for (i = 1; i <= n; i++) print last[keys[i]]
        }' "$timing"
    done
  done
} | column -t -s $'\t'

# --- verdict ----------------------------------------------------------------
# Every gate platform needs a row at GATE_COVERAGE with recall >= MIN_RECALL
# and precision >= MIN_PRECISION.
# gate_value <label> <column>: that column of the label's gate.tsv row.
gate_value() {
  awk -F'\t' -v set="$1" -v c="$2" '
    NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
    $col["label"] == set { print $col[c] }' "$GATE_TSV"
}

VERDICT="PASS"
for platform in "${GATE_PLATFORMS[@]}"; do
  n_rows=0
  for row in "${ROWS[@]}"; do
    IFS=$'\t' read -r sample row_platform coverage <<<"$row"
    [[ "$row_platform" == "$platform" && "$coverage" == "$GATE_COVERAGE" ]] || continue
    n_rows=$((n_rows + 1))
    name="${platform}_${sample}"
    for check in "recall:$MIN_RECALL" "precision:$MIN_PRECISION"; do
      metric="${check%%:*}"
      min="${check#*:}"
      value="$(gate_value "$name" "$metric")"
      if awk -v v="$value" -v m="$min" 'BEGIN { exit !(v != "" && v + 0 >= m + 0) }'; then
        echo "  $name: $metric $value >= $min"
      else
        echo "  $name: $metric ${value:-missing} < $min"
        VERDICT="FAIL"
      fi
    done
  done
  if [[ $n_rows -eq 0 ]]; then
    echo "  $platform: no ${GATE_COVERAGE}x row in $MANIFEST"
    VERDICT="FAIL"
  fi
done

echo
echo "GATE: $VERDICT"
echo "[$(date '+%F %T')] finished"
exit 0
