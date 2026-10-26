#!/usr/bin/bash
# Submit the real-rice regression check for a release candidate.
#
# Runs validation/real_rice/ from the fresh clone made by
# fresh_clone_check.slurm, with that clone's locked `blat` env, so the code and
# tools are exactly what the tag ships. Two jobs:
#   1. SLURM array, one task per sample (run_relocate3.sh).
#   2. After the array succeeds: normalize + compare against the frozen
#      RelocaTE2 calls (run_validation.sh --skip-run).
#
# PATH is set inside --wrap, not with --export: run_relocate3.sh is a login
# script (bash -l), and a login shell resets an exported PATH, so the job would
# fall back to the default pixi env, which has no blat.
#
# Usage, from the RelocaTE3 checkout that ran fresh_clone_check.slurm:
#   bash validation/release/submit_real_rice.sh v0.2.0rc1 [--after <jobid>]
#
# Outputs: results/release-validation/<REF>/clone/validation/real_rice/
#          {results,report,logs}/ ; submission record in
#          results/release-validation/<REF>/real_rice.submission

set -euo pipefail

REF="${1:?usage: submit_real_rice.sh <tag-or-commit> [--after <jobid>]}"
shift
AFTER=""
if [[ "${1:-}" == "--after" ]]; then AFTER="${2:?--after needs a job id}"; fi

BASE_DIR="$(pwd)"
OUTDIR="$BASE_DIR/results/release-validation/${REF}"
CLONE="$OUTDIR/clone"
RR="$CLONE/validation/real_rice"
CONFIG="$RR/config.toml"
ENV_BIN="$CLONE/.pixi/envs/blat/bin"
PLOT_BIN="$BASE_DIR/.pixi/envs/default/bin"   # report-only: plotting stack
RECORD="$OUTDIR/real_rice.submission"

[[ -d "$CLONE/.git" ]] || { echo "ERROR: no fresh clone at $CLONE (run fresh_clone_check.slurm first)" >&2; exit 1; }
[[ -z "$AFTER" && ! -x "$ENV_BIN/relocaTE3" ]] && { echo "ERROR: $ENV_BIN/relocaTE3 missing; pass --after <release-check job>" >&2; exit 1; }
[[ -x "$PLOT_BIN/python3" ]] || { echo "ERROR: plotting env missing: $PLOT_BIN" >&2; exit 1; }
[[ -e "$RECORD" ]] && { echo "ERROR: $RECORD exists; real-rice check already submitted" >&2; exit 1; }
[[ -e "$CONFIG" ]] && { echo "ERROR: $CONFIG exists; move it aside" >&2; exit 1; }

cp "$RR/config.example.toml" "$CONFIG"
SAMPLE_CSV="$(cd "$RR" && python3 -c "from _config import load_config; print(load_config('$CONFIG')['paths']['sample_csv'])")"
N=$(( $(wc -l < "$SAMPLE_CSV") - 1 ))
mkdir -p "$RR/logs"

DEP=()
[[ -n "$AFTER" ]] && DEP=(--dependency="afterok:${AFTER}")

ARRAY_ID=$(sbatch --parsable "${DEP[@]}" \
  -p epyc --cpus-per-task 16 --mem 64gb --time 08:00:00 \
  --array=1-"$N" -J "rt3-rice-${REF}" \
  -o "$RR/logs/relocate3.%A_%a.log" \
  --wrap "export PATH='$ENV_BIN':\$PATH; bash '$RR/run_relocate3.sh' '$CONFIG'")

REPORT_ID=$(sbatch --parsable --dependency="afterok:${ARRAY_ID}" \
  -p short --cpus-per-task 1 --mem 8gb --time 01:00:00 -J "rt3-rice-report-${REF}" \
  -o "$RR/logs/report.%j.log" \
  --wrap "export PATH='$PLOT_BIN':\$PATH; cd '$CLONE' && bash validation/real_rice/run_validation.sh --skip-run --config '$CONFIG'")

{
  echo "submitted : $(date '+%F %T')"
  echo "ref       : $REF ($(git -C "$CLONE" rev-parse HEAD))"
  echo "samples   : $N ($SAMPLE_CSV)"
  echo "after     : ${AFTER:-none}"
  echo "array     : $ARRAY_ID"
  echo "report    : $REPORT_ID"
} | tee "$RECORD"
