#!/usr/bin/bash
# Submit the real-rice regression check for a release candidate.
#
# Runs validation/real_rice/ from the fresh clone made by
# fresh_clone_check.slurm, so the harness and pipeline code are exactly what
# the tag ships. Two jobs, in one of two modes:
#
#   --mode pipeline (default): the Nextflow pipeline with the release image
#     1. Head job: validation/real_rice/run_pipeline.sh <container tag>
#        (-profile hpcc,apptainer; one SLURM job per task) into
#        clone/validation/real_rice/results/.
#     2. After it succeeds: normalize + compare against the frozen RelocaTE2
#        calls (run_validation.sh --mode pipeline --skip-run on that outdir).
#   --mode run-all: relocaTE3 run-all from the clone's locked `blat` env
#     1. SLURM array, one task per sample (run_relocate3.sh).
#     2. After the array succeeds: run_validation.sh --mode run-all --skip-run.
#
# The container tag (pipeline mode, required): the image of the release, e.g.
# 0.3.0rc1 once container.yml has pushed it for the tag, or before tagging the
# release commit's sha-<short> image. Mutable tags (dev, pr-<n>) are refused.
#
# run-all: PATH is set inside --wrap, not with --export: run_relocate3.sh is a
# login script (bash -l), and a login shell resets an exported PATH, so the job
# would fall back to the default pixi env, which has no blat. (run_pipeline.sh
# takes its Python from PYTHON, which a login shell keeps.)
#
# Usage, from the RelocaTE3 checkout that ran fresh_clone_check.slurm:
#   bash validation/release/submit_real_rice.sh <tag-or-commit> --container-tag <tag> [--after <jobid>]
#   bash validation/release/submit_real_rice.sh <tag-or-commit> --mode run-all [--after <jobid>] [--partition <p>]
# e.g. submit_real_rice.sh v0.3.0rc1 --container-tag 0.3.0rc1
# --after waits for the release-check job; --partition sets the run-all
# array's partition (default: epyc; the pipeline head job uses
# run_pipeline.sh's own #SBATCH settings).
#
# Environment:
#   NXF_APPTAINER_CACHEDIR  pipeline image cache, default
#                           <this checkout>/results/apptainer-cache
#   DRY_RUN=1               print the sbatch commands ("+ ..."); submit and write nothing
#
# Outputs: results/release-validation/<REF>/clone/validation/real_rice/
#          {results,report,logs}/ (pipeline mode also results.nextflow/, the
#          launch and work directory); submission record in
#          results/release-validation/<REF>/real_rice.submission

set -euo pipefail

USAGE="usage: submit_real_rice.sh <tag-or-commit> [--mode pipeline|run-all] [--container-tag <tag>] [--after <jobid>] [--partition <p>]"
REF="${1:?$USAGE}"
shift
AFTER=""
PARTITION="epyc"
MODE="pipeline"
TAG=""
DRY_RUN="${DRY_RUN:-0}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --after) AFTER="${2:?--after needs a job id}"; shift 2 ;;
    --partition) PARTITION="${2:?--partition needs a name}"; shift 2 ;;
    --mode) MODE="${2:?--mode needs pipeline or run-all}"; shift 2 ;;
    --container-tag) TAG="${2:?--container-tag needs a tag}"; shift 2 ;;
    *) echo "ERROR: unknown argument: $1" >&2; echo "$USAGE" >&2; exit 2 ;;
  esac
done
case "$MODE" in
  pipeline)
    [[ -n "$TAG" ]] || { echo "ERROR: pipeline mode needs --container-tag <release image tag>, e.g. 0.3.0rc1 or sha-<short>" >&2; echo "$USAGE" >&2; exit 2; }
    [[ "$TAG" =~ ^[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}$ ]] || { echo "ERROR: invalid image tag '$TAG'" >&2; exit 2; }
    if [[ "$TAG" == dev || "$TAG" =~ ^pr-[0-9]+$ ]]; then
      echo "ERROR: container tag '$TAG' is mutable (CI re-pushes it); a release check needs the release image (e.g. 0.3.0rc1) or sha-<short>" >&2
      exit 2
    fi
    ;;
  run-all) [[ -z "$TAG" ]] || { echo "ERROR: --container-tag is pipeline mode only" >&2; exit 2; } ;;
  *) echo "ERROR: --mode must be pipeline or run-all, not '$MODE'" >&2; exit 2 ;;
esac

BASE_DIR="$(pwd)"
OUTDIR="$BASE_DIR/results/release-validation/${REF}"
CLONE="$OUTDIR/clone"
RR="$CLONE/validation/real_rice"
CONFIG="$RR/config.toml"
ENV_BIN="$CLONE/.pixi/envs/blat/bin"
PLOT_BIN="$BASE_DIR/.pixi/envs/default/bin"   # report-only: plotting stack
RECORD="$OUTDIR/real_rice.submission"
CACHE="${NXF_APPTAINER_CACHEDIR:-$BASE_DIR/results/apptainer-cache}"

[[ -d "$CLONE/.git" ]] || { echo "ERROR: no fresh clone at $CLONE (run fresh_clone_check.slurm first)" >&2; exit 1; }
[[ -z "$AFTER" && ! -x "$ENV_BIN/relocaTE3" ]] && { echo "ERROR: $ENV_BIN/relocaTE3 missing; pass --after <release-check job>" >&2; exit 1; }
[[ -x "$PLOT_BIN/python3" ]] || { echo "ERROR: plotting env missing: $PLOT_BIN" >&2; exit 1; }
[[ -e "$RECORD" ]] && { echo "ERROR: $RECORD exists; real-rice check already submitted" >&2; exit 1; }
[[ -e "$CONFIG" ]] && { echo "ERROR: $CONFIG exists; move it aside" >&2; exit 1; }

# A dry run reads the example config in place and writes nothing.
if [[ "$DRY_RUN" == 1 ]]; then
  READ_CONFIG="$RR/config.example.toml"
else
  cp "$RR/config.example.toml" "$CONFIG"
  READ_CONFIG="$CONFIG"
  mkdir -p "$RR/logs"
fi
SAMPLE_CSV="$(cd "$RR" && python3 -c "from _config import load_config; print(load_config('$READ_CONFIG')['paths']['sample_csv'])")"
N=$(( $(wc -l < "$SAMPLE_CSV") - 1 ))

# submit ARGS...: sbatch --parsable ARGS, or print it (to the script's stdout,
# fd 3) under DRY_RUN=1 and return a placeholder job id.
exec 3>&1
submit() {
  if [[ "$DRY_RUN" == 1 ]]; then
    echo "+ sbatch --parsable $*" >&3
    echo "<job-${#SUBMITTED[@]}>"
  else
    sbatch --parsable "$@"
  fi
}
SUBMITTED=()

DEP=()
[[ -n "$AFTER" ]] && DEP=(--dependency="afterok:${AFTER}")

if [[ "$MODE" == pipeline ]]; then
  RUN_ID=$(submit "${DEP[@]}" --chdir "$CLONE" -J "rt3-rice-nf-${REF}" \
    -o "$RR/logs/real_rice_pipeline.%j.log" \
    --export="ALL,CONFIG=$CONFIG,PYTHON=$ENV_BIN/python,NXF_APPTAINER_CACHEDIR=$CACHE" \
    "$RR/run_pipeline.sh" "$TAG" "$RR/results")
  SUBMITTED+=("$RUN_ID")
  printf -v WRAP 'export PATH=%q:$PATH; cd %q && bash validation/real_rice/run_validation.sh --mode pipeline --skip-run --pipeline-outdir %q --config %q' \
    "$PLOT_BIN" "$CLONE" "$RR/results" "$CONFIG"
else
  RUN_ID=$(submit "${DEP[@]}" \
    -p "$PARTITION" --cpus-per-task 16 --mem 64gb --time 08:00:00 \
    --array=1-"$N" -J "rt3-rice-${REF}" \
    -o "$RR/logs/relocate3.%A_%a.log" \
    --wrap "export PATH='$ENV_BIN':\$PATH; bash '$RR/run_relocate3.sh' '$CONFIG'")
  SUBMITTED+=("$RUN_ID")
  printf -v WRAP 'export PATH=%q:$PATH; cd %q && bash validation/real_rice/run_validation.sh --mode run-all --skip-run --config %q' \
    "$PLOT_BIN" "$CLONE" "$CONFIG"
fi

REPORT_ID=$(submit --dependency="afterok:${RUN_ID}" \
  -p short --cpus-per-task 1 --mem 8gb --time 01:00:00 -J "rt3-rice-report-${REF}" \
  -o "$RR/logs/report.%j.log" \
  --wrap "$WRAP")

RECORD_TEXT="$(
  echo "submitted : $(date '+%F %T')"
  echo "ref       : $REF ($(git -C "$CLONE" rev-parse HEAD))"
  echo "mode      : $MODE"
  [[ "$MODE" == pipeline ]] && echo "image tag : $TAG (cache $CACHE)"
  echo "samples   : $N ($SAMPLE_CSV)"
  echo "after     : ${AFTER:-none}"
  [[ "$MODE" == run-all ]] && echo "partition : $PARTITION"
  echo "run       : $RUN_ID"
  echo "report    : $REPORT_ID"
)"
if [[ "$DRY_RUN" == 1 ]]; then
  echo "$RECORD_TEXT"
  echo "dry run   : nothing submitted, $RECORD not written"
else
  echo "$RECORD_TEXT" | tee "$RECORD"
fi
