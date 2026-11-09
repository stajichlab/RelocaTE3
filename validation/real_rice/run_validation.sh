#!/usr/bin/bash -l
# End-to-end real-rice validation:
#   1. Make sure config.toml exists (copy from config.example.toml if not).
#   2. Run RelocaTE3 on every sample, in one of two modes:
#        --mode pipeline (default): the Nextflow pipeline, run_pipeline.sh
#          (-profile hpcc,apptainer, pinned --container-tag), submitted as a
#          head job with `sbatch --wait`. Outputs: --pipeline-outdir, default
#          results/nextflow-real-rice/<tag> under the repository root.
#        --mode run-all: run_relocate3.sh (relocaTE3 run-all) as a SLURM array
#          with `sbatch --wait` (or sequentially with --local). Outputs: the
#          config's paths.relocate3_outdir.
#   3. Normalize both call sets (RelocaTE2 + RelocaTE3) into report/*.tsv.
#   4. Run the comparison and emit summary + venn.
#   Steps 3-4 are the same commands in both modes; in pipeline mode the
#   RelocaTE3 outdir they read is the pipeline outdir (exported as
#   RELOCATE3_VALIDATION_OUTDIR, which _config.py honours; config untouched).
#   This script owns RELOCATE3_VALIDATION_OUTDIR / _REPORT_DIR: a value
#   inherited from the caller is refused (use --pipeline-outdir / --report-dir).
#
# Designed to be launched from the RelocaTE3/ working directory:
#   bash validation/real_rice/run_validation.sh --container-tag sha-2c9a448
#   bash validation/real_rice/run_validation.sh --skip-run --pipeline-outdir DIR   # only re-report
#   bash validation/real_rice/run_validation.sh --mode run-all [--local] [SAMPLE ...]
#   bash validation/real_rice/run_validation.sh --mode run-all --skip-run
#
# Options:
#   --mode pipeline|run-all   default pipeline
#   --container-tag TAG       pipeline mode: image tag (sha-<short> or a release
#                             version; env CONTAINER_TAG). Required unless
#                             --skip-run with --pipeline-outdir. Mutable tags
#                             (dev, pr-<n>) are refused unless ALLOW_MUTABLE_TAG=1.
#   --pipeline-outdir DIR     pipeline mode: pipeline outdir to run into / report on
#   --report-dir DIR          report here instead of the default: the config's
#                             paths.report_dir (run-all) or <paths.report_dir>/<tag>
#                             (pipeline; <tag> = --container-tag, else the
#                             basename of --pipeline-outdir)
#   --config FILE             default validation/real_rice/config.toml
#   --skip-run --skip-nonref --skip-char
#   --local, --force, SAMPLE...   run-all mode only (pipeline mode resumes into
#                                 its outdir and always runs the whole sheet)
#
# Notes:
#   * Both SLURM paths use `sbatch --wait`, so the script blocks until the
#     RelocaTE3 run finishes, then runs normalize + compare.
#   * Re-runs are safe: pipeline mode resumes (-resume) into the same outdir;
#     run-all mode reuses existing per-sample result directories under
#     paths.relocate3_outdir (--force wipes them first).
#   * DRY_RUN=1 prints the commands ("+ ..."), runs nothing and creates nothing
#     (pipeline mode still runs run_pipeline.sh's own DRY_RUN input check).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../.." && pwd)"

MODE="pipeline"        # pipeline | run-all
EXEC="slurm"           # run-all only: slurm | local
SKIP_RUN=0
SKIP_NONREF=0
SKIP_CHAR=0
FORCE=0
CONFIG_ARG=""
TAG="${CONTAINER_TAG:-}"
PIPE_OUT_ARG=""
REPORT_DIR_ARG=""
EXPLICIT_SAMPLES=()    # positional sample names; empty = use the config's sample CSV
LOCAL_SET=0
DRY_RUN="${DRY_RUN:-0}"

die() { echo "ERROR: $*" >&2; exit 2; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)            MODE="${2:?--mode needs pipeline or run-all}"; shift 2 ;;
    --container-tag)   TAG="${2:?--container-tag needs a tag}"; shift 2 ;;
    --pipeline-outdir) PIPE_OUT_ARG="${2:?--pipeline-outdir needs a directory}"; shift 2 ;;
    --report-dir)      REPORT_DIR_ARG="${2:?--report-dir needs a directory}"; shift 2 ;;
    --local)        EXEC="local"; LOCAL_SET=1; shift ;;
    --slurm)        EXEC="slurm"; shift ;;
    --skip-run)     SKIP_RUN=1;   shift ;;
    --skip-nonref)  SKIP_NONREF=1; shift ;;
    --skip-char)    SKIP_CHAR=1;   shift ;;
    --force)        FORCE=1;      shift ;;
    --config)       CONFIG_ARG="${2:?--config needs a file}"; shift 2 ;;
    -h|--help)
      sed -n '2,49p' "$0"; exit 0 ;;
    --)          shift; while [[ $# -gt 0 ]]; do EXPLICIT_SAMPLES+=("$1"); shift; done ;;
    -*)
      die "unknown flag: $1" ;;
    *)
      EXPLICIT_SAMPLES+=("$1"); shift ;;
  esac
done

# The output overrides _config.py reads are set here from the flags only; an
# inherited value would silently redirect the run (sbatch exports the
# environment) and the report.
for var in RELOCATE3_VALIDATION_OUTDIR RELOCATE3_VALIDATION_REPORT_DIR; do
  if [[ -n "${!var:-}" ]]; then
    die "$var is set in the environment (${!var}); unset it and use --pipeline-outdir (pipeline mode) / --report-dir instead"
  fi
done

case "$MODE" in
  pipeline)
    [[ "$LOCAL_SET" == 0 ]] || die "--local is run-all only (the pipeline submits its tasks to SLURM)"
    [[ "$FORCE" == 0 ]] || die "--force is run-all only (pipeline mode resumes; use a new --pipeline-outdir)"
    if [[ ${#EXPLICIT_SAMPLES[@]} -gt 0 && "$SKIP_RUN" == 0 ]]; then
      die "sample names are run-all only (the pipeline runs the whole sample sheet); with --skip-run they scope the report"
    fi
    if [[ -z "$TAG" && ( "$SKIP_RUN" == 0 || -z "$PIPE_OUT_ARG" ) ]]; then
      die "pipeline mode needs --container-tag TAG (or CONTAINER_TAG), e.g. sha-<short> or a release version; --skip-run --pipeline-outdir DIR needs none"
    fi
    ;;
  run-all)
    [[ -z "$PIPE_OUT_ARG" ]] || die "--pipeline-outdir is pipeline mode only"
    ;;
  *) die "--mode must be pipeline or run-all, not '$MODE'" ;;
esac

# Resolve config path. Default lives next to this script; an explicit --config
# argument can be absolute or relative to the caller's CWD.
if [[ -n "$CONFIG_ARG" ]]; then
  if [[ -f "$CONFIG_ARG" ]]; then
    CONFIG="$(cd "$(dirname "$CONFIG_ARG")" && pwd)/$(basename "$CONFIG_ARG")"
  else
    echo "ERROR: --config file not found: $CONFIG_ARG" >&2; exit 1
  fi
else
  CONFIG="$SCRIPT_DIR/config.toml"
  if [[ ! -f "$CONFIG" ]]; then
    if [[ "$DRY_RUN" == 1 ]]; then
      CONFIG="$SCRIPT_DIR/config.example.toml"
    else
      echo "[validate-rice] config.toml not found; copying from config.example.toml"
      cp "$SCRIPT_DIR/config.example.toml" "$CONFIG"
    fi
  fi
fi

# Overrides of the config's output paths (read by _config.py).
if [[ "$MODE" == pipeline ]]; then
  PIPE_OUT="$(realpath -m "${PIPE_OUT_ARG:-$REPO/results/nextflow-real-rice/$TAG}")"
  export RELOCATE3_VALIDATION_OUTDIR="$PIPE_OUT"
fi
if [[ -n "$REPORT_DIR_ARG" ]]; then
  export RELOCATE3_VALIDATION_REPORT_DIR="$(realpath -m "$REPORT_DIR_ARG")"
fi

echo "[validate-rice] date: $(date -Is)"
echo "[validate-rice] mode: $MODE"
echo "[validate-rice] config: $CONFIG"

# Pull a few values out of the config for orchestration.
mapfile -t CFG < <(
  cd "$SCRIPT_DIR"
  python3 -c "
from _config import load_config, load_samples
cfg = load_config('$CONFIG')
print(cfg['paths']['sample_csv'])
print(cfg['paths']['relocate3_outdir'])
print(cfg['paths']['log_dir'])
print(cfg['slurm']['partition'])
print(cfg['slurm']['cpus_per_task'])
print(cfg['slurm']['mem'])
print(cfg['slurm']['time'])
print(','.join(load_samples(cfg['paths']['sample_csv'])))
print(int(bool(cfg['relocate3'].get('characterize', False))))
print(cfg['relocate3'].get('target', 'ALL'))
print(cfg['relocate3'].get('te_name', 'mping'))
print(cfg['paths']['report_dir'])
"
)
[[ ${#CFG[@]} -eq 12 ]] || { echo "ERROR: could not read $CONFIG" >&2; exit 1; }
SAMPLE_CSV="${CFG[0]}"
OUTROOT="${CFG[1]}"
LOG_DIR="${CFG[2]}"
PART="${CFG[3]}"
CPUS="${CFG[4]}"
MEM="${CFG[5]}"
TIMELIM="${CFG[6]}"
DO_CHAR="${CFG[8]}"
TARGET="${CFG[9]}"
TE_NAME="${CFG[10]}"
REPORT_DIR="${CFG[11]}"
# Pipeline mode reports per tag by default, so reports of different images
# (and the run-all report in paths.report_dir) do not overwrite each other.
if [[ "$MODE" == pipeline && -z "$REPORT_DIR_ARG" ]]; then
  REPORT_DIR="$REPORT_DIR/${TAG:-$(basename "$PIPE_OUT")}"
  export RELOCATE3_VALIDATION_REPORT_DIR="$REPORT_DIR"
fi
echo "[validate-rice] relocate3_outdir: $OUTROOT"
echo "[validate-rice] report_dir: $REPORT_DIR"
if [[ "$MODE" == pipeline && -n "$TAG" ]]; then echo "[validate-rice] container_tag: $TAG"; fi
if [[ "$DRY_RUN" == 1 ]]; then echo "[validate-rice] dry run: commands are printed (+), not run"; fi
# Which samples to operate on: explicit positional args win over the CSV.
if [[ ${#EXPLICIT_SAMPLES[@]} -gt 0 ]]; then
  SAMPLES=("${EXPLICIT_SAMPLES[@]}")
  echo "[validate-rice] using ${#SAMPLES[@]} explicit sample(s): ${SAMPLES[*]}"
else
  IFS=',' read -r -a SAMPLES <<< "${CFG[7]}"
fi

show() {
  local w out=()
  for w in "$@"; do
    if [[ "${w}" =~ ^[A-Za-z0-9_./:=,@%+-]+$ ]]; then out+=("${w}"); else out+=("$(printf '%q' "${w}")"); fi
  done
  echo "+ ${out[*]}"
}
# run CMD...: run it, or print it under DRY_RUN=1.
run() {
  if [[ "$DRY_RUN" == 1 ]]; then show "$@"; else "$@"; fi
}
# report CMD...: a normalize/compare step, run from the harness directory.
report() {
  if [[ "$DRY_RUN" == 1 ]]; then show "$@"; else ( cd "$SCRIPT_DIR" && "$@" ); fi
}

if [[ "$DRY_RUN" != 1 ]]; then
  mkdir -p "$LOG_DIR"
  if [[ "$MODE" == run-all ]]; then mkdir -p "$OUTROOT"; fi
fi

if [[ "$FORCE" == "1" ]]; then
  echo "[validate-rice] --force: removing existing per-sample outputs under $OUTROOT"
  for s in "${SAMPLES[@]}"; do run rm -rf "${OUTROOT:?}/${s}"; done
fi

done_file() {  # run-all: the output whose presence marks a sample as done
  if [[ "$DO_CHAR" == "1" ]]; then
    echo "$OUTROOT/$1/results/${TARGET}.${TE_NAME}.all_nonref_insert.characTErized.txt"
  else
    echo "$OUTROOT/$1/results/$1.all_nonref_insert.gff"
  fi
}

# Step 1: run RelocaTE3 on the selected samples.
if [[ "$SKIP_RUN" == "1" ]]; then
  echo "[validate-rice] --skip-run: skipping RelocaTE3 execution"
  if [[ "$MODE" == pipeline && ! -d "$OUTROOT" ]]; then
    echo "ERROR: pipeline outdir not found: $OUTROOT" >&2; exit 1
  fi
elif [[ "$MODE" == pipeline ]]; then
  # Input and tag check first (run_pipeline.sh's own dry run: reads the
  # config, checks every input, refuses mutable tags), then the head job.
  echo "[validate-rice] checking the pipeline launch (run_pipeline.sh DRY_RUN=1)"
  ( cd "$REPO" && DRY_RUN=1 CONFIG="$CONFIG" bash "$SCRIPT_DIR/run_pipeline.sh" "$TAG" "$PIPE_OUT" )
  echo "[validate-rice] submitting the pipeline head job (sbatch --wait)"
  export CONFIG   # read by run_pipeline.sh; sbatch exports the environment
  ( cd "$REPO" && run sbatch --wait -o "$LOG_DIR/real_rice_pipeline.%j.log" \
      "$SCRIPT_DIR/run_pipeline.sh" "$TAG" "$PIPE_OUT" )
elif [[ "$EXEC" == "local" ]]; then
  echo "[validate-rice] running RelocaTE3 locally for ${#SAMPLES[@]} sample(s)"
  for s in "${SAMPLES[@]}"; do
    DONE="$(done_file "$s")"
    if [[ -s "$DONE" ]]; then
      echo "  [skip] $s (existing $(basename "$DONE") found; pass --force to redo)"; continue
    fi
    echo "  [run]  $s"
    if [[ "$DRY_RUN" == 1 ]]; then
      show bash "$SCRIPT_DIR/run_relocate3.sh" "$CONFIG" "$s"
    else
      bash "$SCRIPT_DIR/run_relocate3.sh" "$CONFIG" "$s" \
        > "$LOG_DIR/relocate3.${s}.log" 2>&1
    fi
  done
elif [[ ${#EXPLICIT_SAMPLES[@]} -gt 0 ]]; then
  # Explicit sample list on SLURM: submit one sbatch per sample with --wait
  # so we don't depend on the sample's CSV row index for array indexing.
  echo "[validate-rice] submitting ${#SAMPLES[@]} SLURM job(s) on partition ${PART}"
  for s in "${SAMPLES[@]}"; do
    DONE="$(done_file "$s")"
    if [[ -s "$DONE" ]]; then
      echo "  [skip] $s (existing $(basename "$DONE") found; pass --force to redo)"; continue
    fi
    echo "  [submit] $s"
    run sbatch --wait \
      -p "$PART" --cpus-per-task "$CPUS" --mem "$MEM" --time "$TIMELIM" \
      -J "relocate3_${s}" \
      -o "$LOG_DIR/relocate3.${s}.%j.log" \
      "$SCRIPT_DIR/run_relocate3.sh" "$CONFIG" "$s"
  done
else
  N="${#SAMPLES[@]}"
  echo "[validate-rice] submitting SLURM array 1-${N} on partition ${PART}"
  run sbatch --wait \
    -p "$PART" --cpus-per-task "$CPUS" --mem "$MEM" --time "$TIMELIM" \
    --array=1-"$N" \
    -o "$LOG_DIR/relocate3.%A_%a.log" \
    "$SCRIPT_DIR/run_relocate3.sh" "$CONFIG"
fi

# Step 2: normalize. Restrict to the selected samples so the report reflects
# only what was just run.
SAMPLE_FLAGS=()
if [[ ${#EXPLICIT_SAMPLES[@]} -gt 0 ]]; then
  for s in "${SAMPLES[@]}"; do SAMPLE_FLAGS+=("--samples" "$s"); done
fi

# Stage A: raw non-reference calls.
if [[ "$SKIP_NONREF" == "1" ]]; then
  echo "[validate-rice] --skip-nonref: skipping Stage A (raw non-reference compare)"
else
  echo "[validate-rice] Stage A: normalizing legacy RelocaTE2 non-reference calls"
  report python3 "$SCRIPT_DIR/normalize_relocate2_nonref.py" --config "$CONFIG" "${SAMPLE_FLAGS[@]}"
  echo "[validate-rice] Stage A: normalizing RelocaTE3 non-reference calls"
  report python3 "$SCRIPT_DIR/normalize_relocate3_nonref.py" --config "$CONFIG" "${SAMPLE_FLAGS[@]}"
  echo "[validate-rice] Stage A: comparing non-reference call sets"
  report python3 "$SCRIPT_DIR/compare_nonref.py" --config "$CONFIG"
fi

# Stage B: characterized (genotyped) calls.
if [[ "$SKIP_CHAR" == "1" ]]; then
  echo "[validate-rice] --skip-char: skipping Stage B (characterized compare)"
elif [[ "$DO_CHAR" != "1" ]]; then
  echo "[validate-rice] [relocate3].characterize is false; skipping Stage B"
else
  echo "[validate-rice] Stage B: normalizing legacy RelocaTE2 characterized calls"
  report python3 "$SCRIPT_DIR/normalize_relocate2_char.py" --config "$CONFIG" "${SAMPLE_FLAGS[@]}"
  echo "[validate-rice] Stage B: normalizing RelocaTE3 characterized calls"
  report python3 "$SCRIPT_DIR/normalize_relocate3_char.py" --config "$CONFIG" "${SAMPLE_FLAGS[@]}"
  echo "[validate-rice] Stage B: comparing characterized call sets"
  report python3 "$SCRIPT_DIR/compare_char.py" --config "$CONFIG"
fi

# Print summaries at end.
NONREF_SUM="$REPORT_DIR/nonref/summary.txt"
CHAR_SUM="$REPORT_DIR/characterized/summary.txt"
if [[ "$DRY_RUN" != 1 && -f "$NONREF_SUM" ]]; then
  echo; echo "===== Stage A summary ($NONREF_SUM) ====="; cat "$NONREF_SUM"
fi
if [[ "$DRY_RUN" != 1 && -f "$CHAR_SUM" ]]; then
  echo; echo "===== Stage B summary ($CHAR_SUM) ====="; cat "$CHAR_SUM"
fi

echo "[validate-rice] done."
