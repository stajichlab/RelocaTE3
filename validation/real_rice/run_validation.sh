#!/usr/bin/bash -l
# End-to-end real-rice validation:
#   1. Make sure config.toml exists (copy from config.example.toml if not).
#   2. Run RelocaTE3 on every sample (SLURM array with --wait by default, or
#      sequentially in the current shell when --local is given).
#   3. Normalize both call sets (RelocaTE2 + RelocaTE3) into report/*.tsv.
#   4. Run the comparison and emit summary + venn.
#
# Designed to be launched from the RelocaTE3/ working directory:
#   pixi run validate-rice                 # uses defaults (SLURM array)
#   bash validation/real_rice/run_validation.sh
#   bash validation/real_rice/run_validation.sh --local
#   bash validation/real_rice/run_validation.sh --skip-run        # only re-report
#   bash validation/real_rice/run_validation.sh --config my.toml
#
# Notes:
#   * The SLURM path uses `sbatch --wait` so the script blocks until the array
#     finishes, then runs normalize + compare. `--local` skips SLURM entirely.
#   * Re-runs are safe: existing per-sample result directories under
#     paths.relocate3_outdir are reused; pass --force to wipe them first.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CONFIG_REL="validation/real_rice/config.toml"
EXAMPLE_REL="validation/real_rice/config.example.toml"

MODE="slurm"           # slurm | local
SKIP_RUN=0
FORCE=0
CONFIG_ARG=""
EXPLICIT_SAMPLES=()    # positional sample names; empty = use the config's sample CSV

while [[ $# -gt 0 ]]; do
  case "$1" in
    --local)     MODE="local"; shift ;;
    --slurm)     MODE="slurm"; shift ;;
    --skip-run)  SKIP_RUN=1;   shift ;;
    --force)     FORCE=1;      shift ;;
    --config)    CONFIG_ARG="$2"; shift 2 ;;
    -h|--help)
      sed -n '2,30p' "$0"; exit 0 ;;
    --)          shift; while [[ $# -gt 0 ]]; do EXPLICIT_SAMPLES+=("$1"); shift; done ;;
    -*)
      echo "ERROR: unknown flag: $1" >&2; exit 2 ;;
    *)
      EXPLICIT_SAMPLES+=("$1"); shift ;;
  esac
done

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
    echo "[validate-rice] config.toml not found; copying from config.example.toml"
    cp "$SCRIPT_DIR/config.example.toml" "$CONFIG"
  fi
fi

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
"
)
SAMPLE_CSV="${CFG[0]}"
OUTROOT="${CFG[1]}"
LOG_DIR="${CFG[2]}"
PART="${CFG[3]}"
CPUS="${CFG[4]}"
MEM="${CFG[5]}"
TIMELIM="${CFG[6]}"
# Which samples to operate on: explicit positional args win over the CSV.
if [[ ${#EXPLICIT_SAMPLES[@]} -gt 0 ]]; then
  SAMPLES=("${EXPLICIT_SAMPLES[@]}")
  echo "[validate-rice] using ${#SAMPLES[@]} explicit sample(s): ${SAMPLES[*]}"
else
  IFS=',' read -r -a SAMPLES <<< "${CFG[7]}"
fi

mkdir -p "$OUTROOT" "$LOG_DIR"

if [[ "$FORCE" == "1" ]]; then
  echo "[validate-rice] --force: removing existing per-sample outputs under $OUTROOT"
  for s in "${SAMPLES[@]}"; do rm -rf "${OUTROOT:?}/${s}"; done
fi

# Step 1: run RelocaTE3 on the selected samples.
if [[ "$SKIP_RUN" == "1" ]]; then
  echo "[validate-rice] --skip-run: skipping RelocaTE3 execution"
elif [[ "$MODE" == "local" ]]; then
  echo "[validate-rice] running RelocaTE3 locally for ${#SAMPLES[@]} sample(s)"
  for s in "${SAMPLES[@]}"; do
    if [[ -s "$OUTROOT/$s/results/${s}.all_nonref_insert.gff" ]]; then
      echo "  [skip] $s (existing GFF found; pass --force to redo)"; continue
    fi
    echo "  [run]  $s"
    bash "$SCRIPT_DIR/run_relocate3.sh" "$CONFIG" "$s" \
      > "$LOG_DIR/relocate3.${s}.log" 2>&1
  done
elif [[ ${#EXPLICIT_SAMPLES[@]} -gt 0 ]]; then
  # Explicit sample list on SLURM: submit one sbatch per sample with --wait
  # so we don't depend on the sample's CSV row index for array indexing.
  echo "[validate-rice] submitting ${#SAMPLES[@]} SLURM job(s) on partition ${PART}"
  for s in "${SAMPLES[@]}"; do
    if [[ -s "$OUTROOT/$s/results/${s}.all_nonref_insert.gff" ]]; then
      echo "  [skip] $s (existing GFF found; pass --force to redo)"; continue
    fi
    echo "  [submit] $s"
    sbatch --wait \
      -p "$PART" --cpus-per-task "$CPUS" --mem "$MEM" --time "$TIMELIM" \
      -J "relocate3_${s}" \
      -o "$LOG_DIR/relocate3.${s}.%j.log" \
      "$SCRIPT_DIR/run_relocate3.sh" "$CONFIG" "$s"
  done
else
  N="${#SAMPLES[@]}"
  echo "[validate-rice] submitting SLURM array 1-${N} on partition ${PART}"
  sbatch --wait \
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

echo "[validate-rice] normalizing legacy RelocaTE2 calls"
( cd "$SCRIPT_DIR" && python3 normalize_relocate2.py --config "$CONFIG" "${SAMPLE_FLAGS[@]}" )
echo "[validate-rice] normalizing RelocaTE3 calls"
( cd "$SCRIPT_DIR" && python3 normalize_relocate3.py --config "$CONFIG" "${SAMPLE_FLAGS[@]}" )

# Step 3: compare + report.
echo "[validate-rice] comparing call sets"
( cd "$SCRIPT_DIR" && python3 compare_calls.py --config "$CONFIG" )

echo "[validate-rice] done."
