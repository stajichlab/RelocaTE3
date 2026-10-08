#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=8
#SBATCH --time=96:00:00
# 96h is headroom, not the expected runtime. Re-measured 2026-08-14 across the
# full riceTElib and divergence panels (arrays 27443131/27443164, 126 tasks):
# slowest task 2h01m, and mPing tops out at 55m. Most tasks finish in 10-20 min.
#
# The older note here claimed riceTElib 5x = 8h44m, extrapolating to 15x ~= 26h
# and 30x ~= 52h (array 27262453, 2026-08-06, before the BLAT backend fix).
# Those figures are obsolete and were actively harmful: they are ~35x the
# measured time and led to a panel being scoped down to 6 tasks on the belief
# that a full sweep could not finish overnight. Trust the measurement above.
#SBATCH -o logs/benchmark.%A_%a.log
#
# One array task = one (dataset, caller, sample) tuple. Submit only through
# pipeline/submit_benchmark.sh so DATASET_SELECTION and CONFIG are consistent.
set -euo pipefail

BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"
CONFIG="${CONFIG:-config/benchmark.toml}"
DATASET_SELECTION="${DATASET_SELECTION:?submit via pipeline/submit_benchmark.sh}"
PY=python3.12

echo "[$(date)] benchmark array task starting"
echo "  host      : $(hostname)"
echo "  config    : $CONFIG"
echo "  selection : $DATASET_SELECTION"

eval "$("$PY" pipeline/config_env.py --config "$CONFIG" globals)"

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID required}"
LINE="$(
  "$PY" pipeline/config_env.py --config "$CONFIG" \
    --dataset "$DATASET_SELECTION" tasks |
    sed -n "$((TASK_ID + 1))p"
)"
if [[ -z "$LINE" ]]; then
  echo "ERROR: no task for array index $TASK_ID" >&2
  exit 1
fi
IFS=$'\t' read -r DATASET CALLER SAMPLE COVERAGE REPLICATE R1 R2 <<< "$LINE"
eval "$("$PY" pipeline/config_env.py --config "$CONFIG" dataset-env "$DATASET")"
THREADS="${SLURM_CPUS_PER_TASK:-$THREADS}"

echo "  task_id   : $TASK_ID"
echo "  dataset   : $DATASET ($DATASET_LABEL)"
echo "  caller    : $CALLER"
echo "  sample    : $SAMPLE (cov=$COVERAGE rep=$REPLICATE)"
echo "  R1        : $R1"
echo "  R2        : $R2"
echo "  threads   : $THREADS"

SAMPLE_TRUTH="$DATASET_TRUTH_ROOT/per_sample/${SAMPLE}.tsv"
if [[ ! -s "$SAMPLE_TRUTH" ]]; then
  echo "ERROR: exported sample truth missing or empty: $SAMPLE_TRUTH" >&2
  exit 1
fi
echo "  truth     : $SAMPLE_TRUTH"

OUTDIR="$DATASET_WORK_ROOT/$CALLER/$SAMPLE"
REPORT_DIR="$DATASET_REPORT_ROOT/per_sample/$CALLER/$SAMPLE"
RES_DIR="$DATASET_REPORT_ROOT/resources/$CALLER"
mkdir -p "$RES_DIR" logs

export DATASET SAMPLE R1 R2 OUTDIR THREADS
export REFERENCE TE_LIBRARY REPEATMASKER TE_NAME

eval "$("$PY" pipeline/config_env.py --config "$CONFIG" caller-env "$CALLER")"
# Allowlist of caller-env variables to export. A caller key added to
# ADAPTER_ENV_MAP in config_env.py must ALSO be listed here, or it is silently
# dropped and the adapter falls back to its default -- which is how the
# `-strict` variant first ran without its own flag.
for variable in RT3_REPO TSD_PATTERN RT3_TE_ALIGNER RT3_GENOME_ALIGNER \
                RT3_REQUIRE_BOTH_JUNCTIONS \
                RT2_ALIGNER RT2_SIZE RT2_MISMATCH; do
  if [[ -n "${!variable+x}" ]]; then
    export "${variable?}"
  fi
done
ADAPTER="$("$PY" pipeline/config_env.py --config "$CONFIG" adapter "$CALLER")"
export ADAPTER_DIR="$ADAPTER"

TIME_FILE="$RES_DIR/${SAMPLE}.time-v.txt"
echo "[$(date)] running caller adapter: $ADAPTER/run.sh"
/usr/bin/time -v -o "$TIME_FILE" bash "$ADAPTER/run.sh"

echo "[$(date)] normalizing calls"
"$PY" "$ADAPTER/normalize.py" \
  --outdir "$OUTDIR" \
  --sample "$SAMPLE" \
  --te-name "$TE_NAME" \
  --target ALL

echo "[$(date)] scoring calls"
"$PY" scoring/score_calls.py \
  --truth "$SAMPLE_TRUTH" \
  --calls "$OUTDIR/calls.normalized.tsv" \
  --sample "$SAMPLE" \
  --caller "$CALLER" \
  --window "$MATCH_WINDOW" \
  --coordinate-policy "${SCORING_COORDINATE_POLICY:-legacy}" \
  --outdir "$REPORT_DIR"

echo "[$(date)] parsing resource usage"
"$PY" scoring/parse_time_v.py \
  --time-v "$TIME_FILE" \
  --sample "$SAMPLE" \
  --caller "$CALLER" \
  --coverage "$COVERAGE" \
  --replicate "$REPLICATE" \
  --out "$RES_DIR/${SAMPLE}.tsv"

echo "[$(date)] task $TASK_ID complete: $DATASET / $CALLER / $SAMPLE"
