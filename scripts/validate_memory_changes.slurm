#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH -o logs/relocate3-memory-validation.%j.log
# epyc MaxTime=30-00:00 verified from local SLURM config on 2026-09-23.
# Full baseline adapter took 11h41m; allow headroom for tests and two trims.
set -euo pipefail
BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
PROFILE_DIR="${1:?Pass the prepared candidate profiling directory}"
BASELINE_DIR="${2:?Pass the completed baseline profiling directory}"
[[ -f "$PROFILE_DIR/manifest.json" && -f "$BASELINE_DIR/.profile_complete" ]]
# Exclusive claim: never restart partially completed validation silently.
mkdir "$PROFILE_DIR/validation-started"
export RT3_REPO="$BASE_DIR"
export ADAPTER_DIR="$PROFILE_DIR/snapshot/benchmark/callers/relocate3"
source "$ADAPTER_DIR/env.sh"
export PYTHONPATH
PYTHONPATH="$(realpath "$PROFILE_DIR/snapshot/source/src")"
export PYTHONDONTWRITEBYTECODE=1

echo "[$(date)] job=$SLURM_JOB_ID threads=$THREADS candidate=$PROFILE_DIR baseline=$BASELINE_DIR"
python -c 'import os, pathlib, RelocaTE3; p=pathlib.Path(RelocaTE3.__file__).resolve(); print(p); assert p.is_relative_to(pathlib.Path(os.environ["PYTHONPATH"]).resolve())'
python -m pytest -q tests --junitxml="$PROFILE_DIR/full-tests.xml"
python3.12 "$PROFILE_DIR/snapshot/validate_trim_memory.py" \
  --baseline "$BASELINE_DIR" --candidate "$PROFILE_DIR"
# Full single-sample replay also validates the BLAT streaming change and calls.
python3.12 "$PROFILE_DIR/snapshot/profile_benchmark_memory.py" run --output "$PROFILE_DIR"
echo "[$(date)] Full suite, paired trim profiles, and complete single-sample replay passed"
