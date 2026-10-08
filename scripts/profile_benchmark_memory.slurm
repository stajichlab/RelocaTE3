#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH -o logs/relocate3-memory.%j.log
# 24 h fits epyc MaxTime=30-00:00 in the local SLURM configuration
# inspected 2026-09-21; the baseline sample took 11.2 h with eight threads.
set -euo pipefail
BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
PROFILE_DIR="${1:?Pass a prepared profiling directory relative to the project root}"
[[ -f "$PROFILE_DIR/manifest.json" ]] || { echo "Missing prepared manifest" >&2; exit 1; }
command -v python3.12 >/dev/null
echo "[$(date)] Memory profile starting: job=$SLURM_JOB_ID threads=$THREADS directory=$PROFILE_DIR"
python3.12 "$PROFILE_DIR/snapshot/profile_benchmark_memory.py" run --output "$PROFILE_DIR"
echo "[$(date)] Memory profile and baseline comparison complete"
