#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH -o logs/relocate3-pairing.%A_%a.log
# Validated 2026-09-28: cached cluster config gives epyc MaxTime=30-00:00.
set -euo pipefail
BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
REPLAY_DIR="${1:?Pass the prepared replay directory}"
TASK="${SLURM_ARRAY_TASK_ID:?Submit as an array}"
# Only Python/pysam are used: no aligner, samtools, or bcftools subprocesses.
# The existing project environment supplies the project's Python dependencies.
PYTHON=.pixi/envs/default/bin/python
[[ -x "$PYTHON" && -f "$REPLAY_DIR/manifest.json" ]]
export PYTHONDONTWRITEBYTECODE=1
echo "[$(date)] job=$SLURM_JOB_ID task=$TASK threads=$THREADS replay=$REPLAY_DIR"
"$PYTHON" --version
"$PYTHON" "$REPLAY_DIR/replay_pairing_benchmark.py" run --output "$REPLAY_DIR" --task "$TASK"
echo "[$(date)] paired replay complete"
