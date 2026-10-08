#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH -o logs/family-resolution.%A_%a.log
# Cached cluster config validated epyc MaxTime=30-00:00 on 2026-10-01.
set -euo pipefail
BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
REPLAY_DIR="${1:?Pass the frozen replay directory}"
export PYTHONPATH="$REPLAY_DIR/source/src"
export PYTHONDONTWRITEBYTECODE=1
echo "[$(date -Iseconds)] job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID threads=$THREADS"
.pixi/envs/default/bin/python --version
.pixi/envs/default/bin/python "$REPLAY_DIR/replay_family_resolution.py" run \
  --output "$REPLAY_DIR" --task "${SLURM_ARRAY_TASK_ID:?Submit as array}"
echo "[$(date -Iseconds)] completed"
