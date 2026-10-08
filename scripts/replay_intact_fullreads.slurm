#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH -o logs/intact-fullreads.%A_%a.log
# epyc MaxTime=30 days verified from cached cluster config on 2026-10-02.
set -euo pipefail
BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
REPLAY_DIR="${1:?Pass the frozen replay directory}"
test -f "$REPLAY_DIR/replay_intact_fullreads.py"
export PYTHONPATH="$REPLAY_DIR/source/src"
export PYTHONDONTWRITEBYTECODE=1
echo "[$(date -Iseconds)] job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID threads=$THREADS"
.pixi/envs/default/bin/python "$REPLAY_DIR/replay_intact_fullreads.py" run \
  --output "$REPLAY_DIR" --task "${SLURM_ARRAY_TASK_ID:?Submit as array}"
echo "[$(date -Iseconds)] completed"
