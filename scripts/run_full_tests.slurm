#!/usr/bin/bash -l
#SBATCH -p short
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH -o logs/relocate3-tests.%j.log

set -euo pipefail

BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"

THREADS="${SLURM_CPUS_PER_TASK:-1}"
export PATH="$BASE_DIR/.pixi/envs/blat/bin:$PATH"

for tool in python bwa samtools blat; do
  command -v "$tool" >/dev/null 2>&1 || {
    echo "ERROR: required test tool is unavailable: $tool" >&2
    exit 127
  }
done

echo "[$(date)] RelocaTE3 full test suite starting"
echo "  base_dir : $BASE_DIR"
echo "  host     : $(hostname)"
echo "  threads  : $THREADS"
python --version
python -m pytest --version
bwa 2>&1 | head -3 || true
samtools --version | head -2
blat 2>&1 | head -2 || true

python -m pytest -q tests

echo "[$(date)] RelocaTE3 full test suite complete"
