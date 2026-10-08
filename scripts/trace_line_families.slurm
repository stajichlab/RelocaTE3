#!/usr/bin/bash -l
#SBATCH -p short
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH -o logs/line-family-trace.%j.log
set -euo pipefail
BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
OUTDIR="${1:-results/error-audit/2026-10-01-line-family-trace}"
echo "[$(date -Iseconds)] job=$SLURM_JOB_ID threads=$THREADS output=$OUTDIR"
# Existing project Python supplies pysam; no bioinformatics binaries are invoked.
.pixi/envs/default/bin/python --version
.pixi/envs/default/bin/python scripts/trace_line_families.py \
  --benchmark ../../relocate_benchmark/relocate-benchmark \
  --junctions results/error-audit/2026-10-01-line-junctions/junction_evidence.json \
  --outdir "$OUTDIR"
echo "[$(date -Iseconds)] completed"
