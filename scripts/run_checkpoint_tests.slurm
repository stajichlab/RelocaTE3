#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH -o logs/checkpoint-tests.%j.log
# epyc MaxTime=30 days checked against cached config on 2026-10-05.
set -euo pipefail
BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
# Python libraries/entry point use the existing project environment. BLAT has
# no listed cluster module; other bioinformatics binaries come from modules.
export PATH="$BASE_DIR/.pixi/envs/default/bin:$BASE_DIR/.pixi/envs/blat/bin:$PATH"
module load bwa/0.7.17 bwa-mem2/2.2.1 samtools/1.22.1 minimap2/2.30 \
  bedtools/2.30.0 bowtie2/2.5.3 seqtk/1.4 bcftools/1.22
module list
export PYTHONPATH="$BASE_DIR/src"
export PYTHONDONTWRITEBYTECODE=1
echo "[$(date -Iseconds)] checkpoint tests job=$SLURM_JOB_ID threads=$THREADS"
.pixi/envs/default/bin/python scripts/run_checkpoint_tests.py
echo "[$(date -Iseconds)] completed"
