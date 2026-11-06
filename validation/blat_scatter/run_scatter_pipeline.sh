#!/usr/bin/bash -l
#SBATCH -p intel
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH -J rt3-m4-scatter
#SBATCH -o logs/nextflow_m4_scatter.%j.log
#
# M4 gate: the riceTElib 30x sample through the pipeline with BLAT scatter
# (TE_SPLIT -> TE_SEARCH_CHUNK x n -> TE_MERGE_TRIM), on SLURM (-profile hpcc,pixi).
# Its TE-search outputs are compared with the unsplit baseline
# (run_unsplit_baseline.sh). Genotyping and MultiQC are off: this run is about
# the TE search; the rest of the pipeline is unchanged by M4.
#
# The pixi env loads the worktree's src (editable): do not edit src/ while it runs.
#
# Usage, from the repository root:
#   sbatch validation/blat_scatter/run_scatter_pipeline.sh [EXTRA nextflow params...]
#   e.g. ... --blat_chunks_per_task 4
set -euo pipefail

P="/bigdata/wesslerlab/shared/Rice/Nathan/rice/make_simulated_genome/make_simulation_new/results/riceTElib_benchmark"
SAMPLE="${SAMPLE:-cov30x_rep1}"
OUT="${OUTDIR:-$(pwd)/results/nextflow-m4/scatter/${SAMPLE}}"

echo "date     : $(date -Is)"
echo "commit   : $(git describe --tags --always --dirty)"
echo "outdir   : ${OUT}"
echo "extra    : $*"
[[ -e "${OUT}/nf" ]] && { echo "ERROR: ${OUT}/nf exists; refusing to overwrite" >&2; exit 1; }
mkdir -p "${OUT}"
printf 'sample,fastq_1,fastq_2\n%s,%s,%s\n' "${SAMPLE}" \
  "${P}/reads/${SAMPLE}/R1.fastq.gz" "${P}/reads/${SAMPLE}/R2.fastq.gz" > "${OUT}/samplesheet.csv"
cd "${OUT}"   # own launch dir, so -resume picks this run
pixi run --manifest-path "${OLDPWD}/pixi.toml" --frozen -e nf nextflow run "${OLDPWD}" \
  -profile hpcc,pixi --input "${OUT}/samplesheet.csv" \
  --genome "${P}/genomes/reference.fa" --te_library "${P}/selected_te_library.fa" \
  --genotype false --skip_multiqc \
  --outdir "${OUT}/nf" -work-dir "${OUT}/work" -resume "$@"
echo "done     : $(date -Is)"
