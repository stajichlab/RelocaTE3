#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=08:00:00
#SBATCH -o logs/relocate3.%A_%a.log
#
# Run RelocaTE3 on each rice sample listed in the config's sample CSV.
#
# Usage:
#   # one sample (interactive / single job)
#   ./run_relocate3.sh config.toml <sample_name>
#
#   # whole sheet as a SLURM array (one task per sample)
#   sbatch --array=1-$(($(wc -l < /path/to/samples.csv) - 1)) \
#          ./run_relocate3.sh config.toml
#
# The script expects to be launched from the directory containing it (so the
# default ``logs/`` SBATCH path resolves). Outputs go to paths.relocate3_outdir
# from the config; logs go to paths.log_dir.

set -euo pipefail

CONFIG="${1:-config.toml}"
SAMPLE_ARG="${2:-}"

# Resolve the config to an absolute path first. Under sbatch, BASH_SOURCE
# points to /var/spool/slurmd/<job>/slurm_script (a copy of this script), so
# we anchor SCRIPT_DIR on the config's directory instead — that's the dir
# that holds _config.py and config.toml.
if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: config file not found: $CONFIG" >&2
  echo "       Copy config.example.toml -> config.toml and edit." >&2
  exit 1
fi
CONFIG="$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")"
SCRIPT_DIR="$(dirname "$CONFIG")"
cd "$SCRIPT_DIR"

# Make sure relocaTE3 is on PATH. SLURM compute nodes don't inherit the pixi
# activation from the submitting shell, so we re-exec ourselves inside the
# project's pixi env if needed. Skip if PIXI_INSIDE_RUN is already set, which
# pixi uses to prevent recursive activation.
if [[ -z "${PIXI_INSIDE_RUN:-}" ]] && ! command -v relocaTE3 >/dev/null 2>&1; then
  PIXI_MANIFEST=""
  search_dir="$SCRIPT_DIR"
  while [[ "$search_dir" != "/" && "$search_dir" != "" ]]; do
    if [[ -f "$search_dir/pixi.toml" ]]; then
      PIXI_MANIFEST="$search_dir/pixi.toml"; break
    fi
    search_dir="$(dirname "$search_dir")"
  done
  if [[ -n "$PIXI_MANIFEST" ]] && command -v pixi >/dev/null 2>&1; then
    echo "[run_relocate3] activating pixi env from $PIXI_MANIFEST"
    export PIXI_INSIDE_RUN=1
    exec pixi run --manifest-path "$PIXI_MANIFEST" -- bash "$0" "$@"
  else
    echo "ERROR: relocaTE3 not on PATH and could not locate pixi.toml or pixi binary." >&2
    echo "       Activate the project env (pixi shell) before launching, or install pixi." >&2
    exit 127
  fi
fi

# Pull values out of the TOML using the shared loader.
read_cfg() {
  python3 -c "
from _config import load_config
cfg = load_config('$CONFIG')
print(cfg['paths']['genome'])
print(cfg['paths']['te_library'])
print(cfg['paths']['repeatmasker'])
print(cfg['paths']['fq_dir'])
print(cfg['paths']['sample_csv'])
print(cfg['paths']['relocate3_outdir'])
print(cfg['paths']['log_dir'])
print(cfg['reads']['r1_suffix'])
print(cfg['reads']['r2_suffix'])
print(cfg['relocate3']['aligner'])
print(cfg['relocate3']['mismatch'])
print(cfg['relocate3']['min_match'])
print(cfg['relocate3']['min_trimmed'])
print(cfg['relocate3']['min_mapq'])
print(cfg['relocate3']['te_name'])
print(cfg['relocate3']['target'])
print(cfg['relocate3']['tsd'])
print(cfg['relocate3']['threads'])
print(int(bool(cfg['relocate3']['characterize'])))
"
}

mapfile -t CFG < <(read_cfg)
GENOME="${CFG[0]}"
TE_LIB="${CFG[1]}"
REPMASK="${CFG[2]}"
FQ_DIR="${CFG[3]}"
SAMPLE_CSV="${CFG[4]}"
OUTROOT="${CFG[5]}"
LOG_DIR="${CFG[6]}"
R1_SFX="${CFG[7]}"
R2_SFX="${CFG[8]}"
ALIGNER="${CFG[9]}"
MISMATCH="${CFG[10]}"
MIN_MATCH="${CFG[11]}"
MIN_TRIMMED="${CFG[12]}"
MIN_MAPQ="${CFG[13]}"
TE_NAME="${CFG[14]}"
TARGET="${CFG[15]}"
TSD="${CFG[16]}"
CFG_THREADS="${CFG[17]}"
DO_CHARACTERIZE="${CFG[18]}"

mkdir -p "$OUTROOT" "$LOG_DIR"

# Resolve which sample to process: explicit arg > SLURM_ARRAY_TASK_ID > error.
if [[ -n "$SAMPLE_ARG" ]]; then
  SAMPLE="$SAMPLE_ARG"
elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  # First data row = task 1 (skip header).
  SAMPLE=$(tail -n +2 "$SAMPLE_CSV" | sed -n "${SLURM_ARRAY_TASK_ID}p" | cut -d',' -f1)
else
  echo "ERROR: provide a sample name or run under SLURM with --array=<N>" >&2
  exit 1
fi

if [[ -z "$SAMPLE" ]]; then
  echo "ERROR: empty sample for task ${SLURM_ARRAY_TASK_ID:-?}" >&2
  exit 1
fi

R1="${FQ_DIR}/${SAMPLE}${R1_SFX}"
R2="${FQ_DIR}/${SAMPLE}${R2_SFX}"

for f in "$GENOME" "$TE_LIB" "$REPMASK" "$R1" "$R2"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: required input missing: $f" >&2
    exit 1
  fi
done

THREADS="${SLURM_CPUS_PER_TASK:-$CFG_THREADS}"
SAMPLE_OUTDIR="${OUTROOT}/${SAMPLE}"
mkdir -p "$SAMPLE_OUTDIR"

# Paths emitted by each pipeline step (must match the layout in
# src/RelocaTE3/librelocate.py + align.py + insertions.py).
READ_REPEAT="${SAMPLE_OUTDIR}/te_containing/${SAMPLE}.read_repeat_name.txt"
FLANK5="${SAMPLE_OUTDIR}/flanking/${SAMPLE}.left.flankingReads.fq"
FLANK3="${SAMPLE_OUTDIR}/flanking/${SAMPLE}.right.flankingReads.fq"
GENOME_BAM="${SAMPLE_OUTDIR}/${SAMPLE}.repeat.minimap.sorted.bam"
NONREF_TXT="${SAMPLE_OUTDIR}/results/${TARGET}.${TE_NAME}.all_nonref_insert.txt"

echo "[$(date)] RelocaTE3 validation run"
echo "  sample      : $SAMPLE"
echo "  genome      : $GENOME"
echo "  TE library  : $TE_LIB"
echo "  R1          : $R1"
echo "  R2          : $R2"
echo "  outdir      : $SAMPLE_OUTDIR"
echo "  threads     : $THREADS"
echo "  aligner     : $ALIGNER"
echo "  mismatch    : $MISMATCH"
echo "  match/trim  : ${MIN_MATCH}/${MIN_TRIMMED}"
echo "  te_name/TSD : ${TE_NAME}/${TSD}"
echo "  target      : $TARGET"
echo "  characterize: $DO_CHARACTERIZE"

start=$(date +%s)

# ---------------------------------------------------------------------------
# Step 1: index the reference genome (skipped if .fai + .mmi already exist).
# ---------------------------------------------------------------------------
if [[ ! -s "${GENOME}.fai" || ! -s "${GENOME}.mmi" ]]; then
  echo "[$(date)] index-genome"
  relocaTE3 index-genome -g "$GENOME"
else
  echo "[$(date)] genome already indexed (.fai + .mmi present), skipping index-genome"
fi

# ---------------------------------------------------------------------------
# Steps 2+3: align reads to the TE library and trim TE sequence from them.
# ---------------------------------------------------------------------------
if [[ ! -s "$READ_REPEAT" ]]; then
  echo "[$(date)] run (map + trim)"
  relocaTE3 run \
    --left "$R1" \
    --right "$R2" \
    --te-library "$TE_LIB" \
    --name "$SAMPLE" \
    --outdir "$SAMPLE_OUTDIR" \
    --threads "$THREADS" \
    --aligner "$ALIGNER" \
    --min-match "$MIN_MATCH" \
    --min-trimmed "$MIN_TRIMMED" \
    --mismatch "$MISMATCH"
else
  echo "[$(date)] read_repeat table exists, skipping run (map + trim)"
fi

# ---------------------------------------------------------------------------
# Step 4: align trimmed flanking reads back to the reference genome.
# ---------------------------------------------------------------------------
FLANK_INPUTS=()
[[ -s "$FLANK5" ]] && FLANK_INPUTS+=("$FLANK5")
[[ -s "$FLANK3" ]] && FLANK_INPUTS+=("$FLANK3")
if [[ ${#FLANK_INPUTS[@]} -eq 0 ]]; then
  echo "ERROR: no flanking FASTQs produced under ${SAMPLE_OUTDIR}/flanking/" >&2
  exit 1
fi

if [[ ! -s "$GENOME_BAM" ]]; then
  echo "[$(date)] align-genome (flanking -> genome BAM)"
  relocaTE3 align-genome \
    -g "$GENOME" \
    -f "${FLANK_INPUTS[@]}" \
    -n "$SAMPLE" \
    -o "$SAMPLE_OUTDIR" \
    --threads "$THREADS"
else
  echo "[$(date)] genome-aligned BAM exists, skipping align-genome"
fi

# ---------------------------------------------------------------------------
# Step 5: cluster junction/supporting reads and call non-reference insertions.
# ---------------------------------------------------------------------------
echo "[$(date)] find-insertions"
relocaTE3 find-insertions \
  -b "$GENOME_BAM" \
  --read-repeat "$READ_REPEAT" \
  --tsd "$TSD" \
  --target "$TARGET" \
  --name "$SAMPLE" \
  --outdir "$SAMPLE_OUTDIR" \
  --te-name "$TE_NAME" \
  --reference-ins "$REPMASK" \
  --mismatch "$MISMATCH" \
  --min-mapq "$MIN_MAPQ"

if [[ ! -s "$NONREF_TXT" ]]; then
  echo "WARN: expected non-reference insertion table not produced: $NONREF_TXT" >&2
fi

# ---------------------------------------------------------------------------
# Optional Step 7: characterize zygosity. Needs reads aligned to the genome —
# the legacy validation already has CRAM/BAMs under validation_data/real_rice/aln_input/.
# ---------------------------------------------------------------------------
if [[ "$DO_CHARACTERIZE" == "1" ]]; then
  echo "[$(date)] characterize (step 7) — skipped: needs reads-to-genome BAM not produced here"
  echo "       Pass --characterize false in the config or extend run_relocate3.sh to wire this up."
fi

end=$(date +%s)
echo "[$(date)] done in $((end-start))s"
echo "  insertions table: $NONREF_TXT"
