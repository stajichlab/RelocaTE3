#!/usr/bin/bash -l
#SBATCH -p short
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=02:00:00
#SBATCH -J rt3-release-check
#SBATCH -o logs/release-check.%j.log
#
# Release check for one tag or commit, on a compute node.
#
#   A. Fresh clone at REF, `pixi install -e blat --locked`, full pytest suite.
#   B. Build the sdist and wheel, install each into a clean venv outside the
#      clone, check the version and entry point, then run the Chr3 acceptance
#      tests and a default-aligner measurement against the installed wheel.
#
# Usage, from a RelocaTE3 checkout (outputs go under its results/):
#   sbatch validation/release/fresh_clone_check.slurm v0.2.0rc1 [expected-version]
#
# Outputs: results/release-validation/<REF>/ (clone/, dist/, *.xml, summary.txt,
# .complete). Refuses to reuse a partial attempt; move it aside to rerun.

set -euo pipefail

REF="${1:?usage: sbatch fresh_clone_check.slurm <tag-or-commit> [expected-version]}"
EXPECTED_VERSION="${2:-${REF#v}}"
BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

SRC_REPO="$(cd "$(git rev-parse --git-common-dir)" && pwd)"
OUTDIR="$BASE_DIR/results/release-validation/${REF}"
CLONE="$OUTDIR/clone"
DIST="$OUTDIR/dist"
SUMMARY="$OUTDIR/summary.txt"

# Nothing from the submitting shell may leak a source tree into the checks.
unset PYTHONPATH PYTHONHOME PIXI_INSIDE_RUN PIXI_PROJECT_MANIFEST

if [[ -e "$OUTDIR/.complete" ]]; then
  echo "[$(date)] $OUTDIR/.complete exists; nothing to do"
  exit 0
fi
if [[ -e "$OUTDIR" ]]; then
  echo "ERROR: $OUTDIR exists from an earlier attempt; move it aside to rerun" >&2
  exit 1
fi
command -v pixi >/dev/null 2>&1 || { echo "ERROR: pixi not on PATH" >&2; exit 127; }
mkdir -p "$OUTDIR" "$DIST"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$SUMMARY"; }

log "RelocaTE3 release check"
log "  ref             : $REF"
log "  expected version: $EXPECTED_VERSION"
log "  source repo     : $SRC_REPO"
log "  outdir          : $OUTDIR"
log "  host            : $(hostname)"
log "  threads         : $THREADS"
log "  pixi            : $(pixi --version)"

# --------------------------------------------------------------------------
# A. Fresh clone + locked environment + full test suite
# --------------------------------------------------------------------------
git clone --quiet "$SRC_REPO" "$CLONE"
git -C "$CLONE" checkout --quiet --detach "$REF"
log "  commit          : $(git -C "$CLONE" rev-parse HEAD)"
log "  describe        : $(git -C "$CLONE" describe --tags --always --dirty)"

cd "$CLONE"
log "A1 pixi install -e blat --locked"
pixi install -e blat --locked
ENV_BIN="$CLONE/.pixi/envs/blat/bin"

for tool in python blat bwa samtools seqtk bedtools minimap2 relocaTE3; do
  [[ -x "$ENV_BIN/$tool" ]] || { log "FAIL: $tool missing from the blat env"; exit 1; }
done
log "  relocaTE3 (editable): $("$ENV_BIN/relocaTE3" --version 2>&1)"
log "  python   : $("$ENV_BIN/python" --version 2>&1)"
log "  samtools : $("$ENV_BIN/samtools" --version | head -1)"
log "  bwa      : $("$ENV_BIN/bwa" 2>&1 | grep -m1 Version || true)"
log "  blat     : $("$ENV_BIN/blat" 2>&1 | head -1 || true)"
log "  seqtk    : $("$ENV_BIN/seqtk" 2>&1 | grep -m1 Version || true)"

log "A2 full pytest suite (fresh clone, blat env)"
set +e
pixi run -e blat --locked python -m pytest -ra -q -p no:cacheprovider \
  --junitxml="$OUTDIR/pytest-fresh-clone.xml" tests
A_RC=$?
set -e
log "  pytest exit: $A_RC"
cd "$OUTDIR"

# --------------------------------------------------------------------------
# B. Built artifacts installed outside the clone
# --------------------------------------------------------------------------
log "B1 build sdist + wheel"
"$ENV_BIN/python" -m venv "$OUTDIR/venv-build"
"$OUTDIR/venv-build/bin/python" -m pip install --quiet build
"$OUTDIR/venv-build/bin/python" -m build --outdir "$DIST" "$CLONE"
( cd "$DIST" && sha256sum ./* > SHA256SUMS )
sed 's/^/  /' "$DIST/SHA256SUMS" | tee -a "$SUMMARY"

WHEEL="$(ls "$DIST"/*.whl)"
SDIST="$(ls "$DIST"/*.tar.gz)"

check_install() {  # <label> <artifact>
  local label="$1" artifact="$2" venv="$OUTDIR/venv-$1"
  "$ENV_BIN/python" -m venv "$venv"
  "$venv/bin/python" -m pip install --quiet "$artifact" pytest
  local got where
  got="$("$venv/bin/relocaTE3" --version 2>&1)"
  where="$("$venv/bin/python" -c 'import RelocaTE3, os; print(os.path.dirname(RelocaTE3.__file__))')"
  log "  $label: version=$got module=$where"
  [[ "$got" == "$EXPECTED_VERSION" ]] || { log "FAIL: $label version $got != $EXPECTED_VERSION"; return 1; }
  [[ "$where" == "$venv"/* ]] || { log "FAIL: $label imported RelocaTE3 from outside its venv"; return 1; }
  "$venv/bin/relocaTE3" run-all --help >/dev/null
}

log "B2 install wheel and sdist into clean venvs"
check_install wheel "$WHEEL"
check_install sdist "$SDIST"

# Installed wheel first on PATH; external aligners from the locked env after it.
export PATH="$OUTDIR/venv-wheel/bin:$ENV_BIN:$PATH"
[[ "$(command -v relocaTE3)" == "$OUTDIR/venv-wheel/bin/relocaTE3" ]] \
  || { log "FAIL: relocaTE3 on PATH is not the installed wheel"; exit 1; }

log "B3 acceptance tests against the installed wheel"
cp -r "$CLONE/tests" "$OUTDIR/installed-tests"
set +e
"$OUTDIR/venv-wheel/bin/python" -m pytest -ra -q -p no:cacheprovider \
  --junitxml="$OUTDIR/pytest-installed-wheel.xml" \
  "$OUTDIR/installed-tests/acceptance_test.py"
B_RC=$?
set -e
log "  pytest exit: $B_RC"

log "B4 Chr3 2 Mb measurement, shipped defaults (blat + bwa aln)"
M="$OUTDIR/chr3-defaults"
mkdir -p "$M"
cp "$CLONE/tests/data/sim_genome/MSU7.Chr3_2M.fa" "$M/genome.fa"
D="$CLONE/tests/data"
relocaTE3 run-all \
  -1 "$D/sim_reads/MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz" \
  -2 "$D/sim_reads/MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz" \
  -T "$D/TE_lib/RiceTE.fa" -g "$M/genome.fa" -n HEG4 -o "$M/out" \
  --threads "$THREADS" --mismatch 2 --tsd UNK --te-name repeat \
  --require-both-junctions \
  --repeatmasker "$D/sim_genome/MSU7.Chr3_2M.fa.RepeatMasker.out" --genotype \
  > "$M/run-all.log" 2>&1
"$OUTDIR/venv-wheel/bin/python" - "$D/sim_genome/MSU7.Chr3_2M.ALL.gff" "$M/out/results" \
  > "$M/metrics.tsv" <<'EOF'
import sys
from pathlib import Path

def load(path):
    rows = []
    for line in open(path):
        f = line.split("\t")
        if line.startswith("#") or len(f) < 5:
            continue
        rows.append((f[0], int(f[3]), int(f[4])))
    return rows

def hit(a, b, w=10):
    return a[0] == b[0] and a[1] - w <= b[2] and a[2] + w >= b[1]

truth = load(sys.argv[1])
print("tier\tcalls\trecovered\ttruth\tprecision\tf1")
for name in ("ALL.repeat.all_nonref_insert.gff",
             "ALL.repeat.all_nonref_insert.characTErized.gff"):
    path = Path(sys.argv[2]) / name
    if not path.exists():
        print(f"{name}\tMISSING")
        continue
    calls = load(path)
    rec = sum(any(hit(t, c) for c in calls) for t in truth)
    tp = sum(any(hit(c, t) for t in truth) for c in calls)
    p = tp / len(calls) if calls else 0.0
    r = rec / len(truth)
    f1 = 2 * p * r / (p + r) if p + r else 0.0
    print(f"{name}\t{len(calls)}\t{rec}\t{len(truth)}\t{p:.4f}\t{f1:.4f}")
EOF
sed 's/^/  /' "$M/metrics.tsv" | tee -a "$SUMMARY"

if [[ "$A_RC" -ne 0 || "$B_RC" -ne 0 ]]; then
  log "RESULT: FAIL (fresh-clone pytest=$A_RC, installed-wheel pytest=$B_RC)"
  exit 1
fi
log "RESULT: PASS"
touch "$OUTDIR/.complete"
