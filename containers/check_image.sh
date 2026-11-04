#!/usr/bin/env bash
# Runs INSIDE the RelocaTE3 image (no arguments):
#
#   docker run --rm <image> bash /opt/relocate3/check_image.sh
#   apptainer exec docker://<image> bash /opt/relocate3/check_image.sh
#
# Prints, then asserts:
#   - relocaTE3 --version derives from the git tag recorded at build time
#     (never 0.0.0 or versioningit's default-tag fallback), equals the version
#     build_env.sh installed (relocate3_version), and has no versioningit
#     dirty-tree suffix (+d<YYYYMMDD> / .d<YYYYMMDD>);
#   - blat, bwa, samtools, bcftools, minimap2, seqtk, bedtools, bowtie2 run
#     and report the versions locked in pixi.lock's blat env
#     (expected_versions.tsv, written from `pixi list -e blat --locked` during
#     the build). bwa-mem2 must run and print a version, but it is not
#     compared: bioconda's bwa-mem2 2.3 binary reports 2.2.1;
#   - ps exists (Nextflow task metrics);
#   - python imports RelocaTE3 and pysam, RelocaTE3 non-editable, from the env.
# Exits non-zero if any check fails.
set -uo pipefail

EXPECTED="${EXPECTED_VERSIONS:-/opt/relocate3/meta/expected_versions.tsv}"
# versioningit default-tag in pyproject.toml: the fallback when no tag is found.
DEFAULT_TAG="0.1.0"
FAIL=0

pass() { printf 'PASS  %-10s %s\n' "$1" "$2"; }
fail() { printf 'FAIL  %-10s %s\n' "$1" "$2"; FAIL=1; }

echo "check_image: $(date '+%Y-%m-%d %H:%M:%S') host=$(hostname 2>/dev/null || echo ?) uid=$(id -u)"
echo "check_image: PATH=${PATH}"
echo "check_image: CONDA_PREFIX=${CONDA_PREFIX:-<unset>} expected=${EXPECTED}"

if [[ ! -r "${EXPECTED}" ]]; then
  echo "FAIL  expected versions file not readable: ${EXPECTED}"
  exit 1
fi
expected() { awk -F'\t' -v k="$1" '$1 == k { print $2 }' "${EXPECTED}"; }

# Version found in a tool's own output. `|| true`: several tools print their
# version in a usage message and exit non-zero.
found_version() {
  local tool="$1" out
  case "${tool}" in
    blat)     out="$(blat 2>&1 || true)";               sed -nE 's/.*BLAT v\. *([0-9][^ ]*).*/\1/p' <<<"${out}" | head -1 ;;
    bwa)      out="$(bwa 2>&1 || true)";                sed -nE 's/^Version: *([0-9][^ -]*).*/\1/p' <<<"${out}" | head -1 ;;
    samtools) out="$(samtools --version 2>&1 || true)"; sed -nE 's/^samtools ([0-9][^ ]*).*/\1/p' <<<"${out}" | head -1 ;;
    bcftools) out="$(bcftools --version 2>&1 || true)"; sed -nE 's/^bcftools ([0-9][^ ]*).*/\1/p' <<<"${out}" | head -1 ;;
    # The launcher prints its SIMD choice first; the version is the last line.
    bwa-mem2) out="$(bwa-mem2 version 2>&1 || true)";   grep -E '^[0-9]+\.[0-9][^ ]*$' <<<"${out}" | tail -1 ;;
    minimap2) out="$(minimap2 --version 2>&1 || true)"; sed -nE 's/^([0-9][^ -]*).*/\1/p' <<<"${out}" | head -1 ;;
    seqtk)    out="$(seqtk 2>&1 || true)";              sed -nE 's/^Version: *([0-9][^ -]*).*/\1/p' <<<"${out}" | head -1 ;;
    bedtools) out="$(bedtools --version 2>&1 || true)"; sed -nE 's/^bedtools v([0-9][^ ]*).*/\1/p' <<<"${out}" | head -1 ;;
    bowtie2)  out="$(bowtie2 --version 2>&1 || true)";  sed -nE 's/.* version ([0-9][^ ]*).*/\1/p' <<<"${out}" | head -1 ;;
  esac
}

for tool in blat bwa bwa-mem2 samtools bcftools minimap2 seqtk bedtools bowtie2; do
  want="$(expected "${tool}")"
  if ! path="$(command -v "${tool}")"; then
    fail "${tool}" "not on PATH"
    continue
  fi
  got="$(found_version "${tool}")"
  if [[ "${tool}" == bwa-mem2 && -n "${want}" && -n "${got}" ]]; then
    # Upstream: the 2.3 release still prints 2.2.1 (see header).
    pass "${tool}" "found=${got} locked=${want} (version not compared) (${path})"
  elif [[ -n "${want}" && "${got}" == "${want}" ]]; then
    pass "${tool}" "found=${got} expected=${want} (${path})"
  else
    fail "${tool}" "found=${got:-<none>} expected=${want:-<missing>} (${path})"
  fi
done

# RelocaTE3 version: derived from the tag the build saw.
base="$(expected relocate3_base)"
installed="$(expected relocate3_version)"
if version="$(relocaTE3 --version 2>&1)"; then
  if [[ -z "${base}" ]]; then
    fail relocaTE3 "found=${version}; no relocate3_base in ${EXPECTED}"
  elif [[ "${version}" == "0.0.0" || ( "${base}" != "${DEFAULT_TAG}" && "${version}" == "${DEFAULT_TAG}"* ) ]]; then
    fail relocaTE3 "found=${version}: placeholder/default-tag version"
  elif [[ "${version}" =~ [+.]d[0-9]{8}$ ]]; then
    fail relocaTE3 "found=${version}: dirty-tree build (versioningit date suffix)"
  elif ! [[ "${version}" == "${base}" || "${version}" == "${base}.post"* || "${version}" == "${base}+"* ]]; then
    fail relocaTE3 "found=${version} does not derive from tag_base=${base}"
  elif [[ -z "${installed}" || "${version}" != "${installed}" ]]; then
    fail relocaTE3 "found=${version} != installed at build=${installed:-<missing>}"
  else
    pass relocaTE3 "found=${version} tag_base=${base} installed=${installed}"
  fi
else
  fail relocaTE3 "relocaTE3 --version failed: ${version}"
fi

if path="$(command -v ps)" && ps -p $$ > /dev/null; then
  pass ps "$(ps --version 2>&1 | head -1) (${path})"
else
  fail ps "ps missing or not working"
fi

want_py="$(expected python)"
want_pysam="$(expected pysam)"
if out="$(python - "${want_py}" "${want_pysam}" 2>&1 <<'EOF'
import json
import os
import platform
import sys
from importlib import metadata
from pathlib import Path

import pysam
import RelocaTE3

want_py, want_pysam = sys.argv[1], sys.argv[2]
prefix = Path(os.environ.get("CONDA_PREFIX", sys.prefix)).resolve()
module = Path(RelocaTE3.__file__).resolve()
direct_url = json.loads(metadata.distribution("RelocaTE3").read_text("direct_url.json") or "{}")
problems = []
if prefix not in module.parents:
    problems.append(f"RelocaTE3 imported from {module}, not under {prefix}")
if direct_url.get("dir_info", {}).get("editable"):
    problems.append("RelocaTE3 is an editable install")
if platform.python_version() != want_py:
    problems.append(f"python {platform.python_version()} != locked {want_py}")
if pysam.__version__ != want_pysam:
    problems.append(f"pysam {pysam.__version__} != locked {want_pysam}")
print(f"python={platform.python_version()} pysam={pysam.__version__} RelocaTE3={RelocaTE3.__version__} at {module.parent}")
if problems:
    sys.exit("; ".join(problems))
EOF
)"; then
  pass python "${out}"
else
  fail python "${out}"
fi

if [[ "${FAIL}" -eq 0 ]]; then
  echo "check_image: PASS"
else
  echo "check_image: FAIL"
fi
exit "${FAIL}"
