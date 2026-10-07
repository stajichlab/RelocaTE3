#!/usr/bin/bash -l
# Executed by the process sampler inside the compute allocation.
set -euo pipefail
: "${PROFILE_OUTPUT:?}"
: "${ADAPTER_DIR:?}"
source "$ADAPTER_DIR/env.sh"

# Record the actual resolved tools after benchmark environment activation.
{
  date --iso-8601=seconds
  hostname
  python --version
  for tool in relocaTE3 blat bwa minimap2 samtools seqtk; do
    if command -v "$tool" >/dev/null 2>&1; then
      command -v "$tool"
      sha256sum "$(command -v "$tool")"
    else
      echo "$tool: absent (optional only for seqtk)"
      [[ "$tool" == seqtk ]] || exit 127
    fi
  done
  minimap2 --version
  samtools --version
  bwa 2>&1 | head -3 || true
  blat 2>&1 | head -3 || true
  python -c 'import os, pathlib, RelocaTE3; p = pathlib.Path(RelocaTE3.__file__).resolve(); print("source:", p); print("installed_metadata_version:", RelocaTE3.__version__); assert p.is_relative_to(pathlib.Path(os.environ["PYTHONPATH"]).resolve()), "Not importing frozen source"'
} > "$PROFILE_OUTPUT/tools.txt" 2>&1

exec bash "$ADAPTER_DIR/run.sh"
