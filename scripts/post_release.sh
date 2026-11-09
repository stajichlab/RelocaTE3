#!/usr/bin/env bash
# Move the repository back to a development version after a release tag.
#
#   bash scripts/post_release.sh NEXT_DEV_VERSION [PROJECT_DIR]
#   e.g. bash scripts/post_release.sh 0.3.0rc2.dev0      (after v0.3.0rc1)
#        bash scripts/post_release.sh 0.3.1.dev0         (after v0.3.0)
#
#   NEXT_DEV_VERSION  PEP 440 development version (".devN") that sorts AFTER
#                     the release in the manifest. "0.3.0dev" (= 0.3.0.dev0)
#                     sorts before 0.3.0rc1, so after an rc use the next rc's
#                     dev version.
#   PROJECT_DIR       default: this script's repository
#
# Run it on a new branch from main after the tag is pushed. It edits:
#   nextflow.config       manifest.version = NEXT_DEV_VERSION,
#                         params.container_tag = 'dev'
#   nextflow_schema.json  container_tag default "dev"
#   CHANGELOG.md          a new empty "## Unreleased" section on top
#   docs/source/pipeline.rst  parameter tables, regenerated from the schema
#                         (docs/source/_gen_params.py)
# (tests/nextflow_schema_test.py: a dev manifest defaults to the `dev` image,
# a release manifest to its own version.) The Python package version needs no
# edit: versioningit gives commits after the tag <release>.postN+g<sha>.
#
# Refuses (exit 1, nothing written) unless the manifest is a release whose
# container_tag defaults equal it, and NEXT_DEV_VERSION is a later dev version.
# A rerun with the same NEXT_DEV_VERSION is a no-op. tests/post_release_test.py.
#
# Also refuses before the release tag v<manifest.version> exists (the tag's CI
# build reads the release config), unless ALLOW_UNTAGGED=1.
#
# Env: PYTHON          interpreter >= 3.8 with `packaging` (default: the
#                      pipeline pixi env's python if present, else python3)
#      ALLOW_UNTAGGED  1 = skip the release-tag check
set -euo pipefail

usage() { sed -n '/^#   bash/,/^#   PROJECT_DIR/{s/^# \{0,1\}//;p;}' "${BASH_SOURCE[0]}"; }
if [[ $# -lt 1 || $# -gt 2 || "${1}" == -* ]]; then
  echo "usage: post_release.sh NEXT_DEV_VERSION [PROJECT_DIR]" >&2
  usage >&2
  exit 2
fi
NEXT="$1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${2:-$(dirname "${SCRIPT_DIR}")}"
if [[ -z "${PYTHON:-}" ]]; then
  if [[ -x "${PROJECT_DIR}/.pixi/envs/pipeline/bin/python" ]]; then
    PYTHON="${PROJECT_DIR}/.pixi/envs/pipeline/bin/python"
  else
    PYTHON="python3"
  fi
fi
for f in nextflow.config nextflow_schema.json CHANGELOG.md; do
  [[ -f "${PROJECT_DIR}/${f}" ]] || { echo "ERROR: no ${PROJECT_DIR}/${f}" >&2; exit 1; }
done

echo "post_release: $(date '+%Y-%m-%d %H:%M:%S') next=${NEXT} project=${PROJECT_DIR}"
# The embedded Python uses 3.8+ syntax; check first so an old system python3
# (3.6 on some HPCC nodes) gets this message instead of a SyntaxError.
if ! "${PYTHON}" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 8) else 1)' 2>/dev/null; then
  echo "ERROR: ${PYTHON} is not Python >= 3.8; set PYTHON (e.g. ${PROJECT_DIR}/.pixi/envs/pipeline/bin/python)" >&2
  exit 1
fi
export ALLOW_UNTAGGED="${ALLOW_UNTAGGED:-0}"

"${PYTHON}" - "${NEXT}" "${PROJECT_DIR}" <<'EOF'
from __future__ import annotations

import os
import re
import sys
import tempfile
from pathlib import Path

try:
    from packaging.version import InvalidVersion, Version
except ImportError:
    sys.exit("ERROR: python has no `packaging`; set PYTHON to the pipeline env's python")

next_text, project = sys.argv[1], Path(sys.argv[2])
config_path = project / "nextflow.config"
schema_path = project / "nextflow_schema.json"
log_path = project / "CHANGELOG.md"


def die(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def block_span(text: str, name: str) -> tuple[int, int]:
    """Start/end offsets of the top-level ``name { ... }`` block's body."""
    m = re.search(rf"(?m)^{name}\s*\{{\s*$", text)
    if not m:
        die(f"no top-level {name} {{ }} block in {config_path}")
    end = re.compile(r"(?m)^\}").search(text, m.end())
    return m.end(), end.start()


def sub_in_block(text: str, name: str, key: str, value: str | None = None):
    """Return (current value, new text) for KEY = '...' inside block NAME."""
    start, end = block_span(text, name)
    pat = re.compile(rf"(?m)^(\s*{key}\s*=\s*)'([^']*)'")
    m = pat.search(text, start, end)
    if not m:
        die(f"no {key} = '...' in the {name} block of {config_path}")
    if value is None:
        return m.group(2), text
    return m.group(2), text[: m.start()] + f"{m.group(1)}'{value}'" + text[m.end() :]


SCHEMA_TAG = re.compile(r'("container_tag": \{\s*"type": "string",\s*"default": )"([^"]*)"')

if next_text.startswith("v") or not re.fullmatch(r"[0-9][0-9A-Za-z.+!-]*", next_text):
    die(f"'{next_text}' is not a PEP 440 version (give e.g. 0.3.0rc2.dev0, without 'v')")
try:
    nxt = Version(next_text)
except InvalidVersion:
    die(f"'{next_text}' is not a PEP 440 version (give e.g. 0.3.0rc2.dev0)")
if not nxt.is_devrelease or nxt.local:
    die(f"'{next_text}' is not a development version (.devN, no local part), e.g. 0.3.0rc2.dev0")

config = config_path.read_text()
schema = schema_path.read_text()
changelog = log_path.read_text()
manifest, _ = sub_in_block(config, "manifest", "version")
tag, _ = sub_in_block(config, "params", "container_tag")
m = SCHEMA_TAG.search(schema)
if not m:
    die(f"no container_tag default in {schema_path}")
schema_tag = m.group(2)

if manifest == next_text and tag == "dev" and schema_tag == "dev":
    print(f"already at {next_text} with container_tag 'dev'; nothing to do")
    sys.exit(0)

RESTORE = ("git checkout -- nextflow.config nextflow_schema.json CHANGELOG.md "
           "docs/source/pipeline.rst")
if manifest == next_text and (tag != "dev" or schema_tag != "dev"):
    die(f"half-applied: manifest.version is already {next_text} but the container_tag "
        f"defaults are config '{tag}', schema '{schema_tag}'. Restore the release files "
        f"(in {project}: {RESTORE}) and run this again")

try:
    current = Version(manifest)
except InvalidVersion:
    die(f"manifest.version '{manifest}' is not a PEP 440 version")
if current.is_devrelease or current.is_postrelease or current.local:
    die(f"manifest.version {manifest} is not a release; run this only on the tagged release's "
        f"config (if an earlier run stopped half-way: {RESTORE})")
if tag != manifest or schema_tag != manifest:
    die(
        f"container_tag defaults (config '{tag}', schema '{schema_tag}') differ from "
        f"manifest.version {manifest}; this is not a release config"
    )
if os.environ.get("ALLOW_UNTAGGED") != "1":
    import subprocess

    probe = subprocess.run(
        ["git", "-C", str(project), "rev-parse", "-q", "--verify", f"refs/tags/v{manifest}"],
        capture_output=True, text=True,
    )
    if probe.returncode != 0:
        die(f"release tag v{manifest} not found in {project}; tag and push the release "
            f"first (git fetch --tags if it was pushed elsewhere), or set ALLOW_UNTAGGED=1")
if nxt <= current:
    die(f"{next_text} does not sort after the release {manifest} (PEP 440); e.g. use {current.base_version}rc2.dev0 "
        f"after an rc, or the next patch's .dev0 after a final release")

if str(nxt) != next_text:
    die(f"'{next_text}' is not in PEP 440 normal form; write {nxt}")

_, config = sub_in_block(config, "manifest", "version", next_text)
_, config = sub_in_block(config, "params", "container_tag", "dev")
schema = SCHEMA_TAG.sub(r'\g<1>"dev"', schema, count=1)
if not re.search(r"(?m)^## Unreleased\s*$", changelog):
    first = re.search(r"(?m)^## ", changelog)
    at = first.start() if first else len(changelog)
    changelog = changelog[:at] + "## Unreleased\n\n" + changelog[at:]


def write(path: Path, text: str) -> None:
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    with os.fdopen(fd, "w") as fh:
        fh.write(text)
    os.chmod(tmp, path.stat().st_mode & 0o777)
    os.replace(tmp, path)


for path, text in ((config_path, config), (schema_path, schema), (log_path, changelog)):
    write(path, text)
    print(f"updated {path}")
print(f"manifest.version {manifest} -> {next_text}; container_tag {manifest} -> dev")
EOF

GEN="${PROJECT_DIR}/docs/source/_gen_params.py"
if [[ -f "${GEN}" ]]; then
  "${PYTHON}" "${GEN}"
  "${PYTHON}" "${GEN}" --check
  echo "regenerated the parameter tables of ${PROJECT_DIR}/docs/source/pipeline.rst"
fi

cat <<EOF

Next steps (from ${PROJECT_DIR}):
  git checkout -b chore/post-release-${NEXT}        # if not already on a new branch
  git diff                                         # config, schema, CHANGELOG, pipeline.rst
  pixi run -e pipeline python -m pytest tests -p no:cacheprovider
  pixi run -e nf nf-test test && pixi run -e nf nextflow lint .
  git commit -am "Post-release: development version ${NEXT}"
  open a PR to main; after it merges, main's images are :dev again.
EOF
