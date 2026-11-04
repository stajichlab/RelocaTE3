#!/usr/bin/env bash
# Build-stage commands of containers/Dockerfile, kept in one script so the
# Apptainer dry run (containers/dry_run_recipe.sh) runs exactly the same steps.
#
# Usage: build_env.sh <project_dir> <env_dir>
#   project_dir  holds a copy of the repository's .git (and nothing else needed)
#   env_dir      the absolute prefix the runtime stage copies; must equal
#                <project_dir>/.pixi/envs/blat
#
# Requires: pixi and git on PATH, network access (conda channels + PyPI).
#
# Steps:
#   1. Check out the committed tree from .git, so versioningit sees a clean
#      repository (a partial COPY of tracked files would mark it dirty and add
#      a build-date suffix to the version).
#   2. pixi install --locked -e blat: the environment exactly as in pixi.lock.
#   3. Record the activation (pixi shell-hook) and the locked tool versions,
#      before step 4: shell-hook may re-sync the environment, which would bring
#      the editable install back.
#   4. Replace pixi's editable RelocaTE3 (a .pth pointing at <project_dir>/src)
#      with a normal site-packages install, so the runtime image needs no source.
#   5. Verify the install is non-editable and carries the git-derived version.
set -euo pipefail

PROJECT_DIR="${1:?usage: build_env.sh <project_dir> <env_dir>}"
ENV_DIR="${2:?usage: build_env.sh <project_dir> <env_dir>}"
META_DIR="${PROJECT_DIR}/meta"
PY="${ENV_DIR}/bin/python"

echo "[build_env] $(date '+%Y-%m-%d %H:%M:%S')"
echo "[build_env] project_dir=${PROJECT_DIR} env_dir=${ENV_DIR} pixi=$(pixi --version)"

if [[ "${ENV_DIR}" != "${PROJECT_DIR}/.pixi/envs/blat" ]]; then
  echo "ERROR: env_dir must be ${PROJECT_DIR}/.pixi/envs/blat (got ${ENV_DIR})" >&2
  exit 1
fi
if [[ ! -d "${PROJECT_DIR}/.git" ]]; then
  echo "ERROR: ${PROJECT_DIR}/.git is missing (versioningit needs it)" >&2
  exit 1
fi

cd "${PROJECT_DIR}"

# --- 1. clean checkout --------------------------------------------------------
# The files are owned by whoever built the image layer; trust the directory so
# git does not refuse with "dubious ownership" (git-config(1), safe.directory).
git config --global --add safe.directory "${PROJECT_DIR}"
git reset --quiet --hard HEAD
if [[ -n "$(git status --porcelain --untracked-files=no)" ]]; then
  echo "ERROR: working tree not clean after checkout" >&2
  git status --short >&2
  exit 1
fi
# versioningit falls back to default-tag (0.1.0) when git describe finds no tag,
# e.g. in a shallow clone. Fail instead of shipping a fake version.
if ! TAG="$(git describe --tags --abbrev=0 2>/dev/null)"; then
  echo "ERROR: no git tag reachable from HEAD; clone with full history and tags" >&2
  echo "       (actions/checkout: fetch-depth: 0)" >&2
  exit 1
fi
# versioningit's default tag2version step strips a leading "v".
BASE_VERSION="${TAG#v}"
echo "[build_env] commit=$(git rev-parse HEAD) describe=$(git describe --tags) base_version=${BASE_VERSION}"

# --- 2. locked environment ----------------------------------------------------
pixi install --locked -e blat

# --- 3. activation + locked versions ------------------------------------------
mkdir -p "${META_DIR}"
pixi shell-hook -e blat --locked -s bash > "${META_DIR}/shell-hook.sh"
pixi list -e blat --json --locked > "${META_DIR}/pixi-list-blat.json"
"${PY}" - "${META_DIR}/pixi-list-blat.json" "${BASE_VERSION}" > "${META_DIR}/expected_versions.tsv" <<'EOF'
import json
import sys

packages = json.load(open(sys.argv[1]))
conda = {p["name"]: p["version"] for p in packages if p.get("kind") == "conda"}
for tool in ("blat", "bwa", "bwa-mem2", "samtools", "bcftools", "minimap2", "seqtk", "bedtools", "bowtie2", "pysam", "python"):
    if tool not in conda:
        sys.exit(f"ERROR: {tool} not in the locked blat environment")
    print(f"{tool}\t{conda[tool]}")
print(f"relocate3_base\t{sys.argv[2]}")
EOF
echo "[build_env] expected versions:"
cat "${META_DIR}/expected_versions.tsv"

# --- 4. non-editable RelocaTE3 ------------------------------------------------
# Call the env's python directly, not `pixi run`: pixi would re-sync the env and
# restore the editable install. --no-deps: every dependency is already in the
# locked environment; pip must not change it. The build backend (setuptools +
# versioningit, pyproject [build-system]) runs in pip's isolated build env.
"${PY}" -m pip uninstall -y RelocaTE3
"${PY}" -m pip install --no-deps --no-cache-dir "${PROJECT_DIR}"

# --- 5. verify ------------------------------------------------------------------
"${PY}" - "${ENV_DIR}" "${BASE_VERSION}" <<'EOF'
import json
import sys
from importlib import metadata
from pathlib import Path

env_dir, base = Path(sys.argv[1]), sys.argv[2]
dist = metadata.distribution("RelocaTE3")
direct_url = json.loads(dist.read_text("direct_url.json") or "{}")
if direct_url.get("dir_info", {}).get("editable"):
    sys.exit(f"ERROR: RelocaTE3 is still editable: {direct_url}")
site = env_dir / "lib"
editable_pth = list(site.glob("python3*/site-packages/__editable__*elocate*"))
if editable_pth:
    sys.exit(f"ERROR: editable hook left behind: {editable_pth}")
import RelocaTE3

module = Path(RelocaTE3.__file__).resolve()
if env_dir.resolve() not in module.parents:
    sys.exit(f"ERROR: RelocaTE3 imported from {module}, not from {env_dir}")
version = dist.version
if not (version == base or version.startswith((base + ".post", base + "+"))):
    sys.exit(f"ERROR: version {version} does not derive from tag {base}")
print(f"[build_env] RelocaTE3 {version} non-editable at {module.parent}")
EOF
echo "[build_env] done $(date '+%Y-%m-%d %H:%M:%S')"
