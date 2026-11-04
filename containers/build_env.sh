#!/usr/bin/env bash
# Build-stage commands of containers/Dockerfile, kept in one script so the
# Apptainer dry run (containers/dry_run_recipe.sh) runs exactly the same steps.
#
# Usage: build_env.sh <project_dir> <env_dir>
#        build_env.sh --version-only <project_dir>   (step 1 only; prints
#                     base_version=<v>; tests/release_version_test.py)
#   project_dir  holds a copy of the repository's .git (and nothing else needed)
#   env_dir      the absolute prefix the runtime stage copies; must equal
#                <project_dir>/.pixi/envs/blat
#
# Env: RELEASE_TAG  the git tag being released (v*), set by the Dockerfile's
#      build-arg on tag builds (.github/workflows/container.yml); empty
#      otherwise. versioningit 3.x runs `git describe --long --dirty --always
#      --tags`; with two annotated tags on one commit (an rc and the final
#      release) git picks the newer one, so the rc's build would get the final
#      version. With RELEASE_TAG set, the other tags at HEAD are deleted from
#      THIS .git copy (the image build's own), the RELEASE_TAG must point at
#      HEAD, and the installed version must equal it exactly.
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
#   5. Verify the install is non-editable and carries the git-derived version
#      (exactly RELEASE_TAG's on a release, never versioningit's dirty
#      suffix); record it as relocate3_version for check_image.sh.
set -euo pipefail

VERSION_ONLY=0
if [[ "${1:-}" == "--version-only" ]]; then
  VERSION_ONLY=1
  PROJECT_DIR="${2:?usage: build_env.sh --version-only <project_dir>}"
  ENV_DIR=""
else
  PROJECT_DIR="${1:?usage: build_env.sh <project_dir> <env_dir>}"
  ENV_DIR="${2:?usage: build_env.sh <project_dir> <env_dir>}"
fi
RELEASE_TAG="${RELEASE_TAG:-}"
META_DIR="${PROJECT_DIR}/meta"
PY="${ENV_DIR}/bin/python"

echo "[build_env] $(date '+%Y-%m-%d %H:%M:%S')"
echo "[build_env] project_dir=${PROJECT_DIR} env_dir=${ENV_DIR:-<version only>} release_tag=${RELEASE_TAG:-<none>}"

if [[ "${VERSION_ONLY}" -eq 0 && "${ENV_DIR}" != "${PROJECT_DIR}/.pixi/envs/blat" ]]; then
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
# A release builds exactly its own tag (see RELEASE_TAG in the header).
if [[ -n "${RELEASE_TAG}" ]]; then
  if [[ "$(git rev-parse -q --verify "refs/tags/${RELEASE_TAG}^{commit}" || true)" != "$(git rev-parse HEAD)" ]]; then
    echo "ERROR: RELEASE_TAG ${RELEASE_TAG} does not point at HEAD ($(git rev-parse HEAD))" >&2
    exit 1
  fi
  while read -r other; do
    if [[ -n "${other}" && "${other}" != "${RELEASE_TAG}" ]]; then
      git tag -d "${other}" > /dev/null
      echo "[build_env] RELEASE_TAG=${RELEASE_TAG}: removed tag ${other} (also at HEAD) from this build's .git"
    fi
  done < <(git tag --points-at HEAD)
fi
# The tag versioningit will use: its own describe command. --always prints a
# bare hash when no tag is reachable (e.g. a shallow clone), and versioningit
# then falls back to default-tag (0.1.0): fail instead of shipping that.
DESCRIBE="$(git describe --long --dirty --always --tags)"
if [[ ! "${DESCRIBE}" =~ ^(.+)-[0-9]+-g[0-9a-f]+(-dirty)?$ ]]; then
  echo "ERROR: no git tag reachable from HEAD (describe: ${DESCRIBE}); clone with full" >&2
  echo "       history and tags (actions/checkout: fetch-depth: 0)" >&2
  exit 1
fi
TAG="${BASH_REMATCH[1]}"
# versioningit's default tag2version step strips a leading "v".
BASE_VERSION="${TAG#v}"
echo "[build_env] commit=$(git rev-parse HEAD) describe=${DESCRIBE} base_version=${BASE_VERSION}"
if [[ "${VERSION_ONLY}" -eq 1 ]]; then
  echo "base_version=${BASE_VERSION}"
  exit 0
fi
echo "[build_env] pixi=$(pixi --version)"

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
"${PY}" - "${ENV_DIR}" "${BASE_VERSION}" "${RELEASE_TAG}" "${META_DIR}/expected_versions.tsv" <<'EOF'
import json
import re
import sys
from importlib import metadata
from pathlib import Path

env_dir, base, release_tag, expected = Path(sys.argv[1]), sys.argv[2], sys.argv[3], Path(sys.argv[4])
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
# versioningit's dirty formats end in +d<YYYYMMDD> or .d<YYYYMMDD>.
if re.search(r"[+.]d\d{8}$", version):
    sys.exit(f"ERROR: version {version} is a dirty-tree build")
if release_tag and version != base:
    sys.exit(f"ERROR: release {release_tag}: version {version} != {base}")
with expected.open("a") as fh:
    fh.write(f"relocate3_version\t{version}\n")
print(f"[build_env] RelocaTE3 {version} non-editable at {module.parent}")
EOF
echo "[build_env] done $(date '+%Y-%m-%d %H:%M:%S')"
