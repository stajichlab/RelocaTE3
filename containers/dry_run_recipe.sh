#!/usr/bin/bash -l
#SBATCH -p intel
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH -J relocate3-container-dryrun
#SBATCH -o logs/container_dry_run.%j.log
#
# Dry run of containers/Dockerfile with Apptainer (this cluster has no Docker).
# It replays the recipe's commands, not the Dockerfile syntax (CI checks that):
#
#   build stage   sandbox from docker://<PIXI_IMAGE>; apt git; copy .git;
#                 containers/build_env.sh (the same script the Dockerfile RUNs)
#   runtime stage sandbox from docker://<RUNTIME_IMAGE>; apt procps; copy the
#                 env to the SAME path + meta/ + check_image.sh
#   check         check_image.sh in the runtime sandbox, read-only, as the
#                 calling (non-root) user, with only the Dockerfile's ENV
#                 (PATH, CONDA_PREFIX) on top of a clean environment.
#
# The context is a fresh clone of the project's HEAD (full history + tags, as
# actions/checkout fetch-depth: 0 gives CI), so uncommitted edits are NOT
# tested. Image tags and paths are read from the Dockerfile's ARG lines.
#
# Usage (from the project root):  sbatch containers/dry_run_recipe.sh
# Env: DRYRUN_DIR (default /scratch/$USER/<jobid>/relocate3-dryrun, node-local),
#      KEEP_DRYRUN=1 to keep the sandboxes for inspection.
#
# Root inside the sandboxes: `apptainer exec --fakeroot --writable`. Without
# /etc/subuid entries Apptainer uses a root-mapped user namespace (only uid 0
# exists), so apt's privilege drop to `_apt` cannot work; apt is told to stay
# root with APT::Sandbox::User=root (apt.conf(5)). The Dockerfile needs no such
# option (Docker builds run as real root in the build container).
set -euo pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
DOCKERFILE="${PROJECT_DIR}/containers/Dockerfile"
JOBTAG="${SLURM_JOB_ID:-$$}"
# HPCC creates a node-local /scratch/$USER/$SLURM_JOB_ID per job; /scratch/$USER
# itself is not writable on compute nodes.
DRYRUN_DIR="${DRYRUN_DIR:-/scratch/${USER}/${JOBTAG}/relocate3-dryrun}"
KEEP_DRYRUN="${KEEP_DRYRUN:-0}"

arg() { sed -nE "s/^ARG $1=(\S+)\s*$/\1/p" "${DOCKERFILE}" | head -1; }
PIXI_IMAGE="$(arg PIXI_IMAGE)"
RUNTIME_IMAGE="$(arg RUNTIME_IMAGE)"
IMG_PROJECT_DIR="$(arg PROJECT_DIR)"
ENV_DIR="$(arg ENV_DIR)"

echo "dry_run_recipe: $(date '+%Y-%m-%d %H:%M:%S') host=$(hostname) job=${JOBTAG}"
echo "  project_dir     = ${PROJECT_DIR}"
echo "  commit          = $(git -C "${PROJECT_DIR}" rev-parse HEAD)"
echo "  pixi_image      = ${PIXI_IMAGE}"
echo "  runtime_image   = ${RUNTIME_IMAGE}"
echo "  image project   = ${IMG_PROJECT_DIR}"
echo "  env_dir         = ${ENV_DIR}"
echo "  dryrun_dir      = ${DRYRUN_DIR} (keep=${KEEP_DRYRUN})"

for v in PIXI_IMAGE RUNTIME_IMAGE IMG_PROJECT_DIR ENV_DIR; do
  if [[ -z "${!v}" ]]; then
    echo "ERROR: could not read ${v} from ${DOCKERFILE}" >&2
    exit 1
  fi
done
if [[ -n "$(git -C "${PROJECT_DIR}" status --porcelain --untracked-files=no)" ]]; then
  echo "WARNING: uncommitted changes in ${PROJECT_DIR}; the dry run uses HEAD only" >&2
fi

module load apptainer/1.5.3
apptainer --version

if [[ -e "${DRYRUN_DIR}" ]]; then
  echo "ERROR: ${DRYRUN_DIR} exists; remove it or set DRYRUN_DIR" >&2
  exit 1
fi
mkdir -p "${DRYRUN_DIR}"/{cache,tmp}
export APPTAINER_CACHEDIR="${DRYRUN_DIR}/cache"
export APPTAINER_TMPDIR="${DRYRUN_DIR}/tmp"
cleanup() {
  if [[ "${KEEP_DRYRUN}" != "1" ]]; then
    chmod -R u+w "${DRYRUN_DIR}" 2>/dev/null || true
    rm -rf "${DRYRUN_DIR}"
    echo "cleaned up ${DRYRUN_DIR}"
  fi
}
trap cleanup EXIT

BUILD="${DRYRUN_DIR}/build"
RUNTIME="${DRYRUN_DIR}/runtime"
CTX="${DRYRUN_DIR}/ctx"
APT_ROOT=(-o APT::Sandbox::User=root)
# --no-home: keep $HOME (and its pixi/pip config) out of the container.
ROOT_EXEC=(apptainer exec --fakeroot --writable --no-home --cleanenv)

echo "== context: clone ${PROJECT_DIR} at HEAD"
git clone --quiet --no-hardlinks "${PROJECT_DIR}" "${CTX}"
git -C "${CTX}" checkout --quiet --detach "$(git -C "${PROJECT_DIR}" rev-parse HEAD)"
git -C "${CTX}" describe --tags

echo "== build stage: ${PIXI_IMAGE}"
apptainer build --sandbox "${BUILD}" "docker://${PIXI_IMAGE}"
"${ROOT_EXEC[@]}" "${BUILD}" bash -c \
  "apt-get ${APT_ROOT[*]} update && apt-get ${APT_ROOT[*]} install -y --no-install-recommends git build-essential && rm -rf /var/lib/apt/lists/*"
mkdir -p "${BUILD}${IMG_PROJECT_DIR}"
cp -a "${CTX}/.git" "${BUILD}${IMG_PROJECT_DIR}/.git"              # COPY .git
install -m 0755 "${CTX}/containers/build_env.sh" "${BUILD}/usr/local/bin/build_env.sh"  # COPY build_env.sh
"${ROOT_EXEC[@]}" --pwd "${IMG_PROJECT_DIR}" "${BUILD}" \
  bash /usr/local/bin/build_env.sh "${IMG_PROJECT_DIR}" "${ENV_DIR}"

echo "== runtime stage: ${RUNTIME_IMAGE}"
apptainer build --sandbox "${RUNTIME}" "docker://${RUNTIME_IMAGE}"
"${ROOT_EXEC[@]}" "${RUNTIME}" bash -c \
  "apt-get ${APT_ROOT[*]} update && apt-get ${APT_ROOT[*]} install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*"
mkdir -p "${RUNTIME}${ENV_DIR%/*}"
cp -a "${BUILD}${ENV_DIR}" "${RUNTIME}${ENV_DIR}"                        # COPY --from=build ENV_DIR ENV_DIR
cp -a "${BUILD}${IMG_PROJECT_DIR}/meta" "${RUNTIME}${IMG_PROJECT_DIR}/meta"
cp -a "${BUILD}${IMG_PROJECT_DIR}/containers/check_image.sh" "${RUNTIME}${IMG_PROJECT_DIR}/check_image.sh"
du -sh "${RUNTIME}${ENV_DIR}"
du -sh "${RUNTIME}"

echo "== check_image.sh (runtime sandbox, read-only, uid $(id -u))"
set +e
apptainer exec --no-home --cleanenv \
  --env "PATH=${ENV_DIR}/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
  --env "CONDA_PREFIX=${ENV_DIR}" \
  "${RUNTIME}" bash "${IMG_PROJECT_DIR}/check_image.sh"
RC=$?
set -e

if [[ "${RC}" -eq 0 ]]; then
  echo "DRY RUN: PASS"
else
  echo "DRY RUN: FAIL (check_image.sh exit ${RC})"
fi
echo "dry_run_recipe: done $(date '+%Y-%m-%d %H:%M:%S')"
exit "${RC}"
