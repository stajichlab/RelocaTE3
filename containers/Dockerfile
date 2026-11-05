# RelocaTE3 runtime image: pixi.lock's `blat` environment (the validated tool
# versions) plus a non-editable install of RelocaTE3. Built only by CI.
#
#   docker build -f containers/Dockerfile -t relocate3 .
#   docker run --rm relocate3 bash /opt/relocate3/check_image.sh
#
# The build context must hold the full git history and tags (versioningit;
# actions/checkout with fetch-depth: 0). The image is built from the COMMITTED
# tree in .git (containers/build_env.sh), not from uncommitted edits.
#
# Pattern: pixi's container docs (https://pixi.sh/latest/deployment/container/),
# a multi-stage build with `pixi install --locked`, `pixi shell-hook`, and the
# env copied to the same prefix in a plain runtime image.

# Pixi image tags: https://github.com/prefix-dev/pixi-docker (build matrix in
# .github/workflows/build.yml): <pixi version>-<base>. 0.70.2 is the pixi that
# wrote pixi.lock; `noble` is Ubuntu 24.04.
ARG PIXI_IMAGE=ghcr.io/prefix-dev/pixi:0.70.2-noble
# Same distribution as the pixi image, so the glibc the env was installed
# against is the glibc it runs on.
ARG RUNTIME_IMAGE=ubuntu:24.04
# Conda prefixes are not relocatable: the env must live at this same absolute
# path in both stages. pixi installs to <project>/.pixi/envs/<env>.
ARG PROJECT_DIR=/opt/relocate3
ARG ENV_DIR=/opt/relocate3/.pixi/envs/blat

FROM ${PIXI_IMAGE} AS build
ARG PROJECT_DIR
ARG ENV_DIR
# The git tag being released; empty except on v* tag builds (container.yml).
# build_env.sh then versions the image from exactly this tag.
ARG RELEASE_TAG=
# The pixi image ships only pixi + ca-certificates. Build-stage only: git, for
# versioningit's "git" method (`git describe`). No compiler: every package in
# pixi.lock's blat env is a binary (conda or wheel).
# Keep in sync with containers/dry_run_recipe.sh (tests/container_recipe_test.py).
RUN apt-get update \
    && apt-get install -y --no-install-recommends git \
    && rm -rf /var/lib/apt/lists/*
WORKDIR ${PROJECT_DIR}
COPY .git ${PROJECT_DIR}/.git
COPY containers/build_env.sh /usr/local/bin/build_env.sh
RUN RELEASE_TAG="${RELEASE_TAG}" bash /usr/local/bin/build_env.sh "${PROJECT_DIR}" "${ENV_DIR}"

FROM ${RUNTIME_IMAGE} AS runtime
ARG PROJECT_DIR
ARG ENV_DIR
LABEL org.opencontainers.image.source="https://github.com/stajichlab/RelocaTE3"
LABEL org.opencontainers.image.description="RelocaTE3 with the pixi.lock blat environment"
LABEL org.opencontainers.image.licenses="MIT"
# procps provides `ps`, which Nextflow's task metrics (.command.run) need.
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build ${ENV_DIR} ${ENV_DIR}
# shell-hook.sh (pixi activation, for reference), expected_versions.tsv and the
# raw `pixi list` of the locked env, read by check_image.sh.
COPY --from=build ${PROJECT_DIR}/meta ${PROJECT_DIR}/meta
COPY --from=build ${PROJECT_DIR}/containers/check_image.sh ${PROJECT_DIR}/check_image.sh
# No ENTRYPOINT: Nextflow runs `bash .command.run` and Apptainer ignores it.
# The blat env has no activate.d scripts, so PATH + CONDA_PREFIX is the whole
# activation (see meta/shell-hook.sh).
ENV PATH="${ENV_DIR}/bin:$PATH"
ENV CONDA_PREFIX=${ENV_DIR}
# Apptainer/Singularity bind the host $HOME: keep ~/.local site-packages out of
# the environment's Python.
ENV PYTHONNOUSERSITE=1
