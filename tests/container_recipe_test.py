"""The container recipe stays tied to pixi.lock and to what Nextflow needs.

The image is built only by CI (no Docker on the cluster), so these static checks
guard the recipe's load-bearing lines between CI runs:

- the build context excludes local environments and run/test output;
- the environment comes from ``pixi.lock`` unchanged (``--locked``);
- the pixi base image is pinned, never ``latest``;
- the conda environment sits at the same absolute path in both stages (conda
  prefixes are not relocatable);
- the runtime stage has ``ps`` (procps), which Nextflow's task metrics need,
  and puts the environment on ``PATH`` without an entrypoint.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCKERFILE = ROOT / "containers" / "Dockerfile"
DOCKERIGNORE = ROOT / ".dockerignore"
BUILD_ENV = ROOT / "containers" / "build_env.sh"
CHECK_IMAGE = ROOT / "containers" / "check_image.sh"
DRY_RUN = ROOT / "containers" / "dry_run_recipe.sh"


def _dockerignore_patterns() -> set[str]:
    patterns = set()
    for line in DOCKERIGNORE.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            patterns.add(line.rstrip("/"))
    return patterns


def _dockerfile_args() -> dict[str, str]:
    return dict(re.findall(r"^ARG\s+(\w+)=(\S+)\s*$", DOCKERFILE.read_text(), re.M))


def _stages() -> list[str]:
    """Dockerfile text split at each FROM (stage 0 = build, last = runtime)."""
    return re.split(r"^FROM\s", DOCKERFILE.read_text(), flags=re.M)[1:]


def test_dockerignore_excludes_local_state_but_not_git():
    patterns = _dockerignore_patterns()
    for path in (".pixi", "work", "results", "tests/data", ".nextflow*", ".nf-test", "logs"):
        assert path in patterns, f".dockerignore must exclude {path}"
    # versioningit needs the git history and tags.
    assert not any(p in {".git", ".git/*", "**/.git"} for p in patterns)


def test_pixi_base_image_is_pinned():
    args = _dockerfile_args()
    image = args["PIXI_IMAGE"]
    assert image.startswith("ghcr.io/prefix-dev/pixi:")
    tag = image.split(":", 1)[1]
    assert "latest" not in tag
    assert re.match(r"^\d+\.\d+\.\d+-", tag), f"pixi tag not pinned to a version: {tag}"
    assert "FROM ${PIXI_IMAGE} AS build" in DOCKERFILE.read_text()
    runtime = args["RUNTIME_IMAGE"]
    assert ":" in runtime and "latest" not in runtime


def test_environment_is_installed_locked_and_non_editable():
    script = BUILD_ENV.read_text()
    assert re.search(r"pixi install --locked -e blat\b", script)
    assert "--frozen" not in script
    # Replace the editable install with a normal one.
    assert re.search(r"pip uninstall -y RelocaTE3", script)
    assert re.search(r"pip install --no-deps\b", script)
    assert "pixi shell-hook -e blat" in script
    assert "pixi list -e blat --json --locked" in script


def test_env_copied_to_the_same_path_in_both_stages():
    args = _dockerfile_args()
    env_dir = args["ENV_DIR"]
    assert env_dir.startswith("/") and env_dir.endswith("/.pixi/envs/blat")
    build, *_, runtime = _stages()
    # A global ARG is only visible in a stage that re-declares it.
    assert re.search(r"^ARG ENV_DIR\s*$", build, re.M)
    assert re.search(r"^ARG ENV_DIR\s*$", runtime, re.M)
    assert re.search(r"^COPY --from=build \$\{ENV_DIR\} \$\{ENV_DIR\}\s*$", runtime, re.M)
    # PATH and CONDA_PREFIX point at that same prefix.
    assert re.search(r'^ENV PATH="\$\{ENV_DIR\}/bin:\$PATH"', runtime, re.M)
    assert re.search(r"^ENV CONDA_PREFIX=\$\{ENV_DIR\}", runtime, re.M)


def test_runtime_stage_has_procps_and_no_entrypoint():
    *_, runtime = _stages()
    assert re.search(r"apt-get install -y --no-install-recommends[^\n]*\bprocps\b", runtime)
    assert not re.search(r"^ENTRYPOINT", runtime, re.M)


def test_check_image_and_dry_run_are_executable():
    for script in (CHECK_IMAGE, BUILD_ENV, DRY_RUN):
        assert script.is_file(), script
        assert os.access(script, os.X_OK), f"{script} is not executable"
    dry = DRY_RUN.read_text()
    assert dry.startswith("#!/usr/bin/bash -l")
    assert "set -euo pipefail" in dry
    assert "#SBATCH -p intel" in dry


def _apt_packages(text: str) -> list[str]:
    return [
        m.split()
        for m in re.findall(r"install -y --no-install-recommends ([\w.+\- ]+?)\s*(?:\\|&&)", text)
    ]


def test_dry_run_installs_the_same_apt_packages_as_the_dockerfile():
    build, *_, runtime = _stages()
    dry = DRY_RUN.read_text()
    # build stage: git (versioningit) + a C/C++ toolchain: pixi.lock's blat env
    # carries pybedtools as a PyPI sdist, compiled during `pixi install`.
    assert _apt_packages(build) == [["git", "build-essential"]]
    assert _apt_packages(runtime) == [["procps"]]
    assert _apt_packages(dry) == _apt_packages(build) + _apt_packages(runtime)


# ------------------------------------------------- tools in the locked env


PIXI_LOCK = ROOT / "pixi.lock"


def _locked_conda_versions(env: str, platform: str = "linux-64") -> dict[str, str]:
    """{package: version} of the conda packages pixi.lock pins for ``env``."""
    lines = PIXI_LOCK.read_text().splitlines()
    start = lines.index(f"  {env}:")
    versions: dict[str, str] = {}
    in_platform = False
    for line in lines[start + 1 :]:
        if re.match(r"^ {0,2}\S", line):
            break
        if re.match(r"^      \S+:$", line):
            in_platform = line.strip() == f"{platform}:"
            continue
        m = re.match(r"^      - conda: \S+/([^/]+)-([^-/]+)-[^-/]+\.(?:conda|tar\.bz2)$", line)
        if in_platform and m:
            versions[m.group(1)] = m.group(2)
    return versions


def _build_env_tools() -> list[str]:
    m = re.search(r"for tool in \(([^)]*)\):", BUILD_ENV.read_text())
    assert m, "no expected-tools tuple in build_env.sh"
    return re.findall(r'"([\w-]+)"', m.group(1))


def _check_image_tools() -> list[str]:
    m = re.search(r"^for tool in ([\w -]+); do$", CHECK_IMAGE.read_text(), re.M)
    assert m, "no tool loop in check_image.sh"
    return m.group(1).split()


def test_blat_env_locks_bcftools_in_the_samtools_series():
    # Excision mode (characterize -x, pipeline --excision) shells out to
    # bcftools; it must be in the environment the image is built from.
    locked = _locked_conda_versions("blat")
    assert "samtools" in locked
    assert "bcftools" in locked, "bcftools missing from pixi.lock's blat env"
    series = lambda v: v.split(".")[:2]  # noqa: E731
    assert series(locked["bcftools"]) == series(locked["samtools"]) == series(locked["htslib"])


def test_image_checks_every_locked_tool():
    expected = _build_env_tools()
    checked = _check_image_tools()
    for tool in ("bcftools", "bwa-mem2"):
        assert tool in expected, f"{tool} not in build_env.sh's expected tools"
        assert tool in checked, f"{tool} not checked by check_image.sh"
    # Every expected tool (except the Python packages) is run by check_image.sh.
    assert set(expected) - {"pysam", "python"} == set(checked)
    locked = _locked_conda_versions("blat")
    for tool in checked:
        assert tool in locked, f"{tool} is checked but not locked in the blat env"
