"""A release image carries one version everywhere.

* ``containers/build_env.sh --version-only`` resolves the tag versioningit will
  describe. versioningit runs ``git describe --long --dirty --always --tags``;
  with two annotated tags on one commit (an rc and the final release) git picks
  the newer one, so an rc build could get the final version. ``RELEASE_TAG``
  (the tag being built, passed by .github/workflows/container.yml) wins: the
  other tags at HEAD are deleted from the build's own ``.git`` copy.
* ``containers/check_release.sh`` (release tags only, before push) fails unless
  the tag version, the image's ``relocaTE3 --version``, ``manifest.version``
  and the default ``params.container_tag`` are all equal.
* ``containers/check_image.sh`` compares ``relocaTE3 --version`` with the exact
  version installed at build time and rejects versioningit's dirty suffix.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
BUILD_ENV = ROOT / "containers" / "build_env.sh"
CHECK_RELEASE = ROOT / "containers" / "check_release.sh"
CHECK_IMAGE = ROOT / "containers" / "check_image.sh"
# versioningit 3.x, git method (versioningit/git.py GitRepo.describe).
VERSIONINGIT_DESCRIBE = ["git", "describe", "--long", "--dirty", "--always", "--tags"]


def _git_env(tmp_path: Path) -> dict[str, str]:
    env = dict(os.environ)
    gitconfig = tmp_path / "gitconfig"
    gitconfig.touch()
    env.update(
        GIT_CONFIG_GLOBAL=str(gitconfig),  # build_env.sh adds safe.directory
        GIT_CONFIG_NOSYSTEM="1",
        GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@t", GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@t",
    )
    env.pop("RELEASE_TAG", None)
    return env


def _repo(tmp_path: Path, tags: list[str]) -> tuple[Path, dict[str, str]]:
    """A repo with one commit and annotated ``tags`` on it, oldest first."""
    env = _git_env(tmp_path)
    repo = tmp_path / "project"
    repo.mkdir()
    subprocess.run(["git", "init", "-q"], cwd=repo, env=env, check=True)
    (repo / "f").write_text("x\n")
    subprocess.run(["git", "add", "f"], cwd=repo, env=env, check=True)
    subprocess.run(["git", "commit", "-q", "-m", "c"], cwd=repo, env=env, check=True)
    for day, tag in enumerate(tags, start=1):
        env_t = dict(env, GIT_COMMITTER_DATE=f"2026-11-0{day}T00:00:00")
        subprocess.run(["git", "tag", "-a", "-m", tag, tag], cwd=repo, env=env_t, check=True)
    return repo, env


def _version_only(repo: Path, env: dict[str, str], release_tag: str | None = None):
    if release_tag is not None:
        env = dict(env, RELEASE_TAG=release_tag)
    return subprocess.run(
        ["bash", str(BUILD_ENV), "--version-only", str(repo)],
        env=env, capture_output=True, text=True,
    )


def _describe(repo: Path, env: dict[str, str]) -> str:
    return subprocess.run(VERSIONINGIT_DESCRIBE, cwd=repo, env=env, capture_output=True, text=True, check=True).stdout.strip()


def _base(proc) -> str:
    lines = [l for l in proc.stdout.splitlines() if l.startswith("base_version=")]
    assert lines, proc.stdout + proc.stderr
    return lines[-1].split("=", 1)[1]


def test_rc_and_final_on_one_commit_without_release_tag_follows_versioningit(tmp_path):
    repo, env = _repo(tmp_path, ["v0.3.0", "v0.3.0rc1"])  # rc tagged later: git prefers it
    proc = _version_only(repo, env)
    assert proc.returncode == 0, proc.stderr
    tag = _describe(repo, env).rsplit("-", 2)[0]
    assert _base(proc) == tag.removeprefix("v")


@pytest.mark.parametrize("wanted, other", [("v0.3.0", "v0.3.0rc1"), ("v0.3.0rc1", "v0.3.0")])
def test_release_tag_wins_when_several_tags_point_at_head(tmp_path, wanted, other):
    # `other` is created last, so plain `git describe` would pick it.
    repo, env = _repo(tmp_path, [wanted, other])
    assert _describe(repo, env).startswith(other + "-0-")
    proc = _version_only(repo, env, release_tag=wanted)
    assert proc.returncode == 0, proc.stderr
    assert _base(proc) == wanted.removeprefix("v")
    # What versioningit will run now names the release tag.
    assert _describe(repo, env).startswith(wanted + "-0-")


def test_release_tag_not_at_head_is_an_error(tmp_path):
    repo, env = _repo(tmp_path, ["v0.3.0"])
    subprocess.run(["git", "commit", "-q", "--allow-empty", "-m", "later"], cwd=repo, env=env, check=True)
    proc = _version_only(repo, env, release_tag="v0.3.0")
    assert proc.returncode != 0
    assert "does not point at HEAD" in proc.stderr


def test_no_tag_is_an_error(tmp_path):
    repo, env = _repo(tmp_path, [])
    proc = _version_only(repo, env)
    assert proc.returncode != 0
    assert "no git tag" in proc.stderr


# ----------------------------------------------------------- check_release.sh


def _nf_config(tmp_path: Path, manifest: str, tag: str) -> Path:
    cfg = tmp_path / "nextflow.config"
    cfg.write_text(
        "params {\n    outdir = null\n"
        f"    container_registry = 'ghcr.io/stajichlab/relocate3'\n    container_tag = '{tag}'\n}}\n\n"
        "manifest {\n    name = 'stajichlab/RelocaTE3'\n"
        f"    nextflowVersion = '!>=26.04.0'\n    version = '{manifest}'\n}}\n"
    )
    return cfg


def _check_release(expected: str, image: str, cfg: Path):
    return subprocess.run(
        ["bash", str(CHECK_RELEASE), expected, image, str(cfg)], capture_output=True, text=True
    )


def test_check_release_passes_when_everything_agrees(tmp_path):
    proc = _check_release("0.3.0", "0.3.0", _nf_config(tmp_path, "0.3.0", "0.3.0"))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "RELEASE CHECK: PASS" in proc.stdout


@pytest.mark.parametrize(
    "image, manifest, tag, bad",
    [
        ("0.3.0rc1", "0.3.0", "0.3.0", "relocaTE3 --version"),
        ("0.3.0.post1+gabc", "0.3.0", "0.3.0", "relocaTE3 --version"),
        ("0.3.0", "0.3.0dev", "0.3.0", "manifest.version"),
        ("0.3.0", "0.3.0", "dev", "params.container_tag"),
    ],
)
def test_check_release_names_the_mismatch(tmp_path, image, manifest, tag, bad):
    proc = _check_release("0.3.0", image, _nf_config(tmp_path, manifest, tag))
    assert proc.returncode != 0
    assert "RELEASE CHECK: FAIL" in proc.stdout
    fails = [l for l in proc.stdout.splitlines() if l.startswith("FAIL")]
    assert len(fails) == 1 and bad in fails[0], proc.stdout


def test_check_release_reads_the_real_config():
    # Today's manifest is a dev version, so a release check must fail on it.
    proc = _check_release("0.3.0", "0.3.0", ROOT / "nextflow.config")
    assert proc.returncode != 0
    assert "manifest.version" in proc.stdout


# ------------------------------------------------------------- check_image.sh


def _check_image_version_line(tmp_path: Path, reported: str, installed: str, base: str) -> str:
    bindir = tmp_path / "bin"
    bindir.mkdir()
    stub = bindir / "relocaTE3"
    stub.write_text(f"#!/usr/bin/env bash\necho '{reported}'\n")
    stub.chmod(0o755)
    expected = tmp_path / "expected_versions.tsv"
    expected.write_text(f"relocate3_base\t{base}\nrelocate3_version\t{installed}\n")
    proc = subprocess.run(
        ["bash", str(CHECK_IMAGE)],
        env={"PATH": f"{bindir}:/usr/bin:/bin", "EXPECTED_VERSIONS": str(expected), "HOME": str(tmp_path)},
        capture_output=True, text=True,
    )
    lines = [l for l in proc.stdout.splitlines() if "relocaTE3" in l.split()[1:2]]
    assert len(lines) == 1, proc.stdout
    return lines[0]


def test_check_image_accepts_the_installed_version(tmp_path):
    v = "0.3.0rc1.post2+g1234abc"
    assert _check_image_version_line(tmp_path, v, v, "0.3.0rc1").startswith("PASS")


def test_check_image_rejects_a_version_other_than_the_installed_one(tmp_path):
    line = _check_image_version_line(tmp_path, "0.3.0rc1.post3+g9999999", "0.3.0rc1.post2+g1234abc", "0.3.0rc1")
    assert line.startswith("FAIL")


@pytest.mark.parametrize("v", ["0.3.0rc1+d20261104", "0.3.0rc1.post2+g1234abc.d20261104"])
def test_check_image_rejects_a_dirty_build(tmp_path, v):
    line = _check_image_version_line(tmp_path, v, v, "0.3.0rc1")
    assert line.startswith("FAIL") and "dirty" in line
