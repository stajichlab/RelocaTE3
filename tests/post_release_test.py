"""scripts/post_release.sh: move main back to a development version after a tag.

After ``v0.3.0rc1`` is tagged, ``main`` must not keep claiming 0.3.0rc1 (and its
``:0.3.0rc1`` image) for every later commit. The next development version must
sort *after* the release in PEP 440 (``0.3.0dev`` = ``0.3.0.dev0`` sorts before
``0.3.0rc1``), so the convention is ``<next pre-release>.dev0``, e.g.
``0.3.0rc2.dev0``, with the default ``container_tag`` back to ``dev``
(tests/nextflow_schema_test.py). The script is run on a temp copy of the real
files here; it must never touch the repository.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "post_release.sh"
FILES = (
    "nextflow.config",
    "nextflow_schema.json",
    "CHANGELOG.md",
    "docs/source/pipeline.rst",  # parameter tables generated from the schema
    "docs/source/_gen_params.py",
)
GEN = "docs/source/_gen_params.py"


def _set_release(project: Path, version: str) -> None:
    """Put the copy in the state a release commit has (manifest = tag = version)."""
    cfg = project / "nextflow.config"
    text = cfg.read_text()
    text = re.sub(
        r"(?m)^(\s*version\s*=\s*)'[^']*'", rf"\g<1>'{version}'", text, count=1
    )
    text = re.sub(
        r"(?m)^(\s*container_tag\s*=\s*)'[^']*'", rf"\g<1>'{version}'", text, count=1
    )
    cfg.write_text(text)
    schema = project / "nextflow_schema.json"
    data = schema.read_text()
    data = re.sub(
        r'("container_tag": \{\s*"type": "string",\s*"default": )"[^"]*"',
        rf'\g<1>"{version}"',
        data,
    )
    schema.write_text(data)
    subprocess.run(
        [sys.executable, str(project / GEN)], check=True, capture_output=True
    )
    log = project / "CHANGELOG.md"
    lines = log.read_text().splitlines(keepends=True)
    # A release commit has no Unreleased heading: its section is the version's.
    log.write_text("".join(ln for ln in lines if ln.strip() != "## Unreleased"))


@pytest.fixture
def project(tmp_path: Path) -> Path:
    proj = tmp_path / "proj"
    proj.mkdir()
    for name in FILES:
        (proj / name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / name, proj / name)
    _set_release(proj, "0.3.0rc1")
    _git_tag(proj, "v0.3.0rc1")
    return proj


GIT = ["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "init.defaultBranch=main"]


def _git_tag(project: Path, tag: str) -> None:
    """Make the copy a git repo whose HEAD carries the release tag, as on main."""
    for cmd in (["init", "-q"], ["add", "-A"], ["commit", "-qm", "release"], ["tag", tag]):
        subprocess.run([*GIT, "-C", str(project), *cmd], check=True, capture_output=True)


def _run(*args: str, **extra_env: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["bash", str(SCRIPT), *args],
        capture_output=True,
        text=True,
        check=False,
        env={
            "PATH": f"{Path(sys.executable).parent}:/usr/bin:/bin",
            "PYTHON": sys.executable,
            **extra_env,
        },
    )


def _manifest(project: Path) -> str:
    text = (project / "nextflow.config").read_text()
    block = text[text.index("manifest {") :]
    return re.search(r"(?m)^\s*version\s*=\s*'([^']*)'", block).group(1)


def _container_tag(project: Path) -> str:
    return re.search(
        r"(?m)^\s*container_tag\s*=\s*'([^']*)'",
        (project / "nextflow.config").read_text(),
    ).group(1)


def _schema_tag(project: Path) -> str:
    schema = json.loads((project / "nextflow_schema.json").read_text())
    for group in schema["$defs"].values():
        if "container_tag" in group.get("properties", {}):
            return group["properties"]["container_tag"]["default"]
    raise AssertionError("container_tag not in the schema")


def _snapshot(project: Path) -> dict[str, str]:
    return {name: (project / name).read_text() for name in FILES}


def test_bumps_to_next_dev_version(project: Path):
    before = _snapshot(project)
    res = _run("0.3.0rc2.dev0", str(project))
    assert res.returncode == 0, res.stderr
    assert _manifest(project) == "0.3.0rc2.dev0"
    assert _container_tag(project) == "dev"
    assert _schema_tag(project) == "dev"
    # The schema stays valid JSON and keeps its layout: one line changed.
    after = _snapshot(project)
    diff = [
        (a, b)
        for a, b in zip(
            before["nextflow_schema.json"].splitlines(),
            after["nextflow_schema.json"].splitlines(),
            strict=True,
        )
        if a != b
    ]
    assert diff == [
        (
            '                    "default": "0.3.0rc1",',
            '                    "default": "dev",',
        )
    ]
    cfg_diff = [
        b
        for a, b in zip(
            before["nextflow.config"].splitlines(),
            after["nextflow.config"].splitlines(),
            strict=True,
        )
        if a != b
    ]
    assert [ln.strip() for ln in cfg_diff] == [
        "container_tag = 'dev'",
        "version = '0.3.0rc2.dev0'",
    ]
    # CHANGELOG gets an empty Unreleased section above the release's.
    log = after["CHANGELOG.md"]
    first_release = re.search(r"(?m)^## .*$", before["CHANGELOG.md"]).group(0)
    assert log.index("## Unreleased\n\n" + first_release) < log.index(first_release)
    assert log.replace("## Unreleased\n\n", "", 1) == before["CHANGELOG.md"]
    # The guide's generated parameter table follows the schema.
    gen = subprocess.run(
        [sys.executable, str(project / GEN), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert gen.returncode == 0, gen.stdout + gen.stderr
    assert "``0.3.0rc1``" in before["docs/source/pipeline.rst"]
    assert (
        "``0.3.0rc1``"
        not in after["docs/source/pipeline.rst"]
        .split("BEGIN GENERATED")[1]
        .split("END GENERATED")[0]
    )
    assert "next steps" in res.stdout.lower()
    assert "git checkout -b" in res.stdout


def test_rerun_is_a_no_op(project: Path):
    assert _run("0.3.0rc2.dev0", str(project)).returncode == 0
    before = _snapshot(project)
    res = _run("0.3.0rc2.dev0", str(project))
    assert res.returncode == 0, res.stderr
    assert "already" in res.stdout
    assert _snapshot(project) == before


@pytest.mark.parametrize(
    "version, why",
    [
        ("0.3.0rc2", "development"),  # not a .devN version
        ("0.3.0.dev0", "after"),  # sorts before 0.3.0rc1
        ("0.3.0dev", "after"),  # same as 0.3.0.dev0
        ("0.3.0rc1.dev0", "after"),  # before 0.3.0rc1 too
        ("v0.3.0rc2.dev0", "PEP 440"),  # a tag name, not a version
        ("banana", "PEP 440"),
    ],
)
def test_rejects_bad_versions(project: Path, version: str, why: str):
    before = _snapshot(project)
    res = _run(version, str(project))
    assert res.returncode != 0
    assert why in res.stderr
    assert _snapshot(project) == before


def test_refuses_when_manifest_is_not_a_release(project: Path):
    assert _run("0.3.0rc2.dev0", str(project)).returncode == 0
    before = _snapshot(project)
    res = _run("0.3.0rc3.dev0", str(project))
    assert res.returncode != 0
    assert "not a release" in res.stderr
    assert _snapshot(project) == before


def test_refuses_inconsistent_release_config(project: Path):
    cfg = project / "nextflow.config"
    cfg.write_text(
        re.sub(
            r"(?m)^(\s*container_tag\s*=\s*)'[^']*'",
            r"\g<1>'dev'",
            cfg.read_text(),
            count=1,
        )
    )
    before = _snapshot(project)
    res = _run("0.3.0rc2.dev0", str(project))
    assert res.returncode != 0
    assert "container_tag" in res.stderr
    assert _snapshot(project) == before


def test_usage_without_arguments():
    res = _run()
    assert res.returncode == 2
    assert "usage" in res.stderr.lower()


def test_refuses_before_the_release_tag_exists(project: Path):
    """Run it only after the tag is pushed: otherwise main would leave the
    release config before the tag's CI build reads it."""
    subprocess.run([*GIT, "-C", str(project), "tag", "-d", "v0.3.0rc1"], check=True, capture_output=True)
    before = _snapshot(project)
    res = _run("0.3.0rc2.dev0", str(project))
    assert res.returncode == 1, res.stdout + res.stderr
    assert "v0.3.0rc1" in res.stderr and "tag" in res.stderr
    assert _snapshot(project) == before
    res = _run("0.3.0rc2.dev0", str(project), ALLOW_UNTAGGED="1")
    assert res.returncode == 0, res.stdout + res.stderr
    assert _manifest(project) == "0.3.0rc2.dev0"


def test_half_applied_state_names_the_recovery(project: Path):
    """If a run stopped after writing nextflow.config only, a rerun must say how
    to recover instead of a bare 'not a release'."""
    cfg = project / "nextflow.config"
    cfg.write_text(re.sub(r"(?m)^(\s*version\s*=\s*)'[^']*'", r"\g<1>'0.3.0rc2.dev0'", cfg.read_text(), count=1))
    before = _snapshot(project)
    res = _run("0.3.0rc2.dev0", str(project))
    assert res.returncode == 1, res.stdout + res.stderr
    assert "git checkout --" in res.stderr
    assert _snapshot(project) == before


def test_old_python_gets_a_clear_message(project: Path, tmp_path: Path):
    """System python3 can be 3.6 here; the heredoc uses 3.8+ syntax."""
    fake = tmp_path / "python36"
    fake.write_text("#!/bin/sh\nexit 1\n")
    fake.chmod(0o755)
    before = _snapshot(project)
    res = _run("0.3.0rc2.dev0", str(project), PYTHON=str(fake))
    assert res.returncode == 1
    assert "3.8" in res.stderr
    assert _snapshot(project) == before
