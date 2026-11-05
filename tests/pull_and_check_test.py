"""containers/pull_and_check.sh: pull the published image anonymously, check it.

The anonymous pull is what proves the ghcr package is public; a SIF already in
the cache skips it, so ``--force-pull`` (or ``PULL=1``) pulls anyway. A pull
goes to ``<sif>.tmp`` first: a stale one is removed before pulling, and a
failed pull leaves neither the temporary file nor a changed SIF.

``apptainer`` is a stub here: ``pull DEST IMAGE`` refuses an existing DEST (as
apptainer does without --force) and fails when STUB_PULL_FAIL=1.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "containers" / "pull_and_check.sh"
TAG = "sha-0123456"

STUB = """#!/usr/bin/env bash
echo "apptainer $*" >> "${STUB_LOG}"
case "$1" in
  --version) echo "apptainer version 1.5.3" ;;
  pull)
    [[ -e "$2" ]] && { echo "FATAL: Image file already exists: $2" >&2; exit 255; }
    [[ "${STUB_PULL_FAIL:-0}" == 1 ]] && { echo partial > "$2"; echo "FATAL: pull failed" >&2; exit 255; }
    echo "fresh image" > "$2" ;;
  inspect) echo "org.opencontainers.image.revision: 0123456789abcdef" ;;
  exec) exit 0 ;;
esac
"""


def _run(tmp_path: Path, *args: str, **env: str):
    bindir = tmp_path / "bin"
    bindir.mkdir(exist_ok=True)
    (bindir / "apptainer").write_text(STUB)
    (bindir / "apptainer").chmod(0o755)
    cache = tmp_path / "cache"
    cache.mkdir(exist_ok=True)
    full = {k: v for k, v in os.environ.items() if k not in ("PULL", "STUB_PULL_FAIL")}
    full.update(PATH=f"{bindir}:{os.environ['PATH']}", NXF_APPTAINER_CACHEDIR=str(cache),
                STUB_LOG=str(tmp_path / "calls"), **env)
    proc = subprocess.run(["bash", str(SCRIPT), *args], cwd=tmp_path, env=full,
                          capture_output=True, text=True)
    calls = (tmp_path / "calls").read_text() if (tmp_path / "calls").exists() else ""
    return proc, cache / f"ghcr.io-stajichlab-relocate3-{TAG}.img", calls


def _cached(tmp_path: Path) -> Path:
    sif = tmp_path / "cache" / f"ghcr.io-stajichlab-relocate3-{TAG}.img"
    sif.parent.mkdir(parents=True, exist_ok=True)
    sif.write_text("old image\n")
    return sif


def test_fresh_pull_then_check(tmp_path):
    proc, sif, calls = _run(tmp_path, TAG)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert sif.read_text() == "fresh image\n"
    assert "apptainer pull" in calls and "RESULT: PASS" in proc.stdout
    assert not Path(f"{sif}.tmp").exists()


def test_cached_sif_skips_the_pull_and_says_so(tmp_path):
    _cached(tmp_path)
    proc, sif, calls = _run(tmp_path, TAG)
    assert proc.returncode == 0, proc.stderr
    assert "apptainer pull" not in calls
    assert sif.read_text() == "old image\n"
    assert "--force-pull" in proc.stdout  # tells the user the public check was skipped


@pytest.mark.parametrize("how", [("--force-pull", TAG), ("PULL=1", TAG)])
def test_force_pull_ignores_the_cache(tmp_path, how):
    _cached(tmp_path)
    args, env = ((how[0], how[1]), {}) if how[0].startswith("--") else ((how[1],), {"PULL": "1"})
    proc, sif, calls = _run(tmp_path, *args, **env)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "apptainer pull" in calls
    assert sif.read_text() == "fresh image\n"


def test_stale_tmp_is_removed_before_pulling(tmp_path):
    sif = tmp_path / "cache" / f"ghcr.io-stajichlab-relocate3-{TAG}.img"
    sif.parent.mkdir(parents=True)
    Path(f"{sif}.tmp").write_text("stale\n")
    proc, sif, calls = _run(tmp_path, TAG)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert sif.read_text() == "fresh image\n"


def test_failed_pull_leaves_no_tmp_and_keeps_the_sif(tmp_path):
    _cached(tmp_path)
    proc, sif, calls = _run(tmp_path, "--force-pull", TAG, STUB_PULL_FAIL="1")
    assert proc.returncode != 0
    assert not Path(f"{sif}.tmp").exists()
    assert sif.read_text() == "old image\n"


def test_unknown_option_is_an_error(tmp_path):
    proc, _, _ = _run(tmp_path, "--bogus")
    assert proc.returncode != 0
    assert "Usage" in proc.stderr
