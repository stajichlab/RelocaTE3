"""Shared fixtures.

``uutils_path``: a ``PATH`` whose ``sort`` and ``split`` are uutils coreutils
(the Rust rewrite), the system coreutils on Ubuntu 26.04 and so on GitHub's
``ubuntu-latest`` runner. Two of its differences from GNU broke the suite there:
``split`` writes one empty chunk for empty input (GNU writes none), and
``sort -t=`` is rejected ("separator must be exactly one character long")
up to 0.10.0. The fixture lets those regressions show up on any machine.

The uutils binary is ``$RELOCATE3_UUTILS`` (a path to the multi-call
``coreutils`` binary) if set, otherwise uutils-coreutils 0.10.0 fetched with
``pixi exec`` (cached after the first run). Tests using it are skipped when
neither works, e.g. offline without the pixi cache.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

#: Last uutils release with both behaviours (``sort -t=`` accepted from 0.11.0).
UUTILS_SPEC = "uutils-coreutils==0.10.0"


def _find_uutils() -> str | None:
    explicit = os.environ.get("RELOCATE3_UUTILS")
    if explicit:
        return explicit
    pixi = shutil.which("pixi")
    if pixi is None:
        return None
    try:
        proc = subprocess.run(
            [pixi, "exec", "--spec", UUTILS_SPEC, "--", "sh", "-c", "command -v coreutils"],
            capture_output=True, text=True, timeout=300,
        )  # fmt: skip
    except (OSError, subprocess.TimeoutExpired):
        return None
    lines = proc.stdout.strip().splitlines()
    return lines[-1] if proc.returncode == 0 and lines else None


@pytest.fixture(scope="session")
def uutils_path(tmp_path_factory) -> str:
    """``PATH`` with uutils ``sort`` and ``split`` in front; skip if unavailable."""
    binary = _find_uutils()
    if not binary or not Path(binary).is_file():
        pytest.skip(
            f"uutils coreutils unavailable (set RELOCATE3_UUTILS or allow pixi exec {UUTILS_SPEC})"
        )
    version = subprocess.run(
        [binary, "--version"], capture_output=True, text=True
    ).stdout
    if "multi-call" not in version:
        pytest.skip(
            f"{binary} is not the uutils multi-call coreutils: {version.strip()!r}"
        )
    shim = tmp_path_factory.mktemp("uutils_shim")
    for tool in ("sort", "split"):
        (shim / tool).symlink_to(binary)
    probe = subprocess.run(
        [str(shim / "split"), "--version"], capture_output=True, text=True
    )
    assert probe.returncode == 0 and "split" in probe.stdout, (
        probe.stdout + probe.stderr
    )
    return f"{shim}{os.pathsep}{os.environ.get('PATH', '')}"
