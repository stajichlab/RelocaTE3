"""``tests/nextflow/parity_chr3.sh``: modes, arguments and the commands it runs.

The parity gate itself needs a cluster job (pixi mode) or the published image
(docker/apptainer mode). Here the script runs with ``DRY_RUN=1``, which checks
the arguments and prints each command instead of running it, so the command
lines CI and the cluster run are pinned without running them.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "tests" / "nextflow" / "parity_chr3.sh"
IMAGE = "ghcr.io/stajichlab/relocate3"


def _run(*args: str, dry: bool = True) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if k != "DRY_RUN"}
    if dry:
        env["DRY_RUN"] = "1"
    return subprocess.run(
        ["bash", str(SCRIPT), *args], cwd=ROOT, env=env, capture_output=True, text=True
    )


def _lines(proc: subprocess.CompletedProcess, needle: str) -> list[str]:
    return [line for line in proc.stdout.splitlines() if needle in line]


def test_help_prints_usage():
    proc = _run("--help", dry=False)
    assert proc.returncode == 0
    assert "Usage" in proc.stdout
    for mode in ("pixi", "docker", "apptainer"):
        assert mode in proc.stdout


@pytest.mark.parametrize(
    "args",
    [("podman", "dev"), ("docker",), ("apptainer",), ("docker", "a b"),
     ("docker", "dev", "out", "extra"), ("pixi", "out", "extra")],
)
def test_bad_arguments_fail_with_usage(args, tmp_path):
    proc = _run(*args)
    assert proc.returncode != 0
    assert "Usage" in proc.stderr


def test_pixi_mode_is_the_default(tmp_path):
    out = tmp_path / "out"
    for args in ((), ("pixi",)):
        proc = _run(*args, str(out)) if args else _run()
        assert proc.returncode == 0, proc.stderr
        nf = _lines(proc, "nextflow run")
        assert len(nf) == 1 and "-profile test,pixi" in nf[0]
        assert "--container_tag" not in nf[0]
        run_all = _lines(proc, "relocaTE3 run-all")
        assert len(run_all) == 1
        assert "docker" not in run_all[0] and "apptainer" not in run_all[0]
    assert not out.exists(), "a dry run must not create the output directory"


def test_docker_mode_runs_both_sides_in_the_image(tmp_path):
    out = tmp_path / "out"
    proc = _run("docker", "sha-abc1234", str(out))
    assert proc.returncode == 0, proc.stderr
    image = f"{IMAGE}:sha-abc1234"
    run_all = _lines(proc, "relocaTE3 run-all")
    assert len(run_all) == 1
    assert run_all[0].startswith("docker run --rm -u ")
    assert image in run_all[0]
    nf = _lines(proc, "nextflow run")
    assert len(nf) == 1
    assert "-profile test,docker --container_tag sha-abc1234" in nf[0]
    compare = _lines(proc, "compare_outputs.sh")
    assert len(compare) == 1 and image in compare[0]
    assert not out.exists()


def test_apptainer_mode_runs_both_sides_in_the_image(tmp_path):
    proc = _run("apptainer", "pr-7", str(tmp_path / "out"))
    assert proc.returncode == 0, proc.stderr
    image = f"docker://{IMAGE}:pr-7"
    run_all = _lines(proc, "relocaTE3 run-all")
    assert len(run_all) == 1
    assert run_all[0].startswith("apptainer exec ")
    assert image in run_all[0]
    nf = _lines(proc, "nextflow run")
    assert "-profile test,apptainer --container_tag pr-7" in nf[0]
    compare = _lines(proc, "compare_outputs.sh")
    assert len(compare) == 1 and image in compare[0]


def test_existing_outdir_is_refused(tmp_path):
    proc = _run("pixi", str(tmp_path))
    assert proc.returncode != 0
    assert "refusing to overwrite" in proc.stderr
