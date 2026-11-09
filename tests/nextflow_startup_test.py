"""main.nf starts under every software profile, and ``--help`` exits 0.

``nextflow run -preview`` runs main.nf's workflow body (parameter validation,
the parameter summary log, the sample sheet) and builds the task graph, but
runs no task and pulls no image. The container profiles set
``process.container`` to a single image string; nf-schema 2.7.3's
``paramsSummaryLog`` cast it to a Map and crashed at startup (preview job
31638101).

Each run launches from a temporary directory with a sample sheet of absolute
read paths (nf-schema resolves relative sheet paths against the launch
directory). Needs the pixi ``nf`` environment; skipped without.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
READS = ROOT / "tests/data/sim_reads/MSU7.Chr3_2M.ALL_reads_6X_100_500"
# The default image tag: `dev` on a development manifest, the version on a
# release (tests/nextflow_schema_test.py checks which).
DEFAULT_TAG = re.search(
    r"(?m)^\s*container_tag\s*=\s*'([^']+)'", (ROOT / "nextflow.config").read_text()
).group(1)

pytestmark = pytest.mark.skipif(
    shutil.which("pixi") is None or not (ROOT / ".pixi" / "envs" / "nf").is_dir(),
    reason="needs pixi and the 'nf' environment (pixi install -e nf)",
)

NEXTFLOW = ["pixi", "run", "--manifest-path", str(ROOT / "pixi.toml"), "-e", "nf", "nextflow"]


def _run(tmp_path: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [*NEXTFLOW, "run", str(ROOT), *args],
        cwd=tmp_path, capture_output=True, text=True, timeout=300,
    )


@pytest.mark.parametrize("engine", ["apptainer", "docker", "singularity"])
def test_container_profile_passes_startup(tmp_path, engine):
    sheet = tmp_path / "samplesheet.csv"
    sheet.write_text(f"sample,fastq_1,fastq_2\nchr3,{READS}_1.fq.gz,{READS}_2.fq.gz\n")
    proc = _run(tmp_path, "-profile", f"test,{engine}", "-preview",
                "--input", str(sheet), "--outdir", str(tmp_path / "results"))
    log = proc.stdout + proc.stderr
    assert proc.returncode == 0, log
    # The summary log was printed, with the container image in it.
    assert f"ghcr.io/stajichlab/relocate3:{DEFAULT_TAG}" in log
    assert "Cannot cast" not in log


def test_help_prints_the_parameters_and_exits_0(tmp_path):
    # No --input/--genome/--outdir: help must not run parameter validation.
    proc = _run(tmp_path, "--help")
    log = proc.stdout + proc.stderr
    assert proc.returncode == 0, log
    assert "Missing required parameter" not in log
    for group in ("Input/output options", "Calling options", "Genotyping options",
                  "Container options"):
        assert group in log, group
    for text in ("--input", "--genome", "--te_library", "--te_aligner", "--container_tag"):
        assert text in log, text
    # Help exits before the workflow: no task ran.
    work = tmp_path / "work"
    assert not work.exists() or not any(work.iterdir())


def test_help_for_one_parameter(tmp_path):
    proc = _run(tmp_path, "--help", "te_aligner")
    log = proc.stdout + proc.stderr
    assert proc.returncode == 0, log
    assert "--te_aligner" in log
    # Only that parameter: no parameter groups.
    assert "Calling options" not in log


# ------------------------------------------- test profile, pixi_env, record


def test_test_profile_runs_from_any_directory(tmp_path):
    # `nextflow run stajichlab/RelocaTE3 -profile test,...` launches outside the
    # repository; the test sheet's reads must not depend on the launch dir.
    out = tmp_path / "results"
    proc = _run(tmp_path, "-profile", "test,pixi", "-stub", "--outdir", str(out))
    log = proc.stdout + proc.stderr
    assert proc.returncode == 0, log
    assert (out / "chr3" / "results").is_dir()
    # The run records its software environment (pipeline_info).
    record = (out / "pipeline_info" / "software_environment.txt").read_text()
    assert "profile: test,pixi" in record
    assert "container_engine: none" in record
    assert "container: none" in record
    assert f"pixi_env: {ROOT}/.pixi/envs/pipeline" in record


def test_pixi_profile_needs_relocate3_in_pixi_env(tmp_path):
    proc = _run(tmp_path, "-profile", "test,pixi", "-preview", "--outdir", str(tmp_path / "r"),
                "--pixi_env", str(tmp_path / "no_env"))
    log = proc.stdout + proc.stderr
    assert proc.returncode != 0
    assert "pixi_env" in log and str(tmp_path / "no_env" / "bin" / "relocaTE3") in log


def test_container_profiles_ignore_pixi_env(tmp_path):
    proc = _run(tmp_path, "-profile", "test,apptainer", "-preview", "--outdir", str(tmp_path / "r"),
                "--pixi_env", str(tmp_path / "no_env"))
    assert proc.returncode == 0, proc.stdout + proc.stderr
