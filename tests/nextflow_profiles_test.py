"""Execution profiles resolve to the intended Nextflow settings.

``nextflow config -flat`` shows what a profile sets with its default params.
It takes no ``--param`` options, so CLI overrides are checked by running a
one-line probe workflow with the pipeline's ``nextflow.config``: it prints the
resolved settings from the session config, which is what tasks get.

Needs the pixi ``nf`` environment (Nextflow); skipped without. No container is
pulled and no job is submitted (the probe has no processes).
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
IMAGE = "ghcr.io/stajichlab/relocate3"

pytestmark = pytest.mark.skipif(
    shutil.which("pixi") is None or not (ROOT / ".pixi" / "envs" / "nf").is_dir(),
    reason="needs pixi and the 'nf' environment (pixi install -e nf)",
)

NEXTFLOW = ["pixi", "run", "--manifest-path", str(ROOT / "pixi.toml"), "-e", "nf", "nextflow"]


def _config(profile: str) -> dict[str, str]:
    """``nextflow config -flat -profile <profile>`` as {key: raw value}."""
    proc = subprocess.run(
        [*NEXTFLOW, "config", "-flat", "-profile", profile, str(ROOT)],
        cwd=ROOT, capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    out = {}
    for line in proc.stdout.splitlines():
        key, sep, value = line.partition(" = ")
        if sep:
            out[key.strip()] = value.strip()
    return out


PROBE = """\
workflow {
    def c = workflow.session.config
    println "PROBE container=${c.process.container}"
    println "PROBE queue=${c.process.queue}"
    println "PROBE executor=${c.process.executor}"
}
"""


def _probe(tmp_path: Path, profile: str, *args: str) -> dict[str, str]:
    """Settings a run of the pipeline config resolves, with CLI params ``args``."""
    (tmp_path / "probe.nf").write_text(PROBE)
    proc = subprocess.run(
        [*NEXTFLOW, "-q", "run", "probe.nf", "-c", str(ROOT / "nextflow.config"),
         "-profile", profile, *args],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    out = {}
    for line in proc.stdout.splitlines():
        if line.startswith("PROBE "):
            key, _, value = line[len("PROBE "):].partition("=")
            out[key] = value
    return out


# ----------------------------------------------------------- containers


@pytest.mark.parametrize("engine", ["apptainer", "singularity"])
def test_apptainer_like_profiles(engine):
    cfg = _config(f"test,{engine}")
    assert cfg["process.container"] == f"'{IMAGE}:dev'"
    assert cfg[f"{engine}.enabled"] == "true"
    assert cfg[f"{engine}.autoMounts"] == "true"
    assert cfg[f"{engine}.pullTimeout"] == "'60 min'"
    # Nextflow reads NXF_<ENGINE>_CACHEDIR itself; a cacheDir here would win.
    assert f"{engine}.cacheDir" not in cfg
    for other in {"docker", "apptainer", "singularity"} - {engine}:
        assert cfg[f"{other}.enabled"] == "false"


def test_docker_profile():
    cfg = _config("test,docker")
    assert cfg["process.container"] == f"'{IMAGE}:dev'"
    assert cfg["docker.enabled"] == "true"
    assert cfg["docker.runOptions"] == "'-u $(id -u):$(id -g)'"
    assert cfg["apptainer.enabled"] == "false"
    assert cfg["singularity.enabled"] == "false"


def test_pixi_profile_has_no_container():
    cfg = _config("test,pixi")
    assert "process.container" not in cfg
    for engine in ("docker", "apptainer", "singularity"):
        assert cfg.get(f"{engine}.enabled", "false") == "false"


def test_container_tag_cli_override(tmp_path):
    assert _probe(tmp_path, "apptainer")["container"] == f"{IMAGE}:dev"
    got = _probe(tmp_path, "apptainer", "--container_tag", "pr-7")
    assert got["container"] == f"{IMAGE}:pr-7"
    got = _probe(tmp_path, "docker", "--container_registry", "example.org/rt3",
                 "--container_tag", "sha-abc1234")
    assert got["container"] == "example.org/rt3:sha-abc1234"


# ----------------------------------------------------------------- hpcc


def test_hpcc_profile_submits_to_slurm():
    cfg = _config("hpcc,apptainer")
    assert cfg["process.executor"] == "'slurm'"
    assert cfg["process.queue"] == "'intel'"
    # intel nodes: 64 CPUs, 450000 MB (439 GiB) RAM, 30-day limit.
    assert cfg["process.resourceLimits.cpus"] == "64"
    assert cfg["process.resourceLimits.memory"] == "'430 GB'"
    assert cfg["process.resourceLimits.time"] == "'30d'"
    assert cfg["executor.$slurm.queueSize"] == "50"
    assert cfg["executor.$slurm.submitRateLimit"] == "'10/1min'"
    # The apptainer profile still provides the image; no cacheDir override.
    assert cfg["process.container"] == f"'{IMAGE}:dev'"
    assert "apptainer.cacheDir" not in cfg


def test_slurm_partition_cli_override(tmp_path):
    assert _probe(tmp_path, "hpcc,apptainer") == {
        "container": f"{IMAGE}:dev", "queue": "intel", "executor": "slurm",
    }
    got = _probe(tmp_path, "hpcc,apptainer", "--slurm_partition", "epyc")
    assert got["queue"] == "epyc"


def test_test_profile_stays_local():
    cfg = _config("test,pixi")
    assert "process.executor" not in cfg
    assert "process.queue" not in cfg


def test_nf_test_tasks_fit_a_github_runner():
    """Every nf-test task must fit ubuntu-latest (4 CPUs, 16 GB).

    Module tests run without the `test` profile, so tests/nextflow.config (applied
    by nf-test to every test) caps resources. Without it RELOCATE3_TE_SEARCH
    (process_high, 8 CPUs) failed on CI: "Task requirement exceeds available
    CPUs -- req: 8; avail: 4".
    """
    import re

    text = (ROOT / "tests" / "nextflow.config").read_text()
    m = re.search(r"resourceLimits\s*=\s*\[([^\]]*)\]", text)
    assert m, "tests/nextflow.config must set process.resourceLimits"
    limits = m.group(1)
    cpus = int(re.search(r"cpus:\s*(\d+)", limits).group(1))
    mem_gb = int(re.search(r"memory:\s*(\d+)\.GB", limits).group(1))
    assert cpus <= 4 and mem_gb <= 14, limits


def test_blat_scatter_resources_fit_the_observed_30x_run():
    """Scatter labels sized from the riceTElib 30x trace (2026-11-05).

    TE_SEARCH_CHUNK (8 chunks, 4 BLATs at once) peaked at 0.36 GB, so 2 GB;
    TE_MERGE_TRIM peaked at 3.2 GB using ~1.2 CPUs for 32 min, so 2 CPUs and
    8 GB (8 h), not process_medium's 4 CPUs. Memory and time still grow with
    the attempt, and 137/140 (OOM / SLURM limit) still retry.
    """
    cfg = _config("pixi")
    chunk = "process.'withLabel:process_blat_chunk'"
    assert cfg[f"{chunk}.cpus"] == "4"
    assert cfg[f"{chunk}.memory"] == "{ 2.GB * task.attempt }"
    assert cfg[f"{chunk}.time"] == "{ 4.h * task.attempt }"
    merge = "process.'withLabel:process_merge_trim'"
    assert cfg[f"{merge}.cpus"] == "2"
    assert cfg[f"{merge}.memory"] == "{ 8.GB * task.attempt }"
    assert cfg[f"{merge}.time"] == "{ 8.h * task.attempt }"
    assert "[137, 140]" in cfg["process.errorStrategy"]
    assert cfg["process.maxRetries"] == "2"
    module = (ROOT / "modules" / "local" / "relocate3_te_merge_trim" / "main.nf").read_text()
    assert "label 'process_merge_trim'" in module


# ------------------------------------------- pipeline_info/software_environment.txt

ENV_PROBE = """\
include { softwareEnvironment } from '{root}/subworkflows/local/software_environment.nf'
workflow {
    print softwareEnvironment()
}
"""


def _software_environment(tmp_path: Path, profile: str, *args: str) -> str:
    (tmp_path / "env_probe.nf").write_text(ENV_PROBE.replace("{root}", str(ROOT)))
    proc = subprocess.run(
        [*NEXTFLOW, "-q", "run", "env_probe.nf", "-c", str(ROOT / "nextflow.config"),
         "-profile", profile, *args],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return proc.stdout


@pytest.mark.parametrize("engine", ["apptainer", "docker", "singularity"])
def test_software_environment_records_the_container(tmp_path, engine):
    record = _software_environment(tmp_path, engine, "--container_tag", "sha-abc1234")
    assert f"container_engine: {engine}" in record
    assert f"container: {IMAGE}:sha-abc1234" in record
    assert "digest: not recorded by Nextflow" in record
    assert "pixi_env:" not in record


def test_software_environment_records_the_pixi_env(tmp_path):
    record = _software_environment(tmp_path, "pixi", "--pixi_env", "/some/env")
    assert "container_engine: none" in record
    assert "container: none" in record
    assert "pixi_env: /some/env" in record
