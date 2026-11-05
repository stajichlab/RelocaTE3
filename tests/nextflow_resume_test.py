"""A sample's genotype_bam must not invalidate its other cached tasks.

Nextflow's ``-resume`` hashes every task input, ``val(meta)`` included. If the
genotype BAM path lived in ``meta`` (it did in M1), adding, moving or touching
that BAM changed the hash of every task of the sample, so ``-resume`` reran the
whole sample from the TE search. Only CHARACTERIZE reads the BAM, so only
CHARACTERIZE may rerun.

The pipeline runs in stub mode (each process runs its ``stub:`` block): task
hashing and caching are the same as in a real run, and three runs take well
under a minute. Needs the pixi ``nf`` environment (Nextflow); skipped without.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
READS = ROOT / "tests" / "data" / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500"

pytestmark = pytest.mark.skipif(
    shutil.which("pixi") is None or not (ROOT / ".pixi" / "envs" / "nf").is_dir(),
    reason="needs pixi and the 'nf' environment (pixi install -e nf)",
)


def _run(launch: Path, sheet: Path, resume: bool) -> dict[str, str]:
    """Run the stub pipeline from ``launch``; return {simple process name: status}."""
    cmd = [
        "pixi", "run", "--manifest-path", str(ROOT / "pixi.toml"), "-e", "nf",
        "nextflow", "run", str(ROOT), "-profile", "test,pixi", "-stub",
        "--input", str(sheet), "--outdir", str(launch / "out"),
        "-work-dir", str(launch / "work"),
    ]
    if resume:
        cmd.append("-resume")
    proc = subprocess.run(cmd, cwd=launch, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    trace = launch / "out" / "pipeline_info" / "execution_trace.txt"
    with trace.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    status = {}
    for row in rows:
        simple = row["name"].split(" ")[0].split(":")[-1]
        assert simple not in status, f"{simple} ran twice (one sample expected)"
        status[simple] = row["status"]
    return status


def _sheet(path: Path, bam: Path | None) -> Path:
    path.write_text(
        "sample,fastq_1,fastq_2,genotype_bam\n"
        f"chr3,{READS}_1.fq.gz,{READS}_2.fq.gz,{bam or ''}\n"
    )
    return path


UPSTREAM = (
    "RELOCATE3_INDEX_GENOME",
    "RELOCATE3_TE_SEARCH",
    "RELOCATE3_ALIGN_GENOME",
    "RELOCATE3_FIND_INSERTIONS",
    "RELOCATE3_FIND_REFERENCE",
)


def test_genotype_bam_changes_rerun_only_characterize(tmp_path):
    bam = tmp_path / "given.bam"
    bam.write_text("")
    (tmp_path / "given.bam.bai").write_text("")

    # 1. No genotype_bam: the reads are aligned for genotyping.
    first = _run(tmp_path, _sheet(tmp_path / "plain.csv", None), resume=False)
    assert all(first[p] == "COMPLETED" for p in UPSTREAM), first
    assert first["RELOCATE3_ALIGN_READS"] == "COMPLETED"

    # 2. The same sample, now with a genotype_bam: nothing upstream reruns.
    with_bam = _sheet(tmp_path / "with_bam.csv", bam)
    second = _run(tmp_path, with_bam, resume=True)
    assert {p: second[p] for p in UPSTREAM} == dict.fromkeys(UPSTREAM, "CACHED"), second
    assert "RELOCATE3_ALIGN_READS" not in second
    assert second["RELOCATE3_CHARACTERIZE"] == "COMPLETED"

    # 3. The BAM is touched (new mtime): only CHARACTERIZE reruns.
    st = bam.stat()
    os.utime(bam, (st.st_atime + 60, st.st_mtime + 60))
    third = _run(tmp_path, with_bam, resume=True)
    assert {p: third[p] for p in UPSTREAM} == dict.fromkeys(UPSTREAM, "CACHED"), third
    assert third["RELOCATE3_CHARACTERIZE"] == "COMPLETED"


def test_resume_with_nothing_changed_caches_multiqc(tmp_path):
    """The report's inputs (summary files, versions file, comment) hold nothing
    per-run (no run name or start time), so an unchanged -resume caches it."""
    sheet = _sheet(tmp_path / "plain.csv", None)
    first = _run(tmp_path, sheet, resume=False)
    assert first["MULTIQC"] == "COMPLETED"
    versions = tmp_path / "out" / "pipeline_info" / "software_versions.yml"
    first_versions = versions.read_text()

    second = _run(tmp_path, sheet, resume=True)
    assert set(second.values()) == {"CACHED"}, second
    assert versions.read_text() == first_versions
