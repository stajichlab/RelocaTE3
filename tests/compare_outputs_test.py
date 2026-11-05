"""validation/blat_memory/compare_outputs.sh: the parity gates' file comparison.

Text files byte for byte; BAMs by their records and their @SQ header lines
(@PG lines legitimately differ). A BAM samtools cannot read is a failure, never
"identical": before, two unreadable BAMs both gave empty `samtools view`
output and compared equal. Uses the blat env's real samtools on tiny BAMs.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "validation" / "blat_memory" / "compare_outputs.sh"
PIXI_BIN = ROOT / ".pixi" / "envs" / "blat" / "bin"

pytestmark = pytest.mark.skipif(
    not (PIXI_BIN / "samtools").exists(), reason="needs the blat env's samtools (pixi install -e blat)"
)

RECORD = "r1\t0\tchr1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\n"


def _bam(path: Path, sq_len: int = 100, pg: str = "a", records: str = RECORD) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    sam = f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{sq_len}\n@PG\tID:{pg}\tPN:{pg}\n{records}"
    subprocess.run([str(PIXI_BIN / "samtools"), "view", "-b", "-o", str(path), "-"],
                   input=sam, text=True, check=True)
    return path


def _pair(tmp_path: Path) -> tuple[Path, Path]:
    a, b = tmp_path / "A", tmp_path / "B"
    for d in (a, b):
        (d / "results").mkdir(parents=True)
        (d / "results" / "calls.txt").write_text("same\n")
    return a, b


def _compare(a: Path, b: Path):
    env = {k: v for k, v in os.environ.items() if k != "SLURM_JOB_ID"}
    env["PIXI_BIN"] = str(PIXI_BIN)
    return subprocess.run(["bash", str(SCRIPT), str(a), str(b)], env=env, capture_output=True, text=True)


def test_identical_bams_differing_only_in_pg_pass(tmp_path):
    a, b = _pair(tmp_path)
    _bam(a / "s.bam", pg="run1")
    _bam(b / "s.bam", pg="run2")
    proc = _compare(a, b)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "identical=2 differing=0" in proc.stdout


def test_different_records_fail(tmp_path):
    a, b = _pair(tmp_path)
    _bam(a / "s.bam")
    _bam(b / "s.bam", records=RECORD.replace("\t5\t", "\t6\t"))
    proc = _compare(a, b)
    assert proc.returncode != 0
    assert "DIFF (records) ./s.bam" in proc.stdout


def test_different_sq_header_fails(tmp_path):
    a, b = _pair(tmp_path)
    _bam(a / "s.bam", sq_len=100)
    _bam(b / "s.bam", sq_len=200)
    proc = _compare(a, b)
    assert proc.returncode != 0
    assert "DIFF (@SQ header) ./s.bam" in proc.stdout


@pytest.mark.parametrize("content", [b"", b"not a bam at all\n"], ids=["empty", "corrupt"])
def test_unreadable_bams_on_both_sides_fail(tmp_path, content):
    a, b = _pair(tmp_path)
    for d in (a, b):
        (d / "s.bam").write_bytes(content)
    proc = _compare(a, b)
    assert proc.returncode != 0
    assert "UNREADABLE" in proc.stdout and "s.bam" in proc.stdout


def test_unreadable_bam_on_one_side_fails(tmp_path):
    a, b = _pair(tmp_path)
    _bam(a / "s.bam")
    (b / "s.bam").write_bytes(b"")
    proc = _compare(a, b)
    assert proc.returncode != 0
    assert "UNREADABLE" in proc.stdout
