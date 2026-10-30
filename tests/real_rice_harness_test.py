"""The real-rice harness must run RelocaTE3 the way RelocaTE2 was run.

The 2026-10-26 release check found the harness, not RelocaTE3, behind a 1270-call
gap against the frozen RelocaTE2 calls (docs/2026-10-21-release-v0.2.0-validation.md):

* it compared RelocaTE2's filtered ``ALL.all_nonref_insert.gff`` with
  RelocaTE3's unfiltered ``.all_nonref_insert.txt``;
* it ran RelocaTE3 step by step with ``--tsd "..."``, ``--min-mapq 1`` and
  insert size 500, skipping run-all's mate-anchored alignment and full-read
  false-junction filter. RelocaTE2's real-rice runs used TSD ``UNK``
  (regex.txt) and library size 250.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
HARNESS = REPO / "validation" / "real_rice"

sys.path.insert(0, str(HARNESS))
import normalize_relocate3_nonref  # noqa: E402

GFF_ATTRS = (
    "ID=x;Name=mPing;TSD={tsd};Note=Non-reference, not found in reference;"
    "Right_junction_reads={r};Left_junction_reads={l};"
    "Right_support_reads=0;Left_support_reads=0;"
)


def _txt_row(chrom, start, tsd, r, l):
    return (
        f"mPing\t{tsd}\tS1\t{chrom}\t{start}..{start + 2}\t+\t"
        f"T:{r + l}\tR:{r}\tL:{l}\tST:0\tSR:0\tSL:0\n"
    )


def _gff_row(chrom, start, tsd, r, l):
    attrs = GFF_ATTRS.format(tsd=tsd, r=r, l=l)
    return f"{chrom}\tRelocaTE3\tS1\t{start}\t{start + 2}\t.\t+\t.\t{attrs}\n"


def test_nonref_normalizer_reads_the_same_tier_as_relocate2(tmp_path):
    """RelocaTE2's side is its filtered .gff; RelocaTE3's must be too."""
    results = tmp_path / "r3" / "S1" / "results"
    results.mkdir(parents=True)
    # Unfiltered table: one two-sided call plus a singleton RelocaTE2's
    # clean_false_positive.py would drop from the headline .gff.
    (results / "ALL.mping.all_nonref_insert.txt").write_text(
        _txt_row("Chr1", 1000, "TAA", 3, 4) + _txt_row("Chr1", 5000, "singleton", 1, 0)
    )
    (results / "ALL.mping.all_nonref_insert.gff").write_text(
        _gff_row("Chr1", 1000, "TAA", 3, 4)
    )
    csv_path = tmp_path / "samples.csv"
    csv_path.write_text("Sample\nS1\n")
    config = tmp_path / "config.toml"
    config.write_text(
        "[paths]\n"
        f'sample_csv = "{csv_path}"\n'
        f'relocate3_outdir = "{tmp_path / "r3"}"\n'
        f'report_dir = "{tmp_path / "report"}"\n'
        "[relocate3]\n"
        'te_name = "mping"\n'
        'target = "ALL"\n'
    )
    out = tmp_path / "calls.tsv"

    normalize_relocate3_nonref.main(["--config", str(config), "--out", str(out)])

    rows = out.read_text().splitlines()[1:]
    assert len(rows) == 1, rows
    fields = rows[0].split("\t")
    assert fields[1:4] == ["Chr1", "1000", "1002"]
    assert fields[6] == "TAA"
    assert fields[-1].endswith("ALL.mping.all_nonref_insert.gff")


@pytest.fixture
def harness_run(tmp_path):
    """Run run_relocate3.sh for one sample with a fake relocaTE3 that logs argv."""
    root = tmp_path / "repo"
    rr = root / "validation" / "real_rice"
    shutil.copytree(HARNESS, rr, ignore=shutil.ignore_patterns("__pycache__", "config.toml"))
    te_lib = root / "tests" / "data" / "TE_lib" / "mping.fa"
    te_lib.parent.mkdir(parents=True)
    te_lib.write_text(">mPing\nACGT\n")

    vroot = tmp_path / "validation_data"
    for rel in (
        "ref_genome/MSU_r7.fa",
        "repeatmasker/MSU_r7.fa.out",
        "fq_input/S1_trim_R1.fastq.gz",
        "fq_input/S1_trim_R2.fastq.gz",
        "aln_input/S1/S1.sorted.markdup.cram",
    ):
        p = vroot / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("x\n")
    (vroot / "sample_file").mkdir()
    (vroot / "sample_file" / "samples.csv").write_text("Sample\nS1\n")

    shutil.copy(rr / "config.example.toml", rr / "config.toml")

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    argv_log = tmp_path / "relocaTE3.argv"
    fake = fake_bin / "relocaTE3"
    fake.write_text(f'#!/usr/bin/bash\nprintf "%s\\n" "$@" >> "{argv_log}"\n')
    fake.chmod(0o755)

    env = {
        "PATH": f"{fake_bin}:{Path(sys.executable).parent}:/usr/bin:/bin",
        "HOME": str(tmp_path),
        "RELOCATE3_VALIDATION_ROOT": str(vroot),
    }
    proc = subprocess.run(
        ["bash", str(rr / "run_relocate3.sh"), str(rr / "config.toml"), "S1"],
        env=env,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    calls = argv_log.read_text().split("\n") if argv_log.exists() else []
    return calls, vroot


def _opt(argv, flag):
    return argv[argv.index(flag) + 1]


def test_harness_runs_the_shipped_run_all_pipeline(harness_run):
    argv, vroot = harness_run
    assert argv[0] == "run-all", argv
    # Mate-anchored alignment and the full-read filter need both mates.
    assert _opt(argv, "-1").endswith("S1_trim_R1.fastq.gz")
    assert _opt(argv, "-2").endswith("S1_trim_R2.fastq.gz")


def test_harness_uses_relocate2_real_rice_settings(harness_run):
    argv, _ = harness_run
    assert _opt(argv, "--tsd") == "UNK"
    assert _opt(argv, "--size") == "250"
    assert _opt(argv, "--min-mapq") == "0"
    assert _opt(argv, "--mismatch") == "2"


def test_harness_genotypes_against_the_external_alignment(harness_run):
    argv, vroot = harness_run
    assert "--genotype" in argv
    assert _opt(argv, "--genotype-bam") == str(
        vroot / "aln_input" / "S1" / "S1.sorted.markdup.cram"
    )
