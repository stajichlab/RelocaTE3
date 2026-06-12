"""Tests for single-sample pipeline orchestration and the trim CLI."""

from __future__ import annotations

import os
from pathlib import Path

import pysam

from RelocaTE3.cli import main
from RelocaTE3.pipeline import run_sample
from RelocaTE3.ReadLibrary import ReadLibrary

DATA = Path(__file__).parent / "data"
R1 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"
R2 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"
TELIB = DATA / "mping.fa"
GENOME = DATA / "sim_genome" / "MSU7.Chr3_2M.fa"


def test_run_sample_produces_outputs(tmp_path: Path):
    """run_sample trims, re-aligns, calls insertions, and returns the GFF path."""
    reads = ReadLibrary([str(R1), str(R2)], "HEG4")
    gff = run_sample(reads, str(TELIB), str(GENOME), tmp_path, threads=2)
    assert (tmp_path / "te_containing" / "HEG4.read_repeat_name.txt").exists()
    left = tmp_path / "flanking" / "HEG4.left.flankingReads.fq"
    assert left.exists() and os.path.getsize(left) > 0

    # genome BAM produced and indexed, containing tagged junction reads
    bam = tmp_path / "genome_aln" / "HEG4.genome.bam"
    assert bam.exists() and Path(str(bam) + ".bai").exists()
    with pysam.AlignmentFile(str(bam), "rb") as bf:
        names = [r.query_name for r in bf.fetch() if not r.is_unmapped]
    assert any(":5" in n or ":3" in n for n in names)

    # insertion GFF returned and non-empty
    assert gff.exists() and gff == tmp_path / "results" / "HEG4.all_nonref_insert.gff"
    assert os.path.getsize(gff) > 0


def test_trim_cli(tmp_path: Path):
    """The `trim` subcommand runs end-to-end and exits 0."""
    rc = main(
        [
            "trim",
            "--r1",
            str(R1),
            "--r2",
            str(R2),
            "-t",
            str(TELIB),
            "-o",
            str(tmp_path),
            "--sample",
            "HEG4",
            "-c",
            "2",
        ]
    )
    assert rc == 0
    assert (tmp_path / "flanking" / "HEG4.right.flankingReads.fq").exists()
