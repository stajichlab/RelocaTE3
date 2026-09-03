"""Mate planning must not depend on the input FASTQ naming convention."""

from __future__ import annotations

import shutil
from pathlib import Path

import pysam
import pytest

from RelocaTE3.aligners import canonical_name, get_aligner
from RelocaTE3.genome_align import plan_genome_alignment_inputs
from RelocaTE3.ReadLibrary import ReadLibrary


def _write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as out:
        for name, sequence in records:
            out.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def test_canonical_name_handles_every_convention():
    assert canonical_name("read_500_1/1", "1") == "read_500_1/1"
    assert canonical_name("read_500_1/2", "2") == "read_500_1/2"
    illumina = "A00519:123:HXXXX:1:1101:6768:1000"
    assert canonical_name(illumina, "1") == f"{illumina}/1"
    assert canonical_name(illumina, "2") == f"{illumina}/2"
    assert canonical_name("SRR1234567.1", "1") == "SRR1234567.1/1"
    assert canonical_name(f"{illumina} 1:N:0:ATCACG", "1") == f"{illumina}/1"
    assert canonical_name("read_500_1/1", "2") == "read_500_1/2"


@pytest.mark.parametrize("with_suffix", [False, True])
def test_pair_planner_finds_suffixless_and_suffixed_original_mates(
    tmp_path: Path, with_suffix: bool
):
    """The source file supplies mate identity even when its names are bare."""
    name1 = "readA/1" if with_suffix else "readA"
    name2 = "readA/2" if with_suffix else "readA"
    r1, r2 = tmp_path / "R1.fq", tmp_path / "R2.fq"
    _write_fastq(r1, [(name1, "A" * 40)])
    _write_fastq(r2, [(name2, "G" * 40)])
    flank = tmp_path / "flank.fq"
    _write_fastq(flank, [("readA/1:start:5", "C" * 20)])
    reads = ReadLibrary([str(r1), str(r2)], "sample")

    plan = plan_genome_alignment_inputs(
        [str(flank)],
        {"readA/1"},
        reads,
        tmp_path / "p1.fq",
        tmp_path / "p2.fq",
        tmp_path / "j.fq",
        tmp_path / "s.fq",
    )

    assert (plan.paired, plan.single_junctions, plan.support_reads) == (1, 0, 0)
    assert plan.retag["readA"] == ("readA/1:start:5", "readA/2")
    assert "G" * 40 in (tmp_path / "p2.fq").read_text()


TE_SEQ = (
    "ATGGCCTTAAGGCCATTGACCTTGGAACCTTGGCATTGCAAGGTTCCAAGGTTAACCGGTTAA"
    "CCGGTTAACCGGTTAAGGCCTTAACCGGTTAAGGCCTTAACCGGTTAAGGCCTTAACCGGTTA"
)


@pytest.mark.parametrize("backend_name", ["minimap2", "bwa", "bwamem2", "bowtie2"])
def test_te_library_bams_carry_mate_suffix(tmp_path: Path, backend_name: str):
    binary = {
        "minimap2": "minimap2",
        "bwa": "bwa",
        "bwamem2": "bwa-mem2",
        "bowtie2": "bowtie2",
    }[backend_name]
    if shutil.which(binary) is None:
        pytest.skip(f"{binary} not available")

    te_fa = tmp_path / "te.fa"
    te_fa.write_text(f">mPing\n{TE_SEQ}\n")
    r1, r2 = tmp_path / "sample_R1.fastq", tmp_path / "sample_R2.fastq"
    name = "A00519:1:HXXXX:1:1101:6768:1000"
    _write_fastq(r1, [(name, TE_SEQ[:80])])
    _write_fastq(r2, [(name, TE_SEQ[40:120])])

    reads = ReadLibrary([str(r1), str(r2)], "sample")
    outdir = tmp_path / "out"
    outdir.mkdir()
    bams = get_aligner(backend_name, threads=1).map_te_library(
        reads, str(te_fa), str(outdir)
    )
    seen = []
    for bam in bams:
        with pysam.AlignmentFile(str(bam), "rb") as alignments:
            seen.extend(rec.query_name for rec in alignments.fetch(until_eof=True))
    assert seen
    assert all(name.endswith(("/1", "/2")) for name in seen)


def test_single_end_names_are_left_alone(tmp_path: Path):
    if shutil.which("minimap2") is None:
        pytest.skip("minimap2 not available")
    te_fa = tmp_path / "te.fa"
    te_fa.write_text(f">mPing\n{TE_SEQ}\n")
    r1 = tmp_path / "sample.fastq"
    _write_fastq(r1, [("SRR1234567.1", TE_SEQ[:80])])
    reads = ReadLibrary([str(r1)], "sample")
    outdir = tmp_path / "out"
    outdir.mkdir()
    bams = get_aligner("minimap2", threads=1).map_te_library(
        reads, str(te_fa), str(outdir)
    )
    names = []
    for bam in bams:
        with pysam.AlignmentFile(str(bam), "rb") as alignments:
            names.extend(rec.query_name for rec in alignments.fetch(until_eof=True))
    assert names and all(name == "SRR1234567.1" for name in names)
