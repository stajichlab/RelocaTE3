"""Tests for the RelocaTE2-compatible step-4 alignment-input planner."""

from __future__ import annotations

import errno
import random
import shutil
from pathlib import Path

import pysam
import pytest

from RelocaTE3.aligners import get_aligner
from RelocaTE3.genome_align import (
    align_flanks_anchored,
    align_to_genome,
    collect_junction_fullreads,
    plan_genome_alignment_inputs,
    read_te_hit_names,
    split_mate,
    strip_tag,
)
from RelocaTE3.ReadLibrary import ReadLibrary


def _write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as out:
        for name, sequence in records:
            out.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def _rc(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def test_read_name_helpers_and_all_hit_loader(tmp_path: Path):
    assert strip_tag("read/1:end:5") == "read/1"
    assert strip_tag("read/2:middle") == "read/2"
    assert strip_tag("M00:1:2:3/1") == "M00:1:2:3/1"
    assert split_mate("read/1") == ("read", "1")
    assert split_mate("nomate") == ("nomate", "")

    names = tmp_path / "sample.te_hit_names.txt"
    names.write_text("read/1\nread/2:middle\n\n")
    assert read_te_hit_names(names) == {"read/1", "read/2"}


def test_pair_planner_ports_complete_relocate2_state_machine(tmp_path: Path):
    """Exercise J/J, J/M, J/U, J/N, M/M, M/U, and M/N together."""
    bases = ("jj", "jm", "ju", "jn", "mm", "mu", "mn")
    r1 = tmp_path / "reads_R1.fq"
    r2 = tmp_path / "reads_R2.fq"
    _write_fastq(
        r1,
        [(f"{base}/1", "ACAC" if base == "jn" else "AAAA") for base in bases],
    )
    _write_fastq(
        r2,
        [(f"{base}/2", "GTGT" if base == "mn" else "TTTT") for base in bases],
    )
    reads = ReadLibrary([str(r1), str(r2)], "sample")

    left = tmp_path / "sample.left.flankingReads.fq"
    right = tmp_path / "sample.right.flankingReads.fq"
    _write_fastq(
        left,
        [
            ("jj/1:end:5", "JJJ1"),
            ("jm/1:end:5", "JMJ1"),
            ("ju/1:end:5", "JUJ1"),
            ("mm/1:middle", "MMM1"),
            ("mu/1:middle", "MUM1"),
            ("mn/1:middle", "MNM1"),
        ],
    )
    _write_fastq(
        right,
        [
            ("jj/2:start:3", "JJJ2"),
            ("jm/2:middle", "JMM2"),
            ("jn/2:start:3", "JNJ2"),
            ("mm/2:middle", "MMM2"),
        ],
    )
    all_te_hits = {
        "jj/1",
        "jj/2",
        "jm/1",
        "jm/2",
        "ju/1",
        "ju/2",
        "jn/2",
        "mm/1",
        "mm/2",
        "mu/1",
        "mu/2",
        "mn/1",
    }
    pair1 = tmp_path / "paired_R1.fq"
    pair2 = tmp_path / "paired_R2.fq"
    junctions = tmp_path / "junctions.fq"
    support = tmp_path / "support.fq"

    plan = plan_genome_alignment_inputs(
        [str(left), str(right)],
        all_te_hits,
        reads,
        pair1,
        pair2,
        junctions,
        support,
    )

    assert (plan.paired, plan.single_junctions, plan.support_reads) == (2, 2, 1)
    assert plan.retag == {
        "jj": ("jj/1:end:5", "jj/2:start:3"),
        "jn": ("jn/1", "jn/2:start:3"),
    }
    assert plan.junction_bases == frozenset({"jj", "jm", "ju", "jn"})

    pair1_text, pair2_text = pair1.read_text(), pair2.read_text()
    assert "JJJ1" in pair1_text and "JJJ2" in pair2_text
    # J/N with the junction in R2 must retain original mate order.
    assert "ACAC" in pair1_text and "JNJ2" in pair2_text
    assert "JMJ1" in junctions.read_text() and "JUJ1" in junctions.read_text()
    assert "GTGT" in support.read_text()
    combined = pair1_text + pair2_text + junctions.read_text() + support.read_text()
    assert "MMM1" not in combined and "MUM1" not in combined


def test_pair_planner_rejects_inconsistent_step3_artifacts(tmp_path: Path):
    reads_fq = tmp_path / "reads.fq"
    flank_fq = tmp_path / "flanks.fq"
    _write_fastq(reads_fq, [("read", "AAAA")])
    _write_fastq(flank_fq, [("read:end:5", "AA")])
    reads = ReadLibrary([str(reads_fq)], "sample")

    with pytest.raises(ValueError, match="absent from the all-TE-hit artifact"):
        plan_genome_alignment_inputs(
            [str(flank_fq)],
            set(),
            reads,
            tmp_path / "p1.fq",
            tmp_path / "p2.fq",
            tmp_path / "j.fq",
            tmp_path / "s.fq",
        )


def test_collect_junction_fullreads_uses_tagged_flanks_and_original_records(
    tmp_path: Path,
):
    r1 = tmp_path / "reads_R1.fq"
    r2 = tmp_path / "reads_R2.fq"
    r1.write_text("@keep\nAACCGG\n+\nABCDEF\n@middle\nTTTTTT\n+\nGHIJKL\n")
    r2.write_text("@other\nCCGGTT\n+\nMNOPQR\n")
    reads = ReadLibrary([str(r1), str(r2)], "sample")

    left = tmp_path / "sample.left.flankingReads.fq"
    right = tmp_path / "sample.right.flankingReads.fq"
    left.write_text("@keep/1:end:5\nAAC\n+\nABC\n@middle/1:middle\nTTTTTT\n+\nGHIJKL\n")
    right.write_text("@other/2:start:3\nGTT\n+\nPQR\n")
    out = tmp_path / "junction.fullreads.fq"

    count = collect_junction_fullreads([str(left), str(right)], reads, out)

    with pysam.FastxFile(str(out)) as records:
        observed = {rec.name: (rec.sequence, rec.quality) for rec in records}
    assert count == 2
    assert observed == {
        "keep/1": ("AACCGG", "ABCDEF"),
        "other/2": ("CCGGTT", "MNOPQR"),
    }


def test_collect_junction_fullreads_preserves_pairs_for_mate_anchoring(
    tmp_path: Path,
):
    r1 = tmp_path / "reads_R1.fq"
    r2 = tmp_path / "reads_R2.fq"
    _write_fastq(r1, [("pair/1", "AAAACCCC")])
    _write_fastq(r2, [("pair/2", "GGGGTTTT")])
    reads = ReadLibrary([str(r1), str(r2)], "sample")
    flank = tmp_path / "sample.left.flankingReads.fq"
    _write_fastq(flank, [("pair/2:start:3", "TTTT")])
    out_r1 = tmp_path / "junction.fullreads.R1.fq"
    out_r2 = tmp_path / "junction.fullreads.R2.fq"

    count = collect_junction_fullreads(
        [str(flank)], reads, out_r1, mate_fastq=out_r2
    )

    assert count == 1
    assert "@pair/1\nAAAACCCC" in out_r1.read_text()
    assert "@pair/2\nGGGGTTTT" in out_r2.read_text()


@pytest.mark.skipif(shutil.which("bwa") is None, reason="bwa not available")
@pytest.mark.parametrize("genome_aligner", ["bwa", "bwaaln"])
def test_align_to_genome_anchors_ambiguous_flank_to_mate(
    tmp_path: Path, genome_aligner: str
):
    rng = random.Random(7)
    genome_bases = [rng.choice("ACGT") for _ in range(2200)]
    block = "".join(rng.choice("ACGT") for _ in range(100))
    genome_bases[500:600] = list(block)
    genome_bases[1500:1600] = list(block)
    anchor = "".join(genome_bases[1650:1750])
    genome = tmp_path / "genome.fa"
    genome.write_text(f">chr1\n{''.join(genome_bases)}\n")

    out = tmp_path / "out"
    (out / "flanking").mkdir(parents=True)
    _write_fastq(
        out / "flanking" / "S.left.flankingReads.fq",
        [("r1/1:end:5", block)],
    )
    (out / "te_containing").mkdir()
    (out / "te_containing" / "S.read_repeat_name.txt").write_text("r1/1\tmPing\t+\n")
    (out / "te_containing" / "S.te_hit_names.txt").write_text("r1/1\n")

    r1 = tmp_path / "reads_1.fq"
    r2 = tmp_path / "reads_2.fq"
    _write_fastq(r1, [("r1/1", block)])
    _write_fastq(r2, [("r1/2", _rc(anchor))])
    reads = ReadLibrary([str(r1), str(r2)], "S")

    bam, fullreads_bam = align_to_genome(
        reads, str(genome), out, threads=1, genome_aligner=genome_aligner
    )
    with pysam.AlignmentFile(str(bam), "rb") as alignments:
        flank = next(
            rec
            for rec in alignments.fetch(until_eof=True)
            if rec.query_name == "r1/1:end:5"
        )
    assert flank.is_paired and flank.reference_start > 1000
    assert fullreads_bam is not None
    with pysam.AlignmentFile(str(fullreads_bam), "rb") as fullread_alignments:
        fullread = next(
            rec for rec in fullread_alignments.fetch(until_eof=True) if rec.is_read1
        )
    assert fullread.is_paired and fullread.reference_start > 1000


@pytest.mark.skipif(shutil.which("bwa") is None, reason="bwa not available")
def test_align_genome_subcommand_consumes_all_hit_artifact(tmp_path: Path):
    from RelocaTE3.cli import main

    rng = random.Random(7)
    genome_bases = [rng.choice("ACGT") for _ in range(2200)]
    block = "".join(rng.choice("ACGT") for _ in range(100))
    genome_bases[500:600] = list(block)
    genome_bases[1500:1600] = list(block)
    anchor = "".join(genome_bases[1650:1750])
    genome = tmp_path / "genome.fa"
    genome.write_text(f">chr1\n{''.join(genome_bases)}\n")
    (tmp_path / "flanking").mkdir()
    flank = tmp_path / "flanking" / "S.left.flankingReads.fq"
    _write_fastq(flank, [("r1/1:end:5", block)])
    (tmp_path / "te_containing").mkdir()
    (tmp_path / "te_containing" / "S.te_hit_names.txt").write_text("r1/1\n")
    r1, r2 = tmp_path / "reads_1.fq", tmp_path / "reads_2.fq"
    _write_fastq(r1, [("r1/1", block)])
    _write_fastq(r2, [("r1/2", _rc(anchor))])

    assert (
        main(
            [
                "align-genome",
                "-g",
                str(genome),
                "-f",
                str(flank),
                "-n",
                "S",
                "-o",
                str(tmp_path),
                "--genome-aligner",
                "bwa",
                "--threads",
                "1",
                "-1",
                str(r1),
                "-2",
                str(r2),
            ]
        )
        == 0
    )
    with pysam.AlignmentFile(str(tmp_path / "S.repeat.bwa.sorted.bam"), "rb") as bam:
        record = next(
            rec for rec in bam.fetch(until_eof=True) if ":end:5" in rec.query_name
        )
    assert record.is_paired and record.reference_start > 1000

    fullreads = tmp_path / "genome_aln" / "S.fullreads.genome.bam"
    assert fullreads.is_file() and fullreads.with_suffix(".bam.bai").is_file()
    with pysam.AlignmentFile(str(fullreads), "rb") as bam:
        records = list(bam.fetch(until_eof=True))
    assert len(records) == 2
    junction_fullread = next(rec for rec in records if rec.is_read1)
    assert junction_fullread.is_paired
    assert junction_fullread.query_sequence == block
    assert junction_fullread.reference_start > 1000


@pytest.mark.skipif(shutil.which("bwa") is None, reason="bwa not available")
def test_align_flanks_anchored_moves_single_bam_across_devices(
    tmp_path: Path, monkeypatch
):
    rng = random.Random(1)
    genome = tmp_path / "genome.fa"
    genome.write_text(
        ">chr1\n" + "".join(rng.choice("ACGT") for _ in range(1000)) + "\n"
    )
    flank = tmp_path / "flank.fq"
    _write_fastq(flank, [("read/1:end:5", "ACGTACGTACGTACGT")])
    reads1, reads2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    _write_fastq(reads1, [("read/1", "ACGTACGTACGTACGT")])
    _write_fastq(reads2, [("read/2", "TGCATGCATGCATGCA")])
    reads = ReadLibrary([str(reads1), str(reads2)], "S")
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    out_bam = tmp_path / "out" / "S.repeat.bwa.sorted.bam"
    real_replace = Path.replace

    def _exdev(self, target):
        if Path(target) == out_bam:
            raise OSError(errno.EXDEV, "Invalid cross-device link")
        return real_replace(self, target)

    monkeypatch.setattr(Path, "replace", _exdev)
    result = align_flanks_anchored(
        get_aligner("bwa", 1),
        str(genome),
        [str(flank)],
        {"read/1", "read/2"},
        reads,
        out_bam,
        threads=1,
        tmp=str(scratch),
    )
    assert Path(result).exists() and Path(result).stat().st_size > 0
