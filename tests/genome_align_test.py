"""Unit tests for genome re-alignment helpers (Step 4)."""

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
    build_flank_pairs,
    collect_junction_fullreads,
    recover_support_mates,
    split_mate,
    strip_tag,
)
from RelocaTE3.ReadLibrary import ReadLibrary


def _rc(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

DATA = Path(__file__).parent / "data"
R1 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"
R2 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"


def test_strip_tag():
    assert strip_tag("read_500_1/1:end:5") == "read_500_1/1"
    assert strip_tag("read_500_1/2:start:3") == "read_500_1/2"
    assert strip_tag("read_500_1/1:middle") == "read_500_1/1"
    assert strip_tag("M00:1:2:3/1") == "M00:1:2:3/1"  # real colons preserved


def test_split_mate():
    assert split_mate("read_500_1/1") == ("read_500_1", "1")
    assert split_mate("read_500_1/2") == ("read_500_1", "2")
    assert split_mate("nomate") == ("nomate", "")


def test_recover_support_mates_pulls_only_unmatched_mates(tmp_path: Path):
    reads = ReadLibrary([str(R1), str(R2)], "HEG4")
    # read_500_1/1 matched the TE; its mate /2 did not -> /2 is a supporting read.
    # read_500_2/1 and /2 both matched -> no support read for that pair.
    read_repeat = {
        "read_500_1/1:end:5": ("mPing", "+"),
        "read_500_2/1:end:5": ("mPing", "+"),
        "read_500_2/2:start:3": ("mPing", "-"),
    }
    out_fq = tmp_path / "support.fq"
    n = recover_support_mates(read_repeat, reads, out_fq)
    assert n == 1
    contents = out_fq.read_text()
    assert "read_500_1/2" in contents
    assert "read_500_2/2" not in contents  # mate also matched the TE


def test_recover_support_mates_single_end(tmp_path: Path):
    reads = ReadLibrary([str(R1)], "HEG4")
    out_fq = tmp_path / "support.fq"
    assert (
        recover_support_mates({"read_500_1/1:end:5": ("mPing", "+")}, reads, out_fq)
        == 0
    )


def test_build_flank_pairs_pairs_flank_with_genomic_mate(tmp_path: Path):
    """A junction flank whose mate did NOT match the TE is paired with that mate
    (flank -> R1 keeping its tag; genomic mate -> R2), so the aligner can anchor
    an ambiguous flank to the mate's unique locus."""
    reads = ReadLibrary([str(R1), str(R2)], "HEG4")
    flank_fq = tmp_path / "HEG4.left.flankingReads.fq"
    flank_fq.write_text("@read_500_1/1:end:5\nACGTACGTAC\n+\nIIIIIIIIII\n")
    read_repeat = {"read_500_1/1:end:5": ("mPing", "+")}  # /1 matched; /2 is genomic mate
    r1, r2, se = tmp_path / "r1.fq", tmp_path / "r2.fq", tmp_path / "se.fq"
    n_pair, n_se, retag = build_flank_pairs([str(flank_fq)], read_repeat, reads, r1, r2, se)
    assert (n_pair, n_se) == (1, 0)
    # aligned with neutral matching names; real names restored post-alignment
    assert "@read_500_1/1\n" in r1.read_text()
    assert "@read_500_1/2\n" in r2.read_text()
    assert retag == {"read_500_1": ("read_500_1/1:end:5", "read_500_1/2")}
    assert se.read_text() == ""


def test_build_flank_pairs_se_when_mate_also_matched_te(tmp_path: Path):
    """A junction flank whose mate ALSO matched the TE has no genomic mate to
    anchor it, so it falls back to the single-end file (unchanged behavior)."""
    reads = ReadLibrary([str(R1), str(R2)], "HEG4")
    flank_fq = tmp_path / "HEG4.left.flankingReads.fq"
    flank_fq.write_text("@read_500_2/1:end:5\nACGTACGTAC\n+\nIIIIIIIIII\n")
    read_repeat = {
        "read_500_2/1:end:5": ("mPing", "+"),
        "read_500_2/2:start:3": ("mPing", "-"),  # mate also matched -> no anchor
    }
    r1, r2, se = tmp_path / "r1.fq", tmp_path / "r2.fq", tmp_path / "se.fq"
    n_pair, n_se, retag = build_flank_pairs([str(flank_fq)], read_repeat, reads, r1, r2, se)
    assert (n_pair, n_se) == (0, 1)
    assert retag == {}
    assert "read_500_2/1:end:5" in se.read_text()
    assert r1.read_text() == "" and r2.read_text() == ""


def test_build_flank_pairs_skips_te_internal_middle_reads(tmp_path: Path):
    """``:middle`` reads must never reach the genome aligner (RelocaTE2 parity).

    A ``:middle`` read lies entirely inside the TE, so it carries no flank and
    cannot mark a breakpoint -- it is pure transposon sequence that maps to every
    reference copy of its family. RelocaTE2 drops these before alignment and
    keeps only their genomic mates (``clean_pairs_memory.py``: "...but not reads
    themselve as they are part of repeat").

    RelocaTE3 used to align them: on riceTElib cov30x_rep1, 2,546,333 of
    3,942,639 genome-BAM records (64.6%) were ``:middle`` where RelocaTE2's
    equivalent inputs held 0. In ``_stream_clusters`` each one extends its
    cluster and counts as a supporting read, so they glued unrelated breakpoints
    together and inflated support at every reference TE copy.
    """
    reads = ReadLibrary([str(R1), str(R2)], "HEG4")
    flank_fq = tmp_path / "HEG4.left.flankingReads.fq"
    flank_fq.write_text(
        "@read_500_1/1:end:5\nACGTACGTAC\n+\nIIIIIIIIII\n"      # junction: keep
        "@read_500_3/1:middle\nTTTTTTTTTT\n+\nIIIIIIIIII\n"      # TE-internal: drop
    )
    read_repeat = {
        "read_500_1/1:end:5": ("mPing", "+"),
        "read_500_3/1:middle": ("mPing", "+"),
    }
    r1, r2, se = tmp_path / "r1.fq", tmp_path / "r2.fq", tmp_path / "se.fq"
    n_pair, n_se, retag = build_flank_pairs(
        [str(flank_fq)], read_repeat, reads, r1, r2, se
    )
    written = r1.read_text() + r2.read_text() + se.read_text()
    assert ":middle" not in written
    assert "TTTTTTTTTT" not in written, "the TE-internal sequence must not be aligned"
    # only the junction flank survives, still paired with its genomic mate
    assert (n_pair, n_se) == (1, 0)
    assert retag == {"read_500_1": ("read_500_1/1:end:5", "read_500_1/2")}


def test_middle_read_mate_is_still_recovered_as_support(tmp_path: Path):
    """Dropping the middle read must not drop its genomic mate.

    RelocaTE2 keeps that mate -- it lands in unique genome sequence and brackets
    the insertion, which is exactly what a supporting read is.
    ``recover_support_mates`` reads the same ``read_repeat`` table, so the mate is
    picked up there once the middle read is no longer paired with it.
    """
    reads = ReadLibrary([str(R1), str(R2)], "HEG4")
    read_repeat = {"read_500_3/1:middle": ("mPing", "+")}
    out_fq = tmp_path / "support.fq"
    n = recover_support_mates(read_repeat, reads, out_fq, exclude=set())
    assert n == 1, "the middle read's genomic mate must still be aligned"
    assert "read_500_3/2" in out_fq.read_text()


@pytest.mark.skipif(shutil.which("bwa") is None, reason="bwa not available")
def test_align_to_genome_anchors_ambiguous_flank_to_mate(tmp_path: Path):
    """An ambiguous junction flank (maps to two identical loci) must be placed at
    the locus its uniquely-mapping mate anchors, not scattered to the other copy.

    Genome has a duplicated 100 bp block at pos 500 and pos 1500. The flank equals
    that block (ambiguous). Its mate is unique and sits just past the 1500 copy.
    Single-end, the flank lands on the first (500) copy; paired with its mate it is
    anchored to the true (1500) copy.
    """
    rng = random.Random(7)
    g = [rng.choice("ACGT") for _ in range(2200)]
    block = "".join(rng.choice("ACGT") for _ in range(100))  # the duplicated repeat
    g[500:600] = list(block)
    g[1500:1600] = list(block)
    anchor = "".join(g[1650:1750])  # unique sequence right of the 1500 copy
    genome_seq = "".join(g)

    gdir = tmp_path / "g"
    gdir.mkdir()
    genome = gdir / "g.fa"
    genome.write_text(f">chr1\n{genome_seq}\n")

    out = tmp_path / "out"
    (out / "flanking").mkdir(parents=True)
    (out / "flanking" / "S.left.flankingReads.fq").write_text(
        f"@r1/1:end:5\n{block}\n+\n{'I' * 100}\n"
    )
    (out / "te_containing").mkdir(parents=True)
    (out / "te_containing" / "S.read_repeat_name.txt").write_text("r1/1:end:5\tmPing\t+\n")

    # ReadLibrary whose R2 holds the flank's genomic mate (reverse-comp of the
    # unique anchor, so it maps in proper FR orientation just past pos 1500).
    r1fq = tmp_path / "reads_1.fq"
    r2fq = tmp_path / "reads_2.fq"
    r1fq.write_text(f"@r1/1\n{block}\n+\n{'I' * 100}\n")
    r2fq.write_text(f"@r1/2\n{_rc(anchor)}\n+\n{'I' * 100}\n")
    reads = ReadLibrary([str(r1fq), str(r2fq)], "S")

    bam, _ = align_to_genome(reads, str(genome), out, threads=1, genome_aligner="bwa")
    with pysam.AlignmentFile(str(bam), "rb") as bf:
        flank = next(
            (r for r in bf.fetch(until_eof=True) if r.query_name == "r1/1:end:5"), None
        )
    assert flank is not None and not flank.is_unmapped
    # The fix aligns the flank paired with its mate: single-end it carries no mate
    # (is_paired False) and its placement among identical copies is unanchored.
    assert flank.is_paired, "flank was not aligned paired with its mate"
    # and the unique mate anchors it to the 1500 copy, not the 500 copy
    assert flank.reference_start > 1000, (
        f"flank landed at {flank.reference_start}; expected near 1500 (mate-anchored)"
    )


@pytest.mark.skipif(shutil.which("bwa") is None, reason="bwa not available")
def test_align_genome_subcommand_anchors_with_mate(tmp_path: Path):
    """The `align-genome` CLI subcommand (used by the benchmark) must mate-anchor
    when given the original reads (-1/-2), matching align_to_genome."""
    from RelocaTE3.cli import main

    rng = random.Random(7)
    g = [rng.choice("ACGT") for _ in range(2200)]
    block = "".join(rng.choice("ACGT") for _ in range(100))
    g[500:600] = list(block)
    g[1500:1600] = list(block)
    anchor = "".join(g[1650:1750])
    genome = tmp_path / "g.fa"
    genome.write_text(f">chr1\n{''.join(g)}\n")

    (tmp_path / "flanking").mkdir()
    flank = tmp_path / "flanking" / "S.left.flankingReads.fq"
    flank.write_text(f"@r1/1:end:5\n{block}\n+\n{'I' * 100}\n")
    (tmp_path / "te_containing").mkdir()
    (tmp_path / "te_containing" / "S.read_repeat_name.txt").write_text(
        "r1/1:end:5\tmPing\t+\n"
    )
    r1fq = tmp_path / "reads_1.fq"
    r2fq = tmp_path / "reads_2.fq"
    r1fq.write_text(f"@r1/1\n{block}\n+\n{'I' * 100}\n")
    r2fq.write_text(f"@r1/2\n{_rc(anchor)}\n+\n{'I' * 100}\n")

    rc = main([
        "align-genome", "-g", str(genome), "-f", str(flank),
        "-n", "S", "-o", str(tmp_path), "--genome-aligner", "bwa", "--threads", "1",
        "-1", str(r1fq), "-2", str(r2fq),
    ])
    assert rc == 0
    bam = tmp_path / "S.repeat.bwa.sorted.bam"
    with pysam.AlignmentFile(str(bam), "rb") as bf:
        f = next((r for r in bf.fetch(until_eof=True) if r.query_name == "r1/1:end:5"), None)
    assert f is not None and not f.is_unmapped
    assert f.is_paired, "align-genome subcommand did not mate-anchor the flank"
    assert f.reference_start > 1000  # anchored to the true (1500) copy


@pytest.mark.skipif(shutil.which("bwa") is None, reason="bwa not available")
def test_align_flanks_anchored_moves_single_bam_across_devices(
    tmp_path: Path, monkeypatch
):
    """When only one part BAM is produced (all flanks single-end, e.g. the bwa
    TE-aligner strips mate suffixes), the result must be moved to out_bam even
    when tmp and out_bam are on different filesystems. Path.replace() raises
    EXDEV across devices; align_flanks_anchored must move, not rename.
    """
    rng = random.Random(1)
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\n" + "".join(rng.choice("ACGT") for _ in range(1000)) + "\n")

    flank = tmp_path / "S.left.flankingReads.fq"
    flank.write_text("@read_500_2/1:end:5\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n")
    # mate ALSO matched the TE -> flank goes single-end -> exactly one part BAM
    read_repeat = {
        "read_500_2/1:end:5": ("mPing", "+"),
        "read_500_2/2:start:3": ("mPing", "-"),
    }
    reads = ReadLibrary([str(R1), str(R2)], "S")

    tmp = tmp_path / "scratch"
    tmp.mkdir()
    out_bam = tmp_path / "out" / "S.repeat.bwa.sorted.bam"

    # Simulate a cross-device move only for the final out_bam rename.
    real_replace = Path.replace

    def _exdev(self, target):
        if Path(target) == out_bam:
            raise OSError(errno.EXDEV, "Invalid cross-device link")
        return real_replace(self, target)

    monkeypatch.setattr(Path, "replace", _exdev)

    result = align_flanks_anchored(
        get_aligner("bwa", 1), str(genome), [str(flank)], read_repeat, reads,
        out_bam, threads=1, tmp=str(tmp),
    )
    assert Path(result).exists() and Path(result).stat().st_size > 0


def test_collect_junction_fullreads(tmp_path: Path):
    """Full (untrimmed) sequences are pulled only for 5'/3' junction reads."""
    reads = ReadLibrary([str(R1), str(R2)], "HEG4")
    read_repeat = {
        "read_500_1/1:end:5": ("mPing", "+"),  # junction -> pulled
        "read_500_2/2:start:3": ("mPing", "-"),  # junction -> pulled
        "read_500_3/1:middle": ("mPing", "+"),  # middle -> not pulled
    }
    out_fq = tmp_path / "full.fq"
    n = collect_junction_fullreads(read_repeat, reads, out_fq)
    assert n == 2
    contents = out_fq.read_text()
    assert "read_500_1/1" in contents
    assert "read_500_2/2" in contents
    assert "read_500_3/1" not in contents
