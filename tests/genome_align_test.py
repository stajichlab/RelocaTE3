"""Unit tests for genome re-alignment helpers (Step 4)."""

from __future__ import annotations

import random
import shutil
from pathlib import Path

import pysam
import pytest

from RelocaTE3.genome_align import (
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
