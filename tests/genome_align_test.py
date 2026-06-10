"""Unit tests for genome re-alignment helpers (Step 4)."""

from __future__ import annotations

from pathlib import Path

from RelocaTE3.genome_align import (
    collect_junction_fullreads,
    recover_support_mates,
    split_mate,
    strip_tag,
)
from RelocaTE3.ReadLibrary import ReadLibrary

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
