"""Tests for the depth-based TSD length/sequence inference (R2 parity)."""

from RelocaTE3.insertions import (
    _capture_tsd_from_read,
    _estimate_tsd_length_from_depth,
)


def test_depth_estimator_returns_zero_when_no_reads():
    assert _estimate_tsd_length_from_depth([], breakpoint=100) == 0


def test_depth_estimator_matches_overlap_width():
    # Five reads all overlapping positions 100..102 inclusive (3bp), then tapering.
    reads = [(95, 102), (96, 102), (97, 102), (98, 102), (100, 110)]
    # depth at 100,101,102 == 5; threshold 1.0 yields width 3
    assert _estimate_tsd_length_from_depth(reads, breakpoint=100) == 3


def test_capture_tsd_from_left_read_uses_last_bases():
    # left-junction read: TSD is the last `length` chars of the seq
    assert _capture_tsd_from_read("AAAAATAA", side="left", length=3) == "TAA"


def test_capture_tsd_from_right_read_uses_first_bases():
    # right-junction read: TSD is the first `length` chars of the seq
    assert _capture_tsd_from_read("TTAGGGCC", side="right", length=3) == "TTA"


def test_capture_tsd_short_seq_returns_empty():
    assert _capture_tsd_from_read("AT", side="right", length=3) == ""


from unittest.mock import MagicMock

from RelocaTE3.insertions import _Cluster, _make_insertion
from RelocaTE3.models import JunctionObservation


def _mk_junction(name, side, pos, strand, seq, te_end="5"):
    return JunctionObservation(
        read_name=name,
        side=side,
        position=pos,
        strand=strand,
        te_name="mPing",
        te_end=te_end,
        seq=seq,
    )


def test_make_insertion_one_sided_recovers_tsd_from_depth_and_read():
    chrom = "Chr1"
    cluster = _Cluster(chrom)
    # Three right-junction reads at position 100; all sequences start with "TTA".
    cluster.junctions = [
        _mk_junction(f"r{i}", "right", 100, "+", "TTAGGGCCAAAA") for i in range(3)
    ]
    # Five supporting reads tapering around positions 100..102 so the depth
    # peak (full coverage) is exactly 3bp wide, matching the unit-test fixture.
    cluster.support = [
        ("s1", 95, 102, "+", "NNNNNNNN"),
        ("s2", 96, 102, "+", "NNNNNNN"),
        ("s3", 97, 102, "+", "NNNNNN"),
        ("s4", 98, 102, "+", "NNNNN"),
        ("s5", 100, 110, "+", "NNNNNNNNNNN"),
    ]
    cluster.extend(95, 110)

    genome = MagicMock()  # should NOT be consulted when read capture succeeds
    ins = _make_insertion(
        chrom,
        left_reads=[],
        right_reads=cluster.junctions,
        genome=genome,
        cluster=cluster,
    )
    assert ins.tsd == "TTA"
    genome.fetch.assert_not_called()
