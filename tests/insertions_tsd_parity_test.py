"""Tests for the depth-based TSD length/sequence inference (R2 parity)."""

from unittest.mock import MagicMock

from RelocaTE3.insertions import (
    _Cluster,
    _capture_tsd_from_read,
    _estimate_tsd_length_from_depth,
    _make_insertion,
)
from RelocaTE3.models import JunctionObservation


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
    assert ins.te_name == "mPing"
    assert ins.te_family_support == {"mPing": 3}
    assert ins.te_family_confidence == 1.0
    assert ins.te_family_status == "unique"
    genome.fetch.assert_not_called()


# ---------------------------------------------------------------------------
# RelocaTE2's geometric TSD estimate and its emission gate
# ---------------------------------------------------------------------------
#
# relocaTE_insertionFinder.py:800-818 computes the TSD length twice --
# tsd_finder (depth) and TSD_len_calculate (geometry of the dominant left/right
# breakpoints) -- then:
#
#     if not int(TSD_len) == int(TSD_len_1): TSD_len = int(TSD_len_1)
#     if TSD_len > 0:   # ... otherwise the sub-cluster produces no call
#
# so the geometric value always wins, and a non-positive one means the pairing
# is impossible and RelocaTE2 reports nothing. RelocaTE3 used to fall back to
# the depth estimate in that case and emit the call anyway, which was the
# dominant false-positive source on multi-family libraries (riceTElib
# cov30x_rep1: 944 two-sided calls vs RelocaTE2's 357, from fewer TE reads).


def _two_sided_cluster(left_pos, right_pos):
    cluster = _Cluster("Chr1")
    cluster.junctions = [
        _mk_junction("l1", "left", left_pos, "+", "AAAAATTA"),
        _mk_junction("r1", "right", right_pos, "+", "TTAGGGCC"),
    ]
    for j in cluster.junctions:
        j.gstart, j.gend = right_pos, left_pos
    cluster.support = []
    cluster.extend(min(left_pos, right_pos) - 10, max(left_pos, right_pos) + 10)
    return cluster


def _call(left_pos, right_pos):
    cluster = _two_sided_cluster(left_pos, right_pos)
    left = [j for j in cluster.junctions if j.side == "left"]
    right = [j for j in cluster.junctions if j.side == "right"]
    return _make_insertion("Chr1", left, right, None, cluster)


def test_geometric_tsd_length_is_used_for_a_two_sided_call():
    """left breakpoint 102, right 100 -> a 3bp TSD spanning 100..102."""
    ins = _call(left_pos=102, right_pos=100)
    assert ins is not None
    assert (ins.start, ins.end) == (100, 102)
    assert ins.tsd == "TTA"


def test_short_junction_is_removed_before_fullread_filtering():
    """Reproduce the riceTElib Chr1:22171544 RelocaTE2-only rejection.

    The paired breakpoints imply a 40 bp TSD, but the left flank is only 29 bp.
    RelocaTE2's TSD_check_cluster cannot match its 40-base wildcard against that
    flank, so only the 88 bp right junction reaches the full-read filter.
    """
    cluster = _Cluster("Chr1")
    left = JunctionObservation(
        "cov5x_rep1:clone10:Chr1-37415/1:end:5",
        "left",
        22171583,
        "+",
        "Os1279#DNAauto/CACTA",
        "5",
        "A" * 29,
        22171555,
        22171583,
    )
    right = JunctionObservation(
        "cov5x_rep1:baseline:Chr1-182527/1:start:5",
        "right",
        22171544,
        "+",
        "Os0029#MITE/Stow",
        "5",
        "T" * 88,
        22171544,
        22171631,
    )
    cluster.junctions = [left, right]
    cluster.extend(22171544, 22171631)

    ins = _make_insertion("Chr1", [left], [right], None, cluster)

    assert ins is not None
    assert (ins.start, ins.end, ins.tsd) == (22171544, 22171583, "UNK")
    assert (ins.left_junction_reads, ins.right_junction_reads) == (0, 1)
    assert ins.read_names == [right.read_name]
    assert ins.te_name == "Os0029#MITE/Stow"
    assert ins.strand == "-"


def test_impossible_geometry_emits_no_call():
    """Right breakpoint past the left one: RelocaTE2 emits nothing (:818)."""
    assert _call(left_pos=100, right_pos=140) is None


def test_exactly_adjacent_breakpoints_are_a_zero_length_tsd_and_are_dropped():
    """left == right - 1 gives TSD_len_1 == 0, which fails RelocaTE2's `> 0`."""
    assert _call(left_pos=99, right_pos=100) is None


def test_one_sided_calls_are_unaffected_by_the_gate():
    """TSD_len_calculate needs both sides; one-sided calls keep the depth path."""
    cluster = _Cluster("Chr1")
    cluster.junctions = [
        _mk_junction(f"r{i}", "right", 100, "+", "TTAGGGCCAAAA") for i in range(3)
    ]
    cluster.support = [("s1", 95, 102, "+", "N" * 8)]
    cluster.extend(95, 110)
    ins = _make_insertion(
        "Chr1", [], list(cluster.junctions), None, cluster
    )
    assert ins is not None
