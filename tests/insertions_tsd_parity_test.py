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
