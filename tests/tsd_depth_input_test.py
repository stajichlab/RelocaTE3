"""TSD length must be estimated from junction reads, not supporting mates.

RelocaTE2 builds its depth pileup by calling ``calculate_cluster_depth`` once
per *junction* read in the sub-cluster, and divides by that sub-cluster's
junction-read count (relocaTE_insertionFinder.py:1069-1076, 843-850). The
insight is that junction reads all abut the same breakpoint, so the TSD is the
run of positions where essentially every one of them overlaps.

RelocaTE3 passed ``cluster.support`` instead -- the supporting *mate* reads.
Mates are scattered across the library insert (~500 bp), so no single base is
covered by 60% of them, ``tsd_finder`` returns 0, and the call is reported as
``UNK``. On the Chr3 2 Mb fixture that accounted for 15 of the 16 sites where
RelocaTE2 resolves a concrete TSD and RelocaTE3 does not.
"""

from __future__ import annotations

from RelocaTE3.insertions import _Cluster, _estimate_tsd_length_from_depth
from RelocaTE3.models import JunctionObservation


def _junction(name, side, position, gstart, gend, seq="ACGTACGTAC"):
    return JunctionObservation(
        name, side, position, "+", "mping", "5", seq, gstart, gend
    )


def test_junction_reads_abutting_a_breakpoint_yield_the_tsd_length():
    """Five junction reads whose spans share exactly 3 bases -> TSD length 3."""
    spans = [(100, 102), (100, 102), (100, 102), (100, 102), (100, 102)]
    assert _estimate_tsd_length_from_depth(spans, 100) == 3


def test_scattered_mate_spans_yield_nothing():
    """The old input: mates spread over an insert never reach the depth cutoff."""
    spans = [(100 + 40 * i, 200 + 40 * i) for i in range(10)]
    assert _estimate_tsd_length_from_depth(spans, 300) == 0


def test_lower_thresholds_are_tried_in_order():
    """RelocaTE2 falls back 1.0 -> 0.8 -> 0.6 (insertionFinder.py:801-806).

    Note the fallback only happens when the higher threshold yields *zero*: at
    1.0 a single fully-covered base already returns a length of 1 and wins.
    Here the outlier read shares no base with the rest, so 1.0 yields nothing
    and 0.8 (4 of 5) resolves the 3 bp footprint.
    """
    spans = [(100, 102)] * 4 + [(200, 202)]
    assert _estimate_tsd_length_from_depth(spans, 100) == 3


def test_full_depth_wins_before_the_fallback():
    """A base covered by every read returns length 1 rather than falling back."""
    spans = [(100, 102)] * 4 + [(100, 100)]
    assert _estimate_tsd_length_from_depth(spans, 100) == 1


def test_cluster_junction_spans_are_recorded():
    """JunctionObservation must carry the aligned span, or the pileup is impossible."""
    cluster = _Cluster("Chr3")
    cluster.junctions.append(_junction("r1", "right", 100, 100, 140))
    obs = cluster.junctions[0]
    assert (obs.gstart, obs.gend) == (100, 140)


def test_make_insertion_uses_junction_spans_not_support_spans():
    """End to end: a cluster whose mates are scattered but whose junctions agree.

    Before the fix this produced tsd='UNK' because the scattered mates drove the
    depth estimate to zero.
    """
    from RelocaTE3.insertions import _make_insertion

    cluster = _Cluster("Chr3")
    # Four right-side junction reads sharing a 3 bp footprint at 100..102.
    for i in range(4):
        cluster.junctions.append(_junction(f"j{i}", "right", 100, 100, 102, "TTAGGGCCC"))
    # Mates scattered across an insert: these must not drive the estimate.
    for i in range(12):
        cluster.support.append((f"s{i}", 50 + 30 * i, 150 + 30 * i, "+", "A" * 50))
    cluster.extend(50, 600)

    ins = _make_insertion("Chr3", [], list(cluster.junctions), None, cluster)
    assert ins.tsd != "UNK", "TSD should be resolved from the junction reads"
    assert len(ins.tsd) == 3, f"expected a 3 bp TSD, got {ins.tsd!r}"
