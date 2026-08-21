"""Breakpoint pairing must follow RelocaTE2's sub-clustering exactly.

``_pair_breakpoints`` is a direct port of ``TSD_from_read_depth``
(relocaTE_insertionFinder.py:603-770). RelocaTE2 walks each *left* breakpoint in
ascending order and pairs it with its **nearest** right breakpoint, accepting the
pair only when they are within ``PAIR_MAX_DISTANCE`` (100 bp, :643). Right
breakpoints no left ever claimed are then emitted one-sided (:678).

RelocaTE3 previously used its own scheme -- chain positions into groups by a
25 bp gap, take the most-supported position on each side, then split the pair
again if it looked geometrically impossible. That is not what RelocaTE2 does, and
a chain of <=25 bp steps could span far more than 25 bp, pairing a dominant left
with a distant dominant right.

Note the impossible-pairing case is still handled -- just where RelocaTE2 handles
it. RelocaTE2 pairs the breakpoints here and then rejects the sub-insertion at
emission, because ``TSD_len_calculate`` returns a non-positive length and
``if TSD_len > 0`` (:818) suppresses the call. See ``_make_insertion``.
"""

from __future__ import annotations

from RelocaTE3.insertions import PAIR_MAX_DISTANCE, _pair_breakpoints


def _obs(n: int) -> list[object]:
    """n placeholder junction observations -- only the count is read here."""
    return [object()] * n


def test_nearest_right_wins_not_the_deepest():
    """RelocaTE2 pairs on distance, not support (:610-642)."""
    # 1003 is nearer to 1005 than 1050 is, even though both have equal depth
    assert _pair_breakpoints({1005: _obs(4)}, {1003: _obs(4)}) == [(1005, 1003)]


def test_coincident_breakpoints_pair():
    """A zero-length TSD (left == right) is still one insertion."""
    assert _pair_breakpoints({500: _obs(2)}, {500: _obs(2)}) == [(500, 500)]


def test_lone_pair_beyond_the_window_is_split():
    """":709 -- two junctions are far from each other, might be one end from two
    insertion"."""
    pairs = _pair_breakpoints({1000: _obs(2)}, {2000: _obs(2)})
    assert sorted(pairs, key=lambda p: (p[0] or 0, p[1] or 0)) == [
        (None, 2000),
        (1000, None),
    ]


def test_window_boundary_is_inclusive():
    """`min_dist <= 100` pairs; one bp further splits."""
    d = PAIR_MAX_DISTANCE
    assert _pair_breakpoints({1000: _obs(2)}, {1000 + d: _obs(2)}) == [(1000, 1000 + d)]
    assert len(_pair_breakpoints({1000: _obs(2)}, {1000 + d + 1: _obs(2)})) == 2


def test_one_sided_clusters_are_unchanged():
    assert _pair_breakpoints({700: _obs(3)}, {}) == [(700, None)]
    assert _pair_breakpoints({}, {700: _obs(3)}) == [(None, 700)]


def test_every_left_gets_its_own_sub_insertion():
    """RelocaTE2 emits one sub-cluster per left position, even when two lefts
    share the same nearest right -- it does not collapse them (:610-676)."""
    pairs = _pair_breakpoints({1005: _obs(1), 1006: _obs(7)}, {1003: _obs(4)})
    assert pairs == [(1005, 1003), (1006, 1003)]


def test_unclaimed_right_breakpoints_are_emitted_one_sided():
    """:678 -- rights that no left paired with become their own sub-clusters."""
    # 1003 is nearest to both lefts; 5000 is claimed by nobody
    pairs = _pair_breakpoints({1005: _obs(2), 1006: _obs(2)}, {1003: _obs(3), 5000: _obs(3)})
    assert (None, 5000) in pairs
    assert (1005, 1003) in pairs and (1006, 1003) in pairs


def test_adjacent_insertions_resolve_when_both_are_present():
    """The Chr3 case, with the full breakpoint set RelocaTE2 actually sees.

    truth Chr3:155988  RelocaTE2: 155951..155962 and 155984..155988

    Each left pairs with its own nearest right, so the two insertions separate
    naturally -- no geometric special-case needed.
    """
    left = {155962: _obs(3), 155988: _obs(4)}
    right = {155951: _obs(3), 155984: _obs(5)}
    pairs = _pair_breakpoints(left, right)
    assert pairs == [(155962, 155951), (155988, 155984)]


def test_several_positions_on_one_side_only():
    """:687/:697 -- each becomes its own one-sided sub-insertion."""
    assert _pair_breakpoints({100: _obs(2), 400: _obs(2)}, {}) == [
        (100, None),
        (400, None),
    ]
    assert _pair_breakpoints({}, {100: _obs(2), 400: _obs(2)}) == [
        (None, 100),
        (None, 400),
    ]


# ---------------------------------------------------------------------------
# Reference-TE edge exclusion must not discard two-sided calls
# ---------------------------------------------------------------------------
#
# RelocaTE2 only removes calls near a reference TE boundary when a junction side
# is empty (clean_false_positive.py:82, `Right == 0 or Left == 0`). A call with
# junction reads on both sides is real evidence and is kept even when it abuts a
# reference copy -- transposons insert next to transposons.
#
# RelocaTE3 dropped any call whose edge coincided with a stored reference edge,
# regardless of support. On the Chr3 2 Mb fixture that discarded the mPing
# insertion at 257446..257448 -- 3 left and 6 right junction reads, TSD ACG,
# exactly what RelocaTE2 reports -- because a TEOS1 copy ends at 257444 and the
# loader stores a small window of end positions around it.

from RelocaTE3.insertions import _excluded_by_reference_edge
from RelocaTE3.models import Insertion


def _ins(start, end, left, right):
    return Insertion(
        chrom="Chr3", start=start, end=end, te_name="mPing", strand="+", tsd="ACG",
        left_junction_reads=left, right_junction_reads=right,
    )


EDGES = {"start": {}, "end": {257442: 1, 257443: 1, 257444: 1, 257445: 1}}


def test_two_sided_call_at_a_reference_edge_is_kept():
    """The Chr3:257448 regression: 3 left + 6 right reads, dropped for touching an edge."""
    assert not _excluded_by_reference_edge(_ins(257446, 257448, 3, 6), EDGES)


def test_one_sided_call_at_a_reference_edge_is_still_excluded():
    """RelocaTE2 does remove these -- reads off an intact copy's edge."""
    assert _excluded_by_reference_edge(_ins(257446, 257448, 0, 4), EDGES)
    assert _excluded_by_reference_edge(_ins(257446, 257448, 4, 0), EDGES)


def test_one_sided_call_away_from_any_edge_is_kept():
    assert not _excluded_by_reference_edge(_ins(900000, 900002, 0, 4), EDGES)


def test_start_edge_is_checked_too():
    edges = {"start": {5000: 1}, "end": {}}
    assert _excluded_by_reference_edge(_ins(4998, 5000, 0, 3), edges)
    assert not _excluded_by_reference_edge(_ins(4998, 5000, 3, 3), edges)


# ---------------------------------------------------------------------------
# Cluster-level arbitration (write_output:257-330)
# ---------------------------------------------------------------------------
#
# RelocaTE2 generates many candidates per cluster and then weighs them against
# each other. Porting _pair_breakpoints without this filter raised false
# positives while leaving recall untouched (mPing precision 1.000 -> 0.976,
# riceTElib 5x F1 0.409 -> 0.360): pairing and arbitration are one mechanism.

from RelocaTE3.insertions import MIN_ONE_SIDED_JUNCTIONS, _arbitrate_cluster

NO_EDGES: dict = {"start": {}, "end": {}}


def _cand(left: int, right: int, start: int = 1000) -> Insertion:
    return Insertion(
        chrom="Chr1", start=start, end=start + 2, te_name="mPing", strand="+",
        tsd="TTA", left_junction_reads=left, right_junction_reads=right,
    )


def test_single_two_sided_candidate_is_reported():
    ins = _cand(3, 4)
    assert _arbitrate_cluster([ins], NO_EDGES) == [ins]


def test_single_one_sided_candidate_at_a_reference_edge_is_dropped():
    """:311 else-branch -- a lone one-sided call is checked against existingTE."""
    edges = {"start": {1002: 1}, "end": {}}
    assert _arbitrate_cluster([_cand(0, 4)], edges) == []
    assert _arbitrate_cluster([_cand(0, 4)], NO_EDGES) != []


def test_all_two_sided_candidates_are_kept():
    """:259 -- if every candidate has both junctions, keep them all."""
    cands = [_cand(2, 2, 1000), _cand(3, 1, 2000)]
    assert _arbitrate_cluster(cands, NO_EDGES) == cands


def test_weak_one_sided_dropped_when_a_two_sided_candidate_exists():
    """:271 -- 'if we found both junction for one insertion we should find both
    junction for others too'."""
    good = _cand(2, 2, 1000)
    weak = _cand(0, 1, 2000)
    assert _arbitrate_cluster([good, weak], NO_EDGES) == [good]


def test_well_supported_one_sided_survives_alongside_a_two_sided():
    """:277 -- a one-sided candidate with >= 3 reads is kept anyway."""
    good = _cand(2, 2, 1000)
    strong_one_sided = _cand(0, MIN_ONE_SIDED_JUNCTIONS, 2000)
    kept = _arbitrate_cluster([good, strong_one_sided], NO_EDGES)
    assert kept == [good, strong_one_sided]


def test_no_two_sided_keeps_only_deep_candidates_away_from_edges():
    """:289 -- needs >= 3 junction reads in total and no reference-TE boundary."""
    shallow = _cand(0, 2, 1000)
    deep = _cand(0, 3, 2000)
    assert _arbitrate_cluster([shallow, deep], NO_EDGES) == [deep]
    at_edge = {"start": {2002: 1}, "end": {}}
    assert _arbitrate_cluster([shallow, deep], at_edge) == []


def test_empty_cluster_returns_nothing():
    assert _arbitrate_cluster([], NO_EDGES) == []
