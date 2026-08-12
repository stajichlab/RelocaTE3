"""Adjacent insertions must not be merged into one impossible call.

``_pair_breakpoints`` groups breakpoints within ``SUBCLUSTER_GAP`` and then keeps
the dominant left and dominant right of each group. When two real insertions sit
close together, their *inner* breakpoints fall inside that gap, so the group
holds both -- and the dominant pair ends up being insertion A's left breakpoint
with insertion B's right breakpoint. Both true sites are then lost: one call is
emitted, at neither position.

The pairing is geometrically impossible, which is what makes it detectable
without any extra evidence. A left-junction read marks the *right* edge of the
TSD and a right-junction read marks the *left* edge (see ``_make_insertion``), so
a genuine pair always has ``left_position >= right_position``. When the left
breakpoint lies to the left of the right breakpoint, the two breakpoints cannot
belong to the same insertion.

Observed on the Chr3 2 Mb fixture at two of the three sites RelocaTE2 finds and
RelocaTE3 missed entirely:

    truth Chr3:155988   RelocaTE2: 155951..155962 and 155984..155988
                        RelocaTE3: one call at 155962 (L=155962, R=155984 -> -21)
    truth Chr3:672724   RelocaTE2: 672693..672695 and 672720..672724
                        RelocaTE3: one call at 672695 (L=672695, R=672720 -> -24)
"""

from __future__ import annotations

from RelocaTE3.insertions import _pair_breakpoints


def _obs(n: int) -> list[object]:
    """n placeholder junction observations -- only the count is read here."""
    return [object()] * n


def test_impossible_pairing_is_split_into_two_sub_insertions():
    """L to the left of R cannot be one TSD; emit both breakpoints separately."""
    left = {155962: _obs(3)}
    right = {155984: _obs(5)}
    pairs = _pair_breakpoints(left, right)
    assert (155962, 155984) not in pairs, "must not pair across two insertions"
    assert sorted(pairs, key=lambda p: (p[0] or 0, p[1] or 0)) == [
        (None, 155984),
        (155962, None),
    ]


def test_a_real_tsd_still_pairs():
    """left >= right is the valid arrangement and must be preserved."""
    left = {1005: _obs(4)}
    right = {1003: _obs(4)}
    assert _pair_breakpoints(left, right) == [(1005, 1003)]


def test_coincident_breakpoints_pair():
    """A zero-length TSD (left == right) is still one insertion."""
    left = {500: _obs(2)}
    right = {500: _obs(2)}
    assert _pair_breakpoints(left, right) == [(500, 500)]


def test_distant_breakpoints_remain_separate_groups():
    """Beyond SUBCLUSTER_GAP the grouping already separated them."""
    left = {1000: _obs(2)}
    right = {2000: _obs(2)}
    pairs = _pair_breakpoints(left, right)
    assert sorted(pairs, key=lambda p: (p[0] or 0, p[1] or 0)) == [
        (None, 2000),
        (1000, None),
    ]


def test_one_sided_groups_are_unchanged():
    assert _pair_breakpoints({700: _obs(3)}, {}) == [(700, None)]
    assert _pair_breakpoints({}, {700: _obs(3)}) == [(None, 700)]


def test_dominant_selection_still_applies_within_a_valid_group():
    """Among several left breakpoints the best-supported one wins."""
    left = {1005: _obs(1), 1006: _obs(7)}
    right = {1003: _obs(4)}
    pairs = _pair_breakpoints(left, right)
    assert pairs == [(1006, 1003)]


def test_splitting_does_not_invent_extra_calls():
    """A split yields exactly the two breakpoints, not a third merged one."""
    pairs = _pair_breakpoints({100: _obs(2)}, {120: _obs(2)})
    assert len(pairs) == 2
    assert all((p[0] is None) != (p[1] is None) for p in pairs)


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
