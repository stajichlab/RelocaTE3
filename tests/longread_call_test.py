"""Long-read calls: TE summaries per signature and read, then calls.

Hits and clusters are built by hand; expected values follow from the rules
in ``RelocaTE3.longread.call``, worked by hand, not from the code under test.
"""

from __future__ import annotations

import pytest

from RelocaTE3.longread.call import (
    NESTED_MIN_BP, SignatureTe, read_family, summarize_signature,
)  # fmt: skip
from RelocaTE3.longread.model import NO_HIT, TeHit

MIN_ID = 0.70


def _hit(
    family: str,
    qs: int,
    qe: int,
    qlen: int,
    *,
    strand: str = "+",
    identity: float = 0.95,
    sig: str = "s1",
) -> TeHit:
    te_len = max(qe - qs, 1)
    return TeHit(sig, "r1", family, strand, 0, qe - qs, te_len, qs, qe, qlen,
                 identity)  # fmt: skip


def _none(qlen: int = 440, sig: str = "s1") -> TeHit:
    return TeHit(sig, "r1", NO_HIT, ".", 0, 0, 0, 0, 0, qlen, 0.0)


# --- summarize_signature ----------------------------------------------------


def test_no_hit_row_gives_no_family() -> None:
    s = summarize_signature([_none()], MIN_ID)
    assert s.family == NO_HIT
    assert s.te_strand == "."
    assert s.bases == ()
    assert s.coverage == 0.0
    assert s.nested is False
    assert s.identities == ()


def test_one_hit_family_strand_coverage() -> None:
    s = summarize_signature([_hit("mPing", 5, 435, 440, identity=0.9)], MIN_ID)
    assert (s.family, s.te_strand) == ("mPing", "+")
    assert s.bases == ((("mPing", "+"), 430),)
    assert s.coverage == pytest.approx(430 / 440)
    assert s.nested is False
    assert s.identities == (0.9,)


def test_hit_below_identity_ignored() -> None:
    s = summarize_signature([_hit("mPing", 0, 440, 440, identity=0.69)], MIN_ID)
    assert s.family == NO_HIT
    assert s.coverage == 0.0


def test_low_identity_hit_does_not_count_next_to_good_one() -> None:
    hits = [_hit("Gaijin", 0, 400, 440, identity=0.5),
            _hit("mPing", 0, 200, 440, identity=0.9)]  # fmt: skip
    s = summarize_signature(hits, MIN_ID)
    assert s.family == "mPing"
    assert s.bases == ((("mPing", "+"), 200),)
    assert s.coverage == pytest.approx(200 / 440)


def test_larger_aligned_total_wins() -> None:
    hits = [_hit("mPing", 0, 100, 500), _hit("mPing", 100, 200, 500),
            _hit("Gaijin", 200, 350, 500)]  # fmt: skip
    s = summarize_signature(hits, MIN_ID)
    assert s.family == "mPing"  # 200 vs 150
    assert s.bases == ((("mPing", "+"), 200), (("Gaijin", "+"), 150))


def test_equal_totals_family_name_order() -> None:
    hits = [_hit("mPing", 0, 200, 400), _hit("Gaijin", 200, 400, 400)]
    assert summarize_signature(hits, MIN_ID).family == "Gaijin"


def test_equal_totals_plus_strand_first() -> None:
    hits = [_hit("mPing", 0, 200, 400, strand="-"),
            _hit("mPing", 200, 400, 400, strand="+")]  # fmt: skip
    s = summarize_signature(hits, MIN_ID)
    assert (s.family, s.te_strand) == ("mPing", "+")


def test_coverage_is_union_of_kept_hits() -> None:
    hits = [_hit("mPing", 0, 300, 500), _hit("mPing", 200, 500, 500)]
    assert summarize_signature(hits, MIN_ID).coverage == pytest.approx(1.0)


def test_nested_other_family_in_gap_of_best() -> None:
    hits = [_hit("Retro1", 0, 500, 1509), _hit("mPing", 500, 930, 1509),
            _hit("Retro1", 930, 1509, 1509)]  # fmt: skip
    s = summarize_signature(hits, MIN_ID)
    assert s.family == "Retro1"  # 1079 vs 430
    assert s.nested is True
    assert s.coverage == pytest.approx(1.0)
    assert s.identities == (0.95, 0.95)


def test_other_family_inside_best_span_not_nested() -> None:
    # mPing [400,800) overlaps Retro1 [0,1000) by 400 bp > 10% of 400.
    hits = [_hit("Retro1", 0, 1000, 1000), _hit("mPing", 400, 800, 1000)]
    assert summarize_signature(hits, MIN_ID).nested is False


def test_other_family_overlap_at_ten_percent_still_nested() -> None:
    # mPing [450,950) overlaps Retro1 [0,500) by 50 = 10% of 500: allowed.
    hits = [_hit("Retro1", 0, 500, 1600), _hit("mPing", 450, 950, 1600),
            _hit("Retro1", 950, 1600, 1600)]  # fmt: skip
    assert summarize_signature(hits, MIN_ID).nested is True


def test_short_other_family_hit_never_nested() -> None:
    short = NESTED_MIN_BP - 10
    hits = [_hit("Retro1", 0, 500, 1000), _hit("mPing", 500, 500 + short, 1000),
            _hit("Retro1", 500 + short, 1000, 1000)]  # fmt: skip
    assert summarize_signature(hits, MIN_ID).nested is False


# --- read_family ------------------------------------------------------------


def _st(*bases: tuple[str, str, int]) -> SignatureTe:
    if not bases:
        return SignatureTe(NO_HIT, ".", (), 0.0, False, ())
    items = tuple(((f, s), n) for f, s, n in bases)
    return SignatureTe(bases[0][0], bases[0][1], items, 1.0, False, ())


def test_read_family_sums_over_signatures() -> None:
    sums = [_st(("Gaijin", "+", 300)), _st(("mPing", "+", 200)),
            _st(("mPing", "+", 200), ("Gaijin", "+", 50))]  # fmt: skip
    assert read_family(sums) == ("mPing", "+")  # 400 vs 350


def test_read_family_ties() -> None:
    assert read_family([_st(("mPing", "+", 100)),
                        _st(("Gaijin", "+", 100))]) == ("Gaijin", "+")  # fmt: skip
    assert read_family([_st(("mPing", "-", 100)),
                        _st(("mPing", "+", 100))]) == ("mPing", "+")  # fmt: skip


def test_read_family_none_without_kept_hits() -> None:
    assert read_family([_st(), _st()]) is None
    assert read_family([]) is None
