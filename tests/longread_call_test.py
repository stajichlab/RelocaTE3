"""Long-read calls: TE summaries per signature and read, then calls.

Hits and clusters are built by hand; expected values follow from the rules
in ``RelocaTE3.longread.call``, worked by hand, not from the code under test.
"""

from __future__ import annotations

from itertools import count

import pytest

from RelocaTE3.longread.call import (
    NESTED_MIN_BP, CallParams, SignatureFacts, SignatureTe, call_clusters,
    read_family, summarize_signature,
)  # fmt: skip
from RelocaTE3.longread.model import (
    NO_HIT, Call, CallTier, ClipSide, Cluster, FilterReason, Locus, Signature,
    SignatureKind, TeHit,
)  # fmt: skip

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


# --- call_clusters ----------------------------------------------------------

INS = SignatureKind.INS
SPLIT = SignatureKind.SPLIT_INS
CLIP = SignatureKind.CLIP
P3 = CallParams(min_support=3, min_te_identity=MIN_ID, merge_window=50)
_ids = count()


def _te(
    family: str = "mPing",
    strand: str = "+",
    *,
    bases: int = 430,
    coverage: float = 1.0,
    nested: bool = False,
    identity: float = 0.95,
) -> SignatureTe:
    if family == NO_HIT:
        return SignatureTe(NO_HIT, ".", (), coverage, False, ())
    return SignatureTe(family, strand, (((family, strand), bases),), coverage,
                       nested, (identity,))  # fmt: skip


NONE_TE = _te(NO_HIT, coverage=0.0)

Spec = tuple[str, SignatureFacts, SignatureTe]


def _s(
    read: str,
    pos: int,
    te: SignatureTe,
    kind: SignatureKind = INS,
    *,
    length: int = 430,
    side: ClipSide = ClipSide.NONE,
) -> Spec:
    if kind is CLIP and side is ClipSide.NONE:
        side = ClipSide.RIGHT
    sid = f"{read}:{kind.value}:{pos}:{next(_ids)}"
    return sid, SignatureFacts(read, kind, pos, length, side), te


class World:
    """Clusters plus the facts and TE maps that call_clusters reads."""

    def __init__(self) -> None:
        self.clusters: list[Cluster] = []
        self.facts: dict[str, SignatureFacts] = {}
        self.te: dict[str, SignatureTe] = {}

    def add(self, pos: int, sigs: list[Spec], *, seq: str = "Chr3",
            cid: str | None = None) -> str:  # fmt: skip
        facts = [f for _, f, _ in sigs]
        kinds = [f.kind for f in facts]
        sides = [f.clip_side for f in facts]
        cid = cid or f"ref:{seq}:{pos}"
        self.clusters.append(Cluster(
            cluster_id=cid,
            locus=Locus("ref", seq, pos, pos, "+"),
            n_reads=len({f.read_id for f in facts}),
            n_ins=kinds.count(INS),
            n_split_ins=kinds.count(SPLIT),
            n_clip_left=sides.count(ClipSide.LEFT),
            n_clip_right=sides.count(ClipSide.RIGHT),
            insert_lengths=tuple(sorted(f.insert_length for f in facts
                                        if f.kind is not CLIP)),
            signature_ids=tuple(sid for sid, _, _ in sigs),
        ))  # fmt: skip
        for sid, f, t in sigs:
            self.facts[sid] = f
            self.te[sid] = t
        return cid

    def call(self, params: CallParams = P3) -> list[Call]:
        return call_clusters(self.clusters, self.facts, self.te, params)


def _one(w: World, params: CallParams = P3) -> Call:
    calls = w.call(params)
    assert len(calls) == 1
    return calls[0]


def test_signature_facts_of_signature() -> None:
    sig = Signature("x", CLIP, Locus("ref", "Chr3", 1200, 1200, "+"), 600,
                    ClipSide.LEFT, "r9", 0, 600, 60, "", "", "")  # fmt: skip
    assert SignatureFacts.of(sig) == SignatureFacts("r9", CLIP, 1200, 600,
                                                    ClipSide.LEFT)  # fmt: skip


def test_three_concordant_reads_pass_spanning() -> None:
    w = World()
    cid = w.add(1000, [_s("r1", 1003, _te(), length=432),
                       _s("r2", 998, _te(), length=428),
                       _s("r3", 1000, _te(), length=430)])  # fmt: skip
    c = _one(w)
    assert c.passed
    assert c.filters == ()
    assert c.tier is CallTier.SPANNING
    assert c.call_id == "ref:Chr3:1000"
    assert c.locus == Locus("ref", "Chr3", 1000, 1000, "+")
    assert (c.family, c.te_strand) == ("mPing", "+")
    assert c.family_fraction == 1.0
    assert c.strand_fraction == 1.0
    assert c.te_read_fraction == 1.0
    assert (c.n_reads, c.n_ins, c.n_split_ins, c.n_clip_left, c.n_clip_right) == (
        3, 3, 0, 0, 0)  # fmt: skip
    assert c.median_insert_length == 430
    assert c.median_te_coverage == 1.0
    assert c.median_identity == 0.95
    assert c.nested is False
    assert c.cluster_ids == (cid,)


def test_position_is_lower_median_of_insert_positions() -> None:
    w = World()
    w.add(1000, [_s("r1", 1000, _te()), _s("r2", 1010, _te(), SPLIT),
                 _s("r3", 1020, _te()), _s("r4", 1030, _te()),
                 _s("r5", 900, _te(), CLIP)])  # fmt: skip
    c = _one(w)
    # inserts 1000, 1010, 1020, 1030 -> lower median 1010; the clip is ignored
    assert c.locus.start == 1010
    assert (c.n_ins, c.n_split_ins, c.n_clip_right) == (3, 1, 1)


def test_four_of_five_reads_with_hit_pass() -> None:
    w = World()
    w.add(1000, [_s(f"r{i}", 1000, _te()) for i in range(4)]
          + [_s("r4", 1000, NONE_TE)])  # fmt: skip
    c = _one(w)
    assert c.te_read_fraction == pytest.approx(0.8)
    assert c.family_fraction == 1.0
    assert c.passed


def test_three_of_five_family_discordant() -> None:
    w = World()
    w.add(1000, [_s(f"m{i}", 1000, _te()) for i in range(3)]
          + [_s(f"g{i}", 1000, _te("Gaijin")) for i in range(2)])  # fmt: skip
    c = _one(w)
    assert c.family == "mPing"
    assert c.family_fraction == pytest.approx(0.6)
    assert c.te_read_fraction == 1.0
    assert c.filters == (FilterReason.FAMILY_DISCORDANT,)


def test_three_hit_two_none_few_te_reads() -> None:
    w = World()
    w.add(1000, [_s(f"m{i}", 1000, _te()) for i in range(3)]
          + [_s(f"n{i}", 1000, NONE_TE) for i in range(2)])  # fmt: skip
    c = _one(w)
    assert c.te_read_fraction == pytest.approx(0.6)
    assert c.filters == (FilterReason.FEW_TE_READS,)


def test_low_support() -> None:
    w = World()
    w.add(1000, [_s("r1", 1000, _te())])
    two = CallParams(min_support=2, min_te_identity=MIN_ID, merge_window=50)
    one = CallParams(min_support=1, min_te_identity=MIN_ID, merge_window=50)
    assert _one(w, two).filters == (FilterReason.LOW_SUPPORT,)
    assert _one(w, one).passed


def test_low_te_fraction_only_for_spanning() -> None:
    w = World()
    w.add(1000, [_s(f"r{i}", 1000, _te(coverage=0.5)) for i in range(3)])
    c = _one(w)
    assert c.median_te_coverage == 0.5
    assert c.filters == (FilterReason.LOW_TE_FRACTION,)

    clip = World()
    clip.add(1000, [_s(f"r{i}", 1000, _te(coverage=0.5), CLIP) for i in range(3)])
    c = _one(clip)
    assert c.tier is CallTier.CLIP_ONLY
    assert c.median_te_coverage == 0.5
    assert c.passed


def test_median_te_coverage_is_statistics_median_rounded() -> None:
    w = World()
    covs = [0.9, 0.85, 0.83333333, 0.95]
    w.add(1000, [_s(f"r{i}", 1000, _te(coverage=v)) for i, v in enumerate(covs)])
    # sorted 0.83333333, 0.85, 0.9, 0.95 -> (0.85 + 0.9) / 2
    assert _one(w).median_te_coverage == 0.875


def test_no_hits_clip_only_filters_in_order() -> None:
    w = World()
    w.add(1000, [_s("r1", 1000, NONE_TE, CLIP)])
    two = CallParams(min_support=2, min_te_identity=MIN_ID, merge_window=50)
    c = _one(w, two)
    assert c.filters == (FilterReason.NO_TE_HIT, FilterReason.LOW_SUPPORT)
    assert (c.family, c.te_strand) == (NO_HIT, ".")
    assert (c.family_fraction, c.strand_fraction, c.te_read_fraction) == (
        0.0, 0.0, 0.0)  # fmt: skip
    assert c.median_identity == 0.0
    assert c.nested is False


def test_no_hits_spanning_also_low_te_fraction() -> None:
    w = World()
    w.add(1000, [_s("r1", 1000, NONE_TE)])
    two = CallParams(min_support=2, min_te_identity=MIN_ID, merge_window=50)
    assert _one(w, two).filters == (FilterReason.NO_TE_HIT,
                                    FilterReason.LOW_SUPPORT,
                                    FilterReason.LOW_TE_FRACTION)  # fmt: skip


def test_clip_only_tier_and_position() -> None:
    w = World()
    w.add(1000, [_s("r1", 1010, _te(), CLIP, side=ClipSide.LEFT),
                 _s("r2", 1000, _te(), CLIP), _s("r3", 1020, _te(), CLIP),
                 _s("r4", 1030, _te(), CLIP)])  # fmt: skip
    c = _one(w)
    assert c.tier is CallTier.CLIP_ONLY
    assert c.median_insert_length == 0
    assert c.locus.start == 1010  # lower median of 1000, 1010, 1020, 1030
    assert (c.n_ins, c.n_clip_left, c.n_clip_right) == (0, 1, 3)
    assert c.passed


def test_merge_within_window() -> None:
    w = World()
    a = w.add(1000, [_s("r1", 1000, _te()), _s("r2", 1002, _te())])
    b = w.add(1030, [_s("r2", 1030, _te()), _s("r3", 1032, _te())])
    c = _one(w)
    assert c.cluster_ids == (a, b)
    assert c.n_reads == 3
    assert c.n_ins == 4
    assert c.locus.start == 1002  # lower median of 1000, 1002, 1030, 1032
    assert c.call_id == "ref:Chr3:1002"


def test_no_merge_beyond_window() -> None:
    w = World()
    w.add(1000, [_s(f"a{i}", 1000, _te()) for i in range(3)])
    w.add(1060, [_s(f"b{i}", 1060, _te()) for i in range(3)])
    calls = w.call()
    assert [c.call_id for c in calls] == ["ref:Chr3:1000", "ref:Chr3:1060"]
    assert [c.cluster_ids for c in calls] == [("ref:Chr3:1000",), ("ref:Chr3:1060",)]


def test_merge_chains_on_last_cluster_start() -> None:
    w = World()
    w.add(1000, [_s(f"a{i}", 1000, _te()) for i in range(3)])
    w.add(1040, [_s(f"b{i}", 1040, _te()) for i in range(3)])
    w.add(1080, [_s(f"c{i}", 1080, _te()) for i in range(3)])
    c = _one(w)
    assert c.cluster_ids == ("ref:Chr3:1000", "ref:Chr3:1040", "ref:Chr3:1080")
    assert c.n_reads == 9


def test_same_family_merges_around_other_family() -> None:
    w = World()
    a = w.add(1000, [_s(f"a{i}", 1000, _te()) for i in range(2)])
    g = w.add(1010, [_s(f"g{i}", 1010, _te("Gaijin")) for i in range(3)])
    b = w.add(1030, [_s(f"b{i}", 1030, _te()) for i in range(2)])
    calls = w.call()
    assert [(c.family, c.cluster_ids) for c in calls] == [
        ("mPing", (a, b)), ("Gaijin", (g,))]  # fmt: skip
    assert calls[0].locus.start == 1000  # lower median of 1000 x2, 1030 x2
    assert calls[0].n_reads == 4


def test_different_sequences_never_merge() -> None:
    w = World()
    w.add(1000, [_s(f"c{i}", 1000, _te()) for i in range(3)], seq="Chr3")
    w.add(1000, [_s(f"a{i}", 1000, _te()) for i in range(3)], seq="Chr1")
    assert [c.call_id for c in w.call()] == ["ref:Chr1:1000", "ref:Chr3:1000"]


def test_no_hit_clusters_never_merge() -> None:
    w = World()
    w.add(1000, [_s("r1", 1000, NONE_TE)])
    w.add(1010, [_s("r2", 1010, NONE_TE)])
    calls = w.call()
    assert [c.call_id for c in calls] == ["ref:Chr3:1000", "ref:Chr3:1010"]
    assert all(c.family == NO_HIT for c in calls)


def test_cluster_family_by_read_count_then_bases() -> None:
    # Cluster at 1020: 2 mPing reads (200 bases each) and 2 Gaijin reads
    # (100 each): read counts tie, mPing has more bases, so the cluster is
    # mPing and merges with the mPing cluster at 1000.
    w = World()
    a = w.add(1000, [_s(f"a{i}", 1000, _te()) for i in range(3)])
    b = w.add(1020, [_s("m1", 1020, _te(bases=200)), _s("m2", 1020, _te(bases=200)),
                     _s("g1", 1020, _te("Gaijin", bases=100)),
                     _s("g2", 1020, _te("Gaijin", bases=100))])  # fmt: skip
    c = _one(w)
    assert c.cluster_ids == (a, b)
    assert c.family == "mPing"
    assert c.family_fraction == pytest.approx(5 / 7)


def test_vote_tie_on_reads_goes_to_more_bases() -> None:
    w = World()
    w.add(1000, [_s("m1", 1000, _te(bases=200)), _s("m2", 1000, _te(bases=200)),
                 _s("g1", 1000, _te("Gaijin", bases=100)),
                 _s("g2", 1000, _te("Gaijin", bases=100))])  # fmt: skip
    c = _one(w)
    assert c.family == "mPing"
    assert c.family_fraction == 0.5
    assert c.filters == (FilterReason.FAMILY_DISCORDANT,)


def test_clip_only_cluster_merges_into_spanning() -> None:
    w = World()
    a = w.add(1000, [_s(f"r{i}", 1000, _te()) for i in range(3)])
    b = w.add(1040, [
        _s("c1", 1040, _te(coverage=0.3, identity=0.80), CLIP, side=ClipSide.LEFT),
        _s("c2", 1041, _te(coverage=0.3, identity=0.80), CLIP),
    ])  # fmt: skip
    c = _one(w)
    assert c.cluster_ids == (a, b)
    assert c.tier is CallTier.SPANNING
    assert (c.n_ins, c.n_clip_left, c.n_clip_right, c.n_reads) == (3, 1, 1, 5)
    assert c.locus.start == 1000
    assert c.median_te_coverage == 1.0  # insert signatures only
    assert c.median_identity == 0.95  # 0.80, 0.80, 0.95, 0.95, 0.95
    assert c.passed


def test_strand_tie_goes_to_plus() -> None:
    w = World()
    w.add(1000, [_s("r1", 1000, _te(strand="-")), _s("r2", 1000, _te(strand="+"))])
    two = CallParams(min_support=2, min_te_identity=MIN_ID, merge_window=50)
    c = _one(w, two)
    assert c.te_strand == "+"
    assert c.strand_fraction == 0.5
    assert c.passed


def test_strand_majority_among_top_family_reads() -> None:
    w = World()
    w.add(1000, [_s("r1", 1000, _te(strand="-")), _s("r2", 1000, _te(strand="-")),
                 _s("r3", 1000, _te(strand="+")),
                 _s("g1", 1000, _te("Gaijin", "+"))])  # fmt: skip
    c = _one(w)
    assert (c.family, c.te_strand) == ("mPing", "-")
    assert c.strand_fraction == pytest.approx(2 / 3)
    assert c.family_fraction == pytest.approx(3 / 4)


def test_read_family_uses_only_group_signatures() -> None:
    # r1 has a Gaijin signature in a far cluster; inside the mPing call it
    # is an mPing read.
    w = World()
    w.add(1000, [_s("r1", 1000, _te(bases=100)), _s("r2", 1000, _te()),
                 _s("r3", 1000, _te())])  # fmt: skip
    w.add(
        5000,
        [
            _s("r1", 5000, _te("Gaijin", bases=1000)),
            _s("x1", 5000, _te("Gaijin")),
            _s("x2", 5000, _te("Gaijin")),
        ],
    )
    calls = w.call()
    assert [(c.family, c.family_fraction) for c in calls] == [
        ("mPing", 1.0), ("Gaijin", 1.0)]  # fmt: skip


def test_median_identity_of_top_family_signatures() -> None:
    w = World()
    w.add(1000, [_s("r1", 1000, _te(identity=0.9)), _s("r2", 1000, _te(identity=0.8)),
                 _s("r3", 1000, _te(identity=0.85)),
                 _s("g1", 1000, _te("Gaijin", identity=0.1))])  # fmt: skip
    assert _one(w).median_identity == 0.85


def test_nested_when_half_of_hit_reads_nested() -> None:
    w = World()
    w.add(1000, [_s("r1", 1000, _te(nested=True)), _s("r2", 1000, _te()),
                 _s("r3", 1000, _te(nested=True)), _s("r4", 1000, _te()),
                 _s("n1", 1000, NONE_TE)])  # fmt: skip
    assert _one(w).nested is True
    w = World()
    w.add(1000, [_s("r1", 1000, _te(nested=True)), _s("r2", 1000, _te()),
                 _s("r3", 1000, _te())])  # fmt: skip
    assert _one(w).nested is False


def test_duplicate_positions_get_suffix() -> None:
    w = World()
    w.add(1000, [_s(f"m{i}", 1000, _te()) for i in range(3)], cid="ref:Chr3:1000.1")
    w.add(1000, [_s(f"g{i}", 1000, _te("Gaijin")) for i in range(3)],
          cid="ref:Chr3:1000")  # fmt: skip
    calls = w.call()
    assert [(c.call_id, c.family) for c in calls] == [
        ("ref:Chr3:1000", "Gaijin"), ("ref:Chr3:1000.1", "mPing")]  # fmt: skip


def test_input_order_does_not_matter() -> None:
    w = World()
    w.add(1030, [_s(f"b{i}", 1030, _te()) for i in range(2)])
    w.add(1010, [_s(f"g{i}", 1010, _te("Gaijin")) for i in range(3)])
    w.add(1000, [_s(f"a{i}", 1000, _te()) for i in range(2)])
    calls = w.call()
    assert [c.cluster_ids for c in calls] == [
        ("ref:Chr3:1000", "ref:Chr3:1030"), ("ref:Chr3:1010",)]  # fmt: skip


@pytest.mark.parametrize("drop", ["facts", "te"])
def test_missing_signature_raises(drop: str) -> None:
    w = World()
    w.add(1000, [_s("r1", 1000, _te())])
    sid = w.clusters[0].signature_ids[0]
    del getattr(w, drop)[sid]
    with pytest.raises(ValueError, match=sid):
        w.call()
