"""Long-read data model: loci, references, placements, signatures, clusters."""

from __future__ import annotations

import dataclasses

import pytest

from RelocaTE3.longread.model import (
    NO_HIT, Call, CallTier, ClipSide, Cluster, FilterReason, Locus,
    PlacementRecord, ReadType, Reference, Signature, SignatureKind, TeHit,
)  # fmt: skip


def test_read_type_presets_cover_the_four_types():
    assert {t.value: t.minimap2_preset for t in ReadType} == {
        "ont": "map-ont", "ont-hq": "lr:hq", "hifi": "map-hifi", "clr": "map-pb",
    }  # fmt: skip


@pytest.mark.parametrize(
    ("read_type", "window"),
    [("hifi", 50), ("ont-hq", 50), ("ont", 100), ("clr", 100)],
)
def test_default_cluster_window_by_read_type(read_type, window):
    assert ReadType(read_type).default_cluster_window == window


def test_locus_is_frozen_and_validates():
    locus = Locus("ref", "Chr3", 10, 10, "+")
    with pytest.raises(dataclasses.FrozenInstanceError):
        locus.start = 5  # type: ignore[misc]
    with pytest.raises(ValueError, match="end"):
        Locus("ref", "Chr3", 10, 9, "+")
    with pytest.raises(ValueError, match="strand"):
        Locus("ref", "Chr3", 1, 2, "x")


def test_locus_key_includes_the_reference():
    assert Locus("a", "Chr3", 1, 1, "+").key != Locus("b", "Chr3", 1, 1, "+").key


def test_signature_kind_and_clip_side_values():
    assert {k.value for k in SignatureKind} == {"INS", "SPLIT_INS", "CLIP"}
    assert {s.value for s in ClipSide} == {"left", "right", ""}


def test_records_build():
    locus = Locus("ref", "Chr3", 100, 100, "+")
    rec = PlacementRecord(
        read_id="r1", read_length=1000, query_start=0, query_end=1000,
        locus=Locus("ref", "Chr3", 0, 900, "+"), is_supplementary=False,
        mapq=60, insertions=((100, 300, 400),), left_clip=0, right_clip=0,
        sequence="A" * 1000, has_supplementary=False,
    )  # fmt: skip
    sig = Signature(
        signature_id="s1", kind=SignatureKind.INS, locus=locus, insert_length=100,
        clip_side=ClipSide.NONE, read_id="r1", query_start=300, query_end=400,
        mapq=60, sequence="A" * 100, left_flank="", right_flank="",
    )  # fmt: skip
    cl = Cluster(
        cluster_id="ref:Chr3:100", locus=locus, n_reads=1, n_ins=1, n_split_ins=0,
        n_clip_left=0, n_clip_right=0, insert_lengths=(100,), signature_ids=("s1",),
    )  # fmt: skip
    assert rec.locus.sequence_id == sig.locus.sequence_id == cl.locus.sequence_id
    assert Reference("ref", "g.fa").index is None


@pytest.mark.parametrize(
    "member",
    [*ReadType, *SignatureKind, *ClipSide],
    ids=lambda m: f"{type(m).__name__}.{m.name}",
)
def test_enum_text_form_is_the_value(member):
    assert str(member) == f"{member}" == member.value


def test_clip_side_none_text_is_empty():
    assert str(ClipSide.NONE) == f"{ClipSide.NONE}" == ""


@pytest.mark.parametrize(
    "read_type,support,identity",
    [(ReadType.HIFI, 2, 0.85), (ReadType.ONT_HQ, 2, 0.80),
     (ReadType.ONT, 3, 0.70), (ReadType.CLR, 3, 0.70)],
)  # fmt: skip
def test_read_type_call_defaults(read_type, support, identity) -> None:
    assert read_type.default_min_support == support
    assert read_type.default_min_te_identity == identity


def test_tier_and_reason_text() -> None:
    assert f"{CallTier.CLIP_ONLY}" == "clip_only"
    assert str(CallTier.SPANNING) == "spanning"
    assert [r.value for r in FilterReason] == [
        "no_te_hit", "few_te_reads", "family_discordant", "low_support",
        "low_te_fraction",
    ]  # fmt: skip


def _hit(**kw) -> TeHit:
    base = dict(signature_id="s1", read_id="r1", family="mPing", te_strand="+",
                te_start=0, te_end=430, te_length=430, query_start=5,
                query_end=435, query_length=440, identity=0.99)  # fmt: skip
    base.update(kw)
    return TeHit(**base)


def test_te_hit_checks() -> None:
    assert _hit().aligned == 430
    _hit(family=NO_HIT, te_strand=".", te_start=0, te_end=0, te_length=0,
         query_start=0, query_end=0, identity=0.0)  # fmt: skip
    with pytest.raises(ValueError, match="strand"):
        _hit(te_strand=".")  # a real hit needs + or -
    with pytest.raises(ValueError, match="strand"):
        _hit(family=NO_HIT, te_strand="+")
    with pytest.raises(ValueError, match="query"):
        _hit(query_start=10, query_end=5)


def test_call_passed() -> None:
    loc = Locus("ref", "Chr3", 100, 100, "+")
    common = dict(call_id="c", locus=loc, tier=CallTier.SPANNING,
                  family="mPing", te_strand="+", family_fraction=1.0,
                  strand_fraction=1.0, te_read_fraction=1.0, n_reads=3,
                  n_ins=3, n_split_ins=0, n_clip_left=0, n_clip_right=0,
                  median_insert_length=435, median_te_coverage=0.99,
                  median_identity=0.99, nested=False, cluster_ids=("k",))  # fmt: skip
    assert Call(filters=(), **common).passed
    assert not Call(filters=(FilterReason.LOW_SUPPORT,), **common).passed
    with pytest.raises(ValueError, match="insertion point"):
        Call(filters=(), **{**common, "locus": Locus("ref", "Chr3", 1, 5, "+")})
