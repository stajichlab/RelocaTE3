"""Long-read clustering: Signatures into candidate insertion sites.

Signatures are built by hand; expected clusters follow from the rules in
``RelocaTE3.longread.cluster``, not from the code under test.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from RelocaTE3.longread.cluster import ClusterParams, cluster_signatures
from RelocaTE3.longread.model import (
    ClipSide, Cluster, Locus, ReadType, Signature, SignatureKind,
)  # fmt: skip
from RelocaTE3.longread.tables import read_clusters, write_clusters

REF = "ref"
W50 = ClusterParams(window=50)


def _sig(
    read_id: str,
    pos: int,
    *,
    kind: SignatureKind = SignatureKind.INS,
    length: int = 1000,
    side: ClipSide = ClipSide.NONE,
    seq_id: str = "Chr3",
    ref: str = REF,
    strand: str = "+",
) -> Signature:
    if kind is SignatureKind.CLIP and side is ClipSide.NONE:
        side = ClipSide.RIGHT
    return Signature(
        signature_id=f"{read_id}:{kind.value}@{ref}:{seq_id}:{pos}{strand}",
        kind=kind,
        locus=Locus(ref, seq_id, pos, pos, strand),
        insert_length=length,
        clip_side=side,
        read_id=read_id,
        query_start=0,
        query_end=length,
        mapq=60,
        sequence="",
        left_flank="",
        right_flank="",
    )


def _clip(read_id: str, pos: int, side: ClipSide, **kw: object) -> Signature:
    return _sig(read_id, pos, kind=SignatureKind.CLIP, length=500, side=side,
                **kw)  # type: ignore[arg-type]  # fmt: skip


def _positions(clusters: list[Cluster]) -> list[int]:
    return [c.locus.start for c in clusters]


def test_window_joins_near_and_splits_far() -> None:
    near = cluster_signatures([_sig("r1", 1000), _sig("r2", 1030)], W50)
    assert len(near) == 1
    c = near[0]
    assert c.n_reads == 2
    assert c.n_ins == 2
    assert c.locus == Locus(REF, "Chr3", 1000, 1000, "+")  # lower median
    assert c.cluster_id == "ref:Chr3:1000"
    assert c.insert_lengths == (1000, 1000)

    far = cluster_signatures([_sig("r1", 1000), _sig("r2", 1100)], W50)
    assert _positions(far) == [1000, 1100]
    assert [c.n_reads for c in far] == [1, 1]


def test_single_linkage_chains() -> None:
    sigs = [_sig("r1", 1000), _sig("r2", 1040), _sig("r3", 1080)]
    out = cluster_signatures(sigs, W50)
    assert len(out) == 1
    assert out[0].locus.start == 1040
    assert out[0].n_reads == 3


def test_window_boundary_is_inclusive() -> None:
    assert len(cluster_signatures([_sig("a", 1000), _sig("b", 1050)], W50)) == 1
    assert len(cluster_signatures([_sig("a", 1000), _sig("b", 1051)], W50)) == 2


def test_length_split() -> None:
    same = cluster_signatures(
        [_sig("r1", 1000, length=1000), _sig("r2", 1010, length=1250)], W50
    )
    assert len(same) == 1
    assert same[0].insert_lengths == (1000, 1250)

    split = cluster_signatures(
        [_sig("r1", 1000, length=1000), _sig("r2", 1010, length=1400)], W50
    )
    assert len(split) == 2
    assert sorted(c.insert_lengths for c in split) == [(1000,), (1400,)]


def test_length_split_compares_to_immediately_preceding_length() -> None:
    # 1000 -> 1250 -> 1560: each step <= 30 %, chain stays one subgroup.
    sigs = [_sig("a", 1000, length=1000), _sig("b", 1000, length=1250),
            _sig("c", 1000, length=1560)]  # fmt: skip
    out = cluster_signatures(sigs, W50)
    assert len(out) == 1
    assert out[0].insert_lengths == (1000, 1250, 1560)


def test_same_position_length_subgroups_get_suffixed_ids() -> None:
    sigs = [_sig("a", 1000, length=500), _sig("b", 1000, length=2000)]
    out = cluster_signatures(sigs, W50)
    assert [c.cluster_id for c in out] == ["ref:Chr3:1000", "ref:Chr3:1000.1"]
    # same position: ordered by first signature_id ("a:..." < "b:...")
    assert [c.insert_lengths for c in out] == [(500,), (2000,)]


def test_clips_attach_and_count_by_side() -> None:
    sigs = [
        _sig("r1", 1000),
        _sig("r2", 1005),
        _clip("r3", 1010, ClipSide.RIGHT),
        _clip("r4", 995, ClipSide.LEFT),
        _clip("r5", 5000, ClipSide.RIGHT),
    ]
    out = cluster_signatures(sigs, W50)
    assert _positions(out) == [1000, 5000]
    site, lone = out
    assert (site.n_ins, site.n_split_ins) == (2, 0)
    assert (site.n_clip_left, site.n_clip_right) == (1, 1)
    assert site.n_reads == 4
    assert site.insert_lengths == (1000, 1000)
    assert lone.cluster_id == "ref:Chr3:5000"
    assert lone.insert_lengths == ()
    assert (lone.n_ins, lone.n_split_ins) == (0, 0)
    assert (lone.n_clip_left, lone.n_clip_right) == (0, 1)
    assert lone.n_reads == 1


def test_clip_out_of_window_does_not_attach() -> None:
    out = cluster_signatures([_sig("r1", 1000), _clip("r2", 1051, ClipSide.RIGHT)],
                             W50)  # fmt: skip
    assert _positions(out) == [1000, 1051]
    assert out[1].insert_lengths == ()


def test_clip_only_clusters_single_linkage_and_lower_median() -> None:
    sigs = [_clip("a", 5000, ClipSide.LEFT), _clip("b", 5040, ClipSide.RIGHT),
            _clip("c", 5080, ClipSide.RIGHT), _clip("d", 5090, ClipSide.LEFT)]  # fmt: skip
    out = cluster_signatures(sigs, W50)
    assert len(out) == 1
    assert out[0].locus.start == 5040  # lower median of 4 starts
    assert (out[0].n_clip_left, out[0].n_clip_right) == (2, 2)


def test_clip_goes_to_nearest_subgroup_tie_earlier() -> None:
    # Groups at 1000 and 1100 (separate: 100 > 50); clip at 1050 is 50 from both.
    a, b = _sig("a", 1000), _sig("b", 1100)
    out = cluster_signatures([a, b, _clip("c", 1050, ClipSide.RIGHT)], W50)
    assert _positions(out) == [1000, 1100]
    assert out[0].n_clip_right == 1
    assert out[1].n_clip_right == 0
    # nearer wins over earlier
    out = cluster_signatures([a, b, _clip("c", 1060, ClipSide.RIGHT)], W50)
    assert (out[0].n_clip_right, out[1].n_clip_right) == (0, 1)


def test_clip_tie_same_position_goes_to_smaller_first_signature_id() -> None:
    sigs = [_sig("a", 1000, length=500), _sig("b", 1000, length=2000),
            _clip("c", 1010, ClipSide.RIGHT)]  # fmt: skip
    out = cluster_signatures(sigs, W50)
    assert [c.n_clip_right for c in out] == [1, 0]
    assert out[0].insert_lengths == (500,)


def test_one_read_with_ins_and_clip_counts_once() -> None:
    sigs = [_sig("r1", 1000), _clip("r1", 1000, ClipSide.RIGHT)]
    out = cluster_signatures(sigs, W50)
    assert len(out) == 1
    c = out[0]
    assert c.n_reads == 1
    assert (c.n_ins, c.n_clip_right) == (1, 1)
    assert len(c.signature_ids) == 2
    assert list(c.signature_ids) == sorted(c.signature_ids)


def test_sequences_and_references_never_merge() -> None:
    sigs = [
        _sig("a", 1000, seq_id="Chr3"),
        _sig("b", 1000, seq_id="Chr1"),
        _sig("c", 1000, ref="alt"),
        _clip("d", 1005, ClipSide.RIGHT, seq_id="Chr2"),
    ]
    out = cluster_signatures(sigs, W50)
    assert [c.cluster_id for c in out] == [
        "alt:Chr3:1000", "ref:Chr1:1000", "ref:Chr2:1005", "ref:Chr3:1000",
    ]  # fmt: skip
    assert all(c.n_reads == 1 for c in out)
    assert out[2].insert_lengths == ()


def _mixed() -> list[Signature]:
    sigs = [
        _sig("a", 1000, length=500),
        _sig("b", 1000, length=2000),
        _sig("c", 1008, kind=SignatureKind.SPLIT_INS, length=520),
        _clip("d", 1010, ClipSide.RIGHT),
        _clip("e", 990, ClipSide.LEFT),
        _sig("f", 3000, seq_id="Chr1"),
        _sig("g", 3030, seq_id="Chr1", length=1100),
        _clip("h", 8000, ClipSide.LEFT),
        _clip("i", 8030, ClipSide.RIGHT),
        _sig("j", 1000, ref="alt"),
        _sig("k", 1000, length=4000),
    ]
    return sigs


def test_output_stable_under_shuffles() -> None:
    sigs = _mixed()
    expected = cluster_signatures(sigs, W50)
    assert [c.cluster_id for c in expected] == [
        "alt:Chr3:1000", "ref:Chr1:3000", "ref:Chr3:1000", "ref:Chr3:1000.1",
        "ref:Chr3:1000.2", "ref:Chr3:8000",
    ]  # fmt: skip
    rng = random.Random(7)
    for _ in range(10):
        shuffled = sigs[:]
        rng.shuffle(shuffled)
        assert cluster_signatures(shuffled, W50) == expected
    assert cluster_signatures(iter(sigs), W50) == expected  # any iterable


def test_strand_is_plus_even_when_all_signatures_are_minus() -> None:
    sigs = [_sig("a", 1000, strand="-"), _clip("b", 1010, ClipSide.LEFT, strand="-"),
            _clip("c", 6000, ClipSide.RIGHT, strand="-")]  # fmt: skip
    out = cluster_signatures(sigs, W50)
    assert len(out) == 2
    assert all(c.locus.strand == "+" for c in out)


def test_tsd_offset_ins_and_split_ins_join() -> None:
    sigs = [_sig("a", 1000, length=1000),
            _sig("b", 1008, kind=SignatureKind.SPLIT_INS, length=1010)]  # fmt: skip
    out = cluster_signatures(sigs, W50)
    assert len(out) == 1
    c = out[0]
    assert (c.n_ins, c.n_split_ins) == (1, 1)
    assert c.locus.start == 1000
    assert c.insert_lengths == (1000, 1010)


def test_empty_input() -> None:
    assert cluster_signatures([], W50) == []


def test_duplicate_signature_ids_rejected() -> None:
    s = _sig("a", 1000)
    with pytest.raises(ValueError, match="duplicate signature_id"):
        cluster_signatures([s, s], W50)


def test_round_trips_through_tables(tmp_path: Path) -> None:
    expected = cluster_signatures(_mixed(), W50)
    out = tmp_path / "clusters.tsv.gz"
    n = write_clusters(out, expected, read_type=ReadType.ONT_HQ, reference_id=REF)
    assert n == len(expected)
    meta, rows = read_clusters(out)
    assert meta["reference_id"] == REF
    assert list(rows) == expected
