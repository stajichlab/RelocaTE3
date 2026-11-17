"""Long-read signatures: INS, SPLIT_INS and CLIP from PlacementRecords.

Records are built by hand (no BAM). Each read is assembled from known
pieces, so expected sequences, flanks and coordinates come from the pieces,
not from the code under test.
"""

from __future__ import annotations

import random
from dataclasses import replace
from pathlib import Path

import pytest

from RelocaTE3.longread.model import (
    ClipSide, Locus, PlacementRecord, ReadType, Signature, SignatureKind,
)  # fmt: skip
from RelocaTE3.longread.signatures import (
    SignatureParams, clip_signatures, collect_signatures, ins_signatures,
    split_ins_signatures,
)  # fmt: skip
from RelocaTE3.longread.tables import read_signatures, write_signatures

REF = "msu7"
P = SignatureParams()  # min_ins 50, merge 20, min_clip 100, ref gap 50, flank 200
_COMP = str.maketrans("ACGT", "TGCA")


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _seq(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _rec(
    read_id: str,
    *,
    qs: int,
    qe: int,
    start: int,
    end: int,
    read_length: int,
    sequence: str | None,
    strand: str = "+",
    seq_id: str = "Chr3",
    insertions: tuple[tuple[int, int, int], ...] = (),
    lc: int = 0,
    rc: int = 0,
    lh: int = 0,
    rh: int = 0,
    mapq: int = 60,
    supp: bool = False,
    has_supp: bool = False,
) -> PlacementRecord:
    return PlacementRecord(
        read_id=read_id,
        read_length=read_length,
        query_start=qs,
        query_end=qe,
        locus=Locus(REF, seq_id, start, end, strand),
        is_supplementary=supp,
        mapq=mapq,
        insertions=insertions,
        left_clip=lc,
        right_clip=rc,
        sequence=sequence,
        has_supplementary=has_supp,
        left_hard_clip=lh,
        right_hard_clip=rh,
    )


def _ids(sigs: list[Signature]) -> list[str]:
    return [s.signature_id for s in sigs]


# --- INS ---------------------------------------------------------------------


def test_forward_spanning_insertion() -> None:
    left, insert, right = _seq(300, 1), _seq(200, 2), _seq(300, 3)
    read = left + insert + right
    # 300M200I300M at reference 1000.
    rec = _rec("r1", qs=0, qe=800, start=1000, end=1600, read_length=800,
               sequence=read, insertions=((1300, 300, 500),), mapq=42)  # fmt: skip
    sigs = ins_signatures(rec, P)
    assert len(sigs) == 1
    s = sigs[0]
    assert s.kind is SignatureKind.INS
    assert s.locus == Locus(REF, "Chr3", 1300, 1300, "+")
    assert s.insert_length == 200
    assert s.sequence == insert
    assert s.left_flank == left[-200:]
    assert s.right_flank == right[:200]
    assert (s.query_start, s.query_end) == (300, 500)
    assert read[s.query_start : s.query_end] == insert
    assert s.mapq == 42
    assert s.clip_side is ClipSide.NONE
    assert s.read_id == "r1"
    assert s.signature_id == "r1:INS:300-500@msu7:Chr3:1300+"
    assert clip_signatures(rec, P) == []


def test_reverse_spanning_insertion_with_hard_clip() -> None:
    # Reference-forward picture of the read; the read as sequenced is its
    # reverse complement. The first 30 stored bases are hard-clipped.
    left, insert, right = _seq(300, 4), _seq(200, 5), _seq(300, 6)
    refwd = left + insert + right
    read = _rc(refwd)
    stored = refwd[30:]
    # 30H270M200I300M at reference 1030.
    rec = _rec("r2", qs=0, qe=770, start=1030, end=1600, read_length=800,
               sequence=stored, strand="-", insertions=((1300, 270, 470),),
               lh=30)  # fmt: skip
    (s,) = ins_signatures(rec, P)
    assert s.locus == Locus(REF, "Chr3", 1300, 1300, "-")
    assert s.insert_length == 200
    assert s.sequence == insert  # reference-forward slice
    assert s.left_flank == left[-200:]
    assert s.right_flank == right[:200]
    # Reverse conversion: (L - lh - b, L - lh - a) = (800-30-470, 800-30-270).
    assert (s.query_start, s.query_end) == (300, 500)
    assert _rc(s.sequence) == read[s.query_start : s.query_end]


def test_ins_flanks_truncate_at_read_ends() -> None:
    left, insert, right = _seq(50, 7), _seq(100, 8), _seq(20, 9)
    rec = _rec("r", qs=0, qe=170, start=0, end=70, read_length=170,
               sequence=left + insert + right,
               insertions=((50, 50, 150),))  # fmt: skip
    (s,) = ins_signatures(rec, P)
    assert (s.left_flank, s.right_flank) == (left, right)


def test_ont_pieces_merge_into_one_insertion() -> None:
    a, p1, m, p2, b = (
        _seq(100, 10),
        _seq(30, 11),
        _seq(8, 12),
        _seq(40, 13),
        _seq(100, 14),
    )
    read = a + p1 + m + p2 + b
    # 100M30I8M40I100M at reference 1000.
    rec = _rec("ont", qs=0, qe=len(read), start=1000, end=1208,
               read_length=len(read), sequence=read,
               insertions=((1108, 138, 178), (1100, 100, 130)))  # fmt: skip
    (s,) = ins_signatures(rec, P)
    assert s.locus.start == s.locus.end == 1100
    assert s.insert_length == 30 + 8 + 40
    assert s.sequence == p1 + m + p2
    assert s.left_flank == a
    assert s.right_flank == b
    assert (s.query_start, s.query_end) == (100, 178)


def test_ont_pieces_far_apart_and_short_give_nothing() -> None:
    a, p1, m, p2, b = (
        _seq(100, 15),
        _seq(30, 16),
        _seq(60, 17),
        _seq(40, 18),
        _seq(100, 19),
    )
    read = a + p1 + m + p2 + b
    # 100M30I60M40I100M: pieces 60 bp apart on the reference, each < 50.
    rec = _rec("ont", qs=0, qe=len(read), start=1000, end=1260,
               read_length=len(read), sequence=read,
               insertions=((1100, 100, 130), (1160, 190, 230)))  # fmt: skip
    assert ins_signatures(rec, P) == []


@pytest.mark.parametrize(("gap", "n_sigs"), [(20, 1), (21, 0)])
def test_merge_distance_is_inclusive(gap: int, n_sigs: int) -> None:
    # Two 30 bp pieces: merged they reach 60 + gap >= 50, alone 30 < 50.
    rec = _rec("m", qs=0, qe=260 + gap, start=1000, end=1200 + gap,
               read_length=260 + gap, sequence=_seq(260 + gap, 20),
               insertions=((1100, 100, 130), (1100 + gap, 130 + gap, 160 + gap)),
               )  # fmt: skip
    assert len(ins_signatures(rec, P)) == n_sigs


def test_min_ins_len_is_inclusive() -> None:
    def rec(n: int) -> PlacementRecord:
        return _rec("x", qs=0, qe=200 + n, start=0, end=200,
                    read_length=200 + n, sequence=_seq(200 + n, 21),
                    insertions=((100, 100, 100 + n),))  # fmt: skip

    assert len(ins_signatures(rec(50), P)) == 1
    assert ins_signatures(rec(49), P) == []


def test_nested_inserts_give_two_ins() -> None:
    a, i1, m, i2, b = (
        _seq(100, 22),
        _seq(60, 23),
        _seq(500, 24),
        _seq(70, 25),
        _seq(100, 26),
    )
    read = a + i1 + m + i2 + b
    # 100M60I500M70I100M at reference 1000.
    rec = _rec("n", qs=0, qe=len(read), start=1000, end=1700,
               read_length=len(read), sequence=read,
               insertions=((1100, 100, 160), (1600, 660, 730)))  # fmt: skip
    s1, s2 = ins_signatures(rec, P)
    assert (s1.locus.start, s1.insert_length, s1.sequence) == (1100, 60, i1)
    assert (s2.locus.start, s2.insert_length, s2.sequence) == (1600, 70, i2)
    assert s1.right_flank == m[:200]
    assert s2.left_flank == m[-200:]
    assert len(set(_ids([s1, s2]))) == 2


def test_ins_without_sequence_keeps_length_and_coordinates() -> None:
    rec = _rec("noseq", qs=0, qe=300, start=1000, end=1200, read_length=300,
               sequence=None, insertions=((1100, 100, 200),))  # fmt: skip
    (s,) = ins_signatures(rec, P)
    assert s.insert_length == 100
    assert (s.sequence, s.left_flank, s.right_flank) == ("", "", "")
    assert (s.query_start, s.query_end) == (100, 200)


# --- CLIP --------------------------------------------------------------------


def test_right_clip_500() -> None:
    aligned, clip = _seq(400, 30), _seq(500, 31)
    # 400M500S at reference 1000.
    rec = _rec("c", qs=0, qe=400, start=1000, end=1400, read_length=900,
               sequence=aligned + clip, rc=500, mapq=33)  # fmt: skip
    (s,) = clip_signatures(rec, P)
    assert s.kind is SignatureKind.CLIP
    assert s.clip_side is ClipSide.RIGHT
    assert s.locus == Locus(REF, "Chr3", 1400, 1400, "+")
    assert s.insert_length == 500
    assert s.sequence == clip
    assert s.left_flank == aligned[-200:]
    assert s.right_flank == ""
    assert (s.query_start, s.query_end) == (400, 900)
    assert s.mapq == 33
    assert ins_signatures(rec, P) == []


@pytest.mark.parametrize(("clip_len", "n_sigs"), [(99, 0), (100, 1)])
def test_clip_threshold_inclusive(clip_len: int, n_sigs: int) -> None:
    rec = _rec("c", qs=0, qe=400, start=1000, end=1400,
               read_length=400 + clip_len, sequence=_seq(400 + clip_len, 32),
               rc=clip_len)  # fmt: skip
    assert len(clip_signatures(rec, P)) == n_sigs


def test_left_clip_after_hard_clip_forward() -> None:
    hard, clip, aligned = _seq(30, 33), _seq(150, 34), _seq(400, 35)
    read = hard + clip + aligned
    # 30H150S400M at reference 2000.
    rec = _rec("lc", qs=180, qe=580, start=2000, end=2400, read_length=580,
               sequence=clip + aligned, lc=150, lh=30)  # fmt: skip
    (s,) = clip_signatures(rec, P)
    assert s.clip_side is ClipSide.LEFT
    assert s.locus == Locus(REF, "Chr3", 2000, 2000, "+")
    assert s.insert_length == 150
    assert s.sequence == clip
    assert s.left_flank == ""
    assert s.right_flank == aligned[:200]
    assert (s.query_start, s.query_end) == (30, 180)
    assert read[s.query_start : s.query_end] == clip


def test_left_clip_after_hard_clip_reverse() -> None:
    hard, clip, aligned = _seq(30, 36), _seq(150, 37), _seq(400, 38)
    refwd = hard + clip + aligned
    read = _rc(refwd)
    # Stored 30H150S400M on the reverse strand; read orientation (0, 400).
    rec = _rec("lcr", qs=0, qe=400, start=2000, end=2400, read_length=580,
               sequence=clip + aligned, strand="-", lc=150, lh=30)  # fmt: skip
    (s,) = clip_signatures(rec, P)
    assert s.clip_side is ClipSide.LEFT
    assert s.locus == Locus(REF, "Chr3", 2000, 2000, "-")
    assert s.sequence == clip
    assert (s.query_start, s.query_end) == (400, 550)
    assert _rc(s.sequence) == read[s.query_start : s.query_end]


def test_both_clips_and_right_hard_clip() -> None:
    lclip, aligned, rclip = _seq(120, 39), _seq(300, 40), _seq(110, 41)
    # 5H120S300M110S7H at reference 500; stored length 530, read 542.
    rec = _rec("both", qs=125, qe=425, start=500, end=800, read_length=542,
               sequence=lclip + aligned + rclip, lc=120, rc=110, lh=5,
               rh=7)  # fmt: skip
    left, right = clip_signatures(rec, P)
    assert (left.clip_side, right.clip_side) == (ClipSide.LEFT, ClipSide.RIGHT)
    assert (left.query_start, left.query_end) == (5, 125)
    assert (right.query_start, right.query_end) == (425, 535)
    assert (left.sequence, right.sequence) == (lclip, rclip)
    assert left.right_flank == aligned[:200]
    assert right.left_flank == aligned[-200:]
    assert (left.locus.start, right.locus.start) == (500, 800)
    assert left.signature_id != right.signature_id


@pytest.mark.parametrize(("strand", "interval"), [("+", (425, 535)), ("-", (7, 117))])
def test_right_clip_without_sequence(strand: str, interval: tuple[int, int]) -> None:
    # 5H120S300M110S7H with SEQ '*': stored length from L - lh - rh = 530.
    rec = _rec("ns", qs=0, qe=1, start=500, end=800, read_length=542,
               sequence=None, strand=strand, lc=120, rc=110, lh=5,
               rh=7)  # fmt: skip
    left, right = clip_signatures(rec, P)
    assert right.insert_length == 110
    assert (right.sequence, right.left_flank, right.right_flank) == ("", "", "")
    # Stored [420, 530): forward (5+420, 5+530); reverse (542-5-530, 542-5-420).
    assert (right.query_start, right.query_end) == interval


# --- SPLIT_INS -----------------------------------------------------------------


def _forward_split() -> tuple[str, str, str, str, PlacementRecord, PlacementRecord]:
    a, insert, b = _seq(1000, 50), _seq(800, 51), _seq(1200, 52)
    read = a + insert + b
    rec_a = _rec("s", qs=0, qe=1000, start=1000, end=2000, read_length=3000,
                 sequence=a + insert + b, rc=2000, mapq=60,
                 has_supp=True)  # fmt: skip
    rec_b = _rec("s", qs=1800, qe=3000, start=2010, end=3210, read_length=3000,
                 sequence=b, lh=1800, mapq=30, supp=True,
                 has_supp=True)  # fmt: skip
    return read, a, insert, b, rec_a, rec_b


def test_split_forward() -> None:
    read, a, insert, b, rec_a, rec_b = _forward_split()
    (s,) = split_ins_signatures([rec_b, rec_a], read, P)  # order-free input
    assert s.kind is SignatureKind.SPLIT_INS
    assert s.locus == Locus(REF, "Chr3", 2000, 2000, "+")
    assert s.insert_length == 800
    assert s.sequence == read[1000:1800] == insert
    assert s.left_flank == a[-200:]
    assert s.right_flank == b[:200]
    assert (s.query_start, s.query_end) == (1000, 1800)
    assert s.mapq == 30
    assert s.clip_side is ClipSide.NONE
    assert s.signature_id == "s:SPLIT_INS:1000-1800@msu7:Chr3:2000+"


_ReverseSplit = tuple[str, str, str, str, str, PlacementRecord, PlacementRecord]


def _reverse_split() -> _ReverseSplit:
    # Reference-forward picture: B's bases, the insert, then A's bases.
    b_fw, ins_fw, a_fw = _seq(1200, 60), _seq(800, 61), _seq(1000, 62)
    refwd = b_fw + ins_fw + a_fw
    read = _rc(refwd)  # as sequenced: A first, then the insert, then B
    # A: read (0, 1000), reference [2010, 3010), stored 2000S1000M.
    rec_a = _rec("v", qs=0, qe=1000, start=2010, end=3010, read_length=3000,
                 sequence=refwd, strand="-", lc=2000, mapq=50,
                 has_supp=True)  # fmt: skip
    # B: read (1800, 3000), reference [800, 2000), stored 1200M1800H.
    rec_b = _rec("v", qs=1800, qe=3000, start=800, end=2000, read_length=3000,
                 sequence=refwd[:1200], strand="-", rh=1800, mapq=40,
                 supp=True, has_supp=True)  # fmt: skip
    return read, refwd, b_fw, ins_fw, a_fw, rec_a, rec_b


def test_split_reverse() -> None:
    read, _refwd, b_fw, ins_fw, a_fw, rec_a, rec_b = _reverse_split()
    assert read[0:1000] == _rc(a_fw)  # A really is first in the read
    (s,) = split_ins_signatures([rec_a, rec_b], read, P)
    assert s.locus == Locus(REF, "Chr3", 2000, 2000, "-")
    assert s.insert_length == 800
    assert s.sequence == ins_fw == _rc(read[1000:1800])
    assert s.left_flank == b_fw[-200:]
    assert s.right_flank == a_fw[:200]
    assert (s.query_start, s.query_end) == (1000, 1800)
    assert s.mapq == 40
    # Independent placement on a genome: B sits on [800, 2000), A on
    # [2010, 3010), so the left flank is the reference just left of 2000.
    genome = _seq(800, 63) + b_fw + _seq(10, 64) + a_fw + _seq(500, 65)
    assert genome[1800:2000] == s.left_flank
    assert genome[2010:2210] == s.right_flank


def test_split_without_read_sequence() -> None:
    _read, _a, _insert, _b, rec_a, rec_b = _forward_split()
    (s,) = split_ins_signatures([rec_a, rec_b], None, P)
    assert s.insert_length == 800
    assert (s.sequence, s.left_flank, s.right_flank) == ("", "", "")


def _moved(rec: PlacementRecord, **locus: object) -> PlacementRecord:
    return replace(rec, locus=replace(rec.locus, **locus))


@pytest.mark.parametrize(
    ("b_start", "n_sigs"), [(2060, 0), (2050, 1), (1950, 1), (1940, 0)]
)
def test_split_ref_gap_limit(b_start: int, n_sigs: int) -> None:
    read, _a, _i, _b, rec_a, rec_b = _forward_split()
    rec_b = _moved(rec_b, start=b_start, end=b_start + 1200)
    assert len(split_ins_signatures([rec_a, rec_b], read, P)) == n_sigs


def test_split_needs_same_sequence_and_strand() -> None:
    read, _a, _i, _b, rec_a, rec_b = _forward_split()
    other_chrom = _moved(rec_b, sequence_id="Chr4")
    assert split_ins_signatures([rec_a, other_chrom], read, P) == []
    other_ref = _moved(rec_b, reference_id="other")
    assert split_ins_signatures([rec_a, other_ref], read, P) == []
    other_strand = _moved(rec_b, strand="-")
    assert split_ins_signatures([rec_a, other_strand], read, P) == []


def test_split_read_gap_threshold() -> None:
    read, _a, _i, _b, rec_a, _rec_b = _forward_split()
    rec_b49 = _rec("s", qs=1049, qe=3000, start=2010, end=3961, read_length=3000,
                   sequence=None, lh=1049)  # fmt: skip
    assert split_ins_signatures([rec_a, rec_b49], read, P) == []
    rec_b50 = _rec("s", qs=1050, qe=3000, start=2010, end=3960, read_length=3000,
                   sequence=None, lh=1050)  # fmt: skip
    (s,) = split_ins_signatures([rec_a, rec_b50], read, P)
    assert s.insert_length == 50


# --- collect_signatures --------------------------------------------------------


def test_collect_split_takes_sequence_from_full_primary_forward() -> None:
    read, a, insert, b, rec_a, rec_b = _forward_split()
    stats: dict[str, int] = {}
    sigs = list(collect_signatures([rec_b, rec_a], P, stats))
    split = [s for s in sigs if s.kind is SignatureKind.SPLIT_INS]
    assert len(split) == 1
    assert split[0].sequence == insert
    assert split[0].left_flank == a[-200:]
    assert split[0].right_flank == b[:200]
    # A's 2000 bp right soft clip also gives a CLIP.
    (clip,) = [s for s in sigs if s.kind is SignatureKind.CLIP]
    assert clip.clip_side is ClipSide.RIGHT
    assert (clip.query_start, clip.query_end) == (1000, 3000)
    assert stats == {"ins": 0, "split_ins": 1, "clip": 1, "buffered_reads": 1,
                     "split_ins_no_sequence": 0}  # fmt: skip


def test_collect_split_takes_sequence_from_full_primary_reverse() -> None:
    _read, _refwd, b_fw, ins_fw, a_fw, rec_a, rec_b = _reverse_split()
    sigs = list(collect_signatures([rec_b, rec_a], P))
    (split,) = [s for s in sigs if s.kind is SignatureKind.SPLIT_INS]
    assert split.sequence == ins_fw
    assert split.left_flank == b_fw[-200:]
    assert split.right_flank == a_fw[:200]
    assert split.locus.start == 2000


def test_collect_split_without_full_sequence() -> None:
    _read, _a, _i, _b, rec_a, rec_b = _forward_split()
    no_seq_primary = replace(rec_a, sequence=None)
    stats: dict[str, int] = {}
    sigs = list(collect_signatures([no_seq_primary, rec_b], P, stats))
    (split,) = [s for s in sigs if s.kind is SignatureKind.SPLIT_INS]
    assert (split.sequence, split.left_flank, split.right_flank) == ("", "", "")
    assert split.insert_length == 800
    assert stats["split_ins_no_sequence"] == 1


def test_ids_unique_when_two_records_share_query_span() -> None:
    # A read aligned twice over the same read region (e.g. a duplicated
    # segment): both records carry the same insertion and right clip in read
    # coordinates. Their ids must still differ.
    pre, insert, post, clip = _seq(100, 70), _seq(100, 71), _seq(100, 72), _seq(200, 73)
    read = pre + insert + post + clip
    common = {"qs": 0, "qe": 300, "read_length": 500, "sequence": read,
              "insertions": ((0, 100, 200),), "rc": 200, "has_supp": True}  # fmt: skip
    r1 = _rec(
        "dup", start=1000, end=1200, **{**common, "insertions": ((1100, 100, 200),)}
    )
    r2 = _rec("dup", start=5000, end=5200, supp=True,
              **{**common, "insertions": ((5100, 100, 200),)})  # fmt: skip
    sigs = list(collect_signatures([r1, r2], P))
    assert len(sigs) == 4
    assert len(set(_ids(sigs))) == 4


def _mixed_stream() -> list[PlacementRecord]:
    recs: list[PlacementRecord] = []
    # INS read.
    left, insert, right = _seq(300, 80), _seq(200, 81), _seq(300, 82)
    recs.append(_rec("ins_read", qs=0, qe=800, start=1000, end=1600,
                     read_length=800, sequence=left + insert + right,
                     insertions=((1300, 300, 500),)))  # fmt: skip
    # Clip-only read.
    recs.append(_rec("clip_read", qs=0, qe=400, start=4000, end=4400,
                     read_length=900, sequence=_seq(900, 83), rc=500))  # fmt: skip
    # Forward and reverse split reads (two records each, one with a CLIP).
    recs.extend(_forward_split()[4:])
    recs.extend(_reverse_split()[5:])
    # Nothing to report.
    recs.append(_rec("plain", qs=0, qe=500, start=7000, end=7500,
                     read_length=500, sequence=_seq(500, 84)))  # fmt: skip
    return recs


def test_collect_deterministic_under_shuffle_and_round_trips(tmp_path: Path) -> None:
    base = _mixed_stream()
    stats0: dict[str, int] = {}
    expected = sorted(collect_signatures(base, P, stats0), key=lambda s: s.signature_id)
    # ins_read: 1 INS; clip_read: 1 CLIP; forward split: 1 CLIP (A's right
    # clip) + 1 SPLIT_INS; reverse split: 1 CLIP (A's left clip) + 1 SPLIT_INS.
    assert stats0 == {"ins": 1, "split_ins": 2, "clip": 3, "buffered_reads": 2,
                      "split_ins_no_sequence": 0}  # fmt: skip
    assert len(expected) == 6
    assert len(set(_ids(expected))) == len(expected)
    for seed in range(5):
        shuffled = list(base)
        random.Random(seed).shuffle(shuffled)
        stats: dict[str, int] = {"ins": 10}  # counts are added, not replaced
        got = sorted(
            collect_signatures(shuffled, P, stats), key=lambda s: s.signature_id
        )
        assert got == expected
        assert stats == {**stats0, "ins": 11}

    out = tmp_path / "sigs.tsv.gz"
    n = write_signatures(out, expected, read_type=ReadType.ONT, reference_id=REF)
    assert n == len(expected)
    meta, rows = read_signatures(out)
    assert meta["reference_id"] == REF
    assert list(rows) == expected


# --- SPLIT_INS: nearest collinear partner (TE body aligned elsewhere) ----------


def _te_split_forward() -> tuple[str, str, list[PlacementRecord]]:
    a, insert, b = _seq(1000, 90), _seq(800, 91), _seq(1200, 92)
    read = a + insert + b
    rec_a = _rec("t", qs=0, qe=1000, start=1000, end=2000, read_length=3000,
                 sequence=read, rc=2000, has_supp=True)  # fmt: skip
    # The insert itself aligns to a TE copy elsewhere on the same sequence.
    rec_te = _rec("t", qs=1000, qe=1800, start=50000, end=50800,
                  read_length=3000, sequence=insert, lh=1000, rh=1200,
                  supp=True, has_supp=True)  # fmt: skip
    rec_b = _rec("t", qs=1800, qe=3000, start=2000, end=3200, read_length=3000,
                 sequence=b, lh=1800, supp=True, has_supp=True)  # fmt: skip
    return read, insert, [rec_a, rec_te, rec_b]


def test_split_skips_te_alignment_between_flanks_forward() -> None:
    read, insert, recs = _te_split_forward()
    (s,) = split_ins_signatures(recs, read, P)
    assert s.locus == Locus(REF, "Chr3", 2000, 2000, "+")
    assert s.insert_length == 800
    assert s.sequence == read[1000:1800] == insert
    assert (s.query_start, s.query_end) == (1000, 1800)


def test_split_skips_te_alignment_and_partial_te_alignment() -> None:
    read, insert, recs = _te_split_forward()
    partial = _rec("t", qs=1250, qe=1550, start=60000, end=60300,
                   read_length=3000, sequence=None, lh=1250, rh=1450,
                   supp=True, has_supp=True)  # fmt: skip
    for order in (recs + [partial], [partial] + recs[::-1]):
        (s,) = split_ins_signatures(order, read, P)
        assert (s.locus.start, s.insert_length, s.sequence) == (2000, 800, insert)


def test_split_skips_te_alignment_between_flanks_reverse() -> None:
    b_fw, ins_fw, a_fw = _seq(1200, 93), _seq(800, 94), _seq(1000, 95)
    refwd = b_fw + ins_fw + a_fw
    read = _rc(refwd)  # A, then the insert, then B
    rec_a = _rec("u", qs=0, qe=1000, start=2000, end=3000, read_length=3000,
                 sequence=refwd, strand="-", lc=2000, has_supp=True)  # fmt: skip
    rec_te = _rec("u", qs=1000, qe=1800, start=50000, end=50800,
                  read_length=3000, sequence=ins_fw, strand="-", lh=1200,
                  rh=1000, supp=True, has_supp=True)  # fmt: skip
    rec_b = _rec("u", qs=1800, qe=3000, start=800, end=2000, read_length=3000,
                 sequence=b_fw, strand="-", rh=1800, supp=True,
                 has_supp=True)  # fmt: skip
    (s,) = split_ins_signatures([rec_b, rec_te, rec_a], read, P)
    assert s.locus == Locus(REF, "Chr3", 2000, 2000, "-")
    assert s.insert_length == 800
    assert s.sequence == ins_fw == _rc(read[1000:1800])
    assert s.left_flank == b_fw[-200:]
    assert s.right_flank == a_fw[:200]


def test_split_three_collinear_alignments_give_two() -> None:
    # A, B, C collinear with an insert between each: pairs (A, B), (B, C).
    rec_a = _rec("c3", qs=0, qe=1000, start=1000, end=2000, read_length=3600,
                 sequence=None)  # fmt: skip
    rec_b = _rec("c3", qs=1300, qe=2300, start=2000, end=3000,
                 read_length=3600, sequence=None)  # fmt: skip
    rec_c = _rec("c3", qs=2700, qe=3600, start=3010, end=3910,
                 read_length=3600, sequence=None)  # fmt: skip
    sigs = split_ins_signatures([rec_c, rec_a, rec_b], None, P)
    got = [(s.locus.start, s.query_start, s.query_end, s.insert_length) for s in sigs]
    assert got == [(2000, 1000, 1300, 300), (3000, 2300, 2700, 400)]


# --- CLIP flanks stay inside the aligned bases -----------------------------------


def test_clip_flanks_do_not_reach_the_opposite_clip() -> None:
    lclip, aligned, rclip = _seq(150, 96), _seq(100, 97), _seq(150, 98)
    # 150S100M150S at reference 1000.
    rec = _rec("short", qs=150, qe=250, start=1000, end=1100, read_length=400,
               sequence=lclip + aligned + rclip, lc=150, rc=150)  # fmt: skip
    left, right = clip_signatures(rec, P)
    assert left.right_flank == aligned
    assert right.left_flank == aligned
    assert (left.sequence, right.sequence) == (lclip, rclip)


# --- collect: kept read sequence must have the read length ------------------------


def test_collect_ignores_full_sequence_of_wrong_length() -> None:
    _read, _a, _i, _b, rec_a, rec_b = _forward_split()
    bad = replace(rec_a, sequence=rec_a.sequence[:-1])  # 2999 bases, L = 3000
    stats: dict[str, int] = {}
    sigs = list(collect_signatures([bad, rec_b], P, stats))
    (split,) = [s for s in sigs if s.kind is SignatureKind.SPLIT_INS]
    assert (split.sequence, split.left_flank, split.right_flank) == ("", "", "")
    assert stats["split_ins_no_sequence"] == 1
