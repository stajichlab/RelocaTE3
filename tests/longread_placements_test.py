"""Long-read BAM adapter: pysam alignments to PlacementRecords."""

from __future__ import annotations

import random
from pathlib import Path

import pysam
import pytest

from RelocaTE3.longread.model import Locus
from RelocaTE3.longread.placements import iter_placements, placement_from_segment

CHR3_LEN = 10000
_COMP = str.maketrans("ACGT", "TGCA")


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _header() -> pysam.AlignmentHeader:
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": "Chr3", "LN": CHR3_LEN}]}
    )


def _seq(n: int, seed: int = 1) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _seg(
    name: str,
    cigar: str,
    seq: str | None,
    *,
    start: int = 1000,
    flag: int = 0,
    mapq: int = 60,
    tags: list[tuple[str, object]] | None = None,
    header: pysam.AlignmentHeader | None = None,
) -> pysam.AlignedSegment:
    seg = pysam.AlignedSegment(header or _header())
    seg.query_name = name
    seg.flag = flag
    if flag & 4:
        seg.reference_id = -1
        seg.reference_start = -1
        seg.mapping_quality = 0
        seg.cigarstring = None
    else:
        seg.reference_id = 0
        seg.reference_start = start
        seg.mapping_quality = mapq
        seg.cigarstring = cigar
    if seq is not None:
        seg.query_sequence = seq
        seg.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    if tags:
        seg.set_tags(tags)
    return seg


def _bam(tmp_path: Path, segs: list[pysam.AlignedSegment]) -> Path:
    unsorted = tmp_path / "unsorted.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=_header()) as fh:
        for s in segs:
            fh.write(s)
    out = tmp_path / "reads.bam"
    pysam.sort("-o", str(out), str(unsorted))
    pysam.index(str(out))
    return out


# --- placement_from_segment -------------------------------------------------


def test_forward_spanning_insertion() -> None:
    left, insert, right = _seq(100, 1), _seq(200, 2), _seq(100, 3)
    read = left + insert + right
    rec = placement_from_segment(_seg("r1", "100M200I100M", read), "ref")
    assert rec.read_id == "r1"
    assert rec.insertions == ((1100, 100, 300),)
    # Independent check: the offsets slice out exactly the planted insert.
    ref_pos, qs, qe = rec.insertions[0]
    assert rec.sequence[qs:qe] == insert
    assert (rec.query_start, rec.query_end) == (0, 400)
    assert rec.read_length == 400
    assert (rec.left_clip, rec.right_clip) == (0, 0)
    assert (rec.left_hard_clip, rec.right_hard_clip) == (0, 0)
    assert rec.sequence == read
    assert rec.locus == Locus("ref", "Chr3", 1000, 1200, "+")
    assert rec.mapq == 60
    assert rec.is_supplementary is False
    assert rec.has_supplementary is False


def test_reverse_soft_clip_read_orientation() -> None:
    read = _seq(350, 4)  # as sequenced
    stored = _rc(read)  # BAM stores reverse-strand reads reverse-complemented
    rec = placement_from_segment(_seg("r2", "50S300M", stored, flag=16), "ref")
    assert rec.left_clip == 50
    assert rec.right_clip == 0
    assert rec.read_length == 350
    assert (rec.query_start, rec.query_end) == (0, 300)
    assert rec.locus == Locus("ref", "Chr3", 1000, 1300, "-")
    # Independent check: the aligned stored bases, reverse-complemented, are
    # the read-orientation slice [query_start, query_end) of the read.
    aligned_stored = rec.sequence[50:350]
    assert _rc(aligned_stored) == read[rec.query_start : rec.query_end]


@pytest.mark.parametrize(
    ("flag", "interval"), [(2048, (20, 320)), (2048 | 16, (80, 380))]
)
def test_hard_clipped_supplementary(flag: int, interval: tuple[int, int]) -> None:
    stored = _seq(380, 5)
    seg = _seg(
        "r3",
        "20H300M80S",
        stored,
        flag=flag,
        tags=[("SA", "Chr3,5000,+,300M100S,60,0;")],
    )
    rec = placement_from_segment(seg, "ref")
    assert rec.sequence == stored  # stored bases kept, not None
    assert len(rec.sequence) == 380
    assert rec.left_hard_clip == 20
    assert rec.right_hard_clip == 0
    assert rec.left_clip == 0
    assert rec.right_clip == 80
    assert rec.read_length == 400
    assert (rec.query_start, rec.query_end) == interval
    assert rec.is_supplementary is True
    assert rec.has_supplementary is True


def test_hard_clipped_read_orientation_matches_full_read() -> None:
    # Independent check of the reverse-strand formula: build the full read,
    # hard-clip it, and confirm the reported interval locates the aligned
    # bases in the read as sequenced.
    read = _seq(400, 6)
    full_stored = _rc(read)  # reverse strand, stored orientation
    # CIGAR 20H300M80S: drop the first 20 stored bases.
    stored = full_stored[20:]
    rec = placement_from_segment(_seg("r", "20H300M80S", stored, flag=16), "ref")
    aligned = stored[0:300]
    assert _rc(aligned) == read[rec.query_start : rec.query_end]
    assert (rec.query_start, rec.query_end) == (80, 380)


def test_hard_then_soft_at_both_ends() -> None:
    stored = _seq(5 + 100 + 7, 7)
    rec = placement_from_segment(_seg("r", "10H5S100M7S3H", stored), "ref")
    assert (rec.left_hard_clip, rec.left_clip) == (10, 5)
    assert (rec.right_clip, rec.right_hard_clip) == (7, 3)
    assert rec.read_length == 10 + 5 + 100 + 7 + 3
    assert (rec.query_start, rec.query_end) == (15, 115)
    assert rec.sequence == stored


def test_small_insertions_not_listed_but_shift_offsets() -> None:
    a, small, b, big, c = (
        _seq(50, 8),
        _seq(5, 9),
        _seq(50, 10),
        _seq(30, 11),
        _seq(50, 12),
    )
    read = a + small + b + big + c
    rec = placement_from_segment(_seg("r", "50M5I50M30I50M", read), "ref")
    assert rec.insertions == ((1100, 105, 135),)
    _, qs, qe = rec.insertions[0]
    assert rec.sequence[qs:qe] == big
    assert (rec.query_start, rec.query_end) == (0, len(read))
    # min_ins_piece is inclusive and adjustable.
    rec5 = placement_from_segment(
        _seg("r", "50M5I50M30I50M", read), "ref", min_ins_piece=5
    )
    assert rec5.insertions == ((1050, 50, 55), (1100, 105, 135))


def test_deletion_before_insertion_shifts_reference_pos() -> None:
    a, b, ins, c = _seq(100, 13), _seq(100, 14), _seq(60, 15), _seq(100, 16)
    read = a + b + ins + c
    seg = _seg("r", "100M20D100M60I100M", read)
    rec = placement_from_segment(seg, "ref")
    assert rec.insertions == ((1220, 200, 260),)
    assert rec.sequence[200:260] == ins
    # Independent check against pysam's own pairs: the reference position of
    # the query base right after the insertion is reference_pos.
    pairs = dict(seg.get_aligned_pairs())
    assert pairs[260] == 1220
    assert pairs[199] == 1219
    assert rec.locus.end == seg.reference_end == 1000 + 100 + 20 + 100 + 100


def test_soft_clip_offsets_insertion_into_stored_sequence() -> None:
    clip, a, ins, b = _seq(40, 17), _seq(100, 18), _seq(25, 19), _seq(100, 20)
    stored = clip + a + ins + b
    rec = placement_from_segment(_seg("r", "30H40S100M25I100M", stored), "ref")
    # Offsets index the stored sequence: soft clip counted, hard clip not.
    assert rec.insertions == ((1100, 140, 165),)
    assert rec.sequence[140:165] == ins
    assert (rec.query_start, rec.query_end) == (70, 295)


def test_reverse_strand_insertion_offsets_and_read_orientation() -> None:
    # Reverse version of 30H40S100M25I100M. Build the read as sequenced, then
    # store it as the BAM would: reverse-complemented, first 30 stored bases
    # hard-clipped.
    read = _seq(30 + 40 + 100 + 25 + 100, 23)
    full_stored = _rc(read)
    stored = full_stored[30:]
    ins_stored = full_stored[30 + 140 : 30 + 165]  # the planted I-op bases
    rec = placement_from_segment(_seg("r", "30H40S100M25I100M", stored, flag=16), "ref")
    assert rec.read_length == 295
    assert rec.insertions == ((1100, 140, 165),)
    _, qs, qe = rec.insertions[0]
    assert rec.sequence[qs:qe] == ins_stored
    # Task 6 conversion on a reverse record: read_length - left_hard_clip - q,
    # with start and end swapped.
    r_start = rec.read_length - rec.left_hard_clip - qe
    r_end = rec.read_length - rec.left_hard_clip - qs
    assert (r_start, r_end) == (100, 125)
    assert _rc(rec.sequence[qs:qe]) == read[r_start:r_end]
    # Aligned span, read orientation: stored (70, 295) -> (0, 225).
    assert (rec.query_start, rec.query_end) == (0, 225)


@pytest.mark.parametrize(
    ("cigar", "stored_len"),
    [("10S5S100M", 115), ("50M10S50M", 110), ("50M5H50M", 100), ("100M5H5S", 105)],
)
def test_inner_clips_rejected(cigar: str, stored_len: int) -> None:
    seg = _seg("bad_read", cigar, _seq(stored_len, 24))
    with pytest.raises(ValueError, match=r"bad_read.*" + cigar):
        placement_from_segment(seg, "ref")


def test_missing_cigar_raises_naming_read() -> None:
    seg = _seg("nocig", "100M", _seq(100, 25))
    seg.cigarstring = None
    with pytest.raises(ValueError, match="nocig"):
        placement_from_segment(seg, "ref")


def test_iter_placements_skips_missing_cigar(tmp_path: Path) -> None:
    nocig = _seg("nocig", "100M", _seq(100, 27), start=2000)
    nocig.cigarstring = None
    bam = _bam(tmp_path, [_seg("good", "300M", _seq(300, 26), start=1000), nocig])
    stats: dict[str, int] = {}
    recs = list(iter_placements(bam, "ref", min_mapq=20, stats=stats))
    assert [r.read_id for r in recs] == ["good"]
    assert stats["records"] == 2
    assert stats["no_cigar"] == 1


def test_missing_seq_gives_none_sequence() -> None:
    rec = placement_from_segment(_seg("r", "100M50I100M", None), "ref")
    assert rec.sequence is None
    assert rec.read_length == 250
    assert rec.insertions == ((1100, 100, 150),)


def test_equal_and_mismatch_ops_consume_both() -> None:
    read = _seq(150, 21)
    rec = placement_from_segment(_seg("r", "40=10X20=30I50=", read), "ref")
    assert rec.insertions == ((1070, 70, 100),)
    assert rec.locus.end == 1000 + 120
    assert (rec.query_start, rec.query_end) == (0, 150)


def test_n_op_consumes_reference_only() -> None:
    read = _seq(100 + 40 + 100, 22)
    rec = placement_from_segment(_seg("r", "100M500N100M40I100M", read), "ref")
    assert rec.insertions == ((1700, 200, 240),)


# --- iter_placements ---------------------------------------------------------


def test_iter_placements_filters_and_counts(tmp_path: Path) -> None:
    segs = [
        _seg("keep_fwd", "100M200I100M", _seq(400, 30), start=1000),
        _seg(
            "keep_hc",
            "20H300M80S",
            _seq(380, 31),
            start=3000,
            flag=2048,
            tags=[("SA", "Chr3,100,+,300M100S,60,0;")],
        ),  # fmt: skip
        _seg("secondary", "300M", _seq(300, 32), start=2000, flag=256),
        _seg("unmapped", "", _seq(300, 33), flag=4),
        _seg("lowq", "300M", _seq(300, 34), start=4000, mapq=5),
    ]
    random.Random(0).shuffle(segs)
    bam = _bam(tmp_path, segs)
    stats: dict[str, int] = {}
    recs = list(iter_placements(bam, "ref", min_mapq=20, stats=stats))
    assert [r.read_id for r in recs] == ["keep_fwd", "keep_hc"]
    assert recs[0].insertions == ((1100, 100, 300),)
    assert recs[1].left_hard_clip == 20
    assert stats == {
        "records": 5,
        "skipped_mapq": 1,
        "secondary": 1,
        "unmapped": 1,
        "hard_clipped": 1,
        "no_cigar": 0,
    }


def test_iter_placements_without_stats_accepts_str_path(tmp_path: Path) -> None:
    bam = _bam(tmp_path, [_seg("a", "300M", _seq(300, 40), mapq=20)])
    recs = list(iter_placements(str(bam), "ref", min_mapq=20))
    assert len(recs) == 1  # MAPQ equal to the threshold is kept


def test_iter_placements_reads_cram(tmp_path: Path) -> None:
    genome = _seq(CHR3_LEN, 50)
    fasta = tmp_path / "g.fa"
    fasta.write_text(">Chr3\n" + genome + "\n")
    pysam.faidx(str(fasta))
    read = genome[1000:1100] + _seq(200, 51) + genome[1100:1200]
    cram = tmp_path / "reads.cram"
    with pysam.AlignmentFile(
        str(cram), "wc", header=_header(), reference_filename=str(fasta)
    ) as fh:
        fh.write(_seg("c", "100M200I100M", read))
    recs = list(iter_placements(cram, "ref", min_mapq=20, reference_fasta=fasta))
    assert len(recs) == 1
    assert recs[0].insertions == ((1100, 100, 300),)
    assert recs[0].sequence == read


def test_iter_placements_cram_needs_reference(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="reference_fasta"):
        list(iter_placements(tmp_path / "x.cram", "ref", min_mapq=0))
