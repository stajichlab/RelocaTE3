"""BAM/CRAM adapter: pysam alignments to :class:`PlacementRecord`.

This is the only long-read module that touches pysam alignment objects.
Everything downstream works on :class:`PlacementRecord`, so another format
(e.g. GAF) only needs its own adapter producing the same records.

CIGAR ops: M/=/X consume query and reference, I and S consume query only,
D and N consume reference only, H consumes neither (the bases are absent
from the stored sequence but count toward the read length).
"""

from __future__ import annotations

import os
from collections.abc import Iterator

import pysam

from RelocaTE3.longread.model import Locus, PlacementRecord

# pysam CIGAR operation codes.
_M, _I, _D, _N, _S, _H, _EQ, _X = 0, 1, 2, 3, 4, 5, 7, 8
_QUERY_AND_REF = (_M, _EQ, _X)
_REF_ONLY = (_D, _N)

STATS_KEYS = ("records", "skipped_mapq", "secondary", "unmapped", "hard_clipped")


def _end_clips(cigar: list[tuple[int, int]]) -> tuple[int, int, int, int, int, int]:
    """Return clip lengths and the index range of the inner CIGAR ops.

    Returns:
        ``(left_hard, left_soft, right_soft, right_hard, lo, hi)`` where
        ``cigar[lo:hi]`` are the ops between the end clips.
    """
    lo, hi = 0, len(cigar)
    left_hard = left_soft = right_soft = right_hard = 0
    if lo < hi and cigar[lo][0] == _H:
        left_hard = cigar[lo][1]
        lo += 1
    if lo < hi and cigar[lo][0] == _S:
        left_soft = cigar[lo][1]
        lo += 1
    if hi > lo and cigar[hi - 1][0] == _H:
        right_hard = cigar[hi - 1][1]
        hi -= 1
    if hi > lo and cigar[hi - 1][0] == _S:
        right_soft = cigar[hi - 1][1]
        hi -= 1
    return left_hard, left_soft, right_soft, right_hard, lo, hi


def placement_from_segment(
    seg: pysam.AlignedSegment, reference_id: str, min_ins_piece: int = 10
) -> PlacementRecord:
    """Convert one mapped alignment into a :class:`PlacementRecord`.

    Args:
        seg: A mapped (primary or supplementary) alignment.
        reference_id: Identifier of the reference the BAM was aligned to.
        min_ins_piece: Shortest CIGAR ``I`` op listed in ``insertions``.

    Returns:
        The record. ``query_start``/``query_end`` are in read orientation and
        count hard clips; ``insertions`` holds ``(reference_pos, stored_qstart,
        stored_qend)`` with offsets into ``seg.query_sequence``.
    """
    cigar = seg.cigartuples or []
    left_hard, left_soft, right_soft, right_hard, lo, hi = _end_clips(cigar)

    ref_pos = seg.reference_start
    stored_q = left_soft  # offset into the stored sequence (no hard clips)
    aligned_query = 0
    insertions: list[tuple[int, int, int]] = []
    for op, length in cigar[lo:hi]:
        if op in _QUERY_AND_REF:
            ref_pos += length
            stored_q += length
            aligned_query += length
        elif op == _I:
            if length >= min_ins_piece:
                insertions.append((ref_pos, stored_q, stored_q + length))
            stored_q += length
            aligned_query += length
        elif op in _REF_ONLY:
            ref_pos += length
        elif op == _S:
            # An inner soft clip is invalid SAM; count it as query-only so the
            # offsets still match the stored sequence.
            stored_q += length
        # H inside the CIGAR is invalid SAM and consumes nothing; P likewise.

    read_length = seg.infer_read_length() or 0
    qs = left_hard + left_soft
    qe = qs + aligned_query
    if seg.is_reverse:
        qs, qe = read_length - qe, read_length - qs

    return PlacementRecord(
        read_id=seg.query_name,
        read_length=read_length,
        query_start=qs,
        query_end=qe,
        locus=Locus(
            reference_id,
            seg.reference_name,
            seg.reference_start,
            seg.reference_end,
            "-" if seg.is_reverse else "+",
        ),
        is_supplementary=seg.is_supplementary,
        mapq=seg.mapping_quality,
        insertions=tuple(insertions),
        left_clip=left_soft,
        right_clip=right_soft,
        sequence=seg.query_sequence,
        has_supplementary=seg.has_tag("SA"),
        left_hard_clip=left_hard,
        right_hard_clip=right_hard,
    )


def iter_placements(
    bam_path: str | os.PathLike[str],
    reference_id: str,
    *,
    min_mapq: int,
    stats: dict[str, int] | None = None,
    reference_fasta: str | os.PathLike[str] | None = None,
    min_ins_piece: int = 10,
) -> Iterator[PlacementRecord]:
    """Yield a record for each primary or supplementary alignment in a file.

    Unmapped and secondary records, and records with MAPQ below ``min_mapq``,
    are skipped. Records are read in file order; no index is needed.

    Args:
        bam_path: BAM, or CRAM when the suffix is ``.cram``.
        reference_id: Identifier of the reference the reads were aligned to.
        min_mapq: Lowest MAPQ kept.
        stats: If given, counts ``records`` (all records read),
            ``skipped_mapq``, ``secondary``, ``unmapped`` and
            ``hard_clipped`` (yielded records with any hard clip) are added.
        reference_fasta: Reference FASTA; required to decode a CRAM.
        min_ins_piece: Passed to :func:`placement_from_segment`.

    Yields:
        One :class:`PlacementRecord` per kept alignment.

    Raises:
        ValueError: A CRAM was given without ``reference_fasta``.
    """
    counts = dict.fromkeys(STATS_KEYS, 0)
    path = os.fspath(bam_path)
    if path.lower().endswith(".cram"):
        if reference_fasta is None:
            raise ValueError(f"reading CRAM {path} needs reference_fasta")
        fh = pysam.AlignmentFile(
            path, "rc", reference_filename=os.fspath(reference_fasta)
        )
    else:
        fh = pysam.AlignmentFile(path, "rb")
    try:
        with fh:
            for seg in fh.fetch(until_eof=True):
                counts["records"] += 1
                if seg.is_unmapped:
                    counts["unmapped"] += 1
                    continue
                if seg.is_secondary:
                    counts["secondary"] += 1
                    continue
                if seg.mapping_quality < min_mapq:
                    counts["skipped_mapq"] += 1
                    continue
                rec = placement_from_segment(seg, reference_id, min_ins_piece)
                if rec.left_hard_clip or rec.right_hard_clip:
                    counts["hard_clipped"] += 1
                yield rec
    finally:
        if stats is not None:
            for key, value in counts.items():
                stats[key] = stats.get(key, 0) + value
