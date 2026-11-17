"""Insertion signatures from placement records.

Three kinds of read evidence (:class:`SignatureKind`):

- ``INS``: CIGAR insertions inside one alignment, nearby pieces merged;
- ``SPLIT_INS``: a read gap between consecutive alignments of one read on
  the same sequence and strand, with little reference gap;
- ``CLIP``: a long soft clip at either end of an alignment.

Signature sequences and flanks are reference-forward; query coordinates are
in read orientation and count hard clips. A stored offset range ``[a, b)``
into ``PlacementRecord.sequence`` maps to ``(lh + a, lh + b)`` on a forward
record and ``(L - lh - b, L - lh - a)`` on a reverse one (``L`` the read
length, ``lh`` the left hard clip).

``signature_id`` is ``{read_id}:{kind}:{query_start}-{query_end}`` followed
by ``@{reference_id}:{sequence_id}:{position}{strand}``. The anchor makes ids
unique per read even when two alignments of one read cover the same read
interval, and the id never depends on record order.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass, replace

from RelocaTE3.longread.model import (
    ClipSide, Locus, PlacementRecord, Signature, SignatureKind,
)  # fmt: skip
from RelocaTE3.longread.reference import reverse_complement

STATS_KEYS = ("ins", "split_ins", "clip", "buffered_reads", "split_ins_no_sequence")


@dataclass(frozen=True)
class SignatureParams:
    """Thresholds for signature extraction.

    Attributes:
        min_ins_len: Shortest INS (after merging) or SPLIT_INS read gap kept.
        ins_merge_dist: Largest reference distance between consecutive CIGAR
            insertion pieces that are merged into one INS.
        min_clip_len: Shortest soft clip reported as a CLIP.
        split_max_ref_gap: Largest absolute reference gap between two
            alignments of a read that still counts as a SPLIT_INS.
        flank: Bases of read sequence kept on each side of an insert.
    """

    min_ins_len: int = 50
    ins_merge_dist: int = 20
    min_clip_len: int = 100
    split_max_ref_gap: int = 50
    flank: int = 200


def _signature_id(
    read_id: str, kind: SignatureKind, qs: int, qe: int, locus: Locus
) -> str:
    return (
        f"{read_id}:{kind.value}:{qs}-{qe}"
        f"@{locus.reference_id}:{locus.sequence_id}:{locus.start}{locus.strand}"
    )


def _point(rec_locus: Locus, pos: int, strand: str | None = None) -> Locus:
    return Locus(
        rec_locus.reference_id, rec_locus.sequence_id, pos, pos,
        strand or rec_locus.strand,
    )  # fmt: skip


def _read_interval(rec: PlacementRecord, a: int, b: int) -> tuple[int, int]:
    """Map stored offsets ``[a, b)`` to a read-orientation interval."""
    lh = rec.left_hard_clip
    if rec.locus.strand == "-":
        return rec.read_length - lh - b, rec.read_length - lh - a
    return lh + a, lh + b


def _make(
    rec: PlacementRecord,
    kind: SignatureKind,
    locus: Locus,
    a: int,
    b: int,
    side: ClipSide,
    seq: str,
    left: str,
    right: str,
) -> Signature:
    qs, qe = _read_interval(rec, a, b)
    return Signature(
        signature_id=_signature_id(rec.read_id, kind, qs, qe, locus),
        kind=kind,
        locus=locus,
        insert_length=b - a,
        clip_side=side,
        read_id=rec.read_id,
        query_start=qs,
        query_end=qe,
        mapq=rec.mapq,
        sequence=seq,
        left_flank=left,
        right_flank=right,
    )


def ins_signatures(rec: PlacementRecord, params: SignatureParams) -> list[Signature]:
    """Return INS signatures from one record's CIGAR insertions.

    Pieces are sorted by reference position; consecutive pieces at most
    ``ins_merge_dist`` apart on the reference are merged. A merged group
    spans stored ``[first.qstart, last.qend)`` (the aligned bases between
    pieces included), sits at the first piece's reference position, and is
    kept when that span is at least ``min_ins_len``.

    Args:
        rec: One placement record.
        params: Extraction thresholds.

    Returns:
        INS signatures in reference order. Sequence and flanks are ``""``
        when the record has no sequence.
    """
    pieces = sorted(rec.insertions)
    if not pieces:
        return []
    groups: list[list[tuple[int, int, int]]] = [[pieces[0]]]
    for piece in pieces[1:]:
        if piece[0] - groups[-1][-1][0] <= params.ins_merge_dist:
            groups[-1].append(piece)
        else:
            groups.append([piece])

    seq = rec.sequence
    flank = params.flank
    out = []
    for group in groups:
        pos, qs, qe = group[0][0], group[0][1], group[-1][2]
        if qe - qs < params.min_ins_len:
            continue
        if seq is None:
            ins = left = right = ""
        else:
            ins = seq[qs:qe]
            left = seq[max(0, qs - flank) : qs]
            right = seq[qe : qe + flank]
        out.append(_make(rec, SignatureKind.INS, _point(rec.locus, pos), qs, qe,
                         ClipSide.NONE, ins, left, right))  # fmt: skip
    return out


def clip_signatures(rec: PlacementRecord, params: SignatureParams) -> list[Signature]:
    """Return CLIP signatures for long soft clips at either stored end.

    A LEFT clip sits at ``locus.start``, a RIGHT clip at ``locus.end``. The
    inner flank is the aligned bases next to the clip; the outer flank is
    ``""``. LEFT and RIGHT clips of one record cover disjoint stored
    offsets, so their ids differ.

    Args:
        rec: One placement record.
        params: Extraction thresholds.

    Returns:
        Zero, one or two signatures, LEFT first. Sequence and flanks are
        ``""`` when the record has no sequence; offsets then use the stored
        length ``read_length - left_hard_clip - right_hard_clip``.
    """
    seq = rec.sequence
    flank = params.flank
    n = (
        len(seq) if seq is not None
        else rec.read_length - rec.left_hard_clip - rec.right_hard_clip
    )  # fmt: skip
    out = []
    lc = rec.left_clip
    if lc >= params.min_clip_len:
        clip = seq[:lc] if seq is not None else ""
        right = seq[lc : lc + flank] if seq is not None else ""
        out.append(_make(rec, SignatureKind.CLIP,
                         _point(rec.locus, rec.locus.start), 0, lc,
                         ClipSide.LEFT, clip, "", right))  # fmt: skip
    rc = rec.right_clip
    if rc >= params.min_clip_len:
        clip = seq[n - rc :] if seq is not None else ""
        left = seq[max(0, n - rc - flank) : n - rc] if seq is not None else ""
        out.append(_make(rec, SignatureKind.CLIP,
                         _point(rec.locus, rec.locus.end), n - rc, n,
                         ClipSide.RIGHT, clip, left, ""))  # fmt: skip
    return out


def _split_order(rec: PlacementRecord) -> tuple[int, int, int, int, int]:
    # Total order, so ties in query_start never depend on input order.
    return (rec.query_start, rec.query_end, rec.locus.start, rec.locus.end, rec.mapq)


def split_ins_signatures(
    recs: Iterable[PlacementRecord],
    read_sequence: str | None,
    params: SignatureParams,
) -> list[Signature]:
    """Return SPLIT_INS signatures between alignments of one read.

    Records are grouped by ``(locus.key, locus.strand)`` and ordered by
    read-orientation ``query_start``. For each consecutive pair (A, B) the
    read gap is ``B.query_start - A.query_end``; the reference gap and
    position are ``B.start - A.end`` at ``A.end`` on ``+`` and
    ``A.start - B.end`` at ``B.end`` on ``-``. A signature is emitted when
    ``abs(reference gap) <= split_max_ref_gap`` and the read gap is at least
    ``min_ins_len``.

    Args:
        recs: All kept records of one read (sequences are not needed).
        read_sequence: The full read in read orientation, or None.
        params: Extraction thresholds.

    Returns:
        Signatures in group order, then read order. Sequence and flanks are
        reference-forward, and ``""`` when ``read_sequence`` is None.
    """
    groups: dict[tuple[tuple[str, str], str], list[PlacementRecord]] = {}
    for rec in recs:
        groups.setdefault((rec.locus.key, rec.locus.strand), []).append(rec)

    flank = params.flank
    r = read_sequence
    out = []
    for (_key, strand), group in sorted(groups.items(), key=lambda kv: kv[0]):
        group.sort(key=_split_order)
        for a, b in zip(group, group[1:]):
            gap = b.query_start - a.query_end
            if strand == "+":
                ref_gap, pos = b.locus.start - a.locus.end, a.locus.end
            else:
                ref_gap, pos = a.locus.start - b.locus.end, b.locus.end
            if abs(ref_gap) > params.split_max_ref_gap or gap < params.min_ins_len:
                continue
            s, e = a.query_end, b.query_start
            if r is None:
                ins = left = right = ""
            elif strand == "+":
                ins, left, right = r[s:e], r[max(0, s - flank) : s], r[e : e + flank]
            else:
                ins = reverse_complement(r[s:e])
                left = reverse_complement(r[e : e + flank])
                right = reverse_complement(r[max(0, s - flank) : s])
            locus = _point(a.locus, pos, strand)
            out.append(Signature(
                signature_id=_signature_id(a.read_id, SignatureKind.SPLIT_INS,
                                           s, e, locus),
                kind=SignatureKind.SPLIT_INS,
                locus=locus,
                insert_length=gap,
                clip_side=ClipSide.NONE,
                read_id=a.read_id,
                query_start=s,
                query_end=e,
                mapq=min(a.mapq, b.mapq),
                sequence=ins,
                left_flank=left,
                right_flank=right,
            ))  # fmt: skip
    return out


def collect_signatures(
    placements: Iterable[PlacementRecord],
    params: SignatureParams,
    stats: dict[str, int] | None = None,
) -> Iterator[Signature]:
    """Stream INS, CLIP and SPLIT_INS signatures from placement records.

    Each record's INS and CLIP signatures are yielded as it arrives. Records
    with ``has_supplementary`` are also buffered by read, without their
    sequence; the read-orientation full sequence is kept from the first
    record of the read with no hard clips and a sequence. SPLIT_INS
    signatures follow at the end of the stream, reads in sorted id order.

    Args:
        placements: Placement records, in any order.
        params: Extraction thresholds.
        stats: If given, counts ``ins``, ``split_ins``, ``clip``,
            ``buffered_reads`` and ``split_ins_no_sequence`` (SPLIT_INS
            emitted without a read sequence) are added.

    Yields:
        Signatures; the set yielded does not depend on record order.
    """
    counts = dict.fromkeys(STATS_KEYS, 0)
    buffered: dict[str, list[PlacementRecord]] = {}
    full_seq: dict[str, str] = {}
    try:
        for rec in placements:
            for sig in ins_signatures(rec, params):
                counts["ins"] += 1
                yield sig
            for sig in clip_signatures(rec, params):
                counts["clip"] += 1
                yield sig
            if not rec.has_supplementary:
                continue
            buffered.setdefault(rec.read_id, []).append(
                replace(rec, sequence=None, insertions=())
            )
            if (
                rec.read_id not in full_seq
                and rec.sequence is not None
                and rec.left_hard_clip == 0
                and rec.right_hard_clip == 0
            ):
                full_seq[rec.read_id] = (
                    reverse_complement(rec.sequence)
                    if rec.locus.strand == "-" else rec.sequence
                )  # fmt: skip
        counts["buffered_reads"] = len(buffered)
        for read_id in sorted(buffered):
            read_seq = full_seq.get(read_id)
            for sig in split_ins_signatures(buffered[read_id], read_seq, params):
                counts["split_ins"] += 1
                if read_seq is None:
                    counts["split_ins_no_sequence"] += 1
                yield sig
    finally:
        if stats is not None:
            for key, value in counts.items():
                stats[key] = stats.get(key, 0) + value
