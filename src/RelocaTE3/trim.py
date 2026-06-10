"""Step 3: identify TE-containing reads and trim the TE-matching portion.

This is a faithful port of RelocaTE2's ``relocaTE_trim.py`` (the ``parse_align_bwa``
scoring and the FASTQ trimming branches), restructured around the dataclasses in
``RelocaTE3.models``.

Note on orientation: RelocaTE2 read sequences from the original FASTQ and sliced
them with alignment coordinates. RelocaTE3 instead works from the read sequence as
stored in the BAM record (``pysam`` ``query_sequence`` / ``query_alignment_start``),
which is self-consistent for extracting the flanking sequence to re-map. Exact
byte-for-byte parity with RelocaTE2 read naming is validated at the insertion-calling
stage.
"""

from __future__ import annotations

from pathlib import Path

import pysam

from RelocaTE3.models import JunctionType, TEReadAlignment, TrimmedRead

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of ``seq``."""
    return seq.translate(_COMPLEMENT)[::-1]


def _convert_tags(tags: list) -> dict:
    """Convert a pysam tag list into a dict."""
    return {t[0]: t[1] for t in tags}


def _boundary_score(
    qstart: int, qend: int, qlen: int, tstart: int, tend: int, tlen: int
) -> int:
    """Count how many read/TE ends the alignment reaches (0-4).

    A boundary of 2 is the maximum we expect for a true junction: either one
    read end and one TE end, or both read ends with the TE internal.
    """
    score = 0
    if qstart <= 2:
        score += 1
    if qend + 1 == qlen or qend >= qlen - 3:
        score += 1
    if tstart <= 2:
        score += 1
    if tend + 1 == tlen or tend >= tlen - 3:
        score += 1
    return score


def parse_te_alignments(bamfiles: list[Path]) -> dict[str, TEReadAlignment]:
    """Parse read-to-TE BAM file(s), keeping the best alignment per read.

    Returns a mapping of read name to its best :class:`TEReadAlignment`.
    """
    coord: dict[str, TEReadAlignment] = {}
    for bam in bamfiles:
        fbam = pysam.AlignmentFile(str(bam), "rb")
        refnames = fbam.references
        reflengths = fbam.lengths
        qlen_c = 0
        for record in fbam.fetch(until_eof=True):
            if record.is_unmapped:
                continue
            qname = record.query_name
            qlen = int(record.query_length)
            qstart = int(record.query_alignment_start)
            try:
                qend = int(record.query_alignment_end) - 1
            except (ValueError, TypeError):
                continue
            if qlen == 0:
                # hard-clipped primary record; fall back to last seen length
                qlen = qlen_c
                qend = qlen + qend
            else:
                qlen_c = qlen

            tname = refnames[record.reference_id]
            tlen = int(reflengths[record.reference_id])
            tstart = int(record.reference_start)
            tend = int(record.reference_end) - 1

            tags = _convert_tags(record.tags if record.tags else [])
            match = ins = dele = 0
            for key, length in record.cigartuples or []:
                if key == 0:
                    match += length
                elif key == 1:
                    ins += length
                elif key == 2:
                    dele += length
            mismatch = int(tags.get("NM", 0)) - ins - dele
            match -= mismatch

            strand = "-" if record.is_reverse else "+"
            boundary = _boundary_score(qstart, qend, qlen, tstart, tend, tlen)

            aln = TEReadAlignment(
                read_name=qname,
                read_length=qlen,
                qstart=qstart,
                qend=qend,
                te_name=tname,
                te_length=tlen,
                tstart=tstart,
                tend=tend,
                match=match,
                mismatch=mismatch,
                strand=strand,
                boundary=boundary,
                seq=record.query_sequence or "",
                qual=_qual_string(record),
            )
            existing = coord.get(qname)
            if existing is None or aln.is_better_than(existing):
                coord[qname] = aln
        fbam.close()
    return coord


def _qual_string(record: pysam.AlignedSegment) -> str:
    """Return the ASCII (phred+33) quality string for a record, or ''."""
    quals = record.query_qualities
    if quals is None:
        return ""
    return "".join(chr(q + 33) for q in quals)


def trim_read(
    aln: TEReadAlignment,
    len_cut_match: int = 10,
    len_cut_trim: int = 10,
    mismatch_allowance: int = 0,
) -> TrimmedRead | None:
    """Classify a TE alignment and return the trimmed flanking read, or None.

    Reproduces the three RelocaTE2 cases:

    * 5' junction (``tstart <= 2``): the read crosses the TE 5' end.
    * 3' junction (``tend >= tlen - 3``): the read crosses the TE 3' end.
    * middle (read wholly inside the TE): a supporting read, kept untrimmed.
    """
    seq, qual = aln.seq, aln.qual
    start, end, length = aln.qstart, aln.qend, aln.read_length
    if not seq:
        return None

    read_at_boundary = start <= 2 or end >= length - 3
    match_long_enough = (aln.match + aln.mismatch) >= len_cut_match
    mismatch_ok = aln.mismatch <= mismatch_allowance

    # 5' junction
    if aln.tstart <= 2 and read_at_boundary and match_long_enough and mismatch_ok:
        te_portion = seq[start : end + 1]
        if aln.strand == "-":
            te_portion = reverse_complement(te_portion)
            flank_seq, flank_qual, side = seq[end + 1 :], qual[end + 1 :], "start"
        else:
            flank_seq, flank_qual, side = seq[:start], qual[:start], "end"
        if len(flank_seq) >= len_cut_trim:
            return TrimmedRead(
                read_name=f"{aln.read_name}:{side}:5",
                seq=flank_seq,
                qual=flank_qual,
                te_name=aln.te_name,
                strand=aln.strand,
                junction_type=JunctionType.FIVE_PRIME,
                te_portion=te_portion,
            )
        return None

    # 3' junction
    if (
        aln.tend >= aln.te_length - 3
        and read_at_boundary
        and match_long_enough
        and mismatch_ok
    ):
        te_portion = seq[start : end + 1]
        if aln.strand == "-":
            te_portion = reverse_complement(te_portion)
            flank_seq, flank_qual, side = seq[:start], qual[:start], "end"
        else:
            flank_seq, flank_qual, side = seq[end + 1 :], qual[end + 1 :], "start"
        if len(flank_seq) >= len_cut_trim:
            return TrimmedRead(
                read_name=f"{aln.read_name}:{side}:3",
                seq=flank_seq,
                qual=flank_qual,
                te_name=aln.te_name,
                strand=aln.strand,
                junction_type=JunctionType.THREE_PRIME,
                te_portion=te_portion,
            )
        return None

    # middle: read lies entirely within the TE -> supporting read, untrimmed
    if start <= 2 and end + 1 >= length - 3:
        return TrimmedRead(
            read_name=f"{aln.read_name}:middle",
            seq=seq,
            qual=qual,
            te_name=aln.te_name,
            strand=aln.strand,
            junction_type=JunctionType.MIDDLE,
        )

    return None
