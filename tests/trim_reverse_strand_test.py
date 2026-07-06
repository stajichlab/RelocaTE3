"""Reverse-strand 3' junction: coord-frame translation regression guard.

Builds a 1-record BAM whose alignment cigar/flag mirror the real-world case
that surfaced in the rice validation (see
``plans/2026-07-02-trim-recall-parity.md``):

  - 151-bp read
  - reverse-strand alignment (``flag=16``) at TE position 395 (0-based 394)
  - cigar ``36M115S``; alignment covers TE positions 394..429 (TE length 430)
  - the 115-bp flank lives at the FASTQ-frame START of the read

Before the fix, R3's trim silently drops this read: ``_original_orientation``
returns the FASTQ-frame sequence, but ``_trim_record`` slices it with the
stored-frame ``qstart`` = 0 that pysam reported. ``seq[0:0]`` is empty, so
the ``len(trimmed_seq) >= len_cutoff_l`` gate rejects the record.

Expected after the fix: one flanking record is written, tagged ``:end:3``,
carrying the 115-bp FASTQ-frame flank.
"""

from __future__ import annotations

from pathlib import Path

import pysam

from RelocaTE3.librelocate import RelocaTE

TE_LEN = 430
READ_LEN = 151
MATCH_LEN = 36
FLANK_LEN = READ_LEN - MATCH_LEN  # 115


def _reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def _write_bam(tmp_path: Path) -> Path:
    """One reverse-strand read; TE portion at stored-frame start, flank at end."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"LN": TE_LEN, "SN": "mping"}],
    }
    raw = tmp_path / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as bam:
        r = pysam.AlignedSegment(bam.header)
        r.query_name = "read1"
        r.flag = 16  # mapped, reverse strand
        r.reference_id = 0
        r.reference_start = 394  # 0-based; 1-based POS = 395
        r.mapping_quality = 60
        r.cigartuples = [(0, MATCH_LEN), (4, FLANK_LEN)]  # 36M115S
        # BAM stores the reverse-complemented read: TE-matching bases first,
        # then the reverse-complemented flank. Use distinct bases so we can
        # verify the FASTQ-frame flank exactly.
        r.query_sequence = "A" * MATCH_LEN + "C" * FLANK_LEN
        r.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
        r.set_tag("NM", 0)
        bam.write(r)
    sorted_bam = tmp_path / "syn.bam"
    pysam.sort("-o", str(sorted_bam), str(raw))
    pysam.index(str(sorted_bam))
    return sorted_bam


def test_reverse_strand_three_prime_junction_emits_flank(tmp_path):
    """The 3' junction on a reverse-strand read must emit its 115-bp flank."""
    bam = _write_bam(tmp_path)
    outdir = tmp_path / "out"

    rt = RelocaTE()
    rt.write_trimmed_reads(
        name="syn",
        direction_bams=[("left", bam)],
        outdir=outdir,
        minimum_match_length=10,
        minimum_trimmed_length=10,
        mismatch_allowance=2,
    )

    flank_fq = outdir / "flanking" / "syn.left.flankingReads.fq"
    assert flank_fq.is_file(), f"missing {flank_fq}"
    lines = [line for line in flank_fq.read_text().splitlines() if line]
    # One record = 4 FASTQ lines (@header, seq, +, qual).
    assert len(lines) == 4, lines
    assert lines[0].endswith(":end:3"), lines[0]
    # FASTQ-frame flank: reverse-complement of the last FLANK_LEN bases of the
    # stored seq. Stored last FLANK_LEN bases are all "C" -> FASTQ frame all "G".
    expected_flank = _reverse_complement("C" * FLANK_LEN)
    assert len(lines[1]) == FLANK_LEN, len(lines[1])
    assert lines[1] == expected_flank, lines[1]
