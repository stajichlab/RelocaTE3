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
    """Write one junction plus one admitted but unclassified TE hit."""
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

        weak = pysam.AlignedSegment(bam.header)
        weak.query_name = "weak_hit"
        weak.flag = 0
        weak.reference_id = 0
        weak.reference_start = 100
        weak.mapping_quality = 20
        weak.cigartuples = [(4, 10), (0, 31), (4, 110)]
        weak.query_sequence = "T" * READ_LEN
        weak.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
        weak.set_tag("NM", 0)
        bam.write(weak)
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

    # The weak hit satisfies no junction/middle trim branch. RelocaTE2 still
    # exposes it at the step-3/step-4 seam so it cannot be mistaken for a
    # genome-only mate.
    hit_names = outdir / "te_containing" / "syn.te_hit_names.txt"
    assert set(hit_names.read_text().splitlines()) == {"read1", "weak_hit"}
    containing = outdir / "te_containing" / "syn.left.ContainingReads.fq"
    with pysam.FastxFile(str(containing)) as records:
        assert {record.name for record in records} == {"read1:end:3", "weak_hit"}
    read_repeat = outdir / "te_containing" / "syn.read_repeat_name.txt"
    assert [line.split("\t")[0] for line in read_repeat.read_text().splitlines()] == [
        "read1"
    ]


def test_original_fastq_quality_replaces_synthetic_blat_quality(tmp_path):
    """Step 3 must trim the original quality string, as RelocaTE2 does."""
    bam = _write_bam(tmp_path)
    source = tmp_path / "R1.fastq"
    original_seq = "G" * FLANK_LEN + "T" * MATCH_LEN
    original_qual = "!" * FLANK_LEN + "J" * MATCH_LEN
    source.write_text(
        f"@read1\n{original_seq}\n+\n{original_qual}\n"
        f"@weak_hit\n{'T' * READ_LEN}\n+\n{'#' * READ_LEN}\n"
    )

    outdir = tmp_path / "restored"
    RelocaTE().write_trimmed_reads(
        name="syn",
        direction_bams=[("left", bam)],
        outdir=outdir,
        source_fastqs=[source],
        minimum_match_length=10,
        minimum_trimmed_length=10,
        mismatch_allowance=2,
    )

    flank = outdir / "flanking" / "syn.left.flankingReads.fq"
    with pysam.FastxFile(str(flank)) as records:
        record = next(records)
        assert record.name == "read1:end:3"
        assert record.sequence == "G" * FLANK_LEN
        assert record.quality == "!" * FLANK_LEN

    containing = outdir / "te_containing" / "syn.left.ContainingReads.fq"
    with pysam.FastxFile(str(containing)) as records:
        by_name = {record.name: record for record in records}
    assert by_name["read1:end:3"].quality == original_qual
    assert by_name["weak_hit"].quality == "#" * READ_LEN


def test_quality_restoration_uses_source_mate_for_unsuffixed_names(tmp_path):
    """Illumina-style equal R1/R2 names still resolve to the correct TE BAM."""
    source = tmp_path / "R2.fastq"
    source.write_text("@shared 2:N:0:ACGT\nACGT\n+\n!#%J\n")
    coord = {"shared/2": {"seq": "ACGT", "qual": "IIII"}}

    restored = RelocaTE._restore_original_fastq(coord, source, mate="2")

    assert restored == 1
    assert coord["shared/2"]["seq"] == "ACGT"
    assert coord["shared/2"]["qual"] == "!#%J"
