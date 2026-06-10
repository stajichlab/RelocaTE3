"""Tests for reference TE handling (Steps 0/6)."""

from __future__ import annotations

from pathlib import Path

from RelocaTE3.models import Insertion
from RelocaTE3.reference_te import (
    ReferenceTE,
    filter_reference_overlaps,
    find_reference_insertions,
    parse_repeatmasker,
    write_existing_te_bed,
)

DATA = Path(__file__).parent / "data"
RM = DATA / "sim_genome" / "MSU7.Chr3_2M.fa.RepeatMasker.out"


def test_parse_repeatmasker_real_file():
    tes = parse_repeatmasker(RM)
    assert len(tes) > 100
    # header rows must be skipped; every record has a chromosome + coordinates
    assert all(te.chrom and te.end >= te.start for te in tes)
    # the file has both strands and a mix of intact / fragmented elements
    assert {te.strand for te in tes} == {"+", "-"}
    assert any(te.intact for te in tes)
    # first annotated element is TYPEU4 on Chr3 at 7192-7483 (reverse strand)
    first = tes[0]
    assert first.chrom == "Chr3" and first.start == 7192 and first.strand == "-"


def test_write_existing_te_bed(tmp_path: Path):
    tes = [
        ReferenceTE("Chr3", 100, 200, "mPing", "DNA/MITE", "+", True),
        ReferenceTE("Chr3", 500, 540, "Frag", "DNA", "-", False),
    ]
    bed = tmp_path / "existingTE.bed"
    write_existing_te_bed(tes, bed)
    lines = bed.read_text().splitlines()
    assert lines[0] == "Chr3\t100\t200\tmPing:100-200\t1\t+"
    assert lines[1].endswith("\t0\t-")


def test_filter_reference_overlaps():
    tes = [ReferenceTE("Chr3", 1000, 1200, "mPing", "DNA/MITE", "+", True)]
    inside = Insertion(
        "Chr3", 1100, 1102, "mPing", "+", "TAA"
    )  # same family, overlaps -> drop
    other = Insertion(
        "Chr3", 1100, 1102, "Other", "+", "TAA"
    )  # different family -> keep
    far = Insertion("Chr3", 9000, 9002, "mPing", "+", "TAA")  # no overlap -> keep
    kept = filter_reference_overlaps([inside, other, far], tes)
    assert inside not in kept
    assert other in kept and far in kept


def test_find_reference_insertions(tmp_path: Path):
    """Junction reads at an intact reference TE boundary -> a shared insertion."""
    import pysam

    tes = [ReferenceTE("Chr3", 1000, 1500, "mPing", "DNA/MITE", "+", True)]
    # build a tiny BAM with two junction reads at the TE's left (1000) and right (1500) edges
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "Chr3", "LN": 2_000_000}]}
    bam_path = tmp_path / "g.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for name, pos in (("a/1:end:5", 900), ("b/1:start:5", 1500)):
            rec = pysam.AlignedSegment(bam.header)
            rec.query_name = name
            rec.query_sequence = "A" * 100
            rec.flag = 0
            rec.reference_id = 0
            rec.reference_start = pos
            rec.mapping_quality = 60
            rec.cigartuples = [(0, 100)]
            bam.write(rec)
    pysam.index(str(bam_path))

    read_repeat = {"a/1:end:5": ("mPing", "+"), "b/1:start:5": ("mPing", "+")}
    calls = find_reference_insertions(str(bam_path), read_repeat, tes)
    assert len(calls) == 1
    assert calls[0].te_name == "mPing"
    assert calls[0].note.startswith("Shared")
    assert calls[0].start == 1000 and calls[0].end == 1500
