"""Class-path TSD parity test: wildcard TSD mode emits read-captured 3-mers."""

from pathlib import Path

import pysam

from RelocaTE3.insertions import InsertionFinder


def _write_synthetic_bam(tmp_path: Path) -> Path:
    """Build a 1-chrom BAM with two junction reads framing one TAA TSD at 100..102."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"LN": 1000, "SN": "Chr1"}],
    }
    raw = tmp_path / "syn.raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as bam:
        # Right-junction read: flank covers 1-based 100..109; first 3 bases = TAA.
        r1 = pysam.AlignedSegment(bam.header)
        r1.query_name = "read1:start:5"
        r1.flag = 0  # mapped, forward strand
        r1.reference_id = 0
        r1.reference_start = 99  # 0-based -> 1-based 100
        r1.mapping_quality = 60
        r1.cigartuples = [(0, 10)]  # 10M
        r1.query_sequence = "TAAGGGCCAA"
        r1.query_qualities = pysam.qualitystring_to_array("I" * 10)
        r1.set_tag("NM", 0)
        bam.write(r1)
        # Left-junction read: flank covers 1-based 93..102; last 3 bases = TAA.
        r2 = pysam.AlignedSegment(bam.header)
        r2.query_name = "read2:end:3"
        r2.flag = 0
        r2.reference_id = 0
        r2.reference_start = 92  # 0-based -> 1-based 93
        r2.mapping_quality = 60
        r2.cigartuples = [(0, 10)]
        r2.query_sequence = "AAAAAAATAA"
        r2.query_qualities = pysam.qualitystring_to_array("I" * 10)
        r2.set_tag("NM", 0)
        bam.write(r2)
    sorted_bam = tmp_path / "syn.bam"
    pysam.sort("-o", str(sorted_bam), str(raw))
    pysam.index(str(sorted_bam))
    return sorted_bam


def _write_read_repeat(tmp_path: Path) -> Path:
    p = tmp_path / "read_repeat.txt"
    p.write_text("read1\tmPing\t+\nread2\tmPing\t-\n")
    return p


def test_class_path_wildcard_tsd_captures_read_bases(tmp_path):
    """With tsd='...', InsertionFinder must emit the literal 3 bases each read carries."""
    bam = _write_synthetic_bam(tmp_path)
    read_repeat = _write_read_repeat(tmp_path)
    outdir = tmp_path / "out"
    finder = InsertionFinder(mismatch_allow=2, min_mapq=1)
    out_txt = finder.find_insertions(
        bam_file=bam,
        read_repeat_file=read_repeat,
        tsd="...",  # 3-bp wildcard: capture whatever the junction read shows
        target="Chr1",
        sample="syn",
        outdir=outdir,
        te_name="mPing",
    )
    rows = [
        line.split("\t")
        for line in Path(out_txt).read_text().splitlines()
        if line and not line.lower().startswith("strain")
    ]
    assert len(rows) == 1, rows
    # _emit format (insertions.py:430): TE, TSD, sample, chrom, coord, strand, ...
    assert rows[0][0] == "mPing", rows
    assert rows[0][1] == "TAA", rows
    assert rows[0][3] == "Chr1", rows
    assert rows[0][4] == "100..102", rows


def _write_single_sided_bam(tmp_path: Path) -> Path:
    """Two right-junction reads framing a TAA TSD at Chr1:100..102; no left side."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"LN": 1000, "SN": "Chr1"}],
    }
    raw = tmp_path / "syn1s.raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as bam:
        for i in range(2):
            r = pysam.AlignedSegment(bam.header)
            r.query_name = f"r{i}:start:5"
            r.flag = 0
            r.reference_id = 0
            r.reference_start = 99  # 0-based -> 1-based 100
            r.mapping_quality = 60
            r.cigartuples = [(0, 10)]
            r.query_sequence = "TAAGGGCCAA"
            r.query_qualities = pysam.qualitystring_to_array("I" * 10)
            r.set_tag("NM", 0)
            bam.write(r)
    sorted_bam = tmp_path / "syn1s.bam"
    pysam.sort("-o", str(sorted_bam), str(raw))
    pysam.index(str(sorted_bam))
    return sorted_bam


def test_class_path_emits_captured_tsd_for_single_sided_junction(tmp_path):
    """Single-sided junction in wildcard TSD mode emits top_tsd, not supporting_junction."""
    bam = _write_single_sided_bam(tmp_path)
    repeat = tmp_path / "read_repeat.txt"
    repeat.write_text("r0\tmPing\t+\nr1\tmPing\t+\n")
    outdir = tmp_path / "out"
    finder = InsertionFinder(mismatch_allow=2, min_mapq=1)
    out_txt = finder.find_insertions(
        bam_file=bam,
        read_repeat_file=repeat,
        tsd="...",
        target="Chr1",
        sample="syn",
        outdir=outdir,
        te_name="mPing",
    )
    rows = [
        line.split("\t")
        for line in Path(out_txt).read_text().splitlines()
        if line and not line.lower().startswith("strain")
    ]
    assert len(rows) == 1, rows
    # _emit format: TE, TSD, sample, chrom, coord, strand, T:N, R:N, L:N, ST, SR, SL
    assert rows[0][0] == "mPing"
    assert rows[0][1] == "TAA", rows  # was "supporting_junction" before the fix
    assert rows[0][7] == "R:2"
    assert rows[0][8] == "L:0"
