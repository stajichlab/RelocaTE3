"""Tests for insertion genotyping (Step 7)."""

from __future__ import annotations

from pathlib import Path

import pysam

from RelocaTE3.characterize import (
    characterize_insertions,
    classify_status,
    count_spanners,
)
from RelocaTE3.models import Insertion
from RelocaTE3.pipeline import run_sample
from RelocaTE3.ReadLibrary import ReadLibrary

DATA = Path(__file__).parent / "data"
R1 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"
R2 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"
TELIB = DATA / "mping.fa"
GENOME = DATA / "sim_genome" / "MSU7.Chr3_2M.fa"


def test_classify_status():
    assert classify_status(avg_flankers=5, spanners=0) == "homozygous"
    assert classify_status(avg_flankers=10, spanners=10) == "heterozygous"
    assert classify_status(avg_flankers=1, spanners=50) == "somatic_insertion"
    # many flankers, few spanners -> excision/homozygous bucket
    assert "homozygous" in classify_status(avg_flankers=20, spanners=1)


def _make_bam(tmp_path: Path, reads: list[tuple[int, int, str, int]]) -> Path:
    """Build a tiny indexed BAM from (start0, length, cigar, nm) tuples."""
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "Chr3", "LN": 2_000_000}]}
    bam_path = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for i, (start0, length, cigar, nm) in enumerate(reads):
            rec = pysam.AlignedSegment(bam.header)
            rec.query_name = f"r{i}"
            rec.query_sequence = "A" * length
            rec.flag = 0
            rec.reference_id = 0
            rec.reference_start = start0
            rec.mapping_quality = 60
            rec.cigarstring = cigar
            rec.set_tag("NM", nm)
            bam.write(rec)
    pysam.sort("-o", str(bam_path), str(bam_path))
    pysam.index(str(bam_path))
    return bam_path


def test_count_spanners(tmp_path: Path):
    pos = 1000
    bam_path = _make_bam(
        tmp_path,
        [
            (900, 200, "200M", 0),  # 901..1100 spans cleanly -> spanner
            (900, 200, "200M", 2),  # mismatches -> not a spanner
            (900, 200, "150M50S", 0),  # soft-clipped (multi-op CIGAR) -> not a spanner
            (998, 200, "200M", 0),  # starts after pos-5 -> not spanning
        ],
    )
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        assert count_spanners(bam, "Chr3", pos) == 1


def test_characterize_insertions_homozygous(tmp_path: Path):
    # no spanners at the site -> homozygous
    bam_path = _make_bam(
        tmp_path, [(900, 50, "50M", 0)]
    )  # ends at 950, doesn't span 2000
    ins = Insertion(
        "Chr3",
        1998,
        2000,
        "mPing",
        "+",
        "TAA",
        left_junction_reads=3,
        right_junction_reads=3,
    )
    characterize_insertions([ins], str(bam_path))
    assert ins.spanners == 0
    assert ins.status == "homozygous"


def test_run_sample_with_genotyping(tmp_path: Path):
    """End-to-end with genotyping writes a characterized GFF/TXT."""
    reads = ReadLibrary([str(R1), str(R2)], "HEG4")
    run_sample(reads, str(TELIB), str(GENOME), tmp_path, genotype=True, threads=4)
    char_gff = tmp_path / "results" / "HEG4.all_nonref_insert.characTErized.gff"
    char_txt = tmp_path / "results" / "HEG4.all_nonref_insert.characTErized.txt"
    assert char_gff.exists() and char_txt.exists()
    # at least some insertions were assigned a status
    statuses = [ln.split("\t")[7] for ln in char_txt.read_text().splitlines()[1:]]
    assert statuses
    assert all(s for s in statuses)
