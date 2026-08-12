"""Mate pairing must not depend on the FASTQ's read-name convention.

RelocaTE3 re-derives which reads are mates by string-matching a trailing ``/1``
or ``/2`` on the read name. Reads straight off a modern Illumina instrument put
the mate in a separate field instead (``@INSTRUMENT:... 1:N:0:INDEX``), so the
name is *identical* in R1 and R2; SRA/ENA dumps often carry no mate marker at
all. On those inputs the mate lookup silently found nothing: the run exited 0
and produced plausible output with almost no supporting reads (measured on the
Chr3 fixture: 3955 supporting reads -> 17, and 158/158 calls with support ->
4/154).

The mate is unambiguous from *which file* a read came from, so pairing is keyed
on that rather than on the name.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from RelocaTE3.genome_align import (
    build_flank_pairs,
    canonical_name,
    recover_support_mates,
)
from RelocaTE3.ReadLibrary import ReadLibrary


def _write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


# ---------------------------------------------------------------------------
# canonical_name
# ---------------------------------------------------------------------------


def test_canonical_name_handles_every_convention():
    """All conventions collapse to the same '<base>/<mate>' form."""
    # already suffixed (simulators, the test fixture, the riceTElib benchmark)
    assert canonical_name("read_500_1/1", "1") == "read_500_1/1"
    assert canonical_name("read_500_1/2", "2") == "read_500_1/2"

    # modern Illumina: identical name in both files, mate in the comment field
    illumina = "A00519:123:HXXXX:1:1101:6768:1000"
    assert canonical_name(illumina, "1") == f"{illumina}/1"
    assert canonical_name(illumina, "2") == f"{illumina}/2"

    # SRA / ENA
    assert canonical_name("SRR1234567.1", "1") == "SRR1234567.1/1"

    # a comment field is not part of the name
    assert canonical_name(f"{illumina} 1:N:0:ATCACG", "1") == f"{illumina}/1"

    # the side is authoritative: a stale suffix is corrected, not doubled
    assert canonical_name("read_500_1/1", "2") == "read_500_1/2"


# ---------------------------------------------------------------------------
# mate recovery
# ---------------------------------------------------------------------------


def test_recover_support_mates_without_mate_suffixes(tmp_path: Path):
    """The genomic mate is found even though the FASTQ names carry no /1,/2."""
    r1 = tmp_path / "sample_R1.fastq"
    r2 = tmp_path / "sample_R2.fastq"
    _write_fastq(r1, [("readA", "ACGT" * 10), ("readB", "TTTT" * 10)])
    _write_fastq(r2, [("readA", "GGGG" * 10), ("readB", "CCCC" * 10)])

    reads = ReadLibrary([str(r1), str(r2)], "sample")
    # readA/1 matched a TE; its mate readA/2 did not -> it is a supporting read
    read_repeat = {"readA/1:start:5": ("mPing", "+")}

    out = tmp_path / "support.fq"
    n = recover_support_mates(read_repeat, reads, out)

    assert n == 1, "genomic mate was not recovered from suffix-less FASTQ names"
    assert "GGGG" in out.read_text(), "recovered the wrong mate sequence"


def test_recover_support_mates_still_works_with_suffixes(tmp_path: Path):
    """The existing convention keeps working unchanged."""
    r1 = tmp_path / "sample_R1.fastq"
    r2 = tmp_path / "sample_R2.fastq"
    _write_fastq(r1, [("readA/1", "ACGT" * 10)])
    _write_fastq(r2, [("readA/2", "GGGG" * 10)])

    reads = ReadLibrary([str(r1), str(r2)], "sample")
    out = tmp_path / "support.fq"
    n = recover_support_mates({"readA/1:start:5": ("mPing", "+")}, reads, out)

    assert n == 1
    assert "GGGG" in out.read_text()


def test_mate_that_also_matched_a_te_is_not_a_supporting_read(tmp_path: Path):
    """Both ends matched the TE -> neither is a bracketing genomic mate."""
    r1 = tmp_path / "sample_R1.fastq"
    r2 = tmp_path / "sample_R2.fastq"
    _write_fastq(r1, [("readA", "ACGT" * 10)])
    _write_fastq(r2, [("readA", "GGGG" * 10)])

    reads = ReadLibrary([str(r1), str(r2)], "sample")
    read_repeat = {
        "readA/1:start:5": ("mPing", "+"),
        "readA/2:end:3": ("mPing", "+"),
    }
    out = tmp_path / "support.fq"
    assert recover_support_mates(read_repeat, reads, out) == 0


# ---------------------------------------------------------------------------
# flank pairing (mate anchoring)
# ---------------------------------------------------------------------------


def test_build_flank_pairs_without_mate_suffixes(tmp_path: Path):
    """Junction flanks are anchored to their genomic mate on suffix-less input.

    Mate anchoring is what recovers ambiguous short flanks; losing it was the
    bulk of the damage, since unpaired flanks scatter to equally-scoring loci.
    """
    r1 = tmp_path / "sample_R1.fastq"
    r2 = tmp_path / "sample_R2.fastq"
    _write_fastq(r1, [("readA", "ACGT" * 25)])
    _write_fastq(r2, [("readA", "GGGG" * 25)])

    flank = tmp_path / "flank.fq"
    _write_fastq(flank, [("readA/1:start:5", "ACGT" * 10)])

    reads = ReadLibrary([str(r1), str(r2)], "sample")
    n_pair, n_se, retag = build_flank_pairs(
        [str(flank)],
        {"readA/1:start:5": ("mPing", "+")},
        reads,
        tmp_path / "p1.fq",
        tmp_path / "p2.fq",
        tmp_path / "se.fq",
    )

    assert n_pair == 1, f"flank was not mate-anchored (n_pair={n_pair}, n_se={n_se})"
    assert n_se == 0
    assert retag["readA"] == ("readA/1:start:5", "readA/2")
    assert "GGGG" in (tmp_path / "p2.fq").read_text()


# ---------------------------------------------------------------------------
# TE-library BAMs must carry the mate, whatever the FASTQ looked like
# ---------------------------------------------------------------------------


TE_SEQ = (
    "ATGGCCTTAAGGCCATTGACCTTGGAACCTTGGCATTGCAAGGTTCCAAGGTTAACCGGTTAA"
    "CCGGTTAACCGGTTAAGGCCTTAACCGGTTAAGGCCTTAACCGGTTAAGGCCTTAACCGGTTA"
)


@pytest.mark.parametrize("backend_name", ["minimap2", "bwa", "bwamem2", "bowtie2"])
def test_te_library_bams_carry_mate_suffix(tmp_path: Path, backend_name: str):
    """read_repeat keys must identify the mate even for suffix-less FASTQs.

    The TE stage maps each side to its own BAM, so the side fixes the mate. If a
    backend leaves the name bare, every downstream mate lookup fails and the run
    silently loses all supporting reads.
    """
    import shutil

    import pysam

    from RelocaTE3.aligners import get_aligner

    binary = {"minimap2": "minimap2", "bwa": "bwa", "bwamem2": "bwa-mem2",
              "bowtie2": "bowtie2"}[backend_name]
    if shutil.which(binary) is None:
        pytest.skip(f"{binary} not available")

    te_fa = tmp_path / "te.fa"
    te_fa.write_text(f">mPing\n{TE_SEQ}\n")

    # Illumina-style: identical bare names in both files
    r1 = tmp_path / "sample_R1.fastq"
    r2 = tmp_path / "sample_R2.fastq"
    _write_fastq(r1, [("A00519:1:HXXXX:1:1101:6768:1000", TE_SEQ[:80])])
    _write_fastq(r2, [("A00519:1:HXXXX:1:1101:6768:1000", TE_SEQ[40:120])])

    reads = ReadLibrary([str(r1), str(r2)], "sample")
    outdir = tmp_path / "out"
    outdir.mkdir()

    bams = get_aligner(backend_name, threads=1).map_te_library(
        reads, str(te_fa), str(outdir)
    )
    assert bams, "no TE-library BAMs produced"

    seen = []
    for bam in bams:
        with pysam.AlignmentFile(str(bam), "rb") as fh:
            for rec in fh.fetch(until_eof=True):
                seen.append(rec.query_name)
    assert seen, "no reads aligned to the TE library"
    bare = [n for n in seen if not (n.endswith("/1") or n.endswith("/2"))]
    assert not bare, f"{backend_name} left read names without a mate suffix: {bare[:3]}"


def test_single_end_names_are_left_alone(tmp_path: Path):
    """A single-end library has no mate; do not stamp /1 onto its reads."""
    import shutil

    import pysam

    from RelocaTE3.aligners import get_aligner

    if shutil.which("minimap2") is None:
        pytest.skip("minimap2 not available")

    te_fa = tmp_path / "te.fa"
    te_fa.write_text(f">mPing\n{TE_SEQ}\n")
    r1 = tmp_path / "sample.fastq"
    _write_fastq(r1, [("SRR1234567.1", TE_SEQ[:80])])

    reads = ReadLibrary([str(r1)], "sample")
    outdir = tmp_path / "out"
    outdir.mkdir()
    bams = get_aligner("minimap2", threads=1).map_te_library(
        reads, str(te_fa), str(outdir)
    )
    names = []
    for bam in bams:
        with pysam.AlignmentFile(str(bam), "rb") as fh:
            names.extend(rec.query_name for rec in fh.fetch(until_eof=True))
    assert names and all(n == "SRR1234567.1" for n in names), names
