"""Long-read reference space: the linear FASTA implementation."""

from __future__ import annotations

from pathlib import Path

import pytest

from RelocaTE3.longread.model import Locus, Reference
from RelocaTE3.longread.reference import (
    LinearFastaReferenceSpace, reverse_complement,
)  # fmt: skip


def _write_fasta(path: Path) -> Path:
    # a: 10 bp, last five lowercase (soft-masked); b: 5 bp.
    path.write_text(">a\nACGTAcgtac\n>b\nTTGCA\n")
    return path


@pytest.fixture
def space(tmp_path):
    return LinearFastaReferenceSpace(
        Reference("ref", str(_write_fasta(tmp_path / "g.fa")))
    )


def test_sequences_in_file_order(space):
    assert space.sequences() == [("a", 10), ("b", 5)]


def test_fetch_plus_strand_returns_slice(space):
    assert space.fetch(Locus("ref", "a", 2, 5, "+")) == "GTA"


def test_fetch_minus_strand_returns_reverse_complement(space):
    assert space.fetch(Locus("ref", "a", 2, 5, "-")) == "TAC"


def test_fetch_uppercases_soft_masked_bases(space):
    # pysam returns bases as stored; the reference space uppercases them.
    assert space.fetch(Locus("ref", "a", 4, 8, "+")) == "ACGT"
    assert space.fetch(Locus("ref", "a", 4, 8, "-")) == "ACGT"
    assert space.fetch(Locus("ref", "a", 5, 10, "-")) == "GTACG"


def test_fetch_empty_insertion_point(space):
    assert space.fetch(Locus("ref", "b", 3, 3, "+")) == ""


def test_fetch_other_reference_raises_value_error(space):
    with pytest.raises(ValueError, match="other"):
        space.fetch(Locus("other", "a", 0, 1, "+"))


def test_fetch_unknown_sequence_raises_key_error(space):
    with pytest.raises(KeyError, match="nope"):
        space.fetch(Locus("ref", "nope", 0, 1, "+"))


def test_missing_fasta_names_the_path(tmp_path):
    missing = tmp_path / "absent.fa"
    with pytest.raises(FileNotFoundError, match="absent.fa"):
        LinearFastaReferenceSpace(Reference("ref", str(missing)))


def test_fai_created_when_missing(tmp_path):
    fasta = _write_fasta(tmp_path / "g.fa")
    fai = tmp_path / "g.fa.fai"
    assert not fai.exists()
    LinearFastaReferenceSpace(Reference("ref", str(fasta)))
    assert fai.is_file()


def test_reference_attribute(space):
    assert space.reference == Reference("ref", space.reference.fasta)


def test_reverse_complement_handles_lowercase_and_n():
    assert reverse_complement("ACGTN") == "NACGT"
    assert reverse_complement("acgtn") == "nacgt"
    assert reverse_complement("AcGtN") == "NaCgT"
    assert reverse_complement("") == ""


def test_reverse_complement_iupac_codes():
    pairs = ["RY", "KM", "BV", "DH"]
    for a, b in pairs + [p.lower() for p in pairs]:
        assert reverse_complement(a) == b
        assert reverse_complement(b) == a
    for code in "NSWnsw":
        assert reverse_complement(code) == code
    assert reverse_complement("ARYKMSWBVDHNt") == "aNDHBVWSKMRYT"
