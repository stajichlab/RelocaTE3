"""Unit tests for the TE-trimming logic in RelocaTE3.trim."""

from __future__ import annotations

from RelocaTE3.models import JunctionType, TEReadAlignment
from RelocaTE3.trim import reverse_complement, trim_read


def _aln(**kw) -> TEReadAlignment:
    """Build a TEReadAlignment with sensible defaults for testing."""
    defaults = dict(
        read_name="r1",
        read_length=100,
        qstart=0,
        qend=63,
        te_name="mPing",
        te_length=430,
        tstart=0,
        tend=63,
        match=64,
        mismatch=0,
        strand="+",
        boundary=2,
        seq="A" * 36 + "C" * 64,
        qual="I" * 100,
    )
    defaults.update(kw)
    return TEReadAlignment(**defaults)


def test_reverse_complement():
    assert reverse_complement("ACGTN") == "NACGT"


def test_five_prime_plus_strand_trims_left_flank():
    # read = [36bp genomic flank][64bp into TE 5' end]; aligned region is the TE part
    aln = _aln(qstart=36, qend=99, tstart=0, tend=63, seq="G" * 36 + "C" * 64)
    tr = trim_read(aln, len_cut_match=10, len_cut_trim=10)
    assert tr is not None
    assert tr.junction_type is JunctionType.FIVE_PRIME
    assert tr.read_name == "r1:end:5"
    assert tr.seq == "G" * 36  # the genomic flank, left side


def test_three_prime_plus_strand_trims_right_flank():
    # read = [64bp into TE 3' end][36bp genomic flank]
    aln = _aln(qstart=0, qend=63, tstart=366, tend=429, seq="C" * 64 + "G" * 36)
    tr = trim_read(aln, len_cut_match=10, len_cut_trim=10)
    assert tr is not None
    assert tr.junction_type is JunctionType.THREE_PRIME
    assert tr.read_name == "r1:start:3"
    assert tr.seq == "G" * 36


def test_middle_read_is_supporting_untrimmed():
    # whole read lies inside the TE
    aln = _aln(qstart=0, qend=99, read_length=100, tstart=100, tend=199, seq="C" * 100)
    tr = trim_read(aln)
    assert tr is not None
    assert tr.junction_type is JunctionType.MIDDLE
    assert tr.read_name == "r1:middle"
    assert tr.seq == "C" * 100


def test_short_flank_rejected():
    # only a 3bp flank, below len_cut_trim
    aln = _aln(qstart=3, qend=99, tstart=0, tend=96, seq="G" * 3 + "C" * 97)
    assert trim_read(aln, len_cut_trim=10) is None


def test_too_many_mismatches_rejected():
    aln = _aln(
        qstart=36, qend=99, tstart=0, tend=63, mismatch=5, seq="G" * 36 + "C" * 64
    )
    assert trim_read(aln, mismatch_allowance=0) is None
