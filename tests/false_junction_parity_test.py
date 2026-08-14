"""RelocaTE2's false-junction filter, made reachable and applied to one-sided calls.

RelocaTE2 deletes a candidate when the *untrimmed* versions of its junction reads
map straight through the breakpoint (relocaTE_insertionFinder.py:212-221): if the
whole read aligns across the site, the read never crossed a junction and there is
no insertion there. Its rule is

    left_full >= 0.3 * left_total  AND  right_full >= 0.3 * right_total

Two RelocaTE3 divergences are fixed here.

**It was unreachable.** ``_is_false_junction`` was only ever called by the
module-level ``find_insertions`` used by ``pipeline.run_sample``;
``InsertionFinder`` -- what the CLI and the benchmark actually run -- had no
``fullreads_bam`` parameter at all, so the shipped pipeline applied no
false-junction filter. This is the fifth behavioural difference found hiding in
that split.

**It exempted one-sided calls.** ``if left_total == 0 or right_total == 0:
return False`` skipped exactly the population that is only ~45% accurate.
RelocaTE2 has no such exemption: with a side empty, ``0 >= 0.3*0`` is true, so
the rule reduces to the populated side.

The lookup is region-scoped rather than whole-BAM. ``_load_fullread_spans``
materialises every read name in the file, which is 74.6M reads for the shipped
``original_reads`` BAM and does not complete. RelocaTE2 builds its equivalent
from junction reads only (~30k). Fetching a window around each candidate is the
same information for a few hundred small queries.
"""

from __future__ import annotations

import pytest

from RelocaTE3.insertions import _fullread_false_junction
from RelocaTE3.models import Insertion


class FakeRec:
    def __init__(self, name, start, end):
        self.query_name = name
        self.reference_start = start - 1  # pysam is 0-based
        self.reference_end = end
        self.is_unmapped = False


class FakeBam:
    """Stands in for a pysam AlignmentFile, serving one region."""

    def __init__(self, recs):
        self._recs = recs
        self.calls = []

    def fetch(self, chrom, start, end):
        self.calls.append((chrom, start, end))
        return list(self._recs)


def _ins(start, end, left, right, names):
    return Insertion(
        chrom="Chr1", start=start, end=end, te_name="mping", strand="+", tsd="UNK",
        left_junction_reads=left, right_junction_reads=right, read_names=list(names),
    )


def test_two_sided_call_with_full_reads_through_both_breakpoints_is_false():
    bam = FakeBam([FakeRec("a", 900, 1100), FakeRec("b", 900, 1100)])
    ins = _ins(1000, 1002, 1, 1, ["a:end:5", "b:start:5"])
    assert _fullread_false_junction(bam, ins)


def test_two_sided_call_without_full_reads_is_kept():
    bam = FakeBam([])
    ins = _ins(1000, 1002, 1, 1, ["a:end:5", "b:start:5"])
    assert not _fullread_false_junction(bam, ins)


def test_one_sided_call_is_NOT_exempt():
    """The RelocaTE3 bug: one-sided calls skipped the filter entirely."""
    bam = FakeBam([FakeRec("a", 900, 1100)])
    ins = _ins(1000, 1000, 1, 0, ["a:end:5"])
    assert _fullread_false_junction(bam, ins), "one-sided calls must be filtered too"


def test_one_sided_call_without_full_read_support_survives():
    bam = FakeBam([FakeRec("other", 900, 1100)])
    ins = _ins(1000, 1000, 1, 0, ["a:end:5"])
    assert not _fullread_false_junction(bam, ins)


def test_threshold_is_thirty_percent():
    """1 of 4 reads (25%) is below the bar; 2 of 4 (50%) is above it."""
    names = ["a:end:5", "b:end:5", "c:end:5", "d:end:5"]
    one = FakeBam([FakeRec("a", 900, 1100)])
    assert not _fullread_false_junction(one, _ins(1000, 1000, 4, 0, names))
    two = FakeBam([FakeRec("a", 900, 1100), FakeRec("b", 900, 1100)])
    assert _fullread_false_junction(two, _ins(1000, 1000, 4, 0, names))


def test_read_must_span_with_margin():
    """A full read ending at the breakpoint does not span it."""
    bam = FakeBam([FakeRec("a", 900, 1000)])
    assert not _fullread_false_junction(bam, _ins(1000, 1000, 1, 0, ["a:end:5"]))


def test_lookup_is_region_scoped_not_whole_bam():
    """Must query a window, never iterate the entire file."""
    bam = FakeBam([])
    _fullread_false_junction(bam, _ins(5000, 5002, 1, 1, ["a:end:5", "b:start:5"]))
    assert len(bam.calls) >= 1
    chrom, start, end = bam.calls[0]
    assert chrom == "Chr1"
    assert start >= 0 and end - start < 10000, "window must be bounded"


def test_no_bam_means_no_filtering():
    assert not _fullread_false_junction(None, _ins(1000, 1002, 1, 1, ["a:end:5"]))


def test_call_with_no_junction_reads_is_untouched():
    bam = FakeBam([FakeRec("a", 900, 1100)])
    assert not _fullread_false_junction(bam, _ins(1000, 1002, 0, 0, []))


def test_mate_suffix_is_stripped_when_matching_full_reads():
    """Flank names carry /1 or /2; the untrimmed reads keep the mate in the flag.

    Stripping only the junction tag left `/2` on the key and matched nothing —
    the filter ran and silently did nothing. 0 of 20,000 sampled names in the
    untrimmed BAM carry a mate suffix.
    """
    bam = FakeBam([FakeRec("readX", 900, 1100)])
    ins = _ins(1000, 1000, 1, 0, ["readX/2:end:5"])
    assert _fullread_false_junction(bam, ins)
