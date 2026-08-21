"""RelocaTE2's read-admission gate, which RelocaTE3 had dropped.

``InsertionFinder._passes_quality`` carried the comment "Minimap2-adapted
quality filter (replaces BWA XT/X1/XM/XO logic)" — written when RelocaTE3
targeted minimap2, which does not emit those tags. RelocaTE3 now defaults to
``bwa aln`` for genome placement, so the tags are present on every record and
were simply being ignored.

RelocaTE2's gate (relocaTE_insertionFinder.py:1521-1558):

    properly paired:
        junction read   -> XM <= mismatch AND XO <= 1
        supporting read -> XM <= mismatch
    not properly paired:
        either kind     -> XT == 'U' AND XM <= mismatch AND X1 <= 3

The `XT == 'U'` requirement is the substantive one: a read that is not properly
paired must be *uniquely* mapped. On the mPing benchmark 225 of 639 sampled
reads are not properly paired and 60 carry `XT:A:R` (repeat) — RelocaTE2 rejects
those as evidence, RelocaTE3 accepted them, and they are a prime suspect for
RelocaTE3's one-sided calls being 38% accurate against RelocaTE2's 100%.

The gates apply only when the tags exist, so minimap2/bowtie2 runs are
unaffected.
"""

from __future__ import annotations

import pytest

from RelocaTE3.insertions import InsertionFinder


class FakeRead:
    """Minimal stand-in for a pysam AlignedSegment."""

    def __init__(self, name="r1:start:5", mapq=60, proper=True, tags=None, nm=0):
        self.query_name = name
        self.mapping_quality = mapq
        self.is_proper_pair = proper
        self.cigartuples = [(0, 100)]
        self._tags = dict(tags or {})
        if nm is not None:
            self._tags.setdefault("NM", nm)

    def has_tag(self, tag):
        return tag in self._tags

    def get_tag(self, tag):
        return self._tags[tag]


@pytest.fixture
def finder():
    return InsertionFinder(mismatch_allow=2, min_mapq=1)


# --- not properly paired: uniqueness is required -------------------------


def test_unpaired_multimapping_read_is_rejected(finder):
    """XT:A:R means bwa placed it in a repeat — RelocaTE2 will not use it."""
    assert not finder._passes_quality(
        FakeRead(proper=False, tags={"XT": "R", "X1": 0})
    )


def test_unpaired_unique_read_is_accepted(finder):
    assert finder._passes_quality(FakeRead(proper=False, tags={"XT": "U", "X1": 0}))


def test_unpaired_read_with_many_suboptimal_hits_is_rejected(finder):
    """RelocaTE2 caps X1 (suboptimal hit count) at 3 for unpaired reads."""
    assert not finder._passes_quality(
        FakeRead(proper=False, tags={"XT": "U", "X1": 4})
    )
    assert finder._passes_quality(FakeRead(proper=False, tags={"XT": "U", "X1": 3}))


# --- properly paired: uniqueness NOT required ----------------------------


def test_paired_multimapping_read_is_still_accepted(finder):
    """RelocaTE2 only demands uniqueness when the pair is not proper."""
    assert finder._passes_quality(FakeRead(proper=True, tags={"XT": "R", "X1": 9}))


def test_paired_junction_read_with_two_gaps_is_rejected(finder):
    """XO <= 1 applies to junction reads in a proper pair."""
    assert not finder._passes_quality(
        FakeRead(name="r:start:5", proper=True, tags={"XT": "U", "XO": 2})
    )
    assert finder._passes_quality(
        FakeRead(name="r:start:5", proper=True, tags={"XT": "U", "XO": 1})
    )


def test_paired_supporting_read_is_not_gated_on_gaps(finder):
    """The XO cap is junction-read-only in RelocaTE2."""
    assert finder._passes_quality(
        FakeRead(name="plain_read", proper=True, tags={"XT": "U", "XO": 3})
    )


# --- graceful when the aligner does not emit the tags --------------------


def test_records_without_bwa_tags_are_unaffected(finder):
    """minimap2/bowtie2 emit no XT/X1/XO — behaviour must not change."""
    assert finder._passes_quality(FakeRead(proper=False, tags={}))
    assert finder._passes_quality(FakeRead(proper=True, tags={}))


def test_existing_mismatch_and_mapq_gates_still_apply(finder):
    """``min_mapq`` still gates when a caller opts into it (fixture uses 1)."""
    assert not finder._passes_quality(FakeRead(mapq=0, tags={"XT": "U"}))
    assert not finder._passes_quality(FakeRead(nm=5, tags={"XT": "U"}))


def test_default_admits_mapq_zero_reads_as_relocate2_does():
    """RelocaTE2 has no MAPQ admission gate; the RelocaTE3 default must not either.

    ``relocaTE_insertionFinder.py:1521-1558`` admits a read on XM/XO/XT/X1 only.
    MAPQ is used at :1523,1539 solely to *record* the read as low quality, and
    the call is discarded later (:226-241) only if it rests entirely on such
    reads. RelocaTE3 defaulted ``min_mapq`` to 1, discarding MAPQ-0 reads
    outright. That cost real calls: on the mPing benchmark RelocaTE2 frequently
    resolves the second junction from a single MAPQ-0 read, which RelocaTE3
    never saw, leaving a one-sided call that the both-junctions gate then
    dropped. Restoring the default recovered 4 true calls at cov30x_rep1 with
    no new false positives.
    """
    default = InsertionFinder(mismatch_allow=2)
    assert default.min_mapq == 0
    assert default._passes_quality(FakeRead(mapq=0, tags={"XT": "U"}))
    # ... and it is still recognised as low quality, so it cannot validate a
    # call on its own.
    assert InsertionFinder._is_low_quality(FakeRead(mapq=0, tags={"XT": "U"}))


# ---------------------------------------------------------------------------
# RelocaTE2 deletes calls validated only by low-quality junction reads
# ---------------------------------------------------------------------------
#
# relocaTE_insertionFinder.py:226-241 — after building a candidate, RelocaTE2
# counts how many of its junction reads are NOT in teLowQualityReads and drops
# the call when:
#     total_real == 1 and total_valid == 0   (lone low-quality junction read)
#     left_valid == 0 and right_valid == 0   (no high-quality junction read at all)
#
# RelocaTE3 had no equivalent: a call resting entirely on MAPQ<29 or
# improperly-paired reads was reported as though it were solid.

from RelocaTE3.insertions import _call_validated_by_high_quality
from RelocaTE3.models import JunctionObservation


def _j(side, pos, low):
    return JunctionObservation("r", side, pos, "+", "mping", "5", "ACGT", pos, pos + 50, low)


def test_call_with_a_high_quality_junction_read_is_kept():
    assert _call_validated_by_high_quality([_j("left", 100, False)], [])


def test_call_resting_only_on_low_quality_reads_is_dropped():
    """No high-quality junction read on either side."""
    assert not _call_validated_by_high_quality(
        [_j("left", 100, True), _j("left", 100, True)], [_j("right", 98, True)]
    )


def test_mixed_quality_call_is_kept():
    """One good read is enough — RelocaTE2 requires validity, not purity."""
    assert _call_validated_by_high_quality(
        [_j("left", 100, True)], [_j("right", 98, False)]
    )


def test_lone_low_quality_junction_read_is_dropped():
    """total_real == 1 and total_valid == 0."""
    assert not _call_validated_by_high_quality([_j("left", 100, True)], [])


def test_call_with_no_junction_reads_is_not_gated_here():
    """Support-only clusters take a different path; do not reject them here."""
    assert _call_validated_by_high_quality([], [])
