"""An implausibly long inferred TSD means "no TSD", not a TSD.

Target-site duplications are short: a few bp for most elements, up to ~20 for
the longest. When the depth estimator returns a much larger span it is not
measuring a duplication, it is measuring something else -- and RelocaTE3 was
reporting the genomic sequence over that span as though it were the TSD.

Measured on riceTElib (9 samples, blat/bwa-aln), splitting truth sites by
whether they actually have a TSD:

    site class            R3 detects   R3 reports a TSD   longest true TSD
    has a real TSD             1972    all <= 12 bp                  20 bp
    TSD-less (TSD=NONE)         255    140, of which 114 > 20 bp

So every fabricated TSD over the plausible length is on an element that has no
TSD at all, and no correct call is anywhere near the cap. Reporting UNK instead
loses nothing and stops the caller from inventing sequence.

Note this does not change *detection*: the insertion is still called, only the
TSD field changes. RelocaTE3 detects 255 of these TSD-less sites to RelocaTE2's
10, so it is the caller that has to get this right.
"""

from __future__ import annotations

import pytest

from RelocaTE3.insertions import MAX_PLAUSIBLE_TSD, _resolve_tsd
from RelocaTE3.models import JunctionObservation


def _obs(seq: str, side: str = "right") -> JunctionObservation:
    return JunctionObservation("r", side, 100, "+", "mping", "5", seq, 100, 140)


def test_cap_is_at_least_the_longest_real_tsd_seen():
    """The riceTElib truth tops out at 20 bp; the cap must not cut real TSDs."""
    assert MAX_PLAUSIBLE_TSD >= 20


def test_short_tsd_is_captured_from_the_read():
    tsd = _resolve_tsd([], [_obs("TTAGGGCCCAAA")], "Chr1", 100, 102, 3, None)
    assert tsd == "TTA"


def test_most_supported_junction_tsd_is_selected():
    """Three GCA junctions must outvote one erroneous TCA junction."""
    left = [_obs("AAAAAGCA", side="left")]
    right = [_obs("TCAGGG"), _obs("GCATTT"), _obs("GCACCC")]

    assert _resolve_tsd(left, right, "Chr1", 100, 102, 3, None) == "GCA"


def test_tied_junction_tsds_keep_the_first_right_capture():
    """A tie retains the deterministic right-before-left legacy behavior."""
    left = [_obs("AAAAAGCA", side="left")]
    right = [_obs("TCAGGG")]

    assert _resolve_tsd(left, right, "Chr1", 100, 102, 3, None) == "TCA"


def test_tsd_at_the_cap_is_still_reported():
    length = MAX_PLAUSIBLE_TSD
    seq = "A" * length + "GGGG"
    assert _resolve_tsd([], [_obs(seq)], "Chr1", 100, 100 + length - 1, length, None) == "A" * length


def test_implausibly_long_inference_reports_unk():
    """An 87 bp 'TSD' is not a TSD -- observed on TSD-less elements."""
    length = MAX_PLAUSIBLE_TSD + 1
    seq = "ACGT" * 40
    assert _resolve_tsd([], [_obs(seq)], "Chr1", 100, 100 + length - 1, length, None) == "UNK"


def test_cap_applies_to_the_genome_fallback_too():
    """The fallback path must not smuggle a long span back in."""
    class FakeGenome:
        def fetch(self, chrom, start, end):
            return "A" * (end - start)

    length = MAX_PLAUSIBLE_TSD + 5
    got = _resolve_tsd([], [], "Chr1", 100, 100 + length - 1, length, FakeGenome())
    assert got == "UNK"


@pytest.mark.parametrize("length", [0, -1])
def test_non_positive_length_still_unk(length):
    assert _resolve_tsd([], [_obs("ACGTACGT")], "Chr1", 100, 100, length, None) == "UNK"
