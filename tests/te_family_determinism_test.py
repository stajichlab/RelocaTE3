"""TE-family assignment must not depend on alignment record order.

A read that matches several TE families equally well used to be assigned to
whichever alignment the aligner happened to emit first: ``_is_better`` returned
False on an exact tie, so the first record won. Aligner output order is not
stable across threads for multi-mapping reads, so two identical runs could
label the same call with different families (observed on the Chr3 2 Mb fixture:
``RIRE3`` vs ``mGing`` at Chr3:672695, identical positions and read counts).

These tests pin a deterministic total order so the winner is reproducible.
"""

from __future__ import annotations

import os

from RelocaTE3.librelocate import RelocaTE


def _rec(tname, boundary, match, tstart=0, qstart=0, qend=99):
    """A minimal TE-alignment record of the shape ``_parse_te_bam`` builds."""
    return {
        "tName": tname,
        "boundary": boundary,
        "match": match,
        "tStart": tstart,
        "tEnd": tstart + match,
        "start": qstart,
        "end": qend,
        "len": 100,
        "tLen": 500,
        "mismatch": 0,
        "strand": "+",
        "seq": "A" * 100,
        "qual": "I" * 100,
    }


def _select(records):
    """Reproduce the winner-selection loop used by _parse_te_bam / _merge."""
    best = None
    for rec in records:
        if best is None or RelocaTE._is_better(rec, best):
            best = rec
    return best


def test_equal_scoring_hits_resolve_the_same_way_in_either_order():
    """The core bug: identical scores, so emission order decided the family."""
    a = _rec("RIRE3", boundary=2, match=40)
    b = _rec("mGing", boundary=2, match=40)

    forward = _select([a, b])
    reverse = _select([b, a])

    assert forward["tName"] == reverse["tName"], (
        "TE family depends on record order: "
        f"{forward['tName']} vs {reverse['tName']}"
    )


def test_better_scores_still_win_over_the_tie_break():
    """Determinism must not override the real ranking signals."""
    strong = _rec("zzz_family", boundary=4, match=40)
    weak = _rec("aaa_family", boundary=2, match=40)
    assert _select([weak, strong])["tName"] == "zzz_family"
    assert _select([strong, weak])["tName"] == "zzz_family"

    longer = _rec("zzz_family", boundary=2, match=80)
    shorter = _rec("aaa_family", boundary=2, match=40)
    assert _select([shorter, longer])["match"] == 80
    assert _select([longer, shorter])["match"] == 80


def test_same_family_multiple_loci_resolve_deterministically():
    """Two hits to the SAME family at different target offsets also tie."""
    first = _rec("mPing", boundary=2, match=40, tstart=10)
    second = _rec("mPing", boundary=2, match=40, tstart=900)
    assert _select([first, second])["tStart"] == _select([second, first])["tStart"]


def test_is_better_is_antisymmetric():
    """A strict order: two distinct records cannot each be better than the other."""
    a = _rec("RIRE3", boundary=2, match=40)
    b = _rec("mGing", boundary=2, match=40)
    assert not (RelocaTE._is_better(a, b) and RelocaTE._is_better(b, a))
    # and a record is never better than itself (irreflexive)
    assert not RelocaTE._is_better(a, dict(a))


def test_selection_is_stable_across_many_permutations():
    """Any permutation of the same candidate set yields the same winner."""
    import itertools

    records = [
        _rec("RIRE3", boundary=2, match=40),
        _rec("mGing", boundary=2, match=40),
        _rec("mPing", boundary=2, match=40, tstart=7),
        _rec("Tos17", boundary=2, match=40),
    ]
    winners = {
        (_select(list(p))["tName"], _select(list(p))["tStart"])
        for p in itertools.permutations(records)
    }
    assert len(winners) == 1, f"non-deterministic winner across orderings: {winners}"


# ---------------------------------------------------------------------------
# Cluster-level family vote (second source, downstream of read assignment)
# ---------------------------------------------------------------------------


def test_majority_te_name_breaks_ties_deterministically():
    """A tied family vote must not depend on set iteration order.

    ``max(set(names), key=names.count)`` returns whichever tied name the set
    happens to yield first, and CPython randomises string hashing per process,
    so two runs of the same code on the same data could disagree. Ties are
    resolved by name instead.
    """
    from RelocaTE3.insertions import _majority_te_name

    tied = ["mGing"] * 4 + ["RIRE3"] * 4
    assert _majority_te_name(tied) == "RIRE3"  # lexicographically smallest
    assert _majority_te_name(list(reversed(tied))) == "RIRE3"

    # a real majority still wins regardless of ordering
    clear = ["zzz"] * 5 + ["aaa"] * 2
    assert _majority_te_name(clear) == "zzz"
    assert _majority_te_name(list(reversed(clear))) == "zzz"

    # reads with no family assignment are ignored; empty -> NA
    assert _majority_te_name(["NA", "NA"]) == "NA"
    assert _majority_te_name([]) == "NA"
    assert _majority_te_name(["NA", "mPing", "NA"]) == "mPing"


def test_family_evidence_distinguishes_unique_dominant_and_ambiguous_votes():
    """The primary label stays singular while competing evidence remains visible."""
    from RelocaTE3.insertions import _te_family_evidence

    unique = _te_family_evidence(["mPing", "mPing", "NA"])
    assert unique.primary == "mPing"
    assert unique.support == {"mPing": 2}
    assert unique.confidence == 1.0
    assert unique.status == "unique"

    dominant = _te_family_evidence(["mPing", "mPing", "RIRE3"])
    assert dominant.primary == "mPing"
    assert dominant.support == {"mPing": 2, "RIRE3": 1}
    assert dominant.confidence == 2 / 3
    assert dominant.status == "dominant"

    tied = _te_family_evidence(["mGing", "RIRE3"])
    assert tied.primary == "RIRE3"
    assert tied.support == {"RIRE3": 1, "mGing": 1}
    assert tied.confidence == 0.5
    assert tied.status == "ambiguous"

    plurality = _te_family_evidence(["mPing", "mPing", "RIRE3", "mGing"])
    assert plurality.primary == "mPing"
    assert plurality.confidence == 0.5
    assert plurality.status == "ambiguous"

    missing = _te_family_evidence(["NA", ""])
    assert missing.primary == "NA"
    assert missing.support == {}
    assert missing.confidence == 0.0
    assert missing.status == "unassigned"


def test_majority_te_name_is_stable_across_hash_seeds():
    """The actual failure mode: identical input, different PYTHONHASHSEED."""
    import subprocess
    import sys

    prog = (
        "from RelocaTE3.insertions import _majority_te_name;"
        "print(_majority_te_name(['mGing']*4 + ['RIRE3']*4))"
    )
    results = set()
    for seed in ("0", "1", "2", "3", "4", "5"):
        proc = subprocess.run(
            [sys.executable, "-c", prog],
            capture_output=True,
            text=True,
            env={"PYTHONHASHSEED": seed, "PATH": os.environ.get("PATH", "")},
        )
        assert proc.returncode == 0, proc.stderr
        results.add(proc.stdout.strip())
    assert len(results) == 1, f"family vote varies with hash seed: {results}"
