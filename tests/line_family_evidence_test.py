"""Characterize real LINE ambiguity without requiring original benchmark files.

Archive selections describe the historical behavior; they are not assertions
that its primary family is biologically correct.
"""
import importlib.util
import json
from pathlib import Path

import pytest
from RelocaTE3.librelocate import RelocaTE

root = Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location(
    "line_family_analysis", root / "scripts/analyze_line_family_trace.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


@pytest.fixture(scope="module")
def evidence():
    return json.loads((root / "tests/data/line_family/TE000172.json").read_text())


@pytest.fixture(scope="module")
def comparison(evidence):
    return analysis.analyze(evidence)


def test_original_hit_sets_and_archived_mappings_reproduce(comparison):
    rows, summary = comparison
    assert len(rows) == 11
    assert summary["eligible_alignment_signatures_equal"]
    assert summary["R2_only_eligible_alignments"] == summary["R3_only_eligible_alignments"] == 0
    assert summary["R2_policy_change_reproduces_R3"] == 11
    assert summary["R2_votes"] == {"Os3328#LINE/unknown": 10, "Os0596#LINE/unknown": 1}
    assert summary["R3_votes"] == {"Os3328#LINE/unknown": 2, "Os0596#LINE/unknown": 9}


def test_switched_votes_have_equal_primary_alignment_scores(comparison):
    rows, summary = comparison
    changed = [r for r in rows if r["changed_family"]]
    assert len(changed) == summary["ambiguous_best_family_reads"] == 8
    for row in changed:
        top = json.loads(row["top_psl_records"])
        assert {r["tName"] for r in top} == {"Os3328#LINE/unknown", "Os0596#LINE/unknown"}
        assert len({(r["boundary"], r["match"], r["mismatch"], r["start"], r["end"])
                    for r in top}) == 1
    assert summary["single_best_family_votes"] == {"Os3328#LINE/unknown": 2, "Os0596#LINE/unknown": 1}


def test_deterministic_selection_survives_order_changes(evidence):
    hits = analysis.select_psl(evidence)
    for records in hits.values():
        valid = [r for r in records if r["admitted"]]
        selected = min(valid, key=RelocaTE._match_rank)
        assert selected["tName"] == min(reversed(valid), key=RelocaTE._match_rank)["tName"]
        tied = [r for r in valid if (r["boundary"], r["match"]) ==
                (selected["boundary"], selected["match"])]
        if len({r["tName"] for r in tied}) > 1:
            # R2's first-on-tie policy can change family with the same hit set.
            assert tied[0]["tName"] != tied[-1]["tName"]
