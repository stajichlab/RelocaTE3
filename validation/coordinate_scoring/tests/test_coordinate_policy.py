"""Synthetic matching-contract tests, independent of benchmark outcomes."""

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scoring.score_calls import _interval_distance, _truth_interval, score


def write(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def evaluate(tmp_path, positions, length=20, families=None, policy="tsd-interval",
             call_tsd="UNK", window=10):
    truth = tmp_path / "truth.tsv"
    calls = tmp_path / "calls.tsv"
    write(truth, [dict(event_id="E1", chrom="Chr1", position=100, tsd_length=length,
                       tsd="A" * length if length else "NONE", te_family="TE1",
                       biological_class="somatic_insertion")])
    write(calls, [dict(chrom="Chr1", position=p, te_family=(families or ["TE1"] * len(positions))[i],
                       tsd=call_tsd, strand="-", status="heterozygous", caller="R3", sample="S")
                  for i, p in enumerate(positions)])
    return score(truth, calls, "S", "R3", window, policy)


@pytest.mark.parametrize("position,expected", [(81, 1), (100, 1), (71, 1), (110, 1), (70, 0), (111, 0)])
def test_interval_and_unchanged_tolerance(tmp_path, position, expected):
    _, matches, _, precision = evaluate(tmp_path, [position])
    assert precision["matched_calls"] == expected
    assert matches[0]["position"] == "100"  # original anchor preserved
    assert matches[0]["truth_match_start"] == 81
    assert matches[0]["truth_match_end"] == 100


def test_long_tsd_rescue_not_genotype_or_sequence_credit(tmp_path):
    _, matches, _, precision = evaluate(tmp_path, [81])
    assert precision["matched_calls"] == 1
    assert matches[0]["distance_bp"] == 0
    assert matches[0]["anchor_distance_bp"] == 19
    assert matches[0]["tsd_exact"] == "0"
    assert matches[0]["status_correct"] == "0"
    assert evaluate(tmp_path, [81], policy="legacy")[3]["matched_calls"] == 0


def test_duplicate_calls_are_not_both_credited(tmp_path):
    _, matches, fps, precision = evaluate(tmp_path, [81, 100])
    assert precision["matched_calls"] == 1
    assert len(fps) == 1
    # Original anchor distance is the explicit second tie-breaker.
    assert matches[0]["call_position"] == "100"


def test_wrong_family_stays_false_positive(tmp_path):
    assert evaluate(tmp_path, [81], families=["TE2"])[3]["matched_calls"] == 0


def test_no_tsd_is_a_point_not_length_of_sentinel(tmp_path):
    assert _truth_interval({"position": 100, "tsd": "NONE", "tsd_length": 0}, "tsd-interval") == (100, 100)
    assert evaluate(tmp_path, [89], length=0)[3]["matched_calls"] == 0
    assert evaluate(tmp_path, [90], length=0)[3]["matched_calls"] == 1


def test_old_truth_length_inference_is_conservative():
    assert _truth_interval({"position": 100, "tsd": "ACGT"}, "tsd-interval") == (97, 100)
    for tsd in ("UNK", "UKN", "NONE", "supporting_junction", ""):
        assert _truth_interval({"position": 100, "tsd": tsd}, "tsd-interval") == (100, 100)


@pytest.mark.parametrize("length", [-1, 101, "bad"])
def test_bad_truth_length_fails(length):
    with pytest.raises(ValueError):
        _truth_interval({"position": 100, "tsd_length": length}, "tsd-interval")


def test_one_call_cannot_match_two_truth_events(tmp_path):
    truth, calls = tmp_path / "truth.tsv", tmp_path / "calls.tsv"
    write(truth, [dict(event_id=e, chrom="Chr1", position=p, tsd_length=20,
                      tsd="A"*20, te_family="TE1", biological_class="homozygous")
                  for e, p in (("E1", 100), ("E2", 105))])
    write(calls, [dict(chrom="Chr1", position=95, te_family="TE1", tsd="UNK",
                      strand="+", status="homozygous")])
    _, matches, _, precision = score(truth, calls, "S", "R3", 10, "tsd-interval")
    assert sum(r["matched"] == "1" for r in matches) == 1
    assert precision["matched_calls"] == 1


def test_negative_window_fails(tmp_path):
    with pytest.raises(ValueError):
        evaluate(tmp_path, [100], window=-1)


def test_interval_distance():
    assert [_interval_distance(p, (81, 100)) for p in (80, 81, 90, 100, 101)] == [1, 0, 0, 0, 1]
