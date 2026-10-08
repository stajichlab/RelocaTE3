"""Small-table error-audit coordinate and one-to-one safeguards."""
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "error_audit", Path(__file__).parents[1] / "scripts/audit_benchmark_errors.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def row(position, family="Os_1#LINE", chrom="Chr1", **extra):
    return dict(chrom=chrom, position=str(position), te_family=family, **extra)


def test_interval_endpoints_and_tolerance():
    truth = row(100, truth_match_start="80", truth_match_end="100")
    assert [audit.truth_distance(truth, row(p)) for p in (69, 70, 80, 90, 100, 110, 111)] == [11, 10, 0, 0, 0, 10, 11]


def test_legacy_fallback():
    assert audit.truth_distance(row(100), row(85)) == 15
    assert audit.truth_distance(row(100, truth_match_start="", truth_match_end=""), row(85)) == 15


def test_near_interval_includes_wrong_family_not_wrong_chromosome():
    truth = row(100, truth_match_start="80", truth_match_end="100")
    candidates = [row(70, family="other"), row(111), row(80, chrom="Chr2")]
    actual = audit.near_interval(truth, candidates)
    assert len(actual) == 1 and actual[0]["te_family"] == "other"
    assert "interval_distance_bp" not in candidates[0]


def test_fp_pairing_one_to_one_and_family_normalization():
    assert audit.pair_fps([row(100), row(100)], [row(100, "os1")]) == ({0}, {0})
    assert audit.pair_fps([row(100)], [row(100, "other")]) == (set(), set())


def test_fp_pairing_tolerance_and_empty():
    assert audit.pair_fps([row(100)], [row(110)]) == ({0}, {0})
    assert audit.pair_fps([row(100)], [row(111)]) == (set(), set())
    assert audit.pair_fps([], []) == (set(), set())
