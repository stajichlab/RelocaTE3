"""Real local-cluster characterization and rejected pairing experiment.

These do not require the benchmark or any external aligner. Archived snapshots
retain duplicate behavior as a known limitation, not a desired biological result.
The exclusive-boundary experiment lost 30 detections in the full replay and
must not become production behavior based only on these two local fixtures.
"""

import json
from collections import Counter
from pathlib import Path

import pytest

from RelocaTE3 import insertions as ins
from RelocaTE3.models import JunctionObservation


@pytest.fixture(params=["cov30x_rep1", "cov30x_rep2"])
def locus(request):
    data = json.loads(
        (
            Path(__file__).parent / "data/duplicate_loci" / f"{request.param}.json"
        ).read_text()
    )
    cluster = ins._Cluster(data["chrom"])
    cluster.lo, cluster.hi = data["cluster_lo"], data["cluster_hi"]
    cluster.junctions = [JunctionObservation(**j) for j in data["junctions"]]
    cluster.support = [tuple(r) for r in data["support"]]
    return data, cluster


def breakpoints(cluster):
    return tuple(
        ins._group_by_position([j for j in cluster.junctions if j.side == side])
        for side in ("left", "right")
    )


def test_real_cluster_has_three_breakpoints_and_complete_boundaries(locus):
    data, cluster = locus
    p = data["truth_position"]
    expected = (21, 26) if data["sample"] == "cov30x_rep1" else (32, 32)
    counts = Counter((j.side, j.position) for j in cluster.junctions)
    assert counts == {
        ("left", p): expected[0],
        ("right", p): 1,
        ("right", p + 1): expected[1],
    }
    lo, hi = data["query_0based_half_open"]
    assert cluster.lo > lo + ins.RANGE_ALLOWANCE + 1
    assert cluster.hi < hi - ins.RANGE_ALLOWANCE
    assert all(j.te_name != "NA" for j in cluster.junctions)


@pytest.fixture
def exclusive_pairing(monkeypatch):
    """Reproduce the rejected coordinate experiment only within its tests."""
    original = ins._pair_breakpoints

    def experiment(left, right):
        return [
            (p - 1 if p is not None else None, q)
            for p, q in original({p + 1: reads for p, reads in left.items()}, right)
        ]

    monkeypatch.setattr(ins, "_pair_breakpoints", experiment)


def test_restored_pairing_retains_historical_candidates(locus):
    data, cluster = locus
    p = data["truth_position"]
    assert ins._pair_breakpoints(*breakpoints(cluster)) == [(p, p), (None, p + 1)]
    calls = ins._call_insertions(cluster, None, data["read_repeat"])
    assert sorted(c.start for c in calls) == [p, p + 1]


def test_experiment_removes_archived_duplicate(locus, exclusive_pairing):
    data, cluster = locus
    raw = ins._call_insertions(cluster, None, data["read_repeat"])
    assert len(raw) == 1
    assert raw[0].start == data["truth_position"]
    assert (raw[0].left_junction_reads, raw[0].right_junction_reads) == (0, 1)
    pooled = ins._consolidate_same_start(raw, cluster, data["read_repeat"])
    # Full-read filtering of the changed candidate is validated by HPC replay,
    # not inferred from old candidates' stored filter outcomes.
    qualified = [
        c
        for c in pooled
        if ins._call_validated_by_high_quality(*ins._candidate_junctions(c, cluster))
    ]
    edges = {
        side: {int(k): v for k, v in values.items()}
        for side, values in data["reference_edges"].items()
    }
    chosen = ins._arbitrate_cluster(qualified, edges)
    emitted = [ins._as_supporting_junction(c, cluster) or c for c in chosen]
    assert len(emitted) == (0 if data["sample"] == "cov30x_rep1" else 1)
    if emitted:
        assert (emitted[0].start, emitted[0].end, emitted[0].tsd) == (
            9006303,
            9006305,
            "supporting_junction",
        )
    assert len(data["emitted"]) == 2
    # This is NOT reuse of the same junction read under two call names.
    assert set(data["emitted"][0]["read_names"]).isdisjoint(
        data["emitted"][1]["read_names"]
    )


def test_experiment_preserves_inclusive_output_coordinates(locus, exclusive_pairing):
    data, cluster = locus
    left, right = breakpoints(cluster)
    p = data["truth_position"]
    assert ins._pair_breakpoints(left, right) == [(p, p + 1), (None, p)]
    # The dominant opposed flanks have zero overlap, NOT a one-base TSD.
    assert p - (p + 1) + 1 == 0


def test_boundary_experiment_exposes_sensitivity_cost(locus, exclusive_pairing):
    data, cluster = locus
    candidates = ins._call_insertions(cluster, None, data["read_repeat"])
    assert len(candidates) == 1
    assert candidates[0].left_junction_reads == 0
    assert candidates[0].right_junction_reads == 1
    assert candidates[0].start == data["truth_position"]
    qualified = ins._call_validated_by_high_quality(
        *ins._candidate_junctions(candidates[0], cluster)
    )
    assert qualified == (data["sample"] == "cov30x_rep2")
    if qualified:
        converted = ins._as_supporting_junction(candidates[0], cluster)
        assert converted is not None
        assert converted.tsd == "supporting_junction"
        assert (converted.start, converted.end) == (9006303, 9006305)


def test_experiment_uses_legacy_boundary_distances(locus, exclusive_pairing):
    data, cluster = locus
    p = data["truth_position"]
    assert ins._pair_breakpoints(*breakpoints(cluster)) == [(p, p + 1), (None, p)]


@pytest.mark.parametrize(
    "left,right,expected",
    [
        ({100: [1, 2]}, {100: [3, 4]}, [(100, 100)]),  # genuine one-base TSD
        (
            {103: [1, 2], 153: [3, 4]},
            {100: [5, 6], 150: [7, 8]},
            [(103, 100), (153, 150)],
        ),  # distinct nearby insertions
        ({}, {100: [1, 2, 3]}, [(None, 100)]),  # valid one-sided candidate
        ({100: [1]}, {1: [2]}, [(100, 1)]),  # inclusive distance 99
        ({100: [1]}, {0: [2]}, [(100, 0)]),  # distance exactly 100
        ({100: [1]}, {201: [2]}, [(100, None), (None, 201)]),  # distance 101
        ({100: [1]}, {202: [2]}, [(100, None), (None, 202)]),
    ],
)
def test_pairing_controls(left, right, expected):
    assert ins._pair_breakpoints(left, right) == expected
