"""Acceptance gate: full 14-family run vs the simulation truth.

RelocaTE2's published benchmark on this exact rice Chr3 2 Mb dataset recovered
196 of the 200 simulated insertions (`bedtools window -w 10` against
``MSU7.Chr3_2M.ALL.gff``). This test runs RelocaTE3 with the same full
``RiceTE.fa`` library and asserts comparable recovery, so regressions in
sensitivity are caught. RelocaTE3 uses minimap2 (not blat) and single-end
flank mapping, so it is expected to trail RelocaTE2 somewhat; the thresholds
below leave margin under the observed ~178/200 (~89%) recall at ~90% precision.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from RelocaTE3.pipeline import run_sample
from RelocaTE3.ReadLibrary import ReadLibrary

DATA = Path(__file__).parent / "data"
R1 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"
R2 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"
RICETE = DATA / "RiceTE.fa"
GENOME = DATA / "sim_genome" / "MSU7.Chr3_2M.fa"
TRUTH = DATA / "sim_genome" / "MSU7.Chr3_2M.ALL.gff"

WINDOW = 10  # bp, matching RelocaTE2's `bedtools window -w 10`


def _load_intervals(path: Path) -> list[tuple[str, int, int]]:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 5:
                continue
            rows.append((f[0], int(f[3]), int(f[4])))
    return rows


def _overlaps(a: tuple[str, int, int], b: tuple[str, int, int], window: int) -> bool:
    return a[0] == b[0] and a[1] - window <= b[2] and a[2] + window >= b[1]


@pytest.mark.skipif(not RICETE.exists(), reason="RiceTE.fa not vendored")
def test_acceptance_full_library(tmp_path: Path):
    reads = ReadLibrary([str(R1), str(R2)], "HEG4")
    gff = run_sample(
        reads, str(RICETE), str(GENOME), tmp_path, threads=4, mismatch_allowance=2
    )
    calls = _load_intervals(gff)
    truth = _load_intervals(TRUTH)
    assert len(truth) == 200

    recovered = sum(1 for t in truth if any(_overlaps(t, c, WINDOW) for c in calls))
    true_calls = sum(1 for c in calls if any(_overlaps(c, t, WINDOW) for t in truth))
    precision = true_calls / len(calls) if calls else 0.0

    # RelocaTE2 reference: 196/200. RelocaTE3 (minimap2) target: >= 170/200 recall.
    assert recovered >= 170, f"recall regressed: {recovered}/200"
    assert precision >= 0.85, f"precision regressed: {precision:.2f}"
