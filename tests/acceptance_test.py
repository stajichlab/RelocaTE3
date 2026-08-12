"""Acceptance gate: full 14-family run vs the simulation truth.

RelocaTE2's published benchmark on this exact rice Chr3 2 Mb dataset recovered
196 of the 200 simulated insertions (`bedtools window -w 10` against
``MSU7.Chr3_2M.ALL.gff``). This test runs RelocaTE3 with the same full
``RiceTE.fa`` library and asserts comparable recovery, so regressions in
sensitivity are caught.

This gate pins ``te_aligner="minimap2"`` explicitly even though the shipped
defaults are now blat/bwaaln (matching RelocaTE2). blat cannot live in the
default pixi environment, so a gate on the defaults would simply skip wherever
blat is absent -- including `pixi run test` and CI -- and the thresholds below
are calibrated against minimap2 anyway. `test_run_all_cli_end_to_end` covers
the shipped defaults instead, skipping when blat is unavailable.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from RelocaTE3.pipeline import run_sample
from RelocaTE3.ReadLibrary import ReadLibrary

DATA = Path(__file__).parent / "data"
R1 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"
R2 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"
RICETE = DATA / "TE_lib" / "RiceTE.fa"
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
        reads, str(RICETE), str(GENOME), tmp_path, threads=4, mismatch_allowance=2,
        te_aligner="minimap2", genome_aligner="minimap2",
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


RM_OUT = DATA / "sim_genome" / "MSU7.Chr3_2M.fa.RepeatMasker.out"


@pytest.mark.skipif(not RICETE.exists(), reason="RiceTE.fa not vendored")
@pytest.mark.skipif(
    shutil.which("blat") is None,
    reason="blat not on PATH (it is the default TE aligner; see pixi env 'blat')",
)
def test_run_all_cli_end_to_end(tmp_path: Path):
    """`relocaTE3 run-all` drives the whole pipeline from one command.

    Exercises the real stages (not stubs) including the optional
    find-reference and characterize branches, and asserts every artifact a
    RelocaTE2 user expects is produced.
    """
    from RelocaTE3.cli import main

    genome = tmp_path / "genome.fa"
    genome.write_bytes(GENOME.read_bytes())  # copy: run-all indexes in place

    rc = main(
        [
            "run-all",
            "-1", str(R1),
            "-2", str(R2),
            "-T", str(RICETE),
            "-g", str(genome),
            "-n", "HEG4",
            "-o", str(tmp_path / "out"),
            "--threads", "4",
            "--mismatch", "2",
            "--tsd", "UNK",
            "--te-name", "repeat",
            "--require-both-junctions",
            "--repeatmasker", str(RM_OUT),
            "--genotype",
        ]
    )
    assert rc == 0

    results = tmp_path / "out" / "results"
    nonref = results / "ALL.repeat.all_nonref_insert.txt"
    assert nonref.is_file() and nonref.stat().st_size > 0

    # steps 0/6 -- the RelocaTE2 output that was previously CLI-unreachable
    assert (tmp_path / "out" / "existingTE.bed").is_file()
    assert (results / "HEG4.all_ref_insert.gff").is_file()
    assert (results / "HEG4.all_ref_insert.txt").is_file()

    # step 7 -- genotyped calls
    char = results / "ALL.repeat.all_nonref_insert.characTErized.txt"
    assert char.is_file() and char.stat().st_size > 0

    # the calls must be real: recover a decent share of the 200 simulated sites
    calls = _load_intervals(results / "ALL.repeat.all_nonref_insert.characTErized.gff")
    truth = _load_intervals(TRUTH)
    recovered = sum(1 for t in truth if any(_overlaps(t, c, WINDOW) for c in calls))
    assert recovered > 100, f"only recovered {recovered}/200 insertions"
