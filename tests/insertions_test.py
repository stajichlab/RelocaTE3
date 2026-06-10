"""Tests for insertion calling (Step 5)."""

from __future__ import annotations

import re
from pathlib import Path

import pysam

from RelocaTE3.insertions import (
    _call_insertions,
    _Cluster,
    _junction_info,
    _pair_breakpoints,
)
from RelocaTE3.models import JunctionObservation
from RelocaTE3.pipeline import run_sample
from RelocaTE3.ReadLibrary import ReadLibrary

DATA = Path(__file__).parent / "data"
R1 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"
R2 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"
TELIB = DATA / "mping.fa"
GENOME = DATA / "sim_genome" / "MSU7.Chr3_2M.fa"
TRUTH = DATA / "sim_genome" / "MSU7.Chr3_2M.ALL.gff"


def test_junction_info_breakpoints():
    # + strand: start-flank -> right edge at gstart; end-flank -> left edge at gend
    assert _junction_info("r/1:start:5", "+", 100, 200) == ("right", 100, "5")
    assert _junction_info("r/1:end:5", "+", 100, 200) == ("left", 200, "5")
    # - strand inverts which read end is TE-adjacent
    assert _junction_info("r/1:start:3", "-", 100, 200) == ("left", 200, "3")
    assert _junction_info("r/1:end:3", "-", 100, 200) == ("right", 100, "3")
    # non-junction (middle / untagged)
    assert _junction_info("r/1:middle", "+", 100, 200) is None
    assert _junction_info("r/1", "+", 100, 200) is None


def test_te_orientation():
    def obs(side, te_end):
        return JunctionObservation("r", side, 1, "+", "mPing", te_end)

    assert obs("left", "5").te_orientation == "+"
    assert obs("right", "3").te_orientation == "+"
    assert obs("right", "5").te_orientation == "-"
    assert obs("left", "3").te_orientation == "-"


def test_pair_breakpoints():
    # one left + one right within the TSD window -> a single paired insertion
    assert _pair_breakpoints([100], [98]) == [(100, 98)]
    # far apart -> two one-sided sub-insertions
    assert set(_pair_breakpoints([100], [500])) == {(100, None), (None, 500)}
    # two distinct insertions in one cluster -> two pairs
    pairs = set(_pair_breakpoints([100, 5000], [98, 4998]))
    assert pairs == {(100, 98), (5000, 4998)}


def test_call_insertions_splits_two_sites():
    """A cluster with two well-separated junction pairs yields two insertions."""
    cluster = _Cluster("Chr3")
    cluster.extend(1000, 5002)
    # site 1 near 1000, site 2 near 5000
    cluster.junctions = [
        JunctionObservation("a:end:5", "left", 1002, "+", "mPing", "5"),
        JunctionObservation("b:start:5", "right", 1000, "+", "mPing", "5"),
        JunctionObservation("c:end:5", "left", 5002, "+", "mPing", "5"),
        JunctionObservation("d:start:5", "right", 5000, "+", "mPing", "5"),
    ]
    with pysam.FastaFile(str(GENOME)) as genome:
        calls = _call_insertions(cluster, genome)
    starts = sorted(i.start for i in calls)
    assert len(calls) == 2
    assert starts == [1000, 5000]
    for ins in calls:
        assert ins.left_junction_reads == 1 and ins.right_junction_reads == 1
        assert len(ins.tsd) == 3  # 3 bp TSD from the 2 bp overlap span


def test_find_insertions_recovers_mping(tmp_path: Path):
    """End-to-end: most calls match true mPing insertions within 10 bp."""
    reads = ReadLibrary([str(R1), str(R2)], "HEG4")
    gff = run_sample(reads, str(TELIB), str(GENOME), tmp_path, threads=4)
    assert gff.exists()

    attr_re = re.compile(r"TSD=.*Right_junction_reads=\d+.*Left_support_reads=\d+")
    calls = []
    with open(gff) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            assert attr_re.search(f[8])  # required attributes present
            assert f[6] in ("+", "-")  # TE orientation assigned
            calls.append((f[0], int(f[3]), int(f[4])))

    assert len(calls) >= 10
    truth = _load_truth_mping(TRUTH)
    matched = sum(1 for c in calls if _near(c, truth))
    assert matched >= 0.7 * len(calls)
    assert matched >= 12


def _load_truth_mping(path: Path) -> list[tuple[str, int, int]]:
    coords = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or "ping" not in line.lower():
                continue
            f = line.split("\t")
            coords.append((f[0], int(f[3]), int(f[4])))
    return coords


def _near(
    call: tuple[str, int, int], truths: list[tuple[str, int, int]], window: int = 10
) -> bool:
    chrom, start, end = call
    return any(
        c == chrom and start - window <= e and end + window >= s for c, s, e in truths
    )
