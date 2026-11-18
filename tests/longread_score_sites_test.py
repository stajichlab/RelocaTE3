"""L2 gate scorer: site recall of long-read clusters against a truth table.

``validation/longread_l2/score_sites.py`` counts a truth site as found when a
cluster on the same sequence starts within ``max_dist`` bp of it, and scores
the reference-only control's clusters as pure false positives.
"""

from __future__ import annotations

import sys
from pathlib import Path

from RelocaTE3.longread.model import Cluster, Locus, ReadType
from RelocaTE3.longread.tables import write_clusters

REPO = Path(__file__).resolve().parents[1]
SCORER_DIR = REPO / "validation" / "longread_l2"

sys.path.insert(0, str(SCORER_DIR))
import score_sites  # noqa: E402

COLUMNS = [
    "label", "n_truth", "found", "recall", "found_min2", "recall_min2",
    "n_clusters", "unmatched", "unmatched_min2", "control_clusters",
    "control_clusters_min2", "chance_recall", "chance_recall_min2",
]  # fmt: skip

# Three truth sites on Chr3 (0-based insertion points).
TRUTH_SITES = [("Chr3", 1000), ("Chr3", 5000), ("Chr3", 9000)]


def _truth(tmp_path: Path) -> Path:
    path = tmp_path / "truth_events.tsv"
    lines = ["event_id\tchrom\tposition\tte_family\tstrand\ttsd\ttsd_length\t"
             "biological_class\n"]  # fmt: skip
    for i, (chrom, pos) in enumerate(TRUTH_SITES):
        lines.append(f"e{i}\t{chrom}\t{pos}\tmPing\t+\tTAA\t3\thomozygous\n")
    path.write_text("".join(lines))
    return path


def _clusters(tmp_path: Path, name: str, sites: list[tuple[str, int, int]]) -> Path:
    """Write clusters at ``(chrom, start, n_reads)`` sites."""
    path = tmp_path / f"{name}.clusters.tsv.gz"
    clusters = [
        Cluster(
            cluster_id=f"c{i}",
            locus=Locus("ref", chrom, start, start, "+"),
            n_reads=n_reads,
            n_ins=n_reads,
            n_split_ins=0,
            n_clip_left=0,
            n_clip_right=0,
            insert_lengths=(500,) * n_reads,
            signature_ids=tuple(f"s{i}_{k}" for k in range(n_reads)),
        )  # fmt: skip
        for i, (chrom, start, n_reads) in enumerate(sites)
    ]
    write_clusters(path, clusters, read_type=ReadType.HIFI, reference_id="ref")
    return path


def test_boundary_50_found_51_not(tmp_path):
    truth = _truth(tmp_path)
    clusters = _clusters(tmp_path, "s", [
        ("Chr3", 1050, 1),   # +50: found
        ("Chr3", 4949, 3),   # -51: not found
        ("Chr3", 9000, 2),   # exact
        ("Chr3", 20000, 1),  # far from any site
    ])  # fmt: skip
    row = score_sites.score(truth, clusters, None, 50)
    assert row["n_truth"] == 3
    assert row["found"] == 2
    assert row["recall"] == "0.6667"
    assert row["n_clusters"] == 4
    # 4949 and 20000 are not within 50 bp of any truth site.
    assert row["unmatched"] == 2


def test_other_chromosome_never_matches(tmp_path):
    truth = _truth(tmp_path)
    clusters = _clusters(tmp_path, "s", [
        ("Chr1", 1000, 5), ("Chr1", 5000, 5), ("Chr1", 9000, 5),
    ])  # fmt: skip
    row = score_sites.score(truth, clusters, None, 50)
    assert row["found"] == 0
    assert row["recall"] == "0.0000"
    assert row["unmatched"] == 3


def test_min2_columns_use_only_clusters_with_two_reads(tmp_path):
    truth = _truth(tmp_path)
    clusters = _clusters(tmp_path, "s", [
        ("Chr3", 1000, 1),   # found, single read
        ("Chr3", 5010, 2),   # found, two reads
        ("Chr3", 9000, 4),   # found, four reads
        ("Chr3", 30000, 1),  # unmatched, single read
        ("Chr3", 40000, 2),  # unmatched, two reads
    ])  # fmt: skip
    row = score_sites.score(truth, clusters, None, 50)
    assert row["found"] == 3
    assert row["recall"] == "1.0000"
    assert row["found_min2"] == 2
    assert row["recall_min2"] == "0.6667"
    assert row["n_clusters"] == 5
    assert row["unmatched"] == 2
    assert row["unmatched_min2"] == 1


def test_control_counts_and_chance_recall(tmp_path):
    truth = _truth(tmp_path)
    clusters = _clusters(tmp_path, "s", [("Chr3", 1000, 3)])
    control = _clusters(tmp_path, "k", [
        ("Chr3", 1020, 1),   # near site 1000 by chance, single read
        ("Chr3", 9040, 2),   # near site 9000 by chance, two reads
        ("Chr3", 60000, 2),  # nowhere near a site
    ])  # fmt: skip
    row = score_sites.score(truth, clusters, control, 50)
    assert row["control_clusters"] == 3
    assert row["control_clusters_min2"] == 2
    assert row["chance_recall"] == "0.6667"
    assert row["chance_recall_min2"] == "0.3333"


def test_no_control_leaves_control_fields_empty(tmp_path):
    truth = _truth(tmp_path)
    clusters = _clusters(tmp_path, "s", [("Chr3", 1000, 3)])
    row = score_sites.score(truth, clusters, None, 50)
    for key in ("control_clusters", "control_clusters_min2", "chance_recall",
                "chance_recall_min2"):  # fmt: skip
        assert row[key] == ""


def test_main_prints_header_and_row(tmp_path, capsys):
    truth = _truth(tmp_path)
    clusters = _clusters(tmp_path, "s", [("Chr3", 1000, 3), ("Chr3", 5051, 1)])
    control = _clusters(tmp_path, "k", [("Chr3", 9000, 1)])
    rc = score_sites.main([
        "--truth", str(truth), "--clusters", str(clusters),
        "--control", str(control), "--label", "hifi_cov20x", "--header",
    ])  # fmt: skip
    assert rc == 0
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 2
    assert lines[0].split("\t") == COLUMNS
    assert lines[1].split("\t") == [
        "hifi_cov20x", "3", "1", "0.3333", "1", "0.3333", "2", "1", "0", "1",
        "0", "0.3333", "0.0000",
    ]  # fmt: skip


def test_main_without_header_or_control(tmp_path, capsys):
    truth = _truth(tmp_path)
    clusters = _clusters(tmp_path, "s", [("Chr3", 5000, 2)])
    rc = score_sites.main([
        "--truth", str(truth), "--clusters", str(clusters), "--max-dist", "10",
        "--label", "x",
    ])  # fmt: skip
    assert rc == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines == ["x\t3\t1\t0.3333\t1\t0.3333\t1\t0\t0\t\t\t\t"]
