"""Score long-read clusters against a truth table of insertion sites (L2 gate).

A truth site ``(chrom, position)`` is **found** when any cluster with
``sequence_id == chrom`` has ``abs(cluster.locus.start - position) <=
max_dist``. A sample cluster is **unmatched** when no truth site lies within
``max_dist`` of its start. Every cluster from the reference-only control is a
false positive; **chance recall** is the recall the control clusters alone
reach against the same truth sites. Each metric is also given for clusters
with ``n_reads >= 2`` (``*_min2``).

Inputs:
    --truth     truth_events.tsv (columns include ``chrom`` and ``position``,
                0-based insertion point)
    --clusters  <sample>.clusters.tsv.gz from ``relocaTE3 lr-cluster``
    --control   optional clusters file from the reference-only control reads

Output: one TSV row on stdout (columns in ``COLUMNS``), preceded by the header
line with ``--header``. Control columns are empty without ``--control``.

Usage:
    score_sites.py --truth T --clusters C [--control K] [--max-dist 50] \
        --label hifi_cov20x [--header]
"""

from __future__ import annotations

import argparse
import bisect
import csv
import sys
from collections import defaultdict
from collections.abc import Iterable, Sequence
from pathlib import Path

from RelocaTE3.longread.tables import read_clusters

COLUMNS = (
    "label", "n_truth", "found", "recall", "found_min2", "recall_min2",
    "n_clusters", "unmatched", "unmatched_min2", "control_clusters",
    "control_clusters_min2", "chance_recall", "chance_recall_min2",
)  # fmt: skip

Sites = list[tuple[str, int]]
# (chrom, start, n_reads) for each cluster.
ClusterSites = list[tuple[str, int, int]]


def read_truth(path: str | Path) -> Sites:
    """Return ``(chrom, position)`` for every row of a truth table."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = {"chrom", "position"} - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: truth table lacks columns {sorted(missing)}")
        return [(row["chrom"], int(row["position"])) for row in reader]


def read_cluster_sites(path: str | Path) -> ClusterSites:
    """Return ``(sequence_id, start, n_reads)`` for every cluster in ``path``."""
    _meta, clusters = read_clusters(path)
    return [(c.locus.sequence_id, c.locus.start, c.n_reads) for c in clusters]


def _by_chrom(positions: Iterable[tuple[str, int]]) -> dict[str, list[int]]:
    out: dict[str, list[int]] = defaultdict(list)
    for chrom, pos in positions:
        out[chrom].append(pos)
    for values in out.values():
        values.sort()
    return out


def _near(sorted_positions: Sequence[int], pos: int, max_dist: int) -> bool:
    """True when some value in ``sorted_positions`` is within ``max_dist``."""
    i = bisect.bisect_left(sorted_positions, pos - max_dist)
    return i < len(sorted_positions) and sorted_positions[i] <= pos + max_dist


def _count_near(queries: Sites, targets: Sites, max_dist: int) -> int:
    """Number of ``queries`` with a ``target`` on the same chrom within reach."""
    index = _by_chrom(targets)
    empty: list[int] = []
    return sum(_near(index.get(chrom, empty), pos, max_dist) for chrom, pos in queries)


def _sites(clusters: ClusterSites, min_reads: int) -> Sites:
    return [(chrom, start) for chrom, start, n in clusters if n >= min_reads]


def _ratio(num: int, den: int) -> str:
    return f"{num / den:.4f}" if den else ""


def score(
    truth_path: str | Path,
    clusters_path: str | Path,
    control_path: str | Path | None,
    max_dist: int,
    label: str = "",
) -> dict[str, object]:
    """Score one sample's clusters (and optional control) against the truth.

    Returns:
        A dict keyed by ``COLUMNS``; recall values are 4-decimal strings,
        control fields are ``""`` when ``control_path`` is None.
    """
    if max_dist < 0:
        raise ValueError(f"max_dist must be >= 0, got {max_dist}")
    truth = read_truth(truth_path)
    sample = read_cluster_sites(clusters_path)
    n_truth = len(truth)
    found = _count_near(truth, _sites(sample, 1), max_dist)
    found_min2 = _count_near(truth, _sites(sample, 2), max_dist)
    row: dict[str, object] = {
        "label": label,
        "n_truth": n_truth,
        "found": found,
        "recall": _ratio(found, n_truth),
        "found_min2": found_min2,
        "recall_min2": _ratio(found_min2, n_truth),
        "n_clusters": len(sample),
        "unmatched": len(sample) - _count_near(_sites(sample, 1), truth, max_dist),
        "unmatched_min2": len(_sites(sample, 2))
        - _count_near(_sites(sample, 2), truth, max_dist),
        "control_clusters": "",
        "control_clusters_min2": "",
        "chance_recall": "",
        "chance_recall_min2": "",
    }
    if control_path is not None:
        control = read_cluster_sites(control_path)
        row["control_clusters"] = len(control)
        row["control_clusters_min2"] = len(_sites(control, 2))
        row["chance_recall"] = _ratio(
            _count_near(truth, _sites(control, 1), max_dist), n_truth
        )
        row["chance_recall_min2"] = _ratio(
            _count_near(truth, _sites(control, 2), max_dist), n_truth
        )
    return row


def main(argv: Sequence[str] | None = None) -> int:
    """Print one TSV row (and optionally the header) to stdout."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--truth", required=True, help="truth_events.tsv")
    parser.add_argument("--clusters", required=True, help="sample clusters.tsv.gz")
    parser.add_argument("--control", default=None,
                        help="reference-only control clusters.tsv.gz")  # fmt: skip
    parser.add_argument("--max-dist", type=int, default=50,
                        help="largest |cluster start - truth position| (bp)")  # fmt: skip
    parser.add_argument("--label", required=True, help="row label")
    parser.add_argument("--header", action="store_true", help="print the header")
    args = parser.parse_args(argv)
    row = score(args.truth, args.clusters, args.control, args.max_dist, args.label)
    if args.header:
        sys.stdout.write("\t".join(COLUMNS) + "\n")
    sys.stdout.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
