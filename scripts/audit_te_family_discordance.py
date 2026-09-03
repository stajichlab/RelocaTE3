#!/usr/bin/env python3
"""Audit RelocaTE2 compound TE labels against RelocaTE3 read evidence.

The riceTElib benchmark stores enough intermediate data to separate TE-family
votes from breakpoint-spanning junction reads and bracketing supporting mates.
This script reconstructs that evidence for shared calls whose normalized family
labels differ between RelocaTE2 and RelocaTE3.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path

from RelocaTE3.insertions import (
    _as_supporting_junction,
    _call_insertions,
    _consolidate_same_start,
    _stream_clusters,
)

_JUNCTION_SUFFIX = re.compile(r":(start|end):[53]$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--benchmark-root",
        required=True,
        type=Path,
        help="relocate-benchmark project root",
    )
    parser.add_argument("--dataset", default="ricetelib")
    parser.add_argument("--sample", default="cov5x_rep1")
    parser.add_argument("--r2-caller", default="relocate2")
    parser.add_argument("--r3-caller", default="relocate3-blat-bwaaln")
    parser.add_argument("--chrom", default="Chr1")
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_normalized(path: Path) -> dict[tuple[str, int], dict[str, str]]:
    with path.open() as handle:
        return {
            (row["chrom"], int(row["position"])): row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def read_repeat_table(path: Path) -> dict[str, tuple[str, str]]:
    assignments: dict[str, tuple[str, str]] = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                assignments[fields[0]] = (
                    fields[1],
                    fields[2] if len(fields) >= 3 else "",
                )
    return assignments


def read_r2_evidence(
    reads_list: Path,
    chrom: str,
    targets: set[tuple[str, int]],
) -> dict[tuple[str, int], dict[str, list[str]]]:
    evidence: dict[tuple[str, int], dict[str, list[str]]] = {}
    with reads_list.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4:
                continue
            _te, coordinate, evidence_type, read_names = fields
            observed_chrom, interval = coordinate.split(":", 1)
            key = (observed_chrom, int(interval.split("..", 1)[0]))
            if observed_chrom == chrom and key in targets:
                evidence.setdefault(key, {})[evidence_type] = [
                    name for name in read_names.split(",") if name
                ]
    return evidence


def junction_family(
    read_name: str,
    assignments: dict[str, tuple[str, str]],
) -> str:
    key = _JUNCTION_SUFFIX.sub("", read_name)
    return assignments.get(key, ("NA", ""))[0]


def supporting_family(
    read_name: str,
    assignments: dict[str, tuple[str, str]],
) -> str:
    # This is RelocaTE2's lookup order in insertion_family_supporting().
    for key in (
        read_name,
        f"{read_name}/1",
        f"{read_name}/2",
        f"{read_name}.f",
        f"{read_name}.r",
    ):
        if key in assignments:
            return assignments[key][0]
    return "NA"


def count_string(counts: Counter[str]) -> str:
    return ",".join(
        f"{family}={count}"
        for family, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    )


def primary(counts: Counter[str]) -> str:
    informative = Counter(
        {family: count for family, count in counts.items() if family != "NA"}
    )
    if not informative:
        return "NA"
    return sorted(informative.items(), key=lambda item: (-item[1], item[0]))[0][0]


def reconstruct_r3_evidence(
    bam_path: Path,
    assignments: dict[str, tuple[str, str]],
    targets: set[tuple[str, int]],
    chrom: str,
) -> dict[tuple[str, int], tuple[Counter[str], Counter[str]]]:
    reconstructed: dict[tuple[str, int], tuple[Counter[str], Counter[str]]] = {}
    for cluster in _stream_clusters(str(bam_path), assignments):
        if cluster.chrom != chrom:
            continue
        raw_calls = _call_insertions(
            cluster,
            genome=None,
            read_repeat=assignments,
        )
        for original in _consolidate_same_start(raw_calls, cluster, assignments):
            reported = _as_supporting_junction(original, cluster) or original
            key = (reported.chrom, reported.start)
            if key not in targets:
                continue

            junction_counts = Counter(
                junction_family(name, assignments) for name in original.read_names
            )
            supporting_counts: Counter[str] = Counter()
            for name, gstart, gend, strand, _sequence in cluster.support:
                supports_call = (
                    strand == "+" and gend <= original.start
                ) or (
                    strand == "-" and gstart >= original.end
                )
                if supports_call:
                    supporting_counts[supporting_family(name, assignments)] += 1
            reconstructed[key] = (junction_counts, supporting_counts)
    return reconstructed


def main() -> None:
    args = parse_args()
    root = args.benchmark_root
    run_root = root / "runs" / args.dataset
    r2_root = run_root / args.r2_caller / args.sample
    r3_root = run_root / args.r3_caller / args.sample

    r2_calls = read_normalized(r2_root / "calls.normalized.tsv")
    r3_calls = read_normalized(r3_root / "calls.normalized.tsv")
    targets = {
        key
        for key in r2_calls.keys() & r3_calls.keys()
        if key[0] == args.chrom
        and r2_calls[key]["te_family"] != r3_calls[key]["te_family"]
    }

    r2_raw = r2_root / "raw" / "repeat"
    r2_assignments = read_repeat_table(
        r2_raw / "te_containing_fq" / f"{args.chrom}.read_repeat_name.split.txt"
    )
    r2_evidence = read_r2_evidence(
        r2_raw / "results" / f"{args.chrom}.repeat.reads.list",
        args.chrom,
        targets,
    )

    r3_raw = r3_root / "raw"
    r3_assignments = read_repeat_table(
        r3_raw / "te_containing" / f"{args.sample}.read_repeat_name.txt"
    )
    r3_evidence = reconstruct_r3_evidence(
        r3_raw / f"{args.sample}.repeat.bwaaln.sorted.bam",
        r3_assignments,
        targets,
        args.chrom,
    )

    missing_r2 = sorted(targets - r2_evidence.keys())
    missing_r3 = sorted(targets - r3_evidence.keys())
    if missing_r2 or missing_r3:
        raise RuntimeError(
            f"Failed to reconstruct all targets; R2 missing={missing_r2}, "
            f"R3 missing={missing_r3}"
        )
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)

    columns = [
        "chrom",
        "position",
        "r2_label",
        "r3_primary",
        "r2_junction_support",
        "r2_supporting_support",
        "r3_junction_support",
        "r3_supporting_support",
        "r2_compound_explained",
        "r3_retains_r2_supporting_primary",
        "classification",
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for key in sorted(targets):
            types = r2_evidence[key]
            r2_junction = Counter(
                junction_family(name, r2_assignments)
                for name in types.get("Junction_reads", [])
            )
            r2_supporting = Counter(
                supporting_family(name, r2_assignments)
                for evidence_type in (
                    "Left_supporting_reads",
                    "Right_supporting_reads",
                )
                for name in types.get(evidence_type, [])
            )
            r3_junction, r3_supporting = r3_evidence[key]
            r2_junction_primary = primary(r2_junction)
            r2_supporting_primary = primary(r2_supporting)
            expected_r2_label = f"{r2_junction_primary}/{r2_supporting_primary}"
            writer.writerow(
                {
                    "chrom": key[0],
                    "position": key[1],
                    "r2_label": r2_calls[key]["te_family"],
                    "r3_primary": r3_calls[key]["te_family"],
                    "r2_junction_support": count_string(r2_junction),
                    "r2_supporting_support": count_string(r2_supporting),
                    "r3_junction_support": count_string(r3_junction),
                    "r3_supporting_support": count_string(r3_supporting),
                    "r2_compound_explained": str(
                        r2_calls[key]["te_family"] == expected_r2_label
                    ).lower(),
                    "r3_retains_r2_supporting_primary": str(
                        primary(r3_supporting) == r2_supporting_primary
                    ).lower(),
                    "classification": "junction_vs_supporting_family_disagreement",
                }
            )

    print(f"Audited {len(targets)} discordant shared calls")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
