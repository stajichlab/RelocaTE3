#!/usr/bin/env python3
"""Rank remaining false positives from completed stabilized replay tables."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from audit_benchmark_errors import load, norm, pair_fps, truth_distance, write_table
from summarize_residual_errors import read_raw


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    manifest = json.loads((args.replay / "manifest.json").read_text())
    hashes, annotations, counts, loci = {}, [], Counter(), Counter()
    for index, task in enumerate(manifest["tasks"]):
        folder = args.replay / "tasks" / f"{index:03d}"
        if not (folder / ".complete").exists():
            raise RuntimeError(f"Incomplete task: {index}")
        r2 = load(folder / "relocate2/false_positive_calls.tsv", hashes)
        r3 = load(folder / "tied_evidence/false_positive_calls.tsv", hashes)
        truths = load(folder / "tied_evidence/matches.tsv", hashes)
        paths = list(
            (folder / "tied_evidence/results").glob("*.all_nonref_insert.raw.txt")
        )
        if len(paths) != 1:
            raise ValueError(paths)
        raw = read_raw(paths[0], hashes)
        a, b = pair_fps(r2, r3)
        ds = task["dataset"]
        counts[ds, "shared"] += len(b)
        counts[ds, "R2_only"] += len(r2) - len(a)
        counts[ds, "R3_only"] += len(r3) - len(b)
        for j, row in enumerate(r3):
            choices = [
                r
                for r in raw
                if r["chrom"] == row["chrom"]
                and r["position"] == row["position"]
                and r["te_family"] == row["te_family"]
                and r["tsd"] == row["tsd"]
            ]
            if len(choices) != 1:
                raise ValueError((task["sample"], row, choices))
            near = [
                t
                for t in truths
                if t["chrom"] == row["chrom"] and truth_distance(t, row) <= 10
            ]
            category = (
                "same_family_near_truth"
                if any(norm(t["te_family"]) == norm(row["te_family"]) for t in near)
                else "different_family_near_truth"
                if near
                else "outside_truth_window"
            )
            state = "shared" if j in b else "R3_only"
            if state == "R3_only":
                loci[ds, row["chrom"], row["position"], row["te_family"]] += 1
            annotations.append(
                dict(
                    task=index,
                    dataset=ds,
                    state=state,
                    category=category,
                    **row,
                    raw_evidence=json.dumps(choices[0], sort_keys=True),
                )
            )
    args.output.mkdir(parents=True)
    write_table(args.output / "false_positives.tsv", annotations)
    write_table(
        args.output / "counts.tsv",
        [dict(dataset=d, state=s, count=n) for (d, s), n in sorted(counts.items())],
    )
    write_table(
        args.output / "repeated_R3_only.tsv",
        [
            dict(dataset=d, chrom=c, position=p, family=f, observations=n)
            for (d, c, p, f), n in sorted(
                loci.items(), key=lambda item: (-item[1], item[0])
            )
        ],
    )
    with (args.output / "provenance.json").open("x") as handle:
        json.dump(
            dict(
                input_sha256=hashes,
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            ),
            handle,
            indent=2,
        )
    print(dict(counts))


if __name__ == "__main__":
    main()
