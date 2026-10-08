#!/usr/bin/env python3
"""Summarize completed family replay tables; no alignment scans or caller runs."""

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, records):
    if not records:
        return
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", required=True, type=Path)
    parser.add_argument("--historical", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    manifest = json.loads((args.replay / "manifest.json").read_text())
    for name, digest in manifest["snapshot_sha256"].items():
        if hashlib.sha256((args.replay / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Frozen source changed: {name}")
    grouped = defaultdict(Counter)
    samples, changes, paired, provenance = [], [], [], {}
    historical_agreement = geometry_agreement = 0
    variants = ("relocate2", "historical_r3", "original_evidence", "tied_evidence")

    def read(path):
        provenance[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        return rows(path)

    for index, task in enumerate(manifest["tasks"]):
        folder = args.replay / "tasks" / f"{index:03d}"
        if not (folder / ".complete").is_file():
            raise RuntimeError(f"Incomplete task: {index}")
        comparison = json.loads((folder / "comparison.json").read_text())
        dataset, sample = task["dataset"], task["sample"]
        if (comparison["dataset"], comparison["sample"]) != (dataset, sample):
            raise ValueError(f"Task identity mismatch: {index}")
        historical_agreement += comparison["original_evidence_matches_historical_calls"]
        geometry_agreement += comparison["geometry_and_genotypes_unchanged"]
        coverage = re.search(r"cov(\d+)x", sample).group(1)
        events = {}
        for variant in variants:
            source = (
                args.historical
                / "reports/datasets"
                / dataset
                / "per_sample"
                / "relocate3-blat-bwaaln"
                / sample
                if variant == "historical_r3"
                else folder / variant
            )
            matches = read(source / "matches.tsv")
            precision = read(source / "precision.tsv")[0]
            events[variant] = {row["event_id"]: row for row in matches}
            if len(events[variant]) != len(matches):
                raise ValueError("Duplicate truth event ID")
            tp = sum(int(row["matched"]) for row in matches)
            fp = int(precision["false_positive_calls"])
            if tp != int(precision["matched_calls"]):
                raise ValueError("Match and precision tables disagree")
            samples.append(
                dict(
                    task=index,
                    dataset=dataset,
                    sample=sample,
                    variant=variant,
                    truth=len(matches),
                    TP=tp,
                    FP=fp,
                    FN=len(matches) - tp,
                )
            )
            divergence = matches[0].get("divergence_percent", "") or "NA"
            for row in matches:
                for axis, level in (
                    ("all", "all"),
                    ("coverage", coverage),
                    ("divergence", divergence),
                    ("biological_class", row["biological_class"]),
                    ("te_group", row["te_group"]),
                    (
                        "class_fraction",
                        row["biological_class"] + ":" + row["cellular_fraction"],
                    ),
                ):
                    counts = grouped[dataset, variant, axis, level]
                    counts["truth"] += 1
                    counts["TP"] += int(row["matched"])
                    counts["status_correct"] += int(row.get("status_correct") or 0)
                    counts["tsd_exact"] += int(row.get("tsd_exact") or 0)
            # FP belongs to calls, not truth strata; never repeat FPs per truth class.
            for axis, level in (
                ("all", "all"),
                ("coverage", coverage),
                ("divergence", divergence),
            ):
                grouped[dataset, variant, axis, level]["FP"] += fp
        if any(set(events[v]) != set(events["relocate2"]) for v in variants):
            raise ValueError("Truth event sets disagree")
        for before, after in (
            ("historical_r3", "original_evidence"),
            ("original_evidence", "tied_evidence"),
            ("relocate2", "tied_evidence"),
        ):
            counts = Counter()
            for event, previous in events[before].items():
                current = events[after][event]
                a, b = int(previous["matched"]), int(current["matched"])
                counts[(a, b)] += 1
                if a != b:
                    changes.append(
                        dict(
                            task=index,
                            dataset=dataset,
                            sample=sample,
                            before=before,
                            after=after,
                            event_id=event,
                            change="gain" if b else "loss",
                            te_family=current["te_family"],
                            te_group=current["te_group"],
                            position=current["position"],
                            biological_class=current["biological_class"],
                            cellular_fraction=current["cellular_fraction"],
                        )
                    )
            paired.append(
                dict(
                    dataset=dataset,
                    sample=sample,
                    before=before,
                    after=after,
                    shared=counts[1, 1],
                    gains=counts[0, 1],
                    losses=counts[1, 0],
                )
            )
    output = []
    for (dataset, variant, axis, level), c in sorted(grouped.items()):
        fp = c.get("FP")
        output.append(
            dict(
                dataset=dataset,
                variant=variant,
                axis=axis,
                level=level,
                truth=c["truth"],
                TP=c["TP"],
                FN=c["truth"] - c["TP"],
                FP=fp if fp is not None else "NA",
                recall=c["TP"] / c["truth"],
                precision=c["TP"] / (c["TP"] + fp)
                if fp is not None and c["TP"] + fp
                else "NA",
                status_correct=c["status_correct"],
                tsd_exact=c["tsd_exact"],
            )
        )
    args.output.mkdir(parents=True)
    write(args.output / "summary.tsv", output)
    write(args.output / "samples.tsv", samples)
    write(args.output / "changed_events.tsv", changes)
    write(args.output / "paired.tsv", paired)
    with (args.output / "provenance.json").open("x") as handle:
        json.dump(
            dict(
                recorded=datetime.now().astimezone().isoformat(),
                tasks=len(manifest["tasks"]),
                historical_agreement=historical_agreement,
                geometry_and_genotype_agreement=geometry_agreement,
                coordinate_policy=manifest["coordinate_policy"],
                window_bp=10,
                manifest_sha256=hashlib.sha256(
                    (args.replay / "manifest.json").read_bytes()
                ).hexdigest(),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                input_sha256=provenance,
            ),
            handle,
            indent=2,
        )
    print(
        f"Verified {len(manifest['tasks'])} completed tasks; historical agreement "
        f"{historical_agreement}; geometry/genotype agreement {geometry_agreement}"
    )
    for row in output:
        if row["axis"] == "all":
            print(row)


if __name__ == "__main__":
    main()
