#!/usr/bin/env python3
"""Small-table, paired accuracy audit of the frozen intact-read experiment."""

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime
import hashlib
import json
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--replay", required=True, type=Path)
    p.add_argument("--stabilization", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    manifest = json.loads((args.replay / "manifest.json").read_text())
    hashes = {}
    for name, digest in manifest["snapshot_sha256"].items():
        if hashlib.sha256((args.replay / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Frozen file changed: {name}")

    def read(path):
        data = path.read_bytes()
        hashes[str(path)] = hashlib.sha256(data).hexdigest()
        return list(csv.DictReader(data.decode().splitlines(), delimiter="\t"))

    totals, strata = defaultdict(Counter), defaultdict(Counter)
    changed, samples, paired = [], [], []
    for index, task in enumerate(manifest["tasks"]):
        folder = args.replay / "tasks" / f"{index:03d}"
        if not (folder / ".complete").exists():
            raise ValueError(f"Incomplete task: {index}")
        comparison = json.loads((folder / "comparison.json").read_text())
        assert comparison["baseline_reproduced"]
        assert (comparison["dataset"], comparison["sample"]) == (
            task["dataset"],
            task["sample"],
        )
        ds, sample = task["dataset"], task["sample"]
        events = {}
        for variant in ("relocate2", "baseline", "intact_candidate"):
            source = (
                args.stabilization / "tasks" / f"{index:03d}" / "relocate2"
                if variant == "relocate2"
                else folder / variant
            )
            rows = read(source / "matches.tsv")
            precision = read(source / "precision.tsv")[0]
            events[variant] = {r["event_id"]: r for r in rows}
            assert len(events[variant]) == len(rows)
            tp = sum(int(r["matched"]) for r in rows)
            fp = int(precision["false_positive_calls"])
            assert tp == int(precision["matched_calls"])
            values = dict(truth=len(rows), TP=tp, FP=fp, FN=len(rows) - tp)
            totals[ds, variant].update(values)
            samples.append(dict(dataset=ds, sample=sample, variant=variant, **values))
            for r in rows:
                for axis, level in (
                    ("biological_class", r["biological_class"]),
                    (
                        "class_fraction",
                        r["biological_class"] + ":" + r["cellular_fraction"],
                    ),
                    ("te_group", r["te_group"]),
                    ("divergence", r.get("divergence_percent") or "NA"),
                ):
                    strata[ds, variant, axis, level].update(
                        dict(
                            truth=1,
                            TP=int(r["matched"]),
                            status_correct=int(r.get("status_correct") or 0),
                        )
                    )
        assert (
            set(events["baseline"])
            == set(events["intact_candidate"])
            == set(events["relocate2"])
        )
        for before in ("baseline", "relocate2"):
            c = Counter()
            for e, old in events[before].items():
                new = events["intact_candidate"][e]
                a, b = int(old["matched"]), int(new["matched"])
                c[(a, b)] += 1
                if a != b:
                    changed.append(
                        dict(
                            dataset=ds,
                            sample=sample,
                            task=index,
                            before=before,
                            event_id=e,
                            change="gain" if b else "loss",
                            R2_detected=events["relocate2"][e]["matched"],
                            te_family=old["te_family"],
                            biological_class=old["biological_class"],
                            cellular_fraction=old["cellular_fraction"],
                            chrom=old["chrom"],
                            position=old["position"],
                            baseline_call_position=events["baseline"][e][
                                "call_position"
                            ],
                        )
                    )
            paired.append(
                dict(
                    dataset=ds,
                    sample=sample,
                    before=before,
                    shared=c[1, 1],
                    gains=c[0, 1],
                    losses=c[1, 0],
                )
            )
    args.output.mkdir(parents=True)

    def write(name, rows):
        with (args.output / name).open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    summary = [
        dict(
            dataset=d,
            variant=v,
            **dict(c),
            precision=c["TP"] / (c["TP"] + c["FP"]),
            recall=c["TP"] / c["truth"],
        )
        for (d, v), c in sorted(totals.items())
    ]
    write("summary.tsv", summary)
    write("samples.tsv", samples)
    write("paired.tsv", paired)
    write("changed_events.tsv", changed)
    write(
        "strata.tsv",
        [
            dict(dataset=d, variant=v, axis=a, level=level, **dict(c))
            for (d, v, a, level), c in sorted(strata.items())
        ],
    )
    with (args.output / "provenance.json").open("x") as handle:
        json.dump(
            dict(
                recorded=datetime.now().astimezone().isoformat(),
                completed_tasks=len(manifest["tasks"]),
                input_sha256=hashes,
                manifest_sha256=hashlib.sha256(
                    (args.replay / "manifest.json").read_bytes()
                ).hexdigest(),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            ),
            handle,
            indent=2,
        )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
