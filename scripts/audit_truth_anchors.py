#!/usr/bin/env python3
"""Identify exact long-TSD anchor matches currently counted as FP/FN pairs.

Diagnostic only: never change official reports or claim a complete rescore.
For simulator position p and TSD length L, the leftmost duplicated base is
1-based p-L+1. The caller normalizer currently exports that left endpoint.
"""

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    if args.outdir.exists():
        raise FileExistsError(args.outdir)
    hashes, records, summary = {}, [], []

    def read(path):
        data = path.read_bytes()
        hashes[str(path)] = hashlib.sha256(data).hexdigest()
        return list(csv.DictReader(data.decode().splitlines(), delimiter="\t"))

    def family(value):
        return value.split("#", 1)[0].lower().replace("_", "")

    for dataset in ("mping", "ricetelib", "ricetelib_divergence"):
        for caller in ("relocate2", "relocate3-blat-bwaaln"):
            counts = Counter()
            base = args.benchmark / "reports/datasets" / dataset / "per_sample" / caller
            for sample in sorted(base.iterdir()):
                fps = read(sample / "false_positive_calls.tsv")
                used = set()
                for event in read(sample / "matches.tsv"):
                    length = int(event["tsd_length"])
                    if event["matched"] != "0" or length <= 11:
                        continue
                    left = int(event["position"]) - length + 1
                    choices = [i for i, call in enumerate(fps) if i not in used
                               and call["chrom"] == event["chrom"]
                               and family(call["te_family"]) == family(event["te_family"])
                               and int(call["position"]) == left]
                    if choices:
                        i = choices[0]
                        used.add(i)
                        counts[event["te_group"]] += 1
                        records.append(dict(dataset=dataset, caller=caller, sample=sample.name,
                                            event_id=event["event_id"], chrom=event["chrom"],
                                            truth_position=event["position"], tsd_length=length,
                                            expected_left=left, call_position=fps[i]["position"],
                                            family=event["te_family"], te_group=event["te_group"],
                                            truth_status=event["biological_class"], call_status=fps[i]["status"],
                                            exact_tsd=fps[i]["tsd"] == event["tsd"]))
            summary.append(dict(dataset=dataset, caller=caller, exact_anchor_pairs=sum(counts.values()),
                                LINE=counts["LINE"], SINE=counts["SINE"]))
    args.outdir.mkdir(parents=True, exist_ok=False)
    for name, data in (("exact_anchor_pairs", records), ("summary", summary)):
        with (args.outdir / (name + ".tsv")).open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(data)
    with (args.outdir / "provenance.json").open("x") as handle:
        json.dump({"inputs": hashes, "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   "rule": "Unmatched same-family call exactly at p-L+1; L>11; one-to-one within each sample"}, handle, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
