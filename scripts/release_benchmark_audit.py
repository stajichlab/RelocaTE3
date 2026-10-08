#!/usr/bin/env python3
"""Summarize existing R2/R3 benchmark tables, with explicit denominators.

Small-table analysis only. Never overwrite benchmark scores or existing output.
No statistical noninferiority claim: sample conditions reuse the same truth loci.
"""

import argparse
import csv
import hashlib
import json
import re
import statistics
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

CALLERS = ("relocate2", "relocate3-blat-bwaaln")


def load(path, hashes):
    data = path.read_bytes()
    hashes[str(path)] = hashlib.sha256(data).hexdigest()
    return list(csv.DictReader(data.decode().splitlines(), delimiter="\t"))


def status(value):
    return "somatic" if value.startswith("somatic") else value


def divide(n, d):
    return n / d if d else "NA"


def write(path, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.outdir.exists():
        raise FileExistsError(args.outdir)
    hashes, samples, events, fp_rows, paired, resource_rows = {}, [], [], [], [], []
    for dataset, expected_n in (("mping", 9), ("ricetelib", 9), ("ricetelib_divergence", 54)):
        base = args.benchmark / "reports/datasets" / dataset
        resources = {(r["caller"], r["sample"]): r for r in load(base / "resources.tsv", hashes)}
        per_caller = {}
        for caller in CALLERS:
            paths = sorted((base / "per_sample" / caller).iterdir())
            assert len(paths) == expected_n
            per_caller[caller] = {}
            for path in paths:
                sample = path.name
                coverage = int(re.search(r"cov(\d+)x", sample)[1])
                divergence = int(sample[3:6]) if sample.startswith("div") else 0
                meta = dict(dataset=dataset, caller=caller, sample=sample,
                            coverage=coverage, divergence=divergence)
                matches = load(path / "matches.tsv", hashes)
                fp = load(path / "false_positive_calls.tsv", hashes)
                precision = load(path / "precision.tsv", hashes)[0]
                assert len(matches) == 500
                tp = sum(r["matched"] == "1" for r in matches)
                assert tp == int(precision["matched_calls"])
                assert len(fp) == int(precision["false_positive_calls"])
                assert tp + len(fp) == int(precision["total_calls"])
                row = {**meta, "truth": len(matches), "tp": tp, "fp": len(fp),
                       "fn": len(matches) - tp, "precision": divide(tp, tp + len(fp)),
                       "recall": tp / len(matches), "f1": 2 * tp / (len(matches) + tp + len(fp))}
                resource = resources[(caller, sample)]
                row.update(wall_hours=float(resource["wall_seconds"]) / 3600,
                           rss_gib=float(resource["max_rss_kb"]) / 1048576)
                samples.append(row)
                per_caller[caller][sample] = (row, {m["event_id"]: m for m in matches})
                assert len(per_caller[caller][sample][1]) == 500
                for match in matches:
                    events.append({**meta, **match, "truth_status": status(match["biological_class"]),
                                   "predicted_status": status(match["call_status"]) if match["matched"] == "1" else "missed"})
                fp_rows.extend({**meta, **r, "caller": caller} for r in fp)
        assert per_caller[CALLERS[0]].keys() == per_caller[CALLERS[1]].keys()
        for sample, (a, ma) in per_caller[CALLERS[0]].items():
            b, mb = per_caller[CALLERS[1]][sample]
            assert ma.keys() == mb.keys()
            shared = [k for k in ma if ma[k]["matched"] == mb[k]["matched"] == "1"]
            paired.append(dict(dataset=dataset, sample=sample, coverage=a["coverage"],
                               divergence=a["divergence"], delta_tp=b["tp"] - a["tp"],
                               delta_fp=b["fp"] - a["fp"], delta_f1=b["f1"] - a["f1"],
                               shared_tp=len(shared), lost_tp=a["tp"] - len(shared),
                               gained_tp=b["tp"] - len(shared),
                               shared_same_position=sum(ma[k]["call_position"] == mb[k]["call_position"] for k in shared),
                               shared_same_status=sum(status(ma[k]["call_status"]) == status(mb[k]["call_status"]) for k in shared),
                               wall_ratio=b["wall_hours"] / a["wall_hours"],
                               rss_ratio=b["rss_gib"] / a["rss_gib"]))

    # Detection groups use truth denominators. No subgroup precision is inferred
    # by allocating unmatched calls to a true biological class or cell fraction.
    dimensions = [(), ("coverage",), ("divergence",), ("coverage", "divergence"),
                  ("truth_status",), ("truth_status", "coverage"),
                  ("truth_status", "cellular_fraction", "coverage"),
                  ("divergence", "truth_status", "cellular_fraction", "coverage"),
                  ("te_group",), ("divergence", "te_group")]
    grouped_rows = []
    for dims in dimensions:
        groups = defaultdict(list)
        for r in events:
            groups[(r["dataset"], r["caller"], *(r[d] for d in dims))].append(r)
        for key, group in sorted(groups.items()):
            detected = [r for r in group if r["matched"] == "1"]
            correct = sum(r["truth_status"] == r["predicted_status"] for r in detected)
            exact_tsd = sum(r["tsd_exact"] == "1" for r in detected)
            exact_position = sum(int(r["distance_bp"]) == 0 for r in detected)
            grouped_rows.append({"grouping": "+".join(dims) or "overall",
                                 "dataset": key[0], "caller": key[1], **dict(zip(dims, key[2:])),
                                 "truth": len(group), "tp": len(detected), "fn": len(group) - len(detected),
                                 "recall": len(detected) / len(group), "correct_status": correct,
                                 "status_accuracy_detected": divide(correct, len(detected)),
                                 "correct_status_recall": correct / len(group),
                                 "tsd_exact": exact_tsd, "tsd_accuracy_detected": divide(exact_tsd, len(detected)),
                                 "position_exact": exact_position,
                                 "position_accuracy_detected": divide(exact_position, len(detected))})

    pooled, confusion, somatic = [], [], []
    for dataset in ("mping", "ricetelib", "ricetelib_divergence"):
        for caller in CALLERS:
            selected = [s for s in samples if s["dataset"] == dataset and s["caller"] == caller]
            tp, fp, truth = (sum(s[k] for s in selected) for k in ("tp", "fp", "truth"))
            pooled.append(dict(dataset=dataset, caller=caller, samples=len(selected), tp=tp, fp=fp,
                               fn=truth-tp, precision=divide(tp, tp+fp), recall=tp/truth,
                               f1=2*tp/(truth+tp+fp),
                               median_wall_hours=statistics.median(s["wall_hours"] for s in selected),
                               max_wall_hours=max(s["wall_hours"] for s in selected),
                               median_rss_gib=statistics.median(s["rss_gib"] for s in selected),
                               max_rss_gib=max(s["rss_gib"] for s in selected)))
            es = [e for e in events if e["dataset"] == dataset and e["caller"] == caller]
            counts = Counter((e["truth_status"], e["predicted_status"]) for e in es)
            confusion.extend(dict(dataset=dataset, caller=caller, truth_status=a, predicted_status=b, count=n)
                             for (a, b), n in sorted(counts.items()))
            correct_somatic = counts[("somatic", "somatic")]
            predicted = sum(e["predicted_status"] == "somatic" for e in es)
            unmatched_somatic = sum(status(r["status"]) == "somatic" for r in fp_rows
                                    if r["dataset"] == dataset and r["caller"] == caller)
            somatic.append(dict(dataset=dataset, caller=caller, correctly_labeled_somatic=correct_somatic,
                                matched_calls_labeled_somatic=predicted, unmatched_calls_labeled_somatic=unmatched_somatic,
                                total_calls_labeled_somatic=predicted+unmatched_somatic,
                                somatic_label_precision=divide(correct_somatic, predicted+unmatched_somatic)))
    args.outdir.mkdir(parents=True, exist_ok=False)
    for name, rows in (("pooled", pooled), ("samples", samples), ("paired", paired),
                       ("truth_subgroups", grouped_rows), ("status_confusion", confusion), ("somatic_labels", somatic)):
        write(args.outdir / (name + ".tsv"), rows)
    with (args.outdir / "provenance.json").open("x") as handle:
        json.dump({"recorded_utc": datetime.now(timezone.utc).isoformat(),
                   "input_sha256": hashes, "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   "note": "Historical full benchmark, not the new pairing candidate; repeated loci, descriptive comparisons."}, handle, indent=2)
    (args.outdir / ".complete").touch()
    print(json.dumps(pooled, indent=2))


if __name__ == "__main__":
    main()
