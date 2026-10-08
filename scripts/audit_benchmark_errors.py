#!/usr/bin/env python3
"""Read-only, small-table audit of R2 versus R3 BLAT/bwa-aln errors.

Run from the repository root. This does not scan BAMs, rerun callers, or change
benchmark scores. Caller-to-caller FP pairing uses exact family/position first,
then distance-ordered one-to-one family-compatible pairs within 10 bp.
"""

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def norm(value):
    return value.split("#", 1)[0].lower().replace("_", "")


def load(path, hashes):
    data = path.read_bytes()
    hashes[str(path)] = hashlib.sha256(data).hexdigest()
    return list(csv.DictReader(data.decode().splitlines(), delimiter="\t"))


def nearby(row, candidates, family=False):
    valid = [r for r in candidates if r["chrom"] == row["chrom"]
             and (not family or norm(r["te_family"]) == norm(row["te_family"]))]
    return min(valid, key=lambda r: abs(int(r["position"]) - int(row["position"])),
               default=None)


def describe(row, other):
    if other is None:
        return ""
    return json.dumps({**other, "offset_bp": int(other["position"]) - int(row["position"])},
                      sort_keys=True)


def truth_distance(truth, call):
    """Use the scored truth interval when present; otherwise legacy anchor."""
    start = int(truth.get("truth_match_start") or truth["position"])
    end = int(truth.get("truth_match_end") or truth["position"])
    position = int(call["position"])
    return max(start - position, position - end, 0)


def near_interval(truth, calls):
    return [dict(row, interval_distance_bp=truth_distance(truth, row))
            for row in calls if row["chrom"] == truth["chrom"]
            and truth_distance(truth, row) <= 10]


def pair_fps(left, right):
    edges = sorted((abs(int(a["position"]) - int(b["position"])), i, j)
                   for i, a in enumerate(left) for j, b in enumerate(right)
                   if a["chrom"] == b["chrom"]
                   and norm(a["te_family"]) == norm(b["te_family"])
                   and abs(int(a["position"]) - int(b["position"])) <= 10)
    used_left, used_right = set(), set()
    for _, i, j in edges:
        if i not in used_left and j not in used_right:
            used_left.add(i)
            used_right.add(j)
    return used_left, used_right


def raw_calls(directory, hashes):
    paths = sorted(directory.glob("*.all_nonref_insert.all.txt"))
    if len(paths) != 1:
        raise ValueError(f"Expected one raw call table in {directory}: {paths}")
    path = paths[0]
    data = path.read_bytes()
    hashes[str(path)] = hashlib.sha256(data).hexdigest()
    result = []
    for line in data.decode().splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        start, end = map(int, cols[4].split(".."))
        result.append({"chrom": cols[3], "position": str(start), "end": str(end),
                       "te_family": cols[0], "tsd": cols[1],
                       **dict(c.split(":", 1) for c in cols[6:] if ":" in c)})
    return result


def write_table(path, records):
    fields = list(dict.fromkeys(k for r in records for k in r))
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--raw-benchmark", type=Path,
                        help="Original benchmark containing raw evidence, if reports are rescored")
    args = parser.parse_args()
    if args.outdir.exists():
        raise FileExistsError(f"Refusing existing output: {args.outdir}")
    hashes, events, fps, summaries = {}, [], [], []
    for dataset in ("mping", "ricetelib", "ricetelib_divergence"):
        base = args.benchmark / "reports/datasets" / dataset / "per_sample"
        counts = Counter()
        left_samples = sorted(p.name for p in (base / "relocate2").iterdir() if p.is_dir())
        right_samples = sorted(p.name for p in (base / "relocate3-blat-bwaaln").iterdir() if p.is_dir())
        assert left_samples == right_samples, dataset
        for sample in left_samples:
            tables, calls = [], []
            for caller in ("relocate2", "relocate3-blat-bwaaln"):
                tables.append({name: load(base / caller / sample / (name + ".tsv"), hashes)
                               for name in ("matches", "false_positive_calls", "precision")})
                calls.append(load(args.benchmark / "runs" / dataset / caller / sample
                                  / "calls.normalized.tsv", hashes))
            raw = raw_calls((args.raw_benchmark or args.benchmark) / "runs" / dataset / "relocate3-blat-bwaaln"
                            / sample / "raw/results", hashes)
            a, b = [{r["event_id"]: r for r in t["matches"]} for t in tables]
            assert a.keys() == b.keys() and len(a) == len(tables[0]["matches"])
            for key in a:
                x, y = a[key], b[key]
                for field in ("chrom", "position", "te_family"):
                    assert x[field] == y[field], (dataset, sample, key, field)
                state = ("shared_TP" if x["matched"] == y["matched"] == "1" else
                         "R2_only_TP" if x["matched"] == "1" else
                         "R3_only_TP" if y["matched"] == "1" else "shared_FN")
                counts[state] += 1
                if state in ("R2_only_TP", "R3_only_TP"):
                    events.append({"dataset": dataset, "sample": sample, "category": state,
                                   **x, "R3_matched": y["matched"],
                                   "R3_call_position": y.get("call_position", ""),
                                   "R3_calls_within_truth_interval": json.dumps(near_interval(x, calls[1]), sort_keys=True),
                                   "nearest_R3_call": describe(x, nearby(x, calls[1])),
                                   "nearest_same_family_R3_call": describe(x, nearby(x, calls[1], True)),
                                   "nearest_R3_raw": describe(x, nearby(x, raw))})
            fa, fb = [t["false_positive_calls"] for t in tables]
            ua, ub = pair_fps(fa, fb)
            counts["shared_FP"] += len(ua)
            counts["R2_only_FP"] += len(fa) - len(ua)
            counts["R3_only_FP"] += len(fb) - len(ub)
            for side, rows, used in ((0, fa, ua), (1, fb, ub)):
                for i, row in enumerate(rows):
                    near_truth = min((t for t in a.values() if t["chrom"] == row["chrom"]),
                                     key=lambda t: truth_distance(t, row), default=None)
                    near_other = nearby(row, calls[1 - side])
                    same_family = nearby(row, calls[1 - side], True)
                    annotation = "no_truth_within_10bp"
                    if near_truth and truth_distance(near_truth, row) <= 10:
                        annotation = ("same_family_truth_within_10bp" if
                                      norm(row["te_family"]) == norm(near_truth["te_family"])
                                      else "different_family_truth_within_10bp")
                    fps.append({"dataset": dataset, "sample": sample,
                                "side": "R2" if side == 0 else "R3",
                                "category": "shared_FP" if i in used else "caller_only_FP",
                                **row, "truth_context": annotation,
                                "truth_interval_distance_bp": truth_distance(near_truth, row) if near_truth else "",
                                "nearest_truth": describe(row, near_truth),
                                "nearest_other_call": describe(row, near_other),
                                "nearest_same_family_other_call": describe(row, same_family),
                                "nearest_R3_raw": describe(row, nearby(row, raw))})
            for side, table in enumerate(tables):
                pr = table["precision"][0]
                assert int(pr["false_positive_calls"]) == len(table["false_positive_calls"])
                assert int(pr["matched_calls"]) == sum(r["matched"] == "1" for r in table["matches"])
                assert int(pr["total_calls"]) == len(calls[side])
        summaries.append({"dataset": dataset, "samples": len(left_samples),
                          **{k: counts[k] for k in ("shared_TP", "R2_only_TP", "R3_only_TP",
                                                   "shared_FN", "shared_FP", "R2_only_FP", "R3_only_FP")}})
        print(json.dumps(summaries[-1]), flush=True)
    args.outdir.mkdir(parents=True, exist_ok=False)
    write_table(args.outdir / "summary.tsv", summaries)
    write_table(args.outdir / "discordant_truth_events.tsv", events)
    write_table(args.outdir / "false_positive_audit.tsv", fps)
    with (args.outdir / "provenance.json").open("x") as handle:
        json.dump({"recorded_utc": datetime.now(timezone.utc).isoformat(),
                   "benchmark": str(args.benchmark), "window_bp": 10,
                   "raw_benchmark": str(args.raw_benchmark or args.benchmark),
                   "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   "input_sha256": hashes}, handle, indent=2)
    (args.outdir / ".complete").touch()


if __name__ == "__main__":
    main()
