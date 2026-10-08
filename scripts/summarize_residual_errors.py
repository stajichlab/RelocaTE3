#!/usr/bin/env python3
"""Summarize corrected error tables and diagnostic filter exposure, not new calls."""
import argparse
import json
import hashlib
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from audit_benchmark_errors import load, norm, write_table


def key(row):
    return row["chrom"], int(row["position"]), norm(row["te_family"])


def read_raw(path, hashes):
    data = path.read_bytes()
    hashes[str(path)] = hashlib.sha256(data).hexdigest()
    rows = []
    for line in data.decode().splitlines():
        if not line or line.startswith("#"):
            continue
        c = line.split("\t")
        rows.append(dict(chrom=c[3], position=c[4].split("..")[0],
                         te_family=c[0], tsd=c[1],
                         **dict(s.split(":", 1) for s in c[6:] if ":" in s)))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("audit", "reports", "benchmark", "outdir"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    if args.outdir.exists():
        raise FileExistsError(args.outdir)
    hashes, lost, annotated, exposure = {}, [], [], Counter()
    events = load(args.audit / "discordant_truth_events.tsv", hashes)
    fps = load(args.audit / "false_positive_audit.tsv", hashes)
    for event in events:
        if event["category"] != "R2_only_TP":
            continue
        near = json.loads(event["R3_calls_within_truth_interval"])
        raw = json.loads(event["nearest_R3_raw"] or "{}")
        classification = ("family_mismatch_at_locus" if near else "no_call_in_truth_window")
        lost.append({**event, "diagnosis": classification,
                     "nearby_family_status": raw.get("TE_family_status", "") if near else "",
                     "nearby_family_votes": raw.get("TE_family_support", "") if near else ""})
    caller = "relocate3-blat-bwaaln"
    for dataset in ("mping", "ricetelib", "ricetelib_divergence"):
        base = args.reports / "reports/datasets" / dataset / "per_sample" / caller
        for sample_dir in sorted(base.iterdir()):
            sample = sample_dir.name
            calls = load(args.reports / "runs" / dataset / caller / sample / "calls.normalized.tsv", hashes)
            fp = load(sample_dir / "false_positive_calls.tsv", hashes)
            fp_keys = Counter(key(r) for r in fp)
            directory = args.benchmark / "runs" / dataset / caller / sample / "raw/results"
            # The characterized calls derive from .raw.txt, not always .all.txt.
            paths = sorted(directory.glob("*.all_nonref_insert.raw.txt"))
            assert len(paths) == 1, (dataset, sample, paths)
            raw = read_raw(paths[0], hashes)
            lookup = {}
            for r in raw:
                lookup.setdefault(key(r), []).append(r)
            for call in calls:
                k = key(call)
                choices = [r for r in lookup.get(k, []) if r["tsd"] == call["tsd"]]
                assert len(choices) == 1, (dataset, sample, call, choices)
                r = choices[0]
                state = "FP" if fp_keys[k] else "TP"
                if state == "FP":
                    fp_keys[k] -= 1
                rules = {
                    "all_calls": True,
                    "ambiguous_family": r["TE_family_status"] == "ambiguous",
                    "junction_count_lt3": int(r["T"]) < 3,
                    "one_sided": min(int(r["L"]), int(r["R"])) == 0,
                    "junction_support_family_discordant": r["TE_family_concordance"] == "discordant",
                }
                for rule, applies in rules.items():
                    if applies:
                        exposure[dataset, rule, state] += 1
                if state == "FP":
                    audit = [f for f in fps if f["dataset"] == dataset and f["sample"] == sample
                             and f["side"] == "R3" and key(f) == k]
                    assert len(audit) == 1
                    annotated.append({**audit[0], "exact_raw_evidence": json.dumps(r, sort_keys=True)})
            assert not any(fp_keys.values())
    args.outdir.mkdir(parents=True, exist_ok=False)
    write_table(args.outdir / "lost_events.tsv", lost)
    write_table(args.outdir / "R3_false_positives.tsv", annotated)
    write_table(args.outdir / "filter_exposure.tsv", [
        dict(dataset=d, proposed_exclusion=rule, matched_calls_excluded=exposure[d, rule, "TP"],
             false_positive_calls_excluded=exposure[d, rule, "FP"])
        for d, rule in sorted({(d, rule) for d, rule, _ in exposure})])
    (args.outdir / "provenance.json").write_text(json.dumps({
        "recorded_utc": datetime.now(timezone.utc).isoformat(), "input_sha256": hashes,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "note": "Exclusion exposure only: no new calls or rematching; not a validated filter."}, indent=2))
    (args.outdir / ".complete").touch()
    print(f"Annotated {len(lost)} R2-only events and {len(annotated)} R3 FPs")


if __name__ == "__main__":
    main()
