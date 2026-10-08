#!/usr/bin/env python3
"""Rescore small archived call tables without rerunning or changing callers.

Legacy results must reproduce first. Writes an isolated two-caller report suite;
the benchmark dashboard and historical reports are never overwritten.
"""

import argparse
import csv
import hashlib
import importlib.util
import json
import shutil
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


def load(path, hashes):
    data = path.read_bytes()
    hashes[str(path)] = hashlib.sha256(data).hexdigest()
    return list(csv.DictReader(data.decode().splitlines(), delimiter="\t"))


def same_rows(expected, actual):
    if not expected:
        return not actual
    keys = sorted(expected[0])
    def key(row):
        return tuple(str(row.get(k, "")) for k in keys)
    return Counter(map(key, expected)) == Counter(map(key, actual))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    if args.outdir.exists():
        raise FileExistsError(args.outdir)
    scorer_path = Path("validation/coordinate_scoring/scoring/score_calls.py")
    spec = importlib.util.spec_from_file_location("coordinate_scorer", scorer_path)
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    args.outdir.mkdir(parents=True, exist_ok=False)
    shutil.copy2(scorer_path, args.outdir / "score_calls.py")
    hashes, changes, sample_summary = {}, [], []
    for dataset, n in (("mping", 9), ("ricetelib", 9), ("ricetelib_divergence", 54)):
        aggregate = defaultdict(list)
        source = args.benchmark / "reports/datasets" / dataset
        for caller in ("relocate2", "relocate3-blat-bwaaln"):
            samples = sorted(p for p in (source / "per_sample" / caller).iterdir() if p.is_dir())
            assert len(samples) == n
            for sample_path in samples:
                sample = sample_path.name
                calls = args.benchmark / "runs" / dataset / caller / sample / "calls.normalized.tsv"
                truth = args.benchmark / "truth" / dataset / "per_sample" / (sample + ".tsv")
                load(calls, hashes)
                load(truth, hashes)
                before = scorer.score(truth, calls, sample, caller, 10, "legacy")
                names = ("correctness.tsv", "matches.tsv", "false_positive_calls.tsv", "precision.tsv")
                for name, records in zip(names, (*before[:3], [before[3]])):
                    historical = load(sample_path / name, hashes)
                    if not same_rows(historical, records):
                        raise RuntimeError(f"Legacy reproduction failed: {dataset}/{caller}/{sample}/{name}")
                after = scorer.score(truth, calls, sample, caller, 10, "tsd-interval")
                destination = args.outdir / "reports/datasets" / dataset / "per_sample" / caller / sample
                destination.mkdir(parents=True)
                for name, records in zip(names, (*after[:3], [after[3]])):
                    scorer._write(destination / name, records)
                    if name in ("correctness.tsv", "precision.tsv"):
                        aggregate[name].extend(records)
                (destination / ".complete").touch()
                # The release-audit reader expects the unchanged normalized
                # input calls and resources beside this isolated report tree.
                copied = args.outdir / "runs" / dataset / caller / sample / "calls.normalized.tsv"
                copied.parent.mkdir(parents=True)
                shutil.copy2(calls, copied)
                b = {r["event_id"]: r for r in before[1]}
                a = {r["event_id"]: r for r in after[1]}
                assert a.keys() == b.keys()
                for event in a:
                    old, new = b[event], a[event]
                    if old["matched"] != new["matched"] or any(old.get(k) != new.get(k)
                            for k in ("call_position", "call_tsd", "call_status")):
                        changes.append(dict(dataset=dataset, caller=caller, sample=sample,
                                            event_id=event, te_group=new.get("te_group", ""),
                                            truth_position=new["position"], truth_tsd_length=new.get("tsd_length", ""),
                                            old_matched=old["matched"], new_matched=new["matched"],
                                            old_call_position=old.get("call_position", ""),
                                            new_call_position=new.get("call_position", ""),
                                            old_status=old.get("call_status", ""), new_status=new.get("call_status", "")))
                sample_summary.append(dict(dataset=dataset, caller=caller, sample=sample,
                                           old_tp=before[3]["matched_calls"], new_tp=after[3]["matched_calls"],
                                           old_fp=before[3]["false_positive_calls"], new_fp=after[3]["false_positive_calls"],
                                           lost_events=sum(b[k]["matched"] == "1" and a[k]["matched"] == "0" for k in a)))
        report = args.outdir / "reports/datasets" / dataset
        for name, records in aggregate.items():
            scorer._write(report / name, records)
        resources = load(source / "resources.tsv", hashes)
        scorer._write(report / "resources.tsv", [r for r in resources if r["caller"] in ("relocate2", "relocate3-blat-bwaaln")])
        print(f"Completed {dataset}: {n * 2} legacy-verified sample rescoring runs", flush=True)
    scorer._write(args.outdir / "changed_events.tsv", changes)
    scorer._write(args.outdir / "sample_changes.tsv", sample_summary)
    with (args.outdir / "provenance.json").open("x") as handle:
        json.dump(dict(recorded_utc=datetime.now(timezone.utc).isoformat(), policy="tsd-interval-v1",
                       window=10, input_sha256=hashes,
                       scorer_sha256=hashlib.sha256(scorer_path.read_bytes()).hexdigest(),
                       runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                       samples=144, legacy_reproduction=True,
                       scope="Historical full benchmark, not the new pairing-correction candidate"), handle, indent=2)
    (args.outdir / ".complete").touch()


if __name__ == "__main__":
    main()
