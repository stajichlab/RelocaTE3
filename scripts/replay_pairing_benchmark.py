#!/usr/bin/env python3
"""Freeze and replay step 5/characterization on stored benchmark alignments.

Preparation is lightweight. Run tasks only through SLURM. Each task first
reproduces historical normalized calls using the pre-patch insertion finder;
a mismatch fails the task instead of attributing upstream drift to this patch.
"""

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def dump(path, data):
    with path.open("x") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")


def prepare(args):
    import tomllib

    out = args.output
    if out.exists():
        raise FileExistsError(out)
    benchmark = args.benchmark.resolve()
    config_path = benchmark / "config/benchmark.full-aligners.toml"
    config = tomllib.loads(config_path.read_text())
    tasks = []
    for dataset in ("mping", "ricetelib", "ricetelib_divergence"):
        cfg = config["datasets"][dataset]
        for path in sorted((benchmark / "reports/datasets" / dataset / "per_sample/relocate3-blat-bwaaln").iterdir()):
            sample = path.name
            run = benchmark / "runs" / dataset / "relocate3-blat-bwaaln" / sample
            raw = run / "raw"
            inputs = dict(
                genome_bam=raw / f"{sample}.repeat.bwaaln.sorted.bam",
                read_repeat=raw / "te_containing" / f"{sample}.read_repeat_name.txt",
                fullreads=raw / "genome_aln" / f"{sample}.fullreads.genome.bam",
                characterization_bam=raw / f"{sample}.original_reads.sorted.bam",
                annotation=benchmark / cfg["repeatmasker"],
                truth=benchmark / "truth" / dataset / "per_sample" / (sample + ".tsv"),
                historical_calls=run / "calls.normalized.tsv",
            )
            for name, p in inputs.items():
                if not p.is_file():
                    raise FileNotFoundError(p)
                if name.endswith("bam") or name == "fullreads":
                    if not Path(str(p) + ".bai").is_file():
                        raise FileNotFoundError(str(p) + ".bai")
            tasks.append(dict(dataset=dataset, sample=sample, te_name=cfg["te_name"],
                              inputs={k: str(p) for k, p in inputs.items()},
                              input_stats={k: [p.stat().st_size, p.stat().st_mtime_ns] for k, p in inputs.items()},
                              historical_calls_sha256=sha(inputs["historical_calls"])))
    assert len(tasks) == 72
    out.mkdir(parents=True)
    for variant in ("baseline", "candidate"):
        shutil.copytree("src", out / variant / "src", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    # Before this change insertions.py was clean at HEAD; preserve all other
    # current code identically, including the validated memory improvements.
    baseline = subprocess.check_output(["git", "show", "HEAD:src/RelocaTE3/insertions.py"])
    (out / "baseline/src/RelocaTE3/insertions.py").write_bytes(baseline)
    shutil.copytree("tests", out / "tests", ignore=shutil.ignore_patterns("__pycache__", "results", "*.pyc"))
    for relative in ("lib/calls.py", "scoring/score_calls.py", "config/benchmark.full-aligners.toml"):
        destination = out / "benchmark" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(benchmark / relative, destination)
    shutil.copy2(__file__, out / "replay_pairing_benchmark.py")
    shutil.copy2("scripts/replay_pairing_benchmark.slurm", out / "replay_pairing_benchmark.slurm")
    hashes = {str(p.relative_to(out)): sha(p) for p in out.rglob("*") if p.is_file()}
    dump(out / "manifest.json", {"tasks": tasks, "snapshot_sha256": hashes,
                                 "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
                                 "baseline_note": "HEAD insertion finder; all other source identical to candidate",
                                 "created_unix": time.time()})
    print(f"Prepared {len(tasks)} tasks at {out}")


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def call_key(rows_):
    return Counter(tuple(str(r[k]) for k in ("chrom", "position", "te_family", "tsd", "strand", "status")) for r in rows_)


def run(args):
    root = args.output.resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    for relative, expected in manifest["snapshot_sha256"].items():
        if sha(root / relative) != expected:
            raise RuntimeError(f"Snapshot changed: {relative}")
    task = manifest["tasks"][args.task]
    destination = root / "tasks" / f"{args.task:03d}"
    destination.mkdir(parents=True, exist_ok=False)
    dump(destination / "started.json", {"task": task, "job_id": os.environ.get("SLURM_JOB_ID"), "time": time.time()})
    inp = {k: Path(p) for k, p in task["inputs"].items()}
    for key, path in inp.items():
        if [path.stat().st_size, path.stat().st_mtime_ns] != task["input_stats"][key]:
            raise RuntimeError(f"Input changed since preparation: {path}")
    assert sha(inp["historical_calls"]) == task["historical_calls_sha256"]
    calls_module = load_module("replay_calls", root / "benchmark/lib/calls.py")
    scorer = load_module("replay_score", root / "benchmark/scoring/score_calls.py")
    scores = {}
    for variant in ("baseline", "candidate"):
        work = destination / variant
        work.mkdir()
        env = dict(os.environ, PYTHONPATH=str(root / variant / "src"), PYTHONDONTWRITEBYTECODE="1")
        start = time.monotonic()
        with (work / "pipeline.log").open("x") as log:
            def command(arguments):
                print(json.dumps(arguments), file=log, flush=True)
                subprocess.run([sys.executable, "-m", "RelocaTE3", *arguments], env=env,
                               stdout=log, stderr=subprocess.STDOUT, check=True)
            command(["find-insertions", "-b", str(inp["genome_bam"]), "--read-repeat", str(inp["read_repeat"]),
                     "--tsd", "UNK", "--target", "ALL", "--name", task["sample"], "--outdir", str(work),
                     "--te-name", task["te_name"], "--reference-ins", str(inp["annotation"]),
                     "--mismatch", "2", "--min-mapq", "0", "--require-both-junctions",
                     "--fullreads-bam", str(inp["fullreads"])])
            prefix = "ALL." + task["te_name"] + ".all_nonref_insert"
            command(["characterize", "-s", str(work / "results" / (prefix + ".txt")),
                     "-b", str(inp["characterization_bam"]), "-o", str(work / "results")])
        normalized = list(calls_module.parse_characterized_txt(
            work / "results" / (prefix + ".characTErized.txt"), caller=variant, sample=task["sample"]))
        calls_module.write_normalized(normalized, work / "calls.normalized.tsv")
        if variant == "baseline" and call_key(normalized) != call_key(rows(inp["historical_calls"])):
            raise RuntimeError("Baseline does not reproduce historical calls; candidate comparison blocked")
        summary, matches, fps, precision = scorer.score(inp["truth"], work / "calls.normalized.tsv",
                                                       task["sample"], variant, 10)
        for filename, records in (("correctness.tsv", summary), ("matches.tsv", matches),
                                  ("false_positive_calls.tsv", fps), ("precision.tsv", [precision])):
            scorer._write(work / filename, records)
        scores[variant] = {**precision, "elapsed_seconds": time.monotonic() - start}
    dump(destination / "comparison.json", {"dataset": task["dataset"], "sample": task["sample"],
                                           "baseline_matches_history": True, "scores": scores})
    (destination / ".complete").touch()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "run"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--benchmark", type=Path)
    parser.add_argument("--task", type=int)
    args = parser.parse_args()
    if args.mode == "prepare":
        if args.benchmark is None:
            parser.error("prepare requires --benchmark")
        prepare(args)
    else:
        if args.task is None or not 0 <= args.task < 72:
            parser.error("run requires --task in 0..71")
        run(args)


if __name__ == "__main__":
    main()
