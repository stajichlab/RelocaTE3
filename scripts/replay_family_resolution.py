#!/usr/bin/env python3
"""Freeze and compare original versus tied-family evidence on 72 stored runs.

Run preparation from the repository root; run tasks through SLURM. Both variants
use identical current source, isolating the mapping-evidence change from the
separately pending pairing change. Historical-call agreement is reported.
"""

import argparse
import ast
from collections import Counter
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, value):
    with path.open("x") as handle:
        json.dump(value, handle, indent=2)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare(args):
    import tomllib

    if args.output.exists():
        raise FileExistsError(args.output)
    cached = None
    if args.cached_evidence:
        cached = json.loads((args.cached_evidence / "manifest.json").read_text())
        for name, digest in cached["snapshot_sha256"].items():
            if sha(args.cached_evidence / name) != digest:
                raise RuntimeError(f"Cached replay source changed: {name}")
        if sha(Path("src/RelocaTE3/librelocate.py")) != sha(
            args.cached_evidence / "source/src/RelocaTE3/librelocate.py"
        ):
            raise RuntimeError("TE parser changed; cached evidence must be regenerated")
    config = tomllib.loads(
        (args.benchmark / "config/benchmark.full-aligners.toml").read_text()
    )
    tasks = []
    for dataset in ("mping", "ricetelib", "ricetelib_divergence"):
        cfg = config["datasets"][dataset]
        samples = sorted(
            (
                args.benchmark
                / "reports/datasets"
                / dataset
                / "per_sample/relocate3-blat-bwaaln"
            ).iterdir()
        )
        for directory in samples:
            sample = directory.name
            run = args.benchmark / "runs" / dataset / "relocate3-blat-bwaaln" / sample
            raw = run / "raw"
            inputs = {
                "genome_bam": raw / f"{sample}.repeat.bwaaln.sorted.bam",
                "read_repeat": raw / "te_containing" / f"{sample}.read_repeat_name.txt",
                "fullreads": raw / "genome_aln" / f"{sample}.fullreads.genome.bam",
                "characterization_bam": raw / f"{sample}.original_reads.sorted.bam",
                "left_te_bam": raw / f"{sample}.left.bam",
                "right_te_bam": raw / f"{sample}.right.bam",
                "annotation": args.benchmark / cfg["repeatmasker"],
                "truth": args.benchmark
                / "truth"
                / dataset
                / "per_sample"
                / f"{sample}.tsv",
                "historical_calls": run / "calls.normalized.tsv",
                "R2_calls": args.benchmark
                / "runs"
                / dataset
                / "relocate2"
                / sample
                / "calls.normalized.tsv",
            }
            if cached is not None:
                index = len(tasks)
                previous = cached["tasks"][index]
                directory = args.cached_evidence / "tasks" / f"{index:03d}"
                if (previous["dataset"], previous["sample"]) != (dataset, sample):
                    raise ValueError("Cached task order differs")
                if not (directory / ".complete").is_file():
                    raise RuntimeError(f"Incomplete cached task: {directory}")
                for name, path in inputs.items():
                    if [path.stat().st_size, path.stat().st_mtime_ns] != previous[
                        "input_stats"
                    ][name]:
                        raise RuntimeError(f"Cached evidence input changed: {path}")
                inputs["enriched_mapping"] = directory / "read_repeat.enriched.tsv"
                inputs["evidence_comparison"] = directory / "comparison.json"
            for key, path in inputs.items():
                if not path.is_file():
                    raise FileNotFoundError(path)
                if key.endswith("bam") or key == "fullreads":
                    if not Path(str(path) + ".bai").is_file():
                        raise FileNotFoundError(str(path) + ".bai")
            tasks.append(
                {
                    "dataset": dataset,
                    "sample": sample,
                    "te_name": cfg["te_name"],
                    "inputs": {key: str(path) for key, path in inputs.items()},
                    "input_stats": {
                        key: [path.stat().st_size, path.stat().st_mtime_ns]
                        for key, path in inputs.items()
                    },
                    "small_input_sha256": {
                        key: sha(inputs[key])
                        for key in ("truth", "historical_calls", "R2_calls")
                    },
                }
            )
    if len(tasks) != 72:
        raise ValueError(f"Expected 72 samples, found {len(tasks)}")
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    shutil.copytree(
        "src", out / "source/src", ignore=shutil.ignore_patterns("__pycache__", "*.pyc")
    )
    if args.historical_pairing:
        path = out / "source/src/RelocaTE3/insertions.py"
        current = path.read_text()
        historical = subprocess.check_output(
            ["git", "show", "HEAD:src/RelocaTE3/insertions.py"], text=True
        )

        def function_span(source):
            return next(
                node
                for node in ast.parse(source).body
                if isinstance(node, ast.FunctionDef)
                and node.name == "_pair_breakpoints"
            )

        old, new = function_span(historical), function_span(current)
        lines = current.splitlines(keepends=True)
        lines[new.lineno - 1 : new.end_lineno] = historical.splitlines(keepends=True)[
            old.lineno - 1 : old.end_lineno
        ]
        path.write_text("".join(lines))
    shutil.copy2(args.benchmark / "lib/calls.py", out / "calls.py")
    shutil.copy2(
        "validation/coordinate_scoring/scoring/score_calls.py", out / "score_calls.py"
    )
    shutil.copy2(__file__, out / "replay_family_resolution.py")
    shutil.copy2(
        "scripts/replay_family_resolution.slurm", out / "replay_family_resolution.slurm"
    )
    dump(
        out / "manifest.json",
        {
            "tasks": tasks,
            "snapshot_sha256": {
                str(p.relative_to(out)): sha(p) for p in out.rglob("*") if p.is_file()
            },
            "coordinate_policy": "tsd-interval",
            "created_unix": time.time(),
            "git_head": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], text=True
            ).strip(),
            "comparison": "identical current source, original versus enriched read-family mapping",
            "historical_pairing": args.historical_pairing,
            "cached_evidence_manifest_sha256": sha(
                args.cached_evidence / "manifest.json"
            )
            if cached
            else None,
        },
    )
    print(f"Prepared {len(tasks)} tasks at {out}")


def call_key(rows, family=True):
    fields = ("chrom", "position", "tsd", "strand", "status") + (
        ("te_family",) if family else ()
    )
    return Counter(tuple(str(row[k]) for k in fields) for row in rows)


def run(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Large stored-alignment scans require SLURM")
    import pysam
    from RelocaTE3.librelocate import RelocaTE

    root = args.output
    manifest = json.loads((root / "manifest.json").read_text())
    for relative, digest in manifest["snapshot_sha256"].items():
        if sha(root / relative) != digest:
            raise RuntimeError(f"Frozen source changed: {relative}")
    task = manifest["tasks"][args.task]
    inp = {key: Path(path) for key, path in task["inputs"].items()}
    for key, path in inp.items():
        if [path.stat().st_size, path.stat().st_mtime_ns] != task["input_stats"][key]:
            raise RuntimeError(f"Input changed: {path}")
    for key, digest in task["small_input_sha256"].items():
        if sha(inp[key]) != digest:
            raise RuntimeError(f"Input changed: {key}")
    out = root / "tasks" / f"{args.task:03d}"
    out.mkdir(parents=True, exist_ok=False)
    dump(
        out / "started.json",
        {"job": os.environ["SLURM_JOB_ID"], "task": task, "time": time.time()},
    )
    enriched, ambiguous, mapping_digest = obtain_mapping(inp, out, RelocaTE, pysam)
    calls = load_module("replay_calls", root / "calls.py")
    scorer = load_module("replay_score", root / "score_calls.py")
    replay_calls(
        args,
        manifest,
        task,
        inp,
        out,
        enriched,
        ambiguous,
        mapping_digest,
        calls,
        scorer,
    )


def obtain_mapping(inp, out, RelocaTE, pysam):
    if "enriched_mapping" in inp:
        comparison = json.loads(inp["evidence_comparison"].read_text())
        return (
            inp["enriched_mapping"],
            comparison["ambiguous_mapping_rows"],
            comparison["junction_mapping_sha256"],
        )
    names = set()
    with pysam.AlignmentFile(str(inp["genome_bam"]), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            if re.search(r":(start|end):[53]$", rec.query_name):
                names.add(re.sub(r":(start|end):[53]$", "", rec.query_name))
    hits = {}
    for side in ("left", "right"):
        selected = RelocaTE()._parse_te_bam(inp[side + "_te_bam"], read_names=names)
        for name, rec in selected.items():
            hits[name] = (
                rec if name not in hits else RelocaTE._merge_te_hits(hits[name], rec)
            )
        del selected
    missing = names - hits.keys()
    if missing:
        raise RuntimeError(f"Missing TE hits for {len(missing)} genome junction reads")
    enriched = out / "read_repeat.enriched.tsv"
    observed, ambiguous, digest = set(), 0, hashlib.sha256()
    with inp["read_repeat"].open("rb") as source, enriched.open("x") as dest:
        for line in source:
            digest.update(line)
            fields = line.decode().rstrip("\r\n").split("\t")
            if fields[0] in hits:
                rec = hits[fields[0]]
                if fields[1:3] != [rec["tName"], rec["strand"]]:
                    raise RuntimeError(
                        f"Selected hit differs from archive: {fields[0]}"
                    )
                observed.add(fields[0])
                ambiguous += int(len(rec.get("best_families", ())) > 1)
                RelocaTE._write_family_mapping(dest, fields[0], rec)
            else:
                dest.write(line.decode())
    if names != observed:
        raise RuntimeError("Incomplete junction family mapping")
    del hits, names
    return enriched, ambiguous, digest.hexdigest()


def replay_calls(
    args, manifest, task, inp, out, enriched, ambiguous, mapping_digest, calls, scorer
):
    root = args.output
    scores, normalized = {}, {}
    for variant, mapping in (
        ("original_evidence", inp["read_repeat"]),
        ("tied_evidence", enriched),
    ):
        work = out / variant
        work.mkdir()
        start = time.monotonic()
        env = dict(
            os.environ, PYTHONPATH=str(root / "source/src"), PYTHONDONTWRITEBYTECODE="1"
        )
        with (work / "pipeline.log").open("x") as log:

            def command(arguments):
                subprocess.run(
                    [sys.executable, "-m", "RelocaTE3", *arguments],
                    env=env,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                )

            command(
                [
                    "find-insertions",
                    "-b",
                    str(inp["genome_bam"]),
                    "--read-repeat",
                    str(mapping),
                    "--tsd",
                    "UNK",
                    "--target",
                    "ALL",
                    "--name",
                    task["sample"],
                    "--outdir",
                    str(work),
                    "--te-name",
                    task["te_name"],
                    "--reference-ins",
                    str(inp["annotation"]),
                    "--mismatch",
                    "2",
                    "--min-mapq",
                    "0",
                    "--require-both-junctions",
                    "--fullreads-bam",
                    str(inp["fullreads"]),
                ]
            )
            prefix = "ALL." + task["te_name"] + ".all_nonref_insert"
            command(
                [
                    "characterize",
                    "-s",
                    str(work / "results" / (prefix + ".txt")),
                    "-b",
                    str(inp["characterization_bam"]),
                    "-o",
                    str(work / "results"),
                ]
            )
        normalized[variant] = list(
            calls.parse_characterized_txt(
                work / "results" / (prefix + ".characTErized.txt"),
                caller=variant,
                sample=task["sample"],
            )
        )
        path = work / "calls.normalized.tsv"
        calls.write_normalized(normalized[variant], path)
        correctness, matches, fps, precision = scorer.score(
            inp["truth"],
            path,
            task["sample"],
            variant,
            10,
            coordinate_policy="tsd-interval",
        )
        for filename, rows in (
            ("correctness.tsv", correctness),
            ("matches.tsv", matches),
            ("false_positive_calls.tsv", fps),
            ("precision.tsv", [precision]),
        ):
            scorer._write(work / filename, rows)
        scores[variant] = dict(precision, elapsed_seconds=time.monotonic() - start)
    with inp["historical_calls"].open() as handle:
        historical = list(csv.DictReader(handle, delimiter="\t"))
    if manifest.get("historical_pairing") and call_key(
        normalized["original_evidence"]
    ) != call_key(historical):
        raise RuntimeError(
            "Restored pairing did not reproduce historical calls; attribution blocked"
        )
    r2work = out / "relocate2"
    r2work.mkdir()
    correctness, matches, fps, precision = scorer.score(
        inp["truth"],
        inp["R2_calls"],
        task["sample"],
        "relocate2",
        10,
        coordinate_policy="tsd-interval",
    )
    for filename, rows in (
        ("correctness.tsv", correctness),
        ("matches.tsv", matches),
        ("false_positive_calls.tsv", fps),
        ("precision.tsv", [precision]),
    ):
        scorer._write(r2work / filename, rows)
    scores["relocate2"] = precision
    dump(
        out / "comparison.json",
        {
            "dataset": task["dataset"],
            "sample": task["sample"],
            "scores": scores,
            "junction_mapping_sha256": mapping_digest,
            "ambiguous_mapping_rows": ambiguous,
            "original_evidence_matches_historical_calls": call_key(
                normalized["original_evidence"]
            )
            == call_key(historical),
            "geometry_and_genotypes_unchanged": call_key(
                normalized["original_evidence"], False
            )
            == call_key(normalized["tied_evidence"], False),
        },
    )
    (out / ".complete").touch()
    print(f"Completed {task['dataset']} {task['sample']}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "run"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--benchmark", type=Path)
    parser.add_argument("--task", type=int)
    parser.add_argument(
        "--cached-evidence",
        type=Path,
        help="Reuse mappings from a completed, frozen family replay",
    )
    parser.add_argument(
        "--historical-pairing",
        action="store_true",
        help="Restore only HEAD breakpoint pairing in the frozen snapshot",
    )
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
