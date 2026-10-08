#!/usr/bin/env python3
"""Test the R2 intact-read veto without modifying production calling code.

Preparation freezes the already validated source. Compute tasks reproduce the
stabilized baseline, then add the intact-read veto and rescore against truth.
The existing spanning-read veto is retained: this is an additive candidate,
not a claim to reproduce every R2 full-read filtering detail.
"""

import argparse
from collections import Counter
import contextlib
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, data):
    with path.open("x") as handle:
        json.dump(data, handle, indent=2)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def is_intact(record):
    """R2 read_junction_reads_align: at most ten query bases unaligned.

    Insertions consume query and count as aligned in R2. Extended =/X CIGAR
    operations are treated as M. Neither MAPQ nor proper-pair flag is a gate.
    """
    if record.is_unmapped or record.query_sequence is None:
        return False
    aligned = sum(
        length for op, length in record.cigartuples or () if op in (0, 1, 7, 8)
    )
    return aligned >= len(record.query_sequence) - 10


def veto(ins, intact, key, basekey):
    def found(name):
        k = key(name)
        return k in intact or basekey(k) in intact

    left = sum(found(n) for n in ins.read_names[: ins.left_junction_reads])
    right = sum(found(n) for n in ins.read_names[ins.left_junction_reads :])
    reject = (
        ins.left_junction_reads + ins.right_junction_reads > 0
        and left >= 0.3 * ins.left_junction_reads
        and right >= 0.3 * ins.right_junction_reads
    )
    return reject, left, right


def local_intact(bam, call, module):
    """Use the existing 500-bp candidate fetch window, retaining mate identity.

    Local means overlap with this same-contig, zero-based half-open window;
    it does not require crossing the breakpoint. No window tuning by family.
    """
    if bam is None:
        return set()
    lo = max(0, min(call.start, call.end) - module.FULLREAD_WINDOW)
    hi = max(call.start, call.end) + module.FULLREAD_WINDOW
    try:
        records = bam.fetch(call.chrom, lo, hi)
    except (ValueError, KeyError):
        return set()
    intact = set()
    for record in records:
        if record.query_sequence is None:
            continue
        key = module._fullread_record_key(record)
        if is_intact(record):
            intact.add(key)
        else:
            intact.discard(key)
    return intact


def prepare(args):
    old = json.loads((args.baseline / "manifest.json").read_text())
    for name, digest in old["snapshot_sha256"].items():
        if sha(args.baseline / name) != digest:
            raise ValueError(f"Frozen baseline changed: {name}")
    tasks = []
    for index, previous in enumerate(old["tasks"]):
        folder = args.baseline / "tasks" / f"{index:03d}"
        if not (folder / ".complete").exists():
            raise RuntimeError(f"Baseline incomplete: {index}")
        task = dict(previous)
        task["inputs"] = dict(
            previous["inputs"],
            stabilized_calls=str(folder / "tied_evidence/calls.normalized.tsv"),
        )
        task["input_stats"] = {}
        for key, value in task["inputs"].items():
            p = Path(value)
            stat = [p.stat().st_size, p.stat().st_mtime_ns]
            if key in previous["input_stats"] and stat != previous["input_stats"][key]:
                raise ValueError(f"Baseline input changed: {p}")
            task["input_stats"][key] = stat
        task["small_input_sha256"] = dict(
            previous["small_input_sha256"],
            stabilized_calls=sha(Path(task["inputs"]["stabilized_calls"])),
        )
        tasks.append(task)
    if len(tasks) != 72:
        raise ValueError("Expected 72 samples")
    args.output.mkdir(parents=True, exist_ok=False)
    shutil.copytree(args.baseline / "source", args.output / "source")
    for name in ("calls.py", "score_calls.py"):
        shutil.copy2(args.baseline / name, args.output / name)
    shutil.copy2(__file__, args.output / "replay_intact_fullreads.py")
    shutil.copy2(
        "scripts/replay_intact_fullreads.slurm",
        args.output / "replay_intact_fullreads.slurm",
    )
    dump(
        args.output / "manifest.json",
        dict(
            tasks=tasks,
            scope=args.scope,
            baseline_manifest_sha256=sha(args.baseline / "manifest.json"),
            created_unix=time.time(),
            snapshot_sha256={
                str(p.relative_to(args.output)): sha(p)
                for p in args.output.rglob("*")
                if p.is_file() and "__pycache__" not in p.parts
            },
        ),
    )
    print(f"Prepared {len(tasks)} tasks at {args.output}")


def run(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Full-alignment scans and calling require SLURM")
    import pysam
    from RelocaTE3 import insertions as ins
    from RelocaTE3.cli import main as cli

    root = args.output
    manifest = json.loads((root / "manifest.json").read_text())
    scope = manifest.get("scope", "global")
    if scope not in ("global", "local"):
        raise ValueError(f"Invalid frozen scope: {scope}")
    for name, digest in manifest["snapshot_sha256"].items():
        if sha(root / name) != digest:
            raise ValueError(f"Frozen source changed: {name}")
    task = manifest["tasks"][args.task]
    inp = {k: Path(v) for k, v in task["inputs"].items()}
    for key, p in inp.items():
        if [p.stat().st_size, p.stat().st_mtime_ns] != task["input_stats"][key]:
            raise ValueError(f"Input changed: {p}")
    for key, digest in task["small_input_sha256"].items():
        if sha(inp[key]) != digest:
            raise ValueError(f"Small input changed: {key}")
    out = root / "tasks" / f"{args.task:03d}"
    out.mkdir(parents=True, exist_ok=False)
    dump(
        out / "started.json",
        dict(job=os.environ["SLURM_JOB_ID"], task=task, time=time.time()),
    )
    intact, counts = set(), Counter()
    if scope == "global":
        wanted = set()
        with pysam.AlignmentFile(str(inp["genome_bam"]), "rb") as bam:
            for rec in bam.fetch(until_eof=True):
                if ins._JUNCTION_RE.search(rec.query_name):
                    k = ins._junction_fullread_key(rec.query_name)
                    wanted.update((k, ins._fullread_key(k)))
        with pysam.AlignmentFile(str(inp["fullreads"]), "rb") as bam:
            for rec in bam.fetch(until_eof=True):
                k = ins._fullread_record_key(rec)
                if k not in wanted or rec.query_sequence is None:
                    continue
                counts["selected_records"] += 1
                counts["secondary_or_supplementary"] += int(
                    rec.is_secondary or rec.is_supplementary
                )
                # R2 overwrites the stored per-mate record in input order.
                if is_intact(rec):
                    intact.add(k)
                else:
                    intact.discard(k)
        counts["intact_keys"] = len(intact)
        del wanted
    calls = load_module("intact_calls", root / "calls.py")
    scorer = load_module("intact_score", root / "score_calls.py")
    original = ins._fullread_false_junction
    decisions = []

    def candidate(bam, call):
        old = original(bam, call)
        new, left, right = veto(
            call,
            local_intact(bam, call, ins) if scope == "local" else intact,
            ins._junction_fullread_key,
            ins._fullread_key,
        )
        if new and not old:
            decisions.append(
                dict(
                    chrom=call.chrom,
                    start=call.start,
                    end=call.end,
                    family=call.te_name,
                    tsd=call.tsd,
                    left_total=call.left_junction_reads,
                    right_total=call.right_junction_reads,
                    left_intact=left,
                    right_intact=right,
                    read_names=call.read_names,
                )
            )
        return old or new

    scores = {}
    for variant in ("baseline", "intact_candidate"):
        ins._fullread_false_junction = original if variant == "baseline" else candidate
        work = out / variant
        work.mkdir()
        prefix = f"ALL.{task['te_name']}.all_nonref_insert"
        with (
            (work / "pipeline.log").open("x") as log,
            contextlib.redirect_stdout(log),
            contextlib.redirect_stderr(log),
        ):
            rc = cli(
                [
                    "find-insertions",
                    "-b",
                    str(inp["genome_bam"]),
                    "--read-repeat",
                    str(inp["enriched_mapping"]),
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
            if rc:
                raise RuntimeError(f"find-insertions exit {rc}")
            rc = cli(
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
            if rc:
                raise RuntimeError(f"characterize exit {rc}")
        normalized = list(
            calls.parse_characterized_txt(
                work / "results" / (prefix + ".characTErized.txt"),
                caller=variant,
                sample=task["sample"],
            )
        )
        path = work / "calls.normalized.tsv"
        calls.write_normalized(normalized, path)
        if variant == "baseline":
            with inp["stabilized_calls"].open() as handle:
                expected = list(csv.DictReader(handle, delimiter="\t"))
            fields = ("chrom", "position", "te_family", "tsd", "strand", "status")

            def keys(rows):
                return Counter(tuple(str(row[k]) for k in fields) for row in rows)

            if keys(normalized) != keys(expected):
                raise RuntimeError("Baseline does not reproduce stabilized calls")
        correctness, matches, fps, precision = scorer.score(
            inp["truth"],
            path,
            task["sample"],
            variant,
            10,
            coordinate_policy="tsd-interval",
        )
        for name, rows in (
            ("correctness.tsv", correctness),
            ("matches.tsv", matches),
            ("false_positive_calls.tsv", fps),
            ("precision.tsv", [precision]),
        ):
            scorer._write(work / name, rows)
        scores[variant] = precision
    ins._fullread_false_junction = original
    dump(out / "additional_vetoes.json", decisions)
    dump(
        out / "comparison.json",
        dict(
            dataset=task["dataset"],
            sample=task["sample"],
            scores=scores,
            evidence_counts=dict(counts),
            baseline_reproduced=True,
            scope=scope,
        ),
    )
    (out / ".complete").touch()
    print(f"Completed {task['dataset']} {task['sample']}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "run"))
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--task", type=int)
    parser.add_argument(
        "--scope",
        choices=("global", "local"),
        default="global",
        help="Preparation only: freeze the intact-read evidence scope",
    )
    args = parser.parse_args()
    if args.mode == "prepare":
        if args.baseline is None:
            parser.error("prepare requires --baseline")
        prepare(args)
    else:
        if args.task is None or not 0 <= args.task < 72:
            parser.error("run requires --task in 0..71")
        run(args)


if __name__ == "__main__":
    main()
