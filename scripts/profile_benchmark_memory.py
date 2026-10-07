#!/usr/bin/env python3.12
"""Preserve and replay one benchmark sample with Linux process-memory sampling.

Prepare is lightweight; run is compute-node-only. No calling code is modified.
RSS sums include shared pages more than once; they are not unique memory usage.
Short-lived processes can escape sampling; time -v supplies a kernel high-water
mark for the entire adapter. Python/native allocation tracing is a later step.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys
import time


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_json(path, value):
    with Path(path).open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def capture(command, cwd=None):
    return subprocess.check_output(command, cwd=cwd, text=True).strip()


def table_rows(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return Counter(tuple(sorted(row.items())) for row in reader)


def copy_file(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def prepare(args):
    repo, benchmark, output = (
        Path.cwd(),
        args.benchmark.resolve(),
        args.output.resolve(),
    )
    if not (repo / "src/RelocaTE3/cli.py").is_file():
        raise ValueError("Run from the RelocaTE3 repository root")
    bridge = [sys.executable, "pipeline/config_env.py", "--config", args.config]
    tasks = capture(bridge + ["--dataset", "ricetelib", "tasks"], benchmark)
    selected = [
        line.split("\t")
        for line in tasks.splitlines()
        if line.split("\t")[1:3] == ["relocate3-blat-bwaaln", args.sample]
    ]
    if len(selected) != 1 or selected[0][3] != "30":
        raise ValueError("Expected exactly one 30x riceTElib BLAT/bwa-aln task")
    dataset, caller, sample, coverage, replicate, r1, r2 = selected[0]
    # Read shell assignments through shlex, never execute configuration text.
    settings = {}
    for mode in (["globals"], ["dataset-env", dataset], ["caller-env", caller]):
        for line in capture(bridge + mode, benchmark).splitlines():
            key, value = line.split("=", 1)
            settings[key] = shlex.split(value)[0]
    if Path(settings["RT3_REPO"]).resolve() != repo.resolve():
        raise ValueError(
            "Benchmark configuration points at a different source checkout"
        )
    if (
        settings["RT3_TE_ALIGNER"] != "blat"
        or settings["RT3_GENOME_ALIGNER"] != "bwaaln"
    ):
        raise ValueError("Unexpected aligner settings")
    for key in ("REFERENCE", "TE_LIBRARY", "REPEATMASKER"):
        settings[key] = str((benchmark / settings[key]).resolve())
    old_run = benchmark / settings["DATASET_WORK_ROOT"] / caller / sample
    old_report = (
        benchmark / settings["DATASET_REPORT_ROOT"] / "per_sample" / caller / sample
    )
    truth = benchmark / settings["DATASET_TRUTH_ROOT"] / "per_sample" / f"{sample}.tsv"
    if not (old_run / ".run_complete").is_file():
        raise ValueError("Baseline caller run has no completion sentinel")
    # Refuse all existing destinations, including incomplete previous attempts.
    output.mkdir(parents=True, exist_ok=False)
    baseline, snapshot = output / "baseline", output / "snapshot"
    for filename in ("calls.normalized.tsv", ".run_complete"):
        copy_file(old_run / filename, baseline / filename)
    for path in (old_run / "raw/results").iterdir():
        if path.is_file():
            copy_file(path, baseline / "raw_results" / path.name)
    for path in old_report.glob("*.tsv"):
        copy_file(path, baseline / "score" / path.name)
    copy_file(truth, baseline / "truth.tsv")
    copy_file(benchmark / args.config, snapshot / "benchmark.toml")
    copy_file(
        benchmark / settings["DATASET_REPORT_ROOT"] / "resources.tsv",
        baseline / "panel_resources.tsv",
    )
    for path in (repo / "src").rglob("*.py"):
        copy_file(path, snapshot / "source" / path.relative_to(repo))
    for name in ("pixi.toml", "pixi.lock", "pyproject.toml", "requirements.txt"):
        copy_file(repo / name, snapshot / "source" / name)
    for name in (
        "profile_benchmark_memory.py",
        "memory_profile_adapter.sh",
        "profile_benchmark_memory.slurm",
        "validate_trim_memory.py",
        "validate_memory_changes.slurm",
    ):
        copy_file(repo / "scripts" / name, snapshot / name)
    # Freeze adapter, normalizer, scorer, and their pure-Python dependencies.
    for folder in ("lib", "scoring", "callers/relocate3"):
        for path in (benchmark / folder).glob("*"):
            if path.is_file() and path.suffix in (".py", ".sh", ".toml", ".lock"):
                copy_file(path, snapshot / "benchmark" / folder / path.name)
    adapter = snapshot / "benchmark/callers/relocate3"
    # Reuse benchmark tool installations without copying or updating them.
    for name in ("blat-env", ".pixi"):
        source = benchmark / "callers/relocate3" / name
        if source.exists():
            (adapter / name).symlink_to(source, target_is_directory=True)
    settings.update(
        SAMPLE=sample,
        R1=r1,
        R2=r2,
        TARGET="ALL",
        OUTDIR=str(output / "replay"),
        ADAPTER_DIR=str(adapter),
        PYTHONPATH=str(snapshot / "source/src"),
        PYTHONDONTWRITEBYTECODE="1",
    )
    allowed = (
        "SAMPLE",
        "R1",
        "R2",
        "REFERENCE",
        "TE_LIBRARY",
        "REPEATMASKER",
        "TE_NAME",
        "TSD_PATTERN",
        "TARGET",
        "OUTDIR",
        "ADAPTER_DIR",
        "RT3_REPO",
        "RT3_TE_ALIGNER",
        "RT3_GENOME_ALIGNER",
        "RT3_REQUIRE_BOTH_JUNCTIONS",
        "PYTHONPATH",
        "PYTHONDONTWRITEBYTECODE",
    )
    inputs = [
        Path(settings[key])
        for key in ("R1", "R2", "REFERENCE", "TE_LIBRARY", "REPEATMASKER")
    ]
    # Indexes must already exist; do not permit writes beside external references.
    inputs += [
        Path(settings["REFERENCE"] + suffix)
        for suffix in (".fai", ".mmi", ".bwt", ".ann", ".amb", ".pac", ".sa")
    ]
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_commit": capture(["git", "rev-parse", "HEAD"], repo),
        "source_status": capture(["git", "status", "--short"], repo),
        "benchmark_commit": capture(["git", "rev-parse", "HEAD"], benchmark),
        "benchmark_status": capture(["git", "status", "--short"], benchmark),
        "benchmark_root": str(benchmark),
        "dataset": dataset,
        "caller": caller,
        "sample": sample,
        "coverage": coverage,
        "replicate": replicate,
        "threads": int(settings["THREADS"]),
        "match_window": settings["MATCH_WINDOW"],
        "env": {key: settings[key] for key in allowed},
        "inputs": {
            str(path): {
                "size": path.stat().st_size,
                "mtime_ns": path.stat().st_mtime_ns,
            }
            for path in inputs
        },
    }
    for root in (repo, benchmark):
        name = "source" if root == repo else "benchmark"
        (snapshot / f"{name}.diff").write_text(
            capture(["git", "diff", "HEAD"], root) + "\n"
        )
    metadata["snapshot_sha256"] = {
        str(path.relative_to(output)): digest(path)
        for root in (baseline, snapshot)
        for path in root.rglob("*")
        if path.is_file() and not path.is_symlink()
    }
    write_json(output / "manifest.json", metadata)
    print(f"Prepared immutable call/source baseline: {output}")


def read_process(directory, session):
    """Read one /proc record; process disappearance is normal during sampling."""
    stat = (directory / "stat").read_text().rsplit(")", 1)[1].split()
    if int(stat[3]) != session:
        return None
    status = {}
    for line in (directory / "status").read_text().splitlines():
        key, _, value = line.partition(":")
        if key in ("VmRSS", "VmHWM"):
            status[key] = int(value.split()[0])
    command = (
        (directory / "cmdline")
        .read_bytes()
        .replace(b"\0", b" ")
        .decode(errors="replace")
        .strip()
    )
    return {
        "pid": int(directory.name),
        "ppid": int(stat[1]),
        "start_ticks": int(stat[19]),
        "rss_kib": status.get("VmRSS", 0),
        "hwm_kib": status.get("VmHWM", 0),
        "command": command,
    }


def monitor(command, cwd, env, output, interval=2.0):
    """Run a new process session and sample that session, including grandchildren."""
    output.mkdir(exist_ok=False)
    peaks, maximum_sum = {}, 0
    started = time.monotonic()
    fields = [
        "elapsed_seconds",
        "pid",
        "ppid",
        "start_ticks",
        "rss_kib",
        "hwm_kib",
        "command",
    ]
    with (
        (output / "adapter.log").open("x") as log,
        (output / "process_memory.tsv").open("x") as handle,
    ):
        process = subprocess.Popen(
            command,
            cwd=cwd,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )

        def stop(signum, frame):
            raise InterruptedError(f"Profiler received signal {signum}")

        original = {
            sig: signal.signal(sig, stop) for sig in (signal.SIGTERM, signal.SIGINT)
        }
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        try:
            while True:
                rss_sum = 0
                elapsed = round(time.monotonic() - started, 3)
                for directory in Path("/proc").iterdir():
                    if not directory.name.isdigit():
                        continue
                    try:
                        row = read_process(directory, process.pid)
                    except (OSError, ValueError, IndexError):
                        continue
                    if row is None:
                        continue
                    writer.writerow({"elapsed_seconds": elapsed, **row})
                    rss_sum += row["rss_kib"]
                    key = (row["pid"], row["start_ticks"], row["command"])
                    peak = peaks.setdefault(key, {**row, "first_seen_seconds": elapsed})
                    peak["rss_kib"] = max(peak["rss_kib"], row["rss_kib"])
                    peak["hwm_kib"] = max(peak["hwm_kib"], row["hwm_kib"])
                    peak["last_seen_seconds"] = elapsed
                maximum_sum = max(maximum_sum, rss_sum)
                handle.flush()
                if process.poll() is not None:
                    break
                time.sleep(interval)
        finally:
            # No orphan aligners on cancellation or an error in the sampler.
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
            for sig, handler in original.items():
                signal.signal(sig, handler)
            write_json(
                output / "summary.json",
                {
                    "exit_code": process.returncode,
                    "wall_seconds": time.monotonic() - started,
                    "sampling_interval_seconds": interval,
                    "max_sampled_rss_sum_kib": maximum_sum,
                    "note": "RSS sum double-counts shared pages; sampling can miss short-lived processes. "
                    "VmHWM survives exec and is a process-lifetime, not per-command, maximum.",
                    "processes": sorted(
                        peaks.values(), key=lambda row: row["rss_kib"], reverse=True
                    ),
                },
            )
    return process.returncode


def run(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise ValueError(
            "Substantial replay must run through SLURM, not on a login node"
        )
    output = args.output.resolve()
    metadata = json.loads((output / "manifest.json").read_text())
    if int(os.environ.get("SLURM_CPUS_PER_TASK", "1")) != metadata["threads"]:
        raise ValueError("CPU allocation must match the saved benchmark thread count")
    # Exclusive marker also prevents concurrent duplicate replays.
    with (output / "started.json").open("x") as handle:
        json.dump(
            {
                "job_id": os.environ["SLURM_JOB_ID"],
                "started_utc": datetime.now(timezone.utc).isoformat(),
            },
            handle,
        )
    for path, expected in metadata["snapshot_sha256"].items():
        if digest(output / path) != expected:
            raise ValueError(f"Baseline/snapshot changed: {path}")
    hashes = {}
    for name, expected in metadata["inputs"].items():
        path = Path(name)
        if (
            path.stat().st_size != expected["size"]
            or path.stat().st_mtime_ns != expected["mtime_ns"]
        ):
            raise ValueError(f"Input changed since preparation: {name}")
        hashes[name] = digest(path)
    write_json(output / "input_sha256.json", hashes)
    env = os.environ.copy()
    # Remove inherited caller knobs; the archived adapter supplies its defaults.
    for key in list(env):
        if key.startswith("RT3_") or key in ("TARGET", "TSD_PATTERN", "TE_NAME"):
            del env[key]
    env.update(metadata["env"])
    env.update(THREADS=str(metadata["threads"]), PROFILE_OUTPUT=str(output))
    temp = output / "tmp"
    temp.mkdir()
    env["TMPDIR"] = str(temp)
    snapshot = output / "snapshot"
    command = [
        "/usr/bin/time",
        "-v",
        "-o",
        str(output / "adapter.time-v.txt"),
        "bash",
        str(snapshot / "memory_profile_adapter.sh"),
    ]
    print(
        f"Profiling {metadata['sample']}; log: {output / 'profile/adapter.log'}",
        flush=True,
    )
    code = monitor(command, metadata["benchmark_root"], env, output / "profile")
    if code:
        raise RuntimeError(
            f"Adapter failed with exit {code}; inspect profile/adapter.log"
        )
    code_root = snapshot / "benchmark"
    subprocess.run(
        [
            sys.executable,
            str(code_root / "callers/relocate3/normalize.py"),
            "--outdir",
            str(output / "replay"),
            "--sample",
            metadata["sample"],
            "--te-name",
            env["TE_NAME"],
            "--target",
            "ALL",
        ],
        env=env,
        check=True,
    )
    subprocess.run(
        [
            sys.executable,
            str(code_root / "scoring/score_calls.py"),
            "--truth",
            str(output / "baseline/truth.tsv"),
            "--calls",
            str(output / "replay/calls.normalized.tsv"),
            "--sample",
            metadata["sample"],
            "--caller",
            metadata["caller"],
            "--window",
            metadata["match_window"],
            "--outdir",
            str(output / "score"),
        ],
        env=env,
        check=True,
    )
    checks = {
        "normalized_calls": table_rows(output / "baseline/calls.normalized.tsv")
        == table_rows(output / "replay/calls.normalized.tsv")
    }
    for filename in (
        "matches.tsv",
        "precision.tsv",
        "correctness.tsv",
        "false_positive_calls.tsv",
    ):
        checks[filename] = table_rows(
            output / "baseline/score" / filename
        ) == table_rows(output / "score" / filename)
    write_json(output / "comparison.json", checks)
    if not all(checks.values()):
        raise RuntimeError(
            "Replay differs from saved baseline; investigate comparison.json before optimizing"
        )
    (output / ".profile_complete").touch(exist_ok=False)
    print("Profiling complete; normalized calls and scored tables match baseline.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)
    prep = sub.add_parser(
        "prepare",
        help="Snapshot existing calls/code; no alignment or large-file hashing",
    )
    prep.add_argument(
        "--benchmark",
        type=Path,
        default=Path("../../relocate_benchmark/relocate-benchmark"),
    )
    prep.add_argument("--config", default="config/benchmark.full-aligners.toml")
    prep.add_argument("--sample", default="cov30x_rep1")
    prep.add_argument("--output", required=True, type=Path)
    replay = sub.add_parser("run", help="Replay and profile on a SLURM compute node")
    replay.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        (prepare if args.mode == "prepare" else run)(args)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"ERROR: {error}\n")


if __name__ == "__main__":
    main()
