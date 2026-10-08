#!/usr/bin/env python3.12
"""Compare frozen baseline/candidate trimming from saved BAMs on a compute node."""

import argparse
import json
import os
from pathlib import Path
import shutil

from profile_benchmark_memory import digest, monitor, write_json


def artifact_hashes(raw):
    """Hash all trim artifacts, including empty files and retained weak TE hits."""
    result = {}
    for folder in ("flanking", "te_containing", "te_portions"):
        directory = raw / folder
        if not directory.is_dir():
            raise ValueError(f"Missing trim output directory: {directory}")
        for path in sorted(directory.rglob("*")):
            if path.is_file():
                result[str(path.relative_to(raw))] = digest(path)
    if not result:
        raise ValueError("No trim artifacts found")
    return result


def verify_snapshot(root, manifest):
    for name, expected in manifest["snapshot_sha256"].items():
        if digest(root / name) != expected:
            raise ValueError(f"Snapshot changed: {root / name}")


def validate(baseline, candidate):
    if not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("Run substantial trimming validation through SLURM")
    baseline, candidate = baseline.resolve(), candidate.resolve()
    if not (baseline / ".profile_complete").is_file():
        raise ValueError("Baseline must be a completed, call-identical profiling run")
    old = json.loads((baseline / "manifest.json").read_text())
    new = json.loads((candidate / "manifest.json").read_text())
    for key in ("sample", "caller", "dataset", "threads", "match_window"):
        if old[key] != new[key]:
            raise ValueError(f"Baseline/candidate disagree on {key}")
    for key in ("R1", "R2", "TE_LIBRARY", "REFERENCE", "REPEATMASKER"):
        if old["env"][key] != new["env"][key]:
            raise ValueError(f"Baseline/candidate disagree on input {key}")
    verify_snapshot(baseline, old)
    verify_snapshot(candidate, new)
    output = candidate / "trim-validation"
    output.mkdir(exist_ok=False)
    raw = baseline / "replay/raw"
    sample = old["sample"]
    bams = [raw / f"{sample}.{side}.bam" for side in ("left", "right")]
    fastqs = [Path(old["env"][key]) for key in ("R1", "R2")]
    inputs = bams + [Path(str(bam) + ".bai") for bam in bams] + fastqs
    # Expensive file hashing is deliberately performed only on the compute node.
    hashes = {str(path): digest(path) for path in inputs}
    original_hashes = json.loads((baseline / "input_sha256.json").read_text())
    for fastq in fastqs:
        if hashes[str(fastq)] != original_hashes[str(fastq)]:
            raise ValueError(f"FASTQ differs from completed baseline: {fastq}")
    write_json(output / "inputs.json", hashes)
    expected = artifact_hashes(raw)
    comparisons = {}
    python = shutil.which("python")
    if not python:
        raise ValueError("Activate the benchmark Python environment first")
    for label, source in (("baseline", baseline), ("candidate", candidate)):
        dest = output / label
        dest.mkdir()
        env = os.environ.copy()
        env["PYTHONPATH"] = str(source / "snapshot/source/src")
        env["PYTHONDONTWRITEBYTECODE"] = "1"
        command = [
            "/usr/bin/time",
            "-v",
            "-o",
            str(dest / "time-v.txt"),
            python,
            "-m",
            "RelocaTE3",
            "trim",
            "--bam",
            *(str(path) for path in bams),
            "--fastq",
            *(str(path) for path in fastqs),
            "--name",
            sample,
            "--outdir",
            str(dest / "raw"),
            "--min-match",
            "10",
            "--min-trimmed",
            "10",
            "--mismatch",
            "2",
        ]
        write_json(
            dest / "command.json", {"argv": command, "source": env["PYTHONPATH"]}
        )
        print(f"Starting {label} trim-only profile", flush=True)
        code = monitor(command, Path.cwd(), env, dest / "profile")
        if code:
            raise RuntimeError(f"{label} trimming failed with exit {code}")
        actual = artifact_hashes(dest / "raw")
        comparisons[label] = {
            "identical": actual == expected,
            "files": actual,
            "different_or_missing": sorted(
                key
                for key in expected.keys() | actual.keys()
                if expected.get(key) != actual.get(key)
            ),
        }
        write_json(dest / "comparison.json", comparisons[label])
        if not comparisons[label]["identical"]:
            raise RuntimeError(
                f"{label} trim artifacts differ; inspect {dest / 'comparison.json'}"
            )
    write_json(output / "comparison.json", comparisons)
    (output / ".complete").touch(exist_ok=False)
    print(
        "Both isolated trims match all original trim artifacts byte-for-byte.",
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    args = parser.parse_args()
    try:
        validate(args.baseline, args.candidate)
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(1, f"ERROR: {error}\n")


if __name__ == "__main__":
    main()
