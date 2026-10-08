#!/usr/bin/env python3
"""Run full tests on SLURM, recording source hashes and detecting edits mid-run."""

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def fingerprint():
    paths = [
        p
        for directory in ("src", "tests", "scripts")
        for p in Path(directory).rglob("*")
        if p.is_file()
        and p.suffix in (".py", ".json", ".slurm")
        and "results" not in p.parts
        and "__pycache__" not in p.parts
    ]
    paths += [
        Path(p)
        for p in ("pyproject.toml", "pixi.toml", "pixi.lock")
        if Path(p).is_file()
    ]
    return {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def main():
    job = os.environ.get("SLURM_JOB_ID")
    if not job:
        raise RuntimeError("Full tests require SLURM")
    out = Path("results/checkpoint-tests") / job
    out.mkdir(parents=True, exist_ok=False)
    before = fingerprint()
    tools = {
        name: shutil.which(name)
        for name in (
            "python",
            "relocaTE3",
            "bwa",
            "bwa-mem2",
            "samtools",
            "blat",
            "minimap2",
            "bedtools",
            "bowtie2",
            "seqtk",
            "bcftools",
        )
    }
    if not all(tools.values()):
        raise RuntimeError(f"Missing test executable: {tools}")
    with (out / "provenance.json").open("x") as handle:
        json.dump(
            dict(
                job=job,
                git_head=subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], text=True
                ).strip(),
                source_sha256=before,
                executables=tools,
                python=sys.version,
            ),
            handle,
            indent=2,
        )
    rc = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-ra",
            "-q",
            "tests",
            "--junitxml",
            str(out / "full-tests.xml"),
        ]
    ).returncode
    if fingerprint() != before:
        raise RuntimeError(
            "Source/test inputs changed during full tests; result not a checkpoint gate"
        )
    if rc:
        raise SystemExit(rc)
    (out / ".complete").touch()
    print(f"Full tests passed for recorded source: {out}")


if __name__ == "__main__":
    main()
