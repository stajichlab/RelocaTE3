#!/usr/bin/env python3
"""Bounded indexed extraction of the recurrent Os3912 false-positive locus."""

import argparse
import hashlib
import json
from pathlib import Path

import pysam
from RelocaTE3.insertions import (
    _junction_info,
    _junction_fullread_key,
    _fullread_record_key,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", required=True, type=Path)
    parser.add_argument("--benchmark", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    manifest = json.loads((args.replay / "manifest.json").read_text())
    task = next(
        t
        for t in manifest["tasks"]
        if t["dataset"] == "ricetelib" and t["sample"] == "cov30x_rep1"
    )
    paths = {
        "R3_junctions": Path(task["inputs"]["genome_bam"]),
        "R3_fullreads": Path(task["inputs"]["fullreads"]),
    }
    base = args.benchmark / "runs/ricetelib/relocate2/cov30x_rep1/raw/repeat/bwa_aln"
    paths.update(
        R2_junctions=base / "MSU_r7.repeat.bwa.sorted.bam",
        R2_fullreads=base / "MSU_r7.repeat.fullreads.bwa.sorted.bam",
    )
    result = dict(
        sample=task["sample"],
        chrom="Chr4",
        start=19264325,
        end=19264335,
        query_0based_half_open=[19263825, 19264825],
        alignments={},
    )
    names = set()
    for label, path in paths.items():
        with pysam.AlignmentFile(str(path), "rb") as bam:
            records = list(bam.fetch("Chr4", 19263825, 19264825))
            if len(records) > 2000:
                raise RuntimeError("Bounded query budget exceeded; use SLURM")
            selected = []
            for r in records:
                if "junctions" in label:
                    info = _junction_info(
                        r.query_name,
                        "-" if r.is_reverse else "+",
                        r.reference_start + 1,
                        r.reference_end,
                    )
                    if not info or (info[0], info[1]) not in (
                        ("left", 19264335),
                        ("right", 19264325),
                    ):
                        continue
                elif _fullread_record_key(r) not in names:
                    continue
                selected.append(r)
            if label == "R3_junctions":
                names = {_junction_fullread_key(r.query_name) for r in selected}
            result["alignments"][label] = dict(
                path=str(path),
                size=path.stat().st_size,
                index_sha256=hashlib.sha256(
                    Path(str(path) + ".bai").read_bytes()
                ).hexdigest(),
                header=bam.header.to_dict(),
                sam=[r.to_string() for r in selected],
            )
    result["script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2)
    print({label: len(d["sam"]) for label, d in result["alignments"].items()})


if __name__ == "__main__":
    main()
