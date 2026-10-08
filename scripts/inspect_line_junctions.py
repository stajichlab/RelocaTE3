#!/usr/bin/env python3
"""Bounded indexed inspection of one recurrent LINE family-disagreement locus."""
import argparse
import hashlib
import json
from pathlib import Path
import pysam
from RelocaTE3.insertions import _junction_info


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.outdir.exists():
        raise FileExistsError(args.outdir)
    base = args.benchmark / "runs/ricetelib_divergence"
    sample = "div005_rep03_cov30x"
    r2 = base / "relocate2" / sample / "raw/repeat"
    readlist = r2 / "results/Chr1.repeat.reads.list"
    data = readlist.read_bytes()
    line = [s.split("\t") for s in data.decode().splitlines()
            if "\tChr1:29402710..29402725\tJunction_reads\t" in s]
    assert len(line) == 1
    names = set(line[0][3].split(","))
    paths = {"R2": r2 / "bwa_aln/MSU_r7.repeat.bwa.sorted.bam",
             "R3": base / "relocate3-blat-bwaaln" / sample / "raw" / f"{sample}.repeat.bwaaln.sorted.bam"}
    result = {"sample": sample, "event_id": "TE000172", "R2_reported_junction_names": sorted(names),
              "query_0based_half_open": [29402500, 29402900], "bam": {}}
    for caller, path in paths.items():
        with pysam.AlignmentFile(str(path), "rb") as bam:
            records = list(bam.fetch("Chr1", 29402500, 29402900))
        if len(records) > 1000:
            raise RuntimeError("Bounded query exceeded budget")
        selected = []
        for rec in records:
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            info = _junction_info(rec.query_name, "-" if rec.is_reverse else "+",
                                  rec.reference_start + 1, rec.reference_end)
            if info and ((info[0] == "left" and info[1] == 29402725)
                         or (info[0] == "right" and info[1] == 29402710)):
                selected.append(rec.query_name)
        result["bam"][caller] = {"path": str(path), "size_bytes": path.stat().st_size,
            "index_sha256": hashlib.sha256(Path(str(path) + ".bai").read_bytes()).hexdigest(),
            "breakpoint_junction_names": sorted(selected), "local_sam": [r.to_string() for r in records]}
    a = set(result["bam"]["R2"]["breakpoint_junction_names"])
    b = set(result["bam"]["R3"]["breakpoint_junction_names"])
    result.update(R2_only_names=sorted(a-b), R3_only_names=sorted(b-a), shared_names=len(a & b),
                  R2_breakpoints_match_report=a == names, pysam_version=pysam.__version__,
                  readlist_sha256=hashlib.sha256(data).hexdigest(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.outdir.mkdir(parents=True, exist_ok=False)
    (args.outdir / "junction_evidence.json").write_text(json.dumps(result, indent=2))
    (args.outdir / ".complete").touch()
    print(json.dumps({k: result[k] for k in ("R2_only_names", "R3_only_names", "shared_names", "R2_breakpoints_match_report")}))


if __name__ == "__main__":
    main()
