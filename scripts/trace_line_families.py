#!/usr/bin/env python3
"""Extract archived family mappings and all TE hits for eleven LINE junction reads.

Large sequential scans must run on a compute node. No alignments are rerun.
Outputs are self-contained evidence, not an inferred biological correction.
"""
import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pysam


def select_lines(path, names, column, hashes):
    digest, selected = hashlib.sha256(), []
    with path.open("rb") as handle:
        for number, line in enumerate(handle, 1):
            digest.update(line)
            fields = line.decode().rstrip("\r\n").split("\t")
            if len(fields) > column and fields[column] in names:
                selected.append({"line_number": number, "line": line.decode().rstrip("\r\n")})
    hashes[str(path)] = digest.hexdigest()
    return selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", required=True, type=Path)
    parser.add_argument("--junctions", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    if args.outdir.exists():
        raise FileExistsError(f"Refusing existing output: {args.outdir}")
    junctions = json.loads(args.junctions.read_text())
    names = {re.sub(r":(start|end):[53]$", "", n)
             for n in junctions["R2_reported_junction_names"]}
    if len(names) != 11 or not junctions["R2_breakpoints_match_report"]:
        raise ValueError("Unexpected junction fixture")
    args.outdir.mkdir(parents=True, exist_ok=False)
    sample = junctions["sample"]
    root = args.benchmark / "runs/ricetelib_divergence"
    r2 = root / "relocate2" / sample / "raw/repeat"
    r3 = root / "relocate3-blat-bwaaln" / sample / "raw"
    hashes = {str(args.junctions): hashlib.sha256(args.junctions.read_bytes()).hexdigest()}
    result = {"sample": sample, "event_id": "TE000172", "read_names": sorted(names),
              "R2_mapping": {}, "R2_psl": {}, "R3_mapping": {}, "R3_alignments": {}}
    mapping = r2 / "te_containing_fq/Chr1.read_repeat_name.split.txt"
    result["R2_mapping"][str(mapping)] = select_lines(mapping, names, 0, hashes)
    # Locate original chunk files before scanning only the corresponding PSLs.
    for path in sorted((r2 / "te_containing_fq").glob("*.read_repeat_name.txt")):
        records = select_lines(path, names, 0, hashes)
        if not records:
            continue
        result["R2_mapping"][str(path)] = records
        psl = r2 / "blat_output" / path.name.replace(".read_repeat_name.txt", ".blatout")
        result["R2_psl"][str(psl)] = select_lines(psl, names, 9, hashes)
        print(f"Selected {len(records)} mapping rows from {path.name}", flush=True)
    mapping = r3 / "te_containing" / f"{sample}.read_repeat_name.txt"
    result["R3_mapping"][str(mapping)] = select_lines(mapping, names, 0, hashes)
    bam_provenance = {}
    for side in ("left", "right"):
        path = r3 / f"{sample}.{side}.bam"
        before = path.stat()
        records = []
        with pysam.AlignmentFile(str(path), "rb") as bam:
            header = bam.header.to_dict()
            for rec in bam.fetch(until_eof=True):
                if rec.query_name in names:
                    records.append(rec.to_string())
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise RuntimeError(f"Input changed during extraction: {path}")
        result["R3_alignments"][side] = {"header": header, "sam": records}
        bam_provenance[str(path)] = {"size_bytes": after.st_size, "mtime_ns": after.st_mtime_ns}
        print(f"Selected {len(records)} TE alignments from {side} BAM", flush=True)
    for caller in ("R2", "R3"):
        observed = {row["line"].split("\t")[0]
                    for rows in result[caller + "_mapping"].values() for row in rows}
        if observed != names:
            raise RuntimeError(f"Missing {caller} mappings: {names - observed}")
    psl_names = {row["line"].split("\t")[9]
                 for rows in result["R2_psl"].values() for row in rows}
    sam_names = {line.split("\t")[0]
                 for value in result["R3_alignments"].values() for line in value["sam"]}
    if psl_names != names or sam_names != names:
        raise RuntimeError(f"Missing hits: R2={names - psl_names}, R3={names - sam_names}")
    result["provenance"] = {"recorded_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256": hashes, "bam_metadata": bam_provenance,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "pysam_version": pysam.__version__,
        "note": "Original PSLs retained for R2; R3 temporary PSLs absent, archived TE BAM used. BAMs identified by size/mtime, not full hash."}
    (args.outdir / "family_evidence.json").write_text(json.dumps(result, indent=2))
    (args.outdir / ".complete").touch()
    print("Family evidence extraction complete", flush=True)


if __name__ == "__main__":
    main()
