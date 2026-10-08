#!/usr/bin/env python3
"""Compare extracted LINE hits and reproduce archived per-read TE selections."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import tempfile

import pysam
from RelocaTE3.librelocate import RelocaTE
from RelocaTE3.aligners import _iter_psl_to_sam


def psl_record(line):
    f = line.split("\t")
    if len(f) != 21 or not f[0].isdigit():
        raise ValueError("Expected an extracted PSL alignment")
    admitted = not (int(f[4]) > 1 or int(f[5]) > 3 or int(f[6]) > 1
                    or int(f[7]) > 3 or int(f[17]) >= 3)
    qs, qe, ql, ts, te, tl = map(int, (f[11], f[12], f[10], f[15], f[16], f[14]))
    return dict(qname=f[9], tName=f[13], match=int(f[0]), mismatch=int(f[1]),
                boundary=int(qs <= 2) + int(qe - 1 >= ql - 3)
                + int(ts <= 2) + int(te - 1 >= tl - 3),
                start=qs, end=qe - 1, tStart=ts, tEnd=te - 1,
                len=ql, tLen=tl, strand=f[8], admitted=admitted)


def select_psl(evidence):
    hits = {}
    for rows in evidence["R2_psl"].values():
        for row in rows:
            rec = psl_record(row["line"])
            rec["line_number"] = row["line_number"]
            hits.setdefault(rec["qname"], []).append(rec)
    return hits


def mappings(evidence, caller):
    result = {}
    for source, rows in evidence[caller + "_mapping"].items():
        for row in rows:
            name, family = row["line"].split("\t")[:2]
            if name in result and result[name] != family:
                raise ValueError(f"Inconsistent mapping for {name} in {source}")
            result[name] = family
    return result


def replay_bams(evidence, directory):
    """Call the production parser on 77 original SAM records, with local indexes."""
    merged = {}
    for side, value in evidence["R3_alignments"].items():
        path = directory / (side + ".bam")
        header = pysam.AlignmentHeader.from_dict(value["header"])
        with pysam.AlignmentFile(str(path), "wb", header=header) as out:
            for line in value["sam"]:
                out.write(pysam.AlignedSegment.fromstring(line, header))
        pysam.index(str(path))
        coord = RelocaTE()._parse_te_bam(path)
        for name, rec in coord.items():
            if name not in merged or RelocaTE._is_better(rec, merged[name]):
                merged[name] = rec
    return merged


def analyze(evidence):
    hits = select_psl(evidence)
    r2map, r3map = mappings(evidence, "R2"), mappings(evidence, "R3")
    with tempfile.TemporaryDirectory(prefix="line-family-") as tmp:
        r3replay = replay_bams(evidence, Path(tmp))
    rows = []
    for name in evidence["read_names"]:
        valid = [r for r in hits[name] if r["admitted"]]
        if not valid:
            raise ValueError(f"No eligible R2 hit for {name}")
        # Python max keeps the first equal key, exactly R2's strict replacement.
        legacy = max(valid, key=lambda r: (r["boundary"], r["match"]))
        deterministic = min(valid, key=RelocaTE._match_rank)
        assert legacy["tName"] == r2map[name], name
        assert r3replay[name]["tName"] == r3map[name], name
        top = (legacy["boundary"], legacy["match"])
        families = sorted({r["tName"] for r in valid if (r["boundary"], r["match"]) == top})
        row = dict(read_name=name, R2_family=r2map[name], R3_family=r3map[name],
                   legacy_replay_family=legacy["tName"],
                   deterministic_on_R2_hits=deterministic["tName"],
                   production_R3_replay_family=r3replay[name]["tName"],
                   top_boundary=top[0], top_match=top[1],
                   top_families=json.dumps(families),
                   eligible_R2_hits=len(valid), rejected_R2_hits=len(hits[name])-len(valid),
                   changed_family=int(r2map[name] != r3map[name]),
                   R2_policy_change_reproduces_R3=int(deterministic["tName"] == r3map[name]),
                   top_psl_records=json.dumps([r for r in valid if (r["boundary"], r["match"]) == top], sort_keys=True))
        rows.append(row)
    def signature(line):
        fields = line.split("\t")
        nm = next(s for s in fields[11:] if s.startswith("NM:i:"))
        return tuple(fields[:6]) + (nm,)

    original_psl = [r["line"] for rs in evidence["R2_psl"].values() for r in rs]
    r2signatures = Counter(signature(line) for line in _iter_psl_to_sam(original_psl))
    r3signatures = Counter(signature(line) for value in evidence["R3_alignments"].values() for line in value["sam"])
    summary = {"reads": len(rows), "changed_family": sum(r["changed_family"] for r in rows),
               "R2_policy_change_reproduces_R3": sum(r["R2_policy_change_reproduces_R3"] for r in rows),
               "R2_votes": dict(Counter(r["R2_family"] for r in rows)),
               "R3_votes": dict(Counter(r["R3_family"] for r in rows)),
               "rejected_R2_hits": sum(r["rejected_R2_hits"] for r in rows),
               "eligible_alignment_signatures_equal": r2signatures == r3signatures,
               "R2_only_eligible_alignments": sum((r2signatures-r3signatures).values()),
               "R3_only_eligible_alignments": sum((r3signatures-r2signatures).values()),
               "ambiguous_best_family_reads": sum(len(json.loads(r["top_families"])) > 1 for r in rows),
               "single_best_family_votes": dict(Counter(json.loads(r["top_families"])[0]
                    for r in rows if len(json.loads(r["top_families"])) == 1))}
    return rows, summary


def compact_fixture(evidence, source_hash):
    """Keep original read evidence while omitting unrelated header targets."""
    fixture = dict(sample=evidence["sample"], event_id=evidence["event_id"],
                   read_names=evidence["read_names"], R2_mapping=evidence["R2_mapping"],
                   R3_mapping=evidence["R3_mapping"], R2_psl=evidence["R2_psl"],
                   R3_alignments={}, source_evidence_sha256=source_hash,
                   note="Observed family disagreement; fixture does not assert biological correctness of R3 labels.")
    for side, value in evidence["R3_alignments"].items():
        targets = {line.split("\t")[2] for line in value["sam"]}
        fixture["R3_alignments"][side] = {
            "header": {"HD": value["header"].get("HD", {}),
                       "SQ": [r for r in value["header"]["SQ"] if r["SN"] in targets]},
            "sam": value["sam"]}
    return fixture


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--emit-fixture", action="store_true")
    args = parser.parse_args()
    data = args.evidence.read_bytes()
    evidence = json.loads(data)
    digest = hashlib.sha256(data).hexdigest()
    if args.emit_fixture:
        print(json.dumps(compact_fixture(evidence, digest), indent=2))
        return
    if args.outdir is None:
        parser.error("--outdir is required unless --emit-fixture is used")
    if args.outdir.exists():
        raise FileExistsError(args.outdir)
    rows, summary = analyze(evidence)
    args.outdir.mkdir(parents=True, exist_ok=False)
    with (args.outdir / "per_read.tsv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    (args.outdir / "provenance.json").write_text(json.dumps({
        "evidence_sha256": digest,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "parser_source_sha256": hashlib.sha256(Path(__import__("inspect").getfile(RelocaTE)).read_bytes()).hexdigest(),
        "converter_source_sha256": hashlib.sha256(Path(__import__("inspect").getfile(_iter_psl_to_sam)).read_bytes()).hexdigest(),
        "pysam_version": pysam.__version__,
        "note": "Tiny extracted BAM replay; does not measure full benchmark performance."}, indent=2))
    (args.outdir / ".complete").touch()
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
