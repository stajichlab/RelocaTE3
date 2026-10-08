#!/usr/bin/env python3
"""Extract two bounded, indexed benchmark clusters into portable test fixtures.

Run from the repository root using its existing Python environment. No caller
reruns or whole-BAM scans. Refuse existing outputs and clusters exceeding the
bounded inspection budget. Generated JSON retains actual observations, not
invented reads or coordinate-shifted copies.
"""

import argparse
import hashlib
import json
import tempfile
from dataclasses import asdict
from pathlib import Path

import pysam

from RelocaTE3 import insertions as ins
from RelocaTE3.reference_te import ReferenceTEAnnotator


def reconstruct(benchmark, sample, position):
    raw = benchmark / "runs/ricetelib/relocate3-blat-bwaaln" / sample / "raw"
    path = raw / f"{sample}.repeat.bwaaln.sorted.bam"
    finder = ins.InsertionFinder(mismatch_allow=2)
    with pysam.AlignmentFile(str(path), "rb") as bam, tempfile.TemporaryDirectory() as tmp:
        local = Path(tmp) / "local.bam"
        for radius in (2000, 4000, 8000, 16000):
            lo, hi = position - radius, position + radius
            records = list(bam.fetch("Chr1", lo, hi))
            if len(records) > 10000:
                raise RuntimeError("Inspection budget exceeded; use SLURM")
            with pysam.AlignmentFile(str(local), "wb", template=bam) as out:
                for record in records:
                    out.write(record)
            clusters = list(ins._stream_clusters(str(local), {}, finder._passes_quality))
            candidates = [c for c in clusters if c.lo <= position <= c.hi]
            if len(candidates) != 1:
                raise RuntimeError(f"Expected one cluster at {sample}:{position}")
            cluster = candidates[0]
            # A full allowance of empty space inside the queried interval on
            # each side guarantees that no external accepted read can chain in.
            if cluster.lo > lo + ins.RANGE_ALLOWANCE + 1 and cluster.hi < hi - ins.RANGE_ALLOWANCE:
                break
        else:
            raise RuntimeError("Cluster boundary unresolved; use SLURM")

        names = {ins._strip_junction_tag(j.read_name) for j in cluster.junctions}
        for name, *_ in cluster.support:
            names.update((name, name + "/1", name + "/2", name + ".f", name + ".r"))
        mapping, digest = {}, hashlib.sha256()
        repeat_path = raw / "te_containing" / f"{sample}.read_repeat_name.txt"
        with repeat_path.open("rb") as handle:
            for line in handle:
                digest.update(line)
                cols = line.decode().rstrip().split("\t")
                if cols[0] in names:
                    mapping[cols[0]] = cols[1:3]
        for j in cluster.junctions:
            j.te_name = ins._te_family(mapping, j.read_name)
        if any(j.te_name == "NA" for j in cluster.junctions):
            raise RuntimeError("Missing junction family annotation")

        annotation = benchmark / "cache/reference_annotations/MSU_r7.full.repeatmasker.out"
        edges = ReferenceTEAnnotator.load_existing_te(annotation, "Chr1")["Chr1"]
        edges = {side: {p: value for p, value in values.items()
                        if cluster.lo - 1000 <= p <= cluster.hi + 1000}
                 for side, values in edges.items()}
        raw_candidates = ins._call_insertions(cluster, None, mapping)
        pooled = ins._consolidate_same_start(raw_candidates, cluster, mapping)
        accepted, filters = [], []
        full_path = raw / "genome_aln" / f"{sample}.fullreads.genome.bam"
        with pysam.AlignmentFile(str(full_path), "rb") as full:
            for candidate in pooled:
                false = ins._fullread_false_junction(full, candidate)
                quality = ins._call_validated_by_high_quality(*ins._candidate_junctions(candidate, cluster))
                filters.append({"start": candidate.start, "false_junction": false, "quality": quality})
                if not false and quality:
                    accepted.append(candidate)
        arbitrated = ins._arbitrate_cluster(accepted, edges)
        converted = [ins._as_supporting_junction(c, cluster) or c for c in arbitrated]
        emitted = [c for c in converted if (c.left_junction_reads and c.right_junction_reads)
                   or c.tsd == "supporting_junction"]
        archived = []
        table = raw / "results/ALL.riceTElib.all_nonref_insert.all.txt"
        for line in table.read_text().splitlines():
            cols = line.split("\t")
            start = int(cols[4].split("..")[0])
            if cols[3] == "Chr1" and cluster.lo <= start <= cluster.hi:
                archived.append(line)
        observed = [(c.te_name, c.tsd, f"{c.start}..{c.end}", c.left_junction_reads,
                     c.right_junction_reads) for c in emitted]
        expected = [(c[0], c[1], c[4], int(c[8][2:]), int(c[7][2:]))
                    for c in (line.split("\t") for line in archived)]
        assert sorted(observed) == sorted(expected), (observed, expected)
        r2_path = benchmark / "runs/ricetelib/relocate2" / sample / "raw/repeat/bwa_aln/MSU_r7.repeat.bwa.sorted.bam"
        with pysam.AlignmentFile(str(r2_path), "rb") as r2:
            r2_records = [r.to_string() for r in r2.fetch("Chr1", lo, hi)]
        r2_table = r2_path.parent.parent / "results/Chr1.repeat.all_nonref_insert.txt"
        r2_calls = [line for line in r2_table.read_text().splitlines()
                    if cluster.lo <= int(line.split("\t")[4].split("..")[0]) <= cluster.hi]
        return {
            "sample": sample, "truth_position": position, "chrom": "Chr1",
            "cluster_lo": cluster.lo, "cluster_hi": cluster.hi,
            "query_0based_half_open": [lo, hi], "mismatch_allow": 2,
            "junctions": [asdict(j) for j in cluster.junctions],
            "support": cluster.support, "read_repeat": mapping, "reference_edges": edges,
            "raw_candidates": [asdict(c) for c in raw_candidates],
            "pooled_candidates": [asdict(c) for c in pooled], "filters": filters,
            "arbitrated": [asdict(c) for c in arbitrated],
            "emitted": [asdict(c) for c in emitted], "archived_R3_calls": archived,
            "archived_R2_calls": r2_calls, "R2_window_sam": r2_records,
            "R3_window_sam": [r.to_string() for r in records],
            "provenance": {
                "R3_bam": str(path), "R2_bam": str(r2_path),
                "R3_bam_size": path.stat().st_size,
                "R3_bam_index_sha256": hashlib.sha256(Path(str(path) + ".bai").read_bytes()).hexdigest(),
                "read_repeat_sha256": digest.hexdigest(),
                "R3_calls_sha256": hashlib.sha256(table.read_bytes()).hexdigest(),
                "R2_calls_sha256": hashlib.sha256(r2_table.read_bytes()).hexdigest(),
                "insertions_source_sha256": hashlib.sha256(Path(ins.__file__).read_bytes()).hexdigest(),
                "extractor_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "pysam_version": pysam.__version__,
                "note": "BAM windows retained as SAM; whole BAMs were not hashed. Filters recorded from indexed source queries.",
            },
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.outdir.exists():
        raise FileExistsError(args.outdir)
    data = [reconstruct(args.benchmark, sample, position) for sample, position in
            (("cov30x_rep1", 24775704), ("cov30x_rep2", 9006303))]
    args.outdir.mkdir(parents=True, exist_ok=False)
    for fixture in data:
        with (args.outdir / (fixture["sample"] + ".json")).open("x") as handle:
            json.dump(fixture, handle, indent=2, sort_keys=True)
            handle.write("\n")
        print(fixture["sample"], len(fixture["junctions"]), "junctions;",
              len(fixture["support"]), "supporting reads;",
              len(fixture["emitted"]), "reproduced calls", flush=True)


if __name__ == "__main__":
    main()
