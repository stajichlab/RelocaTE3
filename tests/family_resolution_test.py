"""Independent evidence resolves real read ties while sparse/conflicting reads do not."""

import importlib.util
import io
import json
from pathlib import Path

import pysam
import pytest
from RelocaTE3.family import resolve_family
from RelocaTE3.insertions import (
    InsertionFinder,
    _call_insertions,
    _stream_clusters,
    _mapping_family_evidence,
    _te_family_metadata_columns,
    read_insertions_gff,
    write_insertions_gff,
    _row_to_gff,
)
from RelocaTE3.librelocate import RelocaTE
from RelocaTE3.characterize import Characterizer, write_characterized

root = Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location(
    "line_analysis", root / "scripts/analyze_line_family_trace.py"
)
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def entry(primary, families=(), compatible=True):
    return primary, families, compatible


def test_unique_evidence_resolves_compatible_read_ties():
    evidence = resolve_family(
        [entry("Z"), entry("Z"), entry("A")] + [entry("A", ("A", "Z"))] * 8
    )
    assert evidence.primary == "Z"
    assert evidence.support == {"Z": 2, "A": 1}
    assert evidence.ambiguous_reads == 8
    assert evidence.status == "ambiguous"
    assert evidence.confidence == pytest.approx(2 / 11)
    assert evidence.resolution == "unique_junction_consensus"
    assert evidence.candidate_support == {"Z": 10, "A": 9}


@pytest.mark.parametrize(
    "reads",
    [
        [entry("Z")] + [entry("A", ("A", "Z"))] * 8,
        [entry("A", ("A", "Z"))] * 8,
        [entry("Z")] * 2 + [entry("A")] * 2 + [entry("A", ("A", "Z"))] * 8,
        [entry("Z")] * 2 + [entry("A", ("A", "B"))] * 8,
        [entry("Z")] * 2 + [entry("A", ("A", "Z"), False)] * 8,
    ],
)
def test_sparse_conflicting_incompatible_evidence_does_not_resolve(reads):
    evidence = resolve_family(reads)
    assert evidence.resolution == "unresolved_ties"
    assert evidence.status == "ambiguous"
    assert evidence.primary == "A"


def test_legacy_selected_votes_and_missing_families_are_unchanged():
    evidence = resolve_family([entry("Z")] * 3 + [entry("A"), entry("NA")])
    assert (evidence.primary, evidence.confidence, evidence.status) == (
        "Z",
        0.75,
        "dominant",
    )
    assert evidence.resolution == "selected_votes" and evidence.ambiguous_reads == 0


def test_duplicate_names_cannot_supply_two_independent_votes():
    mapping = {"a": ["Z", "+"], "tie": ["A", "+", ("A", "Z"), True]}
    evidence = _mapping_family_evidence(
        mapping, ["a:end:5", "a:end:5", "tie:start:3"] * 3
    )
    assert evidence.support == {"Z": 1}
    assert evidence.resolution == "unresolved_ties"


def test_unresolved_ties_preserve_pre_deduplication_primary():
    # Minimal reproduction of TE000496: archived support 9:7, distinct
    # selected reads 5:7. Incompatible alternatives must not change the label.
    mapping = {f"a{i}": ["Os1524", "+"] for i in range(4)}
    mapping.update({f"b{i}": ["Os3293", "+"] for i in range(7)})
    mapping["tie"] = ["Os1524", "+", ("Os1524", "Os3293"), False]
    names = list(mapping) + [f"a{i}" for i in range(4)]
    evidence = _mapping_family_evidence(mapping, names)
    assert evidence.primary == "Os1524"
    assert evidence.resolution == "unresolved_ties"
    assert evidence.support == {"Os3293": 7, "Os1524": 4}
    assert evidence.confidence == pytest.approx(4 / 12)
    assert evidence.ambiguous_reads == 1
    assert evidence.candidate_support == {"Os3293": 8, "Os1524": 5}


def test_best_family_sets_survive_hit_order_and_better_score_reset():
    hits = analysis.select_psl(
        json.loads((root / "tests/data/line_family/TE000172.json").read_text())
    )
    records = next(
        [r for r in rs if r["admitted"]]
        for rs in hits.values()
        if len(
            {
                r["tName"]
                for r in rs
                if r["admitted"]
                and (r["boundary"], r["match"])
                == max((h["boundary"], h["match"]) for h in rs if h["admitted"])
            }
        )
        > 1
    )
    results = []
    for order in (records, list(reversed(records))):
        winner = None
        for rec in order:
            rec = dict(rec)
            winner = rec if winner is None else RelocaTE._merge_te_hits(winner, rec)
        results.append(winner)
    assert results[0]["tName"] == results[1]["tName"]
    assert results[0]["best_families"] == results[1]["best_families"]
    assert results[0]["family_ties_compatible"] == results[1]["family_ties_compatible"]
    better = dict(records[0], boundary=5, tName="better")
    selected = RelocaTE._merge_te_hits(results[0], better)
    assert selected["tName"] == "better" and "best_families" not in selected


def test_incompatible_query_geometry_is_retained_as_unresolved():
    hit = dict(
        boundary=2,
        match=40,
        tName="Z",
        start=0,
        end=40,
        strand="+",
        tStart=0,
        tEnd=40,
        tLen=100,
        mismatch=0,
    )
    other = dict(hit, tName="A", start=1, end=41)
    merged = RelocaTE._merge_te_hits(hit, other)
    assert merged["best_families"] == ("A", "Z")
    assert not merged["family_ties_compatible"]


def test_invalid_mapping_evidence_fails_explicitly(tmp_path):
    path = tmp_path / "rr.tsv"
    path.write_text('r\tA\t+\tTE_best_families:{"families":["B"],"compatible":true}\n')
    with pytest.raises(ValueError):
        InsertionFinder._load_read_repeat(path)


def test_real_line_reads_resolve_through_mapping_clustering_and_outputs(tmp_path):
    evidence = json.loads((root / "tests/data/line_family/TE000172.json").read_text())
    alignments = analysis.replay_bams(evidence, tmp_path)
    mapping_path = tmp_path / "rr.tsv"
    out = io.StringIO()
    for name, rec in alignments.items():
        RelocaTE._write_family_mapping(out, name, rec)
    mapping_path.write_text(out.getvalue())
    mapping = InsertionFinder._load_read_repeat(mapping_path)
    assert sum(len(row) > 2 for row in mapping.values()) == 8
    genome = json.loads(
        (root / "tests/data/line_family/TE000172.junctions.json").read_text()
    )
    # A synthetic header supplies the contig; SAM records themselves are archived.
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "Chr1", "LN": 50000000}]}
    )
    records = [
        pysam.AlignedSegment.fromstring(s, header)
        for s in genome["genome_junction_sam"]
    ]
    path = tmp_path / "junctions.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for rec in sorted(records, key=lambda r: r.reference_start):
            bam.write(rec)
    [cluster] = list(_stream_clusters(str(path), mapping))
    [call] = _call_insertions(cluster, None, mapping)
    assert (call.start, call.end, call.tsd) == (29402710, 29402725, "GGAGAGGGAGGTGGCC")
    assert (call.left_junction_reads, call.right_junction_reads) == (2, 9)
    assert call.te_name == "Os3328#LINE/unknown"
    assert call.te_family_support == {
        "Os3328#LINE/unknown": 2,
        "Os0596#LINE/unknown": 1,
    }
    assert call.te_family_ambiguous_reads == 8
    assert call.te_family_confidence == pytest.approx(2 / 11)
    assert call.te_family_resolution == "unique_junction_consensus"
    # Old mappings preserve their previous labels rather than inventing alternatives.
    legacy_mapping = {name: row[:2] for name, row in mapping.items()}
    [legacy_cluster] = list(_stream_clusters(str(path), legacy_mapping))
    [legacy_call] = _call_insertions(legacy_cluster, None, legacy_mapping)
    assert legacy_call.te_name == "Os0596#LINE/unknown"
    assert legacy_call.te_family_resolution == "selected_votes"
    names = [j.read_name for j in cluster.junctions]
    fixed = InsertionFinder._insertion_family_evidence(names, mapping)
    assert fixed.primary == call.te_name and fixed.ambiguous_reads == 8
    gff = tmp_path / "calls.gff"
    write_insertions_gff([call], gff, sample="test")
    [restored] = read_insertions_gff(gff)
    assert restored.te_family_ambiguous_reads == 8
    assert restored.te_family_resolution == call.te_family_resolution
    assert restored.te_family_candidate_support == {
        "Os3328#LINE/unknown": 10,
        "Os0596#LINE/unknown": 9,
    }
    call.status = "homozygous"
    write_characterized(
        [call], tmp_path / "char.gff", tmp_path / "char.txt", sample="test"
    )
    header_names, row = [
        s.split("\t") for s in (tmp_path / "char.txt").read_text().splitlines()
    ]
    assert dict(zip(header_names, row))["TE_family_ambiguous_reads"] == "8"
    metadata = Characterizer._extract_family_metadata(
        [
            "TE_family_ambiguous_reads:8",
            "TE_family_resolution:unique_junction_consensus",
        ]
    )
    assert metadata["family_ambiguous_reads"] == "8"
    assert metadata["family_resolution"] == call.te_family_resolution
    fields = [
        call.te_name,
        call.tsd,
        "test",
        call.chrom,
        f"{call.start}..{call.end}",
        call.strand,
        "T:11",
        "R:9",
        "L:2",
        "ST:0",
        "SR:0",
        "SL:0",
    ]
    tier_gff = _row_to_gff(fields + list(_te_family_metadata_columns(call)), "test")
    assert "TE_family_ambiguous_reads=8;" in tier_gff
    assert "TE_family_resolution=unique_junction_consensus;" in tier_gff
    # Stored-alignment replays retain only requested reads but preserve their evidence.
    left = RelocaTE()._parse_te_bam(tmp_path / "left.bam")
    name = next(iter(left))
    filtered = RelocaTE()._parse_te_bam(tmp_path / "left.bam", read_names={name})
    assert filtered == {name: left[name]}
