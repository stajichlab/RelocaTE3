"""L3 gate scorer: long-read calls against a truth table of insertions.

``validation/longread_l3/score_calls.py`` matches PASS spanning calls to
truth sites one-to-one (greedy by distance, same sequence and family) and
reports recall, precision, F1, strand accuracy and diagnostics.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from RelocaTE3.longread.model import Call, CallTier, FilterReason, Locus, ReadType
from RelocaTE3.longread.tables import write_calls

REPO = Path(__file__).resolve().parents[1]
SCORER_DIR = REPO / "validation" / "longread_l3"

sys.path.insert(0, str(SCORER_DIR))
import score_calls  # noqa: E402

COLUMNS = [
    "label", "n_truth", "n_pass", "tp", "recall", "precision", "f1",
    "tp_position_only", "strand_correct", "strand_accuracy", "n_filtered",
    "clip_only_pass", "clip_only_hits", "control_pass",
    "control_clip_only_pass",
]  # fmt: skip

# Three truth sites on Chr3 (0-based insertion points): (chrom, pos, family, strand).
TRUTH = [
    ("Chr3", 1000, "mPing", "+"),
    ("Chr3", 5000, "Dasheng", "-"),
    ("Chr3", 9000, "RIRE2", "+"),
]

SPAN = CallTier.SPANNING
CLIP = CallTier.CLIP_ONLY
LOW = (FilterReason.LOW_SUPPORT,)


def _truth(tmp_path: Path, sites=TRUTH) -> Path:
    path = tmp_path / "truth_events.tsv"
    lines = ["event_id\tchrom\tposition\tte_family\tstrand\ttsd\ttsd_length\t"
             "biological_class\n"]  # fmt: skip
    for i, (chrom, pos, family, strand) in enumerate(sites):
        lines.append(f"e{i}\t{chrom}\t{pos}\t{family}\t{strand}\tTAA\t3\thomozygous\n")
    path.write_text("".join(lines))
    return path


def _call(
    i: int,
    chrom: str,
    pos: int,
    family: str,
    strand: str = "+",
    tier: CallTier = SPAN,
    filters: tuple[FilterReason, ...] = (),
    reference_id: str = "ref",
) -> Call:
    return Call(
        call_id=f"call{i}", locus=Locus(reference_id, chrom, pos, pos, "+"),
        tier=tier, filters=filters, family=family, te_strand=strand,
        family_fraction=1.0, strand_fraction=1.0, te_read_fraction=1.0,
        n_reads=3, n_ins=3 if tier is SPAN else 0, n_split_ins=0,
        n_clip_left=0 if tier is SPAN else 2, n_clip_right=0 if tier is SPAN else 1,
        median_insert_length=500, median_te_coverage=0.95, median_identity=0.99,
        nested=False, cluster_ids=(f"c{i}",),
    )  # fmt: skip


def _calls(tmp_path: Path, name: str, specs: list[tuple]) -> Path:
    """Write calls from ``(chrom, pos, family[, strand, tier, filters, ref])``."""
    path = tmp_path / f"{name}.calls.tsv.gz"
    calls = [_call(i, *spec) for i, spec in enumerate(specs)]
    write_calls(path, calls, read_type=ReadType.HIFI, reference_id="ref")
    return path


def _perfect(tmp_path: Path) -> Path:
    return _calls(tmp_path, "s", [(c, p, f, s) for c, p, f, s in TRUTH])


def test_perfect_calls(tmp_path):
    row = score_calls.score(_truth(tmp_path), _perfect(tmp_path), None, 50, "x")
    assert row == {
        "label": "x", "n_truth": 3, "n_pass": 3, "tp": 3, "recall": "1.0000",
        "precision": "1.0000", "f1": "1.0000", "tp_position_only": 3,
        "strand_correct": 3, "strand_accuracy": "1.0000", "n_filtered": 0,
        "clip_only_pass": 0, "clip_only_hits": 0, "control_pass": "",
        "control_clip_only_pass": "",
    }  # fmt: skip


def test_wrong_family_counts_only_position(tmp_path):
    calls = _calls(tmp_path, "s", [("Chr3", 1010, "Dasheng")])
    row = score_calls.score(_truth(tmp_path), calls, None, 50)
    assert row["n_pass"] == 1
    assert row["tp"] == 0
    assert row["tp_position_only"] == 1
    assert row["recall"] == "0.0000"
    assert row["precision"] == "0.0000"
    assert row["f1"] == "0.0000"  # 2*0 / (3 + 1)
    assert row["strand_accuracy"] == ""  # no matched pairs


def test_two_calls_near_one_site_one_tp(tmp_path):
    calls = _calls(tmp_path, "s", [("Chr3", 1005, "mPing"), ("Chr3", 1020, "mPing")])
    row = score_calls.score(_truth(tmp_path), calls, None, 50)
    assert row["tp"] == 1
    assert row["precision"] == "0.5000"
    assert row["recall"] == "0.3333"
    assert row["f1"] == "0.4000"  # 2*1 / (3 + 2)


def test_call_between_two_sites_matches_the_nearer(tmp_path):
    truth = _truth(tmp_path, [("Chr3", 1000, "mPing", "+"),
                              ("Chr3", 1060, "mPing", "+")])  # fmt: skip
    # Call 0 at 1040 is 40 bp from site 0 and 20 bp from site 1; call 1 at
    # 1095 reaches only site 1 (35 bp). Pairs by distance: (20, s1, c0),
    # (35, s1, c1), (40, s0, c0): c0 takes s1, then c1 and s0 find their
    # partner used.
    calls = _calls(tmp_path, "s", [("Chr3", 1040, "mPing"), ("Chr3", 1095, "mPing")])
    truth_sites = score_calls.read_truth(truth)
    _, scored = score_calls.read_scored_calls(calls)
    assert score_calls.greedy_match(truth_sites, scored, 50) == [(1, 0)]
    row = score_calls.score(truth, calls, None, 50)
    assert row["tp"] == 1
    assert row["recall"] == "0.5000"
    assert row["precision"] == "0.5000"


def test_boundary_50_matches_51_not(tmp_path):
    calls = _calls(
        tmp_path, "s", [("Chr3", 1050, "mPing"), ("Chr3", 4949, "Dasheng", "-")]
    )
    row = score_calls.score(_truth(tmp_path), calls, None, 50)
    assert row["tp"] == 1
    assert row["tp_position_only"] == 1


def test_other_chromosome_never_matches(tmp_path):
    calls = _calls(tmp_path, "s", [("Chr1", 1000, "mPing")])
    row = score_calls.score(_truth(tmp_path), calls, None, 50)
    assert row["tp"] == 0 and row["tp_position_only"] == 0


def test_filtered_calls_are_not_scored(tmp_path):
    calls = _calls(tmp_path, "s", [
        ("Chr3", 1000, "mPing"),
        ("Chr3", 5000, "Dasheng", "-", SPAN, LOW),
        ("Chr3", 9000, "RIRE2", "+", CLIP, LOW),  # filtered clip-only: neither
    ])  # fmt: skip
    row = score_calls.score(_truth(tmp_path), calls, None, 50)
    assert row["n_pass"] == 1
    assert row["tp"] == 1
    assert row["recall"] == "0.3333"
    assert row["precision"] == "1.0000"
    assert row["f1"] == "0.5000"  # 2*1 / (3 + 1)
    assert row["n_filtered"] == 1
    assert row["clip_only_pass"] == 0


def test_clip_only_pass_at_unmatched_site_is_a_hit(tmp_path):
    calls = _calls(tmp_path, "s", [
        ("Chr3", 1000, "mPing"),                         # spanning tp
        ("Chr3", 5010, "Dasheng", "-", CLIP),            # unmatched site: hit
        ("Chr3", 1000, "mPing", "+", CLIP),              # site already matched
        ("Chr3", 9000, "mPing", "+", CLIP),              # wrong family
        ("Chr3", 9000, "RIRE2", "+", CLIP, LOW),         # not PASS
    ])  # fmt: skip
    row = score_calls.score(_truth(tmp_path), calls, None, 50)
    assert row["n_pass"] == 1
    assert row["tp"] == 1
    assert row["clip_only_pass"] == 3
    assert row["clip_only_hits"] == 1
    assert row["n_filtered"] == 0  # only spanning calls count here


def test_strand_mismatch_lowers_strand_accuracy(tmp_path):
    calls = _calls(tmp_path, "s", [
        ("Chr3", 1000, "mPing", "+"),
        ("Chr3", 5000, "Dasheng", "+"),  # truth is "-"
        ("Chr3", 9000, "RIRE2", "+"),
    ])  # fmt: skip
    row = score_calls.score(_truth(tmp_path), calls, None, 50)
    assert row["tp"] == 3
    assert row["strand_correct"] == 2
    assert row["strand_accuracy"] == "0.6667"


def test_control_columns(tmp_path):
    truth = _truth(tmp_path)
    sample = _perfect(tmp_path)
    row = score_calls.score(truth, sample, None, 50)
    assert row["control_pass"] == "" and row["control_clip_only_pass"] == ""
    control = _calls(tmp_path, "k", [
        ("Chr3", 1000, "mPing"),
        ("Chr3", 70000, "karma"),
        ("Chr3", 80000, "karma", "+", SPAN, LOW),
        ("Chr3", 90000, "karma", "+", CLIP),
        ("Chr3", 95000, "karma", "+", CLIP, LOW),
    ])  # fmt: skip
    row = score_calls.score(truth, sample, control, 50)
    assert row["control_pass"] == 2
    assert row["control_clip_only_pass"] == 1
    assert row["tp"] == 3  # the control does not change the sample's scores


def test_call_with_other_reference_id_raises(tmp_path):
    calls = _calls(tmp_path, "s", [("Chr3", 1000, "mPing", "+", SPAN, (), "other")])
    with pytest.raises(ValueError, match="other"):
        score_calls.score(_truth(tmp_path), calls, None, 50)


def test_control_with_other_reference_id_raises(tmp_path):
    control = _calls(tmp_path, "k", [("Chr3", 1000, "mPing", "+", CLIP, (), "other")])
    with pytest.raises(ValueError, match="other"):
        score_calls.score(_truth(tmp_path), _perfect(tmp_path), control, 50)


def test_per_family_file(tmp_path):
    calls = _calls(tmp_path, "s", [
        ("Chr3", 1000, "mPing"),
        ("Chr3", 5000, "Dasheng", "-"),
        ("Chr3", 7000, "karma"),                       # false positive
        ("Chr3", 9000, "RIRE2", "+", SPAN, LOW),       # filtered
    ])  # fmt: skip
    out = tmp_path / "per_family.tsv"
    score_calls.score(_truth(tmp_path), calls, None, 50, per_family=out)
    assert out.read_text().splitlines() == [
        "family\tn_truth\ttp\trecall\tn_pass\tprecision",
        "Dasheng\t1\t1\t1.0000\t1\t1.0000",
        "RIRE2\t1\t0\t0.0000\t0\t",
        "karma\t0\t0\t\t1\t0.0000",
        "mPing\t1\t1\t1.0000\t1\t1.0000",
    ]  # fmt: skip


def test_main_prints_header_and_row(tmp_path, capsys):
    truth = _truth(tmp_path)
    calls = _calls(tmp_path, "s", [
        ("Chr3", 1000, "mPing", "-"),          # tp, wrong strand
        ("Chr3", 5051, "Dasheng", "-"),        # 51 bp: no match
        ("Chr3", 9000, "RIRE2", "+", CLIP),    # clip-only hit
    ])  # fmt: skip
    control = _calls(tmp_path, "k", [("Chr3", 9000, "RIRE2")])
    per_family = tmp_path / "pf.tsv"
    rc = score_calls.main([
        "--truth", str(truth), "--calls", str(calls), "--control", str(control),
        "--label", "hifi_cov5x", "--header", "--per-family", str(per_family),
    ])  # fmt: skip
    assert rc == 0
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 2
    assert lines[0].split("\t") == COLUMNS
    # n_pass 2, tp 1: recall 1/3, precision 1/2, f1 2/5.
    assert lines[1].split("\t") == [
        "hifi_cov5x", "3", "2", "1", "0.3333", "0.5000", "0.4000", "1", "0",
        "0.0000", "0", "1", "1", "1", "0",
    ]  # fmt: skip
    assert per_family.is_file()


def test_main_without_header_or_control(tmp_path, capsys):
    calls = _calls(tmp_path, "s", [("Chr3", 5000, "Dasheng", "-")])
    rc = score_calls.main(["--truth", str(_truth(tmp_path)), "--calls", str(calls),
                           "--max-dist", "10", "--label", "x"])  # fmt: skip
    assert rc == 0
    assert capsys.readouterr().out.splitlines() == [
        "x\t3\t1\t1\t0.3333\t1.0000\t0.5000\t1\t1\t1.0000\t0\t0\t0\t\t"
    ]


def test_negative_max_dist_rejected(tmp_path):
    with pytest.raises(ValueError, match="max_dist"):
        score_calls.score(_truth(tmp_path), _perfect(tmp_path), None, -1)
