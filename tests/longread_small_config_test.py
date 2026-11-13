"""Tests for the small long-read panel config builder."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "validation" / "longread_small"))

import make_panel_config as mpc  # noqa: E402

DATA = ROOT / "tests" / "data"
GFF = DATA / "sim_genome" / "MSU7.Chr3_2M.ALL.gff"
GENOME = DATA / "sim_genome" / "MSU7.Chr3_2M.fa"
TE_LIB = DATA / "TE_lib" / "RiceTE.fa"


def test_read_truth_gff_parses_one_row(tmp_path):
    gff = tmp_path / "t.gff"
    gff.write_text("Chr3\tMSU7\tmPing\t100\t101\t.\t-\t.\tID=mPing.Chr3.100;TSD=TTA;\n")
    events = mpc.read_truth_gff(gff)
    assert events == [
        mpc.TruthEvent(
            event_id="mPing.Chr3.100",
            chrom="Chr3",
            position=100,
            family="mPing",
            strand="-",
            tsd="TTA",
        )
    ]


def test_read_truth_gff_rejects_missing_tsd(tmp_path):
    gff = tmp_path / "t.gff"
    gff.write_text("Chr3\tMSU7\tmPing\t100\t101\t.\t+\t.\tID=mPing.Chr3.100;\n")
    with pytest.raises(ValueError, match="TSD"):
        mpc.read_truth_gff(gff)


@pytest.mark.parametrize("strand", [".", "", "+-"])
def test_read_truth_gff_rejects_bad_strand(tmp_path, strand):
    gff = tmp_path / "t.gff"
    gff.write_text(
        f"Chr3\tMSU7\tmPing\t100\t101\t.\t{strand}\t.\tID=mPing.Chr3.100;TSD=TTA;\n"
    )
    with pytest.raises(ValueError, match="strand"):
        mpc.read_truth_gff(gff)


@pytest.mark.parametrize(
    ("attrs", "message"),
    [
        ("ID=mPing.Chr3.100;TSD=TTN;", "TSD"),
        ("ID=mPing.Chr3.100;TSD=tta;", "TSD"),
        ("Parent_ID=x;TSD=TTA;", "ID"),
        ("ID=mPing.Chr3.100;XTSD=TTA;", "TSD"),
    ],
)
def test_read_truth_gff_rejects_malformed_attributes(tmp_path, attrs, message):
    gff = tmp_path / "t.gff"
    gff.write_text(f"Chr3\tMSU7\tmPing\t100\t101\t.\t+\t.\t{attrs}\n")
    with pytest.raises(ValueError, match=message):
        mpc.read_truth_gff(gff)


def test_read_truth_gff_rejects_non_integer_position(tmp_path):
    gff = tmp_path / "t.gff"
    gff.write_text("Chr3\tMSU7\tmPing\tabc\t101\t.\t+\t.\tID=m;TSD=TTA;\n")
    with pytest.raises(ValueError, match=r"t\.gff:1"):
        mpc.read_truth_gff(gff)


@pytest.mark.parametrize(
    ("attrs", "message"),
    [
        ("ID=a;ID=b;TSD=TTA;", "ID attribute given 2 times"),
        ("ID=a;TSD=TTA;TSD=GG;", "TSD attribute given 2 times"),
        ("ID= ;TSD=TTA;", "missing ID attribute"),
    ],
)
def test_read_truth_gff_rejects_repeated_or_blank_keys(tmp_path, attrs, message):
    gff = tmp_path / "t.gff"
    gff.write_text(f"Chr3\tMSU7\tmPing\t100\t101\t.\t+\t.\t{attrs}\n")
    with pytest.raises(ValueError, match=message):
        mpc.read_truth_gff(gff)


def test_read_truth_gff_rejects_duplicate_event_id(tmp_path):
    gff = tmp_path / "t.gff"
    row = "Chr3\tMSU7\tmPing\t{}\t101\t.\t+\t.\tID=same;TSD=TTA;\n"
    gff.write_text(row.format(100) + row.format(200))
    with pytest.raises(ValueError, match=r"t\.gff:2: duplicate ID 'same'"):
        mpc.read_truth_gff(gff)


@pytest.mark.parametrize("position", ["0", "-5"])
def test_read_truth_gff_rejects_position_below_one(tmp_path, position):
    gff = tmp_path / "t.gff"
    gff.write_text(f"Chr3\tMSU7\tmPing\t{position}\t101\t.\t+\t.\tID=m;TSD=TTA;\n")
    with pytest.raises(ValueError, match="bad position"):
        mpc.read_truth_gff(gff)


@pytest.mark.parametrize("ncols", [8, 10])
def test_read_truth_gff_requires_exactly_nine_columns(tmp_path, ncols):
    fields = ["Chr3", "MSU7", "mPing", "100", "101", ".", "+", ".", "ID=m;TSD=TTA;"]
    fields = (fields + ["extra"])[:ncols]
    gff = tmp_path / "t.gff"
    gff.write_text("\t".join(fields) + "\n")
    with pytest.raises(ValueError, match=f"expected 9 columns, got {ncols}"):
        mpc.read_truth_gff(gff)


def test_fixture_has_200_events():
    assert len(mpc.read_truth_gff(GFF)) == 200


def test_normalize_te_library_replaces_ambiguity_codes(tmp_path):
    src = tmp_path / "te.fa"
    src.write_text(">A desc\nacgS\nMR\n>B\nACGT\n")
    dst = tmp_path / "out.fa"
    changes = mpc.normalize_te_library(src, dst)
    assert dst.read_text() == ">A\nACGCAA\n>B\nACGT\n"
    assert changes == [("A", 4, "S", "C"), ("A", 5, "M", "A"), ("A", 6, "R", "A")]


def test_normalize_te_library_rejects_n(tmp_path):
    src = tmp_path / "te.fa"
    src.write_text(">A\nACNT\n")
    with pytest.raises(ValueError, match="N"):
        mpc.normalize_te_library(src, tmp_path / "out.fa")


@pytest.mark.parametrize("text", [">\nACGT\n", ">  \nACGT\n"])
def test_read_fasta_rejects_header_without_name(tmp_path, text):
    fasta = tmp_path / "t.fa"
    fasta.write_text(text)
    with pytest.raises(ValueError, match=r"t\.fa:1: header without a name"):
        mpc.read_fasta(fasta)


def test_read_fasta_rejects_empty_file(tmp_path):
    fasta = tmp_path / "t.fa"
    fasta.write_text("")
    with pytest.raises(ValueError, match=r"t\.fa: no FASTA records"):
        mpc.read_fasta(fasta)


def test_build_config_shape():
    events = [
        mpc.TruthEvent("mPing.Chr3.100", "Chr3", 100, "mPing", "-", "TTA"),
        mpc.TruthEvent("karma.Chr3.900", "Chr3", 900, "karma", "+", "A" * 12),
    ]
    config = mpc.build_config(events)
    assert config["schema_version"] == 1
    assert config["backbone"] == "backbone.fa"
    assert config["te_library"] == "te_library.fa"
    assert config["chromosomes"] == ["Chr3"]
    assert config["genomes"] == [
        {"id": "reference", "sample": "Reference", "haplotype": 1, "role": "reference"},
        {"id": "sample_h1", "sample": "Sample", "haplotype": 1, "role": "sample"},
    ]
    assert config["te_rules"] == {
        "karma": {"family": "karma", "tsd_length": 12},
        "mPing": {"family": "mPing", "tsd_length": 3},
    }
    assert config["events"][0] == {
        "id": "mPing.Chr3.100",
        "chromosome": "Chr3",
        "position": 100,
        "te_id": "mPing",
        "strand": "-",
        "carriers": ["sample_h1"],
    }


def test_build_config_rejects_mixed_tsd_length():
    events = [
        mpc.TruthEvent("a", "Chr3", 100, "mPing", "+", "TTA"),
        mpc.TruthEvent("b", "Chr3", 900, "mPing", "+", "TTAA"),
    ]
    with pytest.raises(ValueError, match="mPing"):
        mpc.build_config(events)


def test_fixture_tsd_matches_reference():
    genome = mpc.read_fasta(GENOME)
    for event in mpc.read_truth_gff(GFF):
        start = event.position - len(event.tsd)
        assert genome[event.chrom][start : event.position] == event.tsd, event
