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


def test_fixture_has_200_events():
    assert len(mpc.read_truth_gff(GFF)) == 200
