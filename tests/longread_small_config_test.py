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


def test_write_config_dir_on_fixture(tmp_path):
    out = tmp_path / "config"
    mpc.write_config_dir(GFF, GENOME, TE_LIB, out)
    assert sorted(p.name for p in out.iterdir()) == [
        "backbone.fa",
        "panel.json",
        "te_library.fa",
        "te_library_changes.tsv",
        "truth_events.tsv",
    ]
    import json

    config = json.loads((out / "panel.json").read_text())
    assert len(config["events"]) == 200
    assert len(config["te_rules"]) == 14
    truth = (out / "truth_events.tsv").read_text().splitlines()
    assert truth[0].split("\t") == [
        "event_id",
        "chrom",
        "position",
        "te_family",
        "strand",
        "tsd",
        "tsd_length",
        "biological_class",
    ]
    assert len(truth) == 201
    changes = (out / "te_library_changes.tsv").read_text().splitlines()
    assert len(changes) == 4  # header + the three ambiguity codes
    assert set("".join(mpc.read_fasta(out / "te_library.fa").values())) <= set("ACGT")
    assert not (tmp_path / "config.partial").exists()


def test_write_config_dir_refuses_existing(tmp_path):
    out = tmp_path / "config"
    out.mkdir()
    with pytest.raises(FileExistsError):
        mpc.write_config_dir(GFF, GENOME, TE_LIB, out)


def test_write_config_dir_refuses_existing_partial(tmp_path):
    out = tmp_path / "config"
    (tmp_path / "config.partial").mkdir()
    with pytest.raises(FileExistsError, match="partial"):
        mpc.write_config_dir(GFF, GENOME, TE_LIB, out)
    assert not out.exists()


def _small_inputs(tmp_path, gff_row):
    """Write a 12 bp genome, a one-family TE library and a one-row GFF."""
    genome = tmp_path / "g.fa"
    genome.write_text(">Chr3\nAAACCCGGGTTT\n")
    te_lib = tmp_path / "te.fa"
    te_lib.write_text(">mPing\nACGT\n")
    gff = tmp_path / "t.gff"
    gff.write_text(gff_row)
    return gff, genome, te_lib


def test_write_config_dir_small_inputs_succeed(tmp_path):
    gff, genome, te_lib = _small_inputs(
        tmp_path, "Chr3\tMSU7\tmPing\t6\t7\t.\t+\t.\tID=e1;TSD=CCC;\n"
    )
    out = tmp_path / "config"
    mpc.write_config_dir(gff, genome, te_lib, out)
    assert (out / "panel.json").exists()


@pytest.mark.parametrize(
    ("gff_row", "message"),
    [
        (
            "ChrX\tMSU7\tmPing\t6\t7\t.\t+\t.\tID=e1;TSD=CCC;\n",
            r"e1: chromosome 'ChrX' not in .*g\.fa",
        ),
        (
            "Chr3\tMSU7\tmPing\t2\t3\t.\t+\t.\tID=e1;TSD=CCC;\n",
            r"e1: TSD .* outside Chr3",
        ),
        (
            "Chr3\tMSU7\tmPing\t13\t14\t.\t+\t.\tID=e1;TSD=TTT;\n",
            r"e1: TSD .* outside Chr3",
        ),
        (
            "Chr3\tMSU7\tmPing\t6\t7\t.\t+\t.\tID=e1;TSD=GGG;\n",
            r"t\.gff: e1: TSD 'GGG' does not match the reference 'CCC'",
        ),
        (
            "Chr3\tMSU7\tPong\t6\t7\t.\t+\t.\tID=e1;TSD=CCC;\n",
            r"families not in .*te\.fa: Pong",
        ),
    ],
)
def test_write_config_dir_rejects_bad_input_before_writing(tmp_path, gff_row, message):
    gff, genome, te_lib = _small_inputs(tmp_path, gff_row)
    out = tmp_path / "config"
    with pytest.raises(ValueError, match=message):
        mpc.write_config_dir(gff, genome, te_lib, out)
    assert not out.exists()
    assert not (tmp_path / "config.partial").exists()


def test_read_fasta_rejects_record_without_sequence(tmp_path):
    fasta = tmp_path / "t.fa"
    fasta.write_text(">A\n>B\nACGT\n")
    with pytest.raises(ValueError, match=r"t\.fa: record 'A' has no sequence"):
        mpc.read_fasta(fasta)


def test_normalize_te_library_unsupported_base_names_record_and_position(tmp_path):
    src = tmp_path / "te.fa"
    src.write_text(">A\nAC\nGNT\n")
    with pytest.raises(ValueError) as info:
        mpc.normalize_te_library(src, tmp_path / "out.fa")
    message = str(info.value)
    assert "'A'" in message
    assert "position 4" in message
    assert "'N'" in message
    assert "(N?)" not in message


def test_main_writes_config(tmp_path):
    out = tmp_path / "config"
    argv = [
        "--gff",
        str(GFF),
        "--genome",
        str(GENOME),
        "--te-library",
        str(TE_LIB),
        "--outdir",
        str(out),
    ]
    assert mpc.main(argv) == 0
    assert (out / "panel.json").exists()


def test_write_config_dir_rejects_unsupported_te_base_before_writing(
    tmp_path, monkeypatch
):
    gff, genome, te_lib = _small_inputs(
        tmp_path, "Chr3\tMSU7\tmPing\t6\t7\t.\t+\t.\tID=e1;TSD=CCC;\n"
    )
    te_lib.write_text(">mPing\nACNT\n")
    out = tmp_path / "config"

    def no_write(*args, **kwargs):
        raise AssertionError("wrote output before validation finished")

    monkeypatch.setattr(mpc.shutil, "copyfile", no_write)
    with pytest.raises(
        ValueError, match=r"te\.fa: record 'mPing' position 3: unsupported base 'N'"
    ):
        mpc.write_config_dir(gff, genome, te_lib, out)
    assert not out.exists()
    assert not (tmp_path / "config.partial").exists()


def test_write_config_dir_removes_partial_on_write_error(tmp_path, monkeypatch):
    gff, genome, te_lib = _small_inputs(
        tmp_path, "Chr3\tMSU7\tmPing\t6\t7\t.\t+\t.\tID=e1;TSD=CCC;\n"
    )
    out = tmp_path / "config"

    def disk_full(*args, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(mpc.shutil, "copyfile", disk_full)
    with pytest.raises(OSError, match="disk full"):
        mpc.write_config_dir(gff, genome, te_lib, out)
    assert not out.exists()
    assert not (tmp_path / "config.partial").exists()
