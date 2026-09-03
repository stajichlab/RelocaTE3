"""TE-family evidence survives RelocaTE3's structured TXT/GFF outputs."""

from RelocaTE3.characterize import write_characterized

from RelocaTE3.insertions import (
    _Cluster,
    _count_support,
    read_insertions_gff,
    write_insertions_gff,
    write_insertions_txt,
)
from RelocaTE3.models import Insertion


def _insertion() -> Insertion:
    return Insertion(
        chrom="Chr1",
        start=100,
        end=102,
        te_name="mPing",
        strand="+",
        tsd="TTA",
        left_junction_reads=2,
        right_junction_reads=1,
        te_family_support={"mPing": 2, "RIRE3": 1},
        te_family_confidence=2 / 3,
        te_family_status="dominant",
        te_supporting_family_support={"RIRE3": 2},
        te_supporting_family_confidence=1.0,
        te_supporting_family_status="unique",
        te_family_concordance="discordant",
    )


def test_family_metadata_round_trips_through_gff(tmp_path):
    path = tmp_path / "insertions.gff"
    write_insertions_gff([_insertion()], path, sample="HEG4")

    [observed] = read_insertions_gff(path)
    assert observed.te_name == "mPing"
    assert observed.te_family_support == {"mPing": 2, "RIRE3": 1}
    assert observed.te_family_confidence == 0.666667
    assert observed.te_family_status == "dominant"
    assert observed.te_supporting_family_support == {"RIRE3": 2}
    assert observed.te_supporting_family_confidence == 1.0
    assert observed.te_supporting_family_status == "unique"
    assert observed.te_family_concordance == "discordant"


def test_family_metadata_is_appended_to_structured_txt(tmp_path):
    path = tmp_path / "insertions.txt"
    write_insertions_txt([_insertion()], path)

    header, row = [line.split("\t") for line in path.read_text().splitlines()]
    assert header[-7:] == [
        "TE_family_support",
        "TE_family_confidence",
        "TE_family_status",
        "TE_supporting_family_support",
        "TE_supporting_family_confidence",
        "TE_supporting_family_status",
        "TE_family_concordance",
    ]
    assert row[-7:] == [
        "mPing=2,RIRE3=1",
        "0.666667",
        "dominant",
        "RIRE3=2",
        "1.000000",
        "unique",
        "discordant",
    ]


def test_supporting_family_is_separate_and_cannot_change_primary():
    cluster = _Cluster("Chr1")
    # Two alignments with one read name retain legacy alignment counts but cast
    # only one family vote, so a multimapper cannot inflate family support.
    cluster.support.extend(
        [
            ("pairA", 10, 90, "+", "ACGT"),
            ("pairA", 20, 95, "+", "ACGT"),
            ("pairB", 210, 300, "-", "ACGT"),
        ]
    )
    ins = _insertion()
    _count_support(
        ins,
        cluster,
        {
            "pairA/1": ("RIRE3", "+"),
            "pairB/2": ("RIRE3", "+"),
        },
    )

    assert ins.te_name == "mPing"
    assert ins.left_support_reads == 2
    assert ins.right_support_reads == 1
    assert ins.te_supporting_family_support == {"RIRE3": 2}
    assert ins.te_supporting_family_confidence == 1.0
    assert ins.te_supporting_family_status == "unique"
    assert ins.te_family_concordance == "discordant"


def test_older_gff_without_supporting_family_fields_remains_readable(tmp_path):
    path = tmp_path / "legacy.gff"
    path.write_text(
        "Chr1\tRelocaTE3\tHEG4\t100\t102\t.\t+\t.\t"
        "ID=x;Name=mPing;TSD=TTA;TE_family_support=mPing=2;"
        "TE_family_confidence=1.0;TE_family_status=unique;\n"
    )

    [observed] = read_insertions_gff(path)
    assert observed.te_name == "mPing"
    assert observed.te_supporting_family_support == {}
    assert observed.te_supporting_family_confidence == 0.0
    assert observed.te_supporting_family_status == "unassigned"
    assert observed.te_family_concordance == "unassigned"


def test_supporting_family_metadata_survives_object_characterization(tmp_path):
    ins = _insertion()
    ins.status = "heterozygous"
    ins.spanners = 3
    txt = tmp_path / "characterized.txt"
    gff = tmp_path / "characterized.gff"

    write_characterized([ins], gff, txt, sample="HEG4")

    header, row = [line.split("\t") for line in txt.read_text().splitlines()]
    assert header[-7:] == [
        "TE_family_support",
        "TE_family_confidence",
        "TE_family_status",
        "TE_supporting_family_support",
        "TE_supporting_family_confidence",
        "TE_supporting_family_status",
        "TE_family_concordance",
    ]
    assert row[-4:] == ["RIRE3=2", "1.000000", "unique", "discordant"]
    attributes = gff.read_text()
    assert "TE_supporting_family_support=RIRE3=2;" in attributes
    assert "TE_family_concordance=discordant" in attributes
