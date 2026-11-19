"""Long-read versioned TSV.gz tables: signatures and clusters."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from RelocaTE3.longread.model import (
    ClipSide, Cluster, Locus, ReadType, Signature, SignatureKind,
)  # fmt: skip
from RelocaTE3.longread.tables import (
    CLUSTER_COLUMNS, SIGNATURE_COLUMNS, read_clusters, read_signatures,
    write_clusters, write_signatures,
)  # fmt: skip

REF = "msu7"


def _signatures() -> list[Signature]:
    return [
        Signature(
            "s1", SignatureKind.INS, Locus(REF, "Chr3", 100, 100, "+"), 500,
            ClipSide.NONE, "r1", 10, 510, 60, "ACGT" * 5, "GGCC", "TTAA",
        ),
        Signature(
            "s2", SignatureKind.SPLIT_INS, Locus(REF, "Chr3", 105, 105, "-"),
            800, ClipSide.NONE, "r2", 200, 1000, 0, "", "AC", "",
        ),
        Signature(
            "s3", SignatureKind.CLIP, Locus(REF, "Chr3", 98, 98, "+"), 300,
            ClipSide.LEFT, "r3", 0, 300, 42, "TTTT", "", "GA",
        ),
    ]  # fmt: skip


def _clusters() -> list[Cluster]:
    return [
        Cluster(
            "c1", Locus(REF, "Chr3", 100, 100, "+"), 3, 1, 1, 1, 0,
            (500, 800, 300), ("s1", "s2", "s3"),
        ),
        Cluster(
            "c2", Locus(REF, "Chr3", 5000, 5000, "-"), 1, 0, 0, 0, 1,
            (), ("s9",),
        ),
    ]  # fmt: skip


def _write_raw(path: Path, lines: list[str]) -> Path:
    with gzip.open(path, "wt") as fh:
        fh.write("".join(line + "\n" for line in lines))
    return path


def _schema_line(schema: str, version: str = "v1") -> str:
    return f"#relocate3:{schema}\t{version}\tread_type=ont\treference_id={REF}"


def test_signatures_round_trip_and_metadata(tmp_path):
    path = tmp_path / "sigs.tsv.gz"
    n = write_signatures(
        path, iter(_signatures()), read_type=ReadType.ONT, reference_id=REF
    )
    assert n == 3
    assert path.exists()
    assert not Path(f"{path}.partial").exists()
    meta, rows = read_signatures(path)
    assert meta == {"schema": "lr-signatures", "version": "1",
                    "read_type": "ont", "reference_id": REF}  # fmt: skip
    got = list(rows)
    assert got == _signatures()
    assert got[0].clip_side is ClipSide.NONE
    assert got[2].kind is SignatureKind.CLIP


def test_signature_file_layout(tmp_path):
    path = tmp_path / "sigs.tsv.gz"
    write_signatures(path, _signatures(), read_type=ReadType.HIFI,
                     reference_id=REF)  # fmt: skip
    with gzip.open(path, "rt") as fh:
        lines = fh.read().split("\n")
    assert lines[0] == (
        f"#relocate3:lr-signatures\tv1\tread_type=hifi\treference_id={REF}"
    )
    assert lines[1] == "\t".join(SIGNATURE_COLUMNS)
    assert lines[3].split("\t")[-1] == ""  # empty right_flank kept
    assert lines[-1] == ""  # file ends with a newline


def test_clusters_round_trip(tmp_path):
    path = tmp_path / "clusters.tsv.gz"
    n = write_clusters(path, _clusters(), read_type=ReadType.ONT,
                       reference_id=REF)  # fmt: skip
    assert n == 2
    assert not Path(f"{path}.partial").exists()
    meta, rows = read_clusters(path)
    assert meta["schema"] == "lr-clusters"
    assert meta["version"] == "1"
    got = list(rows)
    assert got == _clusters()
    assert got[1].insert_lengths == ()
    with gzip.open(path, "rt") as fh:
        assert fh.read().split("\n")[1] == "\t".join(CLUSTER_COLUMNS)


def test_existing_output_refused(tmp_path):
    path = tmp_path / "sigs.tsv.gz"
    path.write_text("keep me")
    with pytest.raises(FileExistsError):
        write_signatures(path, _signatures(), read_type=ReadType.ONT,
                         reference_id=REF)  # fmt: skip
    assert path.read_text() == "keep me"
    assert not Path(f"{path}.partial").exists()


def test_wrong_schema_rejected(tmp_path):
    path = tmp_path / "clusters.tsv.gz"
    write_clusters(path, _clusters(), read_type=ReadType.ONT,
                   reference_id=REF)  # fmt: skip
    with pytest.raises(ValueError, match="schema"):
        read_signatures(path)


def test_unknown_version_rejected(tmp_path):
    path = _write_raw(
        tmp_path / "v2.tsv.gz",
        [_schema_line("lr-signatures", "v2"), "\t".join(SIGNATURE_COLUMNS)],
    )
    with pytest.raises(ValueError, match="version"):
        read_signatures(path)


def test_altered_header_rejected(tmp_path):
    cols = list(CLUSTER_COLUMNS)
    cols[3] = "begin"
    path = _write_raw(
        tmp_path / "hdr.tsv.gz", [_schema_line("lr-clusters"), "\t".join(cols)]
    )
    with pytest.raises(ValueError, match="header"):
        read_clusters(path)


def test_short_row_rejected(tmp_path):
    row = "\t".join(["c1", REF, "Chr3", "1", "1", "+", "1"])
    path = _write_raw(
        tmp_path / "short.tsv.gz",
        [_schema_line("lr-clusters"), "\t".join(CLUSTER_COLUMNS), row],
    )
    _, rows = read_clusters(path)
    with pytest.raises(ValueError, match="fields"):
        list(rows)


def test_failed_write_leaves_nothing(tmp_path):
    path = tmp_path / "sigs.tsv.gz"

    def gen():
        yield _signatures()[0]
        raise RuntimeError("boom")

    with pytest.raises(RuntimeError, match="boom"):
        write_signatures(path, gen(), read_type=ReadType.ONT, reference_id=REF)
    assert not path.exists()
    assert not Path(f"{path}.partial").exists()


def test_tab_in_field_rejected(tmp_path):
    path = tmp_path / "sigs.tsv.gz"
    bad = Signature(
        "s\t1", SignatureKind.INS, Locus(REF, "Chr3", 1, 1, "+"), 5,
        ClipSide.NONE, "r1", 0, 5, 60, "ACGTA", "", "",
    )  # fmt: skip
    with pytest.raises(ValueError, match="tab or newline"):
        write_signatures(path, [bad], read_type=ReadType.ONT, reference_id=REF)
    assert not path.exists()
    assert not Path(f"{path}.partial").exists()


@pytest.mark.parametrize("sig_id", ["", "r,1:INS:0"])
def test_ambiguous_signature_id_rejected(tmp_path, sig_id):
    path = tmp_path / "sigs.tsv.gz"
    bad = Signature(
        sig_id, SignatureKind.INS, Locus(REF, "Chr3", 1, 1, "+"), 5,
        ClipSide.NONE, "r1", 0, 5, 60, "ACGTA", "", "",
    )  # fmt: skip
    with pytest.raises(ValueError, match="signature_id"):
        write_signatures(path, [bad], read_type=ReadType.ONT, reference_id=REF)
    assert not path.exists()


@pytest.mark.parametrize("ids", [("",), ("s1", "a,b")])
def test_ambiguous_cluster_signature_ids_rejected(tmp_path, ids):
    path = tmp_path / "clusters.tsv.gz"
    bad = Cluster("c1", Locus(REF, "Chr3", 1, 1, "+"), 1, 1, 0, 0, 0, (5,), ids)
    with pytest.raises(ValueError, match="signature_ids"):
        write_clusters(path, [bad], read_type=ReadType.ONT, reference_id=REF)
    assert not path.exists()


@pytest.mark.parametrize(
    ("col", "value"), [("start", "x"), ("kind", "DEL"), ("strand", ".")]
)
def test_bad_row_value_names_line(tmp_path, col, value):
    path = tmp_path / "sigs.tsv.gz"
    write_signatures(path, _signatures(), read_type=ReadType.ONT,
                     reference_id=REF)  # fmt: skip
    with gzip.open(path, "rt") as fh:
        lines = fh.read().split("\n")
    fields = lines[3].split("\t")
    fields[SIGNATURE_COLUMNS.index(col)] = value
    lines[3] = "\t".join(fields)
    bad = _write_raw(tmp_path / "bad.tsv.gz", lines[:-1])
    _, rows = read_signatures(bad)
    with pytest.raises(ValueError, match=r"bad\.tsv\.gz:4: "):
        list(rows)


def test_none_field_rejected(tmp_path):
    path = tmp_path / "sigs.tsv.gz"
    bad = Signature(
        "s1", SignatureKind.INS, Locus(REF, "Chr3", 1, 1, "+"), 5,
        ClipSide.NONE, "r1", 0, 5, 60, None, "", "",  # type: ignore[arg-type]
    )  # fmt: skip
    with pytest.raises(ValueError, match="None"):
        write_signatures(path, [bad], read_type=ReadType.ONT, reference_id=REF)
    assert not path.exists()
    assert not Path(f"{path}.partial").exists()
