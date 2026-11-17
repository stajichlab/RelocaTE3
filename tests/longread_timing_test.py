"""Long-read per-stage timing table."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

import RelocaTE3
from RelocaTE3.longread.timing import TIMING_COLUMNS, StageCounts, stage_timer

HEADER = (
    "stage", "sample", "start_utc", "wall_s", "max_rss_mb", "records_in",
    "records_out", "relocate3_version", "status",
)  # fmt: skip


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_columns_in_order() -> None:
    assert TIMING_COLUMNS == HEADER


def test_two_stages_append_two_rows_under_one_header(tmp_path: Path) -> None:
    with stage_timer(tmp_path, "lr-align", "s1"):
        pass
    with stage_timer(tmp_path, "lr-signatures", "s1"):
        pass
    lines = (tmp_path / "timing.tsv").read_text().splitlines()
    assert lines[0].split("\t") == list(HEADER)
    assert len(lines) == 3
    rows = _rows(tmp_path / "timing.tsv")
    assert [r["stage"] for r in rows] == ["lr-align", "lr-signatures"]
    assert all(r["status"] == "ok" for r in rows)
    assert all(r["sample"] == "s1" for r in rows)


def test_exception_writes_status_and_propagates(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError, match="boom"):
        with stage_timer(tmp_path, "lr-cluster", "s1"):
            raise RuntimeError("boom")
    (row,) = _rows(tmp_path / "timing.tsv")
    assert row["status"] == "RuntimeError"
    assert row["stage"] == "lr-cluster"


def test_numeric_fields_and_version(tmp_path: Path) -> None:
    with stage_timer(tmp_path, "lr-align", "s1"):
        pass
    (row,) = _rows(tmp_path / "timing.tsv")
    assert float(row["wall_s"]) >= 0.0
    assert len(row["wall_s"].split(".")[1]) == 3
    assert float(row["max_rss_mb"]) > 0.0
    assert len(row["max_rss_mb"].split(".")[1]) == 1
    assert row["relocate3_version"] == RelocaTE3.__version__
    assert row["start_utc"].endswith("+00:00")
    assert "." not in row["start_utc"]


def test_counts_set_in_block_are_written(tmp_path: Path) -> None:
    with stage_timer(tmp_path, "lr-signatures", "s1") as counts:
        assert isinstance(counts, StageCounts)
        assert (counts.records_in, counts.records_out) == (0, 0)
        counts.records_in = 120
        counts.records_out = 7
    (row,) = _rows(tmp_path / "timing.tsv")
    assert (row["records_in"], row["records_out"]) == ("120", "7")


def test_missing_outdir_is_created(tmp_path: Path) -> None:
    outdir = tmp_path / "a" / "b"
    with stage_timer(outdir, "lr-align", "s1"):
        pass
    assert len(_rows(outdir / "timing.tsv")) == 1


def test_header_not_repeated_for_existing_file(tmp_path: Path) -> None:
    with stage_timer(tmp_path, "lr-align", "s1"):
        pass
    with stage_timer(tmp_path, "lr-align", "s2"):
        pass
    text = (tmp_path / "timing.tsv").read_text()
    assert text.count("stage\tsample\t") == 1
    assert [r["sample"] for r in _rows(tmp_path / "timing.tsv")] == ["s1", "s2"]
