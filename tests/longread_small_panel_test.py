"""Shape of the built small long-read panel (skipped when it is absent)."""

from __future__ import annotations

import csv
import os
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PANEL = Path(
    os.environ.get(
        "RELOCATE3_LR_SMALL_PANEL", ROOT / "validation" / "longread_small" / "panel_v1"
    )
)

pytestmark = pytest.mark.skipif(
    not (PANEL / ".complete").is_file(),
    reason="small long-read panel not built (sbatch validation/longread_small/build_panel.slurm)",
)


def test_manifest_lists_four_samples_with_reads():
    with open(PANEL / "panel_manifest.tsv", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert sorted((r["platform"], r["coverage"]) for r in rows) == [
        ("hifi", "20"),
        ("hifi", "5"),
        ("ont-hq", "20"),
        ("ont-hq", "5"),
    ]
    for row in rows:
        for column in ("reads", "control_reads"):
            path = PANEL / row[column]
            assert path.is_file() and path.stat().st_size > 0, path


def test_truth_has_200_homozygous_events():
    with open(PANEL / "truth_events.tsv", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 200
    assert {r["biological_class"] for r in rows} == {"homozygous"}
