"""The five long-read stages end to end on the small panel's hifi 5x sample.

Runs lr-align, lr-signatures, lr-cluster, lr-annotate and lr-call on the
``hifi`` ``cov5x_rep1`` reads of the small long-read panel (200 homozygous
insertions on the 2 Mb Chr3 genome). Skipped when the panel is absent or
minimap2/samtools are not on PATH.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

import pytest

from RelocaTE3.cli import main
from RelocaTE3.longread.model import CallTier
from RelocaTE3.longread.tables import read_calls
from tests.longread_align_test import HAVE_TOOLS

ROOT = Path(__file__).resolve().parents[1]
# An empty RELOCATE3_LR_SMALL_PANEL falls back to the default, not to ".".
PANEL = Path(
    os.environ.get("RELOCATE3_LR_SMALL_PANEL")
    or ROOT / "tests" / "data" / "longread_small"
)
GENOME = ROOT / "tests" / "data" / "sim_genome" / "MSU7.Chr3_2M.fa"
STAGES = ["lr-align", "lr-signatures", "lr-cluster", "lr-annotate", "lr-call"]

pytestmark = [
    pytest.mark.skipif(
        not (PANEL / ".complete").is_file(),
        reason=(
            f"small long-read panel not found at {PANEL}"
            " (ln -s <make_simulation_new>/results/chr3_lr_small"
            " tests/data/longread_small, or set RELOCATE3_LR_SMALL_PANEL)"
        ),
    ),
    pytest.mark.skipif(not HAVE_TOOLS, reason="minimap2/samtools not on PATH"),
]


def _reads(platform: str, sample: str) -> Path:
    with open(PANEL / "panel_manifest.tsv", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    (row,) = [r for r in rows if r["platform"] == platform and r["sample"] == sample]
    return PANEL / row["reads"]


def _library_names(path: Path) -> set[str]:
    with path.open() as handle:
        return {line[1:].split()[0] for line in handle if line.startswith(">")}


def test_hifi_cov5x_five_stages_call_the_panel(tmp_path: Path) -> None:
    library = PANEL / "config" / "te_library.fa"
    names = _library_names(library)
    assert names, f"no FASTA headers in {library}"
    out = tmp_path / "out"
    so = ["--sample", "S", "-o", str(out)]
    sigs = out / "S.signatures.tsv.gz"
    assert main(["lr-align", *so, "--genome", str(GENOME), "--read-type", "hifi",
                 "--reads", str(_reads("hifi", "cov5x_rep1")),
                 "--threads", "2"]) == 0  # fmt: skip
    assert main(["lr-signatures", *so, "--bam", str(out / "S.lr.bam"),
                 "--read-type", "hifi"]) == 0  # fmt: skip
    assert main(["lr-cluster", *so, "--signatures", str(sigs)]) == 0
    assert main(["lr-annotate", *so, "--signatures", str(sigs),
                 "--te-library", str(library), "--threads", "2"]) == 0  # fmt: skip
    assert main(["lr-call", *so, "--clusters", str(out / "S.clusters.tsv.gz"),
                 "--signatures", str(sigs),
                 "--te-hits", str(out / "S.te_hits.tsv.gz")]) == 0  # fmt: skip

    calls_path = out / "S.calls.tsv.gz"
    assert calls_path.is_file()
    _, calls = read_calls(calls_path)
    spanning = [c for c in calls if c.passed and c.tier is CallTier.SPANNING]
    assert len(spanning) >= 150
    named = sum(c.family in names for c in spanning)
    assert named >= 0.9 * len(spanning)

    with (out / "timing.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [r["stage"] for r in rows] == STAGES
    assert all(r["status"] == "ok" for r in rows)
