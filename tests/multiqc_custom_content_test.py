"""MultiQC reads ``relocaTE3 cohort``'s ``cohort.summary_mqc.tsv``.

The summary is MultiQC "custom content" (YAML config in ``#`` lines ahead of a
TSV). This runs the real MultiQC of the pixi env the tests run in (the
``pipeline`` env, the image's environment) on a tiny cohort and checks the
section and its table reach ``multiqc_data/multiqc_data.json``. Skipped when
that env has no ``multiqc`` (e.g. pytest from the ``blat`` env). A ``multiqc``
elsewhere on PATH (e.g. ~/.local/bin) is deliberately not used.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from RelocaTE3.cli import main

MULTIQC = Path(sys.executable).parent / "multiqc"

pytestmark = pytest.mark.skipif(
    not MULTIQC.exists(),
    reason="no multiqc next to the test's Python (run pytest from the pipeline env)",
)

HEADER = (
    "strain\tTE\tTSD\tchromosome.pos\tstrand\tavg_flankers\tspanners\tstatus\t"
    "TE_family_support"
)


def _table(path: Path, sample: str, statuses: list[str]) -> Path:
    rows = [
        f"{sample}\tmPing\tTAA\tChr1:{100 * i}..{100 * i + 2}\t+\t10\t0\t{st}\tmPing=5"
        for i, st in enumerate(statuses, start=1)
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join([HEADER, *rows]) + "\n")
    return path


def test_multiqc_renders_the_cohort_summary_table(tmp_path):
    a = _table(tmp_path / "in" / "A.txt", "A", ["homozygous", "heterozygous", "homozygous/excision_no_footprint"])
    b = _table(tmp_path / "in" / "B.txt", "B", ["somatic"])
    assert main(["cohort", "-i", str(a), str(b), "-o", str(tmp_path / "cohort")]) == 0

    # MultiQC sees only the summary (the call table and matrix are not _mqc).
    mqc_in = tmp_path / "mqc_in"
    mqc_in.mkdir()
    shutil.copy(tmp_path / "cohort" / "cohort.summary_mqc.tsv", mqc_in)
    out = tmp_path / "mqc_out"
    env = {**os.environ, "PYTHONNOUSERSITE": "1"}
    proc = subprocess.run(
        [str(MULTIQC), "--no-ansi", "-o", str(out), str(mqc_in)],
        capture_output=True, text=True, env=env, cwd=tmp_path,
    )  # fmt: skip
    log = proc.stdout + proc.stderr
    assert proc.returncode == 0, log
    assert "relocate3_cohort_summary" in log

    data = json.loads((out / "multiqc_data" / "multiqc_data.json").read_text())
    # pconfig id of the table (cohort.py _MQC_HEADER).
    assert "relocate3_cohort_summary_table" in data["report_plot_data"]
    table = data["report_saved_raw_data"]["multiqc_relocate3_cohort_summary_table"]
    assert set(table) == {"A", "B"}
    assert table["A"]["total_calls"] == 3
    assert table["A"]["homozygous"] == 2
    assert table["A"]["heterozygous"] == 1
    assert table["A"]["excision"] == 1
    assert table["B"]["total_calls"] == 1
    assert (out / "multiqc_report.html").is_file()
