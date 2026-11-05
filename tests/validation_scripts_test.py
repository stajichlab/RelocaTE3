"""Shell scripts under validation/ that M3 added or changed.

* ``validation/real_rice/cohort_report.sh``: cohort + MultiQC report over a
  finished pipeline output directory; run here on a two-sample fake tree.
* ``validation/release/fresh_clone_check.slurm``: the release check's full
  pytest run must use the ``pipeline`` env (blat + MultiQC), or the MultiQC
  tests are skipped; the blat-only checks stay on ``blat``.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from tests.cohort_test import HEADER, _row
from tests.sample_stats_test import _sample

ROOT = Path(__file__).resolve().parent.parent
COHORT_REPORT = ROOT / "validation" / "real_rice" / "cohort_report.sh"
RELEASE_CHECK = ROOT / "validation" / "release" / "fresh_clone_check.slurm"
PIPELINE_BIN = ROOT / ".pixi" / "envs" / "pipeline" / "bin"


# ------------------------------------------------------------ cohort_report.sh


def test_cohort_report_keeps_temp_files_in_the_outdir_and_reads_lines_safely():
    text = COHORT_REPORT.read_text()
    mktemps = re.findall(r"mktemp[^)]*", text)
    assert mktemps and all('-p "${OUT}"' in m for m in mktemps), mktemps
    assert "while IFS= read -r" in text
    assert not re.search(r"for \w+ in \$\(", text), "unquoted word-split loop"


@pytest.mark.skipif(
    not (PIPELINE_BIN / "multiqc").is_file() or shutil.which("bash") is None,
    reason="needs the pipeline env (pixi install -e pipeline)",
)
def test_cohort_report_passes_on_a_two_sample_tree(tmp_path):
    tree = tmp_path / "pipeline_out"
    tables = {}
    for name, n in (("G_1", 2), ("B_10", 1)):
        d = _sample(tree, name=name, te="mping")
        rows = [_row(name, "Chr1", 10 * i + 1, 10 * i + 3) for i in range(n)]
        p = d / "results" / "ALL.mping.all_nonref_insert.characTErized.txt"
        p.write_text("".join(line + "\n" for line in [HEADER, *rows]))
        tables[name] = rows
    out = tmp_path / "report"
    env = {**os.environ, "PIXI_ENV": str(PIPELINE_BIN.parent)}
    proc = subprocess.run(
        ["bash", str(COHORT_REPORT), str(tree), str(out)],
        cwd=ROOT, env=env, capture_output=True, text=True,
    )  # fmt: skip
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "RESULT: PASS" in proc.stdout
    cohort = (out / "cohort" / "cohort.characTErized.txt").read_text().splitlines()
    assert cohort == [HEADER, *tables["B_10"], *tables["G_1"]]
    # the temporary file lived in the outdir and is gone
    assert sorted(p.name for p in out.iterdir()) == ["cohort", "mqc_in", "multiqc"]


# ------------------------------------------------------- fresh_clone_check.slurm


def test_release_check_runs_pytest_in_the_pipeline_env():
    text = RELEASE_CHECK.read_text()
    assert re.search(r"pixi install -e pipeline --locked", text)
    pytest_lines = [
        line for line in text.splitlines()
        if "python -m pytest" in line and "pixi run" in line
    ]  # fmt: skip
    assert pytest_lines and all("-e pipeline" in line for line in pytest_lines), pytest_lines
    # the blat checks and the built-wheel stage still use the blat env
    assert 'ENV_BIN="$CLONE/.pixi/envs/blat/bin"' in text
    assert "pixi install -e blat --locked" in text
