"""Nextflow head jobs ask for realistic SLURM time limits.

A head job only runs Nextflow, which submits and waits for the task jobs. The
10 real-rice samples took 1 h 04 min (job 31830681) and the riceTElib 30x BLAT
scatter run a few hours, yet the scripts asked for 3 days. A long limit is not
free: SLURM holds a job whose limit overlaps a maintenance reservation (job
31830681 pended until its limit was cut to 1 day). Default to at most 1 day;
a bigger run overrides with ``sbatch --time=...``, and ``-resume`` restarts a
head job that hit its limit.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
HEAD_JOBS = [
    "validation/real_rice/run_pipeline.sh",
    "validation/blat_scatter/run_scatter_pipeline.sh",
    "docs/source/pipeline.rst",  # the guide's sbatch template
    "conf/hpcc.config",  # the launch pattern in its header comment
]
# User-facing docs that may carry sbatch templates (checked when they do).
DOC_TEMPLATES = [
    "README.md",
    "docs/source/pipeline.rst",
    "docs/source/installation.rst",
    "conf/hpcc.config",
    "validation/real_rice/README.md",
]
SBATCH_TIME = re.compile(r"#SBATCH[^\n]*?--time[= ](\S+)")
ONE_DAY = 24 * 3600


def _seconds(limit: str) -> int:
    """SLURM time formats: [D-]HH:MM:SS, MM:SS, MM."""
    days = 0
    if "-" in limit:
        d, limit = limit.split("-", 1)
        days = int(d)
    parts = [int(p) for p in limit.split(":")]
    if days:
        parts += [0] * (3 - len(parts))
        h, m, s = parts
    elif len(parts) == 3:
        h, m, s = parts
    elif len(parts) == 2:
        h, m, s = 0, *parts
    else:
        h, m, s = 0, parts[0], 0
    return ((days * 24 + h) * 60 + m) * 60 + s


@pytest.mark.parametrize("path", HEAD_JOBS)
def test_head_job_time_limit_at_most_one_day(path):
    limits = SBATCH_TIME.findall((ROOT / path).read_text())
    assert limits, f"{path}: no #SBATCH --time"
    for limit in limits:
        assert _seconds(limit) <= ONE_DAY, f"{path}: --time={limit} is more than a day"


@pytest.mark.parametrize("path", HEAD_JOBS)
def test_head_job_documents_the_override(path):
    assert "sbatch --time=" in (ROOT / path).read_text()


@pytest.mark.parametrize("path", DOC_TEMPLATES)
def test_doc_templates_time_limit_at_most_one_day(path):
    for limit in SBATCH_TIME.findall((ROOT / path).read_text()):
        assert _seconds(limit) <= ONE_DAY, f"{path}: --time={limit} is more than a day"


@pytest.mark.parametrize("path", DOC_TEMPLATES)
def test_doc_templates_create_the_log_dir_before_sbatch(path):
    """SLURM opens the #SBATCH -o log before the script runs: with a missing
    logs/ the job fails without a log, so a mkdir inside the script is too
    late. A template that logs under logs/ says to create it at submission."""
    text = (ROOT / path).read_text()
    if re.search(r"#SBATCH[^\n]*?(-o|--output)[= ]logs/", text):
        assert "mkdir -p logs && sbatch" in text, f"{path}: submit with mkdir -p logs && sbatch ..."


def test_seconds():
    assert _seconds("1-00:00:00") == ONE_DAY
    assert _seconds("3-00:00:00") == 3 * ONE_DAY
    assert _seconds("02:00:00") == 7200
    assert _seconds("1-12") == 36 * 3600
