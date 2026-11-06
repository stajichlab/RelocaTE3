"""Boolean pipeline settings given on the real command line.

Nextflow >= 26.04 parses scripts with the strict syntax parser, which turns off
the old type detection of command-line parameters: an untyped
``--genotype false`` arrives as the String ``"false"``, which Groovy treats as
true, so ``if (params.genotype)`` aligned the reads anyway (30x job 31766412)
and ``--require_both_junctions false`` would still have passed
``--require-both-junctions``. The typed ``params {}`` block in main.nf converts
them to Boolean.

The nf-tests (tests/main.nf.test) set params in ``params {}`` blocks, which
nf-test hands to Nextflow already typed, so they never saw a command-line
string; these tests run ``nextflow run`` with ``--<param> false`` / ``true`` /
bare ``--<param>`` from a temporary launch directory instead.

Task selection is checked in stub mode (seconds). Flags are checked in the
``.command.sh`` of a real run (stub mode puts the ``stub:`` block there).
Needs the pixi ``nf`` and ``pipeline`` environments; skipped without.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
READS = ROOT / "tests" / "data" / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500"

pytestmark = pytest.mark.skipif(
    shutil.which("pixi") is None or not (ROOT / ".pixi" / "envs" / "nf").is_dir(),
    reason="needs pixi and the 'nf' environment (pixi install -e nf)",
)

NEXTFLOW = ["pixi", "run", "--manifest-path", str(ROOT / "pixi.toml"), "-e", "nf", "nextflow"]


def _sheet(launch: Path) -> Path:
    sheet = launch / "samples.csv"
    sheet.write_text(f"sample,fastq_1,fastq_2\nchr3,{READS}_1.fq.gz,{READS}_2.fq.gz\n")
    return sheet


def _run(launch: Path, *extra: str, stub: bool = True) -> subprocess.CompletedProcess:
    cmd = [*NEXTFLOW, "run", str(ROOT), "-profile", "test,pixi",
           "--input", str(_sheet(launch)), "--outdir", str(launch / "out")]
    if stub:
        cmd.append("-stub")
    cmd.extend(extra)
    return subprocess.run(cmd, cwd=launch, capture_output=True, text=True, timeout=1800)


def _tasks(launch: Path) -> set[str]:
    """Simple names of the processes that ran, from the run's trace."""
    trace = launch / "out" / "pipeline_info" / "execution_trace.txt"
    with trace.open() as fh:
        return {row["name"].split(" ")[0].split(":")[-1]
                for row in csv.DictReader(fh, delimiter="\t")}


def _command(launch: Path, subcommand: str) -> str:
    """The one .command.sh that runs `relocaTE3 <subcommand>`, words joined."""
    found = [p.read_text() for p in (launch / "work").glob("*/*/.command.sh")
             if f"relocaTE3 {subcommand} " in p.read_text()]
    assert len(found) == 1, f"{subcommand}: {len(found)} scripts"
    return " " + " ".join(found[0].replace("\\\n", " ").split()) + " "


# --------------------------------------------------- task selection (stub)

# Every spelling of each value: false, true, and a bare flag (= true).
SPELLINGS = {
    "false": ["--genotype", "false", "--skip_multiqc", "false", "--blat_scatter", "false"],
    "true": ["--genotype", "true", "--skip_multiqc", "true", "--blat_scatter", "true"],
    "bare": ["--genotype", "--skip_multiqc", "--blat_scatter"],
}


@pytest.mark.parametrize("spelling", sorted(SPELLINGS))
def test_boolean_settings_choose_the_tasks(tmp_path, spelling):
    proc = _run(tmp_path, *SPELLINGS[spelling])
    assert proc.returncode == 0, proc.stdout + proc.stderr
    tasks = _tasks(tmp_path)
    on = spelling != "false"
    # genotype: the reads are aligned for genotyping and the sample characterized
    for proc_name in ("RELOCATE3_ALIGN_READS", "RELOCATE3_CHARACTERIZE", "RELOCATE3_COHORT"):
        assert (proc_name in tasks) == on, (spelling, proc_name, tasks)
    # skip_multiqc: no report, no per-sample counts
    for proc_name in ("MULTIQC", "RELOCATE3_SAMPLE_STATS"):
        assert (proc_name in tasks) == (not on), (spelling, proc_name, tasks)
    # blat_scatter: split -> chunks -> merge, instead of one TE_SEARCH task
    assert ("RELOCATE3_TE_SPLIT" in tasks) == on, (spelling, tasks)
    assert ("RELOCATE3_TE_SEARCH" in tasks) == (not on), (spelling, tasks)


def test_params_file_false(tmp_path):
    params = tmp_path / "params.yaml"
    params.write_text("genotype: false\nblat_scatter: false\nskip_multiqc: true\n")
    proc = _run(tmp_path, "-params-file", str(params))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    tasks = _tasks(tmp_path)
    assert "RELOCATE3_ALIGN_READS" not in tasks
    assert "RELOCATE3_TE_SEARCH" in tasks and "RELOCATE3_TE_SPLIT" not in tasks
    assert "MULTIQC" not in tasks


def test_boolean_setting_rejects_other_words(tmp_path):
    proc = _run(tmp_path, "--genotype", "maybe")
    log = proc.stdout + proc.stderr
    assert proc.returncode != 0
    assert "genotype" in log
    assert not (tmp_path / "out" / "pipeline_info" / "execution_trace.txt").exists() \
        or not _tasks(tmp_path)


def test_help_false_runs_the_pipeline(tmp_path):
    proc = _run(tmp_path, "--help", "false", "-preview")
    log = proc.stdout + proc.stderr
    assert proc.returncode == 0, log
    # The parameter summary, not the help text.
    assert "Core Nextflow options" in log
    assert "-profile <pixi|docker" not in log


# -------------------------------------------------- flags (real tiny run)


@pytest.mark.parametrize("args, junctions, excision", [
    (["--require_both_junctions", "false", "--excision", "true"],
     "--no-require-both-junctions", True),
    (["--require_both_junctions", "--excision", "false"],
     "--require-both-junctions", False),
], ids=["false-true", "bare-false"])
def test_boolean_settings_reach_their_flags(tmp_path, args, junctions, excision):
    # Real run (no -stub) without the report, TE search in one task.
    proc = _run(tmp_path, *args, "--skip_multiqc", "true", "--blat_scatter", "false", stub=False)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    find = _command(tmp_path, "find-insertions")
    assert f" {junctions} " in find, find
    other = ({"--require-both-junctions", "--no-require-both-junctions"} - {junctions}).pop()
    assert f" {other} " not in find, find
    characterize = _command(tmp_path, "characterize")
    assert (" --excision " in characterize) == excision, characterize
