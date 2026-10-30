"""``--tsd`` accepts only RelocaTE2's ``UNK`` mode.

RelocaTE2 2.0.1 has no ``--tsd`` option: its driver always writes ``UNK`` into
regex.txt (relocaTE2.py:346), and its fixed-pattern branch cannot run --
``TSD_check`` takes 14 parameters but is called with 13 or 10
(relocaTE_insertionFinder.py:923, 941, 1603, 1611). RelocaTE3's motif
(``TTA``) and wildcard (``...``) modes therefore had no reference to match, so
they were removed; any other value is rejected with an explanation.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from RelocaTE3.cli import build_parser
from RelocaTE3.insertions import InsertionFinder

RUN_ALL = ["run-all", "-1", "r1.fq", "-T", "te.fa", "-g", "g.fa", "-n", "S"]
FIND = ["find-insertions", "-b", "x.bam", "-R", "rr.txt", "-n", "S"]


@pytest.mark.parametrize("base", [RUN_ALL, FIND])
@pytest.mark.parametrize("value", ["UNK", "unk", "UKN", "unknown"])
def test_cli_accepts_unknown_tsd_spellings(base, value):
    assert build_parser().parse_args(base + ["--tsd", value]).tsd == value


@pytest.mark.parametrize("base", [RUN_ALL, FIND])
@pytest.mark.parametrize("value", ["TTA", "...", "TAA"])
def test_cli_rejects_fixed_tsd(base, value, capsys):
    with pytest.raises(SystemExit) as exc:
        build_parser().parse_args(base + ["--tsd", value])
    assert exc.value.code == 2
    err = capsys.readouterr().err
    assert "UNK" in err and "RelocaTE2" in err


@pytest.mark.parametrize("value", ["TTA", "..."])
def test_finder_rejects_fixed_tsd(tmp_path, value):
    with pytest.raises(ValueError, match="UNK"):
        InsertionFinder().find_insertions(
            bam_file=tmp_path / "x.bam",
            read_repeat_file=tmp_path / "rr.txt",
            tsd=value,
            target="ALL",
            sample="S",
            outdir=tmp_path,
        )


@pytest.mark.parametrize(
    "name",
    [
        "_cluster_reads",
        "_assign_cluster",
        "_tsd_check",
        "_merge_offset_starts",
        "_write_output",
        "_write_event_start",
    ],
)
def test_fixed_pattern_path_is_gone(name):
    assert not hasattr(InsertionFinder, name)
