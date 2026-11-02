"""The README's RelocaTE2-defaults notes must agree with the shipped CLI.

Three notes went stale (fixed 2026-11-02): ``--size`` was described as not
implemented, ``--min-mapq`` as defaulting to 1, and ``--require-both-junctions``
as a deliberate divergence from RelocaTE2 (it is RelocaTE2's own
``characterizer.pl:91`` gate).
"""

from __future__ import annotations

import re
from pathlib import Path

from RelocaTE3.cli import build_parser

README = Path(__file__).resolve().parents[1] / "README.md"
FIND = ["find-insertions", "-b", "x.bam", "-R", "rr.txt", "-n", "S"]
RUN_ALL = ["run-all", "-1", "r1.fq", "-T", "te.fa", "-g", "g.fa", "-n", "S"]


def _section(title: str) -> str:
    text = README.read_text()
    start = text.index(f"## {title}")
    end = text.find("\n## ", start + 3)
    return text[start : end if end != -1 else None]


def test_size_note_matches_cli_default():
    notes = _section("RelocaTE2 defaults")
    for argv in (FIND, RUN_ALL):
        size = vars(build_parser().parse_args(argv))["insert_size"]
        assert f"`-s/--size` (default {size})" in notes
    assert "not implemented" not in notes


def test_min_mapq_note_matches_cli_default():
    notes = _section("RelocaTE2 defaults")
    for argv in (FIND, RUN_ALL):
        mapq = vars(build_parser().parse_args(argv))["min_mapq"]
        assert f"`--min-mapq` (default {mapq})" in notes


def test_two_sided_gate_is_not_called_a_divergence():
    text = _section("RelocaTE2 defaults") + _section(
        "One-sided calls and `--require-both-junctions`"
    )
    assert not re.search(r"deliberate(ly)? diverg", text)
    assert "characterizer.pl:91" in text
