"""index-genome builds the chosen genome aligner's index, and run-all always delegates.

In the Nextflow pipeline parallel samples must not each build the genome
aligner's index on demand inside align-genome, so index-genome has to build it
once up front.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from RelocaTE3 import cli
from RelocaTE3.cli import build_parser, main

from tests.main_test import _fake_stage_outputs


def _genome(tmp_path) -> Path:
    """A small random two-contig genome FASTA."""
    rng = random.Random(1)
    fa = tmp_path / "g.fa"
    fa.write_text(
        "".join(
            f">c{i}\n{''.join(rng.choice('ACGT') for _ in range(3000))}\n"
            for i in range(2)
        )
    )
    return fa


def test_index_genome_default_aligner_is_bwaaln():
    ns = build_parser().parse_args(["index-genome", "-g", "g.fa"])
    assert ns.genome_aligner == "bwaaln"


def test_index_genome_rejects_blat(capsys):
    with pytest.raises(SystemExit):
        build_parser().parse_args(
            ["index-genome", "-g", "g.fa", "--genome-aligner", "blat"]
        )
    # rejected as a bad choice, not merely as an unknown option
    assert "invalid choice: 'blat'" in capsys.readouterr().err


def test_index_genome_builds_the_genome_aligner_index(tmp_path):
    fa = _genome(tmp_path)
    assert main(["index-genome", "-g", str(fa), "--genome-aligner", "bwaaln"]) == 0
    for ext in (".fai", ".mmi", ".bwt"):
        assert Path(f"{fa}{ext}").is_file(), ext


def test_run_all_indexes_for_its_genome_aligner(monkeypatch, tmp_path):
    """Even when .fai/.mmi exist, run-all must make sure the genome-aligner index exists."""
    calls: list[tuple[str, object]] = []

    def rec(tag):
        def _inner(args):
            calls.append((tag, args))
            return 0

        return _inner

    for name, tag in [
        ("cmd_index_genome", "index-genome"),
        ("cmd_run", "run"),
        ("cmd_align_genome", "align-genome"),
        ("cmd_find_insertions", "find-insertions"),
        ("cmd_find_reference", "find-reference"),
        ("cmd_characterize", "characterize"),
    ]:
        monkeypatch.setattr(cli, name, rec(tag))

    ref = tmp_path / "ref.fa"
    ref.write_text(">c\nACGT\n")
    Path(f"{ref}.fai").write_text("")
    Path(f"{ref}.mmi").write_bytes(b"")
    _fake_stage_outputs(tmp_path, "S")

    rc = main(
        [
            "run-all",
            "-1",
            "r1.fq",
            "-T",
            "te.fa",
            "-g",
            str(ref),
            "-n",
            "S",
            "-o",
            str(tmp_path),
            "--genome-aligner",
            "bwaaln",
        ]
    )
    assert rc == 0
    tag, args = calls[0]
    assert tag == "index-genome"
    assert args.genome_aligner == "bwaaln"
    assert args.genome_fasta == str(ref)
