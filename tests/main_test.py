from __future__ import annotations

import logging
import shutil
import subprocess
import sys
from typing import Generator

import pytest

from RelocaTE3 import __version__, cli
from RelocaTE3.cli import main


def mockreturn(args):
    """Stand-in command handler that records its Namespace and succeeds."""
    mockreturn.args = args
    return 0


def test_main_version(capsys: pytest.CaptureFixture):
    assert main(["--version"]) == main(["-V"])
    captured: str = capsys.readouterr().out
    captured = captured.split("\n")
    assert captured[0] == captured[1] == __version__


def test_main_map(monkeypatch: Generator):
    monkeypatch.setattr(cli, "cmd_map", mockreturn)
    assert main(["map", "-l", "r1.fq", "-T", "te.fa", "-n", "HEG4"]) == 0
    assert mockreturn.args.name == "HEG4"
    assert mockreturn.args.left == "r1.fq"
    assert mockreturn.args.te_library == "te.fa"


def test_main_map_verbose(monkeypatch: Generator, caplog: pytest.LogCaptureFixture):
    monkeypatch.setattr(cli, "cmd_map", mockreturn)
    with caplog.at_level(logging.DEBUG):
        main(["map", "-l", "r1.fq", "-T", "te.fa", "-n", "HEG4", "-v"])
    assert any(record.levelname == "DEBUG" for record in caplog.records)
    assert "Debug mode enabled." in caplog.text


def test_main_characterize(monkeypatch: Generator):
    monkeypatch.setattr(cli, "cmd_characterize", mockreturn)
    assert main(["characterize", "-s", "sites.txt", "-b", "a.bam", "b.bam"]) == 0
    assert mockreturn.args.sites_file == "sites.txt"
    assert mockreturn.args.bam == ["a.bam", "b.bam"]


def test_module_entry_point_version():
    """`python -m RelocaTE3 --version` runs and prints a version."""
    proc = subprocess.run(
        [sys.executable, "-m", "RelocaTE3", "--version"],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0
    assert proc.stdout.strip()


@pytest.mark.skipif(
    shutil.which("relocaTE3") is None,
    reason="console script not installed on PATH",
)
def test_console_script_version():
    """The installed `relocaTE3 --version` console script runs."""
    proc = subprocess.run(["relocaTE3", "--version"], capture_output=True, text=True)
    assert proc.returncode == 0
    assert proc.stdout.strip()


def test_build_parser_registers_subcommand():
    """build_parser() returns a parser whose subcommands dispatch to cmd_*."""
    from RelocaTE3.cli import build_parser, cmd_index_genome

    parsed = build_parser().parse_args(["index-genome", "-g", "ref.fa"])
    assert parsed.func is cmd_index_genome


def test_te_and_genome_opts_are_parsed_and_shlex_split():
    """--te-opts/--genome-opts reach the aligner backends as argument lists.

    The benchmark drives RelocaTE3 through this CLI, so the per-stage aligner
    option passthrough is only usable if these flags exist and split correctly.

    Values start with "-" (e.g. -minIdentity=80), which argparse would otherwise
    read as another flag, so the "--opt=value" form is required. Callers/adapters
    must use that form.
    """
    from RelocaTE3.cli import build_parser, split_aligner_opts

    parser = build_parser()

    te = parser.parse_args(
        ["run", "-l", "r1.fq", "-T", "te.fa", "-n", "s", "-o", "out",
         "--te-aligner", "blat", "--te-opts=-minIdentity=80"]
    )
    assert split_aligner_opts(te.te_opts) == ["-minIdentity=80"]

    # multi-token values must split on whitespace, not be passed as one arg
    mapped = parser.parse_args(
        ["map", "-l", "r1.fq", "-T", "te.fa", "-n", "s", "-o", "out",
         "--te-opts=-B 4"]
    )
    assert split_aligner_opts(mapped.te_opts) == ["-B", "4"]

    genome = parser.parse_args(
        ["align-genome", "-g", "ref.fa", "-f", "f.fq", "-n", "s", "-o", "out",
         "--genome-aligner", "bwaaln", "--genome-opts=-n 0.10"]
    )
    assert split_aligner_opts(genome.genome_opts) == ["-n", "0.10"]

    # absent flags default to no extra options
    assert split_aligner_opts(None) == []
