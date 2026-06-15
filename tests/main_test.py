from __future__ import annotations

import logging
from typing import Generator

import pytest

from RelocaTE3 import __main__, __version__
from RelocaTE3.__main__ import main


def mockreturn(**kwargs):
    """Stand-in command handler that records its arguments and succeeds."""
    mockreturn.kwargs = kwargs
    return 0


def test_main_version(capsys: pytest.CaptureFixture):
    assert main(["--version"]) == main(["-V"])
    captured: str = capsys.readouterr().out
    captured = captured.split("\n")
    assert captured[0] == captured[1] == __version__


def test_main_map(monkeypatch: Generator):
    monkeypatch.setattr(__main__, "cmd_map", mockreturn)
    assert main(["map", "-l", "r1.fq", "-T", "te.fa", "-n", "HEG4"]) == 0
    assert mockreturn.kwargs["name"] == "HEG4"
    assert mockreturn.kwargs["left"] == "r1.fq"
    assert mockreturn.kwargs["te_library"] == "te.fa"


def test_main_map_verbose(monkeypatch: Generator, caplog: pytest.LogCaptureFixture):
    monkeypatch.setattr(__main__, "cmd_map", mockreturn)
    with caplog.at_level(logging.DEBUG):
        main(["map", "-l", "r1.fq", "-T", "te.fa", "-n", "HEG4", "-v"])
    assert any(record.levelname == "DEBUG" for record in caplog.records)
    assert "Debug mode enabled." in caplog.text


def test_main_characterize(monkeypatch: Generator):
    monkeypatch.setattr(__main__, "cmd_characterize", mockreturn)
    assert main(["characterize", "-s", "sites.txt", "-b", "a.bam", "b.bam"]) == 0
    assert mockreturn.kwargs["sites_file"] == "sites.txt"
    assert mockreturn.kwargs["bam"] == ["a.bam", "b.bam"]
