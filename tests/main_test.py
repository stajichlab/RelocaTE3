from __future__ import annotations

import logging
from typing import Generator

import pytest

from RelocaTE3 import __main__, __version__
from RelocaTE3.__main__ import main


def mockreturn(**kwargs):
    for k, v in kwargs.items():
        print(f"{k}: {v}")


def test_main_version(capsys: pytest.CaptureFixture):
    assert main(["--version"]) == main(["-V"])
    captured: str = capsys.readouterr().out
    captured = captured.split("\n")
    assert captured[0] == captured[1] == __version__


def test_main_align(monkeypatch: Generator):
    monkeypatch.setattr(__main__, "align", mockreturn)
    assert main(["-i", "tests/data/example_input"]) == 0


def test_main_align_verbose(monkeypatch: Generator, caplog: pytest.LogCaptureFixture):
    monkeypatch.setattr(__main__, "align", mockreturn)
    with caplog.at_level(logging.DEBUG):
        main(["-i", "tests/data/example_input", "-v"])
    assert any(record.levelname == "DEBUG" for record in caplog.records)
    assert "Debug mode enabled." in caplog.text
