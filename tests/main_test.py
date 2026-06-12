"""Tests for the RelocaTE3 command-line interface."""

from __future__ import annotations

import logging

import pytest

from RelocaTE3 import __version__
from RelocaTE3.cli import main


def test_main_version(capsys: pytest.CaptureFixture):
    """--version prints the version and exits 0."""
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    captured = capsys.readouterr().out.strip()
    assert captured == __version__


def test_main_no_command_prints_help(capsys: pytest.CaptureFixture):
    """Running with no subcommand prints help and returns 0."""
    assert main([]) == 0
    err = capsys.readouterr().err
    assert "usage" in err.lower()


def test_main_verbose_enables_debug(caplog: pytest.LogCaptureFixture):
    """The global -v flag turns on debug logging."""
    with caplog.at_level(logging.DEBUG):
        main(["-v"])
    assert "Debug mode enabled." in caplog.text
