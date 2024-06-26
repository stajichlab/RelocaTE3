from __future__ import annotations

from pathlib import Path

import pytest

from RelocaTE3.pipeline import align


def test_align(tmp_path: Path):
    with open("tests/data/example_input") as f:
        assert align(f, tmp_path / "output.file") == [1, 4, 9]


def test_align_fileexisterror():
    with pytest.raises(FileExistsError, match="Output file already exists."):
        with open("tests/data/example_input") as f:
            align(f, "examples/example.py")


def test_align_valueerror(tmp_path: Path):
    with pytest.raises(ValueError, match="Input contains non-numeric strings."):
        with open("tests/data/example_input_invalid") as f:
            align(f, tmp_path / "output.file")
