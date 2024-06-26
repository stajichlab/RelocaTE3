"""Pipelines for relocaTE3."""

from __future__ import annotations

from io import TextIOWrapper
from pathlib import Path

from RelocaTE3 import logger


def align(input: TextIOWrapper, output: Path, *args, **kwargs) -> None:
    """Act as an example function that triggered from the main entry point."""
    logger.info(input)

    output = Path(output)
    if output.exists():
        raise FileExistsError("Output file already exists.")

    results = []
    for line in input.read().split("\n"):
        line = line.strip()
        if line:
            try:
                results.append(int(line) * int(line))
            except ValueError:
                raise ValueError("Input contains non-numeric strings.")

    print(results)
    return results
