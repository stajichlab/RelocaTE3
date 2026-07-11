"""Compatibility shim for the RelocaTE3 command-line interface.

The canonical CLI now lives in :mod:`RelocaTE3.__main__`. This module is kept
so that ``from RelocaTE3.cli import main`` continues to work and reaches the
canonical parser.

Compatibility scope: the ``RelocaTE3.cli.main`` entry point is preserved. The
old ``cli.py`` subcommands and flags (e.g. ``--te``, ``--sample``,
``-i/--insertions``, the function-based ``find-reference``, and the
full-pipeline ``run``) are intentionally removed in favor of the single,
validated ``__main__`` interface.
"""

from __future__ import annotations

from RelocaTE3.__main__ import main

__all__ = ["main"]
