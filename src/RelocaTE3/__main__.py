"""Module entry point so ``python -m RelocaTE3`` runs the canonical CLI."""

from __future__ import annotations

from RelocaTE3.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
