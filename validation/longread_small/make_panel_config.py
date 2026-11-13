#!/usr/bin/env python3
"""Build the te-panel config for the small Chr3 long-read panel.

Input: the Chr3 2 Mb fixture (reference FASTA, truth GFF, TE library).
Output: a directory holding everything `simulate-data te-panel` needs.
Standard library only, so it runs in any Python >= 3.10.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class TruthEvent:
    """One simulated insertion, in reference coordinates."""

    event_id: str
    chrom: str
    position: int  # GFF column 4; the TSD is reference[position-L:position]
    family: str
    strand: str
    tsd: str


def _single_attribute(attrs: str, key: str, where: str) -> str | None:
    """Return the one value of ``key`` in a GFF attribute string, or None."""
    values = re.findall(rf"(?:^|;)\s*{key}=([^;]*)", attrs)
    if len(values) > 1:
        raise ValueError(f"{where}: {key} attribute given {len(values)} times")
    return values[0] if values else None


def read_truth_gff(path: Path) -> list[TruthEvent]:
    """Read the fixture's truth GFF.

    Args:
        path: GFF with 9 tab-separated columns per row; column 3 is the TE
            family, column 4 the 1-based position, column 7 the strand, and
            column 9 must carry exactly one ``ID=`` and one ``TSD=``.

    Returns:
        One TruthEvent per data row, in file order.

    Raises:
        ValueError: On a row with the wrong column count, a missing, blank or
            repeated ID, a missing, repeated or non-ACGT TSD, a strand other
            than ``+``/``-``, a position that is not an integer >= 1, or an
            ID already used by an earlier row.
    """
    events = []
    seen: set[str] = set()
    with open(path) as handle:
        for number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            where = f"{path}:{number}"
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{where}: expected 9 columns, got {len(fields)}")
            attrs = fields[8]
            ident = _single_attribute(attrs, "ID", where)
            ident = ident.strip() if ident is not None else ""
            if not ident:
                raise ValueError(f"{where}: missing ID attribute")
            tsd = _single_attribute(attrs, "TSD", where)
            if tsd is None:
                raise ValueError(f"{where}: missing TSD attribute")
            if not re.fullmatch(r"[ACGT]+", tsd):
                raise ValueError(f"{where}: bad TSD {tsd!r}")
            if fields[6] not in ("+", "-"):
                raise ValueError(f"{where}: bad strand {fields[6]!r}")
            try:
                position = int(fields[3])
            except ValueError:
                position = 0
            if position < 1:
                raise ValueError(f"{where}: bad position {fields[3]!r}")
            if ident in seen:
                raise ValueError(f"{where}: duplicate ID {ident!r}")
            seen.add(ident)
            events.append(
                TruthEvent(
                    event_id=ident,
                    chrom=fields[0],
                    position=position,
                    family=fields[2],
                    strand=fields[6],
                    tsd=tsd,
                )
            )
    return events
