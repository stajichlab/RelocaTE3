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


def read_truth_gff(path: Path) -> list[TruthEvent]:
    """Read the fixture's truth GFF."""
    events = []
    with open(path) as handle:
        for number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(f"{path}:{number}: expected 9 columns")
            ident = re.search(r"ID=([^;]+)", fields[8])
            tsd = re.search(r"TSD=([ACGT]+)", fields[8])
            if ident is None or tsd is None:
                raise ValueError(f"{path}:{number}: missing ID or TSD attribute")
            if fields[6] not in ("+", "-"):
                raise ValueError(f"{path}:{number}: bad strand {fields[6]!r}")
            events.append(
                TruthEvent(
                    event_id=ident.group(1),
                    chrom=fields[0],
                    position=int(fields[3]),
                    family=fields[2],
                    strand=fields[6],
                    tsd=tsd.group(1),
                )
            )
    return events
