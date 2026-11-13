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


# First base of each IUPAC code, so the choice is deterministic and recorded.
AMBIGUITY = {
    "R": "A",
    "Y": "C",
    "S": "C",
    "W": "A",
    "K": "G",
    "M": "A",
    "B": "C",
    "D": "A",
    "H": "A",
    "V": "A",
}


def read_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA file.

    Args:
        path: FASTA file; each record is named by the first word of its header.

    Returns:
        ``{record name: uppercase sequence}``, in file order.

    Raises:
        ValueError: On a header with no name, a repeated record name, sequence
            before the first header, or a file with no records.
    """
    records: dict[str, list[str]] = {}
    name = None
    with open(path) as handle:
        for number, line in enumerate(handle, start=1):
            where = f"{path}:{number}"
            if line.startswith(">"):
                words = line[1:].split()
                if not words:
                    raise ValueError(f"{where}: header without a name")
                name = words[0]
                if name in records:
                    raise ValueError(f"{where}: duplicate record {name!r}")
                records[name] = []
            elif name is None:
                raise ValueError(f"{where}: sequence before first header")
            else:
                records[name].append(line.strip().upper())
    if not records:
        raise ValueError(f"{path}: no FASTA records")
    return {key: "".join(parts) for key, parts in records.items()}


def normalize_te_library(src: Path, dst: Path) -> list[tuple[str, int, str, str]]:
    """Write an A/C/G/T-only copy of the TE library.

    Each IUPAC ambiguity code is replaced by the base in ``AMBIGUITY``. Headers
    are cut to the record name and each sequence is written on one line.

    Args:
        src: TE library FASTA.
        dst: Output FASTA; overwritten.

    Returns:
        The changes as (record, 1-based position, original, replacement), in
        file order.

    Raises:
        ValueError: If ``src`` is not valid FASTA (see ``read_fasta``) or holds
            a base that is neither A/C/G/T nor in ``AMBIGUITY``, such as ``N``;
            ``N`` has no defensible single-base replacement.
    """
    changes = []
    lines = []
    for name, sequence in read_fasta(src).items():
        bases = list(sequence)
        for index, base in enumerate(bases):
            if base in "ACGT":
                continue
            if base not in AMBIGUITY:
                raise ValueError(f"{src}: {name} has unsupported base {base!r} (N?)")
            changes.append((name, index + 1, base, AMBIGUITY[base]))
            bases[index] = AMBIGUITY[base]
        lines.append(f">{name}\n{''.join(bases)}\n")
    dst.write_text("".join(lines))
    return changes


def build_config(events: list[TruthEvent]) -> dict:
    """Build the te-panel JSON: one reference genome, one sample haplotype.

    Args:
        events: Truth insertions; every one is carried by ``sample_h1``.

    Returns:
        The te-panel config as a JSON-ready dict; ``te_rules`` is keyed by
        family in sorted order and ``events`` keeps the input order.

    Raises:
        ValueError: If one family has events with different TSD lengths.
    """
    tsd_length: dict[str, int] = {}
    for event in events:
        length = len(event.tsd)
        if tsd_length.setdefault(event.family, length) != length:
            raise ValueError(f"family {event.family} has more than one TSD length")
    chromosomes = sorted({event.chrom for event in events})
    return {
        "schema_version": 1,
        "backbone": "backbone.fa",
        "te_library": "te_library.fa",
        "chromosomes": chromosomes,
        "flank_size": 5,
        "genomes": [
            {
                "id": "reference",
                "sample": "Reference",
                "haplotype": 1,
                "role": "reference",
            },
            {"id": "sample_h1", "sample": "Sample", "haplotype": 1, "role": "sample"},
        ],
        "te_rules": {
            family: {"family": family, "tsd_length": tsd_length[family]}
            for family in sorted(tsd_length)
        },
        "events": [
            {
                "id": event.event_id,
                "chromosome": event.chrom,
                "position": event.position,
                "te_id": event.family,
                "strand": event.strand,
                "carriers": ["sample_h1"],
            }
            for event in events
        ],
    }
