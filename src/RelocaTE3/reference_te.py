"""Steps 0 & 6: reference TE annotation, reference/shared insertions, FP filter.

Ports RelocaTE2's RepeatMasker handling (``existingTE_RM_ALL``) and the
reference-insertion role of ``relocaTE_absenceFinder.py``:

* parse a RepeatMasker ``.out`` into reference TE intervals (with intact-element
  detection) and write the ``existingTE.bed`` used downstream;
* call **reference/shared** insertions — junction-read clusters that coincide with
  an annotated reference TE boundary (the TE is present at that locus in both the
  reference and the resequenced sample);
* drop non-reference calls that overlap a known reference TE (the RelocaTE2 step-7
  ``bedtools intersect -v`` false-positive cleanup).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from RelocaTE3 import logger
from RelocaTE3.insertions import _stream_clusters
from RelocaTE3.models import Insertion

# how close (bp) a junction breakpoint must be to a reference TE edge to match it
BOUNDARY_WINDOW = 5


@dataclass
class ReferenceTE:
    """A transposable element annotated in the reference genome."""

    chrom: str
    start: int  # 1-based
    end: int  # 1-based inclusive
    name: str  # repeat name, e.g. "mPing"
    family: str  # repeat class/family, e.g. "DNA/MITE"
    strand: str  # '+' or '-'
    intact: bool  # full-length copy (both repeat ends present)


def _strip_parens(value: str) -> str:
    """Remove surrounding parentheses from a RepeatMasker position field."""
    return value.replace("(", "").replace(")", "")


def parse_repeatmasker(path: str | Path) -> list[ReferenceTE]:
    """Parse a RepeatMasker ``.out`` file into :class:`ReferenceTE` records.

    Reproduces RelocaTE2's column handling: a full-length (intact) element has
    repeat-start == 1 and repeat-left == 0 (accounting for strand).
    """
    tes: list[ReferenceTE] = []
    with open(path) as fh:
        for line in fh:
            if len(line.strip()) <= 2:
                continue
            unit = re.split(r"\s+", line.rstrip())
            # leading whitespace yields an empty unit[0]; normalize if absent
            if unit[0] != "":
                unit.insert(0, "")
            if len(unit) < 16 or not unit[1].isdigit():
                continue  # header / malformed line
            chrom, begin, end = unit[5], int(unit[6]), int(unit[7])
            rm_strand, name, family = unit[9], unit[10], unit[11]

            intact = False
            if rm_strand == "+":
                left = _strip_parens(unit[14])
                if unit[12] == "1" and left.isdigit() and int(left) == 0:
                    intact = True
                strand = "+"
            else:  # 'C' = reverse complement
                left = _strip_parens(unit[12])
                if unit[14] == "1" and left.isdigit() and int(left) == 0:
                    intact = True
                strand = "-"
            tes.append(ReferenceTE(chrom, begin, end, name, family, strand, intact))
    return tes


def write_existing_te_bed(tes: list[ReferenceTE], path: str | Path) -> None:
    """Write reference TEs as the RelocaTE2 ``existingTE.bed`` (6 columns)."""
    with open(path, "w") as fh:
        for te in tes:
            intact = 1 if te.intact else 0
            fh.write(
                f"{te.chrom}\t{te.start}\t{te.end}\t{te.name}:{te.start}-{te.end}\t{intact}\t{te.strand}\n"
            )


def _intact_boundary_index(tes: list[ReferenceTE]) -> dict[str, list[ReferenceTE]]:
    """Index intact reference TEs by chromosome for boundary lookups."""
    index: dict[str, list[ReferenceTE]] = {}
    for te in tes:
        if te.intact:
            index.setdefault(te.chrom, []).append(te)
    return index


def _matching_reference_te(
    chrom: str, position: int, index: dict[str, list[ReferenceTE]], window: int
) -> ReferenceTE | None:
    """Return an intact reference TE whose start or end is within ``window`` of ``position``."""
    for te in index.get(chrom, []):
        if abs(te.start - position) <= window or abs(te.end - position) <= window:
            return te
    return None


def find_reference_insertions(
    genome_bam: str,
    read_repeat: dict[str, tuple[str, str]],
    reference_tes: list[ReferenceTE],
    window: int = BOUNDARY_WINDOW,
) -> list[Insertion]:
    """Call reference/shared insertions from junction reads at reference TE edges.

    A junction-read cluster whose breakpoint coincides with an intact reference TE
    boundary indicates the element is present at that locus in the sample too
    (shared between reference and sample). Returns insertions sorted by position.
    """
    index = _intact_boundary_index(reference_tes)
    insertions: list[Insertion] = []
    for cluster in _stream_clusters(genome_bam, read_repeat):
        matched: dict[tuple[int, int], list] = {}
        for obs in cluster.junctions:
            te = _matching_reference_te(cluster.chrom, obs.position, index, window)
            if te is not None:
                matched.setdefault((te.start, te.end), []).append((obs, te))
        for (start, end), pairs in matched.items():
            te = pairs[0][1]
            left = sum(1 for obs, _ in pairs if obs.side == "left")
            right = sum(1 for obs, _ in pairs if obs.side == "right")
            insertions.append(
                Insertion(
                    chrom=cluster.chrom,
                    start=start,
                    end=end,
                    te_name=te.name,
                    strand=te.strand,
                    tsd="shared",
                    left_junction_reads=left,
                    right_junction_reads=right,
                    note="Shared, found in reference",
                    read_names=[obs.read_name for obs, _ in pairs],
                )
            )
    insertions.sort(key=lambda i: (i.chrom, i.start, i.end))
    logger.info("Called %d reference/shared insertions", len(insertions))
    return insertions


def filter_reference_overlaps(
    insertions: list[Insertion],
    reference_tes: list[ReferenceTE],
    window: int = BOUNDARY_WINDOW,
) -> list[Insertion]:
    """Drop non-reference calls overlapping a known reference TE of the same family.

    Mirrors RelocaTE2's step-7 ``bedtools intersect -v`` cleanup: a non-reference
    insertion that falls inside an annotated reference TE of the same family is
    almost certainly a mapping artifact, not a novel insertion.
    """
    by_chrom: dict[str, list[ReferenceTE]] = {}
    for te in reference_tes:
        by_chrom.setdefault(te.chrom, []).append(te)

    kept: list[Insertion] = []
    for ins in insertions:
        overlap = any(
            te.start - window <= ins.end
            and te.end + window >= ins.start
            and te.name == ins.te_name
            for te in by_chrom.get(ins.chrom, [])
        )
        if not overlap:
            kept.append(ins)
    return kept


def write_existing_te_bed_from_rm(
    rm_path: str | Path, bed_path: str | Path
) -> list[ReferenceTE]:
    """Parse a RepeatMasker ``.out`` and write ``existingTE.bed``; return the TEs."""
    tes = parse_repeatmasker(rm_path)
    write_existing_te_bed(tes, bed_path)
    return tes
