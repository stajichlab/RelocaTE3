"""Core data structures for RelocaTE3.

These dataclasses replace the deeply nested ``defaultdict`` structures used in
RelocaTE2 to pass state between functions, making the pipeline testable.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class JunctionType(str, Enum):
    """Where a read aligns relative to a transposable element (TE).

    A read that crosses a TE boundary is a *junction* read; one wholly inside
    the TE is a *middle* (supporting) read. The tag value matches the read-name
    suffix convention inherited from RelocaTE2 so downstream clustering can
    recover the junction side from a read name.
    """

    FIVE_PRIME = "5"  # read spans the 5' end of the TE
    THREE_PRIME = "3"  # read spans the 3' end of the TE
    MIDDLE = "middle"  # read lies entirely within the TE (supporting read)
    NONE = "none"  # read did not meet any junction/middle criteria


@dataclass
class TEReadAlignment:
    """Best alignment of a single sequencing read to a TE consensus.

    Coordinates are 0-based; ``qend``/``tend`` are inclusive (matching the
    convention used throughout the original RelocaTE2 trimming code).
    """

    read_name: str
    read_length: int
    qstart: int  # alignment start in the read (0-based)
    qend: int  # alignment end in the read (0-based, inclusive)
    te_name: str
    te_length: int
    tstart: int  # alignment start in the TE (0-based)
    tend: int  # alignment end in the TE (0-based, inclusive)
    match: int
    mismatch: int
    strand: str  # '+' or '-'
    boundary: int  # 0-4: count of read/TE ends reached by the alignment
    seq: str = ""  # read sequence as stored in the BAM record
    qual: str = ""  # read base qualities (ASCII phred+33), may be empty

    def is_better_than(self, other: TEReadAlignment) -> bool:
        """Return True if this alignment should replace ``other`` as the best.

        Prefer more boundaries reached (a read+TE end pair is the strongest
        evidence of a true junction); break ties on longer match length.
        """
        if self.boundary > other.boundary:
            return True
        if self.boundary == other.boundary:
            return self.match > other.match
        return False


@dataclass
class TrimmedRead:
    """A read after the TE-matching portion has been trimmed away.

    ``read_name`` carries the RelocaTE2 junction tag suffix (``:start:5``,
    ``:end:5``, ``:start:3``, ``:end:3`` or ``:middle``) so that downstream
    clustering can recover the junction side from the name alone. ``seq``/``qual``
    are the *flanking* (non-TE) portion to be re-aligned to the genome; for a
    MIDDLE read they are the full read. ``te_portion`` is the TE-matching
    sequence (forward-oriented), kept only for junction reads.
    """

    read_name: str
    seq: str
    qual: str
    te_name: str
    strand: str
    junction_type: JunctionType
    te_portion: str = ""


@dataclass
class JunctionObservation:
    """A junction read mapped to the genome, marking one edge of an insertion.

    ``side`` is 'left' or 'right': a left-junction read's flank lies to the left
    of the insertion (TE attached on its right), a right-junction read's flank to
    the right. ``position`` is the 1-based genomic breakpoint (right edge of the
    TSD for a left read, left edge for a right read).
    """

    read_name: str
    side: str  # 'left' or 'right'
    position: int
    strand: str  # genomic mapped strand of the flank
    te_name: str
    te_end: str = ""  # '5' or '3': which end of the TE this read spans

    @property
    def te_orientation(self) -> str:
        """Infer TE insertion orientation (+/-) from junction side and TE end.

        The TE 5' end lying at the left edge (or 3' end at the right edge) means
        the element is in forward (+) orientation; the opposite implies reverse.
        """
        if (self.te_end == "5" and self.side == "left") or (
            self.te_end == "3" and self.side == "right"
        ):
            return "+"
        if (self.te_end == "5" and self.side == "right") or (
            self.te_end == "3" and self.side == "left"
        ):
            return "-"
        return "+"


@dataclass
class Insertion:
    """A called transposable-element insertion site.

    Coordinates are 1-based inclusive; ``start``/``end`` span the target-site
    duplication (TSD). Counts follow the RelocaTE2 GFF attribute set.
    """

    chrom: str
    start: int
    end: int
    te_name: str
    strand: str
    tsd: str
    left_junction_reads: int = 0
    right_junction_reads: int = 0
    left_support_reads: int = 0
    right_support_reads: int = 0
    note: str = "Non-reference, not found in reference"
    read_names: list = None  # type: ignore[assignment]
    # genotyping (Step 7); populated by characterize.py
    status: str = ""  # homozygous / heterozygous / somatic_insertion / ...
    spanners: int = 0  # reference-allele reads mapping cleanly across the site

    def __post_init__(self) -> None:
        """Default the read-name list to empty."""
        if self.read_names is None:
            self.read_names = []

    @property
    def feature_id(self) -> str:
        """Stable GFF feature id: ``repeat_<chrom>_<start>_<end>``."""
        return f"repeat_{self.chrom}_{self.start}_{self.end}"

    @property
    def avg_flankers(self) -> float:
        """RelocaTE2 ``average_flankers``: total junction reads divided by two."""
        return (self.left_junction_reads + self.right_junction_reads) / 2
