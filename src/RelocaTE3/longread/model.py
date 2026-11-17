"""Long-read data model (milestone L2).

Rules from the long-read proposal ("Designing for multiple references"):
every position names its reference (:class:`Locus`); evidence is separated
from the alignment format (:class:`PlacementRecord` holds no pysam objects);
signatures keep read-side evidence (sequence, flanks, query coordinates).
Coordinates are 0-based half-open; an insertion at ``pos`` has
``start == end == pos``.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class ReadType(str, Enum):
    """Long-read technology, selecting the minimap2 preset."""

    ONT = "ont"
    ONT_HQ = "ont-hq"
    HIFI = "hifi"
    CLR = "clr"

    @property
    def minimap2_preset(self) -> str:
        """The minimap2 ``-x`` preset for this read type."""
        return {"ont": "map-ont", "ont-hq": "lr:hq", "hifi": "map-hifi",
                "clr": "map-pb"}[self.value]  # fmt: skip

    @property
    def default_cluster_window(self) -> int:
        """Breakpoint join distance; noisier reads shift breakpoints more."""
        return 50 if self in (ReadType.HIFI, ReadType.ONT_HQ) else 100


class SignatureKind(str, Enum):
    """Kind of read evidence for an insertion.

    ``INS`` is a CIGAR insertion inside one alignment, ``SPLIT_INS`` an
    insertion implied by the read gap between consecutive alignments of one
    read, ``CLIP`` a long soft-clipped end (with its side).
    """

    INS = "INS"
    SPLIT_INS = "SPLIT_INS"
    CLIP = "CLIP"


class ClipSide(str, Enum):
    """Which end of an alignment a clip signature comes from."""

    LEFT = "left"  # clip before the aligned part (alignment start)
    RIGHT = "right"  # clip after the aligned part (alignment end)
    NONE = ""


@dataclass(frozen=True)
class Locus:
    """A 0-based half-open interval on one sequence of one reference.

    Attributes:
        reference_id: Identifier of the reference (genome) the position is on.
        sequence_id: Sequence (chromosome/contig) name within that reference.
        start: 0-based start, inclusive.
        end: 0-based end, exclusive; equals ``start`` for an insertion point.
        strand: ``"+"`` or ``"-"``.
    """

    reference_id: str
    sequence_id: str
    start: int
    end: int
    strand: str  # "+" or "-"

    def __post_init__(self) -> None:
        """Reject negative or inverted intervals and unknown strands."""
        if self.start < 0 or self.end < self.start:
            raise ValueError(f"bad locus: start={self.start} end={self.end}")
        if self.strand not in ("+", "-"):
            raise ValueError(f"bad locus strand {self.strand!r}")

    @property
    def key(self) -> tuple[str, str]:
        """Grouping key: never the sequence name alone."""
        return (self.reference_id, self.sequence_id)


@dataclass(frozen=True)
class Reference:
    """A reference genome: its identifier, FASTA path and optional index.

    Attributes:
        id: Identifier used as ``Locus.reference_id``.
        fasta: Path to the reference FASTA.
        index: Path to an aligner index (e.g. minimap2 ``.mmi``), if any.
    """

    id: str
    fasta: str
    index: str | None = None


@dataclass(frozen=True)
class PlacementRecord:
    """One primary or supplementary alignment, format-independent.

    ``query_start``/``query_end`` are in read orientation and count hard
    clips. ``insertions`` holds ``(reference_pos, stored_qstart, stored_qend)``
    for CIGAR insertions, with query offsets into ``sequence``.
    ``sequence`` is the stored (reference-forward) read sequence, or ``None``
    when hard clips make it incomplete. ``left_clip``/``right_clip`` are
    soft-clip lengths at the stored left/right ends.
    """

    read_id: str
    read_length: int
    query_start: int
    query_end: int
    locus: Locus
    is_supplementary: bool
    mapq: int
    insertions: tuple[tuple[int, int, int], ...]
    left_clip: int
    right_clip: int
    sequence: str | None
    has_supplementary: bool
    left_hard_clip: int = 0
    right_hard_clip: int = 0


@dataclass(frozen=True)
class Signature:
    """One piece of read evidence for an insertion at ``locus``.

    ``sequence``, ``left_flank`` and ``right_flank`` are in
    reference-forward orientation; ``query_start``/``query_end`` are in
    read orientation. ``clip_side`` is ``ClipSide.NONE`` unless ``kind`` is
    ``SignatureKind.CLIP``.
    """

    signature_id: str
    kind: SignatureKind
    locus: Locus
    insert_length: int
    clip_side: ClipSide
    read_id: str
    query_start: int
    query_end: int
    mapq: int
    sequence: str
    left_flank: str
    right_flank: str


@dataclass(frozen=True)
class Cluster:
    """Signatures joined into one candidate insertion site.

    Counts are per signature kind (and clip side); ``n_reads`` counts
    distinct reads. ``signature_ids`` links back to the member signatures.
    """

    cluster_id: str
    locus: Locus
    n_reads: int
    n_ins: int
    n_split_ins: int
    n_clip_left: int
    n_clip_right: int
    insert_lengths: tuple[int, ...]
    signature_ids: tuple[str, ...]
