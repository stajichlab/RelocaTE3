"""Long-read data model (milestones L2 and L3).

Rules from the long-read proposal ("Designing for multiple references"):
every position names its reference (:class:`Locus`); evidence is separated
from the alignment format (:class:`PlacementRecord` holds no pysam objects);
signatures keep read-side evidence (sequence, flanks, query coordinates).
Coordinates are 0-based half-open; an insertion at ``pos`` has
``start == end == pos``.

L2 types (:class:`PlacementRecord`, :class:`Signature`, :class:`Cluster`)
describe read evidence; L3 types (:class:`TeHit`, :class:`Call`, with
:class:`CallTier` and :class:`FilterReason`) describe TE annotation of that
evidence and the insertion calls made from it.
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

    def __str__(self) -> str:
        """Return the value, so text output is the same on every Python."""
        return self.value

    def __format__(self, spec: str) -> str:
        """Format the value, so f-strings write it on every Python."""
        return format(self.value, spec)

    @property
    def minimap2_preset(self) -> str:
        """The minimap2 ``-x`` preset for this read type."""
        return {"ont": "map-ont", "ont-hq": "lr:hq", "hifi": "map-hifi",
                "clr": "map-pb"}[self.value]  # fmt: skip

    @property
    def default_cluster_window(self) -> int:
        """Breakpoint join distance; noisier reads shift breakpoints more."""
        return 50 if self in (ReadType.HIFI, ReadType.ONT_HQ) else 100

    @property
    def default_min_support(self) -> int:
        """Fewest supporting reads for a PASS call (noisier reads need more)."""
        return 2 if self in (ReadType.HIFI, ReadType.ONT_HQ) else 3

    @property
    def default_min_te_identity(self) -> float:
        """Lowest TE-hit identity used (provisional; see the L3 results note)."""
        return {"hifi": 0.85, "ont-hq": 0.80, "ont": 0.70, "clr": 0.70}[self.value]


class SignatureKind(str, Enum):
    """Kind of read evidence for an insertion.

    ``INS`` is a CIGAR insertion inside one alignment, ``SPLIT_INS`` an
    insertion implied by the read gap between consecutive alignments of one
    read, ``CLIP`` a long soft-clipped end (with its side).
    """

    INS = "INS"
    SPLIT_INS = "SPLIT_INS"
    CLIP = "CLIP"

    def __str__(self) -> str:
        """Return the value, so text output is the same on every Python."""
        return self.value

    def __format__(self, spec: str) -> str:
        """Format the value, so f-strings write it on every Python."""
        return format(self.value, spec)


class ClipSide(str, Enum):
    """Which end of an alignment a clip signature comes from."""

    LEFT = "left"  # clip before the aligned part (alignment start)
    RIGHT = "right"  # clip after the aligned part (alignment end)
    NONE = ""

    def __str__(self) -> str:
        """Return the value, so text output is the same on every Python."""
        return self.value

    def __format__(self, spec: str) -> str:
        """Format the value, so f-strings write it on every Python."""
        return format(self.value, spec)


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

    ``sequence`` is the stored query sequence exactly as in the BAM:
    reference-forward, soft-clipped bases included, hard-clipped bases
    absent. It is ``None`` only when the BAM has no SEQ (``*``), never just
    because of hard clips. ``left_clip``/``right_clip`` are soft-clip
    lengths at the stored left/right ends; ``left_hard_clip``/
    ``right_hard_clip`` say how many read bases are missing from
    ``sequence`` at each stored end. ``sequence`` is the full read only
    when both hard clips are 0.
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


class CallTier(str, Enum):
    """How a call is supported.

    ``SPANNING`` calls have at least one read that spans the insertion
    (``INS`` or ``SPLIT_INS``); ``CLIP_ONLY`` calls rest on clips alone.
    """

    SPANNING = "spanning"
    CLIP_ONLY = "clip_only"

    def __str__(self) -> str:
        """Return the value, so text output is the same on every Python."""
        return self.value

    def __format__(self, spec: str) -> str:
        """Format the value, so f-strings write it on every Python."""
        return format(self.value, spec)


class FilterReason(str, Enum):
    """A check a call failed. Member order is the order reasons are written."""

    NO_TE_HIT = "no_te_hit"
    FEW_TE_READS = "few_te_reads"
    FAMILY_DISCORDANT = "family_discordant"
    LOW_SUPPORT = "low_support"
    LOW_TE_FRACTION = "low_te_fraction"

    def __str__(self) -> str:
        """Return the value, so text output is the same on every Python."""
        return self.value

    def __format__(self, spec: str) -> str:
        """Format the value, so f-strings write it on every Python."""
        return format(self.value, spec)


NO_HIT = "."


@dataclass(frozen=True)
class TeHit:
    """One alignment of a signature's sequence to the TE library.

    A signature with no alignment gets one row with ``family == NO_HIT``,
    ``te_strand == "."`` and zero coordinates, so "looked, nothing found"
    differs from "not looked". ``te_strand`` is the TE's orientation on the
    reference-forward signature sequence. ``query_*`` are on that sequence;
    ``query_length`` is its length. ``identity`` is PAF matches / block length.
    """

    signature_id: str
    read_id: str
    family: str
    te_strand: str
    te_start: int
    te_end: int
    te_length: int
    query_start: int
    query_end: int
    query_length: int
    identity: float

    def __post_init__(self) -> None:
        """Reject a strand that does not fit the hit, or inverted spans."""
        allowed = (".",) if self.family == NO_HIT else ("+", "-")
        if self.te_strand not in allowed:
            raise ValueError(
                f"bad te_strand {self.te_strand!r} for family {self.family!r}"
            )
        if not 0 <= self.query_start <= self.query_end <= self.query_length:
            raise ValueError(f"bad query span {self.query_start}-{self.query_end} "
                             f"of {self.query_length}")  # fmt: skip
        if not 0 <= self.te_start <= self.te_end <= self.te_length:
            raise ValueError(
                f"bad TE span {self.te_start}-{self.te_end} of {self.te_length}"
            )

    @property
    def aligned(self) -> int:
        """Query bases covered by this hit."""
        return self.query_end - self.query_start


@dataclass(frozen=True)
class Call:
    """A TE insertion call made from one or more same-family clusters.

    ``locus`` is the insertion point (``start == end``; strand ``"+"``, as
    for clusters); the TE orientation is ``te_strand`` (``"."`` when no
    family). ``filters`` is empty for a PASS call, else the failed checks in
    :class:`FilterReason` order.
    """

    call_id: str
    locus: Locus
    tier: CallTier
    filters: tuple[FilterReason, ...]
    family: str
    te_strand: str
    family_fraction: float
    strand_fraction: float
    te_read_fraction: float
    n_reads: int
    n_ins: int
    n_split_ins: int
    n_clip_left: int
    n_clip_right: int
    median_insert_length: int
    median_te_coverage: float
    median_identity: float
    nested: bool
    cluster_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        """Reject a locus that is not an insertion point."""
        if self.locus.start != self.locus.end:
            raise ValueError(f"call locus must be an insertion point: {self.locus}")

    @property
    def passed(self) -> bool:
        """True when no filter failed."""
        return not self.filters
