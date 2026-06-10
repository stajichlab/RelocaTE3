"""Step 5: cluster junction/supporting reads and call non-reference insertions.

A restructured port of RelocaTE2's ``relocaTE_insertionFinder.py``. The genome
BAM produced by Step 4 is streamed in coordinate order; junction reads (carrying
the ``:start:5`` / ``:end:3`` ... tags) mark insertion breakpoints. Each cluster
is split into one or more insertions by pairing left/right junction breakpoints
(``TSD_from_read_depth`` logic), the target-site duplication (TSD) is derived from
the overlap, supporting reads are counted, and junctions whose *full* (untrimmed)
read maps cleanly across the breakpoint are filtered as false junctions. Results
are written as GFF3 and a tab-delimited table.
"""

from __future__ import annotations

import re
from pathlib import Path

import pysam

from RelocaTE3 import logger
from RelocaTE3.models import Insertion, JunctionObservation

# read-name junction tag: <name>:(start|end):(5|3)
_JUNCTION_RE = re.compile(r":(start|end):([53])$")
# how far apart (bp) reads may be and still belong to one insertion cluster
RANGE_ALLOWANCE = 1000
# max separation (bp) between a left and right breakpoint to call a shared TSD
TSD_WINDOW = 100
# largest plausible TSD length (bp); wider overlaps are treated as TSD-unknown
MAX_TSD = 20
# a full read extending this far past a breakpoint indicates no insertion
FULLREAD_EXTEND = 10


def _junction_info(
    name: str, strand: str, gstart: int, gend: int
) -> tuple[str, int, str] | None:
    """Return (side, breakpoint, te_end) for a junction read, or None.

    The flank's TE-adjacent edge is the breakpoint. A right-junction read's
    breakpoint is the genomic start (left edge of the TSD); a left-junction
    read's breakpoint is the genomic end (right edge of the TSD).
    """
    m = _JUNCTION_RE.search(name)
    if not m:
        return None
    flank_side, te_end = m.group(1), m.group(2)
    if strand == "+":
        side, pos = ("right", gstart) if flank_side == "start" else ("left", gend)
    else:  # '-' strand inverts which read end is TE-adjacent
        side, pos = ("left", gend) if flank_side == "start" else ("right", gstart)
    return side, pos, te_end


def _te_family(read_repeat: dict[str, tuple[str, str]], read_name: str) -> str:
    """Best-effort TE family name for a junction read."""
    if read_name in read_repeat:
        return read_repeat[read_name][0]
    return "NA"


class _Cluster:
    """Reads grouped within ``RANGE_ALLOWANCE`` bp on one chromosome."""

    def __init__(self, chrom: str):
        self.chrom = chrom
        self.lo: int | None = None
        self.hi: int | None = None
        self.junctions: list[JunctionObservation] = []
        # supporting reads: (name, gstart, gend, strand)
        self.support: list[tuple[str, int, int, str]] = []

    def in_range(self, gstart: int, gend: int) -> bool:
        """True if a read at [gstart, gend] belongs to this cluster."""
        if self.lo is None:
            return True
        return gstart <= self.hi + RANGE_ALLOWANCE and gend >= self.lo - RANGE_ALLOWANCE

    def extend(self, gstart: int, gend: int) -> None:
        """Grow the cluster's coordinate span to include [gstart, gend]."""
        self.lo = gstart if self.lo is None else min(self.lo, gstart)
        self.hi = gend if self.hi is None else max(self.hi, gend)


def _stream_clusters(bam_path: str, read_repeat: dict[str, tuple[str, str]]):
    """Yield :class:`_Cluster` objects by streaming a coordinate-sorted BAM."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        current: _Cluster | None = None
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped:
                continue
            chrom = bam.get_reference_name(rec.reference_id)
            gstart = rec.reference_start + 1
            gend = (
                rec.reference_end
            )  # pysam end is 0-based exclusive == 1-based inclusive
            strand = "-" if rec.is_reverse else "+"
            name = rec.query_name

            if (
                current is None
                or current.chrom != chrom
                or not current.in_range(gstart, gend)
            ):
                if current is not None:
                    yield current
                current = _Cluster(chrom)

            current.extend(gstart, gend)
            info = _junction_info(name, strand, gstart, gend)
            if info is not None:
                side, pos, te_end = info
                current.junctions.append(
                    JunctionObservation(
                        name, side, pos, strand, _te_family(read_repeat, name), te_end
                    )
                )
            else:
                current.support.append((name, gstart, gend, strand))
        if current is not None:
            yield current


def _group_by_position(
    observations: list[JunctionObservation],
) -> dict[int, list[JunctionObservation]]:
    """Group junction observations by their breakpoint position."""
    grouped: dict[int, list[JunctionObservation]] = {}
    for obs in observations:
        grouped.setdefault(obs.position, []).append(obs)
    return grouped


def _pair_breakpoints(
    left_pos: list[int], right_pos: list[int]
) -> list[tuple[int | None, int | None]]:
    """Pair left/right breakpoints into sub-insertions (RelocaTE2 pairing).

    Each returned tuple is (left_position, right_position); either may be None
    for a one-sided junction. Left and right are paired when within ``TSD_WINDOW``.
    """
    nL, nR = len(left_pos), len(right_pos)
    if nL and nR:
        if nL == 1 and nR == 1:
            if abs(left_pos[0] - right_pos[0]) > TSD_WINDOW:
                return [(left_pos[0], None), (None, right_pos[0])]
            return [(left_pos[0], right_pos[0])]
        # greedy nearest-neighbour pairing within the TSD window
        pairs: list[tuple[int | None, int | None]] = []
        remaining_right = sorted(right_pos)
        for lp in sorted(left_pos):
            best = None
            for rp in remaining_right:
                if best is None or abs(rp - lp) < abs(best - lp):
                    best = rp
            if best is not None and abs(best - lp) <= TSD_WINDOW:
                pairs.append((lp, best))
                remaining_right.remove(best)
            else:
                pairs.append((lp, None))
        for rp in remaining_right:
            pairs.append((None, rp))
        return pairs
    if nL:
        return [(lp, None) for lp in sorted(left_pos)]
    return [(None, rp) for rp in sorted(right_pos)]


def _make_insertion(
    chrom: str,
    left_reads: list[JunctionObservation],
    right_reads: list[JunctionObservation],
    genome: pysam.FastaFile,
) -> Insertion:
    """Build an :class:`Insertion` from the left/right junction reads of one site."""
    junctions = left_reads + right_reads
    te_names = [j.te_name for j in junctions if j.te_name != "NA"]
    te_name = max(set(te_names), key=te_names.count) if te_names else "NA"

    orients = [j.te_orientation for j in junctions]
    strand = "+" if orients.count("+") >= orients.count("-") else "-"

    if left_reads and right_reads:
        i_end = left_reads[0].position  # right edge of TSD
        i_start = right_reads[0].position  # left edge of TSD
        if 0 < (i_end - i_start + 1) <= MAX_TSD:
            tsd = _fetch_tsd(genome, chrom, i_start, i_end)
        else:
            i_start = i_end = min(i_start, i_end)
            tsd = "UNK"
    else:
        present = left_reads or right_reads
        i_start = i_end = present[0].position
        tsd = "UNK"

    return Insertion(
        chrom=chrom,
        start=i_start,
        end=i_end,
        te_name=te_name,
        strand=strand,
        tsd=tsd,
        left_junction_reads=len(left_reads),
        right_junction_reads=len(right_reads),
        read_names=[j.read_name for j in junctions],
    )


def _call_insertions(cluster: _Cluster, genome: pysam.FastaFile) -> list[Insertion]:
    """Split a cluster into one or more insertions by pairing breakpoints."""
    left = _group_by_position([j for j in cluster.junctions if j.side == "left"])
    right = _group_by_position([j for j in cluster.junctions if j.side == "right"])
    if not left and not right:
        return []

    insertions: list[Insertion] = []
    for lp, rp in _pair_breakpoints(list(left), list(right)):
        left_reads = left.get(lp, []) if lp is not None else []
        right_reads = right.get(rp, []) if rp is not None else []
        ins = _make_insertion(cluster.chrom, left_reads, right_reads, genome)
        _count_support(ins, cluster)
        insertions.append(ins)
    return insertions


def _fetch_tsd(genome: pysam.FastaFile, chrom: str, start: int, end: int) -> str:
    """Fetch the genomic TSD sequence (1-based inclusive), '.'-padded on failure."""
    try:
        seq = genome.fetch(chrom, start - 1, end)
        return seq.upper() if seq else "." * (end - start + 1)
    except (KeyError, ValueError):
        return "." * (end - start + 1)


def _count_support(ins: Insertion, cluster: _Cluster) -> None:
    """Count bracketing supporting reads (RelocaTE2 ``Supporting_count`` rule)."""
    left = right = 0
    for _name, gstart, gend, strand in cluster.support:
        if strand == "+" and gend <= ins.start:
            left += 1
        elif strand == "-" and gstart >= ins.end:
            right += 1
    ins.left_support_reads = left
    ins.right_support_reads = right


def _load_fullread_spans(
    fullreads_bam: str | None,
) -> dict[str, list[tuple[str, int, int]]]:
    """Map read name -> genome spans for full (untrimmed) junction reads."""
    spans: dict[str, list[tuple[str, int, int]]] = {}
    if not fullreads_bam or not Path(fullreads_bam).exists():
        return spans
    with pysam.AlignmentFile(fullreads_bam, "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped:
                continue
            chrom = bam.get_reference_name(rec.reference_id)
            spans.setdefault(rec.query_name, []).append(
                (chrom, rec.reference_start + 1, rec.reference_end)
            )
    return spans


def _strip_junction_tag(name: str) -> str:
    """Remove the :start:5 / :end:3 junction suffix to recover the read name."""
    return _JUNCTION_RE.sub("", name)


def _is_false_junction(
    ins: Insertion, fullread_spans: dict[str, list[tuple[str, int, int]]]
) -> bool:
    """True if the full reads map across the breakpoint, indicating no insertion.

    Mirrors RelocaTE2: if >=30% of the left junction reads AND >=30% of the right
    junction reads have a full read spanning the breakpoint (with margin), the
    site is a reference locus, not an insertion.
    """
    if not fullread_spans:
        return False
    left_total = ins.left_junction_reads
    right_total = ins.right_junction_reads
    if left_total == 0 or right_total == 0:
        return False

    left_full = right_full = 0
    # read_names are ordered left-reads first, then right-reads
    left_names = ins.read_names[:left_total]
    right_names = ins.read_names[left_total:]
    bp_left, bp_right = ins.start, ins.end
    for tagged in left_names:
        if _maps_through(
            fullread_spans.get(_strip_junction_tag(tagged)), ins.chrom, bp_right
        ):
            left_full += 1
    for tagged in right_names:
        if _maps_through(
            fullread_spans.get(_strip_junction_tag(tagged)), ins.chrom, bp_left
        ):
            right_full += 1
    return left_full >= 0.3 * left_total and right_full >= 0.3 * right_total


def _maps_through(
    spans: list[tuple[str, int, int]] | None, chrom: str, breakpoint: int
) -> bool:
    """True if any full-read span covers the breakpoint with margin on both sides."""
    if not spans:
        return False
    for c, s, e in spans:
        if (
            c == chrom
            and s <= breakpoint - FULLREAD_EXTEND
            and e >= breakpoint + FULLREAD_EXTEND
        ):
            return True
    return False


def find_insertions(
    genome_bam: str,
    read_repeat: dict[str, tuple[str, str]],
    genome_fasta: str,
    fullreads_bam: str | None = None,
    required_junction_reads: int = 1,
) -> list[Insertion]:
    """Call non-reference insertions from the Step-4 genome BAM.

    Keeps clusters with at least ``required_junction_reads`` junction reads on
    either side, dropping false junctions identified via ``fullreads_bam`` (the
    untrimmed junction reads aligned to the genome). Returns insertions sorted by
    chromosome and position.
    """
    fullread_spans = _load_fullread_spans(fullreads_bam)
    insertions: list[Insertion] = []
    n_false = 0
    with pysam.FastaFile(genome_fasta) as genome:
        for cluster in _stream_clusters(genome_bam, read_repeat):
            for ins in _call_insertions(cluster, genome):
                if (
                    ins.left_junction_reads < required_junction_reads
                    and ins.right_junction_reads < required_junction_reads
                ):
                    continue
                if _is_false_junction(ins, fullread_spans):
                    n_false += 1
                    continue
                insertions.append(ins)
    insertions.sort(key=lambda i: (i.chrom, i.start, i.end))
    logger.info(
        "Called %d non-reference insertions (%d false junctions filtered)",
        len(insertions),
        n_false,
    )
    return insertions


def write_insertions_gff(
    insertions: list[Insertion],
    path: str | Path,
    sample: str,
    source: str = "RelocaTE3",
) -> None:
    """Write insertions as GFF3 with the RelocaTE2 attribute set."""
    with open(path, "w") as fh:
        for ins in insertions:
            attrs = (
                f"ID={ins.feature_id};Name={ins.te_name};TSD={ins.tsd};Note={ins.note};"
                f"Right_junction_reads={ins.right_junction_reads};"
                f"Left_junction_reads={ins.left_junction_reads};"
                f"Right_support_reads={ins.right_support_reads};"
                f"Left_support_reads={ins.left_support_reads};"
            )
            fh.write(
                f"{ins.chrom}\t{source}\t{sample}\t{ins.start}\t{ins.end}\t.\t{ins.strand}\t.\t{attrs}\n"
            )


def write_insertions_txt(insertions: list[Insertion], path: str | Path) -> None:
    """Write a tab-delimited insertion summary table."""
    header = [
        "chrom",
        "start",
        "end",
        "TE",
        "strand",
        "TSD",
        "right_junction",
        "left_junction",
        "right_support",
        "left_support",
    ]
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for ins in insertions:
            fh.write(
                "\t".join(
                    str(x)
                    for x in (
                        ins.chrom,
                        ins.start,
                        ins.end,
                        ins.te_name,
                        ins.strand,
                        ins.tsd,
                        ins.right_junction_reads,
                        ins.left_junction_reads,
                        ins.right_support_reads,
                        ins.left_support_reads,
                    )
                )
                + "\n"
            )
