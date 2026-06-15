"""Reference TE annotation helpers for steps 0 and 6.

Includes both the existing-copy annotator and the newer RepeatMasker/reference
insertion helpers used by the pipeline code.
"""

from __future__ import annotations

import re
import subprocess
from collections import defaultdict
from pathlib import Path

from RelocaTE3 import logger

# RepeatMasker / "rm" / ".out" reference annotation files are handled specially.
_RM_HINT = re.compile(r"repeatmasker|rm|\.out", re.IGNORECASE)


class ReferenceTEAnnotator:
    """Locate existing transposon copies in the reference genome."""

    def __init__(self, minimap: str = "minimap2", threads: int = 1, verbose: int = 0):
        """Initialize the annotator.

        Args:
            minimap: path to the ``minimap2`` executable.
            threads: CPU threads for minimap2.
            verbose: verbosity level.
        """
        self.minimap = minimap
        self.threads = threads
        self.verbose = verbose

    def annotate_minimap(
        self,
        te_library: Path,
        genome_fasta: Path,
        outdir: Path,
        min_identity: float = 0.8,
        min_coverage: float = 0.8,
    ) -> Path:
        """Align the TE library to the genome and write a BED of existing copies.

        Args:
            te_library: FASTA of transposon sequences (queries).
            genome_fasta: reference genome FASTA (target).
            outdir: directory to write ``existingTE.bed`` into.
            min_identity: minimum gap-compressed identity (matches / aln block).
            min_coverage: minimum fraction of the TE query covered by the alignment.

        Returns:
            Path to the written ``existingTE.bed``.
        """
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        bed_path = outdir / "existingTE.bed"

        # asm20 tolerates the divergence typical between a TE consensus and its
        # genomic copies; secondary hits keep every copy of a repeat family.
        cmd = [
            self.minimap,
            "-c",
            "-x",
            "asm20",
            "--secondary=yes",
            "-N",
            "100",
            "-p",
            "0.1",
            "-t",
            str(self.threads),
            str(genome_fasta),
            str(te_library),
        ]
        if self.verbose:
            logger.info("Running: %s", " ".join(cmd))
        proc = subprocess.run(cmd, capture_output=True, check=True, text=True)

        rows = []
        for line in proc.stdout.splitlines():
            hit = self._parse_paf_line(line, min_identity, min_coverage)
            if hit is not None:
                rows.append(hit)
        rows.sort(key=lambda r: (r[0], r[1]))
        with open(bed_path, "w") as out:
            for chrom, start, end, name, score, strand in rows:
                out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
        logger.info("Wrote %d existing-TE copies to %s", len(rows), bed_path)
        return bed_path

    @staticmethod
    def _parse_paf_line(line, min_identity, min_coverage):
        """Parse one PAF line into a BED row, or None if below thresholds."""
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 12:
            return None
        qname = cols[0]
        qlen = int(cols[1])
        strand = cols[4]
        tname = cols[5]
        tstart = int(cols[7])
        tend = int(cols[8])
        matches = int(cols[9])
        aln_len = int(cols[10])
        if aln_len == 0 or qlen == 0:
            return None
        identity = matches / aln_len
        coverage = (int(cols[3]) - int(cols[2])) / qlen
        if identity < min_identity or coverage < min_coverage:
            return None
        return (tname, tstart, tend, qname, int(identity * 1000), strand)

    # ------------------------------------------------------------------
    # boundary table consumed by the insertion finder (step 5)
    # ------------------------------------------------------------------
    @classmethod
    def load_existing_te(cls, reference_ins: Path | str, target: str = "ALL") -> dict:
        """Build the ``existingTE_inf`` boundary table from a RM .out or BED file.

        Args:
            reference_ins: a RepeatMasker ``.out`` file or a BED of existing copies.
            target: chromosome to restrict to, or ``"ALL"``.

        Returns:
            Nested dict ``{chrom: {"start": {pos: 1}, "end": {pos: 1}}}``.
        """
        existing = defaultdict(lambda: {"start": {}, "end": {}})
        reference_ins = Path(reference_ins)
        if not reference_ins.exists() or reference_ins.stat().st_size == 0:
            logger.info(
                "Existing TE file does not exist or is empty: %s", reference_ins
            )
            return existing
        if _RM_HINT.search(str(reference_ins)):
            cls._load_repeatmasker(reference_ins, existing, target)
        else:
            cls._load_bed(reference_ins, existing, target)
        return existing

    @staticmethod
    def _record_boundaries(existing, chrom, begin, end):
        """Mark +/- 2 bp windows around a copy's start and end coordinates."""
        for i in range(begin - 2, begin + 3):
            existing[chrom]["start"][i] = 1
        for i in range(end - 2, end + 3):
            existing[chrom]["end"][i] = 1

    @classmethod
    def _load_repeatmasker(cls, infile, existing, target):
        """Parse a RepeatMasker ``.out`` file into the boundary table."""
        with open(infile) as handle:
            for line in handle:
                line = line.rstrip()
                if len(line) <= 2:
                    continue
                unit = re.split(r"\s+", line)
                # normalise so real columns start at index 1 regardless of
                # whether the line had leading whitespace (RM .out usually does)
                if unit[0] != "":
                    unit.insert(0, "")
                # unit[5]=chrom, unit[6]=begin, unit[7]=end, unit[9]=strand(+/C)
                if len(unit) < 10 or not unit[6].isdigit() or not unit[7].isdigit():
                    continue
                chrom = unit[5]
                if target != "ALL" and chrom != target:
                    continue
                if unit[9] in ("+", "C"):
                    cls._record_boundaries(existing, chrom, int(unit[6]), int(unit[7]))

    @classmethod
    def _load_bed(cls, infile, existing, target):
        """Parse a BED of existing copies into the boundary table."""
        with open(infile) as handle:
            for line in handle:
                line = line.rstrip()
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                cols = line.split("\t")
                if len(cols) < 3:
                    continue
                chrom = cols[0]
                if target != "ALL" and chrom != target:
                    continue
                # BED is 0-based half-open; convert to 1-based inclusive boundaries.
                cls._record_boundaries(existing, chrom, int(cols[1]) + 1, int(cols[2]))
from dataclasses import dataclass
from pathlib import Path

from RelocaTE3 import logger
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
    from RelocaTE3.insertions import _stream_clusters

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
