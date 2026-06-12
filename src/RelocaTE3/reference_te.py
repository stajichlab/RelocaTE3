"""Annotate transposon copies already present in the reference genome (RelocaTE2 step 0).

RelocaTE2 used BLAT to align the TE library against the reference genome and
recorded the boundaries of every existing copy so that step 5 can avoid calling
those known/reference insertions as novel. This module performs the same job
with minimap2, and can also ingest a pre-computed RepeatMasker ``.out`` file
(the ``--reference-ins`` path in RelocaTE2).

The boundary table it builds (``existingTE_inf``) maps, per chromosome, the
``start`` and ``end`` coordinates of known TE copies (padded +/- 2 bp) so that a
junction landing on a reference TE edge is recognised and skipped.
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
