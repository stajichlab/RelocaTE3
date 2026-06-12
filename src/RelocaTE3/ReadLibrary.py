"""ReadLibrary for paired-end or single-end read file set."""

from __future__ import annotations

import os
from pathlib import Path

from RelocaTE3.models import JunctionType, TrimmedRead


class ReadLibrary:
    """Represent a sequence library: one (single-end) or two (paired-end) FASTQ files."""

    def __init__(self, fileset: list[str], ind_name: str):
        """Initialize the ReadLibrary.

        Args:
            fileset: list of one or two FASTQ paths (left, optional right).
            ind_name: name of the individual/strain (used as a read-group label).
        """
        if fileset is None or len(fileset) == 0:
            raise ValueError("Fileset needs at least one read file")
        if len(fileset) > 2:
            raise ValueError("Fileset needs to be either one or two files provided")
        self.is_paired = len(fileset) == 2
        self.file_set = list(fileset)
        self.name = ind_name

    def left(self) -> str:
        """Return the left/R1 read file."""
        return self.file_set[0]

    def right(self) -> str | None:
        """Return the right/R2 read file, or None for single-end."""
        return self.file_set[1] if self.is_paired else None


class TrimmedReadLibrary:
    """Reads from one sample after trimming the TE-matching portion.

    Holds the flanking (trimmed) reads that will be re-aligned to the genome,
    organized by direction (0 = left/R1, 1 = right/R2), plus the TE portions and
    read-to-TE assignments needed by downstream steps.
    """

    def __init__(self, ind_name: str = "", is_paired: bool = False):
        """Initialize an empty trimmed library for sample ``ind_name``."""
        self.name = ind_name
        self.is_paired = is_paired
        # flanking (and middle) reads to re-map to the genome, per direction
        self.flanking_reads: list[list[TrimmedRead]] = [[], []]
        # TE-matching portions, kept for 5'/3' junction reads
        self.five_prime: list[TrimmedRead] = []
        self.three_prime: list[TrimmedRead] = []
        # read_name -> (te_name, strand) for every read assigned to a TE
        self.read_repeat: dict[str, tuple[str, str]] = {}

    def add_trimmed_read(self, direction: int, trimmed: TrimmedRead) -> None:
        """Record a trimmed read for ``direction`` (0=left/R1, 1=right/R2)."""
        self.flanking_reads[direction].append(trimmed)
        self.read_repeat[trimmed.read_name] = (trimmed.te_name, trimmed.strand)
        if trimmed.junction_type is JunctionType.FIVE_PRIME and trimmed.te_portion:
            self.five_prime.append(trimmed)
        elif trimmed.junction_type is JunctionType.THREE_PRIME and trimmed.te_portion:
            self.three_prime.append(trimmed)

    @property
    def num_flank_reads(self) -> int:
        """Total number of flanking/middle reads across both directions."""
        return sum(len(d) for d in self.flanking_reads)

    def write_reads(self, outdir: str) -> int:
        """Write flanking FASTQ, TE-portion FASTA, and read-repeat table.

        Returns the number of flanking reads written.
        """
        outdir = Path(outdir)
        flank_dir = outdir / "flanking"
        portion_dir = outdir / "te_portions"
        containing_dir = outdir / "te_containing"
        for d in (flank_dir, portion_dir, containing_dir):
            os.makedirs(d, exist_ok=True)

        written = 0
        directions = ["left", "right"]
        for idx, reads in enumerate(self.flanking_reads):
            if not reads:
                continue
            fq_path = flank_dir / f"{self.name}.{directions[idx]}.flankingReads.fq"
            with open(fq_path, "w") as fh:
                for r in reads:
                    fh.write(f"@{r.read_name}\n{r.seq}\n+\n{r.qual}\n")
                    written += 1

        for label, reads in (
            ("five_prime", self.five_prime),
            ("three_prime", self.three_prime),
        ):
            fa_path = portion_dir / f"{self.name}.{label}.fa"
            with open(fa_path, "w") as fh:
                for r in reads:
                    fh.write(f">{r.read_name}\n{r.te_portion}\n")

        rr_path = containing_dir / f"{self.name}.read_repeat_name.txt"
        with open(rr_path, "w") as fh:
            for read_name, (te_name, strand) in self.read_repeat.items():
                fh.write(f"{read_name}\t{te_name}\t{strand}\n")

        return written
