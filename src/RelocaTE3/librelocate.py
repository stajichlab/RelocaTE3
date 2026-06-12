"""RelocaTE3 driver to genotype individuals for transposon polymorphism."""

from __future__ import annotations

from pathlib import Path

from RelocaTE3 import logger
from RelocaTE3.align import Aligner
from RelocaTE3.ReadLibrary import ReadLibrary, TrimmedReadLibrary
from RelocaTE3.trim import parse_te_alignments, trim_read


class RelocaTE:
    """Process reads to identify transposon insertion and excision sites."""

    def __init__(self, TElib: str = "", threads: int = 1, verbose: int = 0):
        """Initialize the RelocaTE object."""
        self.transposon_library = TElib
        self.cpu_threads = threads
        self.verbose = verbose

    def identify_TE_reads(
        self,
        seqreads: ReadLibrary,
        outdir: Path,
        TE_library: str = "",
        search_tool: str = "minimap2",
        len_cut_match: int = 10,
        len_cut_trim: int = 10,
        mismatch_allowance: int = 0,
    ) -> int:
        """Find reads containing transposon sequence and trim the TE portion.

        Maps the read library to the TE consensus library, trims the
        TE-matching portion from each read, and writes the flanking reads,
        TE portions, and read-to-TE table under ``outdir``.

        Returns the number of flanking (trimmed/supporting) reads written.
        """
        if not TE_library:
            TE_library = self.transposon_library

        if "minimap" not in search_tool.lower():
            raise NotImplementedError(f"Unsupported search tool: {search_tool}")

        alntool = Aligner(self.cpu_threads)
        alntool.verbose = self.verbose > 0
        alntool.index_minimap(TE_library)
        bamfiles = alntool.map_minimap_library(seqreads, str(outdir), TE_library)

        trimlib = self.trim_TE_reads(
            seqreads,
            bamfiles,
            len_cut_match=len_cut_match,
            len_cut_trim=len_cut_trim,
            mismatch_allowance=mismatch_allowance,
        )
        written = trimlib.write_reads(outdir)
        logger.info(
            "%s: %d reads matched a TE, %d flanking reads written",
            seqreads.name,
            len(trimlib.read_repeat),
            written,
        )
        return written

    def trim_TE_reads(
        self,
        reads: ReadLibrary,
        bamfiles: list[Path],
        len_cut_match: int = 10,
        len_cut_trim: int = 10,
        mismatch_allowance: int = 0,
    ) -> TrimmedReadLibrary:
        """Trim TE-matching sequence from reads aligned to the TE library.

        ``bamfiles`` are ordered as produced by :meth:`Aligner.map_minimap_library`
        (left/R1 first, then right/R2 for paired libraries).
        """
        trimlib = TrimmedReadLibrary(reads.name, reads.is_paired)
        for direction, bam in enumerate(bamfiles):
            coords = parse_te_alignments([bam])
            for aln in coords.values():
                trimmed = trim_read(
                    aln,
                    len_cut_match=len_cut_match,
                    len_cut_trim=len_cut_trim,
                    mismatch_allowance=mismatch_allowance,
                )
                if trimmed is not None:
                    trimlib.add_trimmed_read(direction, trimmed)
        return trimlib
