"""RelocaTE3 to genotype individuals for Transposon polymorphism."""

from __future__ import annotations

from pathlib import Path

import pysam
from Bio import SeqIO

from RelocaTE3.align import Aligner
from RelocaTE3.ReadLibrary import ReadLibrary


class RelocaTE:
    """Process reads and mapping to identify transposon insertion and excision sites."""

    cpu_threads = 1         # number of CPU threads to use
    transposon_library = None

    def __init__(self, TElib: str = "", threads: int = 1):
        """Initialize the RelocaTE object."""
        self.transposon_library = TElib
        self.cpu_threads = threads

    def identify_TE_reads(self, seqreads: ReadLibrary, outdir: Path, TE_library: str = "",
                          search_tool: str = "minimap2") -> int:
        """Search for sequence reads containing transposon sequences.

        Args:
            transposon_library (list[Seq]): List of sequences that definE transposon sequences to find in read file.

        Returns:
            int: number of reads (or read pairs) identified as containing transposon sequences.
        """
        if TE_library == "":
            TE_library = self.transposon_library

        alntool = Aligner(self.cpu_threads)
        if "minimap" in search_tool.lower():
            alntool.index_minimap(TE_library)
            bamfiles = alntool.map_minimap_library(TE_library, seqreads, outdir)
            TE_to_readinfo = self.process_TE_BAMfile(bamfiles, TE_library)
            print(TE_to_readinfo)
        print(TE_library)

    def process_TE_BAMfile(self, bamfiles: list[Path], TE_library: str) -> dict:
        """Function to process transposons aligned to reads."""
        TE_to_readinfo = dict()
        # get the lengths of the transposon reference sequences?
        TE_library_sequences = SeqIO.to_dict(TE_library, "fasta")
        print(TE_library_sequences)
        # we need to write code here to figure out the TEs
        # open BAM file with pysam
        # for each read, identify reference (eg mPing or Ping etc) aligned to
        # keep track of orientation and start/end of alignment to the transosable element seq
        for bam in bamfiles:
            readbam = pysam.AlignmentFile(bam, "rb")
            for readaln in readbam.fetch:
                print(readaln)
        return TE_to_readinfo
