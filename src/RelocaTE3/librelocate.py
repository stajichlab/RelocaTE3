"""RelocaTE3 to genotype individuals for Transposon polymorphism."""

from __future__ import annotations

import RelocaTE3.align as align
import RelocaTE3.ReadLibrary as ReadLibrary

# from Bio import SeqIO
# from Bio.Seq import Seq
# possibly support either a list of Seqs or a SeqLibrary fasta file?


class RelocaTE:
    """Process reads and mapping to identify transposon insertion and excision sites."""

    cpu_threads = 1         # number of CPU threads to use

    def __init__(self, threads: int = 1):
        """Initialize the RelocaTE object."""
        self.cpu_threads = threads

    def search_TE_containing_reads(self, seqreads: ReadLibrary, transposon_library: str,
                                   search_tool: str = "minimap2") -> int:
        """Search for sequence reads containing transposon sequences.

        Args:
            transposon_library (list[Seq]): List of sequences that definE transposon sequences to find in read file.

        Returns:
            int: number of reads (or read pairs) identified as containing transposon sequences.
        """
        if "minimap" in search_tool.lower():
            align._index_minimap(transposon_library)
#            align._map_minimap_genome()
        print(transposon_library)
