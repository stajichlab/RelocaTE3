"""RelocaTE3 identify and process transposon sequence containing reads."""

from __future__ import annotations

# from Bio import SeqIO
from Bio.Seq import Seq

# import subprocess
# from multiprocessing import Manager, Pool
# from multiprocessing.pool import ThreadPool


class ReadLibrary():
    """Represent sequence library typically a single or paired-end FASTQ read files."""
    def __init__(self, fileset):
        """Initialize the ReadLibrary."""
        if fileset is not None:
            if len(fileset > 2):
                raise ValueError("Fileset needs to be either one or two files provided")
            self.is_paired = len(fileset) == 2
            self.fileset = fileset

    def search_TE_containing_reads(self, transposon_library: list[Seq], search_tool: str = "minimap2") -> int:
        """Search for sequence reads containing transposon sequences.

        Args:
            transposon_library (list[Seq]): List of sequences that definE transposon sequences to find in read file.

        Returns:
            int: number of reads (or read pairs) identified as containing transposon sequences.
        """
        print(transposon_library)
