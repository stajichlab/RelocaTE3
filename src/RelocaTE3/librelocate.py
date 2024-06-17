"""RelocaTE3 identify and process TE containing reads."""

from __future__ import annotations

import subprocess

from multiprocessing import Manager, Pool
from multiprocessing.pool import ThreadPool

from Bio import SeqIO
from Bio.Seq import Seq

class ReadLibrary(fileset: list[string] = None):
    def __init__(self, fileset):
        if fileset is not None:
            if len(fileset > 2):
                raise ValueError("Fileset needs to be either one or two files provided")
            self.is_paired = len(fileset) == 2
            self.fileset = fileset


    def search_reads(transposon_library: list[Seq]) -> int:
        print(transposon_library)