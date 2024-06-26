"""RelocaTE3 align and map reads to genome or Transposon libraries."""

from __future__ import annotations

import os

import RelocaTE3.ReadLibrary as ReadLibrary

# import subprocess
# import tempfile
# import pysam
# from multiprocessing import Manager, Pool
# from multiprocessing.pool import ThreadPool


class Aligner:
    """Alignment Tool for reads."""
    minimap = "minimap2"


def _map_minimap_library(transposon_library: str, reads: ReadLibrary, outdir: str, thread_count: int = 1) -> str:
    bamfile = os.path.join(outdir, reads.name + ".bam")
    bamfile


def _map_minimap_genome(genome: str, reads: ReadLibrary, outdir: str, thread_count: int = 1) -> str:
    bamfile = os.path.join(outdir, reads.name + ".bam")
    bamfile
