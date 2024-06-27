"""RelocaTE3 align and map reads to genome or Transposon libraries."""

from __future__ import annotations

import os
import subprocess
import warnings
from pathlib import Path

import RelocaTE3.ReadLibrary as ReadLibrary

# import tempfile
# import pysam
# from multiprocessing import Manager, Pool
# from multiprocessing.pool import ThreadPool


class Aligner:
    """Alignment Tool for reads."""
    minimap = "minimap2"
    bwa = "bwa"
    bwamem2 = "bwa-mem2"
    samtools = "samtools"
    verbose = False

    def _index_minimap(self, db: str, indexfile: str = "", force: bool = False) -> int:
        if len(indexfile) == 0:
            indexfile = Path(db+".mmi")
        else:
            indexfile = Path(indexfile)
        db = Path(db)
        if not db.exists():
            raise FileNotFoundError(f"Database file {db} does not exist.")
            return -1
        if indexfile.exists() and force is False:
            warnings.warn(f"minimap2 index file {indexfile} already exists will not recreate without force flag.")
            return 0
        # todo, potentially provide Stderr output if a verbose flag is passed?
        p = subprocess.run(
            [
                self.minimap,
                "-d",
                str(indexfile),
                str(db)
                ],
            stderr=None,
            capture_output=True,
            )
        if self.verbose:
            warnings.warn(p.stderr.decode("utf-8"))
        return 0

    def _map_minimap_library(self, transposon_library: str, reads: ReadLibrary, outdir: str, thread_count: int = 1) -> str:
        bamfile = os.path.join(outdir, reads.name + ".bam")
        bamfile

    def _map_minimap_genome(self, genome: str, reads: ReadLibrary, outdir: str, thread_count: int = 1) -> str:
        bamfile = os.path.join(outdir, reads.name + ".bam")
        bamfile
