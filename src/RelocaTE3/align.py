"""RelocaTE3 align and map reads to genome or Transposon libraries."""

from __future__ import annotations

import os
# import re
import subprocess
import tempfile
import warnings
from pathlib import Path

import pysam

import RelocaTE3.ReadLibrary as ReadLibrary

# from multiprocessing import Manager, Pool
# from multiprocessing.pool import ThreadPool


class Aligner:
    """Alignment Tool for reads."""
    minimap = "minimap2"
    bwa = "bwa"
    bwamem2 = "bwa-mem2"
    samtools = "samtools"
    verbose = False
    # add threads as a state for this? but may use diff thread count for diff tools?

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
            if self.verbose:
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

    def _map_minimap_library(self, transposon_library: str, reads: ReadLibrary, outdir: str, tmpdir: str = "", cpu_threads: int = 1) -> list[str]:
        """Align short reads to transposon library to find those informative for insertions.

        input:
            - transposon_library: string of the fasta sequence library
            - reads: ReadLibrary object
            - outdir: string of where to write the resulting BAM file
            - tmpdir: string of tempfile (SAM) file creation - will use current directory if not
        """
        tmpdirhandle = None
        if len(tmpdir) == 0:
            tmpdirhandle = tempfile.TemporaryDirectory()
            tmpdir = tmpdirhandle.name
        elif not Path(tmpdir).exists():
            os.mkdir(tmpdir)
        # this may not be necessary/performance boost for Transposon library anyways so we might skip this
        # also best practice may be creating index on a SSD scratch volume anyways or loading memory
        # in general these are tiny DBs so it makes little difference I expect.
        index = transposon_library + ".mmi"
        self._index_minimap(transposon_library, str(index))

        temp_sam = os.path.join(tmpdir, "mm.sam")
        temp_bam = os.path.join(tmpdir, "mm.bam")

        # an option here is to run left and right separately as single --sr runs
        # then process the LEFT BAM/SAM file result, keep all mapping reads, AND retrieve the reads from the RIGHT file
        # then process the RIGHT BAM/SAM file and retrieve the LEFT reads that are the paired end of any match
        # do this without duplicating
        read_set = {'left': reads.left()}
        if reads.is_paired:
            read_set['right'] = reads.right()
        bam_files = []
        for direction in read_set:
            read_file = read_set[direction]
            p = subprocess.run(
                [
                    self.minimap,
                    "-t", str(cpu_threads),
                    "-a",
                    "-x",
                    "sr",   # we may need to play with the scoring here to see if this works well enough
                    "-o",
                    temp_sam,
                    str(index),
                    read_file,
                    ],
                stderr=None,
                capture_output=True,
            )
            if self.verbose:
                warnings.warn(p.stderr.decode("utf-8"))
            # potentially have the output to STDOUT and run
            # this as a pipe?
            pysam.sort("-o", temp_bam, temp_sam)
            bamfile = os.path.join(outdir, f"{reads.name}.{direction}.bam")

            p = subprocess.run([
                self.samtools,
                'view',
                '-o',
                bamfile,
                '-F',   # reads that do not match this next bitwise
                '0x4',  # unmapped
                temp_bam
            ])
            bam_files.append(bamfile)
        if tmpdirhandle:
            tmpdirhandle.cleanup()
        return bam_files

    def _map_minimap_genome(self, genome: str, reads: ReadLibrary, outdir: str, tmpdir: str = "", thread_count: int = 1) -> str:
        """Align reads to genome with minimap2 - this may not be best tool so testing."""
        bamfile = os.path.join(outdir, reads.name + ".bam")
        bamfile

    def _map_bwa_genome(self, genome: str, reads: ReadLibrary, outdir: str, thread_count: int = 1) -> str:
        """Align reads to genome with minimap2 - this may not be best tool so testing."""
        bamfile = os.path.join(outdir, reads.name + ".bam")
        bamfile
