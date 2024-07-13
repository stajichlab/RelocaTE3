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
    cpu_threads = 1
    aligntool = "minimap2"  # or BLAT or bwa OR others?
    # add threads as a state for this? but may use diff thread count for diff tools?

    def __init__(self, threads: int = 1, default_aligner: str = "minimap2"):
        """Initialize the aligner object."""
        self.cpu_threads = threads
        self.aligntool = default_aligner

    def index_minimap(self, db: str, indexfile: str = "", force: bool = False) -> int:
        """Create minimap index of target DB."""
        indexfile = Path(f"{db}.mmi") if not indexfile else Path(indexfile)
        db = Path(db)
        if not db.exists():
            raise FileNotFoundError(f"Database file {db} does not exist.")
        if indexfile.exists() and not force:
            if self.verbose:
                warnings.warn(
                    f"minimap2 index file {indexfile} already exists will not recreate without force flag."
                )
            return 0
        # todo, potentially provide Stderr output if a verbose flag is passed?
        p = subprocess.run(
            [self.minimap, "-d", str(indexfile), str(db)],
            stderr=None,
            capture_output=True,
            check=True,
        )
        if self.verbose:
            warnings.warn(p.stderr.decode("utf-8"))
        return 0

    def map_minimap_library(
        self,
        reads: ReadLibrary,
        outdir: str,
        transposon_library: str,
        tmpdir: str = "",
        cpu_threads: int = 0,
    ) -> list[Path]:
        """Align short reads to transposon library to find those informative for insertions.

        input:
            - transposon_library: string of the fasta sequence library
            - reads: ReadLibrary object
            - outdir: string of where to write the resulting BAM file
            - tmpdir: string of tempfile (SAM) file creation - will use current directory if not
        """
        if cpu_threads <= 0:
            cpu_threads = self.cpu_threads
        tmpdirhandle = None
        if not tmpdir:
            tmpdirhandle = tempfile.TemporaryDirectory()
            tmpdir = tmpdirhandle.name
        elif not Path(tmpdir).exists():
            os.mkdir(tmpdir)
        # this may not be necessary/performance boost for Transposon library anyways so we might skip this
        # also best practice may be creating index on a SSD scratch volume anyways or loading memory
        # in general these are tiny DBs so it makes little difference I expect.
        index = f"{transposon_library}.mmi"
        self.index_minimap(transposon_library, str(index))

        temp_sam = os.path.join(tmpdir, "mm.sam")
        temp_bam = os.path.join(tmpdir, "mm.bam")

        # an option here is to run left and right separately as single --sr runs
        # then process the LEFT BAM/SAM file result, keep all mapping reads, AND retrieve the reads from the RIGHT file
        # then process the RIGHT BAM/SAM file and retrieve the LEFT reads that are the paired end of any match
        # do this without duplicating
        read_set = {"left": reads.left()}
        if reads.is_paired:
            read_set["right"] = reads.right()
        bam_files = []
        for direction, read_file in read_set.items():
            p = subprocess.run(
                [
                    self.minimap,
                    "-t",
                    str(cpu_threads),
                    "-a",
                    "-x",
                    "sr",
                    "-o",
                    temp_sam,
                    str(index),
                    read_file,
                ],
                stderr=None,
                capture_output=True,
                check=True,
            )
            if self.verbose:
                warnings.warn(p.stderr.decode("utf-8"))
            # potentially have the output to STDOUT and run
            # this as a pipe?
            pysam.sort("-o", temp_bam, temp_sam)
            bamfile = os.path.join(outdir, f"{reads.name}.{direction}.bam")

            subprocess.run(
                [
                    self.samtools,
                    "view",
                    "-o",
                    bamfile,
                    "-F",  # reads that do not match this next bitwise
                    "0x4",  # unmapped
                    temp_bam,
                ],
                stderr=None,
                check=True,
                capture_output=False,
            )
            self.index_bam(bamfile)

            bam_files.append(Path(bamfile))
        if tmpdirhandle is not None:
            tmpdirhandle.cleanup()
        return bam_files

    def index_bam(self, bamfile: Path) -> bool:
        """Index BAM files."""
        subprocess.run([self.samtools, "index", bamfile])
        # catch errors ...
        return True

    def _map_minimap_genome(
        self,
        genome: str,
        reads: ReadLibrary,
        outdir: str,
        tmpdir: str = "",
        thread_count: int = 1,
    ) -> str:
        """Align reads to genome with minimap2 - this may not be best tool so testing."""
        bamfile = os.path.join(outdir, reads.name + ".bam")
        bamfile

    def _map_bwa_genome(
        self, genome: str, reads: ReadLibrary, outdir: str, thread_count: int = 1
    ) -> str:
        """Align reads to genome with minimap2 - this may not be best tool so testing."""
        bamfile = os.path.join(outdir, reads.name + ".bam")
        bamfile
