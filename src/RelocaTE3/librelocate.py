"""RelocaTE3 to genotype individuals for Transposon polymorphism."""

from __future__ import annotations

import os
import sys
from collections import defaultdict
from pathlib import Path

import pysam

from RelocaTE3.align import Aligner
from RelocaTE3.ReadLibrary import ReadLibrary, TrimmedReadLibrary

# from Bio import SeqIO


class RelocaTE:
    """Process reads and mapping to identify transposon insertion and excision sites."""

    cpu_threads = 1  # number of CPU threads to use
    transposon_library = None
    verbose = 0

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
    ) -> int:
        """Search for sequence reads containing transposon sequences.

        Args:
            transposon_library: Fasta file of sequences that define transposon sequences to find in read file.

        Returns:
            int: number of reads (or read pairs) identified as containing transposon sequences.
        """
        if TE_library == "":
            TE_library = self.transposon_library

        alntool = Aligner(self.cpu_threads)
        if "minimap" in search_tool.lower():
            alntool.index_minimap(TE_library)
            bamfiles = alntool.map_minimap_library(seqreads, outdir, TE_library)
            TE_to_readinfo = self.trim_TE_reads(seqreads, bamfiles)
            print(TE_to_readinfo)

    def trim_TE_reads(
        self,
        reads: ReadLibrary,
        bamfiles: Path[list],
        minimum_match_length: int = 10,
        minimum_trimmed_length: int = 10,
        mismatch_allowance: int = 0,
    ) -> TrimmedReadLibrary:
        """Function to process transposons aligned to reads."""
        trimlib = TrimmedReadLibrary(reads.name)
        coord = defaultdict(lambda: defaultdict(lambda: str))
        # we need to write code here to figure out the TEs
        # open BAM file with pysam
        # for each read, identify reference (eg mPing or Ping etc) aligned to
        # keep track of orientation and start/end of alignment to the transposable element seq
        for bam in bamfiles:
            if not os.path.exists(str(bam) + ".bai"):
                alntool = Aligner()
                alntool.index_bam(bam)
            fbam = pysam.AlignmentFile(bam, "rb")
            refnames = fbam.references
            reflengths = fbam.lengths
            qlen_c = 0
            for TE_ref in refnames:
                # each TE in the DB matched would be a sep entry / query / analysis
                query = fbam.fetch(reference=TE_ref, until_eof=True)
                for record in query:
                    qName = record.query_name
                    qLen = int(record.query_length)
                    qStart = int(record.query_alignment_start)
                    qEnd = 0
                    try:
                        qEnd = int(record.query_alignment_end) - 1
                    except ValueError:
                        continue
                    if qLen == 0:
                        qLen = qlen_c
                        qEnd = qLen + qEnd
                    else:
                        qlen_c = qLen

                    tName = refnames[record.reference_id]
                    tLen = int(reflengths[record.reference_id])
                    tStart = int(record.reference_start)
                    tEnd = int(record.reference_end) - 1
                    # match and mismatch
                    tag = record.tags if record.tags else []
                    tags = self._convert_tag(tag)
                    match = 0
                    ins0 = 0
                    del0 = 0
                    for key, length in record.cigartuples:
                        # print key, length
                        if int(key) == 0:
                            match += length
                        elif int(key) == 1:
                            ins0 += length
                        elif int(key) == 2:
                            del0 += length
                    mismatch = int(tags["NM"]) - int(ins0) - int(del0)
                    match = match - mismatch
                    # strand, flag is 0 is read if read is unpaired and mapped to plus strand
                    strand = ""
                    flag = record.flag
                    if int(flag) == 0:
                        strand = "+"
                    else:
                        strand = "-" if record.is_reverse else "+"
                    # update data
                    # boundary = 1 if int(qStart) == 0 or int(qEnd) + 1 == int(qLen) else 0
                    addRecord = 0
                    boundary = 0
                    boundary_qry_left = 0
                    boundary_tar_left = 0
                    boundary_qry_right = 0
                    boundary_tar_right = 0
                    if int(qStart) == 0 or int(qStart) <= 2:
                        boundary_qry_left = 1
                    if int(qEnd) + 1 == int(qLen) or int(qEnd) >= int(qLen) - 3:
                        boundary_qry_right = 1
                    if int(tStart) == 0 or int(tStart) <= 2:
                        boundary_tar_left = 1
                    if int(tEnd) + 1 == int(tLen) or int(tEnd) >= int(tLen) - 3:
                        boundary_tar_right = 1
                    # max boundary should be 2: 1. match one read end and one repeat end; 2. match two read end and internal of repeat
                    # we expect more boundary and compare match length when having equal number of boundary
                    boundary = (
                        boundary_qry_left
                        + boundary_tar_left
                        + boundary_qry_right
                        + boundary_tar_right
                    )
                    # change to use warnings?
                    if self.verbose >= 3:
                        print(
                            qName,
                            qLen,
                            qStart,
                            qEnd,
                            tName,
                            tLen,
                            tStart,
                            tEnd,
                            match,
                            mismatch,
                            boundary,
                            file=sys.stderr,
                        )
                    if qName in coord:
                        # keep the best match to TE
                        if int(boundary) > int(coord[qName]["boundary"]):
                            addRecord = 1
                        elif int(boundary) == int(coord[qName]["boundary"]):
                            if int(match) > int(coord[qName]["match"]):
                                addRecord = 1
                            else:
                                addRecord = 0
                        else:
                            addRecord = 0
                    else:
                        addRecord = 1
                    if addRecord == 1:
                        coord[qName]["match"] = match
                        coord[qName]["len"] = qLen
                        coord[qName]["start"] = qStart
                        coord[qName]["end"] = qEnd
                        coord[qName]["tLen"] = tLen
                        coord[qName]["mismatch"] = mismatch
                        coord[qName]["strand"] = strand
                        coord[qName]["tName"] = tName
                        coord[qName]["tStart"] = tStart
                        coord[qName]["tEnd"] = tEnd
                        coord[qName]["boundary"] = boundary
            fbam.close()
        trimlib.trimmed_coordinates(coord)
        return trimlib

    def _convert_tag(self, tag: list):
        """Convert SAM tags."""
        tags = {}
        for t in tag:
            tags[t[0]] = t[1]
        return tags

    def _update_coord(self, header1, header, coord):
        coord[header]["start"] = int(coord[header1]["start"])
        coord[header]["len"] = int(coord[header1]["len"])
        coord[header]["end"] = int(coord[header1]["end"])
        coord[header]["tName"] = coord[header1]["tName"]
        coord[header]["tStart"] = int(coord[header1]["tStart"])
        coord[header]["tEnd"] = int(coord[header1]["tEnd"])
        coord[header]["tLen"] = int(coord[header1]["tLen"])
        coord[header]["mismatch"] = int(coord[header1]["mismatch"])
        coord[header]["match"] = int(coord[header1]["match"])
        coord[header]["strand"] = coord[header1]["strand"]
        coord[header]["boundary"] = int(coord[header1]["boundary"])

        del coord[header1]["start"]
        del coord[header1]["len"]
        del coord[header1]["end"]
        del coord[header1]["tName"]
        del coord[header1]["tStart"]
        del coord[header1]["tEnd"]
        del coord[header1]["tLen"]
        del coord[header1]["mismatch"]
        del coord[header1]["match"]
        del coord[header1]["strand"]
        del coord[header1]["boundary"]
