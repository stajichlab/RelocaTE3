"""RelocaTE3 to genotype individuals for Transposon polymorphism.

Implements RelocaTE2 step 3: align reads to a transposon library, identify the
reads that contain TE sequence, and trim the TE portion off so that the
remaining *flanking* sequence can be aligned back to the genome (step 4) to
locate insertions (step 5).

The trim step emits, per read that overlaps a TE:

- a **flanking FASTQ** read whose name carries a junction suffix
  ``:start:5`` / ``:end:5`` / ``:start:3`` / ``:end:3`` (or ``:middle`` for a
  read fully inside a TE) encoding which read-end was trimmed and which TE end
  it abutted. These suffixes are what the insertion finder (step 5) keys on.
- a line in the **read_repeat_name** table mapping the *original* (unsuffixed)
  read name to the TE family and strand. Step 5 strips the junction suffix to
  recover this original name, so the table must use the unsuffixed name.
- the **TE-matched portion** as a FASTA entry (five_prime / three_prime).
- the full read in the **ContainingReads** FASTQ.

Sequence/orientation handling mirrors RelocaTE2's ``parse_align_bwa`` exactly:
read-relative ``start``/``end`` come from the aligned orientation
(``query_alignment_start``/``end``) while the sequence is taken in the original
read orientation (the reverse complement of the BAM SEQ for reverse alignments).
Reconstructing both from the BAM alone avoids re-reading the FASTQ and keeps the
trim output consistent with the step-5 ``TSD_check`` port.
"""

from __future__ import annotations

import os
from collections import defaultdict
from pathlib import Path

import pysam

from RelocaTE3 import logger
from RelocaTE3.align import Aligner
from RelocaTE3.ReadLibrary import ReadLibrary, TrimmedReadLibrary

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


class RelocaTE:
    """Process reads to identify transposon insertion and excision sites."""

    cpu_threads = 1  # number of CPU threads to use
    transposon_library = None
    verbose = 0

    def __init__(self, TElib: str = "", threads: int = 1, verbose: int = 0):
        """Initialize the RelocaTE object."""
        self.transposon_library = TElib
        self.cpu_threads = threads
        self.verbose = verbose

    # ------------------------------------------------------------------
    # orchestration
    # ------------------------------------------------------------------
    def identify_TE_reads(
        self,
        seqreads: ReadLibrary,
        outdir: Path,
        TE_library: str = "",
        search_tool: str = "minimap2",
        minimum_match_length: int = 10,
        minimum_trimmed_length: int = 10,
        mismatch_allowance: int = 0,
    ) -> int:
        """Align reads to the TE library, then trim and write flanking reads.

        Returns:
            int: number of flanking (trimmed) reads written.
        """
        if not TE_library:
            TE_library = self.transposon_library
        if "minimap" not in search_tool.lower():
            raise NotImplementedError(f"search_tool {search_tool!r} is not supported")

        alntool = Aligner(self.cpu_threads)
        alntool.index_minimap(TE_library)
        bamfiles = alntool.map_minimap_library(seqreads, outdir, TE_library)

        directions = self._bam_directions(bamfiles)
        return self.write_trimmed_reads(
            seqreads.name,
            list(zip(directions, bamfiles)),
            outdir,
            minimum_match_length=minimum_match_length,
            minimum_trimmed_length=minimum_trimmed_length,
            mismatch_allowance=mismatch_allowance,
        )

    @staticmethod
    def _bam_directions(bamfiles: list[Path]) -> list[str]:
        """Infer the read direction (left/right/...) from each BAM filename."""
        directions = []
        for i, bam in enumerate(bamfiles):
            stem = Path(bam).name
            if ".left." in stem:
                directions.append("left")
            elif ".right." in stem:
                directions.append("right")
            else:
                directions.append(f"set{i}")
        return directions

    # ------------------------------------------------------------------
    # writing trimmed / flanking reads (the core of step 3)
    # ------------------------------------------------------------------
    def write_trimmed_reads(
        self,
        name: str,
        direction_bams: list[tuple[str, Path]],
        outdir: Path,
        minimum_match_length: int = 10,
        minimum_trimmed_length: int = 10,
        mismatch_allowance: int = 0,
    ) -> int:
        """Trim TE sequence from reads and write the step-3 output files.

        Args:
            name: sample/individual name (output prefix).
            direction_bams: list of ``(direction, bam_path)`` pairs; each BAM is a
                TE-library alignment of one read file.
            outdir: base output directory.
            minimum_match_length: minimum TE match length (RelocaTE2 ``len_cut_match``).
            minimum_trimmed_length: minimum retained flank length (``len_cut_trim``).
            mismatch_allowance: maximum mismatches in the TE alignment.

        Returns:
            Number of flanking (trimmed) reads written across all directions.
        """
        outdir = Path(outdir)
        flank_dir = outdir / "flanking"
        contain_dir = outdir / "te_containing"
        portion_dir = outdir / "te_portions"
        for d in (flank_dir, contain_dir, portion_dir):
            d.mkdir(parents=True, exist_ok=True)

        read_repeat_path = contain_dir / f"{name}.read_repeat_name.txt"
        five_path = portion_dir / f"{name}.five_prime.fa"
        three_path = portion_dir / f"{name}.three_prime.fa"

        flank_written = 0
        # read_repeat / TE-portion files aggregate across directions (append)
        with open(read_repeat_path, "w") as rr_out, open(
            five_path, "w"
        ) as te5_out, open(three_path, "w") as te3_out:
            for direction, bam in direction_bams:
                coord = self._parse_te_bam(
                    Path(bam), mismatch_allowance=mismatch_allowance
                )
                flank_path = flank_dir / f"{name}.{direction}.flankingReads.fq"
                contain_path = contain_dir / f"{name}.{direction}.ContainingReads.fq"
                flank_written += self._write_direction(
                    coord,
                    flank_path,
                    contain_path,
                    rr_out,
                    te5_out,
                    te3_out,
                    minimum_match_length,
                    minimum_trimmed_length,
                    mismatch_allowance,
                )
        logger.info(
            "Wrote %d flanking reads; read_repeat table at %s",
            flank_written,
            read_repeat_path,
        )
        return flank_written

    def _write_direction(
        self,
        coord,
        flank_path,
        contain_path,
        rr_out,
        te5_out,
        te3_out,
        len_cutoff_m,
        len_cutoff_l,
        mismatch_allowance,
    ) -> int:
        """Apply the trim branch logic for one read file and write its outputs."""
        flank_written = 0
        with open(flank_path, "w") as flank_out, open(contain_path, "w") as contain_out:
            for rl_name, rec in coord.items():
                trimmed = self._trim_record(
                    rl_name,
                    rec,
                    rr_out,
                    te5_out,
                    te3_out,
                    len_cutoff_m,
                    len_cutoff_l,
                    mismatch_allowance,
                )
                if trimmed is None:
                    continue
                header, trimmed_seq, trimmed_qual = trimmed
                if len(trimmed_seq) >= len_cutoff_l:
                    flank_out.write(f"@{header}\n{trimmed_seq}\n+\n{trimmed_qual}\n")
                    flank_written += 1
                # every read that overlapped a TE is written to ContainingReads
                contain_out.write(f"@{header}\n{rec['seq']}\n+\n{rec['qual']}\n")
        return flank_written

    def _trim_record(
        self,
        rl_name,
        rec,
        rr_out,
        te5_out,
        te3_out,
        len_cutoff_m,
        len_cutoff_l,
        mismatch_allowance,
    ):
        """Classify one TE-overlapping read and return ``(header, seq, qual)``.

        Faithful port of the three RelocaTE2 trim branches (5' end of TE,
        3' end of TE, fully internal). Returns ``None`` when the read does not
        meet any branch's criteria.
        """
        start = int(rec["start"])
        end = int(rec["end"])
        length = int(rec["len"])
        t_name = rec["tName"]
        t_start = int(rec["tStart"])
        t_end = int(rec["tEnd"])
        t_len = int(rec["tLen"])
        mismatch = int(rec["mismatch"])
        match = int(rec["match"])
        strand = rec["strand"]
        seq = rec["seq"]
        qual = rec["qual"]

        match_span = match + mismatch
        flank_len = length - match_span
        ends_align = start <= 2 or end >= length - 3
        passes = (
            match_span >= len_cutoff_m
            and flank_len >= len_cutoff_l
            and mismatch <= mismatch_allowance
        )

        header = None
        trimmed_seq = trimmed_qual = ""
        te_subseq = seq[start : end + 1]

        if t_start <= 2 and ends_align and passes:
            # read overlaps the 5' end of the TE
            if strand == "-":
                te_subseq = reverse_complement(te_subseq)
                trimmed_seq, trimmed_qual = seq[end + 1 :], qual[end + 1 :]
                header = f"{rl_name}:start:5"
            else:
                trimmed_seq, trimmed_qual = seq[0:start], qual[0:start]
                header = f"{rl_name}:end:5"
            if len(trimmed_seq) >= len_cutoff_l:
                rr_out.write(f"{rl_name}\t{t_name}\t{strand}\n")
                te5_out.write(
                    self._te_fasta(
                        header, start, end, t_name, t_start, t_end, mismatch, te_subseq
                    )
                )
        elif t_end >= (t_len - 3) and ends_align and passes:
            # read overlaps the 3' end of the TE
            if strand == "-":
                te_subseq = reverse_complement(te_subseq)
                trimmed_seq, trimmed_qual = seq[0:start], qual[0:start]
                header = f"{rl_name}:end:3"
            else:
                trimmed_seq, trimmed_qual = seq[end + 1 :], qual[end + 1 :]
                header = f"{rl_name}:start:3"
            if len(trimmed_seq) >= len_cutoff_l:
                rr_out.write(f"{rl_name}\t{t_name}\t{strand}\n")
                te3_out.write(
                    self._te_fasta(
                        header, start, end, t_name, t_start, t_end, mismatch, te_subseq
                    )
                )
        elif start <= 2 and end + 1 >= length - 3:
            # read lies fully inside the TE: keep it whole, label :middle
            trimmed_seq, trimmed_qual = seq, qual
            header = f"{rl_name}:middle"
            rr_out.write(f"{rl_name}\t{t_name}\t{strand}\n")

        if header is None:
            return None
        return header, trimmed_seq, trimmed_qual

    @staticmethod
    def _te_fasta(
        header, start, end, t_name, t_start, t_end, mismatch, te_subseq
    ) -> str:
        """Format a FASTA entry for the TE-matched portion of a read."""
        return (
            f">{header} {start + 1}..{end + 1} matches {t_name}:{t_start + 1}.."
            f"{t_end + 1} mismatches:{mismatch}\n{te_subseq}\n"
        )

    # ------------------------------------------------------------------
    # alignment parsing
    # ------------------------------------------------------------------
    def trim_TE_reads(
        self,
        reads: ReadLibrary,
        bamfiles: list[Path],
        minimum_match_length: int = 10,
        minimum_trimmed_length: int = 10,
        mismatch_allowance: int = 0,
    ) -> TrimmedReadLibrary:
        """Parse TE-library BAM(s) into best-match coordinates per read.

        Returns a :class:`TrimmedReadLibrary` whose ``trimmed_coordinates`` holds
        the merged best TE match for every read across the supplied BAMs.
        """
        trimlib = TrimmedReadLibrary(ind_name=reads.name)
        merged: dict = defaultdict(dict)
        for bam in bamfiles:
            coord = self._parse_te_bam(Path(bam), mismatch_allowance=mismatch_allowance)
            for qname, rec in coord.items():
                if qname not in merged or self._is_better(rec, merged[qname]):
                    merged[qname] = rec
        trimlib.trimmed_coordinates = merged
        return trimlib

    @staticmethod
    def _is_better(new_rec: dict, old_rec: dict) -> bool:
        """Best-match preference: more boundary contact, then longer match."""
        if int(new_rec["boundary"]) > int(old_rec["boundary"]):
            return True
        if int(new_rec["boundary"]) == int(old_rec["boundary"]):
            return int(new_rec["match"]) > int(old_rec["match"])
        return False

    def _parse_te_bam(self, bam: Path, mismatch_allowance: int = 0) -> dict:
        """Parse one TE-library BAM into ``{read_name: best-match record}``.

        Each record stores the read-relative match coordinates, the TE target
        coordinates, the mismatch/match counts, strand, a boundary score used to
        pick the best match, and the read sequence/qualities reconstructed in
        original (FASTQ) orientation.
        """
        coord: dict = defaultdict(dict)
        if not os.path.exists(f"{bam}.bai"):
            Aligner().index_bam(bam)
        fbam = pysam.AlignmentFile(str(bam), "rb")
        try:
            refnames = fbam.references
            reflengths = fbam.lengths
            qlen_c = 0
            for record in fbam.fetch(until_eof=True):
                if record.is_unmapped:
                    continue
                qname = record.query_name
                qlen = int(record.query_length)
                qstart = int(record.query_alignment_start)
                try:
                    qend = int(record.query_alignment_end) - 1
                except (ValueError, TypeError):
                    continue
                if qlen == 0:
                    qlen = qlen_c
                    qend = qlen + qend
                else:
                    qlen_c = qlen

                tname = refnames[record.reference_id]
                tlen = int(reflengths[record.reference_id])
                tstart = int(record.reference_start)
                tend = int(record.reference_end) - 1

                tags = self._convert_tag(record.tags if record.tags else [])
                match = ins0 = del0 = 0
                for key, length in record.cigartuples or []:
                    if key == 0:
                        match += length
                    elif key == 1:
                        ins0 += length
                    elif key == 2:
                        del0 += length
                mismatch = int(tags.get("NM", 0)) - ins0 - del0
                match = match - mismatch

                strand = (
                    "+" if record.flag == 0 else ("-" if record.is_reverse else "+")
                )

                boundary = (
                    (1 if qstart <= 2 else 0)
                    + (1 if tstart <= 2 else 0)
                    + (1 if qend + 1 == qlen or qend >= qlen - 3 else 0)
                    + (1 if tend + 1 == tlen or tend >= tlen - 3 else 0)
                )

                seq, qual = self._original_orientation(record)
                new_rec = {
                    "match": match,
                    "len": qlen,
                    "start": qstart,
                    "end": qend,
                    "tLen": tlen,
                    "mismatch": mismatch,
                    "strand": strand,
                    "tName": tname,
                    "tStart": tstart,
                    "tEnd": tend,
                    "boundary": boundary,
                    "seq": seq,
                    "qual": qual,
                }
                if qname not in coord or self._is_better(new_rec, coord[qname]):
                    coord[qname] = new_rec
        finally:
            fbam.close()
        return coord

    @staticmethod
    def _original_orientation(record) -> tuple[str, str]:
        """Reconstruct the read sequence and qualities in original FASTQ orientation."""
        seq = record.query_sequence or ""
        quals = record.query_qualities
        qual_str = (
            pysam.qualities_to_qualitystring(quals)
            if quals is not None
            else "I" * len(seq)
        )
        if record.is_reverse:
            seq = reverse_complement(seq)
            qual_str = qual_str[::-1]
        return seq, qual_str

    @staticmethod
    def _convert_tag(tag: list) -> dict:
        """Convert SAM tags list into a dict."""
        return {t[0]: t[1] for t in tag}
