"""Find non-reference transposon insertions from genome-aligned flanking reads (RelocaTE2 step 5).

This is a first-cut port of ``relocaTE_insertionFinder.py``. It implements the
core algorithm for the case where the TSD (target site duplication) motif is
known: trimmed flanking reads are aligned to the reference genome (step 4),
clustered by genomic position, and each cluster that has converging left/right
junctions defines a candidate insertion. The output ``*.all_nonref_insert.txt``
table is the input consumed by :mod:`RelocaTE3.characterize` (step 7).

Junction reads are recognised by the name suffix the trim step attaches:
``<read>:start:5`` / ``<read>:end:5`` / ``<read>:start:3`` / ``<read>:end:3``
where ``start``/``end`` is which end of the read the TE was trimmed from and
``5``/``3`` is which end of the TE the flank abutted.

Adaptation from RelocaTE2: the original quality filtering relied on BWA-specific
SAM tags (``XT``/``X1``/``XM``/``XO``). minimap2 does not emit those, so reads
are filtered here on MAPQ (uniqueness) and the ``NM`` edit distance minus
indels (mismatch count) instead.

Not yet ported (documented TODO): TSD-unknown inference from read depth
(``TSD_from_read_depth``), the supporting-read-only ``all_nonref_supporting``
output, and reference-insertion subtraction via bedtools.
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import pysam

from RelocaTE3 import logger
from RelocaTE3.reference_te import ReferenceTEAnnotator

# junction-read name suffix: "<read>:start|end:5|3"
_JUNCTION = re.compile(r"(.*):(start|end):([53])")
_RANGE_ALLOWANCE = 1000  # max gap (bp) before a read starts a new cluster


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


class InsertionFinder:
    """Cluster genome-aligned flanking reads into non-reference insertion calls."""

    def __init__(self, mismatch_allow: int = 0, min_mapq: int = 1, verbose: int = 0):
        """Initialize the finder.

        Args:
            mismatch_allow: maximum read/genome mismatches (excluding indels).
            min_mapq: minimum MAPQ for a read to be considered uniquely mapped.
            verbose: verbosity level.
        """
        self.mismatch_allow = mismatch_allow
        self.min_mapq = min_mapq
        self.verbose = verbose

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------
    def find_insertions(
        self,
        bam_file: Path,
        read_repeat_file: Path,
        tsd: str,
        target: str,
        sample: str,
        outdir: Path,
        te_name: str = "repeat",
        reference_ins: Path | None = None,
    ) -> Path:
        """Find non-reference insertions and write the ``all_nonref_insert`` table.

        Args:
            bam_file: coordinate-sorted BAM of flanking reads aligned to the genome.
            read_repeat_file: ``read_repeat_name`` table mapping read -> (TE, strand).
            tsd: the TSD motif (e.g. ``"TTA"``). TSD-unknown mode is not yet supported.
            target: chromosome to analyze, or ``"ALL"``.
            sample: sample/experiment name (the ``exper`` column).
            outdir: directory for the ``results/`` output.
            te_name: TE label used in output filenames.
            reference_ins: optional RepeatMasker ``.out``/BED of existing copies to skip.

        Returns:
            Path to the written ``*.all_nonref_insert.txt`` file.
        """
        if re.search(r"UNK|UKN|unknown", tsd, re.IGNORECASE):
            raise NotImplementedError(
                "TSD-unknown (read-depth) inference is not yet ported; provide a TSD motif."
            )

        read_repeat = self._load_read_repeat(read_repeat_file)
        existing_te = (
            ReferenceTEAnnotator.load_existing_te(reference_ins, target)
            if reference_ins
            else defaultdict(lambda: {"start": {}, "end": {}})
        )

        # teInsertions[event][tsd_start][tsd_seq] -> counts; reads kept in parallel
        te_insertions: dict = defaultdict(
            lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        )
        te_insertions_reads: dict = defaultdict(
            lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        )
        cluster_chrom: dict[int, str] = {}

        self._cluster_reads(
            bam_file,
            target,
            tsd,
            read_repeat,
            existing_te,
            te_insertions,
            te_insertions_reads,
            cluster_chrom,
        )

        result_dir = Path(outdir) / "results"
        result_dir.mkdir(parents=True, exist_ok=True)
        out_txt = result_dir / f"{target}.{te_name}.all_nonref_insert.txt"
        self._write_output(
            out_txt,
            sample,
            target,
            read_repeat,
            te_insertions,
            te_insertions_reads,
            cluster_chrom,
        )
        return out_txt

    # ------------------------------------------------------------------
    # inputs
    # ------------------------------------------------------------------
    @staticmethod
    def _load_read_repeat(read_repeat_file: Path) -> dict:
        """Load read -> [TE_name, strand] from the trim step's mapping table."""
        data: dict[str, list[str]] = {}
        path = Path(read_repeat_file)
        if not path.exists():
            logger.warning("read_repeat_name file not found: %s", path)
            return data
        with open(path) as handle:
            for line in handle:
                line = line.rstrip()
                if len(line) <= 2:
                    continue
                unit = line.split("\t")
                # name -> [TE_name, strand]; extra columns ignored
                data[unit[0]] = (
                    unit[1:3]
                    if len(unit) >= 3
                    else [unit[1] if len(unit) > 1 else "NA", ""]
                )
        return data

    # ------------------------------------------------------------------
    # clustering
    # ------------------------------------------------------------------
    def _cluster_reads(
        self,
        bam_file,
        target,
        tsd,
        read_repeat,
        existing_te,
        te_insertions,
        te_insertions_reads,
        cluster_chrom,
    ):
        """Walk the BAM, group reads into positional clusters and score junctions."""
        ref = None if target == "ALL" else target
        bam = pysam.AlignmentFile(str(bam_file), "rb")
        try:
            rnames = bam.references
            bin_ins = [0]
            count = 0
            for record in bam.fetch(reference=ref, until_eof=True):
                if record.is_unmapped:
                    continue
                if not self._passes_quality(record):
                    continue
                name = record.query_name
                start = int(record.reference_start) + 1  # 1-based
                end = int(record.reference_end) + 1
                seq = record.query_sequence or ""
                chro = rnames[record.reference_id]
                strand = "-" if record.is_reverse else "+"

                bin_ins, count = self._assign_cluster(
                    bin_ins,
                    count,
                    start,
                    end,
                    name,
                    seq,
                    chro,
                    strand,
                    tsd,
                    read_repeat,
                    existing_te,
                    te_insertions,
                    te_insertions_reads,
                )
                cluster_chrom[count] = chro
        finally:
            bam.close()

    def _passes_quality(self, record) -> bool:
        """Minimap2-adapted quality filter (replaces BWA XT/X1/XM/XO logic)."""
        if record.mapping_quality < self.min_mapq:
            return False
        mismatch = self._mismatch_count(record)
        return mismatch is None or mismatch <= self.mismatch_allow

    @staticmethod
    def _mismatch_count(record):
        """Edit distance (NM) minus indel bases, mirroring librelocate's logic."""
        if not record.has_tag("NM"):
            return None
        ins = sum(length for op, length in (record.cigartuples or []) if op == 1)
        dele = sum(length for op, length in (record.cigartuples or []) if op == 2)
        return int(record.get_tag("NM")) - ins - dele

    def _assign_cluster(
        self,
        bin_ins,
        count,
        start,
        end,
        name,
        seq,
        chro,
        strand,
        tsd,
        read_repeat,
        existing_te,
        te_insertions,
        te_insertions_reads,
    ):
        """Place a read in the current cluster or open a new one, then score it."""
        padded_start = bin_ins[0] - _RANGE_ALLOWANCE
        padded_end = bin_ins[-1] + _RANGE_ALLOWANCE
        in_range = (padded_start <= start <= padded_end) or (
            padded_start <= end <= padded_end
        )
        if in_range:
            bin_ins.extend([start, end])
            bin_ins.sort()
        else:
            count += 1
            bin_ins = [start, end]

        match = _JUNCTION.search(name)
        if match:
            real_name = match.group(1)
            self._tsd_check(
                count,
                seq,
                chro,
                start,
                end,
                real_name,
                read_repeat,
                name,
                tsd,
                strand,
                existing_te,
                te_insertions,
                te_insertions_reads,
            )
        return bin_ins, count

    def _tsd_check(
        self,
        event,
        seq,
        chro,
        start,
        end,
        real_name,
        read_repeat,
        name,
        tsd,
        strand,
        existing_te,
        te_insertions,
        te_insertions_reads,
    ):
        """Faithful port of RelocaTE2 ``TSD_check`` (known-TSD path).

        Determines, for a junction read, which boundary (left/right) it marks,
        the TE orientation, the TSD position and sequence, then records the
        junction unless it lands on a known reference-TE edge.
        """
        rev_com = _reverse_complement(seq)
        r5 = re.compile(r"start:[53]$")
        r3 = re.compile(r"end:[53]$")
        r5_tsd = re.compile(rf"^({tsd})")
        r3_tsd = re.compile(rf"({tsd})$")

        result = 0
        pos = ""
        te_orient = 0
        tsd_start = 0
        tsd_seq = ""

        # start: TE trimmed from start of read; 5/3: which TE end the flank abuts
        if strand == "+":
            if r5.search(name) and (r5_tsd.search(seq) or r3_tsd.search(rev_com)):
                result = 1
                m = r5_tsd.search(seq)
                tsd_seq = m.group(1) if m else "UNK"
                pos = "right"
                te_orient = "-" if name[-1] == "5" else "+"
                tsd_start = start
            elif r3.search(name) and (r5_tsd.search(rev_com) or r3_tsd.search(seq)):
                result = 1
                m = r3_tsd.search(seq)
                tsd_seq = m.group(1) if m else "UNK"
                pos = "left"
                te_orient = "+" if name[-1] == "5" else "-"
                tsd_start = end - len(tsd)
        elif strand == "-":
            if r5.search(name) and (r5_tsd.search(rev_com) or r3_tsd.search(seq)):
                result = 1
                m = r3_tsd.search(seq)
                tsd_seq = m.group(1) if m else "UNK"
                pos = "left"
                te_orient = "+" if name[-1] == "5" else "-"
                tsd_start = end - len(tsd)
            elif r3.search(name) and (r5_tsd.search(seq) or r3_tsd.search(rev_com)):
                result = 1
                m = r5_tsd.search(seq)
                tsd_seq = m.group(1) if m else "UNK"
                pos = "right"
                te_orient = "-" if name[-1] == "5" else "+"
                tsd_start = start

        if not (result and te_orient):
            return

        tir1_end = end if pos == "left" else 0
        tir2_end = (start - 1) if pos == "right" else 0
        # skip junctions that land on a known reference-TE boundary
        if tir1_end > 0 and tir1_end in existing_te[chro]["start"]:
            return
        if tir2_end > 0 and tir2_end in existing_te[chro]["end"]:
            return

        bucket = te_insertions[event][tsd_start][tsd_seq]
        bucket["count"] += 1
        bucket[pos] += 1
        bucket[te_orient] += 1
        te_insertions_reads[event][tsd_start][tsd_seq]["read"].append(name)
        if pos == "left":
            te_insertions_reads[event][tsd_start][tsd_seq]["left_read"].append(name)
        else:
            te_insertions_reads[event][tsd_start][tsd_seq]["right_read"].append(name)

    # ------------------------------------------------------------------
    # output
    # ------------------------------------------------------------------
    def _write_output(
        self,
        out_txt,
        sample,
        target,
        read_repeat,
        te_insertions,
        te_insertions_reads,
        cluster_chrom,
    ):
        """Write the ``all_nonref_insert`` table consumed by characterize (step 7)."""
        with open(out_txt, "w") as out:
            for event in sorted(te_insertions, key=int):
                for tsd_start in sorted(te_insertions[event], key=int):
                    self._write_event_start(
                        out,
                        event,
                        tsd_start,
                        sample,
                        target,
                        read_repeat,
                        te_insertions,
                        te_insertions_reads,
                    )
        logger.info("Wrote insertions table %s", out_txt)

    def _write_event_start(
        self,
        out,
        event,
        tsd_start,
        sample,
        target,
        read_repeat,
        te_insertions,
        te_insertions_reads,
    ):
        """Emit one insertion row, picking the dominant TSD and TE orientation."""
        total_count = left_count = right_count = 0
        fwd = rev = 0
        tsd_count: dict[str, int] = {}
        reads: list[str] = []
        for found_tsd in sorted(te_insertions[event][tsd_start]):
            b = te_insertions[event][tsd_start][found_tsd]
            total_count += b["count"]
            left_count += b["left"]
            right_count += b["right"]
            fwd += b["+"]
            rev += b["-"]
            tsd_count[found_tsd] = b["count"]
            reads.extend(te_insertions_reads[event][tsd_start][found_tsd]["read"])

        if not tsd_count:
            return
        top_tsd = max(tsd_count.items(), key=lambda kv: kv[1])[0]
        te_orient = "+" if fwd > rev else "-"
        repeat_family = self._insertion_family(reads, read_repeat)

        # coordinate range: TSD spans [tsd_start, tsd_start + len(top_tsd) - 1]
        coor_start = tsd_start
        coor = tsd_start + max(len(top_tsd) - 1, 0)

        # status mirrors RelocaTE2: a true junction needs both left and right reads
        if left_count > 0 and right_count > 0:
            tsd_field = top_tsd
        elif total_count == 1:
            tsd_field = "singleton"
        else:
            tsd_field = "supporting_junction"

        out.write(
            f"{repeat_family}\t{tsd_field}\t{sample}\t{target}\t{coor_start}..{coor}\t"
            f"{te_orient}\tT:{total_count}\tR:{right_count}\tL:{left_count}\t"
            f"ST:0\tSR:0\tSL:0\n"
        )

    @staticmethod
    def _insertion_family(reads, read_repeat) -> str:
        """Pick the dominant TE family among a cluster's junction reads."""
        family: dict[str, int] = defaultdict(int)
        for read in reads:
            m = _JUNCTION.search(read)
            real = m.group(1) if m else None
            if real and real in read_repeat:
                family[read_repeat[real][0]] += 1
        if not family:
            return ""
        return max(family.items(), key=lambda kv: kv[1])[0]
