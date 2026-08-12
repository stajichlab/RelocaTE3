"""Step 5: call non-reference insertions from genome-aligned flanking reads.

Includes both the original class-based port and the newer function-based helpers
used by the pipeline code.
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import pysam

from RelocaTE3 import logger

# junction-read name suffix: "<read>:start|end:5|3"
_JUNCTION = re.compile(r"(.*):(start|end):([53])")
_RANGE_ALLOWANCE = 1000  # max gap (bp) before a read starts a new cluster


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


class InsertionFinder:
    """Cluster genome-aligned flanking reads into non-reference insertion calls."""

    def __init__(
        self,
        mismatch_allow: int = 0,
        min_mapq: int = 1,
        verbose: int = 0,
        require_both_junctions: bool = False,
        insert_size: int = 500,
    ):
        """Initialize the finder.

        Args:
            mismatch_allow: maximum read/genome mismatches (excluding indels).
            min_mapq: minimum MAPQ for a read to be considered uniquely mapped.
            verbose: verbosity level.
            require_both_junctions: when True, emit only insertions supported by
                both a left and a right junction read (drop single-sided calls),
                matching RelocaTE2, which drops one-sided ``insufficient_data``/
                ``singleton`` sites. Single-sided calls cannot resolve a TSD (they
                report ``UNK``) and are the bulk of RelocaTE3's false positives.
        """
        self.mismatch_allow = mismatch_allow
        self.min_mapq = min_mapq
        self.verbose = verbose
        self.require_both_junctions = require_both_junctions
        self.insert_size = insert_size

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
            return self._find_insertions_unknown_tsd(
                bam_file, read_repeat_file, target, sample, outdir, te_name,
                reference_ins,
            )

        read_repeat = self._load_read_repeat(read_repeat_file)
        if reference_ins:
            from RelocaTE3.reference_te import ReferenceTEAnnotator

            existing_te = ReferenceTEAnnotator.load_existing_te(reference_ins, target)
        else:
            existing_te = defaultdict(lambda: {"start": {}, "end": {}})

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

        # Collapse tsd_start sub-buckets split by the fixed-length wildcard TSD so
        # one insertion is not emitted as two adjacent single-sided calls.
        self._merge_offset_starts(te_insertions, te_insertions_reads)

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

    def _find_insertions_unknown_tsd(
        self,
        bam_file: Path,
        read_repeat_file: Path,
        target: str,
        sample: str,
        outdir: Path,
        te_name: str,
        reference_ins: Path | None,
    ) -> Path:
        """Call variable-length TSDs using the RelocaTE2 ``UNK`` strategy.

        The live CLI historically used only the fixed-pattern class path.  The
        function-based path below preserves its output contract while using the
        already ported breakpoint/depth inference helpers for each cluster.
        """
        read_repeat = self._load_read_repeat(read_repeat_file)
        if reference_ins:
            from RelocaTE3.reference_te import ReferenceTEAnnotator

            existing_te = ReferenceTEAnnotator.load_existing_te(reference_ins, target)
        else:
            existing_te = defaultdict(lambda: {"start": {}, "end": {}})

        result_dir = Path(outdir) / "results"
        result_dir.mkdir(parents=True, exist_ok=True)
        out_txt = result_dir / f"{target}.{te_name}.all_nonref_insert.txt"
        support_only: list[Insertion] = []
        with open(out_txt, "w") as out:
            for cluster in _stream_clusters(
                str(bam_file), read_repeat, quality_filter=self._passes_quality
            ):
                if target != "ALL" and cluster.chrom != target:
                    continue
                wrote_any = False
                for ins in _call_insertions(cluster, genome=None):
                    # RelocaTE2 parity: drop single-sided (one-junction) calls,
                    # which cannot resolve a TSD and dominate the false positives.
                    if self.require_both_junctions and (
                        ins.left_junction_reads == 0
                        or ins.right_junction_reads == 0
                    ):
                        continue
                    # Match the fixed-TSD path's reference-TE edge exclusion.
                    edges = existing_te[cluster.chrom]
                    if ins.end in edges["start"] or ins.start - 1 in edges["end"]:
                        continue
                    out.write(
                        f"{ins.te_name}\t{ins.tsd}\t{sample}\t{ins.chrom}\t"
                        f"{ins.start}..{ins.end}\t{ins.strand}\t"
                        f"T:{ins.left_junction_reads + ins.right_junction_reads}\t"
                        f"R:{ins.right_junction_reads}\tL:{ins.left_junction_reads}\t"
                        f"ST:0\tSR:{ins.right_support_reads}\tSL:{ins.left_support_reads}\n"
                    )
                    wrote_any = True
                # RelocaTE2 calls a site from mates alone only when no junction
                # read produced one, and files it separately (NONSUP).
                if not wrote_any:
                    call = call_support_only(cluster, insert_size=self.insert_size)
                    if call is not None:
                        support_only.append(call)
        logger.info("Wrote variable-length TSD insertions table %s", out_txt)
        write_supporting_reads(result_dir, target, te_name, sample, support_only)
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
            prev_chro: str | None = None
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

                # Force a new cluster when the chromosome changes so that
                # cluster_chrom is never overwritten and clusters don't span
                # contigs that share a coordinate range.
                if prev_chro is not None and chro != prev_chro:
                    count += 1
                    bin_ins = [start, end]

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
                prev_chro = chro
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
    # sub-cluster reconciliation
    # ------------------------------------------------------------------
    def _merge_offset_starts(self, te_insertions, te_insertions_reads):
        """Merge tsd_start sub-buckets split by the fixed-length wildcard TSD.

        The right junction records ``tsd_start = start`` (the TSD left edge) while
        the left junction records ``tsd_start = end - len(tsd)``. When the true TSD
        is longer than the wildcard pattern the two sides diverge by
        ``true_len - wildcard_len`` (1-2 bp), fragmenting one insertion into two
        single-sided rows. Within each cluster, collapse tsd_start entries whose
        TSD coordinate spans overlap (gap < captured TSD length) into a single
        canonical start, pooling counts and reads. Distinct insertions have
        non-overlapping TSDs (at least a TSD length apart) and are never merged.

        Mutates ``te_insertions`` and ``te_insertions_reads`` in place.
        """
        for event in list(te_insertions.keys()):
            starts = te_insertions[event]
            if len(starts) < 2:
                continue
            # Captured TSD length for this cluster (wildcard captures are uniform);
            # a length < 2 leaves nothing that could overlap, so skip.
            tsd_len = 0
            for pos_map in starts.values():
                for seq in pos_map:
                    tsd_len = max(tsd_len, len(seq))
            if tsd_len < 2:
                continue

            totals = {st: sum(b["count"] for b in starts[st].values()) for st in starts}

            # Chain consecutive starts whose spans overlap (gap < tsd_len).
            ordered = sorted(starts.keys(), key=int)
            groups: list[list[int]] = []
            current = [ordered[0]]
            for st in ordered[1:]:
                if int(st) - int(current[-1]) < tsd_len:
                    current.append(st)
                else:
                    groups.append(current)
                    current = [st]
            groups.append(current)

            for group in groups:
                if len(group) < 2:
                    continue
                # Canonical start = highest total count; tie -> smallest coordinate.
                canonical = group[0]
                best = (totals[canonical], -int(canonical))
                for st in group[1:]:
                    cand = (totals[st], -int(st))
                    if cand > best:
                        best = cand
                        canonical = st
                for st in group:
                    if st == canonical:
                        continue
                    for seq, bucket in te_insertions[event][st].items():
                        dst = te_insertions[event][canonical][seq]
                        for key, val in bucket.items():
                            dst[key] += val
                        src_reads = te_insertions_reads[event][st][seq]
                        dst_reads = te_insertions_reads[event][canonical][seq]
                        for key, val in src_reads.items():
                            dst_reads[key].extend(val)
                    del te_insertions[event][st]
                    if st in te_insertions_reads.get(event, {}):
                        del te_insertions_reads[event][st]

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
                chrom = cluster_chrom.get(event, target)
                for tsd_start in sorted(te_insertions[event], key=int):
                    self._write_event_start(
                        out,
                        event,
                        tsd_start,
                        sample,
                        chrom,
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
        chrom,
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

        # Emit the read-captured TSD whenever one was inferred. Mirrors
        # RelocaTE2's behavior post-TSD_from_read_depth: that pass synthesizes
        # the missing-side junction reads from the depth pileup, lifting
        # single-sided clusters into the both-sided emission. We achieve the
        # same emission directly when the wildcard TSD mode (validation
        # config tsd="...") has filled top_tsd with a real read-captured
        # 3-mer. Literal-TSD callers keep the original behavior because
        # top_tsd will be "UNK" when no read matched the motif.
        real_capture = bool(top_tsd) and top_tsd not in {"UNK", "UKN"}
        if real_capture and left_count > 0 and right_count > 0:
            tsd_field = top_tsd
        elif real_capture and total_count > 1:
            # single-sided junction with multiple reads; trust the wildcard capture
            tsd_field = top_tsd
        elif total_count == 1:
            tsd_field = "singleton"
        else:
            tsd_field = "supporting_junction"

        out.write(
            f"{repeat_family}\t{tsd_field}\t{sample}\t{chrom}\t{coor_start}..{coor}\t"
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


from RelocaTE3.models import Insertion, JunctionObservation

# read-name junction tag: <name>:(start|end):(5|3)
_JUNCTION_RE = re.compile(r":(start|end):([53])$")
# how far apart (bp) reads may be and still belong to one insertion cluster
RANGE_ALLOWANCE = 1000
# max gap (bp) between breakpoints of the SAME insertion. Above a TSD (+ read
# jitter), below the spacing between distinct insertions, so nearby breakpoints
# group into one sub-insertion while distinct sites stay separate.
SUBCLUSTER_GAP = 25
# a full read extending this far past a breakpoint indicates no insertion
FULLREAD_EXTEND = 10
# minimum bracketing reads per strand for a support-only (no junction) call
MIN_SUPPORT_ONLY = 2


def _junction_info(
    name: str, strand: str, gstart: int, gend: int
) -> tuple[str, int, str] | None:
    """Return (side, breakpoint, te_end) for a junction read, or None.

    The flank's TE-adjacent edge is the breakpoint. A right-junction read's
    breakpoint is the genomic start (left edge of the TSD); a left-junction
    read's breakpoint is the genomic end (right edge of the TSD).
    """
    m = _JUNCTION_RE.search(name)
    if not m:
        return None
    flank_side, te_end = m.group(1), m.group(2)
    if strand == "+":
        side, pos = ("right", gstart) if flank_side == "start" else ("left", gend)
    else:  # '-' strand inverts which read end is TE-adjacent
        side, pos = ("left", gend) if flank_side == "start" else ("right", gstart)
    return side, pos, te_end


def _te_family(read_repeat: dict[str, tuple[str, str]], read_name: str) -> str:
    """Best-effort TE family name for a junction read."""
    # The genome-aligned flank carries the junction suffix, while the trim
    # step's read_repeat_name table is keyed by the original untagged name.
    real_name = _JUNCTION_RE.sub("", read_name)
    if real_name in read_repeat:
        return read_repeat[real_name][0]
    return "NA"


class _Cluster:
    """Reads grouped within ``RANGE_ALLOWANCE`` bp on one chromosome."""

    def __init__(self, chrom: str):
        self.chrom = chrom
        self.lo: int | None = None
        self.hi: int | None = None
        self.junctions: list[JunctionObservation] = []
        # supporting reads: (name, gstart, gend, strand, seq)
        self.support: list[tuple[str, int, int, str, str]] = []

    def in_range(self, gstart: int, gend: int) -> bool:
        """True if a read at [gstart, gend] belongs to this cluster."""
        if self.lo is None:
            return True
        return gstart <= self.hi + RANGE_ALLOWANCE and gend >= self.lo - RANGE_ALLOWANCE

    def extend(self, gstart: int, gend: int) -> None:
        """Grow the cluster's coordinate span to include [gstart, gend]."""
        self.lo = gstart if self.lo is None else min(self.lo, gstart)
        self.hi = gend if self.hi is None else max(self.hi, gend)


def _stream_clusters(
    bam_path: str,
    read_repeat: dict[str, tuple[str, str]],
    quality_filter=None,
):
    """Yield :class:`_Cluster` objects by streaming a coordinate-sorted BAM."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        current: _Cluster | None = None
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped:
                continue
            if quality_filter is not None and not quality_filter(rec):
                continue
            chrom = bam.get_reference_name(rec.reference_id)
            gstart = rec.reference_start + 1
            gend = (
                rec.reference_end
            )  # pysam end is 0-based exclusive == 1-based inclusive
            strand = "-" if rec.is_reverse else "+"
            name = rec.query_name
            seq = rec.query_sequence or ""

            if (
                current is None
                or current.chrom != chrom
                or not current.in_range(gstart, gend)
            ):
                if current is not None:
                    yield current
                current = _Cluster(chrom)

            current.extend(gstart, gend)
            info = _junction_info(name, strand, gstart, gend)
            if info is not None:
                side, pos, te_end = info
                current.junctions.append(
                    JunctionObservation(
                        name,
                        side,
                        pos,
                        strand,
                        _te_family(read_repeat, name),
                        te_end,
                        seq,
                        gstart,
                        gend,
                    )
                )
            else:
                current.support.append((name, gstart, gend, strand, seq))
        if current is not None:
            yield current


def _group_by_position(
    observations: list[JunctionObservation],
) -> dict[int, list[JunctionObservation]]:
    """Group junction observations by their breakpoint position."""
    grouped: dict[int, list[JunctionObservation]] = {}
    for obs in observations:
        grouped.setdefault(obs.position, []).append(obs)
    return grouped


def _pair_breakpoints(
    left: dict[int, list], right: dict[int, list]
) -> list[tuple[int | None, int | None]]:
    """Pair left/right breakpoints into sub-insertions (RelocaTE2-style).

    ``left``/``right`` map each breakpoint position to its supporting junction
    reads (so ``len`` is the read support). Breakpoints are grouped into
    sub-insertions by proximity (one insertion's TSD-adjacent breakpoints fall
    within ``SUBCLUSTER_GAP``; distinct insertions are farther apart), then within
    each sub-cluster the **most-supported** (dominant) left and right breakpoint is
    selected -- mirroring RelocaTE2, which keys reads by breakpoint and resolves
    the TSD from the dominant boundary. Each returned tuple is
    ``(left_position, right_position)``; either may be ``None`` for a one-sided
    sub-cluster.

    A group whose dominant breakpoints are ordered ``left < right`` cannot be one
    insertion (the TSD geometry requires ``left >= right``); it is two adjacent
    insertions caught in the same group, and is split back into two one-sided
    sub-insertions rather than merged into a call at neither position.

    This replaces greedy nearest-neighbour pairing within a 100 bp window, which
    paired a lone left with a distant lone right and reported the ~100 bp gap as a
    (runaway) TSD -- the dominant source of RelocaTE3's false positives.
    """
    if not left and not right:
        return []
    marks = sorted([(p, "L") for p in left] + [(p, "R") for p in right])
    groups: list[list[tuple[int, str]]] = [[]]
    for pos, side in marks:
        if groups[-1] and pos - groups[-1][-1][0] > SUBCLUSTER_GAP:
            groups.append([])
        groups[-1].append((pos, side))

    pairs: list[tuple[int | None, int | None]] = []
    for group in groups:
        lpos = [p for p, s in group if s == "L"]
        rpos = [p for p, s in group if s == "R"]
        # dominant = most reads; ties -> the position closer to the other boundary
        # (left is the right edge of the TSD, right is the left edge) for a
        # sensible small TSD span.
        lp = max(lpos, key=lambda p: (len(left[p]), p)) if lpos else None
        rp = max(rpos, key=lambda p: (len(right[p]), -p)) if rpos else None
        if lp is not None and rp is not None and lp < rp:
            # Geometrically impossible as a single insertion: a left-junction
            # read marks the TSD's right edge and a right-junction read its left
            # edge (see _make_insertion), so a genuine pair has lp >= rp. lp < rp
            # means the group swept up two adjacent insertions whose *inner*
            # breakpoints happen to sit within SUBCLUSTER_GAP; pairing them
            # yields a call at neither site and loses both. Emit each breakpoint
            # as its own sub-insertion instead.
            pairs.append((lp, None))
            pairs.append((None, rp))
            continue
        pairs.append((lp, rp))
    return pairs


def _resolve_tsd(
    left_reads: list[JunctionObservation],
    right_reads: list[JunctionObservation],
    chrom: str,
    i_start: int,
    i_end: int,
    tsd_len: int,
    genome: pysam.FastaFile | None,
) -> str:
    """Capture TSD literally from a junction read; fall back to the genome.

    Mirrors RelocaTE2's read-derived TSD reporting (TSD_check_cluster). The
    genome fetch is only used when no read has the bases (e.g. supporting-only
    insertions).
    """
    if tsd_len <= 0:
        return "UNK"
    for obs in right_reads:
        captured = _capture_tsd_from_read(obs.seq, "right", tsd_len)
        if captured:
            return captured
    for obs in left_reads:
        captured = _capture_tsd_from_read(obs.seq, "left", tsd_len)
        if captured:
            return captured
    if genome is None:
        return "UNK"
    fetched = _fetch_tsd(genome, chrom, i_start, i_end)
    return fetched or "UNK"


def _majority_te_name(te_names: list[str]) -> str:
    """Pick the TE family for a cluster by majority vote, breaking ties by name.

    Reads tagged ``NA`` (no family assignment) do not vote; an empty vote yields
    ``"NA"``.

    The tie-break is load-bearing. This was previously
    ``max(set(names), key=names.count)``, which returns whichever tied name the
    set iterator happened to yield first -- and CPython randomises string
    hashing per process, so the same cluster could be labelled ``RIRE3`` on one
    run and ``mGing`` on the next, with byte-identical positions and read
    counts. Ranking by (-count, name) makes the winner reproducible. Which name
    wins among equals is arbitrary; that it is always the same one is the point.
    """
    votes = [name for name in te_names if name != "NA"]
    if not votes:
        return "NA"
    return min(set(votes), key=lambda name: (-votes.count(name), name))


def _make_insertion(
    chrom: str,
    left_reads: list[JunctionObservation],
    right_reads: list[JunctionObservation],
    genome: pysam.FastaFile | None,
    cluster: "_Cluster",
) -> Insertion:
    """Build an :class:`Insertion` from the left/right junction reads of one site.

    TSD inference follows RelocaTE2's read-depth path: when both breakpoints
    exist, the TSD width is the breakpoint overlap; otherwise it's estimated
    from supporting-read depth pileups (``_estimate_tsd_length_from_depth``).
    The TSD bases are captured from a junction read (R2 parity) and only fall
    back to the reference genome when no read sequence is available.
    """
    junctions = left_reads + right_reads
    te_name = _majority_te_name([j.te_name for j in junctions])

    orients = [j.te_orientation for j in junctions]
    strand = "+" if orients.count("+") >= orients.count("-") else "-"

    # RelocaTE2 estimates TSD length from a depth pileup over the *junction*
    # reads of the cluster, dividing by their count
    # (relocaTE_insertionFinder.py:1069-1076 feeding tsd_finder at :843).
    # Junction reads all abut the same breakpoint, so the TSD is the run of
    # positions nearly all of them cover. Supporting mates are spread across the
    # library insert and never reach the cutoff -- using them returned 0 for 15
    # of the 16 Chr3 sites where RelocaTE2 resolves a TSD and we reported UNK.
    spans = [(j.gstart, j.gend) for j in cluster.junctions if j.gend >= j.gstart > 0]
    if not spans:  # support-only clusters keep the previous behaviour
        spans = [(s, e) for _n, s, e, _strand, _seq in cluster.support]

    if left_reads and right_reads:
        i_end = left_reads[0].position  # right edge of TSD
        i_start = right_reads[0].position  # left edge of TSD
        overlap = i_end - i_start + 1
        if overlap > 0:
            tsd_len = overlap
        else:
            tsd_len = _estimate_tsd_length_from_depth(spans, min(i_start, i_end))
            i_start = i_end = min(i_start, i_end)
    else:
        present = left_reads or right_reads
        bp = present[0].position
        tsd_len = _estimate_tsd_length_from_depth(spans, bp)
        if tsd_len > 0 and right_reads:
            i_start, i_end = bp, bp + tsd_len - 1
        elif tsd_len > 0:
            i_start, i_end = bp - tsd_len + 1, bp
        else:
            i_start = i_end = bp

    tsd = _resolve_tsd(left_reads, right_reads, chrom, i_start, i_end, tsd_len, genome)

    return Insertion(
        chrom=chrom,
        start=i_start,
        end=i_end,
        te_name=te_name,
        strand=strand,
        tsd=tsd,
        left_junction_reads=len(left_reads),
        right_junction_reads=len(right_reads),
        read_names=[j.read_name for j in junctions],
    )


def _call_insertions(
    cluster: _Cluster, genome: pysam.FastaFile | None
) -> list[Insertion]:
    """Split a cluster into one or more insertions by pairing breakpoints."""
    left = _group_by_position([j for j in cluster.junctions if j.side == "left"])
    right = _group_by_position([j for j in cluster.junctions if j.side == "right"])
    if not left and not right:
        return []

    insertions: list[Insertion] = []
    for lp, rp in _pair_breakpoints(left, right):
        left_reads = left.get(lp, []) if lp is not None else []
        right_reads = right.get(rp, []) if rp is not None else []
        ins = _make_insertion(cluster.chrom, left_reads, right_reads, genome, cluster)
        _count_support(ins, cluster)
        insertions.append(ins)
    return insertions


def _fetch_tsd(genome: pysam.FastaFile, chrom: str, start: int, end: int) -> str:
    """Fetch the genomic TSD sequence (1-based inclusive), '.'-padded on failure."""
    try:
        seq = genome.fetch(chrom, start - 1, end)
        return seq.upper() if seq else "." * (end - start + 1)
    except (KeyError, ValueError):
        return "." * (end - start + 1)


def _estimate_tsd_length_from_depth(
    spans: list[tuple[int, int]],
    breakpoint: int,
    thresholds: tuple[float, ...] = (1.0, 0.8, 0.6),
) -> int:
    """Estimate TSD length from read-depth overlap near ``breakpoint``.

    Port of RelocaTE2 ``tsd_finder`` (relocaTE_insertionFinder.py:843). Builds a
    per-base depth pileup from ``spans`` (1-based inclusive ``(start, end)``
    tuples), then for each fractional threshold (in order) counts the positions
    whose depth >= ``threshold * len(spans)``. Returns the first non-zero
    length, or 0 if none qualify. Positions need not be contiguous, matching
    ``tsd_finder``; and because the first non-zero result wins, a single base
    covered by every read yields 1 rather than falling to a looser threshold.

    ``spans`` must be the cluster's *junction*-read spans -- see the call site
    in ``_make_insertion``.

    The ``breakpoint`` argument is reserved for future locality refinement; the
    R2 reference implementation also passes a candidate position but does not
    use it to bound the depth window.
    """
    del breakpoint  # currently unused; matches R2 signature shape
    if not spans:
        return 0
    depth: dict[int, int] = {}
    for s, e in spans:
        for p in range(s, e + 1):
            depth[p] = depth.get(p, 0) + 1
    total = len(spans)
    for frac in thresholds:
        cutoff = frac * total
        length = sum(1 for d in depth.values() if d >= cutoff)
        if length:
            return length
    return 0


def _capture_tsd_from_read(seq: str, side: str, length: int) -> str:
    """Return the literal TSD characters from a junction read.

    Mirrors RelocaTE2 ``TSD_check_cluster`` (relocaTE_insertionFinder.py:1249):
    a *right*-side junction read carries the TSD at the start of the read, a
    *left*-side read at the end. Returns ``""`` when the read is too short or
    ``length <= 0``.
    """
    if length <= 0 or len(seq) < length:
        return ""
    return (seq[:length] if side == "right" else seq[-length:]).upper()


def _count_support(ins: Insertion, cluster: _Cluster) -> None:
    """Count bracketing supporting reads (RelocaTE2 ``Supporting_count`` rule)."""
    left = right = 0
    for _name, gstart, gend, strand, _seq in cluster.support:
        if strand == "+" and gend <= ins.start:
            left += 1
        elif strand == "-" and gstart >= ins.end:
            right += 1
    ins.left_support_reads = left
    ins.right_support_reads = right


#: RelocaTE2's -s/--size: sequencing library insert size (relocaTE2.py:198).
DEFAULT_INSERT_SIZE = 500

#: RelocaTE2 extends a one-sided supporting cluster by insert size * (1 + 0.2),
#: its comment reading "insertion size of library * (1 + sd of library)"
#: (relocaTE_insertionFinder.py:446).
_INSERT_SD_FACTOR = 1.2


def call_support_only(
    cluster: _Cluster, insert_size: int = DEFAULT_INSERT_SIZE
) -> Insertion | None:
    """Call an insertion from supporting (mate) reads alone -- RelocaTE2's NONSUP.

    RelocaTE2 writes these to a separate ``all_nonref_supporting`` file rather
    than into ``all_nonref_insert`` (relocaTE_insertionFinder.py:431-459): they
    carry no junction reads at all (T:0 R:0 L:0), so they are much weaker
    evidence and are kept apart from the main call set. RelocaTE3 follows that
    split exactly -- these never enter the insertion tiers.

    Three cases, matching RelocaTE2:

    * **both strands** -- forward mates bracket the site on the left, reverse
      mates on the right; the insertion lies in the gap between the innermost
      reads. Overlapping brackets are ambiguous and rejected (:440).
    * **forward only** -- the site starts at the rightmost forward end and runs
      one library insert onwards (:446).
    * **reverse only** -- the mirror image (:455).

    Unlike :func:`_call_support_only`, no minimum read count is imposed, because
    RelocaTE2 imposes none; the separate output file is what marks these as
    provisional.
    """
    plus_ends = [gend for _n, _s, gend, strand, _q in cluster.support if strand == "+"]
    minus_starts = [
        gstart for _n, gstart, _e, strand, _q in cluster.support if strand == "-"
    ]
    if not plus_ends and not minus_starts:
        return None

    span = int(float(insert_size) * _INSERT_SD_FACTOR)
    if plus_ends and minus_starts:
        ins_start = max(plus_ends)
        ins_end = min(minus_starts)
        if ins_start > ins_end:
            return None  # mates overlap: cannot bracket a gap
    elif plus_ends:
        ins_start = max(plus_ends)
        ins_end = ins_start + span
    else:
        ins_end = min(minus_starts)
        ins_start = max(1, ins_end - span)  # stay on the contig

    return Insertion(
        chrom=cluster.chrom,
        start=ins_start,
        end=ins_end,
        te_name="NA",
        strand="+",
        tsd="supporting_reads",
        left_support_reads=len(plus_ends),
        right_support_reads=len(minus_starts),
        note="Non-reference, supporting reads only",
    )


def _call_support_only(
    cluster: _Cluster, min_support: int = MIN_SUPPORT_ONLY
) -> Insertion | None:
    """Call an insertion from supporting reads alone, when junctions failed to map.

    Mirrors RelocaTE2's both-strand ``supporting_reads`` path: forward-strand mates
    bracket the insertion on the left, reverse-strand mates on the right, and the
    site lies in the gap between the innermost reads (``get_boundary``). Requires
    ``min_support`` reads on each strand. Lower confidence than a junction call —
    short junction flanks often don't map uniquely, but the paired-end mates do.
    """
    plus_ends = [
        gend for _n, _s, gend, strand, _seq in cluster.support if strand == "+"
    ]
    minus_starts = [
        gstart for _n, gstart, _e, strand, _seq in cluster.support if strand == "-"
    ]
    if len(plus_ends) < min_support or len(minus_starts) < min_support:
        return None
    ins_start = max(plus_ends)  # rightmost extent of left-bracketing reads
    ins_end = min(minus_starts)  # leftmost extent of right-bracketing reads
    if ins_start > ins_end:
        return None  # reads overlap: ambiguous, not a clean bracket
    return Insertion(
        chrom=cluster.chrom,
        start=ins_start,
        end=ins_end,
        te_name="NA",
        strand="+",
        tsd="supporting_reads",
        left_support_reads=len(plus_ends),
        right_support_reads=len(minus_starts),
        note="Non-reference, supporting reads only",
    )


def _load_fullread_spans(
    fullreads_bam: str | None,
) -> dict[str, list[tuple[str, int, int]]]:
    """Map read name -> genome spans for full (untrimmed) junction reads."""
    spans: dict[str, list[tuple[str, int, int]]] = {}
    if not fullreads_bam or not Path(fullreads_bam).exists():
        return spans
    with pysam.AlignmentFile(fullreads_bam, "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped:
                continue
            chrom = bam.get_reference_name(rec.reference_id)
            spans.setdefault(rec.query_name, []).append(
                (chrom, rec.reference_start + 1, rec.reference_end)
            )
    return spans


def _strip_junction_tag(name: str) -> str:
    """Remove the :start:5 / :end:3 junction suffix to recover the read name."""
    return _JUNCTION_RE.sub("", name)


def _is_false_junction(
    ins: Insertion, fullread_spans: dict[str, list[tuple[str, int, int]]]
) -> bool:
    """True if the full reads map across the breakpoint, indicating no insertion.

    Mirrors RelocaTE2: if >=30% of the left junction reads AND >=30% of the right
    junction reads have a full read spanning the breakpoint (with margin), the
    site is a reference locus, not an insertion.
    """
    if not fullread_spans:
        return False
    left_total = ins.left_junction_reads
    right_total = ins.right_junction_reads
    if left_total == 0 or right_total == 0:
        return False

    left_full = right_full = 0
    # read_names are ordered left-reads first, then right-reads
    left_names = ins.read_names[:left_total]
    right_names = ins.read_names[left_total:]
    bp_left, bp_right = ins.start, ins.end
    for tagged in left_names:
        if _maps_through(
            fullread_spans.get(_strip_junction_tag(tagged)), ins.chrom, bp_right
        ):
            left_full += 1
    for tagged in right_names:
        if _maps_through(
            fullread_spans.get(_strip_junction_tag(tagged)), ins.chrom, bp_left
        ):
            right_full += 1
    return left_full >= 0.3 * left_total and right_full >= 0.3 * right_total


def _maps_through(
    spans: list[tuple[str, int, int]] | None, chrom: str, breakpoint: int
) -> bool:
    """True if any full-read span covers the breakpoint with margin on both sides."""
    if not spans:
        return False
    for c, s, e in spans:
        if (
            c == chrom
            and s <= breakpoint - FULLREAD_EXTEND
            and e >= breakpoint + FULLREAD_EXTEND
        ):
            return True
    return False


def find_insertions(
    genome_bam: str,
    read_repeat: dict[str, tuple[str, str]],
    genome_fasta: str,
    fullreads_bam: str | None = None,
    required_junction_reads: int = 1,
    include_support_only: bool = True,
) -> list[Insertion]:
    """Call non-reference insertions from the Step-4 genome BAM.

    Keeps clusters with at least ``required_junction_reads`` junction reads on
    either side, dropping false junctions identified via ``fullreads_bam`` (the
    untrimmed junction reads aligned to the genome). When ``include_support_only``
    is set, clusters with no mapped junction reads but two-sided paired-end support
    are also called (lower confidence; recovers sites whose short junction flanks
    failed to map). Returns insertions sorted by chromosome and position.
    """
    fullread_spans = _load_fullread_spans(fullreads_bam)
    insertions: list[Insertion] = []
    n_false = 0
    n_support_only = 0
    with pysam.FastaFile(genome_fasta) as genome:
        for cluster in _stream_clusters(genome_bam, read_repeat):
            calls = []
            for ins in _call_insertions(cluster, genome):
                if (
                    ins.left_junction_reads < required_junction_reads
                    and ins.right_junction_reads < required_junction_reads
                ):
                    continue
                if _is_false_junction(ins, fullread_spans):
                    n_false += 1
                    continue
                calls.append(ins)
            if not calls and include_support_only:
                support_call = _call_support_only(cluster)
                if support_call is not None:
                    calls.append(support_call)
                    n_support_only += 1
            insertions.extend(calls)
    insertions.sort(key=lambda i: (i.chrom, i.start, i.end))
    logger.info(
        "Called %d non-reference insertions (%d junction-based, %d support-only, %d false junctions filtered)",
        len(insertions),
        len(insertions) - n_support_only,
        n_support_only,
        n_false,
    )
    return insertions


def write_insertions_gff(
    insertions: list[Insertion],
    path: str | Path,
    sample: str,
    source: str = "RelocaTE3",
) -> None:
    """Write insertions as GFF3 with the RelocaTE2 attribute set."""
    with open(path, "w") as fh:
        for ins in insertions:
            attrs = (
                f"ID={ins.feature_id};Name={ins.te_name};TSD={ins.tsd};Note={ins.note};"
                f"Right_junction_reads={ins.right_junction_reads};"
                f"Left_junction_reads={ins.left_junction_reads};"
                f"Right_support_reads={ins.right_support_reads};"
                f"Left_support_reads={ins.left_support_reads};"
            )
            fh.write(
                f"{ins.chrom}\t{source}\t{sample}\t{ins.start}\t{ins.end}\t.\t{ins.strand}\t.\t{attrs}\n"
            )


def _gff_attr(attrs: str, key: str, default: str = "") -> str:
    """Extract ``key=value`` from a GFF attribute column."""
    m = re.search(rf"{key}=([^;]*)", attrs)
    return m.group(1) if m else default


def read_insertions_gff(path: str | Path) -> list[Insertion]:
    """Parse a RelocaTE3 insertion GFF back into :class:`Insertion` records."""
    insertions: list[Insertion] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = f[8]
            insertions.append(
                Insertion(
                    chrom=f[0],
                    start=int(f[3]),
                    end=int(f[4]),
                    te_name=_gff_attr(attrs, "Name", "NA"),
                    strand=f[6],
                    tsd=_gff_attr(attrs, "TSD", "UNK"),
                    left_junction_reads=int(
                        _gff_attr(attrs, "Left_junction_reads", "0")
                    ),
                    right_junction_reads=int(
                        _gff_attr(attrs, "Right_junction_reads", "0")
                    ),
                    left_support_reads=int(_gff_attr(attrs, "Left_support_reads", "0")),
                    right_support_reads=int(
                        _gff_attr(attrs, "Right_support_reads", "0")
                    ),
                    note=_gff_attr(attrs, "Note", ""),
                )
            )
    return insertions


def write_insertions_txt(insertions: list[Insertion], path: str | Path) -> None:
    """Write a tab-delimited insertion summary table."""
    header = [
        "chrom",
        "start",
        "end",
        "TE",
        "strand",
        "TSD",
        "right_junction",
        "left_junction",
        "right_support",
        "left_support",
    ]
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for ins in insertions:
            fh.write(
                "\t".join(
                    str(x)
                    for x in (
                        ins.chrom,
                        ins.start,
                        ins.end,
                        ins.te_name,
                        ins.strand,
                        ins.tsd,
                        ins.right_junction_reads,
                        ins.left_junction_reads,
                        ins.right_support_reads,
                        ins.left_support_reads,
                    )
                )
                + "\n"
            )


def write_supporting_reads(
    result_dir: Path | str,
    target: str,
    te_name: str,
    sample: str,
    insertions: list[Insertion],
) -> Path:
    """Write ``<target>.<TE>.all_nonref_supporting.{txt,gff}`` (RelocaTE2 NONSUP).

    Always written, even when empty, so downstream steps never have to test for
    the file's existence -- RelocaTE2 opens it unconditionally too
    (relocaTE_insertionFinder.py:151).
    """
    result_dir = Path(result_dir)
    result_dir.mkdir(parents=True, exist_ok=True)
    txt_path = result_dir / f"{target}.{te_name}.all_nonref_supporting.txt"
    gff_path = result_dir / f"{target}.{te_name}.all_nonref_supporting.gff"
    with open(txt_path, "w") as txt:
        for ins in insertions:
            txt.write(
                f"{ins.te_name}\t{ins.tsd}\t{sample}\t{ins.chrom}\t"
                f"{ins.start}..{ins.end}\t{ins.strand}\tT:0\tR:0\tL:0\t"
                f"ST:{ins.left_support_reads + ins.right_support_reads}\t"
                f"SR:{ins.right_support_reads}\tSL:{ins.left_support_reads}\n"
            )
    write_insertions_gff(insertions, gff_path, sample)
    logger.info(
        "Wrote %d supporting-reads-only insertions to %s",
        len(insertions),
        txt_path,
    )
    return txt_path


# ---------------------------------------------------------------------------
# Tiered output (RelocaTE2 clean_false_positive.py)
# ---------------------------------------------------------------------------

#: Call classes RelocaTE2 strips from its headline output. RelocaTE3 records the
#: class in the TSD column, exactly where RelocaTE2 does.
LOW_CONFIDENCE_CLASSES = frozenset(
    {"singleton", "insufficient_data", "supporting_reads"}
)


def _table_stem(table: Path) -> str:
    """``.../ALL.mping.all_nonref_insert.txt`` -> ``.../ALL.mping.all_nonref_insert``."""
    return str(table)[: -len(table.suffix)] if table.suffix else str(table)


def _row_class(fields: list[str]) -> str:
    """The TSD column, which doubles as the call class for low-confidence calls."""
    return fields[1] if len(fields) > 1 else ""


def _junction_counts(fields: list[str]) -> tuple[int, int] | None:
    """``(right, left)`` junction-read counts, or None if unparseable."""
    if len(fields) < 9:
        return None
    try:
        return (
            int(fields[7].removeprefix("R:")),
            int(fields[8].removeprefix("L:")),
        )
    except ValueError:
        return None


def _has_empty_side(fields: list[str]) -> bool:
    """True when one junction side has no reads at all.

    This is the test RelocaTE2's *boundary* filter uses
    (clean_false_positive.py:82, ``Right == 0 or Left == 0``). It is broader
    than the high_conf rule below -- do not conflate them.
    """
    counts = _junction_counts(fields)
    if counts is None:
        return False
    return (counts[0] == 0) != (counts[1] == 0)


def _is_single_read_one_sided(fields: list[str]) -> bool:
    """True for RelocaTE2's high_conf exclusions: exactly one read against zero.

    clean_false_positive.py:108 greps out the literal
    ``Right_junction_reads=1;Left_junction_reads=0`` and its mirror -- and
    nothing else. A one-sided call backed by several junction reads survives
    into high_conf. Treating every zero-sided call as low confidence is stricter
    than RelocaTE2 and discards real insertions: on the Chr3 2 Mb fixture it
    removes 16 calls RelocaTE2 keeps, costing ~0.05 recall for no precision
    gain.
    """
    counts = _junction_counts(fields)
    if counts is None:
        return False
    return counts in {(1, 0), (0, 1)}


def _load_te_boundaries(reference_ins: Path | str) -> dict[str, set[int]]:
    """``{chrom: {boundary positions}}`` from a RepeatMasker .out or BED."""
    from RelocaTE3.reference_te import ReferenceTEAnnotator

    table = ReferenceTEAnnotator.load_existing_te(reference_ins, "ALL")
    return {
        chrom: set(sides["start"]) | set(sides["end"])
        for chrom, sides in table.items()
    }


def _at_te_boundary(
    fields: list[str], boundaries: dict[str, set[int]], distance: int
) -> bool:
    """True when either endpoint sits within ``distance`` bp of a TE boundary.

    RelocaTE2 compares the call's start and end against the reference copy's
    start and end (clean_false_positive.py:84-91, four ``elif`` arms).
    """
    if len(fields) < 5:
        return False
    positions = boundaries.get(fields[3])
    if not positions:
        return False
    span = fields[4]
    start_s, _, end_s = span.partition("..")
    try:
        ends = (int(start_s), int(end_s or start_s))
    except ValueError:
        return False
    return any(
        abs(end - boundary) <= distance for end in ends for boundary in positions
    )


def _row_to_gff(fields: list[str], sample: str, source: str = "RelocaTE3") -> str | None:
    """Render one table row as GFF3 with RelocaTE2's attribute set."""
    if len(fields) < 12:
        return None
    te_name, tsd, _exper, chrom, span, strand = fields[:6]
    start, _, end = span.partition("..")
    right_j = fields[7].removeprefix("R:")
    left_j = fields[8].removeprefix("L:")
    right_s = fields[10].removeprefix("SR:")
    left_s = fields[11].removeprefix("SL:")
    attrs = (
        f"ID={chrom}.{start}.spanners;Name={te_name};TSD={tsd};"
        f"Note=Non-reference, not found in reference;"
        f"Right_junction_reads={right_j};Left_junction_reads={left_j};"
        f"Right_support_reads={right_s};Left_support_reads={left_s};"
    )
    return f"{chrom}\t{source}\t{sample}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}"


def write_insertion_tiers(
    table: Path | str,
    sample: str,
    reference_ins: Path | str | None = None,
    distance: int = 3,
) -> dict[str, Path]:
    """Derive RelocaTE2's tiered call sets from a written all_nonref_insert table.

    RelocaTE2 publishes several files, not one (``clean_false_positive.py``).
    Users comparing RelocaTE3 against a RelocaTE2 run were comparing RelocaTE3's
    unfiltered calls against RelocaTE2's filtered ones. This writes, alongside
    the table already produced:

    ``<stem>.raw.{txt,gff}``
        Every call RelocaTE3 made.
    ``<stem>.all.{txt,gff}``
        Minus one-sided calls whose breakpoint sits within ``distance`` bp of a
        reference TE boundary (clean_false_positive.py:82-91). Reads from an
        intact reference copy's edge mimic a novel junction; a *two-sided* call
        at a boundary is genuine and is kept. Requires ``reference_ins``;
        without it this tier equals ``.raw``.
    ``<stem>.gff``
        RelocaTE2's headline set: minus ``singleton`` / ``insufficient_data`` /
        ``supporting_reads`` (clean_false_positive.py:107).
    ``<stem>.high_conf.{txt,gff}``
        Additionally minus calls with exactly one junction read on one side and
        none on the other (clean_false_positive.py:108). Note this is narrower
        than "one-sided": a call with several junction reads on a single side
        survives, as it does in RelocaTE2.

    **The plain ``<stem>.txt`` is left exactly as written.** RelocaTE2 runs
    clean_false_positive.py on the GFF only (relocaTE2.py:704) and concatenates
    the table unfiltered (:703), then genotypes that unfiltered table with
    characterizer.pl (:707). Filtering the table here would quietly drop calls
    from genotyping that RelocaTE2 genotypes.

    Args:
        table: the ``all_nonref_insert.txt`` just written.
        sample: sample name, for the GFF's third column.
        reference_ins: RepeatMasker ``.out`` or BED of reference TE copies. The
            boundary filter is skipped when omitted, matching RelocaTE2, which
            bails out when its overlap file comes back empty
            (clean_false_positive.py:64).
        distance: RelocaTE2's ``-d/--distance``, default 3.

    Returns:
        Mapping of tier name (``raw``/``all``/``final``/``high_conf``) to its
        GFF path.
    """
    table = Path(table)
    stem = _table_stem(table)
    rows = [
        line.split("\t")
        for line in table.read_text().splitlines()
        if line.strip()
    ]

    raw = rows
    if reference_ins:
        boundaries = _load_te_boundaries(reference_ins)
        every = [
            r
            for r in raw
            if not (_has_empty_side(r) and _at_te_boundary(r, boundaries, distance))
        ]
    else:
        every = list(raw)
    headline = [r for r in every if _row_class(r) not in LOW_CONFIDENCE_CLASSES]
    high_conf = [r for r in headline if not _is_single_read_one_sided(r)]

    written: dict[str, Path] = {}
    for name, suffix, subset in (
        ("raw", ".raw", raw),
        ("all", ".all", every),
        ("final", "", headline),
        ("high_conf", ".high_conf", high_conf),
    ):
        gff_path = Path(f"{stem}{suffix}.gff")
        with open(gff_path, "w") as fh:
            for fields in subset:
                line = _row_to_gff(fields, sample)
                if line:
                    fh.write(line + "\n")
        written[name] = gff_path
        # Companion tables for the suffixed tiers only. The plain <stem>.txt is
        # the caller's own output and stays unfiltered (see the docstring).
        if suffix:
            with open(f"{stem}{suffix}.txt", "w") as fh:
                for fields in subset:
                    fh.write("\t".join(fields) + "\n")

    logger.info(
        "Insertion tiers: %d raw, %d all, %d final, %d high_conf",
        len(raw),
        len(every),
        len(headline),
        len(high_conf),
    )
    return written
