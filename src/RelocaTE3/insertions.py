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

        # R2 reports the read-captured TSD whenever one was inferred, including
        # for single-sided junctions; the "supporting_junction" / "singleton"
        # sentinels only kick in when no usable TSD string was captured.
        if top_tsd and top_tsd not in {"UNK", "UKN"}:
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
# max separation (bp) between a left and right breakpoint to call a shared TSD
TSD_WINDOW = 100
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
    if read_name in read_repeat:
        return read_repeat[read_name][0]
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


def _stream_clusters(bam_path: str, read_repeat: dict[str, tuple[str, str]]):
    """Yield :class:`_Cluster` objects by streaming a coordinate-sorted BAM."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        current: _Cluster | None = None
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped:
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
    left_pos: list[int], right_pos: list[int]
) -> list[tuple[int | None, int | None]]:
    """Pair left/right breakpoints into sub-insertions (RelocaTE2 pairing).

    Each returned tuple is (left_position, right_position); either may be None
    for a one-sided junction. Left and right are paired when within ``TSD_WINDOW``.
    """
    nL, nR = len(left_pos), len(right_pos)
    if nL and nR:
        if nL == 1 and nR == 1:
            if abs(left_pos[0] - right_pos[0]) > TSD_WINDOW:
                return [(left_pos[0], None), (None, right_pos[0])]
            return [(left_pos[0], right_pos[0])]
        # greedy nearest-neighbour pairing within the TSD window
        pairs: list[tuple[int | None, int | None]] = []
        remaining_right = sorted(right_pos)
        for lp in sorted(left_pos):
            best = None
            for rp in remaining_right:
                if best is None or abs(rp - lp) < abs(best - lp):
                    best = rp
            if best is not None and abs(best - lp) <= TSD_WINDOW:
                pairs.append((lp, best))
                remaining_right.remove(best)
            else:
                pairs.append((lp, None))
        for rp in remaining_right:
            pairs.append((None, rp))
        return pairs
    if nL:
        return [(lp, None) for lp in sorted(left_pos)]
    return [(None, rp) for rp in sorted(right_pos)]


def _resolve_tsd(
    left_reads: list[JunctionObservation],
    right_reads: list[JunctionObservation],
    chrom: str,
    i_start: int,
    i_end: int,
    tsd_len: int,
    genome: pysam.FastaFile,
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
    fetched = _fetch_tsd(genome, chrom, i_start, i_end)
    return fetched or "UNK"


def _make_insertion(
    chrom: str,
    left_reads: list[JunctionObservation],
    right_reads: list[JunctionObservation],
    genome: pysam.FastaFile,
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
    te_names = [j.te_name for j in junctions if j.te_name != "NA"]
    te_name = max(set(te_names), key=te_names.count) if te_names else "NA"

    orients = [j.te_orientation for j in junctions]
    strand = "+" if orients.count("+") >= orients.count("-") else "-"

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


def _call_insertions(cluster: _Cluster, genome: pysam.FastaFile) -> list[Insertion]:
    """Split a cluster into one or more insertions by pairing breakpoints."""
    left = _group_by_position([j for j in cluster.junctions if j.side == "left"])
    right = _group_by_position([j for j in cluster.junctions if j.side == "right"])
    if not left and not right:
        return []

    insertions: list[Insertion] = []
    for lp, rp in _pair_breakpoints(list(left), list(right)):
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
    tuples), then for each fractional threshold (in order) counts contiguous
    positions whose depth >= ``threshold * len(spans)``. Returns the first
    non-zero length, or 0 if none qualify.

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
