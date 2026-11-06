"""Pluggable aligner backends for RelocaTE3.

RelocaTE3 has two alignment stages: TE-library search (``map_te_library``) and
genome re-alignment (``map_genome``). Each backend wraps one aligner binary and
returns the contract every downstream step assumes: a coordinate-sorted, indexed
BAM containing mapped reads only, with read names preserved and (where the
aligner provides it) an ``NM`` tag.

Backends are stateless subprocess wrappers with explicit file I/O so they map
cleanly onto Nextflow processes and keep the aligner choice out of the Python
parsing that is slated for a Rust rewrite. ``minimap2`` (the default) delegates
to the original, acceptance-tested :class:`RelocaTE3.align.Aligner` so its
behavior is preserved exactly.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import uuid
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass
from pathlib import Path

import pysam


@dataclass(frozen=True)
class TEChunk:
    """TE-library hits for one contiguous block of query reads.

    ``bam`` is a coordinate-sorted, indexed BAM holding the hits of exactly the
    ``n_reads`` consecutive records of the side's read file that this chunk
    covers (records ``[sum of earlier n_reads, + n_reads)``). Because a read's
    hits never leave its chunk, per-read best-hit selection over one chunk is
    the same as over the whole side.
    """

    bam: Path
    n_reads: int


@dataclass(frozen=True)
class ChunkedTEAlignment:
    """One read side's TE-library alignment, kept as its per-chunk parts.

    ``bam`` is the merged side BAM (the usual ``map_te_library`` contract);
    ``chunks`` are its parts in read-file order; ``chunk_dir`` holds them and
    can be removed once the trim step has consumed them.
    """

    bam: Path
    chunks: tuple[TEChunk, ...]
    chunk_dir: Path


@dataclass(frozen=True)
class TESide:
    """One read side's TE-search chunks, ready for :func:`~RelocaTE3.librelocate.merge_trim`.

    ``direction`` is ``left`` or ``right`` (the side BAM is
    ``<name>.<direction>.bam``); ``read_file`` is the FASTQ the chunks were cut
    from, in order; ``chunks`` are the chunk BAMs in that order (possibly none,
    for an empty read file). ``chunk_dir``, when set, is a directory holding
    only this side's chunk BAMs, which the caller may remove after the merge.
    """

    direction: str
    read_file: Path
    chunks: tuple[TEChunk, ...]
    chunk_dir: Path | None = None


#: Query reads per BLAT chunk, matching RelocaTE2's ``fastq_split -s``
#: (relocaTE2.py:91). BLAT is single-threaded, so chunking is the only way to
#: use more than one core (or, with the scatter subcommands, more than one job).
BLAT_CHUNK_SEQS = 200_000


@contextmanager
def _tmp(tmpdir):
    """Yield a scratch dir, creating a self-cleaning temp dir when none given."""
    if tmpdir:
        Path(tmpdir).mkdir(parents=True, exist_ok=True)
        yield str(tmpdir)
    else:
        handle = tempfile.TemporaryDirectory()
        try:
            yield handle.name
        finally:
            handle.cleanup()


def _sam_to_sorted_mapped_bam(sam: str, out_bam, threads: int) -> Path:
    """Shared contract tail: drop unmapped reads, coordinate-sort, index.

    Secondary/supplementary alignments (multi-mappers, needed for multi-copy TE
    families) are retained -- only the unmapped flag (0x4) is filtered.
    """
    out_bam = Path(out_bam)
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    mapped = str(out_bam) + ".mapped.tmp.bam"
    subprocess.run(
        ["samtools", "view", "-b", "-F", "0x4", "-o", mapped, sam], check=True
    )
    pysam.sort("-@", str(threads), "-o", str(out_bam), mapped)
    pysam.index(str(out_bam))
    os.unlink(mapped)
    return out_bam


def _concat(fastqs, dest: str) -> str:
    """Concatenate FASTQ files into ``dest`` (mirrors Aligner.map_reads_to_genome)."""
    with open(dest, "wb") as out:
        for f in fastqs:
            with open(f, "rb") as src:
                shutil.copyfileobj(src, out)
    return dest


def _count_fasta_records(path) -> int:
    """Count ``>`` header lines in a FASTA file without loading it."""
    count = 0
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                count += 1
    return count


def _merge_sorted_bams(bams, out_bam, sam_header: str, threads: int, tmpdir) -> Path:
    """``samtools merge`` coordinate-sorted BAMs into ``out_bam`` and index it.

    Streams the inputs, so memory does not grow with the number of records.
    With no inputs, writes a header-only BAM so the output contract holds.
    """
    out_bam = Path(out_bam)
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    if bams:
        listing = Path(tmpdir) / f"{out_bam.name}.merge.list"
        listing.write_text("".join(f"{b}\n" for b in bams))
        subprocess.run(
            [
                "samtools", "merge", "-f", "--no-PG", "-@", str(max(0, int(threads) - 1)),
                "-b", str(listing), "-o", str(out_bam),
            ],
            check=True,
        )
        listing.unlink()
    else:
        header_sam = Path(tmpdir) / f"{out_bam.name}.header.sam"
        header_sam.write_text(sam_header.replace("SO:unsorted", "SO:coordinate"))
        subprocess.run(
            ["samtools", "view", "--no-PG", "-b", "-o", str(out_bam), str(header_sam)],
            check=True,
        )
        header_sam.unlink()
    subprocess.run(["samtools", "index", str(out_bam)], check=True)
    return out_bam


def _fai_is_complete(fasta: Path, fai: Path) -> bool:
    """Whether ``fai`` is a whole ``samtools faidx`` index of ``fasta``.

    Every line must have the five integer-valued columns, and the last record
    must end where the FASTA does (only whitespace after it): an index cut
    short at a line boundary lacks the last records, so ``>`` follows.
    """
    try:
        text = fai.read_text()
    except OSError:
        return False
    if not text.endswith("\n"):
        return False
    try:
        for line in text.splitlines():
            name, length, offset, line_bases, line_width = line.split("\t")
            length, offset = int(length), int(offset)
            line_bases, line_width = int(line_bases), int(line_width)
    except ValueError:
        return False
    full_lines, rest = divmod(length, line_bases) if line_bases else (0, 0)
    end = offset + full_lines * line_width + rest
    with open(fasta, "rb") as fh:
        fh.seek(end)
        while block := fh.read(1 << 16):
            if block.strip():
                return False
    return fasta.stat().st_size >= end


def ensure_fasta_index(fasta) -> Path:
    """Make sure ``<fasta>.fai`` exists and is complete; return its path.

    htslib (``pysam.FastaFile``) writes a missing ``.fai`` in place, so a job
    opening the library while another is still writing it reads a truncated
    index -- with several ``te-search-chunk`` jobs on one shared TE library
    that gives chunk BAMs a short ``@SQ`` list. Here the index is built under
    a temporary name in the same directory and renamed into place, so readers
    see either no index or a whole one; concurrent builders each rename an
    identical file. A complete existing index is left alone; an incomplete
    one is rebuilt.
    """
    fasta = Path(fasta)
    fai = Path(f"{fasta}.fai")
    if _fai_is_complete(fasta, fai):
        return fai
    # a name no other builder uses; samtools creates it (so the umask applies,
    # as for an index htslib writes itself)
    tmp = str(fai.with_name(f".{fai.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp"))
    try:
        pysam.faidx(str(fasta), "--fai-idx", tmp)
        os.replace(tmp, fai)
    finally:
        Path(tmp).unlink(missing_ok=True)
    return fai


def te_sam_header(te_library) -> str:
    """SAM header (``@HD`` + one ``@SQ`` per TE) for BLAT hits on ``te_library``.

    It depends only on the library, so every chunk search (in any process)
    and the merge of an empty side compute the same header independently.
    The library's ``.fai`` is created atomically (:func:`ensure_fasta_index`)
    so concurrent calls never read a partial index.
    """
    ensure_fasta_index(te_library)
    with pysam.FastaFile(str(te_library)) as fa:
        return "@HD\tVN:1.6\tSO:unsorted\n" + "".join(
            f"@SQ\tSN:{r}\tLN:{fa.get_reference_length(r)}\n" for r in fa.references
        )


def sq_lines(sam_header: str) -> list[str]:
    """The ``@SQ`` lines of a SAM header text, in order."""
    return [line for line in sam_header.splitlines() if line.startswith("@SQ")]


def check_chunk_headers(bams, sam_header: str) -> None:
    """Refuse chunk BAMs whose ``@SQ`` list is not ``sam_header``'s.

    ``samtools merge`` silently unions differing dictionaries, so a chunk
    searched against another (or a half-indexed) TE library would otherwise
    slip into the side BAM.

    Raises:
        ValueError: naming the first chunk BAM whose ``@SQ`` lines differ.
    """
    want = sq_lines(sam_header)
    for bam in bams:
        with pysam.AlignmentFile(str(bam), "rb") as fh:
            got = sq_lines(str(fh.header))
        if got != want:
            raise ValueError(
                f"chunk BAM {bam} has {len(got)} @SQ line(s) that do not match the "
                f"TE library's {len(want)}; it was searched against another TE "
                "library (or a partly written .fai of this one). Re-run "
                "te-search-chunk for it with the same --te-library"
            )


_CHUNK_INDEX_DIGITS = 5  # split -a 5: chunk.<stem>.00000.fa ... 99999


def split_te_query(
    read_file, outdir, stem: str, chunk_reads: int = BLAT_CHUNK_SEQS
) -> list[tuple[Path, int]]:
    """Convert one read side to FASTA and cut it into BLAT query chunks.

    Step 1 of the BLAT TE search (then :func:`search_te_chunk` per chunk, then
    :func:`RelocaTE3.librelocate.merge_trim`). Chunks are written to ``outdir``
    as ``chunk.<stem>.<NNNNN>.fa`` (zero-padded index, file order), each
    holding the next ``chunk_reads`` records of ``read_file``.

    BLAT reads FASTA, not FASTQ (it treats a FASTQ's ``@name`` lines as a list
    of filenames), so the query has to be converted. RelocaTE2 uses ``seqtk``
    for this (relocaTE2.py step 2) and so do we: doing it in Python dominates
    the TE-search stage on a full-size library -- roughly 4M reads per 35 min,
    hours for one 30x side before BLAT even starts. ``seqtk seq -A -l 0`` emits
    exactly two lines per record and preserves the ``/1``, ``/2`` mate suffix,
    so plain ``split -l`` cuts on record boundaries. Falls back to a pure-Python
    conversion (same chunk files) when seqtk is absent.

    Args:
        read_file: FASTQ or FASTA input, optionally gzipped.
        outdir: directory for the chunk files (created if needed).
        stem: name fragment making chunk names unique, e.g. ``<sample>.left``.
        chunk_reads: records per chunk (default :data:`BLAT_CHUNK_SEQS`). The
            outputs of the whole TE search do not depend on it.

    Returns:
        ``(chunk FASTA path, number of reads in it)`` per chunk, in file order;
        empty for an empty read file.

    Raises:
        FileExistsError: if ``outdir`` already holds chunks for ``stem`` (they
            would be mixed with the new ones).
    """
    chunk_reads = int(chunk_reads)
    if chunk_reads < 1:
        raise ValueError(f"chunk_reads must be at least 1, got {chunk_reads}")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = f"chunk.{stem}."
    pattern = re.compile(re.escape(prefix) + r"\d{%d}\.fa$" % _CHUNK_INDEX_DIGITS)

    def listing() -> list[Path]:
        return sorted(p for p in outdir.iterdir() if pattern.match(p.name))

    stale = listing()
    if stale:
        raise FileExistsError(
            f"{outdir} already holds {len(stale)} query chunk(s) for {stem} "
            f"({stale[0].name}, ...); remove them or use an empty directory"
        )

    if shutil.which("seqtk") is not None:
        query_fa = outdir / f".query.{stem}.fa"
        try:
            with open(query_fa, "w") as out:
                subprocess.run(
                    ["seqtk", "seq", "-A", "-l", "0", str(read_file)],
                    stdout=out,
                    check=True,
                )
            subprocess.run(
                [
                    "split",
                    "-l",
                    str(2 * chunk_reads),
                    "-d",
                    "-a",
                    str(_CHUNK_INDEX_DIGITS),
                    "--additional-suffix=.fa",
                    str(query_fa),
                    str(outdir / prefix),
                ],
                check=True,
            )
        finally:
            query_fa.unlink(missing_ok=True)
        return [(path, _count_fasta_records(path)) for path in listing()]

    chunks: list[tuple[Path, int]] = []
    out = None
    try:
        with pysam.FastxFile(str(read_file)) as fx:
            for i, rec in enumerate(fx):
                if i % chunk_reads == 0:
                    if out is not None:
                        out.close()
                    path = outdir / f"{prefix}{len(chunks):0{_CHUNK_INDEX_DIGITS}d}.fa"
                    chunks.append((path, 0))
                    out = open(path, "w")
                out.write(f">{rec.name}\n{rec.sequence}\n")
                chunks[-1] = (chunks[-1][0], chunks[-1][1] + 1)
    finally:
        if out is not None:
            out.close()
    return chunks


_MATE_SUFFIX_RE = re.compile(r"^(.*)/([12])$")


def canonical_name(raw_name: str, mate: str) -> str:
    """Return ``<base>/<mate>``, whatever convention the FASTQ used.

    RelocaTE3 identifies mates by a trailing ``/1``/``/2``. Plenty of real data
    does not use that: reads off a modern Illumina instrument carry the mate in
    a separate field (``@A00519:... 1:N:0:ATCACG``), so the *name* is identical
    in R1 and R2; SRA/ENA dumps often carry no mate marker at all. Matching such
    names against ``<base>/2`` silently found nothing, so no genomic mates were
    recovered and no junction flank was mate-anchored -- with no error, just
    quietly worse calls.

    Which file a read came from settles the question, so ``mate`` ('1' or '2')
    is authoritative: any comment field and any existing suffix are stripped
    before it is applied. An empty ``mate`` just normalises the base name.
    """
    base = raw_name.split(None, 1)[0] if raw_name else raw_name
    m = _MATE_SUFFIX_RE.match(base)
    if m:
        base = m.group(1)
    return f"{base}/{mate}" if mate else base


def _restore_mate_suffix(bam, mate: str) -> None:
    """Rewrite read names in ``bam`` to the canonical ``<base>/{mate}`` form.

    The TE stage maps each side (R1/R2) to its own BAM, so the side fixes the
    mate. Some aligners strip a trailing mate suffix (bwa mem) and some inputs
    never had one (Illumina, SRA), so the name alone cannot be trusted --
    without this, downstream flank-pairing finds no genomic mate for any read.
    Rewrites in place (coordinate order unchanged) and re-indexes.
    """
    bam = Path(bam)
    tmp = bam.with_suffix(".matefix.bam")
    with pysam.AlignmentFile(str(bam), "rb") as inp:
        with pysam.AlignmentFile(str(tmp), "wb", template=inp) as out:
            for rec in inp.fetch(until_eof=True):
                rec.query_name = canonical_name(rec.query_name, mate)
                out.write(rec)
    tmp.replace(bam)
    pysam.index(str(bam))


def canonicalize_te_bams(bams: list, reads) -> list:
    """Stamp the mate onto every TE-library BAM name, per side.

    ``bams`` are the per-side BAMs in ``[left, right]`` order as returned by the
    backends' ``map_te_library``. Single-end libraries are left alone: there is
    no mate to recover, and adding ``/1`` would only churn names.
    """
    if not getattr(reads, "is_paired", False):
        return bams
    for index, bam in enumerate(bams[:2]):
        _restore_mate_suffix(bam, "1" if index == 0 else "2")
    return bams


class AlignerBackend(ABC):
    """One aligner binary, exposing the two RelocaTE3 alignment stages."""

    name: str = ""

    def __init__(self, threads: int = 1, te_opts=(), genome_opts=()):
        """Store the thread count and per-stage extra aligner options.

        ``te_opts``/``genome_opts`` are appended verbatim to the TE-search and
        genome-placement command lines. They exist so sensitivity knobs can be
        swept from configuration (blat ``-minIdentity``, minimap2 ``-B``,
        bowtie2 ``--very-sensitive-local``, bwa tuning) without code changes.
        The two stages are kept separate because they want opposite settings:
        TE search benefits from permissive, multi-mapping alignment while
        genome placement needs near-unique hits to preserve precision.
        """
        self.threads = threads
        self.te_opts = tuple(te_opts or ())
        self.genome_opts = tuple(genome_opts or ())

    def _threads(self, threads):
        return threads if threads and threads > 0 else self.threads

    def _stage_opts(self, stage: str) -> list[str]:
        """Extra options for ``stage`` ('te' or 'genome')."""
        return list(self.te_opts if stage == "te" else self.genome_opts)

    @abstractmethod
    def index(self, reference, *, force: bool = False) -> None:
        """Build whatever index this aligner needs for ``reference``."""

    @abstractmethod
    def map_te_library(
        self, reads, te_library, outdir, *, threads=None, tmpdir=None
    ) -> list[Path]:
        """Map reads to the TE library; return per-side ``{name}.{side}.bam`` list."""

    @abstractmethod
    def map_genome(
        self,
        genome,
        fastqs,
        out_bam,
        *,
        paired: bool = False,
        threads=None,
        tmpdir=None,
    ) -> Path:
        """Map (flanking/support) reads to the genome; return the sorted BAM."""


class MinimapBackend(AlignerBackend):
    """minimap2 backend -- delegates to the original Aligner (behavior-preserving)."""

    name = "minimap2"

    def __init__(self, threads: int = 1, te_opts=(), genome_opts=()):
        """Wrap a :class:`RelocaTE3.align.Aligner` configured for minimap2."""
        super().__init__(threads, te_opts=te_opts, genome_opts=genome_opts)
        from RelocaTE3.align import Aligner

        self._aln = Aligner(threads)

    def index(self, reference, *, force: bool = False) -> None:
        """Create the minimap2 ``.mmi`` index."""
        self._aln.index_minimap(str(reference), force=force)

    def map_te_library(
        self, reads, te_library, outdir, *, threads=None, tmpdir=None
    ) -> list[Path]:
        """Delegate to ``Aligner.map_minimap_library`` (per-side BAMs)."""
        bams = self._aln.map_minimap_library(
            reads,
            str(outdir),
            str(te_library),
            tmpdir=str(tmpdir) if tmpdir else "",
            cpu_threads=self._threads(threads),
            extra_opts=self._stage_opts("te"),
        )
        return canonicalize_te_bams(bams, reads)

    def map_genome(
        self,
        genome,
        fastqs,
        out_bam,
        *,
        paired: bool = False,
        threads=None,
        tmpdir=None,
    ) -> Path:
        """Map to the genome. Paired runs align R1/R2 with ``minimap2 -ax sr`` so a
        uniquely-mapping mate anchors an ambiguous junction flank; otherwise every
        FASTQ is aligned single-end (concatenated, behavior-preserving)."""
        if paired and len(fastqs) >= 2:
            t = self._threads(threads)
            with _tmp(tmpdir) as td:
                sam = os.path.join(td, "genome.pe.sam")
                # -k 11 -w 5 matches the single-end junction-flank params (short
                # 10-15 bp flanks need tight seeds); -a for SAM output.
                cmd = [
                    "minimap2", "-t", str(t), "-a", "-x", "sr",
                    "-k", "11", "-w", "5",
                    *self._stage_opts("genome"),
                    str(genome), str(fastqs[0]), str(fastqs[1]),
                ]
                with open(sam, "w") as fh:
                    subprocess.run(cmd, stdout=fh, check=True)
                return _sam_to_sorted_mapped_bam(sam, out_bam, t)
        return self._aln.map_reads_to_genome(
            str(genome),
            [str(f) for f in fastqs],
            str(out_bam),
            tmpdir=str(tmpdir) if tmpdir else "",
            cpu_threads=self._threads(threads),
            extra_opts=self._stage_opts("genome"),
        )


class BwaBackend(AlignerBackend):
    """bwa-mem backend. TE stage keeps all alignments (-a) for multi-copy TEs."""

    name = "bwa"
    binary = "bwa"
    _index_sentinel = ".bwt"

    def index(self, reference, *, force: bool = False) -> None:
        """Build the bwa index unless it already exists."""
        if force or not Path(f"{reference}{self._index_sentinel}").exists():
            subprocess.run([self.binary, "index", str(reference)], check=True)

    def _mem_cmd(self, reference, read_files, *, threads=1, extra=(), stage="te"):
        """Build the ``bwa mem`` command for ``stage`` ('te' or 'genome')."""
        return [
            self.binary,
            "mem",
            "-t",
            str(threads),
            *extra,
            *self._stage_opts(stage),
            str(reference),
            *[str(f) for f in read_files],
        ]

    def _mem(
        self, reference, read_files, out_bam, threads, tmpdir, extra=(), stage="te"
    ):
        sam = os.path.join(tmpdir, "bwa.sam")
        with open(sam, "w") as fh:
            subprocess.run(
                self._mem_cmd(
                    reference, read_files, threads=threads, extra=extra, stage=stage
                ),
                stdout=fh,
                check=True,
            )
        return _sam_to_sorted_mapped_bam(sam, out_bam, threads)

    def map_te_library(
        self, reads, te_library, outdir, *, threads=None, tmpdir=None
    ) -> list[Path]:
        """Map each read side to the TE library (SE, -a to keep multi-mappers)."""
        threads = self._threads(threads)
        self.index(te_library)
        read_set = {"left": reads.left()}
        if reads.is_paired:
            read_set["right"] = reads.right()
        bams = []
        with _tmp(tmpdir) as td:
            for side, read_file in read_set.items():
                out_bam = Path(outdir) / f"{reads.name}.{side}.bam"
                self._mem(te_library, [read_file], out_bam, threads, td, extra=["-a"])
                bams.append(out_bam)
        return canonicalize_te_bams(bams, reads)

    def map_genome(
        self,
        genome,
        fastqs,
        out_bam,
        *,
        paired: bool = False,
        threads=None,
        tmpdir=None,
    ) -> Path:
        """Map flanking/support reads to the genome (SE concat, or R1/R2 if paired)."""
        threads = self._threads(threads)
        self.index(genome)
        with _tmp(tmpdir) as td:
            if paired and len(fastqs) >= 2:
                return self._mem(
                    genome, fastqs[:2], out_bam, threads, td, stage="genome"
                )
            combined = _concat(fastqs, os.path.join(td, "reads.fq"))
            return self._mem(
                genome, [combined], out_bam, threads, td, stage="genome"
            )


class BwaMem2Backend(BwaBackend):
    """bwa-mem2 backend -- same algorithm/output as bwa-mem, different binary."""

    name = "bwamem2"
    binary = "bwa-mem2"
    _index_sentinel = ".bwt.2bit.64"


class BwaAlnBackend(BwaBackend):
    """bwa ``aln``/``samse``/``sampe`` backend for the genome stage (RelocaTE2
    parity).

    ``bwa mem`` (and minimap2 ``-ax sr``) impose a seed-length floor (mem's
    ``-k`` defaults to 19), so they silently drop the short (<~19 bp) junction
    flanks that decide low-coverage / low-VAF insertions. ``bwa aln`` is the
    backtracking short-read algorithm with no such floor and is exactly what
    RelocaTE2 uses for flank->genome placement; using it here makes the genome
    step a true RelocaTE2 parity. It is a subcommand of the same modern bwa
    binary -- not an older version -- and its SAM is piped through the shared
    sorted/mapped/indexed-BAM helper (no RelocaTE2-style intermediate shell
    scripts). ``aln`` takes no tuning flags here, matching RelocaTE2's defaults.

    TE-library search (``map_te_library``) is inherited from :class:`BwaBackend`
    (bwa mem); the parity variant uses BLAT for that stage, so this backend is
    used only for ``map_genome``.
    """

    name = "bwaaln"

    def _aln_cmd(self, reference, fastq, threads):
        """Build the ``bwa aln`` command (genome stage).

        Extra genome options (e.g. ``-n 0.10`` to widen the edit-distance budget
        for divergent flanks, or ``-l``/``-k`` seed tuning) are appended here.
        """
        return [
            self.binary,
            "aln",
            "-t",
            str(threads),
            *self._stage_opts("genome"),
            str(reference),
            str(fastq),
        ]

    def _aln_sai(self, reference, fastq, sai, threads):
        with open(sai, "wb") as fh:
            subprocess.run(
                self._aln_cmd(reference, fastq, threads),
                stdout=fh,
                check=True,
            )

    def _aln_to_sam(self, reference, fastqs, sam, threads, tmpdir, paired):
        if paired and len(fastqs) >= 2:
            sai1 = os.path.join(tmpdir, "r1.sai")
            sai2 = os.path.join(tmpdir, "r2.sai")
            self._aln_sai(reference, fastqs[0], sai1, threads)
            self._aln_sai(reference, fastqs[1], sai2, threads)
            cmd = [
                self.binary, "sampe", str(reference),
                sai1, sai2, str(fastqs[0]), str(fastqs[1]),
            ]
        else:
            combined = _concat(fastqs, os.path.join(tmpdir, "reads.fq"))
            sai = os.path.join(tmpdir, "se.sai")
            self._aln_sai(reference, combined, sai, threads)
            cmd = [self.binary, "samse", str(reference), sai, combined]
        with open(sam, "w") as fh:
            subprocess.run(cmd, stdout=fh, check=True)

    def map_genome(
        self,
        genome,
        fastqs,
        out_bam,
        *,
        paired: bool = False,
        threads=None,
        tmpdir=None,
    ) -> Path:
        """Map flanking/support reads to the genome with bwa ``aln`` -- paired
        (``sampe``, so a flank is anchored by its mate) or single-end (``samse``,
        concatenated). Short flanks that ``bwa mem`` drops are placed here."""
        threads = self._threads(threads)
        self.index(genome)
        with _tmp(tmpdir) as td:
            sam = os.path.join(td, "bwaaln.sam")
            self._aln_to_sam(genome, fastqs, sam, threads, td, paired)
            return _sam_to_sorted_mapped_bam(sam, out_bam, threads)


class Bowtie2Backend(AlignerBackend):
    """bowtie2 backend. TE stage uses -k to keep multi-mappers for multi-copy TEs."""

    name = "bowtie2"

    def index(self, reference, *, force: bool = False) -> None:
        """Build the bowtie2 index unless it already exists."""
        if force or not Path(f"{reference}.1.bt2").exists():
            subprocess.run(
                ["bowtie2-build", "--quiet", str(reference), str(reference)], check=True
            )

    def _bt2_cmd(self, reference, read_args, sam, threads, stage="te"):
        """Build the bowtie2 command for ``stage`` ('te' or 'genome').

        Extra options (e.g. ``--very-sensitive-local``/``-N 1`` for divergent TE
        copies) are appended after ``read_args`` so they override the presets.
        """
        return [
            "bowtie2",
            "-p",
            str(threads),
            "-x",
            str(reference),
            "-S",
            str(sam),
            *read_args,
            *self._stage_opts(stage),
        ]

    def _run(self, reference, read_args, out_bam, threads, tmpdir, stage="te"):
        sam = os.path.join(tmpdir, "bt2.sam")
        subprocess.run(
            self._bt2_cmd(reference, read_args, sam, threads, stage=stage),
            check=True,
        )
        return _sam_to_sorted_mapped_bam(sam, out_bam, threads)

    def map_te_library(
        self, reads, te_library, outdir, *, threads=None, tmpdir=None
    ) -> list[Path]:
        """Map each read side to the TE library (SE, -k 20 to keep multi-mappers).

        Uses ``--local`` so bowtie2 soft-clips TE-junction reads (part TE, part
        genomic flank) instead of dropping them: in the default end-to-end mode a
        read must align over its full length, so junction reads -- the ones that
        carry the non-reference insertion signal -- align 0 times and no flank is
        recovered. ``--local`` matches the soft-clipping behavior of the bwa-mem
        and minimap2 backends.
        """
        threads = self._threads(threads)
        self.index(te_library)
        read_set = {"left": reads.left()}
        if reads.is_paired:
            read_set["right"] = reads.right()
        bams = []
        with _tmp(tmpdir) as td:
            for side, read_file in read_set.items():
                out_bam = Path(outdir) / f"{reads.name}.{side}.bam"
                self._run(
                    te_library,
                    ["-k", "20", "--local", "-U", str(read_file)],
                    out_bam,
                    threads,
                    td,
                )
                bams.append(out_bam)
        return canonicalize_te_bams(bams, reads)

    def map_genome(
        self,
        genome,
        fastqs,
        out_bam,
        *,
        paired: bool = False,
        threads=None,
        tmpdir=None,
    ) -> Path:
        """Map flanking/support reads to the genome (SE concat, or -1/-2 if paired)."""
        threads = self._threads(threads)
        self.index(genome)
        with _tmp(tmpdir) as td:
            if paired and len(fastqs) >= 2:
                args = ["-1", str(fastqs[0]), "-2", str(fastqs[1])]
            else:
                combined = _concat(fastqs, os.path.join(td, "reads.fq"))
                args = ["-U", combined]
            return self._run(genome, args, out_bam, threads, td, stage="genome")


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def psl_to_sam(psl_lines, sequences=None):
    """Return SAM lines as a list, preserving the existing converter interface.

    Production file writing uses ``_iter_psl_to_sam`` directly to avoid
    retaining every converted alignment. Both paths share conversion logic.
    """
    return list(_iter_psl_to_sam(psl_lines, sequences=sequences))


def _iter_psl_to_sam(psl_lines, sequences=None):
    """Yield SAM lines in input order, retaining only the current alignment.

    Handles the common short-read case (single or multi-block, ungapped or with
    small gaps). CIGAR is soft-clip + M/I/D from block coordinates; ``NM`` is the
    PSL mismatch + inserted-base counts. ``-`` strand sets flag 0x10.

    Before conversion, reject gap-complex alignments using RelocaTE2's exact
    PSL admission thresholds: at most one query and target gap opening, at most
    three inserted query and target bases, and fewer than three alignment
    blocks. Rejected alignments therefore cannot participate in downstream
    per-read best-hit selection.

    ``sequences`` (a ``{qname: forward_seq}`` map) fills the SAM SEQ field
    (reverse-complemented for ``-`` strand, per SAM convention) -- required so the
    trim step can emit flanking reads, since BLAT's PSL carries no sequence. When
    omitted, SEQ is ``*``.
    """
    for line in psl_lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 21 or not f[0].isdigit():
            continue  # skip PSL header / blank lines
        mismatches = int(f[1])
        q_num_ins = int(f[4])
        q_base_ins = int(f[5])
        t_num_ins = int(f[6])
        t_base_ins = int(f[7])
        block_count = int(f[17])
        if (
            q_num_ins > 1
            or q_base_ins > 3
            or t_num_ins > 1
            or t_base_ins > 3
            or block_count >= 3
        ):
            continue
        strand = f[8]
        qname, qsize = f[9], int(f[10])
        qstart, qend = int(f[11]), int(f[12])
        tname = f[13]
        block_sizes = [int(x) for x in f[18].rstrip(",").split(",") if x]
        q_starts = [int(x) for x in f[19].rstrip(",").split(",") if x]
        t_starts = [int(x) for x in f[20].rstrip(",").split(",") if x]

        cigar = []
        lead_clip = qstart if strand == "+" else qsize - qend
        if lead_clip:
            cigar.append((lead_clip, "S"))
        for i, bsize in enumerate(block_sizes):
            if i > 0:
                q_gap = q_starts[i] - (q_starts[i - 1] + block_sizes[i - 1])
                t_gap = t_starts[i] - (t_starts[i - 1] + block_sizes[i - 1])
                if q_gap:
                    cigar.append((q_gap, "I"))
                if t_gap:
                    cigar.append((t_gap, "D"))
            cigar.append((bsize, "M"))
        tail_clip = qsize - qend if strand == "+" else qstart
        if tail_clip:
            cigar.append((tail_clip, "S"))

        cigar_str = "".join(f"{n}{op}" for n, op in cigar if n > 0)
        flag = 16 if strand == "-" else 0
        nm = mismatches + q_base_ins + t_base_ins
        pos = t_starts[0] + 1  # SAM is 1-based
        seq_field = "*"
        if sequences is not None:
            s = sequences.get(qname)
            if s:
                seq_field = _revcomp(s) if strand == "-" else s
        yield (
            f"{qname}\t{flag}\t{tname}\t{pos}\t255\t{cigar_str}\t*\t0\t0\t{seq_field}\t*\tNM:i:{nm}"
        )


class BlatBackend(AlignerBackend):
    """BLAT backend -- TE-library search only (PSL output, single-end).

    BLAT does not emit SAM or support paired-end mapping, so ``map_genome``
    raises. Provide ``blat`` on PATH (it is intentionally not pinned in
    ``pixi.toml``; see the note there).
    """

    name = "blat"

    def index(self, reference, *, force: bool = False) -> None:
        """BLAT indexes at run time; nothing to pre-build."""

    @staticmethod
    def _require_blat() -> None:
        """Fail early, and usefully, when ``blat`` is not installed.

        blat is the default TE aligner (matching RelocaTE2) but cannot be pinned
        in the main pixi environment -- the bioconda build hard-pins
        zlib 1.2.11 and ships linux-64 only. So the common first-run failure is
        a missing binary, and without this it surfaced as a bare
        FileNotFoundError from subprocess deep inside the trim step, which reads
        like a RelocaTE3 bug rather than a missing dependency.
        """
        if shutil.which("blat") is None:
            raise RuntimeError(
                "blat is not on PATH. It is the default TE aligner (matching "
                "RelocaTE2) but is not part of the default pixi environment "
                "because its bioconda build conflicts with the plotting stack "
                "and is linux-64 only.\n"
                "  Install it:  pixi run -e blat true   (then use `pixi run -e "
                "blat relocaTE3 ...`)\n"
                "  Or pick another aligner, e.g. --te-aligner bowtie2; on the "
                "rice benchmark every alternative finds far fewer insertions "
                "than blat (about half; see README, 'Choosing an aligner')"
            )

    #: RelocaTE2's BLAT sensitivity settings (relocaTE2.py:545). BLAT's own
    #: defaults are minScore=30, tileSize=11.
    RT2_BLAT_OPTS = ("-minScore=10", "-tileSize=7")

    def _blat_cmd(self, te_library, query_fa, psl):
        # Runs BLAT with RelocaTE2's parameters, which is the whole point of
        # defaulting to BLAT at all: at BLAT's own defaults the TE search simply
        # does not find a large share of the junction reads RelocaTE2 finds.
        # Measured on mPing cov30x_rep1, over the 19 junction reads RelocaTE2
        # used at sites RelocaTE3 missed entirely: BLAT defaults find 0 of 19,
        # -minScore=10 -tileSize=7 finds 19 of 19.
        #
        # History, because this was reverted once (02e68eb) and the reasoning
        # was sound at the time: benchmarking on riceTElib showed these params
        # raising recall but collapsing precision with depth (0.635 / 0.406 /
        # 0.231 at 5x / 15x / 30x). The extra calls were overwhelmingly
        # one-sided junctions -- and RelocaTE3 was then admitting every
        # one-sided call at characterize, because characterize.py ported
        # RelocaTE2's gate (characterizer.pl:91) with `or` where RelocaTE2 has
        # `and`. With that gate corrected, one-sided calls are dropped the way
        # RelocaTE2 drops them, which is precisely the "filter" the old comment
        # said these settings needed before they were safe. Re-measure both
        # panels if you touch either of these two things again; they interact.
        #
        # --te-opts is appended after, so a run can still override any of this.
        # The TE library is the database and the reads are the query, so db
        # precedes query.
        return [
            "blat",
            str(te_library),
            str(query_fa),
            *self.RT2_BLAT_OPTS,
            *self._stage_opts("te"),
            "-noHead",
            "-out=psl",
            str(psl),
        ]

    #: Query sequences per BLAT chunk (:data:`BLAT_CHUNK_SEQS`); an instance
    #: or class override changes the chunking of the in-process search.
    BLAT_CHUNK_SEQS = BLAT_CHUNK_SEQS

    @staticmethod
    def _query_sequences(query_fa, psl, tmpdir) -> dict[str, str]:
        """Load sequences for just the reads BLAT matched.

        BLAT's PSL carries no sequence, so the PSL->SAM conversion needs the
        query sequences to fill SEQ. Retain only matched reads, which can still
        number in the millions for a complex TE library. Stream the selected
        FASTA without retaining an additional whole-output string/line list.
        The seqtk path spools to a temporary file so subprocess failure is
        checked before parsing, without buffering its stdout in Python memory.
        """
        names: set[str] = set()
        with open(psl) as ph:
            for line in ph:
                fields = line.split("\t")
                if len(fields) > 9:
                    names.add(fields[9])
        if not names:
            return {}

        seqs: dict[str, str] = {}
        with ExitStack() as stack:
            if shutil.which("seqtk") is not None:
                names_file = os.path.join(tmpdir, "matched.names")
                with open(names_file, "w") as fh:
                    fh.writelines(f"{name}\n" for name in names)
                stream = stack.enter_context(
                    tempfile.TemporaryFile(mode="w+t", dir=tmpdir)
                )
                subprocess.run(
                    ["seqtk", "subseq", str(query_fa), names_file],
                    stdout=stream,
                    check=True,
                )
                stream.seek(0)
            else:
                stream = stack.enter_context(open(query_fa))
            name = None
            for line in stream:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                elif name is not None and name in names:
                    seqs[name] = seqs.get(name, "") + line.strip()
        return seqs

    #: ``samtools sort -m`` for one chunk. Chunks are small (200k reads), and
    #: up to ``threads`` sorts run at once, so keep each one modest.
    CHUNK_SORT_MEM = "256M"

    def _blat_chunk(self, te_library, chunk_fa, chunk_bam, sam_header, mate, workdir):
        """BLAT one query chunk and turn its hits into a sorted, indexed BAM.

        This is RelocaTE2's per-chunk unit (relocaTE2.py ``--split``): only
        this chunk's matched sequences are loaded to fill SAM SEQ. Read names
        get their canonical ``/1``/``/2`` here (``mate``), exactly what
        :func:`canonicalize_te_bams` would do to the merged BAM, so the chunk
        BAMs are usable by the trim step directly.
        """
        work = Path(workdir)
        work.mkdir(parents=True, exist_ok=True)
        psl = work / "aln.psl"
        subprocess.run(self._blat_cmd(te_library, chunk_fa, psl), check=True)
        seqs = self._query_sequences(chunk_fa, psl, str(work))
        sam = work / "aln.sam"
        with open(sam, "w") as fh:
            fh.write(sam_header)
            with open(psl) as ph:
                for rec in _iter_psl_to_sam(ph, sequences=seqs):
                    if mate:
                        qname, rest = rec.split("\t", 1)
                        rec = f"{canonical_name(qname, mate)}\t{rest}"
                    fh.write(rec + "\n")
        del seqs
        mapped = work / "mapped.bam"
        subprocess.run(
            ["samtools", "view", "--no-PG", "-b", "-F", "0x4", "-o", str(mapped), str(sam)],
            check=True,
        )
        subprocess.run(
            [
                "samtools", "sort", "--no-PG", "-m", self.CHUNK_SORT_MEM,
                "-T", str(work / "sort"), "-o", str(chunk_bam), str(mapped),
            ],
            check=True,
        )
        subprocess.run(["samtools", "index", str(chunk_bam)], check=True)
        n_reads = _count_fasta_records(chunk_fa)
        shutil.rmtree(work)
        return TEChunk(bam=Path(chunk_bam), n_reads=n_reads)

    def _search_side(
        self, te_library, read_file, chunk_dir, stem, threads, tmpdir, mate=None
    ) -> list[TEChunk]:
        """Split one read side and search its chunks concurrently (steps 1-2).

        :func:`split_te_query` cuts the side into ``BLAT_CHUNK_SEQS``-read
        chunk FASTAs under ``tmpdir``; each is searched on its own
        (:meth:`_blat_chunk`, the body of :func:`search_te_chunk`) into
        ``chunk_dir/chunk.<NNNNN>.bam``, ``threads`` at a time. BLAT is
        single-threaded and RT2_BLAT_OPTS make it far slower per base, so one
        process over a whole library is not viable: a single mPing 5x side ran
        >13 min unchunked while RelocaTE2 does the entire sample in about that.
        RelocaTE2 splits into 200k-read chunks and runs a process pool
        (relocaTE2.py:91,110); this is the same strategy, carried through
        PSL->BAM conversion so no step holds more than one chunk's hits or
        sequences. The chunk FASTAs are removed once searched.

        Returns the chunks in read-file order.
        """
        chunk_dir = Path(chunk_dir)
        chunk_dir.mkdir(parents=True, exist_ok=True)
        queries = split_te_query(
            read_file, Path(tmpdir) / f"query.{stem}", stem, self.BLAT_CHUNK_SEQS
        )
        sam_header = te_sam_header(te_library)
        parts: list[TEChunk] = []
        if queries:
            workers = max(1, min(int(threads or 1), len(queries)))
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = [
                    pool.submit(
                        self._blat_chunk,
                        te_library,
                        chunk,
                        chunk_dir / f"chunk.{i:05d}.bam",
                        sam_header,
                        mate,
                        os.path.join(tmpdir, f"work.{stem}.{i:05d}"),
                    )
                    for i, (chunk, _n) in enumerate(queries)
                ]
                parts = [fut.result() for fut in futures]  # re-raise the first failure
        shutil.rmtree(Path(tmpdir) / f"query.{stem}", ignore_errors=True)
        return parts

    def _blat_side(
        self, te_library, read_file, out_bam, threads, tmpdir, mate=None, chunk_dir=None
    ) -> ChunkedTEAlignment:
        """Search one read side in chunks; keep the chunk BAMs and merge them.

        The side BAM is ``samtools merge`` of the chunk BAMs. Both sort and merge
        are stable (ties keep input order, then file order), so it has the same
        records in the same order as sorting the concatenated chunks at once.
        """
        stem = Path(out_bam).stem
        chunk_dir = Path(chunk_dir) if chunk_dir else Path(out_bam).with_suffix(".te_chunks")
        if chunk_dir.exists():
            shutil.rmtree(chunk_dir)  # stale parts from an interrupted run
        parts = self._search_side(
            te_library, read_file, chunk_dir, stem, threads, tmpdir, mate=mate
        )
        _merge_sorted_bams(
            [p.bam for p in parts], out_bam, te_sam_header(te_library), threads, tmpdir
        )
        return ChunkedTEAlignment(bam=Path(out_bam), chunks=tuple(parts), chunk_dir=chunk_dir)

    @staticmethod
    def _sides(reads):
        """``(direction, read file, mate)`` per side, in ``[left, right]`` order.

        Paired sides get canonical ``/1``, ``/2`` names (as
        :func:`canonicalize_te_bams`); single-end names are left alone.
        """
        sides = [("left", reads.left())]
        if reads.is_paired:
            sides.append(("right", reads.right()))
        return [
            (side, read_file, str(index + 1) if reads.is_paired else None)
            for index, (side, read_file) in enumerate(sides)
        ]

    def search_te_chunks(
        self, reads, te_library, outdir, *, threads=None, tmpdir=None
    ) -> list[TESide]:
        """Steps 1-2 of the BLAT TE search for every read side, in-process.

        Returns one :class:`TESide` per side, ``[left, right]``, whose chunk
        BAMs are in ``<outdir>/<name>.<side>.te_chunks/`` (``chunk_dir``). Hand
        them to :func:`RelocaTE3.librelocate.merge_trim` (step 3), then remove
        ``chunk_dir``. This is what ``relocaTE3 run`` does; the scatter
        subcommands run the same steps in separate processes.
        """
        self._require_blat()
        threads = self._threads(threads)
        found = []
        with _tmp(tmpdir) as td:
            for side, read_file, mate in self._sides(reads):
                chunk_dir = Path(outdir) / f"{reads.name}.{side}.te_chunks"
                if chunk_dir.exists():
                    shutil.rmtree(chunk_dir)  # stale parts from an interrupted run
                parts = self._search_side(
                    te_library, read_file, chunk_dir, f"{reads.name}.{side}",
                    threads, td, mate=mate,
                )
                found.append(
                    TESide(
                        direction=side,
                        read_file=Path(read_file),
                        chunks=tuple(parts),
                        chunk_dir=chunk_dir,
                    )
                )
        return found

    def map_te_library_chunked(
        self, reads, te_library, outdir, *, threads=None, tmpdir=None
    ) -> list[ChunkedTEAlignment]:
        """Map each read side via BLAT and keep its per-chunk BAMs.

        Returns one :class:`ChunkedTEAlignment` per side in ``[left, right]``
        order. The merged side BAMs are the same files :meth:`map_te_library`
        returns; the chunk BAMs (``<outdir>/<name>.<side>.te_chunks/``) let the
        trim step work one chunk at a time. The caller removes ``chunk_dir``
        when done with it.
        """
        self._require_blat()
        threads = self._threads(threads)
        alignments = []
        with _tmp(tmpdir) as td:
            for side, read_file, mate in self._sides(reads):
                out_bam = Path(outdir) / f"{reads.name}.{side}.bam"
                alignments.append(
                    self._blat_side(
                        te_library,
                        read_file,
                        out_bam,
                        threads,
                        td,
                        mate=mate,
                        chunk_dir=Path(outdir) / f"{reads.name}.{side}.te_chunks",
                    )
                )
        return alignments

    def map_te_library(
        self, reads, te_library, outdir, *, threads=None, tmpdir=None
    ) -> list[Path]:
        """Map each read side to the TE library via BLAT (PSL -> SAM)."""
        alignments = self.map_te_library_chunked(
            reads, te_library, outdir, threads=threads, tmpdir=tmpdir
        )
        for alignment in alignments:
            shutil.rmtree(alignment.chunk_dir, ignore_errors=True)
        return [alignment.bam for alignment in alignments]

    def map_genome(self, *args, **kwargs) -> Path:
        """Not supported -- BLAT is TE-search only."""
        raise NotImplementedError(
            "blat genome alignment is not supported; use "
            "--genome-aligner {minimap2,bwa,bwamem2,bowtie2}"
        )


def search_te_chunk(
    te_library, chunk_fa, chunk_bam, mate=None, te_opts=(), *, workdir=None, sam_header=None
) -> TEChunk:
    """BLAT one query chunk against the TE library into a sorted, indexed BAM.

    Step 2 of the BLAT TE search (after :func:`split_te_query`, before
    :func:`RelocaTE3.librelocate.merge_trim`); exactly what the in-process
    search runs per chunk (:meth:`BlatBackend._blat_chunk`): RelocaTE2's BLAT
    options (:attr:`BlatBackend.RT2_BLAT_OPTS`) followed by ``te_opts``, hits
    converted to SAM with this chunk's sequences, read names canonicalised to
    ``<base>/<mate>`` when ``mate`` is ``'1'``/``'2'`` (paired sides; ``None``
    for single-end), unmapped dropped, ``samtools sort``, index.

    Args:
        te_library: TE library FASTA (the BLAT database).
        chunk_fa: one chunk FASTA from :func:`split_te_query`.
        chunk_bam: output BAM path (``.bai`` is written next to it).
        mate: ``'1'``, ``'2'`` or ``None``.
        te_opts: extra BLAT options (``--te-opts``), appended after RT2's.
        workdir: scratch directory, removed afterwards (default
            ``<chunk_bam>.work``).
        sam_header: precomputed :func:`te_sam_header` of ``te_library``.

    Returns:
        :class:`TEChunk` with the BAM and the chunk's read count.
    """
    backend = BlatBackend(te_opts=te_opts)
    backend._require_blat()
    chunk_bam = Path(chunk_bam)
    chunk_bam.parent.mkdir(parents=True, exist_ok=True)
    return backend._blat_chunk(
        te_library,
        chunk_fa,
        chunk_bam,
        sam_header if sam_header is not None else te_sam_header(te_library),
        mate,
        workdir if workdir is not None else f"{chunk_bam}.work",
    )


BACKENDS: dict[str, type[AlignerBackend]] = {
    "minimap2": MinimapBackend,
    "bwa": BwaBackend,
    "bwamem2": BwaMem2Backend,
    "bwaaln": BwaAlnBackend,
    "bowtie2": Bowtie2Backend,
    "blat": BlatBackend,
}

TE_ALIGNERS = tuple(BACKENDS)  # all backends can do TE-library search
GENOME_ALIGNERS = tuple(n for n in BACKENDS if n != "blat")  # blat has no genome stage


def get_aligner(
    name: str, threads: int = 1, te_opts=(), genome_opts=()
) -> AlignerBackend:
    """Instantiate the aligner backend registered under ``name``.

    ``te_opts``/``genome_opts`` are extra command-line options appended to the
    TE-search and genome-placement invocations respectively (see
    :class:`AlignerBackend`).
    """
    try:
        backend = BACKENDS[name]
    except KeyError:
        raise ValueError(
            f"unknown aligner {name!r}; choices: {sorted(BACKENDS)}"
        ) from None
    return backend(threads, te_opts=te_opts, genome_opts=genome_opts)
