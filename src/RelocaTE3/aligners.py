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
import shutil
import subprocess
import tempfile
from abc import ABC, abstractmethod
from contextlib import contextmanager
from pathlib import Path

import pysam



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


class AlignerBackend(ABC):
    """One aligner binary, exposing the two RelocaTE3 alignment stages."""

    name: str = ""

    def __init__(self, threads: int = 1):
        """Store the default thread count for this backend."""
        self.threads = threads

    def _threads(self, threads):
        return threads if threads and threads > 0 else self.threads

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

    def __init__(self, threads: int = 1):
        """Wrap a :class:`RelocaTE3.align.Aligner` configured for minimap2."""
        super().__init__(threads)
        from RelocaTE3.align import Aligner

        self._aln = Aligner(threads)

    def index(self, reference, *, force: bool = False) -> None:
        """Create the minimap2 ``.mmi`` index."""
        self._aln.index_minimap(str(reference), force=force)

    def map_te_library(
        self, reads, te_library, outdir, *, threads=None, tmpdir=None
    ) -> list[Path]:
        """Delegate to ``Aligner.map_minimap_library`` (per-side BAMs)."""
        return self._aln.map_minimap_library(
            reads,
            str(outdir),
            str(te_library),
            tmpdir=str(tmpdir) if tmpdir else "",
            cpu_threads=self._threads(threads),
        )

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
        """Delegate to ``Aligner.map_reads_to_genome`` (SE, concatenated)."""
        return self._aln.map_reads_to_genome(
            str(genome),
            [str(f) for f in fastqs],
            str(out_bam),
            tmpdir=str(tmpdir) if tmpdir else "",
            cpu_threads=self._threads(threads),
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

    def _mem(self, reference, read_files, out_bam, threads, tmpdir, extra=()):
        sam = os.path.join(tmpdir, "bwa.sam")
        with open(sam, "w") as fh:
            subprocess.run(
                [
                    self.binary,
                    "mem",
                    "-t",
                    str(threads),
                    *extra,
                    str(reference),
                    *[str(f) for f in read_files],
                ],
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
        return bams

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
                return self._mem(genome, fastqs[:2], out_bam, threads, td)
            combined = _concat(fastqs, os.path.join(td, "reads.fq"))
            return self._mem(genome, [combined], out_bam, threads, td)


class BwaMem2Backend(BwaBackend):
    """bwa-mem2 backend -- same algorithm/output as bwa-mem, different binary."""

    name = "bwamem2"
    binary = "bwa-mem2"
    _index_sentinel = ".bwt.2bit.64"


class Bowtie2Backend(AlignerBackend):
    """bowtie2 backend. TE stage uses -k to keep multi-mappers for multi-copy TEs."""

    name = "bowtie2"

    def index(self, reference, *, force: bool = False) -> None:
        """Build the bowtie2 index unless it already exists."""
        if force or not Path(f"{reference}.1.bt2").exists():
            subprocess.run(
                ["bowtie2-build", "--quiet", str(reference), str(reference)], check=True
            )

    def _run(self, reference, read_args, out_bam, threads, tmpdir):
        sam = os.path.join(tmpdir, "bt2.sam")
        subprocess.run(
            [
                "bowtie2",
                "-p",
                str(threads),
                "-x",
                str(reference),
                "-S",
                sam,
                *read_args,
            ],
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
        return bams

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
            return self._run(genome, args, out_bam, threads, td)


def psl_to_sam(psl_lines):
    """Convert headerless BLAT PSL records to SAM alignment lines.

    Handles the common short-read case (single or multi-block, ungapped or with
    small gaps). CIGAR is soft-clip + M/I/D from block coordinates; ``NM`` is the
    PSL mismatch + inserted-base counts. ``-`` strand sets flag 0x10 (the read is
    reported reverse-complemented by BLAT, so the sequence is emitted as ``*``).
    """
    out = []
    for line in psl_lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 21 or not f[0].isdigit():
            continue  # skip PSL header / blank lines
        mismatches = int(f[1])
        q_base_ins = int(f[5])
        t_base_ins = int(f[7])
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
        out.append(
            f"{qname}\t{flag}\t{tname}\t{pos}\t255\t{cigar_str}\t*\t0\t0\t*\t*\tNM:i:{nm}"
        )
    return out


class BlatBackend(AlignerBackend):
    """BLAT backend -- TE-library search only (PSL output, single-end).

    BLAT does not emit SAM or support paired-end mapping, so ``map_genome``
    raises. Provide ``blat`` on PATH (it is intentionally not pinned in
    ``pixi.toml``; see the note there).
    """

    name = "blat"

    def index(self, reference, *, force: bool = False) -> None:
        """BLAT indexes at run time; nothing to pre-build."""

    def _blat_side(self, te_library, read_file, out_bam, threads, tmpdir):
        # BLAT wants FASTA/FASTQ query; write a header for SAM assembly.
        psl = os.path.join(tmpdir, "aln.psl")
        subprocess.run(
            ["blat", str(te_library), str(read_file), "-noHead", "-out=psl", psl],
            check=True,
        )
        ref_lengths = {
            r: fa.get_reference_length(r)
            for fa in [pysam.FastaFile(str(te_library))]
            for r in fa.references
        }
        sam = os.path.join(tmpdir, "aln.sam")
        with open(sam, "w") as fh:
            fh.write("@HD\tVN:1.6\tSO:unsorted\n")
            for r, ln in ref_lengths.items():
                fh.write(f"@SQ\tSN:{r}\tLN:{ln}\n")
            with open(psl) as ph:
                for rec in psl_to_sam(ph):
                    fh.write(rec + "\n")
        return _sam_to_sorted_mapped_bam(sam, out_bam, threads)

    def map_te_library(
        self, reads, te_library, outdir, *, threads=None, tmpdir=None
    ) -> list[Path]:
        """Map each read side to the TE library via BLAT (PSL -> SAM)."""
        threads = self._threads(threads)
        read_set = {"left": reads.left()}
        if reads.is_paired:
            read_set["right"] = reads.right()
        bams = []
        with _tmp(tmpdir) as td:
            for side, read_file in read_set.items():
                out_bam = Path(outdir) / f"{reads.name}.{side}.bam"
                self._blat_side(te_library, read_file, out_bam, threads, td)
                bams.append(out_bam)
        return bams

    def map_genome(self, *args, **kwargs) -> Path:
        """Not supported -- BLAT is TE-search only."""
        raise NotImplementedError(
            "blat genome alignment is not supported; use "
            "--genome-aligner {minimap2,bwa,bwamem2,bowtie2}"
        )


BACKENDS: dict[str, type[AlignerBackend]] = {
    "minimap2": MinimapBackend,
    "bwa": BwaBackend,
    "bwamem2": BwaMem2Backend,
    "bowtie2": Bowtie2Backend,
    "blat": BlatBackend,
}

TE_ALIGNERS = tuple(BACKENDS)  # all backends can do TE-library search
GENOME_ALIGNERS = tuple(n for n in BACKENDS if n != "blat")  # blat has no genome stage


def get_aligner(name: str, threads: int = 1) -> AlignerBackend:
    """Instantiate the aligner backend registered under ``name``."""
    try:
        backend = BACKENDS[name]
    except KeyError:
        raise ValueError(
            f"unknown aligner {name!r}; choices: {sorted(BACKENDS)}"
        ) from None
    return backend(threads)
