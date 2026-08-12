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
    """Convert headerless BLAT PSL records to SAM alignment lines.

    Handles the common short-read case (single or multi-block, ungapped or with
    small gaps). CIGAR is soft-clip + M/I/D from block coordinates; ``NM`` is the
    PSL mismatch + inserted-base counts. ``-`` strand sets flag 0x10.

    ``sequences`` (a ``{qname: forward_seq}`` map) fills the SAM SEQ field
    (reverse-complemented for ``-`` strand, per SAM convention) -- required so the
    trim step can emit flanking reads, since BLAT's PSL carries no sequence. When
    omitted, SEQ is ``*``.
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
        seq_field = "*"
        if sequences is not None:
            s = sequences.get(qname)
            if s:
                seq_field = _revcomp(s) if strand == "-" else s
        out.append(
            f"{qname}\t{flag}\t{tname}\t{pos}\t255\t{cigar_str}\t*\t0\t0\t{seq_field}\t*\tNM:i:{nm}"
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

    def _blat_cmd(self, te_library, query_fa, psl):
        # Runs BLAT at its defaults (minScore=30, tileSize=11). We briefly
        # hardcoded RelocaTE2's -minScore=10 -tileSize=7 here (relocaTE2.py:545)
        # on the theory that BLAT sensitivity was the recall gap. Benchmarking on
        # riceTElib (9 samples, 3 coverages) disproved it: those params raise
        # recall by ~0.05 but collapse precision as coverage rises -- 0.635 /
        # 0.406 / 0.231 at 5x / 15x / 30x, versus RelocaTE2's steady ~0.83 --
        # because false positives grow superlinearly with depth (286 -> 1222 ->
        # 3320) while RelocaTE2's stay flat. The extra calls are short-TSD
        # (median 3 bp) LINE/SINE junctions that RelocaTE2 discards downstream
        # and we do not, so the permissive settings are only safe once that
        # filter exists. At BLAT's defaults RelocaTE3 already beat RelocaTE2 on
        # F1 at every coverage (0.444 / 0.645 / 0.729 vs 0.438 / 0.620 / 0.696).
        # See relocate-benchmark docs/2026-08-10-ricetelib-parity-results.md.
        #
        # The sensitized values remain available per-run via --te-opts, which is
        # appended below (before the output path, so it can override anything
        # defaulted here). The TE library is the database and the reads are the
        # query, so db precedes query.
        return [
            "blat",
            str(te_library),
            str(query_fa),
            *self._stage_opts("te"),
            "-noHead",
            "-out=psl",
            str(psl),
        ]

    def _blat_side(self, te_library, read_file, out_bam, threads, tmpdir):
        # BLAT reads FASTA, not FASTQ (it treats a FASTQ's "@name" lines as a list
        # of filenames). Read libraries are FASTQ, so convert the query to FASTA
        # first; FastxFile handles FASTQ/FASTA (and gzip) and preserves the /1,/2
        # mate suffix in the name.
        query_fa = os.path.join(tmpdir, f"query.{Path(out_bam).stem}.fa")
        seqs: dict[str, str] = {}
        with pysam.FastxFile(str(read_file)) as fx, open(query_fa, "w") as out:
            for rec in fx:
                out.write(f">{rec.name}\n{rec.sequence}\n")
                seqs[rec.name] = rec.sequence
        psl = os.path.join(tmpdir, "aln.psl")
        subprocess.run(self._blat_cmd(te_library, query_fa, psl), check=True)
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
                for rec in psl_to_sam(ph, sequences=seqs):
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
        return canonicalize_te_bams(bams, reads)

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
