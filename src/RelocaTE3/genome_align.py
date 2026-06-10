"""Step 4: re-align trimmed flanking reads (and their mates) to the genome.

This replaces RelocaTE2's ``relocaTE_align.py`` + ``clean_pairs_memory.py``. The
trimmed flanking reads (precise junction breakpoints) and the genomic *mates* of
TE-containing reads (paired-end support) are mapped to the reference genome and
merged into a single coordinate-sorted BAM consumed by the insertion finder.

Reads are mapped single-end: the insertion finder distinguishes junction reads
from supporting reads by the read-name tag (``:start:5`` etc.) and genomic
proximity, not by BAM proper-pair flags.
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path

import pysam

from RelocaTE3 import logger
from RelocaTE3.align import Aligner
from RelocaTE3.ReadLibrary import ReadLibrary

# matches the junction-tag suffixes appended during trimming
_TAG_RE = re.compile(r":(?:start|end):[53]$|:middle$")
# matches only the junction (5'/3') tags, not :middle
_JUNCTION_TAG_RE = re.compile(r":(?:start|end):[53]$")
# splits a trailing /1 or /2 mate designator
_MATE_RE = re.compile(r"^(.*)/([12])$")


def strip_tag(name: str) -> str:
    """Remove a RelocaTE3 junction tag (``:start:5`` / ``:middle`` ...) from a name."""
    return _TAG_RE.sub("", name)


def split_mate(name: str) -> tuple[str, str]:
    """Split a read name into (base, mate) where mate is '1', '2', or ''."""
    m = _MATE_RE.match(name)
    if m:
        return m.group(1), m.group(2)
    return name, ""


def read_read_repeat(path: Path) -> dict[str, tuple[str, str]]:
    """Load a ``read_repeat_name.txt`` table written by the trim step."""
    read_repeat: dict[str, tuple[str, str]] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                read_repeat[parts[0]] = (parts[1], parts[2])
    return read_repeat


def _fetch_reads_by_name(
    reads: ReadLibrary, needed: dict[str, set[str]], out_fastq: Path
) -> int:
    """Stream the original FASTQ(s) and write reads whose names are in ``needed``.

    ``needed`` maps mate-end ('1'/'2') to the set of read names (including the
    ``/1`` or ``/2`` suffix) to pull from that file. Returns the count written.
    """
    written = 0
    file_for_end = {"1": reads.left(), "2": reads.right()}
    with open(out_fastq, "w") as out:
        for end, names in needed.items():
            if not names or not file_for_end.get(end):
                continue
            with pysam.FastxFile(file_for_end[end]) as fx:
                for rec in fx:
                    if rec.name in names:
                        qual = (
                            rec.quality
                            if rec.quality is not None
                            else "I" * len(rec.sequence)
                        )
                        out.write(f"@{rec.name}\n{rec.sequence}\n+\n{qual}\n")
                        written += 1
    return written


def recover_support_mates(
    read_repeat: dict[str, tuple[str, str]],
    reads: ReadLibrary,
    out_fastq: Path,
) -> int:
    """Write the genomic mates of TE-containing reads to ``out_fastq``.

    A supporting read is the mate of a TE-containing read when that mate did not
    itself match the TE (i.e. it lands in unique genome sequence and brackets the
    insertion). Returns the number of supporting reads written.
    """
    if not reads.is_paired:
        out_fastq.write_text("")
        return 0

    # (base, mate-end) pairs that matched a TE
    te_ends: set[tuple[str, str]] = set()
    for tagged in read_repeat:
        base, mate = split_mate(strip_tag(tagged))
        te_ends.add((base, mate))

    # mate names we need to pull, grouped by which original file holds them
    needed: dict[str, set[str]] = {"1": set(), "2": set()}
    for base, mate in te_ends:
        if mate not in ("1", "2"):
            continue
        other = "2" if mate == "1" else "1"
        if (base, other) not in te_ends:
            needed[other].add(f"{base}/{other}")
    return _fetch_reads_by_name(reads, needed, out_fastq)


def collect_junction_fullreads(
    read_repeat: dict[str, tuple[str, str]],
    reads: ReadLibrary,
    out_fastq: Path,
) -> int:
    """Write the full (untrimmed) sequences of junction reads to ``out_fastq``.

    These are re-aligned to the genome so the insertion finder can drop false
    junctions whose full read maps cleanly across the breakpoint. Returns the
    count written.
    """
    needed: dict[str, set[str]] = {"1": set(), "2": set()}
    for tagged in read_repeat:
        if not _JUNCTION_TAG_RE.search(tagged):
            continue  # only 5'/3' junction reads, not :middle
        full_name = strip_tag(tagged)  # e.g. read_500_470/1
        _base, mate = split_mate(full_name)
        if mate in ("1", "2"):
            needed[mate].add(full_name)
    return _fetch_reads_by_name(reads, needed, out_fastq)


def align_to_genome(
    reads: ReadLibrary,
    genome: str,
    outdir: str | Path,
    aligner: Aligner | None = None,
    threads: int = 1,
) -> tuple[Path, Path | None]:
    """Map trimmed flanking reads + supporting mates to ``genome``.

    Expects the trim step to have populated ``<outdir>/flanking`` and
    ``<outdir>/te_containing/<sample>.read_repeat_name.txt``. Also aligns the full
    (untrimmed) junction reads to support false-junction filtering. Returns
    ``(genome_bam, fullreads_bam)``; the fullreads BAM is None if there are no
    junction reads.
    """
    outdir = Path(outdir)
    sample = reads.name
    aligner = aligner or Aligner(threads)

    flanking_dir = outdir / "flanking"
    flanking_files = sorted(
        str(p) for p in flanking_dir.glob(f"{sample}.*.flankingReads.fq")
    )
    if not flanking_files:
        raise FileNotFoundError(
            f"No flanking reads found in {flanking_dir}; run the trim step first."
        )

    genome_dir = outdir / "genome_aln"
    genome_dir.mkdir(parents=True, exist_ok=True)

    rr_path = outdir / "te_containing" / f"{sample}.read_repeat_name.txt"
    read_repeat = read_read_repeat(rr_path) if rr_path.exists() else {}

    fastq_inputs = list(flanking_files)
    fullreads_bam: Path | None = None
    with tempfile.TemporaryDirectory() as tmp:
        support_fq = Path(tmp) / f"{sample}.support.fq"
        n_support = recover_support_mates(read_repeat, reads, support_fq)
        if n_support > 0:
            fastq_inputs.append(str(support_fq))
        logger.info("%s: %d supporting mate reads recovered", sample, n_support)

        outbam = genome_dir / f"{sample}.genome.bam"
        aligner.map_reads_to_genome(
            genome, fastq_inputs, str(outbam), tmpdir=tmp, cpu_threads=threads
        )

        # full (untrimmed) junction reads for false-junction filtering
        full_fq = Path(tmp) / f"{sample}.fullreads.fq"
        n_full = collect_junction_fullreads(read_repeat, reads, full_fq)
        if n_full > 0:
            fullreads_bam = genome_dir / f"{sample}.fullreads.genome.bam"
            aligner.map_reads_to_genome(
                genome,
                [str(full_fq)],
                str(fullreads_bam),
                tmpdir=tmp,
                cpu_threads=threads,
            )
        logger.info(
            "%s: %d full junction reads aligned for false-junction filtering",
            sample,
            n_full,
        )

    return outbam, fullreads_bam
