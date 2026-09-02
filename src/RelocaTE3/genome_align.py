"""Step 4: re-align trimmed flanking reads (and their mates) to the genome.

This replaces RelocaTE2's ``relocaTE_align.py`` + ``clean_pairs_memory.py``. The
trimmed flanking reads (precise junction breakpoints) and the genomic *mates* of
TE-containing reads (paired-end support) are mapped to the reference genome and
merged into a single coordinate-sorted BAM consumed by the insertion finder.

RelocaTE2's paired-read state machine is preserved: some junctions are mapped
with a genomic mate, two junction mates are mapped together, and the remaining
junction/support evidence is mapped single-end.
"""

from __future__ import annotations

import re
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pysam

from RelocaTE3 import logger
from RelocaTE3.aligners import canonical_name, get_aligner
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


def read_te_hit_names(path: Path) -> set[str]:
    """Load names of all reads with a selected, aligner-admitted TE hit."""
    with open(path) as fh:
        return {strip_tag(line.strip()) for line in fh if line.strip()}


def _fetch_reads_by_name(
    reads: ReadLibrary, needed: dict[str, set[str]], out_fastq: Path
) -> int:
    """Stream the original FASTQ(s) and write reads whose names are in ``needed``.

    ``needed`` maps mate-end ('1'/'2', or '' for single-end input) to the set of
    read names to pull from that file. Returns the count written.
    """
    written = 0
    file_for_end = {"": reads.left(), "1": reads.left(), "2": reads.right()}
    with open(out_fastq, "w") as out:
        for end, names in needed.items():
            if not names or not file_for_end.get(end):
                continue
            with pysam.FastxFile(file_for_end[end]) as fx:
                for rec in fx:
                    # `end` says which file this is, so it fixes the mate even
                    # when the FASTQ name carries no /1,/2 of its own.
                    name = canonical_name(rec.name, end)
                    if name in names:
                        qual = (
                            rec.quality
                            if rec.quality is not None
                            else "I" * len(rec.sequence)
                        )
                        out.write(f"@{name}\n{rec.sequence}\n+\n{qual}\n")
                        written += 1
    return written


@dataclass(frozen=True)
class _FastqRecord:
    """One FASTQ record retained by the step-4 planner."""

    name: str
    sequence: str
    quality: str


@dataclass(frozen=True)
class GenomeAlignmentPlan:
    """Observable summary of the RelocaTE2-compatible alignment inputs."""

    paired: int
    single_junctions: int
    support_reads: int
    retag: dict[str, tuple[str, str]]
    junction_bases: frozenset[str]


def _write_fastq_record(out, name: str, rec: _FastqRecord) -> None:
    out.write(f"@{name}\n{rec.sequence}\n+\n{rec.quality}\n")


def _read_flanking_records(
    flanking_files: list[str], all_te_hits: set[str]
) -> dict[str, dict[str, _FastqRecord]]:
    """Load classified flanks by pair base and validate the step-3 contract."""
    pairs: dict[str, dict[str, _FastqRecord]] = {}
    for fq in flanking_files:
        with pysam.FastxFile(fq) as fx:
            for rec in fx:
                original = strip_tag(rec.name)
                if original not in all_te_hits:
                    raise ValueError(
                        f"Classified read {rec.name!r} is absent from the all-TE-hit "
                        "artifact; rerun step 3 with this RelocaTE3 version"
                    )
                base, mate = split_mate(original)
                if mate not in ("1", "2"):
                    mate = ""
                    base = original
                by_end = pairs.setdefault(base, {})
                if mate in by_end:
                    raise ValueError(f"Duplicate classified read end: {rec.name}")
                quality = rec.quality or "I" * len(rec.sequence)
                by_end[mate] = _FastqRecord(rec.name, rec.sequence, quality)
    return pairs


def _fetch_original_records(
    reads: ReadLibrary, needed: dict[str, set[str]]
) -> dict[str, _FastqRecord]:
    """Fetch requested original reads in one streaming pass per mate file."""
    found: dict[str, _FastqRecord] = {}
    file_for_end = {"1": reads.left(), "2": reads.right()}
    for end, names in needed.items():
        read_file = file_for_end.get(end)
        if not names or not read_file:
            continue
        with pysam.FastxFile(read_file) as fx:
            for rec in fx:
                name = canonical_name(rec.name, end)
                if name in names:
                    quality = rec.quality or "I" * len(rec.sequence)
                    found[name] = _FastqRecord(name, rec.sequence, quality)
    return found


def _record_state(rec: _FastqRecord | None, name: str, all_te_hits: set[str]) -> str:
    """Return J, M, U, or N for one mate end."""
    if rec is not None:
        return "J" if _JUNCTION_TAG_RE.search(rec.name) else "M"
    return "U" if name in all_te_hits else "N"


def plan_genome_alignment_inputs(
    flanking_files: list[str],
    all_te_hits: set[str],
    reads: ReadLibrary,
    paired_r1_out: Path,
    paired_r2_out: Path,
    junction_se_out: Path,
    support_se_out: Path,
) -> GenomeAlignmentPlan:
    """Emit step-4 FASTQs using RelocaTE2's complete mate-state machine.

    The classified flanking FASTQs distinguish junction (J) and middle (M)
    reads. ``all_te_hits`` additionally identifies admitted but unclassified TE
    hits (U); only absence from that set is a genomic/no-TE-hit mate (N).
    Paired FASTQs always preserve the original R1/R2 order.
    """
    pairs = _read_flanking_records(flanking_files, all_te_hits)
    needed: dict[str, set[str]] = {"1": set(), "2": set()}

    if reads.is_paired:
        for base, records in pairs.items():
            states = {
                end: _record_state(records.get(end), f"{base}/{end}", all_te_hits)
                for end in ("1", "2")
            }
            for end, other in (("1", "2"), ("2", "1")):
                if states[end] in ("J", "M") and states[other] == "N":
                    needed[other].add(f"{base}/{other}")
    originals = _fetch_original_records(reads, needed)

    n_paired = n_junction = n_support = 0
    retag: dict[str, tuple[str, str]] = {}
    junction_bases: set[str] = set()
    with (
        open(paired_r1_out, "w") as paired_r1,
        open(paired_r2_out, "w") as paired_r2,
        open(junction_se_out, "w") as junction_se,
        open(support_se_out, "w") as support_se,
    ):
        for base in sorted(pairs):
            records = pairs[base]
            # A single-end input has no mate-state decision: keep junctions and
            # discard TE-internal reads, as RelocaTE2 does.
            if not reads.is_paired or "" in records:
                for rec in records.values():
                    if _JUNCTION_TAG_RE.search(rec.name):
                        _write_fastq_record(junction_se, rec.name, rec)
                        n_junction += 1
                        junction_bases.add(base)
                continue

            r1 = records.get("1")
            r2 = records.get("2")
            states = {
                "1": _record_state(r1, f"{base}/1", all_te_hits),
                "2": _record_state(r2, f"{base}/2", all_te_hits),
            }
            junction_ends = [end for end in ("1", "2") if states[end] == "J"]
            if junction_ends:
                junction_bases.add(base)

            if states == {"1": "J", "2": "J"}:
                assert r1 is not None and r2 is not None
                _write_fastq_record(paired_r1, f"{base}/1", r1)
                _write_fastq_record(paired_r2, f"{base}/2", r2)
                retag[base] = (r1.name, r2.name)
                n_paired += 1
                continue

            if len(junction_ends) == 1:
                junction_end = junction_ends[0]
                other = "2" if junction_end == "1" else "1"
                junction = records[junction_end]
                mate_name = f"{base}/{other}"
                mate = originals.get(mate_name) if states[other] == "N" else None
                if mate is None:
                    # J/M and J/U suppress the TE-hit mate; a missing N mate also
                    # falls back to the R2 single-junction behavior.
                    _write_fastq_record(junction_se, junction.name, junction)
                    n_junction += 1
                    continue

                ordered = {
                    junction_end: junction,
                    other: mate,
                }
                _write_fastq_record(paired_r1, f"{base}/1", ordered["1"])
                _write_fastq_record(paired_r2, f"{base}/2", ordered["2"])
                retag[base] = (ordered["1"].name, ordered["2"].name)
                n_paired += 1
                continue

            # With no junction, only M/N contributes: discard M and emit the
            # original N mate as single-end support. M/M and M/U emit nothing.
            for middle_end, other in (("1", "2"), ("2", "1")):
                if states[middle_end] == "M" and states[other] == "N":
                    mate = originals.get(f"{base}/{other}")
                    if mate is not None:
                        _write_fastq_record(support_se, mate.name, mate)
                        n_support += 1

    return GenomeAlignmentPlan(
        paired=n_paired,
        single_junctions=n_junction,
        support_reads=n_support,
        retag=retag,
        junction_bases=frozenset(junction_bases),
    )


def collect_junction_fullreads(
    flanking_files: list[str],
    reads: ReadLibrary,
    out_fastq: Path,
    mate_fastq: Path | None = None,
) -> int:
    """Write the full (untrimmed) sequences of junction reads to ``out_fastq``.

    These are re-aligned to the genome so the insertion finder can drop false
    junctions whose full read maps cleanly across the breakpoint. Junction
    membership comes from the tagged flanking FASTQs, not
    ``read_repeat_name.txt``: that table deliberately stores untagged original
    names for step 5.

    When ``mate_fastq`` is supplied for paired input, both original ends of
    every junction-containing pair are written in R1/R2 order. Mapping those
    files as a pair lets a unique genomic mate anchor a repetitive full
    junction read, matching RelocaTE2's ``matched.fullreads.bwa.mates`` path.
    Returns the number of junction full reads selected (not the number of mate
    records written).
    """
    needed: dict[str, set[str]] = {"": set(), "1": set(), "2": set()}
    for flank_path in flanking_files:
        with pysam.FastxFile(str(flank_path)) as records:
            for record in records:
                if not _JUNCTION_TAG_RE.search(record.name):
                    continue  # only 5'/3' junction reads, not :middle
                full_name = strip_tag(record.name)  # e.g. read_500_470/1
                _base, mate = split_mate(full_name)
                if mate in ("1", "2"):
                    needed[mate].add(full_name)
                elif reads.is_paired:
                    raise ValueError(
                        f"Paired junction read lacks a canonical mate suffix: "
                        f"{record.name!r}"
                    )
                else:
                    needed[""].add(full_name)

    if reads.is_paired and mate_fastq is not None:
        junction_count = len(needed["1"]) + len(needed["2"])
        pair_bases = {
            split_mate(name)[0] for end in ("1", "2") for name in needed[end]
        }
        pair_names = {
            "1": {f"{base}/1" for base in pair_bases},
            "2": {f"{base}/2" for base in pair_bases},
        }
        originals = _fetch_original_records(reads, pair_names)
        missing = sorted(
            name for names in pair_names.values() for name in names
            if name not in originals
        )
        if missing:
            preview = ", ".join(missing[:3])
            raise ValueError(
                f"Could not recover both original ends for {len(missing)} "
                f"junction mate records (first: {preview})"
            )
        with open(out_fastq, "w") as r1_out, open(mate_fastq, "w") as r2_out:
            for base in sorted(pair_bases):
                _write_fastq_record(r1_out, f"{base}/1", originals[f"{base}/1"])
                _write_fastq_record(r2_out, f"{base}/2", originals[f"{base}/2"])
        return junction_count

    return _fetch_reads_by_name(reads, needed, out_fastq)


def align_junction_fullreads(
    backend,
    genome: str,
    flanking_files: list[str],
    reads: ReadLibrary,
    out_bam: Path,
    threads: int,
    tmp: str,
) -> tuple[Path | None, int]:
    """Align original junction reads, mate-anchored when input is paired."""
    full_fq = Path(tmp) / f"{reads.name}.fullreads.R1.fq"
    mate_fq = (
        Path(tmp) / f"{reads.name}.fullreads.R2.fq" if reads.is_paired else None
    )
    n_full = collect_junction_fullreads(
        flanking_files, reads, full_fq, mate_fastq=mate_fq
    )
    fullreads_bam: Path | None = None
    if n_full > 0:
        fastqs = [str(full_fq)]
        if mate_fq is not None:
            fastqs.append(str(mate_fq))
        fullreads_bam = backend.map_genome(
            genome,
            fastqs,
            str(out_bam),
            paired=mate_fq is not None,
            threads=threads,
            tmpdir=tmp,
        )
    logger.info(
        "%s: %d full junction reads aligned for false-junction filtering",
        reads.name,
        n_full,
    )
    return fullreads_bam, n_full


def _restore_pair_names(
    in_bam: Path, out_bam: Path, retag: dict[str, tuple[str, str]]
) -> None:
    """Restore real R1/R2 names after alignment of neutral-name pair FASTQs."""
    with (
        pysam.AlignmentFile(str(in_bam), "rb") as inb,
        pysam.AlignmentFile(str(out_bam), "wb", template=inb) as outb,
    ):
        for rec in inb.fetch(until_eof=True):
            names = retag.get(split_mate(rec.query_name)[0])
            if names is not None:
                rec.query_name = names[0] if rec.is_read1 else names[1]
            outb.write(rec)


def align_flanks_anchored(
    backend,
    genome: str,
    flanking_files: list[str],
    all_te_hits: set[str],
    reads: ReadLibrary,
    out_bam,
    threads: int,
    tmp: str,
) -> Path:
    """Align RelocaTE2-compatible junction/support inputs to ``genome``.

    The paired and single-end inputs are planned from classified flanks plus the
    larger all-TE-hit population. Writes the merged coordinate-sorted BAM to
    ``out_bam``.

    Shared by :func:`align_to_genome` (pipeline path) and the ``align-genome`` CLI
    subcommand so both get mate-anchoring.
    """
    tmp_path = Path(tmp)
    sample = reads.name
    r1_fq = tmp_path / f"{sample}.flank_R1.fq"
    r2_fq = tmp_path / f"{sample}.flank_R2.fq"
    junction_fq = tmp_path / f"{sample}.junction_se.fq"
    support_fq = tmp_path / f"{sample}.support.fq"
    plan = plan_genome_alignment_inputs(
        flanking_files,
        all_te_hits,
        reads,
        r1_fq,
        r2_fq,
        junction_fq,
        support_fq,
    )
    logger.info(
        "%s: %d paired inputs, %d single-end junctions, %d support mates",
        sample,
        plan.paired,
        plan.single_junctions,
        plan.support_reads,
    )

    part_bams: list[str] = []
    if plan.paired > 0:
        raw_bam = tmp_path / f"{sample}.paired.raw.bam"
        backend.map_genome(
            genome,
            [str(r1_fq), str(r2_fq)],
            str(raw_bam),
            paired=True,
            threads=threads,
            tmpdir=tmp,
        )
        paired_bam = tmp_path / f"{sample}.paired.bam"
        _restore_pair_names(raw_bam, paired_bam, plan.retag)
        part_bams.append(str(paired_bam))
    se_inputs = [
        str(path)
        for path, count in (
            (junction_fq, plan.single_junctions),
            (support_fq, plan.support_reads),
        )
        if count > 0
    ]
    if se_inputs:
        se_bam = tmp_path / f"{sample}.se.bam"
        backend.map_genome(
            genome,
            se_inputs,
            str(se_bam),
            paired=False,
            threads=threads,
            tmpdir=tmp,
        )
        part_bams.append(str(se_bam))

    out_bam = Path(out_bam)
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    if not part_bams:
        raise RuntimeError("No junction or support reads remained after mate planning")
    if len(part_bams) == 1:
        # shutil.move (not Path.replace/os.rename) so the tmp -> out_bam move
        # works across filesystems (scratch tmp vs. network run dir): a rename
        # across devices raises OSError EXDEV.
        shutil.move(part_bams[0], str(out_bam))
    else:
        pysam.merge("-f", str(out_bam), *part_bams)
    pysam.index(str(out_bam))
    return out_bam


def align_to_genome(
    reads: ReadLibrary,
    genome: str,
    outdir: str | Path,
    threads: int = 1,
    genome_aligner: str = "bwaaln",
    genome_opts=(),
) -> tuple[Path, Path | None]:
    """Map trimmed flanking reads + supporting mates to ``genome``.

    Expects the trim step to have populated ``<outdir>/flanking`` plus the
    ``read_repeat_name`` and ``te_hit_names`` tables under ``te_containing``.
    Also aligns the full (untrimmed) junction reads to support false-junction
    filtering. Returns ``(genome_bam, fullreads_bam)``; the fullreads BAM is None
    if there are no junction reads.
    """
    outdir = Path(outdir)
    sample = reads.name
    backend = get_aligner(genome_aligner, threads, genome_opts=genome_opts)

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

    te_hits_path = outdir / "te_containing" / f"{sample}.te_hit_names.txt"
    if not te_hits_path.exists():
        raise FileNotFoundError(
            f"All-TE-hit table not found: {te_hits_path}; rerun step 3 with this "
            "RelocaTE3 version"
        )
    all_te_hits = read_te_hit_names(te_hits_path)

    fullreads_bam: Path | None = None
    outbam = genome_dir / f"{sample}.genome.bam"
    with tempfile.TemporaryDirectory() as tmp:
        align_flanks_anchored(
            backend, genome, flanking_files, all_te_hits, reads, outbam, threads, tmp
        )

        fullreads_bam, _ = align_junction_fullreads(
            backend,
            genome,
            flanking_files,
            reads,
            genome_dir / f"{sample}.fullreads.genome.bam",
            threads,
            tmp,
        )

    return outbam, fullreads_bam
