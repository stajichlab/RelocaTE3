"""Long-read subcommands: lr-align, lr-signatures and lr-cluster.

Each stage reads the previous stage's file and writes one file into
``--outdir``, named after ``--sample``, and appends one row to
``<outdir>/timing.tsv``:

- ``lr-align``: reads (minimap2 | samtools sort) or an existing BAM/CRAM
  (checked, then symlinked) -> ``<sample>.lr.bam`` (+ index), or
  ``<sample>.lr.cram`` (+ ``.crai``) for a CRAM.
- ``lr-signatures``: BAM/CRAM -> ``<sample>.signatures.tsv.gz``.
- ``lr-cluster``: signatures -> ``<sample>.clusters.tsv.gz``.

``--read-type`` exists only here; ``run-all``/``run-batch`` are unchanged.
Errors are raised as exceptions; :func:`RelocaTE3.cli.main` logs them and
exits 1. No stage overwrites an existing output.
"""

from __future__ import annotations

import argparse
import logging
import os
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import TypeVar

import pysam

from RelocaTE3.longread.align import (
    _BAM_INDEX_SUFFIXES, _CRAM_INDEX_SUFFIXES, align_reads, check_alignment_input,
)  # fmt: skip
from RelocaTE3.longread.cluster import ClusterParams, cluster_signatures
from RelocaTE3.longread.model import Reference, ReadType
from RelocaTE3.longread.placements import iter_placements
from RelocaTE3.longread.reference import LinearFastaReferenceSpace
from RelocaTE3.longread.signatures import SignatureParams, collect_signatures
from RelocaTE3.longread.tables import read_signatures, write_clusters, write_signatures
from RelocaTE3.longread.timing import StageCounts, stage_timer

log = logging.getLogger(__name__)

T = TypeVar("T")

READ_TYPES = [rt.value for rt in ReadType]


# --- helpers ----------------------------------------------------------------


def _require_file(path: Path, what: str) -> Path:
    if path.is_symlink() and not path.exists():
        raise FileNotFoundError(
            f"{what}: dangling symlink (target missing): {path} -> {os.readlink(path)}"
        )
    if not path.is_file():
        raise FileNotFoundError(f"{what} not found: {path}")
    return path


def _require_genome_for_cram(aln: Path, genome: str | Path | None) -> None:
    if aln.suffix.lower() == ".cram" and genome is None:
        raise ValueError(f"CRAM input needs --genome (the reference FASTA): {aln}")


def _refuse_existing(*paths: Path) -> None:
    """Raise if any path exists; a dangling symlink counts as existing."""
    for path in paths:
        if os.path.lexists(path):
            raise FileExistsError(f"refusing to overwrite {path}")


def _counted(items: Iterable[T], counts: StageCounts) -> Iterator[T]:
    """Yield ``items``, adding one to ``counts.records_in`` per item."""
    for item in items:
        counts.records_in += 1
        yield item


def _find_index(aln: Path, suffixes: tuple[str, ...]) -> tuple[Path, str]:
    """Return the index next to ``aln`` (as ``check_alignment_input`` finds it)."""
    for suffix in suffixes:
        for cand in (Path(str(aln) + suffix), aln.with_suffix(suffix)):
            if cand.is_file():
                return cand, suffix
    raise ValueError(f"{aln} has no {'/'.join(suffixes)} index")


def _mapped_records(aln: Path, genome: Path) -> int:
    """Mapped-record count from the index, or 0 if the index has none."""
    try:
        with pysam.AlignmentFile(str(aln), reference_filename=str(genome)) as fh:
            return int(fh.mapped)
    except (ValueError, OSError) as err:  # e.g. a CRAM index holds no counts
        log.debug("no mapped-record count for %s: %s", aln, err)
        return 0


# --- lr-align ---------------------------------------------------------------


def _link_alignment(bam: Path, genome: Path, reference_id: str, outdir: Path,
                    sample: str) -> Path:  # fmt: skip
    """Check a user BAM/CRAM and symlink it and its index into ``outdir``."""
    is_cram = bam.suffix.lower() == ".cram"
    space = LinearFastaReferenceSpace(Reference(reference_id, str(genome)))
    check_alignment_input(bam, space, reference_fasta=genome if is_cram else None)
    index, suffix = _find_index(
        bam, _CRAM_INDEX_SUFFIXES if is_cram else _BAM_INDEX_SUFFIXES
    )
    link = outdir / f"{sample}.lr.{'cram' if is_cram else 'bam'}"
    link_index = Path(str(link) + suffix)
    _refuse_existing(link, link_index)
    link.symlink_to(bam.resolve())
    try:
        link_index.symlink_to(index.resolve())
    except BaseException:
        link.unlink()
        raise
    log.info("linked %s -> %s (index %s)", link, bam.resolve(), link_index)
    return link


def cmd_lr_align(args: argparse.Namespace) -> int:
    """Align long reads, or check and link an existing alignment."""
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    genome = Path(args.genome)
    read_type = ReadType(args.read_type)
    with stage_timer(outdir, "lr-align", args.sample) as counts:
        _require_file(genome, "genome FASTA")
        if args.bam is not None:
            bam = _require_file(Path(args.bam), "alignment file")
            _require_genome_for_cram(bam, genome)
            counts.records_in = _mapped_records(bam, genome)
            out = _link_alignment(bam, genome, args.reference_id, outdir, args.sample)
        else:
            reads = [_require_file(Path(r), "read file") for r in args.reads]
            out = outdir / f"{args.sample}.lr.bam"
            _refuse_existing(out, Path(f"{out}.bai"))
            log.info("aligning %d read file(s) as %s with %d thread(s)",
                     len(reads), read_type, args.threads)  # fmt: skip
            align_reads(genome, reads, read_type, out, args.threads)
        counts.records_out = _mapped_records(out, genome)
    log.info("lr-align wrote %s (%d mapped records)", out, counts.records_out)
    return 0


# --- lr-signatures ----------------------------------------------------------


def cmd_lr_signatures(args: argparse.Namespace) -> int:
    """Extract insertion signatures from a long-read alignment."""
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"{args.sample}.signatures.tsv.gz"
    read_type = ReadType(args.read_type)
    params = SignatureParams(
        min_ins_len=args.min_ins_len, ins_merge_dist=args.ins_merge_dist,
        min_clip_len=args.min_clip_len, split_max_ref_gap=args.split_max_ref_gap,
    )  # fmt: skip
    placement_stats: dict[str, int] = {}
    signature_stats: dict[str, int] = {}
    with stage_timer(outdir, "lr-signatures", args.sample) as counts:
        bam = _require_file(Path(args.bam), "alignment file")
        _require_genome_for_cram(bam, args.genome)
        if args.genome is not None:
            _require_file(Path(args.genome), "genome FASTA")
        _refuse_existing(out)
        log.info("signature parameters: min_mapq=%d %s", args.min_mapq, params)
        placements = iter_placements(
            bam, args.reference_id, min_mapq=args.min_mapq, stats=placement_stats,
            reference_fasta=args.genome,
        )  # fmt: skip
        signatures = collect_signatures(
            _counted(placements, counts), params, signature_stats
        )
        counts.records_out = write_signatures(
            out, signatures, read_type=read_type, reference_id=args.reference_id
        )
    log.info("placement stats: %s", placement_stats)
    log.info("signature stats: %s", signature_stats)
    log.info("lr-signatures wrote %s (%d signatures from %d placements)",
             out, counts.records_out, counts.records_in)  # fmt: skip
    return 0


# --- lr-cluster -------------------------------------------------------------


def cmd_lr_cluster(args: argparse.Namespace) -> int:
    """Cluster insertion signatures into candidate sites."""
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"{args.sample}.clusters.tsv.gz"
    with stage_timer(outdir, "lr-cluster", args.sample) as counts:
        sig_path = _require_file(Path(args.signatures), "signatures file")
        _refuse_existing(out)
        meta, signatures = read_signatures(sig_path)
        try:
            read_type = ReadType(meta["read_type"])
        except ValueError:
            raise ValueError(
                f"{sig_path}: unknown read_type {meta['read_type']!r} "
                f"(expected one of {', '.join(READ_TYPES)})"
            ) from None
        window = args.cluster_window
        if window is None:
            window = read_type.default_cluster_window
            log.info("cluster window %d bp (default for %s)", window, read_type)
        else:
            log.info("cluster window %d bp (--cluster-window)", window)
        if window < 1:
            raise ValueError(f"--cluster-window must be at least 1, got {window}")
        params = ClusterParams(window=window, length_tolerance=args.length_tolerance)
        clusters = cluster_signatures(_counted(signatures, counts), params)
        counts.records_out = write_clusters(
            out, clusters, read_type=read_type, reference_id=meta["reference_id"]
        )
    log.info("lr-cluster wrote %s (%d clusters from %d signatures)",
             out, counts.records_out, counts.records_in)  # fmt: skip
    return 0


# --- parsers ----------------------------------------------------------------


def _add_sample_outdir(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--sample", required=True,
                        help="Sample name; prefixes every output file")  # fmt: skip
    parser.add_argument("-o", "--outdir", required=True,
                        help="Output directory (created if missing)")  # fmt: skip


def _add_read_type(parser: argparse.ArgumentParser, help: str) -> None:
    parser.add_argument("--read-type", required=True, choices=READ_TYPES, help=help)


def _add_reference_id(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--reference-id", default="ref",
                        help="Identifier of the reference, recorded on every "
                        "position")  # fmt: skip


def _menu_lr_align(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    _add_sample_outdir(parser)
    parser.add_argument("--genome", required=True,
                        help="Reference genome FASTA (also the CRAM reference)")  # fmt: skip
    _add_read_type(parser, "Long-read technology; selects the minimap2 preset")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--reads", nargs="+", metavar="FASTQ",
                        help="Long-read FASTQ file(s) to align")  # fmt: skip
    source.add_argument("--bam", metavar="BAM",
                        help="Existing coordinate-sorted, indexed BAM/CRAM "
                        "against --genome; checked and symlinked, not copied")  # fmt: skip
    parser.add_argument("--threads", type=int, default=1,
                        help="Threads for minimap2 and samtools sort")  # fmt: skip
    _add_reference_id(parser)
    parser.set_defaults(func=cmd_lr_align)
    return parser


def _menu_lr_signatures(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    _add_sample_outdir(parser)
    parser.add_argument("--bam", required=True,
                        help="Long-read BAM/CRAM (e.g. from lr-align)")  # fmt: skip
    _add_read_type(parser, "Long-read technology, recorded in the output")
    _add_reference_id(parser)
    parser.add_argument("--genome", default=None,
                        help="Reference genome FASTA; needed only for a CRAM")  # fmt: skip
    parser.add_argument("--min-mapq", type=int, default=20,
                        help="Lowest MAPQ of an alignment used")  # fmt: skip
    parser.add_argument("--min-ins-len", type=int, default=50,
                        help="Shortest insertion (INS or split-read gap) kept")  # fmt: skip
    parser.add_argument("--ins-merge-dist", type=int, default=20,
                        help="Largest reference distance between CIGAR "
                        "insertion pieces merged into one")  # fmt: skip
    parser.add_argument("--min-clip-len", type=int, default=100,
                        help="Shortest soft clip reported")  # fmt: skip
    parser.add_argument("--split-max-ref-gap", type=int, default=50,
                        help="Largest reference gap between two alignments of "
                        "a read that still counts as a split-read insertion")  # fmt: skip
    parser.set_defaults(func=cmd_lr_signatures)
    return parser


def _menu_lr_cluster(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    _add_sample_outdir(parser)
    parser.add_argument("--signatures", required=True,
                        help="Signatures file from lr-signatures")  # fmt: skip
    parser.add_argument("--cluster-window", type=int, default=None,
                        help="Largest distance (bp) joining signatures (default: "
                        "50 for hifi/ont-hq, 100 for ont/clr, from the file's "
                        "read type)")  # fmt: skip
    parser.add_argument("--length-tolerance", type=float, default=0.30,
                        help="Relative insert-length step that starts a new "
                        "cluster")  # fmt: skip
    parser.set_defaults(func=cmd_lr_cluster)
    return parser


def register_longread_commands(subparsers: argparse._SubParsersAction) -> None:
    """Add the lr-align, lr-signatures and lr-cluster subcommands."""
    # Imported here: RelocaTE3.cli imports this module while it is loading.
    from RelocaTE3.cli import CustomHelpFormatter, _add_common_args

    menus = [
        (_menu_lr_align, "lr-align",
         "Long reads, 1/3: align reads (or link a BAM/CRAM) to the genome",
         "Align long reads with minimap2 (preset from --read-type) into a "
         "sorted, indexed <outdir>/<sample>.lr.bam, or check an existing "
         "coordinate-sorted, indexed BAM/CRAM against --genome and symlink it "
         "(and its index) as <outdir>/<sample>.lr.bam or .lr.cram. Appends a "
         "row to <outdir>/timing.tsv."),
        (_menu_lr_signatures, "lr-signatures",
         "Long reads, 2/3: extract insertion signatures from the alignment",
         "Turn primary and supplementary alignments into insertion signatures "
         "(CIGAR insertions, split-read gaps, long soft clips) in "
         "<outdir>/<sample>.signatures.tsv.gz. Appends a row to "
         "<outdir>/timing.tsv."),
        (_menu_lr_cluster, "lr-cluster",
         "Long reads, 3/3: cluster signatures into candidate insertion sites",
         "Cluster the signatures of one sample into candidate insertion sites "
         "in <outdir>/<sample>.clusters.tsv.gz, keeping the input's read type "
         "and reference id. Appends a row to <outdir>/timing.tsv."),
    ]  # fmt: skip
    for menu, name, help, description in menus:
        parser = subparsers.add_parser(
            name, formatter_class=CustomHelpFormatter, help=help,
            description=description,
        )  # fmt: skip
        menu(parser)
        _add_common_args(parser)
