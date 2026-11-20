"""Long-read subcommands: lr-align, lr-signatures, lr-cluster, lr-annotate, lr-call.

Each stage reads the previous stage's file and writes one file into
``--outdir``, named after ``--sample``, and appends one row to
``<outdir>/timing.tsv``:

- ``lr-align``: reads (minimap2 | samtools sort) or an existing BAM/CRAM
  (checked, then symlinked) -> ``<sample>.lr.bam`` (+ index), or
  ``<sample>.lr.cram`` (+ ``.crai``) for a CRAM, plus the sidecar
  ``<sample>.lr.meta.tsv`` (``reference_id``, ``read_type``, ``genome``).
- ``lr-signatures``: BAM/CRAM -> ``<sample>.signatures.tsv.gz``. If the
  sidecar sits beside the BAM/CRAM (same name, last suffix replaced by
  ``.meta.tsv``), its reference id and read type must match the options.
- ``lr-cluster``: signatures -> ``<sample>.clusters.tsv.gz``.
- ``lr-annotate``: signatures + TE library (one minimap2 run) ->
  ``<sample>.te_hits.tsv.gz``.
- ``lr-call``: clusters + signatures + TE hits -> ``<sample>.calls.tsv.gz``.
  The three inputs must share read type and reference id, and every
  signature must have at least one TE-hit row.

The table outputs keep their input's read type and reference id.

``--read-type`` exists only here; ``run-all``/``run-batch`` are unchanged.
Errors are raised as exceptions; :func:`RelocaTE3.cli.main` logs them and
exits 1. No stage overwrites an existing output.
"""

from __future__ import annotations

import argparse
import logging
import os
from collections import Counter
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import TypeVar

import pysam

from RelocaTE3.longread.align import align_reads, check_alignment_input, find_index
from RelocaTE3.longread.annotate import TE_MINIMAP2_ARGS, annotate_signatures
from RelocaTE3.longread.call import (
    CallParams, SignatureFacts, SignatureTe, call_clusters, summarize_signature,
)  # fmt: skip
from RelocaTE3.longread.cluster import ClusterParams, cluster_signatures
from RelocaTE3.longread.model import Reference, ReadType, TeHit
from RelocaTE3.longread.placements import iter_placements
from RelocaTE3.longread.reference import LinearFastaReferenceSpace
from RelocaTE3.longread.signatures import SignatureParams, collect_signatures
from RelocaTE3.longread.tables import (
    read_clusters, read_signatures, read_te_hits, write_calls, write_clusters,
    write_signatures, write_te_hits,
)  # fmt: skip
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


def _read_type_of(meta: dict[str, str], path: Path) -> ReadType:
    """The read type recorded in a table's schema line."""
    try:
        return ReadType(meta["read_type"])
    except ValueError:
        raise ValueError(
            f"{path}: unknown read_type {meta['read_type']!r} "
            f"(expected one of {', '.join(READ_TYPES)})"
        ) from None


def _mapped_records(aln: Path, genome: Path) -> int:
    """Mapped-record count from the index, or 0 if the index has none."""
    try:
        with pysam.AlignmentFile(str(aln), reference_filename=str(genome)) as fh:
            return int(fh.mapped)
    except (ValueError, OSError) as err:  # e.g. a CRAM index holds no counts
        log.debug("no mapped-record count for %s: %s", aln, err)
        return 0


# --- reference-id sidecar ----------------------------------------------------

_META_KEYS = ("reference_id", "read_type", "genome")


def _meta_path(alignment: Path) -> Path:
    """Sidecar of an alignment: ``x.lr.bam`` / ``x.lr.cram`` -> ``x.lr.meta.tsv``."""
    return alignment.with_suffix(".meta.tsv")


def _meta_values(reference_id: str, read_type: ReadType,
                 genome: Path) -> dict[str, str]:  # fmt: skip
    """Sidecar values, checked to hold no tab or newline."""
    values = {"reference_id": reference_id, "read_type": read_type.value,
              "genome": str(genome.resolve())}  # fmt: skip
    for key, value in values.items():
        if any(c in value for c in "\t\n\r"):
            raise ValueError(f"sidecar {key} contains a tab or newline: {value!r}")
    return values


def _write_meta(path: Path, values: dict[str, str]) -> None:
    """Write checked sidecar values via ``<path>.partial`` and a rename."""
    partial = Path(f"{path}.partial")
    try:
        with partial.open("w") as fh:
            fh.write("key\tvalue\n")
            for key in _META_KEYS:
                fh.write(f"{key}\t{values[key]}\n")
        partial.replace(path)
    except BaseException:
        partial.unlink(missing_ok=True)
        raise
    log.info("wrote %s", path)


def _read_meta(path: Path) -> dict[str, str]:
    """Read a sidecar into a dict; the header must be ``key<TAB>value``."""
    with path.open() as fh:
        lines = fh.read().splitlines()
    if not lines or lines[0] != "key\tvalue":
        raise ValueError(f"{path}: not a long-read sidecar (header 'key\\tvalue')")
    meta: dict[str, str] = {}
    for n, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != 2:
            raise ValueError(f"{path}:{n}: expected 2 fields, got {len(fields)}")
        if fields[0] in meta:
            raise ValueError(f"{path}:{n}: duplicate key {fields[0]!r}")
        meta[fields[0]] = fields[1]
    return meta


def _check_meta(alignment: Path, reference_id: str, read_type: ReadType) -> None:
    """Compare the alignment's sidecar, if any, with the given options."""
    path = _meta_path(alignment)
    if not os.path.lexists(path):  # a dangling link or directory fails in _read_meta
        log.info("no sidecar %s (a user BAM/CRAM); --reference-id %r and "
                 "--read-type %s are not checked", path, reference_id,
                 read_type)  # fmt: skip
        return
    meta = _read_meta(path)
    for key, given, option in (("reference_id", reference_id, "--reference-id"),
                               ("read_type", read_type.value, "--read-type")):  # fmt: skip
        if key not in meta:
            raise ValueError(f"{path}: sidecar has no {key}")
        if meta[key] != given:
            raise ValueError(f"{option} {given!r} differs from {key} "
                             f"{meta[key]!r} recorded by lr-align in {path}")  # fmt: skip
    log.info("sidecar %s matches (reference_id %s, read_type %s)", path,
             reference_id, read_type)  # fmt: skip


# --- lr-align ---------------------------------------------------------------


def _link_alignment(bam: Path, genome: Path, reference_id: str, outdir: Path,
                    sample: str) -> Path:  # fmt: skip
    """Check a user BAM/CRAM and symlink it and its index into ``outdir``."""
    is_cram = bam.suffix.lower() == ".cram"
    space = LinearFastaReferenceSpace(Reference(reference_id, str(genome)))
    check_alignment_input(bam, space, reference_fasta=genome if is_cram else None)
    index = find_index(bam)
    if index is None:  # check_alignment_input already requires one
        raise ValueError(f"{bam} has no index")
    suffix = index.suffix
    link = outdir / f"{sample}.lr.{'cram' if is_cram else 'bam'}"
    link_index = Path(str(link) + suffix)
    _refuse_existing(link, link_index, _meta_path(link))
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
        meta_values = _meta_values(args.reference_id, read_type, genome)
        _require_file(genome, "genome FASTA")
        if args.bam is not None:
            bam = _require_file(Path(args.bam), "alignment file")
            _require_genome_for_cram(bam, genome)
            _check_meta(bam, args.reference_id, read_type)
            counts.records_in = _mapped_records(bam, genome)
            out = _link_alignment(bam, genome, args.reference_id, outdir, args.sample)
        else:
            reads = [_require_file(Path(r), "read file") for r in args.reads]
            out = outdir / f"{args.sample}.lr.bam"
            _refuse_existing(out, Path(f"{out}.bai"), _meta_path(out))
            log.info("aligning %d read file(s) as %s with %d thread(s)",
                     len(reads), read_type, args.threads)  # fmt: skip
            align_reads(genome, reads, read_type, out, args.threads)
        counts.records_out = _mapped_records(out, genome)
        _write_meta(_meta_path(out), meta_values)
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
        _check_meta(bam, args.reference_id, read_type)
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
        read_type = _read_type_of(meta, sig_path)
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


# --- lr-annotate ------------------------------------------------------------


def _require_fasta(path: Path) -> None:
    """Raise unless the first non-blank line of ``path`` starts with ``>``."""
    with path.open(errors="replace") as fh:
        for line in fh:
            if line.strip():
                if line.startswith(">"):
                    return
                raise ValueError(f"TE library is not FASTA (first line does not "
                                 f"start with '>'): {path}")  # fmt: skip
    raise ValueError(f"TE library is empty (no FASTA record): {path}")


def cmd_lr_annotate(args: argparse.Namespace) -> int:
    """Align every signature's sequence to the TE library."""
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"{args.sample}.te_hits.tsv.gz"
    stats: dict[str, int] = {}
    with stage_timer(outdir, "lr-annotate", args.sample) as counts:
        sig_path = _require_file(Path(args.signatures), "signatures file")
        library = _require_file(Path(args.te_library), "TE library FASTA")
        _require_fasta(library)
        _refuse_existing(out)
        meta, signatures = read_signatures(sig_path)
        read_type = _read_type_of(meta, sig_path)
        log.info("TE minimap2 flags: %s -t %d", " ".join(TE_MINIMAP2_ARGS),
                 args.threads)  # fmt: skip
        hits = annotate_signatures(
            _counted(signatures, counts), library, threads=args.threads,
            workdir=outdir, stats=stats,
        )  # fmt: skip
        counts.records_out = write_te_hits(
            out, hits, read_type=read_type, reference_id=meta["reference_id"]
        )
    log.info("annotation stats: %s", stats)
    if stats.get("queries", 0) > 0 and stats.get("with_hit", 0) == 0:
        log.warning("no signature hit the TE library %s (%d queries); check the "
                    "library", library, stats["queries"])  # fmt: skip
    log.info("lr-annotate wrote %s (%d hit rows for %d signatures)",
             out, counts.records_out, counts.records_in)  # fmt: skip
    return 0


# --- lr-call ----------------------------------------------------------------


def _check_same_origin(tables: list[tuple[Path, dict[str, str]]]) -> None:
    """Raise unless every table has the first one's read type and reference id."""
    first_path, first = tables[0]
    for path, meta in tables[1:]:
        for key in ("read_type", "reference_id"):
            if meta[key] != first[key]:
                raise ValueError(f"{key} differs between inputs: {first[key]!r} "
                                 f"in {first_path}, {meta[key]!r} in {path}")  # fmt: skip


def _summaries(hits: Iterable[TeHit], facts: dict[str, SignatureFacts],
               min_identity: float, path: Path) -> dict[str, SignatureTe]:  # fmt: skip
    """Group hits by signature id and summarize each signature."""
    by_sig: dict[str, list[TeHit]] = {}
    for hit in hits:
        if hit.signature_id not in facts:
            raise ValueError(f"{path}: hit for unknown signature "
                             f"{hit.signature_id!r}")  # fmt: skip
        by_sig.setdefault(hit.signature_id, []).append(hit)
    missing = [sid for sid in facts if sid not in by_sig]
    if missing:
        shown = ", ".join(repr(sid) for sid in missing[:5])
        raise ValueError(f"{len(missing)} signature(s) not annotated in {path} "
                         f"(no TE-hit row): {shown}")  # fmt: skip
    return {sid: summarize_signature(h, min_identity) for sid, h in by_sig.items()}


def cmd_lr_call(args: argparse.Namespace) -> int:
    """Call TE insertion sites from clusters, signatures and TE hits."""
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"{args.sample}.calls.tsv.gz"
    with stage_timer(outdir, "lr-call", args.sample) as counts:
        clu_path = _require_file(Path(args.clusters), "clusters file")
        sig_path = _require_file(Path(args.signatures), "signatures file")
        hit_path = _require_file(Path(args.te_hits), "TE hits file")
        _refuse_existing(out)
        clu_meta, clusters = read_clusters(clu_path)
        sig_meta, signatures = read_signatures(sig_path)
        hit_meta, hits = read_te_hits(hit_path)
        _check_same_origin([(sig_path, sig_meta), (clu_path, clu_meta),
                            (hit_path, hit_meta)])  # fmt: skip
        read_type = _read_type_of(sig_meta, sig_path)
        params = CallParams(
            min_support=(read_type.default_min_support if args.min_support is None
                         else args.min_support),
            min_te_identity=(read_type.default_min_te_identity
                             if args.min_te_identity is None
                             else args.min_te_identity),
            merge_window=(read_type.default_cluster_window
                          if args.merge_window is None else args.merge_window),
            min_te_reads_fraction=args.min_te_reads_fraction,
            min_family_fraction=args.min_family_fraction,
            min_te_coverage=args.min_te_coverage,
        )  # fmt: skip
        log.info("call parameters (%s): %s", read_type, params)
        facts = {s.signature_id: SignatureFacts.of(s) for s in signatures}
        te = _summaries(hits, facts, params.min_te_identity, hit_path)
        calls = call_clusters(_counted(clusters, counts), facts, te, params)
        counts.records_out = write_calls(
            out, calls, read_type=read_type, reference_id=sig_meta["reference_id"]
        )
    by_tier = Counter((str(c.tier), "PASS" if c.passed else "filtered") for c in calls)
    log.info("calls by tier: %s",
             ", ".join(f"{t} {p} {n}" for (t, p), n in sorted(by_tier.items())))  # fmt: skip
    log.info("lr-call wrote %s (%d calls from %d clusters, %d signatures)",
             out, counts.records_out, counts.records_in, len(facts))  # fmt: skip
    return 0


# --- parsers ----------------------------------------------------------------


def _add_sample_outdir(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--sample", required=True,
                        help="Sample name; prefixes every output file")  # fmt: skip
    parser.add_argument("-o", "--outdir", required=True,
                        help="Output directory (created if missing)")  # fmt: skip


def _add_read_type(parser: argparse.ArgumentParser, help_text: str) -> None:
    parser.add_argument("--read-type", required=True, choices=READ_TYPES,
                        help=help_text)  # fmt: skip


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


def _menu_lr_annotate(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    _add_sample_outdir(parser)
    parser.add_argument("--signatures", required=True,
                        help="Signatures file from lr-signatures")  # fmt: skip
    parser.add_argument("--te-library", required=True,
                        help="TE library FASTA; a hit's family is the first "
                        "word of the header")  # fmt: skip
    parser.add_argument("--threads", type=int, default=1,
                        help="Threads for minimap2")  # fmt: skip
    parser.set_defaults(func=cmd_lr_annotate)
    return parser


def _menu_lr_call(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    _add_sample_outdir(parser)
    parser.add_argument("--clusters", required=True,
                        help="Clusters file from lr-cluster")  # fmt: skip
    parser.add_argument("--signatures", required=True,
                        help="Signatures file the clusters were made from")  # fmt: skip
    parser.add_argument("--te-hits", required=True,
                        help="TE hits file from lr-annotate")  # fmt: skip
    parser.add_argument("--min-support", type=int, default=None,
                        help="Fewest distinct reads for a PASS call (default: 2 "
                        "for hifi/ont-hq, 3 for ont/clr)")  # fmt: skip
    parser.add_argument("--min-te-identity", type=float, default=None,
                        help="Lowest TE-hit identity used (default: 0.85 hifi, "
                        "0.80 ont-hq, 0.70 ont/clr)")  # fmt: skip
    parser.add_argument("--merge-window", type=int, default=None,
                        help="Largest distance (bp) merging same-family "
                        "clusters (default: the read type's cluster window)")  # fmt: skip
    parser.add_argument("--min-te-reads-fraction", type=float, default=0.80,
                        help="Lowest fraction of reads with a TE hit")  # fmt: skip
    parser.add_argument("--min-family-fraction", type=float, default=0.80,
                        help="Lowest fraction of hit reads of the top family")  # fmt: skip
    parser.add_argument("--min-te-coverage", type=float, default=0.80,
                        help="Lowest median TE coverage of a spanning call")  # fmt: skip
    parser.set_defaults(func=cmd_lr_call)
    return parser


def register_longread_commands(subparsers: argparse._SubParsersAction) -> None:
    """Add the lr-align, lr-signatures, lr-cluster, lr-annotate and lr-call subcommands."""
    # Function-local to avoid a circular import (RelocaTE3.cli imports this
    # module); _add_common_args stays private because RelocaTE3/cli.py is not
    # edited beyond registering these commands.
    from RelocaTE3.cli import CustomHelpFormatter, _add_common_args

    menus = [
        (_menu_lr_align, "lr-align",
         "Long reads, 1/5: align reads (or link a BAM/CRAM) to the genome",
         "Align long reads with minimap2 (preset from --read-type) into a "
         "sorted, indexed <outdir>/<sample>.lr.bam, or check an existing "
         "coordinate-sorted, indexed BAM/CRAM against --genome and symlink it "
         "(and its index) as <outdir>/<sample>.lr.bam or .lr.cram. Appends a "
         "row to <outdir>/timing.tsv."),
        (_menu_lr_signatures, "lr-signatures",
         "Long reads, 2/5: extract insertion signatures from the alignment",
         "Turn primary and supplementary alignments into insertion signatures "
         "(CIGAR insertions, split-read gaps, long soft clips) in "
         "<outdir>/<sample>.signatures.tsv.gz. Appends a row to "
         "<outdir>/timing.tsv."),
        (_menu_lr_cluster, "lr-cluster",
         "Long reads, 3/5: cluster signatures into candidate insertion sites",
         "Cluster the signatures of one sample into candidate insertion sites "
         "in <outdir>/<sample>.clusters.tsv.gz, keeping the input's read type "
         "and reference id. Appends a row to <outdir>/timing.tsv."),
        (_menu_lr_annotate, "lr-annotate",
         "Long reads, 4/5: annotate signatures against the TE library",
         "Align every signature's sequence to the TE library with one "
         "minimap2 run and write all hits (one no-hit row for a signature "
         "without a hit) to <outdir>/<sample>.te_hits.tsv.gz, keeping the "
         "input's read type and reference id. Appends a row to "
         "<outdir>/timing.tsv."),
        (_menu_lr_call, "lr-call",
         "Long reads, 5/5: call TE insertion sites from clusters",
         "Merge same-family clusters, vote the family and strand, and filter "
         "the calls into <outdir>/<sample>.calls.tsv.gz (PASS or the failed "
         "checks). The clusters, signatures and TE hits must share read "
         "type and reference id. Appends a row to <outdir>/timing.tsv."),
    ]  # fmt: skip
    for menu, name, help, description in menus:
        parser = subparsers.add_parser(
            name, formatter_class=CustomHelpFormatter, help=help,
            description=description,
        )  # fmt: skip
        menu(parser)
        _add_common_args(parser)
