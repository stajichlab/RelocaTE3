"""RelocaTE3 command-line interface (canonical implementation)."""

from __future__ import annotations

import argparse
import shlex
import re
import sys
import textwrap
from pathlib import Path

from RelocaTE3 import __author__, __entry_points__, __version__, logger
from RelocaTE3.aligners import GENOME_ALIGNERS, TE_ALIGNERS


class CustomHelpFormatter(argparse.HelpFormatter):
    """HelpFormatter that have customized function for text filling, line splitting and default parameter showing."""

    def _fill_text(self, text, width, indent):
        text = [
            self._whitespace_matcher.sub(" ", line).strip()
            for line in text.split("\n\n")
            if line != ""
        ]
        return "\n\n".join([textwrap.fill(line, width) for line in text])

    def _split_lines(self, text, width):
        text = [
            self._whitespace_matcher.sub(" ", line).strip()
            for line in text.split("\n")
            if line != ""
        ]
        formatted_text = []
        for line in text:
            formatted_text.extend(textwrap.wrap(line, width))
        return formatted_text

    def _get_help_string(self, action):
        help = action.help
        pattern = r"\(default: .+\)"
        if re.search(pattern, action.help) is None and action.default not in [
            argparse.SUPPRESS,
            None,
            False,
        ]:
            defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
            if action.option_strings or action.nargs in defaulting_nargs:
                help += " (default: %(default)s)"
        return help


def _add_common_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add --verbose and --version flags shared by all subcommands."""
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose mode for debug"
    )
    parser.add_argument("-V", "--version", action="version", version=__version__)
    return parser


def _menu_map(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Arguments for the 'map' subcommand."""
    parser.add_argument(
        "-l",
        "--left",
        "--r1",
        required=True,
        metavar="R1",
        help="Left/R1 read file (FASTQ)",
    )
    parser.add_argument(
        "-r",
        "--right",
        "--r2",
        metavar="R2",
        help="Right/R2 read file for paired-end (FASTQ)",
    )
    parser.add_argument(
        "-T",
        "--te-library",
        required=True,
        dest="te_library",
        metavar="FASTA",
        help="Transposon library FASTA file",
    )
    parser.add_argument(
        "-n",
        "--name",
        required=True,
        help="Sample/individual name (used as output prefix)",
    )
    parser.add_argument("-o", "--outdir", default=".", help="Output directory")
    parser.add_argument(
        "--threads", type=int, default=1, help="CPU threads for alignment"
    )
    parser.add_argument(
        "--te-aligner",
        "--aligner",
        dest="te_aligner",
        default="minimap2",
        choices=list(TE_ALIGNERS),
        help="Aligner for TE-library search (--aligner is a deprecated alias)",
    )
    parser.add_argument(
        "--te-opts",
        dest="te_opts",
        default="",
        help=(
            "Extra options passed verbatim to the TE-search aligner, e.g. "
            "\"-minIdentity=80\" (blat) or \"-B 4\" (minimap2)"
        ),
    )
    _add_common_args(parser)
    parser.set_defaults(func=cmd_map)
    return parser


def _menu_trim(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Arguments for the 'trim' subcommand."""
    parser.add_argument(
        "-b",
        "--bam",
        nargs="+",
        required=True,
        metavar="BAM",
        help="BAM file(s) from TE library alignment (output of 'map')",
    )
    parser.add_argument("-n", "--name", required=True, help="Sample/individual name")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory")
    parser.add_argument(
        "--min-match",
        type=int,
        default=10,
        dest="minimum_match_length",
        help="Minimum alignment match length to TE",
    )
    parser.add_argument(
        "--min-trimmed",
        type=int,
        default=10,
        dest="minimum_trimmed_length",
        help="Minimum trimmed flanking sequence length to retain",
    )
    parser.add_argument(
        "--mismatch",
        type=int,
        default=0,
        dest="mismatch_allowance",
        help="Allowed mismatches in TE alignment",
    )
    _add_common_args(parser)
    parser.set_defaults(func=cmd_trim)
    return parser


def _menu_run(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Arguments for the 'run' subcommand (full pipeline)."""
    parser.add_argument(
        "-l",
        "--left",
        "--r1",
        required=True,
        metavar="R1",
        help="Left/R1 read file (FASTQ)",
    )
    parser.add_argument(
        "-r",
        "--right",
        "--r2",
        metavar="R2",
        help="Right/R2 read file for paired-end (FASTQ)",
    )
    parser.add_argument(
        "-T",
        "--te-library",
        required=True,
        dest="te_library",
        metavar="FASTA",
        help="Transposon library FASTA file",
    )
    parser.add_argument(
        "-n",
        "--name",
        required=True,
        help="Sample/individual name (used as output prefix)",
    )
    parser.add_argument("-o", "--outdir", default=".", help="Output directory")
    parser.add_argument(
        "--threads", type=int, default=1, help="CPU threads for alignment"
    )
    parser.add_argument(
        "--te-aligner",
        "--aligner",
        dest="te_aligner",
        default="minimap2",
        choices=list(TE_ALIGNERS),
        help="Aligner for TE-library search (--aligner is a deprecated alias)",
    )
    parser.add_argument(
        "--te-opts",
        dest="te_opts",
        default="",
        help=(
            "Extra options passed verbatim to the TE-search aligner, e.g. "
            "\"-minIdentity=80\" (blat) or \"-B 4\" (minimap2)"
        ),
    )
    parser.add_argument(
        "--min-match",
        type=int,
        default=10,
        dest="minimum_match_length",
        help="Minimum alignment match length to TE",
    )
    parser.add_argument(
        "--min-trimmed",
        type=int,
        default=10,
        dest="minimum_trimmed_length",
        help="Minimum trimmed flanking sequence length to retain",
    )
    parser.add_argument(
        "--mismatch",
        type=int,
        default=0,
        dest="mismatch_allowance",
        help="Allowed mismatches in TE alignment",
    )
    _add_common_args(parser)
    parser.set_defaults(func=cmd_run)
    return parser


def _menu_characterize(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Arguments for the 'characterize' subcommand."""
    parser.add_argument(
        "-s",
        "--sites-file",
        required=True,
        dest="sites_file",
        metavar="TXT",
        help="RelocaTE non-reference insertion table (e.g. SAMPLE.mping.all_nonref.txt)",
    )
    parser.add_argument(
        "-b",
        "--bam",
        nargs="+",
        required=True,
        metavar="BAM_OR_CRAM",
        dest="bam",
        help="BAM or CRAM file(s) of original reads aligned to the reference genome (before TE trimming)",
    )
    parser.add_argument(
        "-g",
        "--genome-fasta",
        dest="genome_fasta",
        metavar="FASTA",
        help="Reference genome FASTA (required for CRAM input or --excision)",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default=None,
        help="Output directory (default: directory of the sites file)",
    )
    parser.add_argument(
        "-x",
        "--excision",
        action="store_true",
        help="Also search for excision events that leave a footprint",
    )
    parser.add_argument(
        "--samtools", default="samtools", help="Path to samtools executable"
    )
    parser.add_argument(
        "--bcftools", default="bcftools", help="Path to bcftools executable"
    )
    _add_common_args(parser)
    parser.set_defaults(func=cmd_characterize)
    return parser


def _menu_annotate_ref(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Arguments for the 'annotate-ref' subcommand (step 0)."""
    parser.add_argument(
        "-T",
        "--te-library",
        required=True,
        dest="te_library",
        metavar="FASTA",
        help="Transposon library FASTA file",
    )
    parser.add_argument(
        "-g",
        "--genome-fasta",
        required=True,
        dest="genome_fasta",
        metavar="FASTA",
        help="Reference genome FASTA file",
    )
    parser.add_argument("-o", "--outdir", default=".", help="Output directory")
    parser.add_argument(
        "--threads", type=int, default=1, help="CPU threads for minimap2"
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.8,
        dest="min_identity",
        help="Minimum alignment identity to keep an existing-TE copy",
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.8,
        dest="min_coverage",
        help="Minimum TE query coverage to keep an existing-TE copy",
    )
    _add_common_args(parser)
    parser.set_defaults(func=cmd_annotate_ref)
    return parser


def _menu_run_all(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Arguments for the 'run-all' subcommand (whole pipeline, one command)."""
    parser.add_argument(
        "-1", "--left", "--r1", required=True, dest="left", metavar="R1",
        help="Left/R1 read file (FASTQ, may be gzipped)",
    )
    parser.add_argument(
        "-2", "--right", "--r2", dest="right", metavar="R2",
        help="Right/R2 read file for paired-end (FASTQ, may be gzipped)",
    )
    parser.add_argument(
        "-T", "--te-library", required=True, dest="te_library", metavar="FASTA",
        help="Transposon library FASTA file",
    )
    parser.add_argument(
        "-g", "--genome-fasta", required=True, dest="genome_fasta", metavar="FASTA",
        help="Reference genome FASTA file (indexed automatically if needed)",
    )
    parser.add_argument("-n", "--name", required=True, help="Sample/individual name")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory")
    parser.add_argument(
        "--threads", type=int, default=1, help="CPU threads for alignment"
    )

    te = parser.add_argument_group("TE-library search (steps 2-3)")
    te.add_argument(
        "--te-aligner", dest="te_aligner", default="minimap2",
        choices=["minimap2", "bwa", "bwamem2", "bwaaln", "bowtie2", "blat"],
        help="Aligner for TE-library search",
    )
    te.add_argument(
        "--te-opts", dest="te_opts", default="",
        help='Extra options passed verbatim to the TE-search aligner. Values '
        'start with "-", so the --te-opts=VALUE form is required',
    )
    te.add_argument(
        "--min-match", type=int, default=10, dest="minimum_match_length",
        help="Minimum alignment match length to TE",
    )
    te.add_argument(
        "--min-trimmed", type=int, default=10, dest="minimum_trimmed_length",
        help="Minimum trimmed flanking sequence length to retain",
    )
    te.add_argument(
        "--mismatch", type=int, default=0, dest="mismatch_allowance",
        help="Allowed mismatches (TE alignment and read/genome comparison)",
    )

    gen = parser.add_argument_group("Genome placement (step 4)")
    gen.add_argument(
        "--genome-aligner", dest="genome_aligner", default="minimap2",
        choices=["minimap2", "bwa", "bwamem2", "bwaaln", "bowtie2"],
        help="Aligner for genome re-alignment (blat is not supported here)",
    )
    gen.add_argument(
        "--genome-opts", dest="genome_opts", default="",
        help="Extra options passed verbatim to the genome aligner "
        "(use the --genome-opts=VALUE form)",
    )

    call = parser.add_argument_group("Insertion calling (step 5)")
    call.add_argument(
        "--tsd", default="UNK",
        help="TSD motif (e.g. TTA), fixed wildcard (e.g. ...), or UNK to infer each TSD",
    )
    call.add_argument(
        "-c", "--target", default="ALL", help="Chromosome to analyze, or ALL"
    )
    call.add_argument(
        "--te-name", dest="te_name", default="repeat",
        help="TE label used in output filenames",
    )
    call.add_argument(
        "--min-mapq", type=int, default=1, dest="min_mapq",
        help="Minimum MAPQ for a uniquely-mapped read",
    )
    call.add_argument(
        "--require-both-junctions", action="store_true", dest="require_both_junctions",
        help="Emit only insertions with both a left and right junction read "
        "(RelocaTE2 parity, fewer false positives)",
    )

    ref = parser.add_argument_group("Reference TEs (steps 0/6)")
    ref.add_argument(
        "--repeatmasker", metavar="OUT", default=None,
        help="RepeatMasker .out for the reference. Filters novel calls that land "
        "on a known reference TE, and additionally emits all_ref_insert.{gff,txt}",
    )

    gt = parser.add_argument_group("Genotyping (step 7)")
    gt.add_argument(
        "--genotype", action="store_true",
        help="Classify each insertion as homozygous/heterozygous/somatic. Aligns "
        "the original reads to the genome unless --genotype-bam is given",
    )
    gt.add_argument(
        "--genotype-bam", dest="genotype_bam", default=None, metavar="BAM",
        help="Existing BAM/CRAM of the original reads aligned to the genome, "
        "reused for genotyping instead of realigning (implies --genotype)",
    )
    gt.add_argument(
        "-x", "--excision", action="store_true",
        help="Also search for excision events that leave a footprint",
    )
    gt.add_argument(
        "--bcftools", default="bcftools", help="Path to bcftools executable"
    )

    _add_common_args(parser)
    parser.set_defaults(func=cmd_run_all)
    return parser


def _menu_find_reference(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Arguments for the 'find-reference' subcommand (steps 0/6)."""
    parser.add_argument(
        "-b",
        "--bam",
        required=True,
        help="Sorted BAM of flanking reads aligned to the genome (output of align-genome)",
    )
    parser.add_argument(
        "-R",
        "--read-repeat",
        required=True,
        dest="read_repeat",
        metavar="TXT",
        help="read_repeat_name table from the trim step (read -> TE mapping)",
    )
    parser.add_argument(
        "--repeatmasker",
        required=True,
        metavar="OUT",
        help="RepeatMasker .out annotation of reference TE copies",
    )
    parser.add_argument("-n", "--name", required=True, help="Sample/experiment name")
    parser.add_argument(
        "-o", "--outdir", default=".", help="Output directory (writes results/ subdir)"
    )
    _add_common_args(parser)
    parser.set_defaults(func=cmd_find_reference)
    return parser


def _menu_index_genome(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Arguments for the 'index-genome' subcommand (step 1)."""
    parser.add_argument(
        "-g",
        "--genome-fasta",
        required=True,
        dest="genome_fasta",
        metavar="FASTA",
        help="Reference genome FASTA file",
    )
    parser.add_argument(
        "--force", action="store_true", help="Recreate indexes even if present"
    )
    _add_common_args(parser)
    parser.set_defaults(func=cmd_index_genome)
    return parser


def _menu_align_genome(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Arguments for the 'align-genome' subcommand (step 4)."""
    parser.add_argument(
        "-g",
        "--genome-fasta",
        required=True,
        dest="genome_fasta",
        metavar="FASTA",
        help="Reference genome FASTA file",
    )
    parser.add_argument(
        "-f",
        "--fastq",
        nargs="+",
        required=True,
        metavar="FASTQ",
        dest="fastq",
        help="Trimmed flanking-read FASTQ file(s) from the trim step",
    )
    parser.add_argument(
        "-n", "--name", required=True, help="Sample/individual name (output prefix)"
    )
    parser.add_argument("-o", "--outdir", default=".", help="Output directory")
    parser.add_argument(
        "-1",
        "--left",
        "--r1",
        default=None,
        metavar="R1",
        help="Original left/R1 reads (FASTQ). When given (with --right), junction "
        "flanks are aligned paired with their genomic mate so ambiguous flanks are "
        "anchored to the true insertion (reads read_repeat from --outdir).",
    )
    parser.add_argument(
        "-2",
        "--right",
        "--r2",
        default=None,
        metavar="R2",
        help="Original right/R2 reads (FASTQ), the mate file for --left.",
    )
    parser.add_argument(
        "--paired",
        action="store_true",
        help="Align the first two FASTQs as a read pair (legacy; prefer -1/-2 for "
        "mate-anchored flank alignment)",
    )
    parser.add_argument(
        "--genome-aligner",
        dest="genome_aligner",
        default="minimap2",
        choices=list(GENOME_ALIGNERS),
        help="Aligner for genome re-alignment (blat is not supported here)",
    )
    parser.add_argument(
        "--genome-opts",
        dest="genome_opts",
        default="",
        help=(
            "Extra options passed verbatim to the genome aligner, e.g. "
            "\"-n 0.10\" (bwa aln edit-distance budget)"
        ),
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="CPU threads for alignment"
    )
    _add_common_args(parser)
    parser.set_defaults(func=cmd_align_genome)
    return parser


def _menu_find_insertions(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Arguments for the 'find-insertions' subcommand (step 5)."""
    parser.add_argument(
        "-b",
        "--bam",
        required=True,
        metavar="BAM",
        help="Sorted BAM of flanking reads aligned to the genome (output of align-genome)",
    )
    parser.add_argument(
        "-R",
        "--read-repeat",
        required=True,
        dest="read_repeat",
        metavar="TXT",
        help="read_repeat_name table from the trim step (read -> TE mapping)",
    )
    parser.add_argument(
        "--tsd",
        required=True,
        help="TSD motif (e.g. TTA), fixed wildcard (e.g. ...), or UNK to infer each TSD",
    )
    parser.add_argument(
        "-c", "--target", default="ALL", help="Chromosome to analyze, or ALL"
    )
    parser.add_argument("-n", "--name", required=True, help="Sample/experiment name")
    parser.add_argument(
        "-o", "--outdir", default=".", help="Output directory (writes results/ subdir)"
    )
    parser.add_argument(
        "--te-name",
        default="repeat",
        dest="te_name",
        help="TE label used in output filenames",
    )
    parser.add_argument(
        "--reference-ins",
        dest="reference_ins",
        metavar="FILE",
        help="Existing-TE RepeatMasker .out or BED to skip reference insertions",
    )
    parser.add_argument(
        "--mismatch",
        type=int,
        default=0,
        dest="mismatch_allow",
        help="Allowed read/genome mismatches (excluding indels)",
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=1,
        dest="min_mapq",
        help="Minimum MAPQ for a uniquely-mapped read",
    )
    parser.add_argument(
        "--require-both-junctions",
        action="store_true",
        dest="require_both_junctions",
        help="Emit only insertions with both a left and right junction read "
        "(drop single-sided calls; RelocaTE2 parity, fewer false positives)",
    )
    _add_common_args(parser)
    parser.set_defaults(func=cmd_find_insertions)
    return parser


# ---------------------------------------------------------------------------
# Command handlers
# ---------------------------------------------------------------------------


def cmd_map(args: argparse.Namespace) -> int:
    """Align reads to TE library and write BAM files."""
    from RelocaTE3.ReadLibrary import ReadLibrary

    fileset = [args.left] + ([args.right] if args.right else [])
    reads = ReadLibrary(fileset, args.name)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    from RelocaTE3.aligners import get_aligner

    backend = get_aligner(
        args.te_aligner,
        args.threads,
        te_opts=split_aligner_opts(getattr(args, "te_opts", "")),
    )
    backend.index(args.te_library)
    bamfiles = backend.map_te_library(reads, args.te_library, out)
    logger.info("%d BAM file(s) written to %s", len(bamfiles), args.outdir)
    return 0


def cmd_trim(args: argparse.Namespace) -> int:
    """Trim TE sequence from TE-library BAM reads and emit junction-named flanking FASTQs."""
    from RelocaTE3.librelocate import RelocaTE

    bam_paths = [Path(b) for b in args.bam]
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    relocate = RelocaTE(verbose=int(args.verbose))
    directions = relocate._bam_directions(bam_paths)
    flank_written = relocate.write_trimmed_reads(
        args.name,
        list(zip(directions, bam_paths)),
        out,
        minimum_match_length=args.minimum_match_length,
        minimum_trimmed_length=args.minimum_trimmed_length,
        mismatch_allowance=args.mismatch_allowance,
    )
    logger.info("Wrote %d flanking read(s) to %s", flank_written, out / "flanking")
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    """Identify TE-containing reads and generate flanking reads (map + trim).

    Maps reads to the TE library then trims the TE-matching portion. This is
    steps 2-3 (TE-read identification and flank generation), not the complete
    insertion-calling pipeline.
    """
    from RelocaTE3.librelocate import RelocaTE
    from RelocaTE3.ReadLibrary import ReadLibrary

    fileset = [args.left] + ([args.right] if args.right else [])
    reads = ReadLibrary(fileset, args.name)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    relocate = RelocaTE(
        TElib=args.te_library, threads=args.threads, verbose=int(args.verbose)
    )
    n = relocate.identify_TE_reads(
        reads,
        out,
        te_aligner=args.te_aligner,
        te_opts=split_aligner_opts(getattr(args, "te_opts", "")),
        len_cut_match=args.minimum_match_length,
        len_cut_trim=args.minimum_trimmed_length,
        mismatch_allowance=args.mismatch_allowance,
    )
    logger.info("%d read(s) written", n)
    return 0


def cmd_characterize(args: argparse.Namespace) -> int:
    """Characterize insertion sites as homozygous, heterozygous, somatic or excised."""
    from RelocaTE3.characterize import Characterizer

    characterizer = Characterizer(
        samtools=args.samtools, bcftools=args.bcftools, verbose=int(args.verbose)
    )
    bam_paths = [Path(b) for b in args.bam]
    txt_path, gff_path = characterizer.characterize(
        sites_file=Path(args.sites_file),
        bam_files=bam_paths,
        genome_fasta=Path(args.genome_fasta) if args.genome_fasta else None,
        outdir=Path(args.outdir) if args.outdir else None,
        excision=args.excision,
    )
    logger.info("Characterization written to %s and %s", txt_path, gff_path)
    return 0


def cmd_annotate_ref(args: argparse.Namespace) -> int:
    """Annotate existing/reference TE copies in the genome with minimap2 (step 0)."""
    from RelocaTE3.reference_te import ReferenceTEAnnotator

    annotator = ReferenceTEAnnotator(threads=args.threads, verbose=int(args.verbose))
    bed = annotator.annotate_minimap(
        te_library=Path(args.te_library),
        genome_fasta=Path(args.genome_fasta),
        outdir=Path(args.outdir),
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
    )
    logger.info("Existing-TE annotation written to %s", bed)
    return 0


def cmd_run_all(args: argparse.Namespace) -> int:
    """Run the whole pipeline for one sample with a single command.

    This is deliberately a thin orchestrator over the *same* ``cmd_*`` handlers
    the staged workflow uses (``index-genome`` -> ``run`` -> ``align-genome`` ->
    ``find-insertions`` -> optionally ``find-reference`` and ``characterize``),
    rather than a second implementation.

    That matters: RelocaTE3 has two non-equivalent insertion callers -- the
    CLI's ``InsertionFinder`` (the path validated against RelocaTE2 on the
    benchmark) and ``pipeline.run_sample``'s module-level ``find_insertions``.
    Wrapping the latter would give a one-command mode that quietly produces
    different numbers from the documented workflow, so we compose the former.
    """
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    verbose = args.verbose

    def _step(label, handler, /, **kwargs):
        # positional-only: the steps below pass name=/outdir=/bam= through as
        # kwargs, which would otherwise collide with these parameter names.
        logger.info("run-all: %s", label)
        rc = handler(argparse.Namespace(verbose=verbose, **kwargs))
        if rc != 0:
            raise RuntimeError(f"run-all: step '{label}' failed with exit code {rc}")

    # Step 1 -- index the genome unless it is already indexed. The staged
    # workflow expects prebuilt indexes; a one-command mode should not.
    genome = Path(args.genome_fasta)
    indexed = (
        genome.with_suffix(genome.suffix + ".fai").exists()
        and genome.with_suffix(genome.suffix + ".mmi").exists()
    )
    if indexed:
        logger.info("run-all: genome already indexed, skipping index-genome")
    else:
        _step(
            "index-genome",
            cmd_index_genome,
            genome_fasta=args.genome_fasta,
            force=False,
        )

    # Steps 2-3 -- TE-library search and trimming.
    _step(
        "run",
        cmd_run,
        left=args.left,
        right=args.right,
        te_library=args.te_library,
        name=args.name,
        outdir=str(outdir),
        threads=args.threads,
        te_aligner=args.te_aligner,
        te_opts=args.te_opts,
        minimum_match_length=args.minimum_match_length,
        minimum_trimmed_length=args.minimum_trimmed_length,
        mismatch_allowance=args.mismatch_allowance,
    )

    read_repeat = outdir / "te_containing" / f"{args.name}.read_repeat_name.txt"
    if not read_repeat.is_file() or read_repeat.stat().st_size == 0:
        raise RuntimeError(
            f"run-all: no TE-containing reads were found (missing/empty {read_repeat}). "
            "Check that the TE library matches the reads."
        )

    flanking = sorted(
        str(p) for p in (outdir / "flanking").glob(f"{args.name}.*.flankingReads.fq")
    )
    if not flanking:
        raise RuntimeError(
            f"run-all: no flanking FASTQs were produced under {outdir / 'flanking'}"
        )

    # Step 4 -- place flanks on the genome. Passing the original reads enables
    # mate-anchored alignment, which is what recovers ambiguous short flanks.
    _step(
        "align-genome",
        cmd_align_genome,
        genome_fasta=args.genome_fasta,
        fastq=flanking,
        name=args.name,
        outdir=str(outdir),
        left=args.left,
        right=args.right,
        paired=False,
        genome_aligner=args.genome_aligner,
        genome_opts=args.genome_opts,
        threads=args.threads,
    )

    genome_bam = _find_genome_bam(outdir, args.name)

    # Step 5 -- cluster junctions into non-reference insertion calls.
    _step(
        "find-insertions",
        cmd_find_insertions,
        bam=str(genome_bam),
        read_repeat=str(read_repeat),
        tsd=args.tsd,
        target=args.target,
        name=args.name,
        outdir=str(outdir),
        te_name=args.te_name,
        reference_ins=args.repeatmasker,
        mismatch_allow=args.mismatch_allowance,
        min_mapq=args.min_mapq,
        require_both_junctions=args.require_both_junctions,
    )

    # Steps 0/6 -- reference/shared calls (RelocaTE2's all_ref_insert.*).
    if args.repeatmasker:
        _step(
            "find-reference",
            cmd_find_reference,
            bam=str(genome_bam),
            read_repeat=str(read_repeat),
            repeatmasker=args.repeatmasker,
            name=args.name,
            outdir=str(outdir),
        )

    # Step 7 -- zygosity, which needs the ORIGINAL (untrimmed) reads on the
    # genome. Build that BAM unless the user supplied one.
    if args.genotype or args.genotype_bam:
        sites = outdir / "results" / f"{args.target}.{args.te_name}.all_nonref_insert.txt"
        if not sites.is_file():
            raise RuntimeError(f"run-all: expected insertion table not found: {sites}")

        genotype_bam = args.genotype_bam
        if not genotype_bam:
            from RelocaTE3.align import Aligner
            from RelocaTE3.ReadLibrary import ReadLibrary

            fileset = [args.left] + ([args.right] if args.right else [])
            reads = ReadLibrary(fileset, args.name)
            genotype_bam = str(outdir / "genome_aln" / f"{args.name}.reads.genome.bam")
            logger.info("run-all: aligning original reads to the genome for genotyping")
            Aligner(args.threads).map_library_to_genome(
                args.genome_fasta, reads, genotype_bam, cpu_threads=args.threads
            )

        _step(
            "characterize",
            cmd_characterize,
            sites_file=str(sites),
            bam=[genotype_bam],
            genome_fasta=args.genome_fasta,
            outdir=str(outdir / "results"),
            excision=args.excision,
            samtools="samtools",
            bcftools=args.bcftools,
        )

    logger.info("run-all: complete; results under %s", outdir / "results")
    return 0


def _find_genome_bam(outdir: Path, sample: str) -> Path:
    """Locate the genome BAM written by align-genome.

    The filename embeds the aligner (e.g. ``.repeat.bwaaln.sorted.bam``), so it
    is discovered rather than reconstructed.
    """
    matches = sorted(outdir.glob(f"{sample}.repeat.*.sorted.bam"))
    if not matches:
        raise RuntimeError(
            f"run-all: align-genome produced no BAM matching "
            f"{sample}.repeat.*.sorted.bam in {outdir}"
        )
    return matches[0]


def cmd_find_reference(args: argparse.Namespace) -> int:
    """Call reference/shared TE insertions from junction reads (steps 0/6).

    Distinct from ``annotate-ref``: that only records where reference TE copies
    ARE (an ``existingTE.bed`` used to filter novel calls). This command asks
    which of those copies are also supported by reads in *this sample*, and
    emits them as RelocaTE2's ``all_ref_insert.{gff,txt}``.
    """
    from RelocaTE3.genome_align import read_read_repeat
    from RelocaTE3.insertions import write_insertions_gff, write_insertions_txt
    from RelocaTE3.reference_te import (
        find_reference_insertions,
        write_existing_te_bed_from_rm,
    )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    result_dir = outdir / "results"
    result_dir.mkdir(parents=True, exist_ok=True)

    bed_path = outdir / "existingTE.bed"
    reference_tes = write_existing_te_bed_from_rm(args.repeatmasker, bed_path)
    logger.info("Wrote %d reference TE copies to %s", len(reference_tes), bed_path)

    read_repeat = read_read_repeat(Path(args.read_repeat))
    insertions = find_reference_insertions(
        str(args.bam), read_repeat, reference_tes
    )

    gff_path = result_dir / f"{args.name}.all_ref_insert.gff"
    txt_path = result_dir / f"{args.name}.all_ref_insert.txt"
    write_insertions_gff(insertions, gff_path, args.name)
    write_insertions_txt(insertions, txt_path)
    logger.info(
        "Wrote %d reference/shared insertions to %s and %s",
        len(insertions),
        txt_path,
        gff_path,
    )
    return 0


def cmd_index_genome(args: argparse.Namespace) -> int:
    """Index/format the reference genome (samtools faidx + minimap2 index) (step 1)."""
    from RelocaTE3.align import Aligner

    aln = Aligner()
    aln.verbose = bool(args.verbose)
    aln.index_genome(args.genome_fasta, force=args.force)
    logger.info("Indexed genome %s", args.genome_fasta)
    return 0


def cmd_align_genome(args: argparse.Namespace) -> int:
    """Align trimmed flanking reads to the reference genome (step 4)."""
    if args.left:
        # Mate-anchored path: pair each junction flank with its genomic mate so
        # ambiguous (multi-mapping) flanks are anchored to the true insertion.
        import tempfile

        from RelocaTE3.aligners import get_aligner
        from RelocaTE3.genome_align import (
            align_flanks_anchored,
            read_read_repeat,
        )
        from RelocaTE3.ReadLibrary import ReadLibrary

        reads = ReadLibrary(
            [args.left] + ([args.right] if args.right else []), args.name
        )
        outdir = Path(args.outdir)
        rr_path = outdir / "te_containing" / f"{args.name}.read_repeat_name.txt"
        read_repeat = read_read_repeat(rr_path) if rr_path.exists() else {}
        suffix = "minimap" if args.genome_aligner == "minimap2" else args.genome_aligner
        out_bam = outdir / f"{args.name}.repeat.{suffix}.sorted.bam"
        backend = get_aligner(
            args.genome_aligner,
            args.threads,
            genome_opts=split_aligner_opts(getattr(args, "genome_opts", "")),
        )
        backend.index(args.genome_fasta)
        with tempfile.TemporaryDirectory() as tmp:
            bam = align_flanks_anchored(
                backend,
                args.genome_fasta,
                list(args.fastq),
                read_repeat,
                reads,
                out_bam,
                args.threads,
                tmp,
            )
        logger.info("Genome-aligned BAM (mate-anchored) written to %s", bam)
        return 0
    if args.genome_aligner == "minimap2":
        # Preserve the original minimap2 path byte-for-byte (tuned flags, paired
        # handling, and the {name}.repeat.minimap.sorted.bam output name).
        from RelocaTE3.align import Aligner

        aln = Aligner(threads=args.threads)
        aln.verbose = bool(args.verbose)
        bam = aln.map_genome_minimap(
            genome=args.genome_fasta,
            fastqs=args.fastq,
            name=args.name,
            outdir=args.outdir,
            paired=args.paired,
        )
    else:
        from RelocaTE3.aligners import get_aligner

        out_bam = (
            Path(args.outdir) / f"{args.name}.repeat.{args.genome_aligner}.sorted.bam"
        )
        backend = get_aligner(
            args.genome_aligner,
            args.threads,
            genome_opts=split_aligner_opts(getattr(args, "genome_opts", "")),
        )
        backend.index(args.genome_fasta)
        bam = backend.map_genome(
            args.genome_fasta,
            args.fastq,
            out_bam,
            paired=args.paired,
            threads=args.threads,
        )
    logger.info("Genome-aligned BAM written to %s", bam)
    return 0


def cmd_find_insertions(args: argparse.Namespace) -> int:
    """Find non-reference insertions from genome-aligned flanking reads (step 5)."""
    from RelocaTE3.insertions import InsertionFinder

    finder = InsertionFinder(
        mismatch_allow=args.mismatch_allow,
        min_mapq=args.min_mapq,
        verbose=int(args.verbose),
        require_both_junctions=getattr(args, "require_both_junctions", False),
    )
    out_txt = finder.find_insertions(
        bam_file=Path(args.bam),
        read_repeat_file=Path(args.read_repeat),
        tsd=args.tsd,
        target=args.target,
        sample=args.name,
        outdir=Path(args.outdir),
        te_name=args.te_name,
        reference_ins=Path(args.reference_ins) if args.reference_ins else None,
    )
    logger.info("Non-reference insertions written to %s", out_txt)
    return 0


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------



def split_aligner_opts(value) -> list[str]:
    """Split a --te-opts/--genome-opts string into argv tokens.

    The value is passed through verbatim to the aligner, so it is split with
    shell-like quoting rather than a naive ``str.split``. ``None``/empty yields
    no extra options.
    """
    if not value:
        return []
    if isinstance(value, (list, tuple)):
        return [str(v) for v in value]
    return shlex.split(value)


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level argument parser with all subcommands registered."""
    prog = __entry_points__.get(__name__, "relocaTE3")

    parser = argparse.ArgumentParser(
        prog=prog,
        formatter_class=CustomHelpFormatter,
        description=main.__doc__,
        epilog=f"Written by {__author__}",
    )
    parser.add_argument("-V", "--version", action="version", version=__version__)

    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")
    subparsers.required = True

    _menu_map(
        subparsers.add_parser(
            "map",
            formatter_class=CustomHelpFormatter,
            help="Align reads to TE library (produces BAM files)",
            description="Align FASTQ reads to a transposon sequence library using minimap2 or bwa.",
        )
    )
    _menu_trim(
        subparsers.add_parser(
            "trim",
            formatter_class=CustomHelpFormatter,
            help="Trim TE sequences from TE-library-aligned BAM files",
            description="Process BAM files from 'map' to identify and trim transposon sequences from reads.",
        )
    )
    _menu_run(
        subparsers.add_parser(
            "run",
            formatter_class=CustomHelpFormatter,
            help="Identify TE-containing reads and generate flanking reads (map + trim)",
            description="Map reads to the TE library then trim the TE-matching "
            "portion, emitting flanking reads and the read_repeat_name table. "
            "This is TE-read identification and flank generation (steps 2-3), "
            "NOT the complete insertion-calling pipeline.",
        )
    )
    _menu_annotate_ref(
        subparsers.add_parser(
            "annotate-ref",
            formatter_class=CustomHelpFormatter,
            help="Annotate existing TE copies in the reference genome (step 0)",
            description="Align the TE library to the reference genome with minimap2 to record existing "
            "(reference) transposon copies, so novel-insertion calling can skip them.",
        )
    )
    _menu_run_all(
        subparsers.add_parser(
            "run-all",
            formatter_class=CustomHelpFormatter,
            help="Run the whole pipeline for one sample (one command)",
            description="Run every step for a single sample: index-genome -> run "
            "(map + trim) -> align-genome -> find-insertions, plus find-reference "
            "when --repeatmasker is given and characterize when --genotype is. "
            "Dispatches the same handlers as the staged subcommands, so results "
            "match the staged workflow exactly.",
        )
    )
    _menu_find_reference(
        subparsers.add_parser(
            "find-reference",
            formatter_class=CustomHelpFormatter,
            help="Call reference/shared TE insertions present in this sample (steps 0/6)",
            description="Call TE copies that are annotated in the reference AND supported by "
            "junction reads in this sample, emitting RelocaTE2's all_ref_insert.{gff,txt}. "
            "Unlike 'annotate-ref' (which only records where reference TEs are, to filter "
            "novel calls), this reports them as insertion calls.",
        )
    )
    _menu_index_genome(
        subparsers.add_parser(
            "index-genome",
            formatter_class=CustomHelpFormatter,
            help="Index/format the reference genome (step 1)",
            description="Create samtools (.fai) and minimap2 (.mmi) indexes for the reference genome.",
        )
    )
    _menu_align_genome(
        subparsers.add_parser(
            "align-genome",
            formatter_class=CustomHelpFormatter,
            help="Align trimmed flanking reads to the reference genome (step 4)",
            description="Align trimmed flanking-read FASTQs to the reference genome with minimap2, "
            "producing a sorted, indexed BAM for insertion finding.",
        )
    )
    _menu_find_insertions(
        subparsers.add_parser(
            "find-insertions",
            formatter_class=CustomHelpFormatter,
            help="Find non-reference TE insertions from genome-aligned flanking reads (step 5)",
            description="Cluster genome-aligned flanking reads into candidate non-reference insertion "
            "sites and write the all_nonref_insert table consumed by 'characterize'.",
        )
    )
    _menu_characterize(
        subparsers.add_parser(
            "characterize",
            formatter_class=CustomHelpFormatter,
            help="Characterize insertion sites (homozygous/heterozygous/somatic/excision)",
            description="Characterize RelocaTE non-reference insertion sites using read support from "
            "genome-aligned BAM files, classifying each site by zygosity and excision status.",
        )
    )

    return parser


def main(args: list[str] | None = None) -> int:
    """Tool for identifying Transposable transposition from WGS data by comparison to a reference genome."""
    parser = build_parser()

    try:
        cli_args = args or sys.argv[1:]
        if not cli_args:
            parser.print_help(sys.stderr)
            raise SystemExit(0)

        parsed = parser.parse_args(cli_args)

        if getattr(parsed, "verbose", False):
            logger.setLevel("DEBUG")
            for handler in logger.handlers:
                handler.setLevel("DEBUG")
            logger.debug("Debug mode enabled.")

        return parsed.func(parsed)

    except KeyboardInterrupt:
        logger.warning("Terminated by user.")
        return 1

    except SystemExit as err:
        if err.code != 0:
            logger.error(err)
            return 1

    except Exception as err:
        logger.error(err)
        return 1

    return 0
