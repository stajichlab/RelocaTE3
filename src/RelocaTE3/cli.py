"""RelocaTE3 command-line interface.

Subcommands are thin drivers over the library so the same code paths can be run
on a laptop, a single HPC node, or scattered by a workflow engine.
"""

from __future__ import annotations

import argparse
import re
import sys
import textwrap
from pathlib import Path

from RelocaTE3 import __author__, __version__, logger
from RelocaTE3.characterize import characterize_insertions, write_characterized
from RelocaTE3.genome_align import align_to_genome, read_read_repeat
from RelocaTE3.insertions import (
    find_insertions,
    read_insertions_gff,
    write_insertions_gff,
    write_insertions_txt,
)
from RelocaTE3.librelocate import RelocaTE
from RelocaTE3.pipeline import run_sample
from RelocaTE3.ReadLibrary import ReadLibrary
from RelocaTE3.reference_te import (
    find_reference_insertions,
    write_existing_te_bed_from_rm,
)


class CustomHelpFormatter(argparse.HelpFormatter):
    """HelpFormatter with text filling, line splitting, and default display."""

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
        help_text = action.help
        pattern = r"\(default: .+\)"
        if re.search(pattern, action.help or "") is None and action.default not in [
            argparse.SUPPRESS,
            None,
            False,
        ]:
            defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
            if action.option_strings or action.nargs in defaulting_nargs:
                help_text += " (default: %(default)s)"
        return help_text


def _add_trim_parser(subparsers: argparse._SubParsersAction) -> None:
    """Register the ``trim`` subcommand (Step 3: TE read identification)."""
    p = subparsers.add_parser(
        "trim",
        help="Map reads to a TE library and trim the TE-matching portion",
        formatter_class=CustomHelpFormatter,
    )
    p.add_argument(
        "-l",
        "--left",
        "--r1",
        dest="left",
        required=True,
        help="Left/R1 read file (FASTQ[.gz])",
    )
    p.add_argument(
        "-r",
        "--right",
        "--r2",
        dest="right",
        help="Right/R2 read file (FASTQ[.gz]) for paired-end",
    )
    p.add_argument(
        "-t", "--te", dest="te_library", required=True, help="TE/repeat consensus FASTA"
    )
    p.add_argument("-o", "--outdir", required=True, help="Output directory")
    p.add_argument(
        "--sample",
        default="sample",
        help="Sample/strain name for output and read groups",
    )
    p.add_argument("-c", "--threads", type=int, default=1, help="Number of CPU threads")
    p.add_argument(
        "--len-cut-match", type=int, default=10, help="Minimum read-to-TE match length"
    )
    p.add_argument(
        "--len-cut-trim",
        type=int,
        default=10,
        help="Minimum trimmed flank length to keep",
    )
    p.add_argument(
        "--mismatch",
        type=int,
        default=0,
        help="Mismatches allowed in read-to-TE alignment",
    )
    p.set_defaults(func=_run_trim)


def _run_trim(args: argparse.Namespace) -> int:
    """Execute the ``trim`` subcommand."""
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fileset = [args.left] + ([args.right] if args.right else [])
    reads = ReadLibrary(fileset, args.sample)
    relocate = RelocaTE(
        TElib=args.te_library, threads=args.threads, verbose=1 if args.verbose else 0
    )
    relocate.identify_TE_reads(
        reads,
        outdir,
        TE_library=args.te_library,
        len_cut_match=args.len_cut_match,
        len_cut_trim=args.len_cut_trim,
        mismatch_allowance=args.mismatch,
    )
    return 0


def _add_align_genome_parser(subparsers: argparse._SubParsersAction) -> None:
    """Register the ``align-genome`` subcommand (Step 4: re-align to genome)."""
    p = subparsers.add_parser(
        "align-genome",
        help="Re-align trimmed flanking reads and supporting mates to the genome",
        formatter_class=CustomHelpFormatter,
    )
    p.add_argument(
        "-l",
        "--left",
        "--r1",
        dest="left",
        required=True,
        help="Original left/R1 read file",
    )
    p.add_argument(
        "-r",
        "--right",
        "--r2",
        dest="right",
        help="Original right/R2 read file (paired-end)",
    )
    p.add_argument("-g", "--genome", required=True, help="Reference genome FASTA")
    p.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory (must contain trim-step output)",
    )
    p.add_argument("--sample", default="sample", help="Sample/strain name")
    p.add_argument("-c", "--threads", type=int, default=1, help="Number of CPU threads")
    p.set_defaults(func=_run_align_genome)


def _run_align_genome(args: argparse.Namespace) -> int:
    """Execute the ``align-genome`` subcommand."""
    fileset = [args.left] + ([args.right] if args.right else [])
    reads = ReadLibrary(fileset, args.sample)
    bam, fullreads = align_to_genome(
        reads, args.genome, args.outdir, threads=args.threads
    )
    logger.info("Genome alignment written to %s (fullreads: %s)", bam, fullreads)
    return 0


def _add_find_insertions_parser(subparsers: argparse._SubParsersAction) -> None:
    """Register the ``find-insertions`` subcommand (Step 5: call insertions)."""
    p = subparsers.add_parser(
        "find-insertions",
        help="Cluster junction/supporting reads and call non-reference insertions",
        formatter_class=CustomHelpFormatter,
    )
    p.add_argument(
        "-b", "--bam", required=True, help="Genome BAM from the align-genome step"
    )
    p.add_argument("-g", "--genome", required=True, help="Reference genome FASTA")
    p.add_argument(
        "--fullreads-bam",
        help="Full-read genome BAM from align-genome (for false-junction filtering)",
    )
    p.add_argument(
        "--read-repeat",
        required=True,
        help="read_repeat_name.txt table written by the trim step",
    )
    p.add_argument("-o", "--outdir", required=True, help="Output directory")
    p.add_argument("--sample", default="sample", help="Sample/strain name")
    p.add_argument(
        "--required-junction-reads",
        type=int,
        default=1,
        help="Minimum junction reads on a side to call an insertion",
    )
    p.set_defaults(func=_run_find_insertions)


def _run_find_insertions(args: argparse.Namespace) -> int:
    """Execute the ``find-insertions`` subcommand."""
    read_repeat = read_read_repeat(Path(args.read_repeat))
    insertions = find_insertions(
        args.bam,
        read_repeat,
        args.genome,
        fullreads_bam=args.fullreads_bam,
        required_junction_reads=args.required_junction_reads,
    )
    results_dir = Path(args.outdir) / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    gff = results_dir / f"{args.sample}.all_nonref_insert.gff"
    write_insertions_gff(insertions, gff, args.sample)
    write_insertions_txt(
        insertions, results_dir / f"{args.sample}.all_nonref_insert.txt"
    )
    logger.info("Wrote %d insertions to %s", len(insertions), gff)
    return 0


def _add_run_parser(subparsers: argparse._SubParsersAction) -> None:
    """Register the ``run`` subcommand (full single-sample pipeline)."""
    p = subparsers.add_parser(
        "run",
        help="Run the full single-sample pipeline (trim + genome re-alignment)",
        formatter_class=CustomHelpFormatter,
    )
    p.add_argument(
        "-l",
        "--left",
        "--r1",
        dest="left",
        required=True,
        help="Left/R1 read file (FASTQ[.gz])",
    )
    p.add_argument(
        "-r",
        "--right",
        "--r2",
        dest="right",
        help="Right/R2 read file (FASTQ[.gz]) for paired-end",
    )
    p.add_argument(
        "-t", "--te", dest="te_library", required=True, help="TE/repeat consensus FASTA"
    )
    p.add_argument("-g", "--genome", required=True, help="Reference genome FASTA")
    p.add_argument("-o", "--outdir", required=True, help="Output directory")
    p.add_argument(
        "--repeatmasker",
        help="RepeatMasker .out for the genome (enables reference/shared insertions + FP filter)",
    )
    p.add_argument(
        "--genotype",
        action="store_true",
        help="Align original reads to the genome and classify zygosity (Step 7)",
    )
    p.add_argument("--sample", default="sample", help="Sample/strain name")
    p.add_argument("-c", "--threads", type=int, default=1, help="Number of CPU threads")
    p.add_argument(
        "--len-cut-match", type=int, default=10, help="Minimum read-to-TE match length"
    )
    p.add_argument(
        "--len-cut-trim",
        type=int,
        default=10,
        help="Minimum trimmed flank length to keep",
    )
    p.add_argument(
        "--mismatch",
        type=int,
        default=0,
        help="Mismatches allowed in read-to-TE alignment",
    )
    p.set_defaults(func=_run_pipeline)


def _run_pipeline(args: argparse.Namespace) -> int:
    """Execute the ``run`` subcommand."""
    fileset = [args.left] + ([args.right] if args.right else [])
    reads = ReadLibrary(fileset, args.sample)
    gff = run_sample(
        reads,
        args.te_library,
        args.genome,
        args.outdir,
        repeatmasker=args.repeatmasker,
        genotype=args.genotype,
        threads=args.threads,
        len_cut_match=args.len_cut_match,
        len_cut_trim=args.len_cut_trim,
        mismatch_allowance=args.mismatch,
        verbose=1 if args.verbose else 0,
    )
    logger.info("Pipeline complete; non-reference insertions at %s", gff)
    return 0


def _add_characterize_parser(subparsers: argparse._SubParsersAction) -> None:
    """Register the ``characterize`` subcommand (Step 7: genotyping)."""
    p = subparsers.add_parser(
        "characterize",
        help="Genotype insertions (homo/heterozygous/somatic) from a reads-to-genome BAM",
        formatter_class=CustomHelpFormatter,
    )
    p.add_argument(
        "-i", "--insertions", required=True, help="Non-reference insertion GFF"
    )
    p.add_argument(
        "-b",
        "--reads-bam",
        required=True,
        help="Original reads aligned to the genome (sorted+indexed BAM)",
    )
    p.add_argument("-o", "--outdir", required=True, help="Output directory")
    p.add_argument("--sample", default="sample", help="Sample/strain name")
    p.set_defaults(func=_run_characterize)


def _run_characterize(args: argparse.Namespace) -> int:
    """Execute the ``characterize`` subcommand."""
    insertions = read_insertions_gff(args.insertions)
    characterize_insertions(insertions, args.reads_bam)
    results_dir = Path(args.outdir) / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    write_characterized(
        insertions,
        results_dir / f"{args.sample}.all_nonref_insert.characTErized.gff",
        results_dir / f"{args.sample}.all_nonref_insert.characTErized.txt",
        args.sample,
    )
    logger.info("Genotyping complete for %d insertions", len(insertions))
    return 0


def _add_find_reference_parser(subparsers: argparse._SubParsersAction) -> None:
    """Register the ``find-reference`` subcommand (Steps 0/6: reference TEs)."""
    p = subparsers.add_parser(
        "find-reference",
        help="Build existingTE.bed and call reference/shared insertions",
        formatter_class=CustomHelpFormatter,
    )
    p.add_argument(
        "-b", "--bam", required=True, help="Genome BAM from the align-genome step"
    )
    p.add_argument(
        "--repeatmasker", required=True, help="RepeatMasker .out for the genome"
    )
    p.add_argument(
        "--read-repeat",
        required=True,
        help="read_repeat_name.txt table written by the trim step",
    )
    p.add_argument("-o", "--outdir", required=True, help="Output directory")
    p.add_argument("--sample", default="sample", help="Sample/strain name")
    p.set_defaults(func=_run_find_reference)


def _run_find_reference(args: argparse.Namespace) -> int:
    """Execute the ``find-reference`` subcommand."""
    outdir = Path(args.outdir)
    results_dir = outdir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    reference_tes = write_existing_te_bed_from_rm(
        args.repeatmasker, outdir / "existingTE.bed"
    )
    read_repeat = read_read_repeat(Path(args.read_repeat))
    ref_insertions = find_reference_insertions(args.bam, read_repeat, reference_tes)
    gff = results_dir / f"{args.sample}.all_ref_insert.gff"
    write_insertions_gff(ref_insertions, gff, args.sample)
    write_insertions_txt(
        ref_insertions, results_dir / f"{args.sample}.all_ref_insert.txt"
    )
    logger.info("Wrote %d reference/shared insertions to %s", len(ref_insertions), gff)
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Construct the top-level argument parser."""
    parser = argparse.ArgumentParser(
        prog="relocaTE3",
        formatter_class=CustomHelpFormatter,
        description="Identify transposable element transposition from resequencing data.",
        epilog=f"Written by {__author__}",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose/debug logging"
    )
    parser.add_argument("-V", "--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", metavar="<command>")
    _add_trim_parser(subparsers)
    _add_align_genome_parser(subparsers)
    _add_find_insertions_parser(subparsers)
    _add_find_reference_parser(subparsers)
    _add_characterize_parser(subparsers)
    _add_run_parser(subparsers)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point for the ``relocaTE3`` console script."""
    parser = build_parser()
    argv = argv if argv is not None else sys.argv[1:]
    args = parser.parse_args(argv)

    if getattr(args, "verbose", False):
        logger.setLevel("DEBUG")
        for handler in logger.handlers:
            handler.setLevel("DEBUG")
        logger.debug("Debug mode enabled.")

    if not getattr(args, "command", None):
        parser.print_help(sys.stderr)
        return 0

    try:
        return args.func(args)
    except KeyboardInterrupt:
        logger.warning("Terminated by user.")
        return 1
    except Exception as err:  # noqa: BLE001 - top-level CLI guard
        logger.error(err)
        return 1
