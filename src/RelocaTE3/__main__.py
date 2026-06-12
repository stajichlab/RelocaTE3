"""RelocaTE3 Interface and CLI logic."""

from __future__ import annotations

import argparse
import re
import sys
import textwrap
from pathlib import Path

from RelocaTE3 import __author__, __entry_points__, __version__, logger


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
        "--aligner",
        default="minimap2",
        choices=["minimap2", "bwa"],
        help="Alignment tool",
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
        "--aligner",
        default="minimap2",
        choices=["minimap2", "bwa"],
        help="Alignment tool",
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
        metavar="BAM",
        dest="bam",
        help="BAM file(s) of original reads aligned to the reference genome (before TE trimming)",
    )
    parser.add_argument(
        "-g",
        "--genome-fasta",
        dest="genome_fasta",
        metavar="FASTA",
        help="Reference genome FASTA (required with --excision)",
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
        "--paired",
        action="store_true",
        help="Align the first two FASTQs as a read pair",
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="CPU threads for minimap2"
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
        help="TSD motif (e.g. TTA); TSD-unknown mode is not yet supported",
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
    _add_common_args(parser)
    parser.set_defaults(func=cmd_find_insertions)
    return parser


# ---------------------------------------------------------------------------
# Command handlers
# ---------------------------------------------------------------------------


def cmd_map(
    left, right, te_library, name, outdir, threads, aligner, verbose, **kwargs
) -> int:
    """Align reads to TE library and write BAM files."""
    from RelocaTE3.align import Aligner
    from RelocaTE3.ReadLibrary import ReadLibrary

    fileset = [left] + ([right] if right else [])
    reads = ReadLibrary(fileset, name)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    aln = Aligner(threads=threads, default_aligner=aligner)
    bamfiles = aln.map_minimap_library(reads, out, te_library)
    logger.info("%d BAM file(s) written to %s", len(bamfiles), outdir)
    return 0


def cmd_trim(
    bam,
    name,
    minimum_match_length,
    minimum_trimmed_length,
    mismatch_allowance,
    outdir,
    verbose,
    **kwargs,
) -> int:
    """Trim TE sequence from TE-library BAM reads and emit junction-named flanking FASTQs."""
    from RelocaTE3.librelocate import RelocaTE

    bam_paths = [Path(b) for b in bam]
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    relocate = RelocaTE(verbose=int(verbose))
    directions = relocate._bam_directions(bam_paths)
    flank_written = relocate.write_trimmed_reads(
        name,
        list(zip(directions, bam_paths)),
        out,
        minimum_match_length=minimum_match_length,
        minimum_trimmed_length=minimum_trimmed_length,
        mismatch_allowance=mismatch_allowance,
    )
    logger.info("Wrote %d flanking read(s) to %s", flank_written, out / "flanking")
    return 0


def cmd_run(
    left,
    right,
    te_library,
    name,
    outdir,
    threads,
    aligner,
    minimum_match_length,
    minimum_trimmed_length,
    mismatch_allowance,
    verbose,
    **kwargs,
) -> int:
    """Run full pipeline: align reads to TE library then trim TE sequences."""
    from RelocaTE3.librelocate import RelocaTE
    from RelocaTE3.ReadLibrary import ReadLibrary

    fileset = [left] + ([right] if right else [])
    reads = ReadLibrary(fileset, name)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    relocate = RelocaTE(TElib=te_library, threads=threads, verbose=int(verbose))
    n = relocate.identify_TE_reads(reads, out, search_tool=aligner)
    logger.info("%d read(s) written", n)
    return 0


def cmd_characterize(
    sites_file,
    bam,
    genome_fasta,
    outdir,
    excision,
    samtools,
    bcftools,
    verbose,
    **kwargs,
) -> int:
    """Characterize insertion sites as homozygous, heterozygous, somatic or excised."""
    from RelocaTE3.characterize import Characterizer

    characterizer = Characterizer(
        samtools=samtools, bcftools=bcftools, verbose=int(verbose)
    )
    bam_paths = [Path(b) for b in bam]
    txt_path, gff_path = characterizer.characterize(
        sites_file=Path(sites_file),
        bam_files=bam_paths,
        genome_fasta=Path(genome_fasta) if genome_fasta else None,
        outdir=Path(outdir) if outdir else None,
        excision=excision,
    )
    logger.info("Characterization written to %s and %s", txt_path, gff_path)
    return 0


def cmd_annotate_ref(
    te_library,
    genome_fasta,
    outdir,
    threads,
    min_identity,
    min_coverage,
    verbose,
    **kwargs,
) -> int:
    """Annotate existing/reference TE copies in the genome with minimap2 (step 0)."""
    from RelocaTE3.reference_te import ReferenceTEAnnotator

    annotator = ReferenceTEAnnotator(threads=threads, verbose=int(verbose))
    bed = annotator.annotate_minimap(
        te_library=Path(te_library),
        genome_fasta=Path(genome_fasta),
        outdir=Path(outdir),
        min_identity=min_identity,
        min_coverage=min_coverage,
    )
    logger.info("Existing-TE annotation written to %s", bed)
    return 0


def cmd_index_genome(genome_fasta, force, verbose, **kwargs) -> int:
    """Index/format the reference genome (samtools faidx + minimap2 index) (step 1)."""
    from RelocaTE3.align import Aligner

    aln = Aligner()
    aln.verbose = bool(verbose)
    aln.index_genome(genome_fasta, force=force)
    logger.info("Indexed genome %s", genome_fasta)
    return 0


def cmd_align_genome(
    genome_fasta, fastq, name, outdir, paired, threads, verbose, **kwargs
) -> int:
    """Align trimmed flanking reads to the reference genome (step 4)."""
    from RelocaTE3.align import Aligner

    aln = Aligner(threads=threads)
    aln.verbose = bool(verbose)
    bam = aln.map_genome_minimap(
        genome=genome_fasta,
        fastqs=fastq,
        name=name,
        outdir=outdir,
        paired=paired,
    )
    logger.info("Genome-aligned BAM written to %s", bam)
    return 0


def cmd_find_insertions(
    bam,
    read_repeat,
    tsd,
    target,
    name,
    outdir,
    te_name,
    reference_ins,
    mismatch_allow,
    min_mapq,
    verbose,
    **kwargs,
) -> int:
    """Find non-reference insertions from genome-aligned flanking reads (step 5)."""
    from RelocaTE3.insertions import InsertionFinder

    finder = InsertionFinder(
        mismatch_allow=mismatch_allow, min_mapq=min_mapq, verbose=int(verbose)
    )
    out_txt = finder.find_insertions(
        bam_file=Path(bam),
        read_repeat_file=Path(read_repeat),
        tsd=tsd,
        target=target,
        sample=name,
        outdir=Path(outdir),
        te_name=te_name,
        reference_ins=Path(reference_ins) if reference_ins else None,
    )
    logger.info("Non-reference insertions written to %s", out_txt)
    return 0


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main(args: list[str] | None = None) -> int:
    """Tool for identifying Transposable transposition from WGS data by comparison to a reference genome."""
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
            help="Run full pipeline: map + trim",
            description="Run the complete RelocaTE3 pipeline: align reads to the TE library then trim TE sequences.",
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

        parsed.func(**vars(parsed))

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


if __name__ == "__main__":
    main()
