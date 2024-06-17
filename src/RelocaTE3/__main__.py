"""RelocaTE3 Interface and CLI logic."""
from __future__ import annotations

import argparse
import logging
import re
import sys
import textwrap
from typing import Callable

from RelocaTE3 import __author__, __entry_points__, __version__


class CustomHelpFormatter(argparse.HelpFormatter):
    """HelpFormatter that have customized function for text filling, line splitting and default parameter showing."""

    def _fill_text(self, text, width, indent):
        text = [self._whitespace_matcher.sub(" ", line).strip() for line in text.split("\n\n") if line != ""]
        return "\n\n".join([textwrap.fill(line, width) for line in text])

    def _split_lines(self, text, width):
        text = [self._whitespace_matcher.sub(" ", line).strip() for line in text.split("\n") if line != ""]
        formatted_text = []
        [formatted_text.extend(textwrap.wrap(line, width)) for line in text]
        # The textwrap module is used only for formatting help.
        # Delay its import for speeding up the common usage of argparse.
        return formatted_text

    def _get_help_string(self, action):
        helpstr = action.help
        pattern = r"\(default: .+\)"
        if (re.search(pattern, action.help) is None and action.default not in [argparse.SUPPRESS, None, False]):
            defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
            if action.option_strings or action.nargs in defaulting_nargs:
                helpstr += " (default: %(default)s)"
        return helpstr


def args_parser(
    parser_func: Callable[[argparse.ArgumentParser], argparse.ArgumentParser],
    args: list[str] | None,
    *,
    prog: str | None = None,
    description: str | None = None,
    epilog: str | None = None,
):
    """Preset menu structure for entry-point scripts."""
    try:
        logging.basicConfig(format="%(asctime)s %(name)s %(levelname)s %(message)s", level="INFO")
        logger = logging.getLogger(prog)

        parser = argparse.ArgumentParser(
            prog=prog,
            formatter_class=CustomHelpFormatter,
            description=description,
            epilog=epilog,
        )
        parser = parser_func(parser)

        args = args or sys.argv[1:]
        if not args or "help" in args:
            parser.print_help(sys.stderr)
            raise SystemExit(0)
        args = parser.parse_args(args)

        if args.verbose:
            logger.setLevel("DEBUG")
            for handler in logger.handlers:
                handler.setLevel("DEBUG")
            logging.debug("Debug mode enabled.")

        args.func(**vars(args))

    except KeyboardInterrupt:
        logging.warning("Terminated by user.")
        return 1

    except SystemExit as err:
        if err.code != 0:
            logging.error(err)
            return 1

    return 0


def _menu(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Menu for this entry point."""
    # do we want to have multiple libraries processed for each strain in a run?
    #
"""
  -b BAM, --bam BAM     Name of BAM file of reads mapped reference genome
  -t TE_FASTA, --te_fasta TE_FASTA
                        Name of fasta sequence of repeat element
  -d FQ_DIR, --fq_dir FQ_DIR
                        Name of directory of input fastq sequence data
  -g GENOME_FASTA, --genome_fasta GENOME_FASTA
                        Name of fasta file of reference genome sequence
  -r REFERENCE_INS, --reference_ins REFERENCE_INS
                        Name of RepeatMasker TE annotation of reference genome
  -o OUTDIR, --outdir OUTDIR
                        Name of output directory where to put temperary and
                        final results
  -s SIZE, --size SIZE  Insert size of sequence library, default = 500
  -c CPU, --cpu CPU     Number of CPUs to use for multiplex, default = 1
  -1 MATE_1_ID, --mate_1_id MATE_1_ID
                        string define paired-end read1, default = "_1"
  -2 MATE_2_ID, --mate_2_id MATE_2_ID
                        string define paired-end read2, default = "_2"
  -u UNPAIRED_ID, --unpaired_id UNPAIRED_ID
                        string defining single-end reads, default = ".unPaired"
  --sample SAMPLE       string defining sample name which will present in output
                        GFF, default = "not_given"
  --aligner ALIGNER     aligner used to map reads to repeat elements,
                        default=blat
  --len_cut_match LEN_CUT_MATCH
                        length cutoff threshold for match between reads and
                        repeat elements. Large value will lead to less
                        sensitive but more accuracy, default = 10
  --len_cut_trim LEN_CUT_TRIM
                        length cutoff threshold for trimed reads after
                        trimming repeat sequence from reads. Large value will
                        lead to less sensitive but more accuracy, default = 10
  --mismatch MISMATCH   Number of mismatches allowed for matches between reads
                        and repeat elements, default = 2
  --mismatch_junction MISMATCH_JUNCTION
                        Number of mismatches allowed for matches between
                        junction reads and repeat elements, default = 2
  --step STEP           Number to control steps of pipeline, default =
                        "1234567"
  --dry_run             write shell scripts only while this script excute
  --run                 run while this script excute
  --split               split fastq into 1 million reads chunks to run blat/bwa jobs
  -v VERBOSE, --verbose VERBOSE
                        verbose grade to print out information in all scripts:
                        range from 0 to 4, default = 2 """
    parser.add_argument("-i", "--input", type=argparse.FileType("r"),
                        default=sys.stdin,
                        help="Input file (or by stdin) (default: stdin)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose mode for debug")
    parser.add_argument("-V", "--version", action="version",
                        version=__version__)

    return parser


def main(args: list[str] | None = None) -> int:
    """Tool for identifying Transposable transposition from resequencing data by comparison to a reference genome"""
    return args_parser(_menu, args, prog=__entry_points__[__name__],
                        description=main.__doc__,
                        epilog=f"Written by {__author__}")


if __name__ == "__main__":
    main()
