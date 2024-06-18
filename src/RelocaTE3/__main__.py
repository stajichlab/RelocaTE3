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
        for line in text:
            formatted_text.extend(textwrap.wrap(line, width))
        # The textwrap module is used only for formatting help.
        # Delay its import for speeding up the common usage of argparse.
        return formatted_text

    def _get_help_string(self, action):
        help = action.help
        pattern = r"\(default: .+\)"
        if (re.search(pattern, action.help) is None and action.default not in [argparse.SUPPRESS, None, False]):
            defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
            if action.option_strings or action.nargs in defaulting_nargs:
                help += " (default: %(default)s)"
        return help


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
    parser.add_argument("-i", "--input", type=argparse.FileType("r"),
                        default=sys.stdin,
                        help="Input file (or by stdin) (default: stdin)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose mode for debug")
    parser.add_argument("-V", "--version", action="version",
                        version=__version__)

    return parser


def main(args: list[str] | None = None) -> int:
    """Tool for identifying Transposable transposition from WGS data by comparison to a reference genome."""
    return args_parser(_menu, args, prog=__entry_points__[__name__], description=main.__doc__,
                       epilog=f"Written by {__author__}")


if __name__ == "__main__":
    main()
