"""Every shipped default must match RelocaTE2's, so an untuned run reproduces it.

RelocaTE3 v1 aims to be a drop-in for RelocaTE2: a user who passes only inputs
should get RelocaTE2's behaviour. These values come from RelocaTE2's own source,
cited per assertion, not from the RelocaTE3 benchmark's tuned configuration.
"""

from __future__ import annotations

import inspect

import pytest

from RelocaTE3.cli import build_parser


def _defaults(argv: list[str]) -> dict:
    """Parse a subcommand with only its required args and return the namespace."""
    return vars(build_parser().parse_args(argv))


REQUIRED_RUN = ["run", "-l", "r1.fq", "-T", "te.fa", "-n", "S"]
REQUIRED_RUN_ALL = ["run-all", "-1", "r1.fq", "-T", "te.fa", "-g", "g.fa", "-n", "S"]
REQUIRED_TRIM = ["trim", "-b", "x.bam", "-n", "S"]
REQUIRED_FIND = ["find-insertions", "-b", "x.bam", "-R", "rr.txt", "-n", "S"]


@pytest.mark.parametrize(
    "argv,key",
    [
        (REQUIRED_RUN, "mismatch_allowance"),
        (REQUIRED_RUN_ALL, "mismatch_allowance"),
        (REQUIRED_TRIM, "mismatch_allowance"),
        (REQUIRED_FIND, "mismatch_allow"),
    ],
)
def test_mismatch_defaults_to_two(argv, key):
    """relocaTE2.py:207-208 -- --mismatch and --mismatch_junction both default 2."""
    assert _defaults(argv)[key] == 2


def test_te_aligner_defaults_to_blat():
    """relocaTE2.py:204 -- --aligner default='blat'."""
    assert _defaults(REQUIRED_RUN)["te_aligner"] == "blat"
    assert _defaults(REQUIRED_RUN_ALL)["te_aligner"] == "blat"


def test_genome_aligner_defaults_to_bwaaln():
    """RelocaTE2 places flanking reads with `bwa aln` (relocaTE_align.py)."""
    assert _defaults(REQUIRED_RUN_ALL)["genome_aligner"] == "bwaaln"


def test_length_cutoffs_default_to_ten():
    """relocaTE2.py:205-206 -- len_cut_match and len_cut_trim both default 10."""
    ns = _defaults(REQUIRED_RUN_ALL)
    assert ns["minimum_match_length"] == 10
    assert ns["minimum_trimmed_length"] == 10


def test_tsd_defaults_to_unk():
    """relocaTE2.py:346 writes regex.txt with a hardcoded UNK TSD field.

    RelocaTE2 never asks the user for a TSD, so neither should RelocaTE3 --
    find-insertions must default it rather than require it.
    """
    assert _defaults(REQUIRED_RUN_ALL)["tsd"] == "UNK"
    assert _defaults(REQUIRED_FIND)["tsd"] == "UNK"


def test_both_junctions_are_not_required():
    """relocaTE_insertionFinder.py:365 keeps a site on `l_count >= 1 OR r_count >= 1`.

    RelocaTE2 emits single-junction insertions, so requiring both would silently
    drop calls RelocaTE2 reports. --require-both-junctions must stay opt-in.
    """
    assert _defaults(REQUIRED_RUN_ALL)["require_both_junctions"] is False
    assert _defaults(REQUIRED_FIND)["require_both_junctions"] is False


def test_required_junction_reads_defaults_to_one():
    """relocaTE_insertionFinder.py:1732-1734 -- required_(left|right)_reads = 1."""
    from RelocaTE3.insertions import find_insertions

    sig = inspect.signature(find_insertions)
    assert sig.parameters["required_junction_reads"].default == 1


def test_library_entry_points_agree_with_the_cli():
    """pipeline.run_sample is a second front door; it must not drift from the CLI."""
    from RelocaTE3.pipeline import run_sample

    params = inspect.signature(run_sample).parameters
    assert params["mismatch_allowance"].default == 2
    assert params["te_aligner"].default == "blat"
    assert params["genome_aligner"].default == "bwaaln"
    assert params["required_junction_reads"].default == 1
