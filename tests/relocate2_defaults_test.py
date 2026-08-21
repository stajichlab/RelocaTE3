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


def test_both_junctions_are_required_by_default():
    """Deliberate divergence from RelocaTE2, justified by measurement.

    RelocaTE2 keeps a site on `l_count >= 1 OR r_count >= 1`
    (relocaTE_insertionFinder.py:365) and its one-sided calls are reliable --
    on mPing all 53 of them are true. RelocaTE3's are not: 86 one-sided calls,
    only 33 true (38%). Matching RelocaTE2's *policy* without matching its
    candidate quality is far worse than diverging from it:

        riceTElib F1   5x     15x    30x
          one-sided allowed   0.170  0.155  0.144
          both required       0.448  0.643  0.722
          RelocaTE2           0.439  0.622  0.701

    Allowing one-sided calls degrades *with* coverage on a multi-family library
    (more depth, more noise), so the failure is worst exactly where users have
    the best data. Requiring both junctions is what makes RelocaTE3 beat
    RelocaTE2 on whole-library runs; opt out with --no-require-both-junctions
    when sensitivity matters more than precision (single-element studies).
    """
    assert _defaults(REQUIRED_RUN_ALL)["require_both_junctions"] is True
    assert _defaults(REQUIRED_FIND)["require_both_junctions"] is True


def test_one_sided_calls_can_be_re_enabled():
    """The RelocaTE2-equivalent policy stays reachable, for single-element work."""
    ns = _defaults(REQUIRED_RUN_ALL + ["--no-require-both-junctions"])
    assert ns["require_both_junctions"] is False
    ns = _defaults(REQUIRED_FIND + ["--no-require-both-junctions"])
    assert ns["require_both_junctions"] is False


def test_explicit_flag_still_works():
    """Existing scripts passing --require-both-junctions must not break."""
    assert _defaults(REQUIRED_FIND + ["--require-both-junctions"])["require_both_junctions"] is True


def test_min_mapq_defaults_to_no_gate():
    """RelocaTE2 has no MAPQ admission gate, so neither may RelocaTE3's default.

    ``relocaTE_insertionFinder.py:1521-1558`` admits reads on XM/XO/XT/X1 only;
    MAPQ is used at :1523,1539 purely to mark a read low quality, and calls
    resting solely on low-quality reads are dropped later (:226-241). RelocaTE3
    defaulted to 1, silently discarding MAPQ-0 reads. On the mPing benchmark
    those are frequently the single read carrying the second junction, so the
    call collapsed to one-sided and was dropped: restoring the RelocaTE2 default
    recovered 4 true calls at cov30x_rep1 with no new false positives.
    """
    assert _defaults(REQUIRED_RUN_ALL)["min_mapq"] == 0
    assert _defaults(REQUIRED_FIND)["min_mapq"] == 0


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
