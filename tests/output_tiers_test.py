"""Tiered non-reference output, matching RelocaTE2's clean_false_positive.py.

RelocaTE2 does not publish one call set, it publishes four, produced by
``clean_false_positive.py`` from the raw calls:

    ALL.all_nonref_insert.raw.gff        everything
    ALL.all_nonref_insert.all.gff        minus calls near an existing reference TE
    ALL.all_nonref_insert.gff            minus singleton/insufficient_data/
                                         supporting_reads  <- RelocaTE2's headline
    ALL.all_nonref_insert.high_conf.gff  minus one-sided junction calls

RelocaTE3 emitted only the equivalent of ``.all``, so users comparing against a
RelocaTE2 run were comparing an unfiltered set to a filtered one. On the Chr3
2 Mb fixture that is the whole precision difference: every one of RelocaTE3's
21 false positives is a one-sided call, and the high_conf tier scores 178/178
(precision 1.000).

The exact filters are transcribed from clean_false_positive.py:107-108.

Note which file gets filtered. RelocaTE2 runs clean_false_positive.py on the
**GFF only** (relocaTE2.py:704); the .txt is concatenated unfiltered
(relocaTE2.py:703) and characterizer.pl genotypes that unfiltered table
(relocaTE2.py:707). So the plain .txt must keep every call -- filtering it would
silently change what `characterize` sees.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from RelocaTE3.insertions import write_insertion_tiers

# TE, TSD-or-class, sample, chrom, start..end, strand, T:, R:, L:, ST:, SR:, SL:
TWO_SIDED = "mping\tTTA\tHEG4\tChr3\t100..102\t+\tT:6\tR:3\tL:3\tST:0\tSR:9\tSL:8"
ONE_SIDED_R = "mping\tUNK\tHEG4\tChr3\t200..200\t-\tT:2\tR:2\tL:0\tST:0\tSR:5\tSL:0"
ONE_SIDED_L = "mping\tUNK\tHEG4\tChr3\t300..300\t+\tT:1\tR:0\tL:1\tST:0\tSR:0\tSL:4"
SUPPORT_ONLY = (
    "NA\tsupporting_reads\tHEG4\tChr3\t400..500\t+\tT:0\tR:0\tL:0\tST:9\tSR:5\tSL:4"
)
SINGLETON = "mping\tsingleton\tHEG4\tChr3\t600..600\t+\tT:1\tR:1\tL:0\tST:0\tSR:1\tSL:0"
INSUFFICIENT = (
    "mping\tinsufficient_data\tHEG4\tChr3\t700..700\t+\tT:1\tR:0\tL:1\tST:0\tSR:0\tSL:1"
)

ALL_ROWS = [
    TWO_SIDED,
    ONE_SIDED_R,
    ONE_SIDED_L,
    SUPPORT_ONLY,
    SINGLETON,
    INSUFFICIENT,
]


@pytest.fixture
def written(tmp_path: Path) -> Path:
    table = tmp_path / "ALL.mping.all_nonref_insert.txt"
    table.write_text("\n".join(ALL_ROWS) + "\n")
    write_insertion_tiers(table, sample="HEG4")
    return table


def _tier(table: Path, suffix: str) -> Path:
    """Sibling tier path. The filename is dotted, so with_suffix() is unsafe."""
    stem = str(table)[: -len(".txt")]
    return Path(f"{stem}{suffix}")


def _rows(path: Path) -> list[str]:
    return [ln for ln in path.read_text().splitlines() if ln.strip()]


def _starts(path: Path) -> set[str]:
    return {ln.split("\t")[4] for ln in _rows(path)}


def _starts_gff(path: Path) -> set[str]:
    return {ln.split("\t")[3] for ln in _rows(path)}


def test_all_tier_keeps_every_call(written: Path):
    """`.all` is the unfiltered set, matching the plain table."""
    assert len(_rows(_tier(written, ".all.gff"))) == len(ALL_ROWS)
    assert len(_rows(_tier(written, ".all.txt"))) == len(ALL_ROWS)


def test_headline_gff_drops_the_three_low_confidence_classes(written: Path):
    """clean_false_positive.py:107 greps out singleton/insufficient_data/supporting_reads."""
    assert _starts_gff(_tier(written, ".gff")) == {"100", "200", "300"}


def test_plain_table_is_left_unfiltered(written: Path):
    """RelocaTE2 filters the GFF, not the table; characterize reads the table.

    relocaTE2.py:703 concatenates the .txt without filtering and :707 feeds it
    straight to characterizer.pl. Filtering it here would drop calls from
    genotyping that RelocaTE2 genotypes.
    """
    assert len(_rows(written)) == len(ALL_ROWS)


def test_high_conf_additionally_drops_one_sided_calls(written: Path):
    """clean_false_positive.py:108 removes Right=1;Left=0 and Right=0;Left=1."""
    assert _starts_gff(_tier(written, ".high_conf.gff")) == {"100"}


def test_every_tier_has_a_matching_gff(written: Path):
    """RelocaTE2's headline output is a GFF; find-insertions wrote none at all."""
    for suffix in (".txt", ".gff", ".all.txt", ".all.gff",
                   ".high_conf.txt", ".high_conf.gff"):
        assert _tier(written, suffix).is_file(), f"missing {suffix}"


def test_gff_carries_the_relocate2_attributes(written: Path):
    """Downstream RelocaTE2 tooling greps these attribute keys."""
    gff = _tier(written, ".gff")
    first = _rows(gff)[0]
    for key in ("TSD=", "Name=", "Right_junction_reads=", "Left_junction_reads=",
                "Right_support_reads=", "Left_support_reads="):
        assert key in first, f"{key} missing from GFF attributes"


def test_table_keeps_its_name_and_contents(tmp_path: Path):
    """characterize/run-all consume the unsuffixed .txt -- it must not move or shrink."""
    table = tmp_path / "ALL.mping.all_nonref_insert.txt"
    table.write_text("\n".join(ALL_ROWS) + "\n")
    before = table.read_text()
    write_insertion_tiers(table, sample="HEG4")
    assert table.read_text() == before


def test_tiers_are_nested(written: Path):
    """high_conf subset of headline subset of all -- a tier can only remove."""
    high = _starts_gff(_tier(written, ".high_conf.gff"))
    head = _starts_gff(_tier(written, ".gff"))
    every = _starts_gff(_tier(written, ".all.gff"))
    assert high <= head <= every


def test_empty_input_still_writes_every_tier(tmp_path: Path):
    """A sample with no calls must not leave downstream steps missing files."""
    table = tmp_path / "ALL.mping.all_nonref_insert.txt"
    table.write_text("")
    write_insertion_tiers(table, sample="HEG4")
    for suffix in (".gff", ".all.txt", ".all.gff", ".high_conf.txt", ".high_conf.gff"):
        assert _tier(table, suffix).is_file()
        assert _rows(_tier(table, suffix)) == []
