"""Per-sample counts for the MultiQC summary report.

``relocaTE3 sample-stats`` reads files a sample's run already wrote, in the
``run-all`` layout of one sample directory, and writes one MultiQC
custom-content row (``<id>.stats_mqc.tsv``). It only counts lines; it never
re-derives or filters calls.

==========================  ===================================================  ==============
column                      file (under the sample directory)                    count
==========================  ===================================================  ==============
``te_containing_reads``     ``te_containing/<id>.read_repeat_name.txt``          lines
``flanking_reads``          ``flanking/<id>.*.flankingReads.fq`` (all of them)   lines / 4
``nonref_raw``              ``results/<target>.<te>.all_nonref_insert.raw.txt``  lines
``nonref_all``              ``results/<target>.<te>.all_nonref_insert.all.txt``  lines
``nonref_headline``         ``results/<target>.<te>.all_nonref_insert.gff``      lines
``nonref_high_conf``        ``results/<target>.<te>.all_nonref_insert.high_conf.txt``  lines
``characterized``           ``results/<target>.<te>.all_nonref_insert.characTErized.txt``  lines - header
``reference``               ``results/<id>.all_ref_insert.txt``                  lines - header
==========================  ===================================================  ==============

Notes:

- The non-reference tier tables and GFFs have no header line; the
  characterized table and ``all_ref_insert.txt`` have exactly one.
- The headline tier (RelocaTE2's ``clean_false_positive.py`` set) exists only
  as ``<stem>.gff``: the plain ``<stem>.txt`` is the unfiltered table
  (``insertions.write_tiers``), so it would repeat ``nonref_raw``.
- ``read_repeat_name.txt`` has one line per TE-containing read record (a pair
  contributes one line per mate that hit the TE library).
- Blank lines and GFF ``#`` comment lines are not counted.
- ``characterized`` and ``reference`` are optional: their files exist only
  with genotyping and with a RepeatMasker annotation. A missing optional file
  leaves its column out (MultiQC shows an empty cell), never a zero.
  Any other missing file is an error.
"""

from __future__ import annotations

import os
from pathlib import Path

COLUMNS = (
    "te_containing_reads",
    "flanking_reads",
    "nonref_raw",
    "nonref_all",
    "nonref_headline",
    "nonref_high_conf",
    "characterized",
    "reference",
)
OPTIONAL = ("characterized", "reference")

# Column -> (title, description) for the MultiQC table header.
_HEADERS = {
    "te_containing_reads": ("TE reads", "Reads that hit the TE library (lines of read_repeat_name.txt)"),
    "flanking_reads": ("Flanking reads", "Trimmed flanking reads aligned to the genome (FASTQ records, all mates)"),
    "nonref_raw": ("Non-ref raw", "Non-reference calls, raw tier (every call)"),
    "nonref_all": ("Non-ref all", "Non-reference calls, all tier (minus one-sided calls at reference TE edges)"),
    "nonref_headline": ("Non-ref headline", "Non-reference calls, headline tier (RelocaTE2's set: minus singleton, insufficient_data and supporting_reads)"),
    "nonref_high_conf": ("Non-ref high-conf", "Non-reference calls, high_conf tier (minus single-read one-sided calls)"),
    "characterized": ("Characterized", "Non-reference calls with a zygosity call (characTErized table)"),
    "reference": ("Reference", "Reference (shared) TE calls (all_ref_insert.txt)"),
}  # fmt: skip

# MultiQC custom content: YAML in '# ' comment lines ahead of the TSV
# (https://github.com/MultiQC/MultiQC/blob/main/docs/markdown/custom_content/index.md,
# "MultiQC-specific data file").
# Every sample writes its own file with the same id; MultiQC merges them into
# one table.
_MQC_HEADER = [
    "# id: 'relocate3_sample_stats'",
    "# section_name: 'RelocaTE3 per-sample counts'",
    "# description: 'Reads and calls per sample, counted from the files each sample step wrote. Tiers: raw >= all >= headline >= high_conf. Empty cells: no genotyping (characterized) or no RepeatMasker annotation (reference).'",
    "# format: 'tsv'",
    "# plot_type: 'table'",
    "# pconfig:",
    "#     id: 'relocate3_sample_stats_table'",
    "#     title: 'RelocaTE3: reads and calls per sample'",
    "# headers:",
]

def _yaml_quote(text: str) -> str:
    """A YAML single-quoted scalar (a ' inside is written twice)."""
    return "'" + text.replace("'", "''") + "'"


for _col in COLUMNS:
    _title, _desc = _HEADERS[_col]
    _MQC_HEADER += [
        f"#     {_col}:",
        f"#         title: {_yaml_quote(_title)}",
        f"#         description: {_yaml_quote(_desc)}",
        "#         format: '{:,.0f}'",
    ]


class SampleStatsError(ValueError):
    """A file the counts need is missing or malformed."""


def _count_lines(path: Path, skip_comments: bool = False) -> int:
    n = 0
    with open(path, "rb") as fh:
        for line in fh:
            if not line.strip():
                continue
            if skip_comments and line.startswith(b"#"):
                continue
            n += 1
    return n


def _require(path: Path) -> Path:
    if not path.is_file():
        raise SampleStatsError(f"sample-stats: required file not found: {path}")
    return path


def _fastq_reads(path: Path) -> int:
    n = 0
    with open(path, "rb") as fh:
        for n, _ in enumerate(fh, start=1):
            pass
    if n % 4:
        raise SampleStatsError(
            f"sample-stats: {path} has {n} lines, not a multiple of 4 (truncated FASTQ?)"
        )
    return n // 4


def sample_stats(
    sample_dir: str | Path, name: str, target: str = "ALL", te_name: str = "repeat"
) -> dict[str, int]:
    """Counts for one sample directory (run-all layout), in :data:`COLUMNS` order."""
    d = Path(sample_dir)
    stem = d / "results" / f"{target}.{te_name}.all_nonref_insert"

    flanking = sorted((d / "flanking").glob(f"{name}.*.flankingReads.fq"))
    if not flanking:
        raise SampleStatsError(
            f"sample-stats: no {name}.*.flankingReads.fq in {d / 'flanking'}"
        )

    stats: dict[str, int] = {
        "te_containing_reads": _count_lines(
            _require(d / "te_containing" / f"{name}.read_repeat_name.txt")
        ),
        "flanking_reads": sum(_fastq_reads(f) for f in flanking),
        "nonref_raw": _count_lines(_require(Path(f"{stem}.raw.txt"))),
        "nonref_all": _count_lines(_require(Path(f"{stem}.all.txt"))),
        "nonref_headline": _count_lines(_require(Path(f"{stem}.gff")), skip_comments=True),
        "nonref_high_conf": _count_lines(_require(Path(f"{stem}.high_conf.txt"))),
    }
    optional = {
        "characterized": Path(f"{stem}.characTErized.txt"),
        "reference": d / "results" / f"{name}.all_ref_insert.txt",
    }
    for col, path in optional.items():
        if path.is_file():
            # One header line, then one row per call.
            stats[col] = max(_count_lines(path) - 1, 0)
    return stats


def stats_lines(name: str, stats: dict[str, int]) -> list[str]:
    columns = [c for c in COLUMNS if c in stats]
    return [
        *_MQC_HEADER,
        "\t".join(["Sample", *columns]),
        "\t".join([name, *(str(stats[c]) for c in columns)]),
    ]


def write_sample_stats(
    sample_dir: str | Path,
    name: str,
    output: str | Path,
    target: str = "ALL",
    te_name: str = "repeat",
) -> Path:
    """Count, then write ``output`` atomically (temp file + rename)."""
    lines = stats_lines(name, sample_stats(sample_dir, name, target, te_name))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_name(f".{output.name}.tmp")
    try:
        tmp.write_text("".join(line + "\n" for line in lines))
        os.replace(tmp, output)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    return output
