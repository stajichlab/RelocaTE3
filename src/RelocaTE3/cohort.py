"""Combine per-sample characterized tables into cohort-level outputs.

``relocaTE3 cohort`` reads the ``*.all_nonref_insert.characTErized.txt`` tables
that ``characterize`` writes (one per sample) and writes three files:

``cohort.characTErized.txt``
    The per-sample header once, then every input data row byte-for-byte,
    ordered by sample name and then by the row's position in its own file.

``cohort.genotype_matrix.tsv``
    A plain TSV (no comment lines, so any TSV reader works). Columns are
    ``chromosome start end strand TE`` followed by one column per sample
    (sorted). One row per *exact* site key ``(chromosome, start, end, strand,
    TE)`` parsed from ``chromosome.pos`` (``Chr1:1132975..1132977``), ``strand``
    and ``TE``. A cell is that sample's ``status`` at the site, or ``.`` when
    the sample has no call there. Nearby sites are NOT merged; deciding when
    two close sites are the same insertion is deferred to the VCF design.
    Rows sort by chromosome in natural order (see :func:`natural_key`), then
    start and end numerically, then strand and TE as text.

``cohort.summary_mqc.tsv``
    A MultiQC custom-content table (one row per sample). ``total_calls`` is the
    number of data rows. Each status is split on ``/``; the first part is the
    genotype class (``homozygous``, ``heterozygous``, ``somatic_insertion``,
    ``other``; any other first part gets its own column after ``excision``).
    ``excision`` counts statuses containing the text ``excision`` anywhere
    (``.../excision_no_footprint`` and ``.../excision_with_footprint``), so a
    ``homozygous/excision_no_footprint`` call counts once in ``homozygous``
    and once in ``excision``.

Inputs are given as plain paths, where the sample name is the single value of
the file's ``strain`` column, or as ``NAME=PATH``, where every ``strain``
value must equal ``NAME``. A header-only file (a sample with no characterized
calls) must use ``NAME=PATH``: it then appears in the matrix as all ``.`` and
in the summary with zero counts.

All inputs are validated before anything is written. The command fails when
headers differ, two inputs resolve to the same sample, a file mixes strains,
coordinates cannot be parsed, or one sample has two rows with the same site
key (never observed in real data; refusing beats silently picking one).
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path

from RelocaTE3 import logger

COHORT_TABLE = "cohort.characTErized.txt"
GENOTYPE_MATRIX = "cohort.genotype_matrix.tsv"
SUMMARY_MQC = "cohort.summary_mqc.tsv"

REQUIRED_COLUMNS = ("strain", "TE", "chromosome.pos", "strand", "status")
SITE_COLUMNS = ("chromosome", "start", "end", "strand", "TE")
# The classes RelocaTE3's status ladder emits (characterize.classify_status);
# they are always reported, even at zero, so summary columns are stable.
GENOTYPE_CLASSES = ("homozygous", "heterozygous", "somatic_insertion", "other")
NO_CALL = "."

_COORD = re.compile(r"^(?P<chrom>.+):(?P<start>\d+)\.\.(?P<end>\d+)$")
_DIGITS = re.compile(r"(\d+)")

# MultiQC custom content: YAML config in '# ' comment lines ahead of the TSV.
# Format: https://github.com/MultiQC/MultiQC/blob/main/docs/markdown/custom_content/index.md
# ("MultiQC-specific data file"; first table column = sample name).
_MQC_HEADER = """\
# id: 'relocate3_cohort_summary'
# section_name: 'RelocaTE3 cohort summary'
# description: 'Characterized non-reference insertions per sample. Genotype class is the part of the status before the first slash; excision counts statuses that mention excision.'
# format: 'tsv'
# plot_type: 'table'
# pconfig:
#     id: 'relocate3_cohort_summary_table'
#     title: 'RelocaTE3: characterized insertions per sample'
"""


class CohortError(ValueError):
    """Invalid or inconsistent cohort input."""


@dataclass
class SampleTable:
    """One sample's characterized table, rows kept verbatim."""

    sample: str
    path: Path
    header: str
    rows: list[str] = field(default_factory=list)
    sites: dict[tuple, str] = field(default_factory=dict)  # site key -> status
    statuses: list[str] = field(default_factory=list)


@dataclass
class Cohort:
    header: str
    tables: list[SampleTable]  # sorted by sample

    @property
    def samples(self) -> list[str]:
        return [t.sample for t in self.tables]


def natural_key(name: str) -> tuple:
    """Natural sort key for chromosome names: Chr1 < Chr2 < Chr10.

    The name is split into alternating text and digit runs; text compares as
    text (case-sensitive) and digit runs as integers. Names without digits sort
    by their text among the others (``Chr10`` < ``ChrSy`` < ``ChrUn`` <
    ``scaffold_9``). Ties such as ``Chr01``/``Chr1`` fall back to the raw name.
    """
    parts = _DIGITS.split(name)
    return tuple(int(p) if i % 2 else p for i, p in enumerate(parts)), name


def status_class(status: str) -> tuple[str, bool]:
    """``(genotype class, is_excision)`` for a characterize ``status`` value."""
    return status.split("/", 1)[0], "excision" in status


def parse_input_spec(spec: str) -> tuple[str | None, Path]:
    """``NAME=PATH`` -> (NAME, PATH); a plain path -> (None, PATH).

    An argument that is an existing file is always a plain path, so file names
    containing ``=`` still work.
    """
    if "=" in spec and not Path(spec).is_file():
        name, path = spec.split("=", 1)
        if not name:
            raise CohortError(f"empty sample name in input {spec!r}")
        return name, Path(path)
    return None, Path(spec)


def _read_table(spec: str) -> SampleTable:
    name, path = parse_input_spec(spec)
    if not path.is_file():
        raise CohortError(f"characterized table not found: {path}")
    with open(path) as fh:
        lines = fh.read().splitlines()
    if not lines:
        raise CohortError(f"{path}: file is empty (expected at least a header)")

    header, rows = lines[0], [line for line in lines[1:] if line != ""]
    columns = header.split("\t")
    missing = [c for c in REQUIRED_COLUMNS if c not in columns]
    if missing:
        raise CohortError(f"{path}: header lacks column(s) {', '.join(missing)}")
    idx = {c: columns.index(c) for c in REQUIRED_COLUMNS}

    split_rows = []
    strains: list[str] = []
    for lineno, line in enumerate(rows, start=2):
        cols = line.split("\t")
        if len(cols) != len(columns):
            raise CohortError(
                f"{path}:{lineno}: {len(cols)} fields, header has {len(columns)}"
            )
        split_rows.append(cols)
        if cols[idx["strain"]] not in strains:
            strains.append(cols[idx["strain"]])
    if len(strains) > 1:
        raise CohortError(
            f"{path}: more than one strain in a per-sample table ({', '.join(strains)})"
        )
    if name is None:
        if not strains:
            raise CohortError(
                f"{path}: no data rows, so the sample name cannot be read from the "
                "strain column; pass it as NAME=PATH"
            )
        name = strains[0]
    elif strains and strains[0] != name:
        raise CohortError(
            f"{path}: given as sample {name!r} but its strain column says {strains[0]!r}"
        )
    if not name or any(c in name for c in "\t\r\n") or name in SITE_COLUMNS:
        raise CohortError(f"{path}: unusable sample name {name!r}")

    table = SampleTable(sample=name, path=path, header=header, rows=rows)
    first_line: dict[tuple, int] = {}
    for lineno, cols in enumerate(split_rows, start=2):
        coord = _COORD.match(cols[idx["chromosome.pos"]])
        if not coord:
            raise CohortError(
                f"{path}:{lineno}: cannot parse chromosome.pos "
                f"{cols[idx['chromosome.pos']]!r} (expected CHROM:START..END)"
            )
        key = (
            coord["chrom"],
            int(coord["start"]),
            int(coord["end"]),
            cols[idx["strand"]],
            cols[idx["TE"]],
        )
        status = cols[idx["status"]]
        if key in table.sites:
            raise CohortError(
                f"duplicate site in sample {name!r}: {key[0]}:{key[1]}..{key[2]} "
                f"strand {key[3]} TE {key[4]} at {path} lines {first_line[key]} and "
                f"{lineno}; the genotype matrix needs one call per site and sample"
            )
        table.sites[key] = status
        first_line[key] = lineno
        table.statuses.append(status)
    return table


def build_cohort(specs: list[str]) -> Cohort:
    """Read and validate every input; nothing is written."""
    if not specs:
        raise CohortError("no characterized tables given")
    tables = [_read_table(spec) for spec in specs]

    header = tables[0].header
    for t in tables[1:]:
        if t.header != header:
            raise CohortError(
                f"header of {t.path} differs from header of {tables[0].path}"
            )
    by_sample: dict[str, SampleTable] = {}
    for t in tables:
        if t.sample in by_sample:
            raise CohortError(
                f"sample {t.sample!r} is claimed by more than one input: "
                f"{by_sample[t.sample].path} and {t.path}"
            )
        by_sample[t.sample] = t
    return Cohort(header=header, tables=[by_sample[s] for s in sorted(by_sample)])


def _site_sort_key(key: tuple) -> tuple:
    chrom, start, end, strand, te = key
    return natural_key(chrom), start, end, strand, te


def cohort_table_lines(cohort: Cohort) -> list[str]:
    return [cohort.header] + [row for t in cohort.tables for row in t.rows]


def genotype_matrix_lines(cohort: Cohort) -> list[str]:
    keys = sorted({k for t in cohort.tables for k in t.sites}, key=_site_sort_key)
    lines = ["\t".join([*SITE_COLUMNS, *cohort.samples])]
    for key in keys:
        cells = [t.sites.get(key, NO_CALL) for t in cohort.tables]
        lines.append("\t".join([str(v) for v in key] + cells))
    return lines


def summary_lines(cohort: Cohort) -> list[str]:
    counts = []
    extra: set[str] = set()
    for t in cohort.tables:
        c: dict[str, int] = {"total_calls": len(t.statuses), "excision": 0}
        for status in t.statuses:
            klass, excision = status_class(status)
            c[klass] = c.get(klass, 0) + 1
            c["excision"] += int(excision)
            if klass not in GENOTYPE_CLASSES:
                extra.add(klass)
        counts.append(c)
    columns = ["total_calls", *GENOTYPE_CLASSES, "excision", *sorted(extra)]
    lines = _MQC_HEADER.splitlines() + ["\t".join(["Sample", *columns])]
    for t, c in zip(cohort.tables, counts):
        lines.append("\t".join([t.sample] + [str(c.get(col, 0)) for col in columns]))
    return lines


def _write_atomic(path: Path, lines: list[str]) -> None:
    """Write via a hidden temp file + rename, so readers never see a partial file."""
    tmp = path.with_name(f".{path.name}.tmp")
    try:
        with open(tmp, "w") as fh:
            fh.write("".join(line + "\n" for line in lines))
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def write_cohort(specs: list[str], outdir: str | Path) -> dict[str, Path]:
    """Validate ``specs`` then write the three cohort files into ``outdir``."""
    cohort = build_cohort(specs)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        COHORT_TABLE: cohort_table_lines(cohort),
        GENOTYPE_MATRIX: genotype_matrix_lines(cohort),
        SUMMARY_MQC: summary_lines(cohort),
    }
    paths = {}
    for name, lines in outputs.items():
        path = outdir / name
        if path.exists():
            logger.info("cohort: replacing existing %s", path)
        _write_atomic(path, lines)
        paths[name] = path
    n_rows = sum(len(t.rows) for t in cohort.tables)
    logger.info(
        "cohort: %d sample(s), %d call(s), %d site(s) -> %s",
        len(cohort.tables),
        n_rows,
        len(outputs[GENOTYPE_MATRIX]) - 1,
        outdir,
    )
    return paths
