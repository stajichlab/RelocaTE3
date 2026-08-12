"""Resolve a batch of samples from a sample sheet or a directory of FASTQs.

RelocaTE2 accepted ``--fq_dir``; RelocaTE3 took one sample per invocation, so
running a cohort meant hand-writing a loop. Both inputs here resolve to the same
explicit :class:`Sample` list, which ``relocaTE3 run-batch`` then feeds through
the ordinary single-sample pipeline one sample at a time.

Sheet columns follow the schema documented for the planned Nextflow work in
``plans/FEATURES.md`` -- ``sample_id, r1_fq, r2_fq`` plus optional per-row
``te_library`` / ``reference_genome`` / ``repeatmasker`` overrides -- so a sheet
written for one entry point works for the other.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

# FASTQ extensions we recognise when scanning a directory
_FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

# trailing mate markers: sample_R1 / sample_1 / sample.R1 / sample-R1
_MATE_MARKER = re.compile(r"(?P<sep>[._-])(?:R)?(?P<mate>[12])$", re.IGNORECASE)

# accepted column spellings -> canonical field
_COLUMNS = {
    "sample_id": "name",
    "sample": "name",
    "name": "name",
    "r1_fq": "r1",
    "r1": "r1",
    "fastq_1": "r1",
    "r2_fq": "r2",
    "r2": "r2",
    "fastq_2": "r2",
    "te_library": "te_library",
    "reference_genome": "genome",
    "genome": "genome",
    "repeatmasker": "repeatmasker",
    "repeatmasker_out": "repeatmasker",
}


@dataclass(frozen=True)
class Sample:
    """One sample to run: a name, its reads, and any per-sample overrides."""

    name: str
    r1: str
    r2: str | None = None
    te_library: str | None = None
    genome: str | None = None
    repeatmasker: str | None = None

    @property
    def is_paired(self) -> bool:
        """True when a second read file was supplied."""
        return bool(self.r2)


def _strip_fastq_suffix(filename: str) -> str | None:
    """Return ``filename`` without its FASTQ extension, or None if not a FASTQ."""
    lowered = filename.lower()
    for suffix in _FASTQ_SUFFIXES:
        if lowered.endswith(suffix):
            return filename[: -len(suffix)]
    return None


def _require_file(path: str, sample: str, field: str) -> str:
    resolved = Path(path).expanduser()
    if not resolved.is_file():
        raise ValueError(
            f"sample '{sample}': {field} does not exist or is not a file: {path}"
        )
    return str(resolved)


def read_sample_sheet(path: str | Path) -> list[Sample]:
    """Parse a CSV/TSV sample sheet into :class:`Sample` records.

    The delimiter is inferred from the file extension (``.tsv``/``.tab`` are tab
    separated, anything else comma). Blank lines and ``#`` comments are skipped.
    Every referenced read file must exist -- a batch that dies twenty samples in
    on a typo is worse than one that refuses to start.
    """
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"sample sheet not found: {path}")

    delimiter = "\t" if path.suffix.lower() in (".tsv", ".tab") else ","
    rows: list[dict[str, str]] = []
    with open(path, newline="") as handle:
        lines = [
            line
            for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        ]
        if not lines:
            raise ValueError(f"sample sheet is empty: {path}")
        reader = csv.DictReader(lines, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"sample sheet has no header row: {path}")

        mapping = {}
        for column in reader.fieldnames:
            key = (column or "").strip().lower()
            if key in _COLUMNS:
                mapping[column] = _COLUMNS[key]
        if "name" not in mapping.values():
            raise ValueError(
                f"sample sheet {path} needs a sample_id column "
                f"(got: {', '.join(reader.fieldnames)})"
            )
        if "r1" not in mapping.values():
            raise ValueError(
                f"sample sheet {path} needs an r1_fq column "
                f"(got: {', '.join(reader.fieldnames)})"
            )
        for row in reader:
            rows.append(
                {
                    field: (row.get(column) or "").strip()
                    for column, field in mapping.items()
                }
            )

    samples: list[Sample] = []
    seen: set[str] = set()
    for row in rows:
        name = row.get("name", "")
        if not name:
            continue
        if name in seen:
            raise ValueError(f"sample sheet {path}: duplicate sample id '{name}'")
        seen.add(name)
        if not row.get("r1"):
            raise ValueError(f"sample '{name}': r1_fq is empty")
        samples.append(
            Sample(
                name=name,
                r1=_require_file(row["r1"], name, "r1_fq"),
                r2=_require_file(row["r2"], name, "r2_fq") if row.get("r2") else None,
                te_library=(
                    _require_file(row["te_library"], name, "te_library")
                    if row.get("te_library")
                    else None
                ),
                genome=(
                    _require_file(row["genome"], name, "reference_genome")
                    if row.get("genome")
                    else None
                ),
                repeatmasker=(
                    _require_file(row["repeatmasker"], name, "repeatmasker")
                    if row.get("repeatmasker")
                    else None
                ),
            )
        )
    if not samples:
        raise ValueError(f"sample sheet {path} lists no samples")
    return samples


def discover_fastq_dir(directory: str | Path) -> list[Sample]:
    """Pair the FASTQs in ``directory`` into samples (RelocaTE2's ``--fq_dir``).

    Recognises ``sample_R1``/``sample_R2`` and ``sample_1``/``sample_2`` (also
    with ``.`` or ``-`` before the marker). A FASTQ with no mate marker at all is
    taken as single-end; a file that *does* carry a marker but whose partner is
    missing raises, because a half-present pair is much more often a staging
    mistake than a deliberate single-end run.
    """
    directory = Path(directory).expanduser()
    if not directory.is_dir():
        raise ValueError(f"--fq-dir is not a directory: {directory}")

    stems: dict[str, dict[str, str]] = {}
    singles: dict[str, str] = {}
    for entry in sorted(directory.iterdir()):
        if not entry.is_file():
            continue
        stem = _strip_fastq_suffix(entry.name)
        if stem is None:
            continue
        match = _MATE_MARKER.search(stem)
        if match:
            base = stem[: match.start()]
            stems.setdefault(base, {})[match.group("mate")] = str(entry)
        else:
            singles[stem] = str(entry)

    if not stems and not singles:
        raise ValueError(f"no FASTQ files found in {directory}")

    samples: list[Sample] = []
    for base, mates in stems.items():
        if "1" not in mates:
            raise ValueError(
                f"sample '{base}' in {directory} has an R2 file but no R1"
            )
        if "2" not in mates:
            raise ValueError(
                f"sample '{base}' in {directory} has an R1 file but no R2. "
                "Move it aside, or pass it explicitly with run-all --left for a "
                "genuine single-end run."
            )
        samples.append(Sample(name=base, r1=mates["1"], r2=mates["2"]))
    samples.extend(Sample(name=base, r1=path) for base, path in singles.items())

    duplicates = {s.name for s in samples}
    if len(duplicates) != len(samples):
        raise ValueError(f"duplicate sample names discovered in {directory}")
    return sorted(samples, key=lambda s: s.name)
