#!/usr/bin/env python3
"""Check a finished small long-read panel and write qc/panel_check.tsv.

Checks: manifest paths exist with completion markers; every FASTQ is a
complete gzip stream; delivered coverage is within tolerance; the truth
table has 200 rows; and the fixture's original short reads fit the rebuilt
sample genome (few clipped reads) far better than the reference.
Needs minimap2 on PATH for the last check.

``<panel>/.complete`` is written only when every check passes. A panel that
already has it is refused; delete the marker to recheck.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import subprocess
import sys
import tempfile
import zlib
from collections.abc import Iterable
from pathlib import Path

COVERAGE_TOLERANCE = 0.25  # HiFi loses about 8% of molecules in ccs
MAX_CLIP_RATIO = 0.05
MIN_CLIP = 20
EXPECTED_EVENTS = 200
MANIFEST_COLUMNS = ("sample", "platform", "coverage", "reads", "control_reads")
REPORT_HEADER = (
    "platform",
    "sample",
    "kind",
    "reads",
    "bases",
    "coverage",
    "n50",
    "max_length",
    "status",
)
STDERR_TAIL_LINES = 20


def fastq_stats(path: Path) -> dict[str, int]:
    """Read a whole gzipped FASTQ and summarise its read lengths.

    Reading to the end also proves the gzip stream is complete.

    Args:
        path: Gzipped four-line FASTQ.

    Returns:
        ``reads``, ``bases``, ``n50`` and ``max_length``.

    Raises:
        ValueError: If a record is truncated or malformed, its sequence and
            quality lengths differ, or the file holds no reads.
        OSError: If the file is missing or not valid gzip.
        EOFError: If the gzip stream is cut short.
    """
    lengths = []
    with gzip.open(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip("\n")
            plus = handle.readline()
            quality = handle.readline().rstrip("\n")
            if not quality or not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"truncated or malformed FASTQ: {path}")
            if len(sequence) != len(quality):
                raise ValueError(
                    f"sequence and quality length differ ({len(sequence)} vs "
                    f"{len(quality)}) for {header.split()[0]} in {path}"
                )
            lengths.append(len(sequence))
    if not lengths:
        raise ValueError(f"no reads in FASTQ: {path}")
    lengths.sort(reverse=True)
    total = sum(lengths)
    running = 0
    n50 = 0
    for length in lengths:
        running += length
        if running * 2 >= total:
            n50 = length
            break
    return {"reads": len(lengths), "bases": total, "n50": n50, "max_length": lengths[0]}


def count_clipped_reads(sam_lines: Iterable[str], min_clip: int) -> int:
    """Count mapped primary records with a soft clip of at least ``min_clip``.

    Args:
        sam_lines: SAM lines without trailing newlines; headers are skipped.
        min_clip: Smallest soft clip (bases, either end) that counts.

    Returns:
        The number of clipped primary alignments.
    """
    count = 0
    for line in sam_lines:
        if line.startswith("@"):
            continue
        fields = line.split("\t")
        flag = int(fields[1])
        if flag & 0x904:  # unmapped, secondary or supplementary
            continue
        clips = [int(n) for n in re.findall(r"(\d+)S", fields[5])]
        if clips and max(clips) >= min_clip:
            count += 1
    return count


def coverage_ok(bases: int, genome_length: int, requested: float) -> bool:
    """Decide whether delivered coverage is within tolerance of the request.

    Args:
        bases: Bases delivered.
        genome_length: Length of the genome the reads were simulated from.
        requested: Requested depth.

    Returns:
        True if ``|bases / genome_length - requested|`` is at most
        ``COVERAGE_TOLERANCE * requested``.
    """
    delivered = bases / genome_length
    return abs(delivered - requested) <= COVERAGE_TOLERANCE * requested


def clipped_short_reads(genome: Path, r1: Path, r2: Path, threads: int) -> int:
    """Map paired short reads to ``genome`` and count clipped reads.

    Args:
        genome: FASTA to map to (indexed on the fly by minimap2).
        r1: Read-1 FASTQ.
        r2: Read-2 FASTQ.
        threads: minimap2 threads.

    Returns:
        Primary alignments soft-clipped by at least ``MIN_CLIP`` bases.

    Raises:
        RuntimeError: If minimap2 exits non-zero; the message holds the
            tail of its stderr.
    """
    with tempfile.TemporaryFile(mode="w+") as stderr:
        process = subprocess.Popen(
            [
                "minimap2",
                "-ax",
                "sr",
                "-t",
                str(threads),
                str(genome),
                str(r1),
                str(r2),
            ],
            stdout=subprocess.PIPE,
            stderr=stderr,
            text=True,
        )
        assert process.stdout is not None
        count = count_clipped_reads(
            (line.rstrip("\n") for line in process.stdout), MIN_CLIP
        )
        if process.wait() != 0:
            stderr.seek(0)
            tail = "".join(stderr.readlines()[-STDERR_TAIL_LINES:])
            raise RuntimeError(
                f"minimap2 failed (exit {process.returncode}) mapping {r1} and "
                f"{r2} to {genome}; stderr tail:\n{tail}"
            )
    return count


def fasta_length(path: Path) -> int:
    """Total sequence length of a FASTA.

    Args:
        path: Uncompressed FASTA.

    Returns:
        Bases over all records.

    Raises:
        OSError: If the file cannot be read.
    """
    with open(path) as handle:
        return sum(len(line.strip()) for line in handle if not line.startswith(">"))


def read_manifest(path: Path) -> list[dict[str, str]]:
    """Read ``panel_manifest.tsv``.

    Args:
        path: Manifest written by ``simulate_reads.py``.

    Returns:
        One dict per row.

    Raises:
        OSError: If the file cannot be read.
        ValueError: If a column in ``MANIFEST_COLUMNS`` is missing or the
            manifest has no rows.
    """
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [c for c in MANIFEST_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path} lacks columns: {', '.join(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path} has no rows")
    return rows


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point.

    Every problem found is printed as a ``FAIL:`` line; missing or
    unreadable inputs are failures, not tracebacks.

    Args:
        argv: Arguments without the program name; ``None`` reads ``sys.argv``.

    Returns:
        0 if every check passed (``<panel>/.complete`` is then written),
        1 otherwise or if the panel was already checked.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--short-r1", type=Path, required=True)
    parser.add_argument("--short-r2", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args(argv)

    marker = args.panel / ".complete"
    if marker.exists():
        print(
            f"ERROR: panel already checked ({marker} exists); delete the marker "
            f"to recheck",
            file=sys.stderr,
        )
        return 1

    failures = []
    rows = []

    genome_length = None
    try:
        genome_length = fasta_length(args.reference)
        if genome_length == 0:
            failures.append(f"reference has no sequence: {args.reference}")
            genome_length = None
    except OSError as error:
        failures.append(f"cannot read reference {args.reference}: {error}")

    truth_path = args.panel / "truth_events.tsv"
    try:
        events = len(truth_path.read_text().splitlines()) - 1
        if events != EXPECTED_EVENTS:
            failures.append(
                f"{truth_path} has {events} events, expected {EXPECTED_EVENTS}"
            )
    except OSError as error:
        failures.append(f"cannot read truth table {truth_path}: {error}")

    manifest_path = args.panel / "panel_manifest.tsv"
    try:
        manifest = read_manifest(manifest_path)
    except (OSError, ValueError) as error:
        failures.append(f"bad manifest {manifest_path}: {error}")
        manifest = []
    for entry in manifest:
        for column in ("reads", "control_reads"):
            fastq = args.panel / entry[column]
            if not (fastq.parent / ".complete").is_file():
                failures.append(f"no completion marker: {fastq.parent}")
                continue
            try:
                stats = fastq_stats(fastq)
            except (OSError, EOFError, ValueError, zlib.error) as error:
                failures.append(f"unreadable FASTQ {fastq}: {error}")
                continue
            try:
                requested = float(entry["coverage"])
            except ValueError:
                failures.append(f"bad coverage {entry['coverage']!r} for {fastq}")
                continue
            if genome_length is None:
                ok = False
                delivered = "NA"
            else:
                ok = coverage_ok(stats["bases"], genome_length, requested)
                delivered = f"{stats['bases'] / genome_length:.2f}"
                if not ok:
                    failures.append(
                        f"coverage off target: {fastq} ({delivered}x delivered, "
                        f"{requested:g}x requested)"
                    )
            rows.append(
                (
                    entry["platform"],
                    entry["sample"],
                    column,
                    stats["reads"],
                    stats["bases"],
                    delivered,
                    stats["n50"],
                    stats["max_length"],
                    "ok" if ok else "FAIL",
                )
            )

    clipped_reference = clipped_sample = "NA"
    sample_genome = args.panel / "truth" / "genomes" / "sample_h1.fa"
    inputs = (args.reference, sample_genome, args.short_r1, args.short_r2)
    absent = [str(path) for path in inputs if not path.is_file()]
    if absent:
        failures.append(f"short-read fit not checked; missing: {', '.join(absent)}")
    else:
        try:
            clipped_reference = clipped_short_reads(
                args.reference, args.short_r1, args.short_r2, args.threads
            )
            clipped_sample = clipped_short_reads(
                sample_genome, args.short_r1, args.short_r2, args.threads
            )
        except (OSError, RuntimeError) as error:
            failures.append(f"short-read fit not checked: {error}")
        else:
            if clipped_sample > MAX_CLIP_RATIO * clipped_reference:
                failures.append(
                    f"rebuilt genome does not explain the original short reads: "
                    f"{clipped_sample} clipped reads vs {clipped_reference} on "
                    f"the reference (max ratio {MAX_CLIP_RATIO})"
                )

    qc = args.panel / "qc"
    if args.panel.is_dir():
        qc.mkdir(exist_ok=True)
        with open(qc / "panel_check.tsv", "w") as handle:
            handle.write("\t".join(REPORT_HEADER) + "\n")
            for row in rows:
                handle.write("\t".join(str(value) for value in row) + "\n")
            handle.write(f"# clipped_short_reads_reference\t{clipped_reference}\n")
            handle.write(f"# clipped_short_reads_sample\t{clipped_sample}\n")
    else:
        failures.append(f"panel directory not found: {args.panel}")

    for failure in failures:
        print(f"FAIL: {failure}", file=sys.stderr)
    if failures:
        print(f"RESULT: FAIL ({len(failures)} problems)", file=sys.stderr)
        return 1
    marker.touch()
    print(f"RESULT: PASS ({len(rows)} FASTQs checked)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
