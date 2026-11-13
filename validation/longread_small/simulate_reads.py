#!/usr/bin/env python3
"""Simulate ONT-HQ and HiFi reads for the small Chr3 long-read panel.

Run inside the simulate_data toolkit's pixi environment (PBSIM3 + ccs). The
toolkit is imported only inside `simulate`, so the job list, input checks,
manifest and skip rules can be tested without it.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import hashlib
import io
import sys
import tempfile
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

MANIFEST_HEADER = (
    "sample",
    "platform",
    "coverage",
    "replicate",
    "reads",
    "control_reads",
)

SIMULATION_HEADER = (
    "platform",
    "coverage",
    "seed",
    "genome",
    "genome_bytes",
    "genome_sha256",
    "reads",
    "ccs_yield",
)

# Fields a completed directory must match before a rerun may skip it.
MATCHED_FIELDS = ("platform", "coverage", "seed", "genome_sha256")

# The toolkit's PLATFORMS keys this driver supports; checked without importing it.
SUPPORTED_PLATFORMS = ("ont-hq", "hifi")


@dataclass(frozen=True)
class Job:
    """One platform x coverage x replicate sample, plus its reference control."""

    platform: str
    coverage: int
    replicate: int
    seed: int
    control_seed: int

    @property
    def sample(self) -> str:
        """Sample name, ``cov<coverage>x_rep<replicate>``."""
        return f"cov{self.coverage}x_rep{self.replicate}"

    @property
    def reads(self) -> str:
        """Sample FASTQ path, relative to the panel root."""
        return f"reads/{self.platform}/{self.sample}/reads.fastq.gz"

    @property
    def control_reads(self) -> str:
        """Reference-control FASTQ path, relative to the panel root."""
        return f"reads/{self.platform}/{self.sample}_reference_control/reads.fastq.gz"


def validate_inputs(
    platforms: list[str],
    coverages: list[int],
    replicates: int,
    genomes: list[Path],
) -> None:
    """Check command-line inputs before any simulation starts.

    Args:
        platforms: Platform keys; each must be in ``SUPPORTED_PLATFORMS``.
        coverages: Target depths; each must be an integer >= 1.
        replicates: Replicates per platform and coverage; must be >= 1.
        genomes: FASTA files the reads are simulated from.

    Raises:
        ValueError: On an unknown or repeated platform, a repeated coverage,
            a coverage below 1, or fewer than 1 replicate.
        FileNotFoundError: If a genome is not an existing file.
    """
    seen_platforms: set[str] = set()
    for platform in platforms:
        if platform not in SUPPORTED_PLATFORMS:
            raise ValueError(
                f"unknown platform {platform!r}; expected one of "
                f"{', '.join(SUPPORTED_PLATFORMS)}"
            )
        if platform in seen_platforms:
            raise ValueError(f"repeated platform {platform!r}")
        seen_platforms.add(platform)
    seen_coverages: set[int] = set()
    for coverage in coverages:
        if coverage < 1:
            raise ValueError(f"coverage must be >= 1, got {coverage}")
        if coverage in seen_coverages:
            raise ValueError(f"repeated coverage {coverage}")
        seen_coverages.add(coverage)
    if replicates < 1:
        raise ValueError(f"replicates must be >= 1, got {replicates}")
    for genome in genomes:
        if not Path(genome).is_file():
            raise FileNotFoundError(f"genome FASTA not found: {genome}")


def jobs(
    platforms: list[str], coverages: list[int], replicates: int, seed: int
) -> list[Job]:
    """List jobs in platform, coverage, replicate order with distinct seeds.

    Args:
        platforms: Platform keys, in output order.
        coverages: Target depths, in output order.
        replicates: Replicates per platform and coverage.
        seed: Base seed; job ``i`` gets ``seed + 2i`` for the sample and
            ``seed + 2i + 1`` for its control.

    Returns:
        One Job per platform x coverage x replicate.
    """
    out = []
    for platform in platforms:
        for coverage in coverages:
            for replicate in range(1, replicates + 1):
                index = len(out)
                out.append(
                    Job(
                        platform,
                        coverage,
                        replicate,
                        seed + 2 * index,
                        seed + 2 * index + 1,
                    )
                )
    return out


def write_manifest(job_list: list[Job], path: Path) -> None:
    """Write the manifest with the full long-read panel's columns.

    Args:
        job_list: Jobs, one row each, in order.
        path: Output TSV; overwritten (it is derived only from ``job_list``).
    """
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(MANIFEST_HEADER)
        for job in job_list:
            writer.writerow(
                (
                    job.sample,
                    job.platform,
                    job.coverage,
                    job.replicate,
                    job.reads,
                    job.control_reads,
                )
            )


def needs_run(directory: Path) -> bool:
    """Decide whether a read directory still has to be simulated.

    Output without a completion marker is never overwritten or deleted.

    Args:
        directory: One sample or control read directory.

    Returns:
        True if the directory is absent or empty; False if it holds a
        ``.complete`` marker.

    Raises:
        RuntimeError: If the directory holds files but no ``.complete``
            marker; the message lists them and the command to move them aside.
    """
    if not directory.exists() or not any(directory.iterdir()):
        return True
    if (directory / ".complete").is_file():
        return False
    found = ", ".join(sorted(entry.name for entry in directory.iterdir()))
    raise RuntimeError(
        f"incomplete output in {directory} (found: {found}); a previous run "
        f"failed or was killed. Move it aside and rerun: "
        f"mv {directory} {directory}.failed"
    )


def genome_identity(genome: Path) -> tuple[int, str]:
    """Return a genome file's size and content digest.

    Args:
        genome: FASTA file.

    Returns:
        ``(size in bytes, SHA-256 hex digest)``.
    """
    digest = hashlib.sha256()
    with open(genome, "rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return Path(genome).stat().st_size, digest.hexdigest()


def write_simulation_record(
    directory: Path,
    platform_key: str,
    coverage: float,
    seed: int,
    genome: Path,
    reads: int,
    ccs_yield: float | None,
) -> None:
    """Write ``simulation.tsv``: what a read directory was simulated from.

    Args:
        directory: The read directory.
        platform_key: Platform key.
        coverage: Target depth.
        seed: PBSIM3 seed.
        genome: Source FASTA, recorded as given plus its size and SHA-256.
        reads: Number of reads written.
        ccs_yield: Mean ccs ZMW yield, or None when ccs was not run.
    """
    size, digest = genome_identity(genome)
    row = (
        platform_key,
        f"{coverage:g}",
        str(seed),
        str(genome),
        str(size),
        digest,
        str(reads),
        "NA" if ccs_yield is None else f"{ccs_yield:.4f}",
    )
    with open(directory / "simulation.tsv", "w") as handle:
        handle.write("\t".join(SIMULATION_HEADER) + "\n")
        handle.write("\t".join(row) + "\n")


def check_complete(
    directory: Path,
    platform_key: str,
    coverage: float,
    seed: int,
    genome: Path,
) -> None:
    """Check that a completed read directory matches the current job.

    Nothing is overwritten or deleted, whatever the outcome.

    Args:
        directory: A read directory holding a ``.complete`` marker.
        platform_key: Platform the job asks for.
        coverage: Depth the job asks for.
        seed: Seed the job asks for.
        genome: FASTA the job would simulate from; compared by SHA-256.

    Raises:
        RuntimeError: If ``simulation.tsv`` is missing, is not one header row
            plus one record with ``SIMULATION_HEADER`` columns, or differs
            from the job in any of ``MATCHED_FIELDS``.
    """
    path = directory / "simulation.tsv"
    if not path.is_file():
        raise RuntimeError(
            f"{directory} is marked complete but has no simulation.tsv; "
            f"cannot tell what it holds. Move it aside: mv {directory} "
            f"{directory}.unverified"
        )
    lines = path.read_text().splitlines()
    if len(lines) != 2 or tuple(lines[0].split("\t")) != SIMULATION_HEADER:
        raise RuntimeError(
            f"unreadable {path}: expected a header {SIMULATION_HEADER} and one "
            f"record. Move the directory aside: mv {directory} "
            f"{directory}.unverified"
        )
    values = lines[1].split("\t")
    if len(values) != len(SIMULATION_HEADER):
        raise RuntimeError(
            f"unreadable {path}: record has {len(values)} fields, expected "
            f"{len(SIMULATION_HEADER)}"
        )
    found = dict(zip(SIMULATION_HEADER, values, strict=True))
    wanted = {
        "platform": platform_key,
        "coverage": f"{coverage:g}",
        "seed": str(seed),
        "genome_sha256": genome_identity(genome)[1],
    }
    differing = [
        f"{field} (completed {found[field]}, requested {wanted[field]})"
        for field in MATCHED_FIELDS
        if found[field] != wanted[field]
    ]
    if differing:
        raise RuntimeError(
            f"{directory} was completed for a different job: "
            f"{'; '.join(differing)}. Rerun with the original arguments or "
            f"use a new --panel; it is never overwritten"
        )


@contextlib.contextmanager
def open_pooled(path: Path, compresslevel: int = 6) -> Iterator[io.TextIOWrapper]:
    """Open a gzip text stream whose bytes depend only on what is written.

    The gzip header gets mtime 0 and no file name, so the same reads give a
    byte-identical file.

    Args:
        path: Output file; overwritten.
        compresslevel: gzip level.

    Yields:
        A text stream writing gzip-compressed bytes to ``path``.
    """
    with open(path, "wb") as raw:
        with gzip.GzipFile(
            filename="", mode="wb", compresslevel=compresslevel, fileobj=raw, mtime=0
        ) as compressed:
            with io.TextIOWrapper(compressed, encoding="ascii", newline="\n") as text:
                yield text


def simulate(
    platform_key: str,
    genome: Path,
    coverage: float,
    seed: int,
    directory: Path,
    name_prefix: str,
    threads: int,
    min_ccs_yield: float,
) -> int:
    """Simulate one pooled FASTQ into ``directory``.

    Reads are pooled into ``reads.fastq.gz.partial``, renamed to
    ``reads.fastq.gz`` when done, then ``simulation.tsv`` (see
    ``write_simulation_record``) and ``.complete`` are written. A failure leaves the ``.partial`` file in place so
    ``needs_run`` refuses the directory until it is moved aside.

    Args:
        platform_key: Key into the toolkit's ``PLATFORMS``.
        genome: FASTA to simulate from.
        coverage: Target depth.
        seed: PBSIM3 seed.
        directory: Output directory; created if absent.
        name_prefix: Prefix for every read name (must not contain spaces).
        threads: Threads for ccs (HiFi only).
        min_ccs_yield: Fail if ccs passes fewer than this fraction of ZMWs.

    Returns:
        The number of reads written.

    Raises:
        RuntimeError: If PBSIM3 or ccs fails, ccs yield is below
            ``min_ccs_yield``, or no reads were simulated.
        FileNotFoundError: If PBSIM3 output is missing.
    """
    from simulate_data.longread import (
        GZIP_LEVEL,
        PLATFORMS,
        build_pbsim_command,
        pool_fastq_gz,
        run_ccs,
    )

    # Private toolkit helper, verified against simulate_data commit 1e19882
    # (working tree dirty); the toolkit's own drivers import it the same way.
    from simulate_data.modules.te_benchmark_panel import _pbsim_contig_outputs
    from simulate_data.utils import run_command

    platform = PLATFORMS[platform_key]
    directory.mkdir(parents=True, exist_ok=True)
    partial = directory / "reads.fastq.gz.partial"
    count = 0
    yields = []
    with tempfile.TemporaryDirectory(dir=directory) as tmp:
        prefix = Path(tmp) / "sim"
        run_command(build_pbsim_command(platform, genome, coverage, prefix, seed))
        with open_pooled(partial, GZIP_LEVEL) as pooled:
            for _ref, _maf, payload in _pbsim_contig_outputs(prefix):
                if payload.suffix == ".bam":
                    fastq = Path(str(payload)[: -len(".bam")] + ".fastq.gz")
                    stats = run_ccs(payload, fastq, threads, min_ccs_yield)
                    yields.append(stats["yield_fraction"])
                else:
                    fastq = payload
                count += pool_fastq_gz(fastq, pooled, name_prefix)
    if count == 0:
        raise RuntimeError(f"no reads simulated for {directory}")
    partial.rename(directory / "reads.fastq.gz")
    mean_yield = sum(yields) / len(yields) if yields else None
    write_simulation_record(
        directory, platform_key, coverage, seed, genome, count, mean_yield
    )
    (directory / ".complete").touch()
    return count


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point.

    Args:
        argv: Arguments without the program name; ``None`` reads ``sys.argv``.

    Returns:
        The process exit status, 0 on success.

    Raises:
        ValueError: On invalid inputs (see ``validate_inputs``).
        FileNotFoundError: If a genome FASTA is missing.
        RuntimeError: On incomplete earlier output, completed output from a
            different job (see ``check_complete``), or a failed simulation.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True, help="panel root")
    parser.add_argument("--sample-genome", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--platform", action="append", default=None)
    parser.add_argument("--coverage", action="append", type=int, default=None)
    parser.add_argument("--replicates", type=int, default=1)
    parser.add_argument("--seed", type=int, default=916)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--min-ccs-yield", type=float, default=0.80)
    args = parser.parse_args(argv)

    platforms = args.platform or ["ont-hq", "hifi"]
    coverages = args.coverage or [5, 20]
    validate_inputs(
        platforms, coverages, args.replicates, [args.sample_genome, args.reference]
    )
    job_list = jobs(platforms, coverages, args.replicates, args.seed)
    for job in job_list:
        for genome, relative, seed, label in (
            (args.sample_genome, job.reads, job.seed, "sample"),
            (args.reference, job.control_reads, job.control_seed, "control"),
        ):
            directory = (args.panel / relative).parent
            if not needs_run(directory):
                check_complete(directory, job.platform, job.coverage, seed, genome)
                print(f"skip (complete): {directory}")
                continue
            count = simulate(
                job.platform,
                genome,
                job.coverage,
                seed,
                directory,
                f"{job.sample}:{label}:s{seed}",
                args.threads,
                args.min_ccs_yield,
            )
            print(f"wrote {count} reads: {directory}")
    write_manifest(job_list, args.panel / "panel_manifest.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
