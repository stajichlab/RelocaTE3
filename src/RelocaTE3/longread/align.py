"""Long-read alignment to the genome and checks on user-supplied BAM/CRAM.

Reads are aligned with ``minimap2`` using the preset of their
:class:`~RelocaTE3.longread.model.ReadType` and piped straight into
``samtools sort``; the result is indexed and renamed into place only when
every step succeeded, so a final BAM is never half-written.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import IO

import pysam

from RelocaTE3.longread.model import ReadType
from RelocaTE3.longread.reference import ReferenceSpace

_STDERR_TAIL = 20
_INDEX_SUFFIXES = (".bai", ".csi", ".crai")


def _check_threads(threads: int) -> None:
    if threads < 1:
        raise ValueError(f"threads must be at least 1, got {threads}")


def minimap2_command(
    genome: Path, reads: list[Path], read_type: ReadType, threads: int
) -> list[str]:
    """Build the minimap2 command that aligns long reads to the genome.

    Args:
        genome: Genome FASTA (or minimap2 ``.mmi`` index).
        reads: Read files (FASTQ/FASTA, optionally gzipped).
        read_type: Selects the ``-x`` preset.
        threads: minimap2 threads (``-t``); at least 1.

    Returns:
        The argument list, SAM to stdout with soft-clipped supplementaries
        (``-Y``), ``MD`` tags and no secondary alignments.

    Raises:
        ValueError: If ``threads`` is below 1.
    """
    _check_threads(threads)
    return [
        "minimap2", "-a", "-x", read_type.minimap2_preset, "-Y", "--MD",
        "--secondary=no", "-t", str(threads), str(genome), *map(str, reads),
    ]  # fmt: skip


def _tail(fh: IO[bytes]) -> str:
    fh.seek(0)
    lines = fh.read().decode(errors="replace").splitlines()
    return "\n".join(lines[-_STDERR_TAIL:])


def _partial_path(out_bam: Path) -> Path:
    stem = out_bam.name[:-4] if out_bam.name.endswith(".bam") else out_bam.name
    return out_bam.with_name(stem + ".partial.bam")


def align_reads(
    genome: Path, reads: list[Path], read_type: ReadType, out_bam: Path, threads: int
) -> Path:
    """Align long reads to a coordinate-sorted, indexed BAM.

    Runs ``minimap2 ... | samtools sort -@ threads -o <out>.partial.bam -``,
    indexes the partial BAM, then renames the BAM and its ``.bai`` to
    ``out_bam`` and ``out_bam.bai``.

    Args:
        genome: Genome FASTA (or minimap2 ``.mmi`` index).
        reads: Read files.
        read_type: Selects the minimap2 preset.
        out_bam: Final BAM path; must not exist.
        threads: Threads for minimap2 and samtools sort; at least 1.

    Returns:
        ``out_bam``.

    Raises:
        ValueError: If ``threads`` is below 1.
        FileExistsError: If ``out_bam`` or its index already exists.
        RuntimeError: If minimap2 or samtools is not on PATH, or a step
            fails (the message holds the last lines of that tool's stderr).
    """
    _check_threads(threads)
    out_bam = Path(out_bam)
    out_bai = Path(str(out_bam) + ".bai")
    for path in (out_bam, out_bai):
        if path.exists():
            raise FileExistsError(f"refusing to overwrite {path}")
    for tool in ("minimap2", "samtools"):
        if shutil.which(tool) is None:
            raise RuntimeError(f"{tool} not found on PATH")

    partial = _partial_path(out_bam)
    partial_bai = Path(str(partial) + ".bai")

    def cleanup() -> None:
        for p in (partial, partial_bai):
            p.unlink(missing_ok=True)

    cleanup()  # stale partials from an interrupted run
    mm2 = minimap2_command(genome, reads, read_type, threads)
    sort = ["samtools", "sort", "-@", str(threads), "-o", str(partial), "-"]
    # stderr goes to files, not pipes, so a chatty tool cannot deadlock.
    with tempfile.TemporaryFile() as mm2_err, tempfile.TemporaryFile() as sort_err:
        p_mm2 = subprocess.Popen(mm2, stdout=subprocess.PIPE, stderr=mm2_err)
        try:
            p_sort = subprocess.Popen(sort, stdin=p_mm2.stdout, stderr=sort_err)
        except BaseException:
            p_mm2.kill()
            p_mm2.wait()
            raise
        finally:
            assert p_mm2.stdout is not None
            p_mm2.stdout.close()  # sort holds the read end; mm2 sees SIGPIPE
        rc_sort = p_sort.wait()
        rc_mm2 = p_mm2.wait()
        if rc_mm2 != 0 or rc_sort != 0:
            cleanup()
            # minimap2 failing makes sort fail too; report the root cause.
            cmd, rc, err = (
                (mm2, rc_mm2, mm2_err) if rc_mm2 != 0 else (sort, rc_sort, sort_err)
            )
            raise RuntimeError(
                f"{cmd[0]} {'sort ' if cmd is sort else ''}failed (exit {rc}): "
                f"{' '.join(cmd)}\n{_tail(err)}"
            )

    idx = subprocess.run(
        ["samtools", "index", str(partial)], capture_output=True, text=True
    )
    if idx.returncode != 0:
        cleanup()
        tail = "\n".join(idx.stderr.splitlines()[-_STDERR_TAIL:])
        raise RuntimeError(f"samtools index failed (exit {idx.returncode})\n{tail}")
    os.replace(partial_bai, out_bai)
    os.replace(partial, out_bam)
    return out_bam


def _has_index(path: Path) -> bool:
    candidates = [Path(str(path) + s) for s in _INDEX_SUFFIXES]
    candidates += [path.with_suffix(s) for s in _INDEX_SUFFIXES]
    return any(c.is_file() for c in candidates)


def check_alignment_input(
    bam: Path, space: ReferenceSpace, reference_fasta: Path | None = None
) -> None:
    """Check that a user-supplied BAM/CRAM can be used as long-read input.

    Args:
        bam: BAM, or CRAM when the suffix is ``.cram``.
        space: Reference the alignments must be against.
        reference_fasta: Reference FASTA; required to open a CRAM.

    Raises:
        ValueError: If a CRAM has no ``reference_fasta``; the header lacks
            ``@HD SO:coordinate`` (fix: ``samtools sort``); there is no
            ``.bai``/``.csi``/``.crai`` index next to it (fix:
            ``samtools index``); or the ``@SQ`` names and lengths differ from
            ``space.sequences()`` as sets (fix: realign, or
            ``samtools reheader`` if only names differ).
    """
    bam = Path(bam)
    if bam.suffix.lower() == ".cram":
        if reference_fasta is None:
            raise ValueError(f"reading CRAM {bam} needs reference_fasta")
        fh = pysam.AlignmentFile(
            str(bam), "rc", reference_filename=str(reference_fasta)
        )
    else:
        fh = pysam.AlignmentFile(str(bam), "rb")
    with fh:
        header = fh.header.to_dict()

    order = header.get("HD", {}).get("SO")
    if order != "coordinate":
        raise ValueError(
            f"{bam} is not coordinate-sorted (@HD SO:{order}); "
            f"run: samtools sort -o sorted.bam {bam}"
        )
    if not _has_index(bam):
        raise ValueError(
            f"{bam} has no .bai/.csi/.crai index; run: samtools index {bam}"
        )

    got = {sq["SN"]: int(sq["LN"]) for sq in header.get("SQ", [])}
    want = dict(space.sequences())
    if got == want:
        return
    problems = []
    missing = sorted(set(want) - set(got))
    extra = sorted(set(got) - set(want))
    if missing:
        problems.append(f"missing from @SQ: {', '.join(missing)}")
    if extra:
        problems.append(f"extra in @SQ: {', '.join(extra)}")
    lengths = [
        f"{n} (@SQ {got[n]}, reference {want[n]})"
        for n in sorted(set(got) & set(want))
        if got[n] != want[n]
    ]
    if lengths:
        problems.append(f"length differs: {', '.join(lengths)}")
    raise ValueError(
        f"{bam} @SQ lines do not match reference {space.reference.id!r}: "
        + "; ".join(problems)
        + "; realign to this reference (minimap2 | samtools sort), or fix"
        " names with samtools reheader"
    )
