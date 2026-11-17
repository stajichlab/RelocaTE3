"""Per-stage timing table for the long-read stages.

Each stage (lr-align, lr-signatures, lr-cluster) wraps its work in
:func:`stage_timer`, which appends one row to ``<outdir>/timing.tsv`` when the
block finishes, whether it succeeds or raises. The header is written only
when the file is new or empty. Each row goes out in a single ``write()`` on a
file opened in append mode, so rows from separate stages do not interleave.
"""

from __future__ import annotations

import resource
import sys
import time
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import RelocaTE3

TIMING_FILE = "timing.tsv"
TIMING_COLUMNS = (
    "stage", "sample", "start_utc", "wall_s", "max_rss_mb", "records_in",
    "records_out", "relocate3_version", "status",
)  # fmt: skip


@dataclass
class StageCounts:
    """Record counts a stage fills in inside a :func:`stage_timer` block.

    Attributes:
        records_in: Records the stage read.
        records_out: Records the stage wrote.
    """

    records_in: int = 0
    records_out: int = 0


def _max_rss_mb() -> float:
    """Return the peak RSS of this process or its children, in MiB.

    ``ru_maxrss`` is KiB on Linux and bytes on macOS. The larger of the
    self and children peaks is reported.
    """
    peak = max(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
    )
    divisor = 1024 * 1024 if sys.platform == "darwin" else 1024
    return peak / divisor


@contextmanager
def stage_timer(outdir: Path, stage: str, sample: str) -> Iterator[StageCounts]:
    """Append one row to ``<outdir>/timing.tsv`` when the block finishes.

    Columns: stage, sample, start_utc, wall_s, max_rss_mb (self and
    children, from resource.getrusage, KiB on Linux), records_in,
    records_out, relocate3_version, status ("ok" or the exception class).
    The header is written only when the file is new. On an exception the
    row is still written and the exception re-raised.

    Args:
        outdir: Directory holding ``timing.tsv``; created if missing.
        stage: Stage name, e.g. ``lr-align``.
        sample: Sample name.

    Yields:
        A :class:`StageCounts` the block fills in.
    """
    counts = StageCounts()
    start_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    t0 = time.perf_counter()
    status = "ok"
    try:
        yield counts
    except BaseException as exc:
        status = type(exc).__name__
        raise
    finally:
        wall_s = time.perf_counter() - t0
        fields = (
            stage, sample, start_utc, f"{wall_s:.3f}", f"{_max_rss_mb():.1f}",
            str(counts.records_in), str(counts.records_out),
            RelocaTE3.__version__, status,
        )  # fmt: skip
        row = "\t".join(fields) + "\n"
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        with (outdir / TIMING_FILE).open("a") as fh:
            if fh.tell() == 0:
                row = "\t".join(TIMING_COLUMNS) + "\n" + row
            fh.write(row)
