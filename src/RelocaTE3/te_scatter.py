"""TE-search chunk manifest: how the BLAT search is spread over jobs.

The BLAT TE search runs as three steps (``relocaTE3 run`` runs them in one
process; the ``te-split`` / ``te-search-chunk`` / ``te-merge-trim``
subcommands run them in separate ones):

1. :func:`RelocaTE3.aligners.split_te_query` per read side -> chunk FASTAs;
2. :func:`RelocaTE3.aligners.search_te_chunk` per chunk -> chunk BAM;
3. :func:`RelocaTE3.librelocate.merge_trim` -> the files ``relocaTE3 run``
   writes.

The split step records its chunks in a manifest so the merge step can put the
chunk BAMs back in order and refuse an incomplete or inconsistent set.

Manifest format (TSV, one header line, one row per chunk, sides in
``left``, ``right`` order, chunks in index order)::

    side   index  chunk                    n_reads  mate  read_file
    left   0      chunk.S.left.00000.fa    200000   1     S_R1.fq.gz

- ``side``: ``left`` (R1 / single-end) or ``right`` (R2);
- ``index``: 0-based chunk position in the side's read file;
- ``chunk``: chunk FASTA file name (no directory); the chunk BAM's name is
  :func:`chunk_bam_name` of it, e.g. ``chunk.S.left.00000.bam``;
- ``n_reads``: records in the chunk; the side's read file must hold exactly
  Σ ``n_reads`` records;
- ``mate``: ``1`` / ``2`` for paired sides (read names get ``/1``, ``/2``),
  ``.`` for single-end;
- ``read_file``: base name of the read file the chunks came from (recorded for
  provenance; not enforced, since staged inputs may be renamed).

A side with no rows had no reads (empty read file).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from RelocaTE3 import logger
from RelocaTE3.aligners import TEChunk, TESide

MANIFEST_NAME = "manifest.tsv"
MANIFEST_COLUMNS = ("side", "index", "chunk", "n_reads", "mate", "read_file")
SIDES = ("left", "right")
_NO_MATE = "."


@dataclass(frozen=True)
class ManifestRow:
    """One chunk of one read side (a manifest line)."""

    side: str
    index: int
    chunk: str
    n_reads: int
    mate: str | None
    read_file: str


def chunk_bam_name(chunk_fasta) -> str:
    """The chunk BAM's file name for a chunk FASTA: ``.fa`` -> ``.bam``."""
    return Path(chunk_fasta).with_suffix(".bam").name


def write_manifest(path, rows) -> Path:
    """Write ``rows`` (:class:`ManifestRow`) as a manifest TSV at ``path``."""
    path = Path(path)
    lines = ["\t".join(MANIFEST_COLUMNS)]
    for r in rows:
        mate = r.mate if r.mate else _NO_MATE
        lines.append(f"{r.side}\t{r.index}\t{r.chunk}\t{r.n_reads}\t{mate}\t{r.read_file}")
    path.write_text("\n".join(lines) + "\n")
    return path


def read_manifest(path) -> list[ManifestRow]:
    """Parse a manifest TSV written by :func:`write_manifest`.

    Raises:
        ValueError: on a wrong header or a malformed line.
    """
    path = Path(path)
    lines = path.read_text().splitlines()
    if not lines or tuple(lines[0].split("\t")) != MANIFEST_COLUMNS:
        raise ValueError(
            f"{path} is not a TE chunk manifest: expected the header "
            f"{chr(9).join(MANIFEST_COLUMNS)!r}"
        )
    rows = []
    for number, line in enumerate(lines[1:], start=2):
        if not line:
            continue
        fields = line.split("\t")
        try:
            if len(fields) != len(MANIFEST_COLUMNS):
                raise ValueError(f"{len(fields)} columns")
            side, index, chunk, n_reads, mate, read_file = fields
            if side not in SIDES:
                raise ValueError(f"side {side!r} is not left or right")
            if mate not in ("1", "2", _NO_MATE):
                raise ValueError(f"mate {mate!r} is not 1, 2 or {_NO_MATE}")
            row = ManifestRow(
                side, int(index), chunk, int(n_reads),
                None if mate == _NO_MATE else mate, read_file,
            )
            if row.index < 0 or row.n_reads < 0 or not chunk:
                raise ValueError("negative index/n_reads or empty chunk name")
        except ValueError as err:
            raise ValueError(f"{path} line {number}: {err}") from None
        rows.append(row)
    return rows


def sides_from_manifest(rows, chunk_bams, read_files) -> list[TESide]:
    """Match chunk BAMs to the manifest and return the sides for ``merge_trim``.

    Args:
        rows: the manifest (:func:`read_manifest`).
        chunk_bams: every chunk BAM, any order; matched to rows by file name
            (:func:`chunk_bam_name`).
        read_files: the sample's read files, ``[R1]`` or ``[R1, R2]``.

    Returns:
        :class:`TESide` per read file, chunks in index order.

    Raises:
        ValueError: on a manifest side without a read file, chunk indices that
            are not exactly ``0..k-1`` per side, a wrong mate, or chunk BAMs
            that are missing, not in the manifest, or given more than once.
    """
    read_files = [Path(f) for f in read_files]
    if not 1 <= len(read_files) <= 2:
        raise ValueError(f"expected one or two read files, got {len(read_files)}")
    paired = len(read_files) == 2
    directions = SIDES[: len(read_files)]

    by_side: dict[str, list[ManifestRow]] = {d: [] for d in directions}
    for row in rows:
        if row.side not in by_side:
            raise ValueError(
                f"the manifest has chunks for side {row.side!r} but only one read "
                "file was given (single-end); give -r/--right for paired reads"
            )
        by_side[row.side].append(row)

    names = Counter(Path(b).name for b in chunk_bams)
    repeated = sorted(n for n, c in names.items() if c > 1)
    if repeated:
        raise ValueError(f"chunk BAM(s) given more than once: {', '.join(repeated)}")
    bam_paths = {Path(b).name: Path(b) for b in chunk_bams}

    expected: set[str] = set()
    sides = []
    for number, direction in enumerate(directions):
        side_rows = sorted(by_side[direction], key=lambda r: r.index)
        indices = [r.index for r in side_rows]
        if indices != list(range(len(side_rows))):
            raise ValueError(
                f"chunk indices for side {direction!r} must be "
                f"0..{len(side_rows) - 1} with no gaps or repeats; got {indices}"
            )
        want_mate = str(number + 1) if paired else None
        for r in side_rows:
            if r.mate != want_mate:
                raise ValueError(
                    f"mate {r.mate or _NO_MATE!r} for side {direction!r} "
                    f"(chunk {r.chunk}); expected {want_mate or _NO_MATE!r}"
                )
            if r.read_file != read_files[number].name:
                logger.warning(
                    "manifest side %s came from %s; merging with %s",
                    direction, r.read_file, read_files[number],
                )
        missing = [chunk_bam_name(r.chunk) for r in side_rows
                   if chunk_bam_name(r.chunk) not in bam_paths]
        if missing:
            raise ValueError(
                f"missing chunk BAM(s) for side {direction!r}: {', '.join(missing)} "
                f"({len(missing)} of {len(side_rows)}); every chunk in the "
                "manifest must be searched (te-search-chunk) before te-merge-trim"
            )
        expected.update(chunk_bam_name(r.chunk) for r in side_rows)
        sides.append(
            TESide(
                direction=direction,
                read_file=read_files[number],
                chunks=tuple(
                    TEChunk(bam=bam_paths[chunk_bam_name(r.chunk)], n_reads=r.n_reads)
                    for r in side_rows
                ),
            )
        )
    extra = sorted(set(bam_paths) - expected)
    if extra:
        raise ValueError(f"chunk BAM(s) not in the manifest: {', '.join(extra)}")
    return sides
