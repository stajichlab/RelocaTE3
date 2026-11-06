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
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from RelocaTE3 import logger
from RelocaTE3.aligners import (
    BLAT_CHUNK_SEQS,
    TEChunk,
    TESide,
    search_te_chunk,
    split_te_query,
    te_sam_header,
)

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
        ValueError: on a wrong header, a malformed line or a chunk file name
            given on more than one line.
    """
    path = Path(path)
    lines = path.read_text().splitlines()
    if not lines or tuple(lines[0].split("\t")) != MANIFEST_COLUMNS:
        raise ValueError(
            f"{path} is not a TE chunk manifest: expected the header "
            f"{chr(9).join(MANIFEST_COLUMNS)!r}"
        )
    rows = []
    seen: dict[str, int] = {}
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
            if chunk in seen:
                raise ValueError(
                    f"chunk {chunk} named more than once (first on line {seen[chunk]})"
                )
            seen[chunk] = number
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
        ValueError: on a manifest side without a read file, a chunk named
            twice, chunk indices that
            are not exactly ``0..k-1`` per side, a wrong mate, or chunk BAMs
            that are missing, not in the manifest, or given more than once.
    """
    read_files = [Path(f) for f in read_files]
    if not 1 <= len(read_files) <= 2:
        raise ValueError(f"expected one or two read files, got {len(read_files)}")
    paired = len(read_files) == 2
    directions = SIDES[: len(read_files)]

    chunk_names = Counter(row.chunk for row in rows)
    repeated_chunks = sorted(n for n, c in chunk_names.items() if c > 1)
    if repeated_chunks:
        raise ValueError(
            f"chunk(s) {', '.join(repeated_chunks)} named more than once in the manifest"
        )

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


# ---------------------------------------------------------------------------
# the three steps as the subcommands run them
# ---------------------------------------------------------------------------


def chunk_dir_for(outdir, name: str) -> Path:
    """``<outdir>/<name>.te_chunks``: where :func:`split_reads` writes."""
    return Path(outdir) / f"{name}.te_chunks"


def split_reads(read_files, name: str, outdir, chunk_reads: int = BLAT_CHUNK_SEQS) -> Path:
    """``te-split``: cut each read side into chunk FASTAs and write the manifest.

    Writes ``<outdir>/<name>.te_chunks/chunk.<name>.<side>.NNNNN.fa`` and
    ``manifest.tsv`` there; returns the manifest path.

    Raises:
        FileExistsError: if ``<name>.te_chunks`` already exists.
    """
    read_files = [Path(f) for f in read_files]
    if not 1 <= len(read_files) <= 2:
        raise ValueError(f"expected one or two read files, got {len(read_files)}")
    chunk_dir = chunk_dir_for(outdir, name)
    if chunk_dir.exists():
        raise FileExistsError(
            f"{chunk_dir} already exists; remove it or choose another --outdir"
        )
    chunk_dir.mkdir(parents=True)
    paired = len(read_files) == 2
    rows = []
    for number, (side, read_file) in enumerate(zip(SIDES, read_files)):
        chunks = split_te_query(read_file, chunk_dir, f"{name}.{side}", chunk_reads)
        rows.extend(
            ManifestRow(
                side, index, Path(fa).name, n_reads,
                str(number + 1) if paired else None, read_file.name,
            )
            for index, (fa, n_reads) in enumerate(chunks)
        )
        logger.info("%s: %d chunk(s) of up to %d reads", read_file, len(chunks), chunk_reads)
    return write_manifest(chunk_dir / MANIFEST_NAME, rows)


def search_chunks(te_library, manifest, chunk_fastas, outdir, te_opts=(), threads: int = 1):
    """``te-search-chunk``: BLAT each chunk FASTA into ``<outdir>/<chunk>.bam``.

    Each chunk's mate comes from the manifest (so read names are canonicalised
    exactly as in ``relocaTE3 run``); ``threads`` chunks run at once.

    Returns:
        The chunk BAM paths, in the order of ``chunk_fastas``.

    Raises:
        ValueError: if a chunk is not in the manifest, is given twice, or
            holds a different number of reads than the manifest says.
    """
    rows = {r.chunk: r for r in read_manifest(manifest)}
    names = [Path(fa).name for fa in chunk_fastas]
    unknown = [n for n in names if n not in rows]
    if unknown:
        raise ValueError(f"chunk(s) not in the manifest {manifest}: {', '.join(unknown)}")
    repeated = sorted(n for n, c in Counter(names).items() if c > 1)
    if repeated:
        raise ValueError(f"chunk(s) given more than once: {', '.join(repeated)}")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sam_header = te_sam_header(te_library)  # once, before the threads start

    def one(fa) -> TEChunk:
        row = rows[Path(fa).name]
        found = search_te_chunk(
            te_library, fa, outdir / chunk_bam_name(fa), mate=row.mate,
            te_opts=te_opts, sam_header=sam_header,
        )
        if found.n_reads != row.n_reads:
            raise ValueError(
                f"{fa} holds {found.n_reads} reads; the manifest says {row.n_reads}"
            )
        return found

    if not chunk_fastas:
        return []
    with ThreadPoolExecutor(max_workers=max(1, min(int(threads or 1), len(chunk_fastas)))) as pool:
        return [chunk.bam for chunk in pool.map(one, chunk_fastas)]
