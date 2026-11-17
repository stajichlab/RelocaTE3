"""Versioned TSV.gz tables that join the long-read stages.

Each file is gzip text. Line 1 is the schema line, tab-separated:
``#relocate3:<schema>``, ``v<version>``, then ``key=value`` metadata
(``read_type``, ``reference_id``). Line 2 is the column header; rows follow.
A stage refuses a file of the wrong schema or version, or with a changed
header, so stages cannot silently read each other's output.

Writes go to ``<path>.partial`` and are renamed on success; on any error the
partial file is removed. An existing ``<path>`` is never overwritten.
Fields are joined with tabs and must not contain tabs or newlines.
"""

from __future__ import annotations

import gzip
import os
from collections.abc import Callable, Iterable, Iterator, Sequence
from pathlib import Path
from typing import TextIO, TypeVar

from RelocaTE3.longread.model import (
    ClipSide, Cluster, Locus, ReadType, Signature, SignatureKind,
)  # fmt: skip

T = TypeVar("T")

SIGNATURE_SCHEMA = "lr-signatures"
CLUSTER_SCHEMA = "lr-clusters"
VERSION = 1
_MAGIC = "#relocate3:"
_META_KEYS = ("read_type", "reference_id")

SIGNATURE_COLUMNS = (
    "signature_id", "kind", "reference_id", "sequence_id", "start", "end",
    "strand", "insert_length", "clip_side", "read_id", "query_start",
    "query_end", "mapq", "sequence", "left_flank", "right_flank",
)  # fmt: skip
CLUSTER_COLUMNS = (
    "cluster_id", "reference_id", "sequence_id", "start", "end", "strand",
    "n_reads", "n_ins", "n_split_ins", "n_clip_left", "n_clip_right",
    "insert_lengths", "signature_ids",
)  # fmt: skip


def _join_line(fields: Sequence[object]) -> str:
    texts = [str(f) for f in fields]
    for text in texts:
        if "\t" in text or "\n" in text or "\r" in text:
            raise ValueError(f"field contains a tab or newline: {text!r}")
    return "\t".join(texts) + "\n"


def _write_table(
    path: str | os.PathLike[str],
    schema: str,
    columns: tuple[str, ...],
    items: Iterable[T],
    to_fields: Callable[[T], Sequence[object]],
    read_type: ReadType,
    reference_id: str,
) -> int:
    out = Path(path)
    if out.exists():
        raise FileExistsError(f"output exists, not overwriting: {out}")
    partial = Path(f"{out}.partial")
    n = 0
    try:
        with gzip.open(partial, "wt", encoding="utf-8", newline="") as fh:
            fh.write(_join_line((
                f"{_MAGIC}{schema}", f"v{VERSION}",
                f"read_type={read_type}", f"reference_id={reference_id}",
            )))  # fmt: skip
            fh.write(_join_line(columns))
            for item in items:
                fh.write(_join_line(to_fields(item)))
                n += 1
        if out.exists():
            raise FileExistsError(f"output appeared while writing: {out}")
        os.replace(partial, out)
    except BaseException:
        partial.unlink(missing_ok=True)
        raise
    return n


def _read_table(
    path: str | os.PathLike[str],
    schema: str,
    columns: tuple[str, ...],
    from_fields: Callable[[list[str]], T],
) -> tuple[dict[str, str], Iterator[T]]:
    fh = gzip.open(path, "rt", encoding="utf-8")
    try:
        meta = _parse_schema_line(fh.readline(), schema, path)
        header = tuple(fh.readline().rstrip("\n").split("\t"))
        if header != columns:
            raise ValueError(
                f"{path}: header differs from {schema} v{VERSION}: "
                f"expected {columns}, got {header}"
            )
    except BaseException:
        fh.close()
        raise
    return meta, _iter_rows(fh, path, len(columns), from_fields)


def _parse_schema_line(line: str, schema: str, path: object) -> dict[str, str]:
    parts = line.rstrip("\n").split("\t")
    if not parts[0].startswith(_MAGIC):
        raise ValueError(f"{path}: not a relocate3 table (no schema line)")
    found = parts[0][len(_MAGIC) :]
    if found != schema:
        raise ValueError(f"{path}: wrong schema {found!r}, expected {schema!r}")
    version = parts[1] if len(parts) > 1 else ""
    if version != f"v{VERSION}":
        raise ValueError(
            f"{path}: unknown {schema} version {version!r}, expected 'v{VERSION}'"
        )
    meta = {"schema": found, "version": str(VERSION)}
    for item in parts[2:]:
        key, sep, value = item.partition("=")
        if not sep:
            raise ValueError(f"{path}: bad schema metadata {item!r}")
        meta[key] = value
    missing = [k for k in _META_KEYS if k not in meta]
    if missing:
        raise ValueError(f"{path}: schema line lacks metadata {missing}")
    return meta


def _iter_rows(
    fh: TextIO,
    path: object,
    n_columns: int,
    from_fields: Callable[[list[str]], T],
) -> Iterator[T]:
    with fh:
        for lineno, line in enumerate(fh, start=3):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != n_columns:
                raise ValueError(
                    f"{path}:{lineno}: expected {n_columns} fields, got {len(fields)}"
                )
            yield from_fields(fields)


def _split_list(text: str) -> list[str]:
    return text.split(",") if text else []


def _signature_fields(s: Signature) -> tuple[object, ...]:
    loc = s.locus
    return (
        s.signature_id, s.kind, loc.reference_id, loc.sequence_id, loc.start,
        loc.end, loc.strand, s.insert_length, s.clip_side, s.read_id,
        s.query_start, s.query_end, s.mapq, s.sequence, s.left_flank,
        s.right_flank,
    )  # fmt: skip


def _signature_from(f: list[str]) -> Signature:
    return Signature(
        signature_id=f[0],
        kind=SignatureKind(f[1]),
        locus=Locus(f[2], f[3], int(f[4]), int(f[5]), f[6]),
        insert_length=int(f[7]),
        clip_side=ClipSide(f[8]),
        read_id=f[9],
        query_start=int(f[10]),
        query_end=int(f[11]),
        mapq=int(f[12]),
        sequence=f[13],
        left_flank=f[14],
        right_flank=f[15],
    )


def _cluster_fields(c: Cluster) -> tuple[object, ...]:
    loc = c.locus
    return (
        c.cluster_id, loc.reference_id, loc.sequence_id, loc.start, loc.end,
        loc.strand, c.n_reads, c.n_ins, c.n_split_ins, c.n_clip_left,
        c.n_clip_right, ",".join(str(n) for n in c.insert_lengths),
        ",".join(c.signature_ids),
    )  # fmt: skip


def _cluster_from(f: list[str]) -> Cluster:
    return Cluster(
        cluster_id=f[0],
        locus=Locus(f[1], f[2], int(f[3]), int(f[4]), f[5]),
        n_reads=int(f[6]),
        n_ins=int(f[7]),
        n_split_ins=int(f[8]),
        n_clip_left=int(f[9]),
        n_clip_right=int(f[10]),
        insert_lengths=tuple(int(n) for n in _split_list(f[11])),
        signature_ids=tuple(_split_list(f[12])),
    )


def write_signatures(
    path: str | os.PathLike[str],
    signatures: Iterable[Signature],
    *,
    read_type: ReadType,
    reference_id: str,
) -> int:
    """Write signatures to a new ``lr-signatures`` v1 TSV.gz file.

    Args:
        path: Output path; must not exist.
        signatures: Signatures to write, consumed once.
        read_type: Read technology, recorded in the schema line.
        reference_id: Reference identifier, recorded in the schema line.

    Returns:
        Number of rows written.

    Raises:
        FileExistsError: If ``path`` already exists.
        ValueError: If a field contains a tab or newline.
    """
    return _write_table(path, SIGNATURE_SCHEMA, SIGNATURE_COLUMNS, signatures,
                        _signature_fields, read_type, reference_id)  # fmt: skip


def read_signatures(
    path: str | os.PathLike[str],
) -> tuple[dict[str, str], Iterator[Signature]]:
    """Read an ``lr-signatures`` v1 TSV.gz file.

    The schema line and header are checked before returning; rows are
    streamed by the returned iterator.

    Args:
        path: Input path.

    Returns:
        The metadata (``schema``, ``version``, ``read_type``,
        ``reference_id``) and an iterator over the signatures.

    Raises:
        ValueError: On a wrong schema, unknown version or changed header;
            the iterator raises it on a row with the wrong field count.
    """
    return _read_table(path, SIGNATURE_SCHEMA, SIGNATURE_COLUMNS,
                       _signature_from)  # fmt: skip


def write_clusters(
    path: str | os.PathLike[str],
    clusters: Iterable[Cluster],
    *,
    read_type: ReadType,
    reference_id: str,
) -> int:
    """Write clusters to a new ``lr-clusters`` v1 TSV.gz file.

    Args:
        path: Output path; must not exist.
        clusters: Clusters to write, consumed once.
        read_type: Read technology, recorded in the schema line.
        reference_id: Reference identifier, recorded in the schema line.

    Returns:
        Number of rows written.

    Raises:
        FileExistsError: If ``path`` already exists.
        ValueError: If a field contains a tab or newline.
    """
    return _write_table(path, CLUSTER_SCHEMA, CLUSTER_COLUMNS, clusters,
                        _cluster_fields, read_type, reference_id)  # fmt: skip


def read_clusters(
    path: str | os.PathLike[str],
) -> tuple[dict[str, str], Iterator[Cluster]]:
    """Read an ``lr-clusters`` v1 TSV.gz file.

    The schema line and header are checked before returning; rows are
    streamed by the returned iterator.

    Args:
        path: Input path.

    Returns:
        The metadata (``schema``, ``version``, ``read_type``,
        ``reference_id``) and an iterator over the clusters.

    Raises:
        ValueError: On a wrong schema, unknown version or changed header;
            the iterator raises it on a row with the wrong field count.
    """
    return _read_table(path, CLUSTER_SCHEMA, CLUSTER_COLUMNS, _cluster_from)
