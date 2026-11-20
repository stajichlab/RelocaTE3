"""TE annotation of long-read signatures (milestone L3, step 4).

Every signature's sequence (inserted bases; clipped bases for CLIP) is
aligned to the TE library with one minimap2 run per sample. Flags were
measured on the small Chr3 panel (plans/2026-11-19-long-read-l3-plan.md):
the HiFi preset misses 145 bp TEs entirely; ``map-ont`` with k=11, w=5 hits
every family at 1-10 % error and nothing in random sequence; ``-r 100,100``
splits nested inserts. The same flags serve every read type, because the
query is a short insert, not a read.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from collections.abc import Iterable, Mapping
from pathlib import Path

from RelocaTE3.longread.model import NO_HIT, Signature, TeHit

_STDERR_TAIL = 20

TE_MINIMAP2_ARGS = ("-c", "-x", "map-ont", "-k", "11", "-w", "5", "-r",
                    "100,100", "--secondary=no")  # fmt: skip


def te_minimap2_command(library: Path, queries: Path, threads: int) -> list[str]:
    """Build the minimap2 command that aligns signature sequences to TEs.

    Args:
        library: TE library FASTA.
        queries: FASTA of signature sequences, named by signature id.
        threads: minimap2 threads (``-t``); at least 1.

    Returns:
        The argument list; PAF with CIGAR (``-c``) goes to stdout.

    Raises:
        ValueError: If ``threads`` is below 1.
    """
    if threads < 1:
        raise ValueError(f"threads must be at least 1, got {threads}")
    return ["minimap2", *TE_MINIMAP2_ARGS, "-t", str(threads), str(library),
            str(queries)]  # fmt: skip


def parse_paf_line(line: str, read_of: Mapping[str, str]) -> TeHit | None:
    """Turn one PAF line into a :class:`TeHit`.

    Args:
        line: A PAF line; the query name is a signature id.
        read_of: Read id of each signature id.

    Returns:
        The hit, or None for a secondary alignment (``tp:A:S``).

    Raises:
        KeyError: If the query name is not in ``read_of``.
    """
    f = line.rstrip("\n").split("\t")
    if any(tag == "tp:A:S" for tag in f[12:]):
        return None
    nmatch, block = int(f[9]), int(f[10])
    return TeHit(
        signature_id=f[0], read_id=read_of[f[0]], family=f[5], te_strand=f[4],
        te_start=int(f[7]), te_end=int(f[8]), te_length=int(f[6]),
        query_start=int(f[2]), query_end=int(f[3]), query_length=int(f[1]),
        identity=round(nmatch / block, 4) if block else 0.0,
    )  # fmt: skip


def _no_hit(sig_id: str, read_id: str, length: int) -> TeHit:
    return TeHit(sig_id, read_id, NO_HIT, ".", 0, 0, 0, 0, 0, length, 0.0)


def _tail(path: Path) -> str:
    lines = path.read_text(errors="replace").splitlines()
    return "\n".join(lines[-_STDERR_TAIL:])


def annotate_signatures(
    signatures: Iterable[Signature],
    library: Path,
    *,
    threads: int,
    workdir: Path,
    stats: dict[str, int] | None = None,
) -> list[TeHit]:
    """Align every signature's sequence to the TE library.

    Signatures with an empty sequence are not sent to minimap2. Temporary
    files live in a hidden directory under ``workdir`` that is removed on
    success and on failure.

    Args:
        signatures: Signatures with unique ``signature_id``.
        library: TE library FASTA.
        threads: minimap2 threads; at least 1.
        workdir: Existing directory for the temporary files.
        stats: If given, filled with ``signatures``, ``queries`` (sequences
            sent to minimap2), ``with_hit``, ``without_hit`` and ``hits``.

    Returns:
        Per signature, in input order, its hits sorted by
        ``(query_start, query_end, family)``, or one ``NO_HIT`` row.

    Raises:
        FileNotFoundError: If ``library`` is not a file.
        ValueError: If a signature id repeats, or ``threads`` is below 1.
        RuntimeError: If minimap2 is not on PATH (checked only when there
            is a query) or exits non-zero (the message names the command and
            holds the last lines of its stderr).
    """
    library = Path(library)
    if not library.is_file():
        raise FileNotFoundError(f"TE library not found: {library}")
    if threads < 1:
        raise ValueError(f"threads must be at least 1, got {threads}")

    order: list[tuple[str, str, int]] = []
    read_of: dict[str, str] = {}
    queries: list[tuple[str, str]] = []
    for sig in signatures:
        sid = sig.signature_id
        if sid in read_of:
            raise ValueError(f"duplicate signature_id {sid!r}")
        read_of[sid] = sig.read_id
        order.append((sid, sig.read_id, len(sig.sequence)))
        if sig.sequence:
            queries.append((sid, sig.sequence))

    by_sig: dict[str, list[TeHit]] = {}
    if queries:
        if shutil.which("minimap2") is None:
            raise RuntimeError("minimap2 not found on PATH")
        with tempfile.TemporaryDirectory(dir=workdir, prefix=".lr-annotate.") as tmp:
            tmpdir = Path(tmp)
            fasta = tmpdir / "queries.fa"
            with fasta.open("w") as fh:
                for sid, seq in queries:
                    fh.write(f">{sid}\n{seq}\n")
            paf = tmpdir / "hits.paf"
            err = tmpdir / "minimap2.err"
            cmd = te_minimap2_command(library, fasta, threads)
            with paf.open("wb") as out, err.open("wb") as errfh:
                rc = subprocess.run(cmd, stdout=out, stderr=errfh).returncode
            if rc != 0:
                raise RuntimeError(
                    f"minimap2 failed (exit {rc}): {' '.join(cmd)}\n{_tail(err)}"
                )
            with paf.open() as fh:
                for line in fh:
                    if not line.strip():
                        continue
                    hit = parse_paf_line(line, read_of)
                    if hit is not None:
                        by_sig.setdefault(hit.signature_id, []).append(hit)

    result: list[TeHit] = []
    with_hit = 0
    n_hits = 0
    for sid, read_id, length in order:
        group = by_sig.get(sid)
        if group:
            group.sort(key=lambda h: (h.query_start, h.query_end, h.family))
            result.extend(group)
            with_hit += 1
            n_hits += len(group)
        else:
            result.append(_no_hit(sid, read_id, length))
    if stats is not None:
        stats.update(signatures=len(order), queries=len(queries),
                     with_hit=with_hit, without_hit=len(order) - with_hit,
                     hits=n_hits)  # fmt: skip
    return result
