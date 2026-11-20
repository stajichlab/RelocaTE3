"""Long-read calls: TE summaries of signatures and reads.

Signature summary (:func:`summarize_signature`): hits with ``identity <
min_identity`` are dropped. Aligned query bases are summed per
``(family, te_strand)``; the best pair has the most bases (ties: family
name, then strand, ``"+"`` before ``"-"``). Coverage is the union of all
kept hits over the query length. The signature is nested when a hit of
another family, at least ``NESTED_MIN_BP`` long, overlaps the best family's
hits by at most ``NESTED_MAX_OVERLAP`` of its length.

Read summary (:func:`read_family`): bases per ``(family, strand)`` summed
over the read's signatures; the top pair, ranked as above.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from RelocaTE3.longread.model import NO_HIT, TeHit

NESTED_MIN_BP = 50
NESTED_MAX_OVERLAP = 0.10


@dataclass(frozen=True)
class SignatureTe:
    """TE summary of one signature after the identity filter.

    ``bases`` maps (family, strand) to aligned query bases; ``identities``
    holds the identity of each kept hit of the best family.

    Attributes:
        family: Best family, or ``NO_HIT`` when no hit was kept.
        te_strand: Best family's strand, or ``"."`` when no hit was kept.
        bases: ``((family, strand), bases)`` pairs, best first.
        coverage: Fraction of the query covered by the union of kept hits.
        nested: True when another family sits outside the best family's hits.
        identities: Identity of each kept hit of the best family.
    """

    family: str  # NO_HIT when no kept hit
    te_strand: str  # "." when no kept hit
    bases: tuple[tuple[tuple[str, str], int], ...]
    coverage: float
    nested: bool
    identities: tuple[float, ...]


def _rank(item: tuple[tuple[str, str], int]) -> tuple[int, str, str]:
    (family, strand), n = item
    return (-n, family, strand)


def _union(spans: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[list[int]] = []
    for a, b in sorted(spans):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def _overlap(span: tuple[int, int], union: list[tuple[int, int]]) -> int:
    a, b = span
    return sum(max(0, min(b, y) - max(a, x)) for x, y in union)


def summarize_signature(hits: Sequence[TeHit], min_identity: float) -> SignatureTe:
    """Summarize one signature's TE hits.

    Args:
        hits: All hits of one signature (or its single ``NO_HIT`` row).
        min_identity: Hits with a lower identity are ignored.

    Returns:
        The signature's best family, strand, coverage and nesting.
    """
    kept = [h for h in hits if h.family != NO_HIT and h.identity >= min_identity]
    if not kept:
        return SignatureTe(NO_HIT, ".", (), 0.0, False, ())
    bases: Counter[tuple[str, str]] = Counter()
    for h in kept:
        bases[(h.family, h.te_strand)] += h.aligned
    ranked = sorted(bases.items(), key=_rank)
    family, strand = ranked[0][0]
    length = kept[0].query_length
    covered = sum(b - a for a, b in _union((h.query_start, h.query_end) for h in kept))
    best_union = _union(
        (h.query_start, h.query_end) for h in kept if h.family == family
    )
    nested = any(
        h.family != family and h.aligned >= NESTED_MIN_BP
        and _overlap((h.query_start, h.query_end), best_union)
        <= NESTED_MAX_OVERLAP * h.aligned
        for h in kept
    )  # fmt: skip
    return SignatureTe(
        family, strand, tuple(ranked), covered / length if length else 0.0, nested,
        tuple(h.identity for h in kept if h.family == family),
    )  # fmt: skip


def read_family(summaries: Iterable[SignatureTe]) -> tuple[str, str] | None:
    """Top ``(family, strand)`` of one read over its signature summaries.

    Args:
        summaries: Summaries of the read's signatures.

    Returns:
        The pair with the most summed bases (ties: family name, then
        strand), or None when no signature has a kept hit.
    """
    total: Counter[tuple[str, str]] = Counter()
    for s in summaries:
        for key, n in s.bases:
            total[key] += n
    if not total:
        return None
    return sorted(total.items(), key=_rank)[0][0]
