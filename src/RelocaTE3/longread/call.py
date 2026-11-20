"""Long-read calls: TE summaries, then calls from clusters.

Signature summary (:func:`summarize_signature`): hits with ``identity <
min_identity`` are dropped. Aligned query bases are summed per
``(family, te_strand)``; the best pair has the most bases (ties: family
name, then strand, ``"+"`` before ``"-"``). Coverage is the union of all
kept hits over the query length. The signature is nested when a hit of
another family, at least ``NESTED_MIN_BP`` long, overlaps the best family's
hits by at most ``NESTED_MAX_OVERLAP`` of its length.

Read summary (:func:`read_family`): bases per ``(family, strand)`` summed
over the read's signatures; the top pair, ranked as above.

Calls (:func:`call_clusters`):

1. Inputs: clusters; for every signature id, a :class:`SignatureFacts` and a
   :class:`SignatureTe`. A cluster signature missing from either map raises
   ``ValueError`` naming it.
2. Cluster family: :func:`read_family` per read of the cluster; the cluster
   family is the family held by most reads (ties: more summed bases of that
   family over those reads, then name). ``NO_HIT`` if no read has one.
3. Merge: clusters are walked in ``(reference_id, sequence_id, start,
   cluster_id)`` order. A cluster of family F joins the open group of
   ``(key, F)`` when its start minus that group's last start is
   ``<= merge_window``; otherwise it opens a new group for ``(key, F)``.
   ``NO_HIT`` clusters are always their own group.
4. Vote, on the group's union of signatures: reads are the distinct read
   ids; each read's family is :func:`read_family` of its signatures in the
   group only; hit reads are reads with a family. ``te_read_fraction`` =
   hit reads / reads; the top family is chosen as in rule 2 over hit reads;
   ``family_fraction`` = its read count / hit reads; the strand is the
   majority strand among reads of the top family (tie: ``"+"``) and
   ``strand_fraction`` = that count / those reads. With no hit reads:
   family ``NO_HIT``, strand ``"."``, all three fractions 0.0.
5. Tier: ``SPANNING`` if the group has an INS or SPLIT_INS signature, else
   ``CLIP_ONLY``.
6. Position: lower median of insert-signature positions (``SPANNING``),
   else of all signature positions.
7. Medians: ``median_insert_length`` is the lower median of insert-signature
   lengths (``SPANNING``), else 0. ``median_te_coverage`` is the median
   coverage over insert signatures (``SPANNING``) or all signatures
   (``CLIP_ONLY``), rounded to 4 decimals. ``median_identity`` is the median
   identity over all group signatures whose family is the top family,
   rounded to 4 decimals, 0.0 if none. ``nested``: at least half of the hit
   reads have a nested signature (False with no hit reads).
8. Filters, each checked independently, in :class:`FilterReason` order:
   ``NO_TE_HIT`` (no hit reads); ``FEW_TE_READS`` (hit reads and
   ``te_read_fraction < min_te_reads_fraction``); ``FAMILY_DISCORDANT``
   (hit reads and ``family_fraction < min_family_fraction``);
   ``LOW_SUPPORT`` (``n_reads < min_support``); ``LOW_TE_FRACTION``
   (``SPANNING`` and ``median_te_coverage < min_te_coverage``).
9. Output is ordered by ``(reference_id, sequence_id, position, first
   cluster_id)``; ``call_id`` is ``{reference_id}:{sequence_id}:{position}``
   with ``.1``, ``.2`` ... on repeats, as for clusters.

Signature counts (``n_ins`` ... ``n_clip_right``) are over the group's
union of signatures. ``min_te_identity`` is applied when ``te`` is built
with :func:`summarize_signature`; :class:`CallParams` carries it so one
object describes a run.
"""

from __future__ import annotations

import statistics
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace

from RelocaTE3.longread.model import (
    NO_HIT, Call, CallTier, ClipSide, Cluster, FilterReason, Locus, Signature,
    SignatureKind, TeHit,
)  # fmt: skip

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


@dataclass(frozen=True)
class CallParams:
    """Calling parameters.

    Attributes:
        min_support: Fewest distinct reads for a PASS call.
        min_te_identity: Lowest TE-hit identity used (applied when
            summarizing signatures).
        merge_window: Largest start distance (bp) for merging same-family
            clusters.
        min_te_reads_fraction: Lowest fraction of reads with a TE hit.
        min_family_fraction: Lowest fraction of hit reads of the top family.
        min_te_coverage: Lowest median TE coverage of spanning calls.
    """

    min_support: int
    min_te_identity: float
    merge_window: int
    min_te_reads_fraction: float = 0.80
    min_family_fraction: float = 0.80
    min_te_coverage: float = 0.80


@dataclass(frozen=True)
class SignatureFacts:
    """The fields of a :class:`Signature` that calling needs.

    Attributes:
        read_id: Read the signature comes from.
        kind: Signature kind.
        position: Signature start on the reference.
        insert_length: Inserted (or clipped) length.
        clip_side: Clip side; ``ClipSide.NONE`` unless a clip.
    """

    read_id: str
    kind: SignatureKind
    position: int
    insert_length: int
    clip_side: ClipSide

    @classmethod
    def of(cls, s: Signature) -> SignatureFacts:
        """Build from a full signature.

        Args:
            s: The signature.

        Returns:
            Its calling facts.
        """
        return cls(s.read_id, s.kind, s.locus.start, s.insert_length, s.clip_side)


@dataclass(frozen=True)
class _Vote:
    n_reads: int
    hit_reads: dict[str, tuple[str, str]]
    family: str
    strand: str
    family_fraction: float
    strand_fraction: float
    te_read_fraction: float


def _lower_median(values: list[int]) -> int:
    ordered = sorted(values)
    return ordered[(len(ordered) - 1) // 2]


def _vote(sig_ids: Sequence[str], facts: Mapping[str, SignatureFacts],
          te: Mapping[str, SignatureTe]) -> _Vote:  # fmt: skip
    by_read: dict[str, list[SignatureTe]] = {}
    for sid in sig_ids:
        by_read.setdefault(facts[sid].read_id, []).append(te[sid])
    hit_reads: dict[str, tuple[str, str]] = {}
    family_bases: Counter[str] = Counter()
    for read, summaries in by_read.items():
        top = read_family(summaries)
        if top is None:
            continue
        hit_reads[read] = top
        family_bases[top[0]] += sum(n for s in summaries for (f, _), n in s.bases
                                    if f == top[0])  # fmt: skip
    if not hit_reads:
        return _Vote(len(by_read), {}, NO_HIT, ".", 0.0, 0.0, 0.0)
    counts = Counter(f for f, _ in hit_reads.values())
    family = min(counts, key=lambda f: (-counts[f], -family_bases[f], f))
    strands = Counter(s for f, s in hit_reads.values() if f == family)
    strand = "+" if strands["+"] >= strands["-"] else "-"
    return _Vote(
        len(by_read), hit_reads, family, strand,
        counts[family] / len(hit_reads), strands[strand] / counts[family],
        len(hit_reads) / len(by_read),
    )  # fmt: skip


def _group_clusters(
    clusters: list[Cluster], facts: Mapping[str, SignatureFacts],
    te: Mapping[str, SignatureTe], window: int,
) -> list[list[Cluster]]:  # fmt: skip
    groups: list[list[Cluster]] = []
    open_group: dict[tuple[str, str, str], list[Cluster]] = {}
    order = sorted(clusters, key=lambda c: (*c.locus.key, c.locus.start, c.cluster_id))
    for c in order:
        family = _vote(c.signature_ids, facts, te).family
        slot = (*c.locus.key, family)
        group = open_group.get(slot)
        if (family != NO_HIT and group is not None
                and c.locus.start - group[-1].locus.start <= window):  # fmt: skip
            group.append(c)
            continue
        group = [c]
        groups.append(group)
        if family != NO_HIT:
            open_group[slot] = group
    return groups


def _make_call(group: list[Cluster], facts: Mapping[str, SignatureFacts],
               te: Mapping[str, SignatureTe], params: CallParams) -> Call:  # fmt: skip
    sig_ids = list(dict.fromkeys(sid for c in group for sid in c.signature_ids))
    vote = _vote(sig_ids, facts, te)
    inserts = [sid for sid in sig_ids if facts[sid].kind is not SignatureKind.CLIP]
    spanning = bool(inserts)
    basis = inserts if spanning else sig_ids
    position = _lower_median([facts[sid].position for sid in basis])
    coverage = round(statistics.median(te[sid].coverage for sid in basis), 4)
    identities = [i for sid in sig_ids if te[sid].family == vote.family
                  for i in te[sid].identities]  # fmt: skip
    nested_reads = {facts[sid].read_id for sid in sig_ids if te[sid].nested}
    n_nested = sum(1 for r in vote.hit_reads if r in nested_reads)
    hits = len(vote.hit_reads)
    checks = {
        FilterReason.NO_TE_HIT: hits == 0,
        FilterReason.FEW_TE_READS: hits > 0
        and vote.te_read_fraction < params.min_te_reads_fraction,
        FilterReason.FAMILY_DISCORDANT: hits > 0
        and vote.family_fraction < params.min_family_fraction,
        FilterReason.LOW_SUPPORT: vote.n_reads < params.min_support,
        FilterReason.LOW_TE_FRACTION: spanning and coverage < params.min_te_coverage,
    }
    kinds = [facts[sid].kind for sid in sig_ids]
    sides = [facts[sid].clip_side for sid in sig_ids]
    key = group[0].locus.key
    return Call(
        call_id="",
        locus=Locus(key[0], key[1], position, position, "+"),
        tier=CallTier.SPANNING if spanning else CallTier.CLIP_ONLY,
        filters=tuple(r for r in FilterReason if checks[r]),
        family=vote.family,
        te_strand=vote.strand,
        family_fraction=vote.family_fraction,
        strand_fraction=vote.strand_fraction,
        te_read_fraction=vote.te_read_fraction,
        n_reads=vote.n_reads,
        n_ins=kinds.count(SignatureKind.INS),
        n_split_ins=kinds.count(SignatureKind.SPLIT_INS),
        n_clip_left=sides.count(ClipSide.LEFT),
        n_clip_right=sides.count(ClipSide.RIGHT),
        median_insert_length=(_lower_median([facts[s].insert_length for s in inserts])
                              if spanning else 0),
        median_te_coverage=coverage,
        median_identity=round(statistics.median(identities), 4) if identities else 0.0,
        nested=hits > 0 and 2 * n_nested >= hits,
        cluster_ids=tuple(c.cluster_id for c in group),
    )  # fmt: skip


def call_clusters(
    clusters: Iterable[Cluster],
    facts: Mapping[str, SignatureFacts],
    te: Mapping[str, SignatureTe],
    params: CallParams,
) -> list[Call]:
    """Merge same-family clusters into calls and filter them.

    Args:
        clusters: Clusters from ``cluster_signatures``, in any order.
        facts: Calling facts for every signature id named by a cluster.
        te: TE summary for every signature id named by a cluster.
        params: Merge window and filter thresholds.

    Returns:
        Calls ordered by reference, sequence, position and first cluster
        id, with unique ``call_id``s.

    Raises:
        ValueError: If a cluster names a signature missing from ``facts``
            or ``te``.
    """
    clusters = list(clusters)
    for c in clusters:
        for sid in c.signature_ids:
            for name, table in (("facts", facts), ("te", te)):
                if sid not in table:
                    raise ValueError(f"cluster {c.cluster_id!r}: signature "
                                     f"{sid!r} missing from {name}")  # fmt: skip
    calls = [
        _make_call(g, facts, te, params)
        for g in _group_clusters(clusters, facts, te, params.merge_window)
    ]
    calls.sort(key=lambda c: (*c.locus.key, c.locus.start, c.cluster_ids[0]))
    out: list[Call] = []
    used: dict[str, int] = {}
    for call in calls:
        locus = call.locus
        base = f"{locus.reference_id}:{locus.sequence_id}:{locus.start}"
        n = used.get(base, 0)
        used[base] = n + 1
        out.append(replace(call, call_id=base if n == 0 else f"{base}.{n}"))
    return out
