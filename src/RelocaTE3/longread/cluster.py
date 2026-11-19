"""Long-read clustering: Signatures into candidate insertion sites.

Rules, applied per ``locus.key`` (``(reference_id, sequence_id)``, keys in
sorted order; different keys never merge):

1. INS and SPLIT_INS, sorted by ``(start, signature_id)``, are joined by
   single linkage: a signature joins the current group when its start minus
   the previous start is ``<= window``.
2. Each group, sorted by ``(insert_length, signature_id)``, is split into
   subgroups where a length exceeds the immediately preceding one by more
   than ``length_tolerance`` (``l > prev * (1 + tol)``).
3. Each CLIP joins the subgroup whose position (lower median of member
   starts) is nearest, if within ``window``; ties go to the earlier
   position, then to the subgroup with the smaller first ``signature_id``.
   Unattached clips are single-linkage grouped into clip-only clusters.
4. A cluster's position is the lower median of its INS/SPLIT_INS starts (of
   its clip starts for clip-only clusters); its strand is ``+``.
5. Clusters are ordered by ``(reference_id, sequence_id, position, first
   signature_id)``; ``cluster_id`` is ``{reference_id}:{sequence_id}:{pos}``
   with ``.1``, ``.2`` ... on the second and later clusters sharing that id.

minimap2 left-aligns gaps, so INS and LEFT clips sit ``TSD length`` before
the true site while SPLIT_INS and RIGHT clips sit at it; the window absorbs
that offset. The output does not depend on input order.
"""

from __future__ import annotations

from bisect import bisect_left
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, replace

from RelocaTE3.longread.model import (
    ClipSide, Cluster, Locus, Signature, SignatureKind,
)  # fmt: skip


@dataclass(frozen=True)
class ClusterParams:
    """Clustering parameters.

    Attributes:
        window: Maximum distance (bp) for single linkage and clip attachment.
        length_tolerance: Relative length step that starts a new subgroup.
    """

    window: int
    length_tolerance: float = 0.30


def _by_pos(s: Signature) -> tuple[int, str]:
    return (s.locus.start, s.signature_id)


def _lower_median(values: list[int]) -> int:
    ordered = sorted(values)
    return ordered[(len(ordered) - 1) // 2]


def _single_linkage(sigs: list[Signature], window: int) -> list[list[Signature]]:
    groups: list[list[Signature]] = []
    for s in sorted(sigs, key=_by_pos):
        if groups and s.locus.start - groups[-1][-1].locus.start <= window:
            groups[-1].append(s)
        else:
            groups.append([s])
    return groups


def _split_by_length(group: list[Signature], tol: float) -> list[list[Signature]]:
    subgroups: list[list[Signature]] = []
    for s in sorted(group, key=lambda s: (s.insert_length, s.signature_id)):
        if subgroups and s.insert_length <= subgroups[-1][-1].insert_length * (1 + tol):
            subgroups[-1].append(s)
        else:
            subgroups.append([s])
    return subgroups


def _record(key: tuple[str, str], insertions: list[Signature],
            clips: list[Signature]) -> Cluster:  # fmt: skip
    members = insertions + clips
    pos = _lower_median([s.locus.start for s in (insertions or clips)])
    kinds = [s.kind for s in insertions]
    sides = [s.clip_side for s in clips]
    return Cluster(
        cluster_id="",
        locus=Locus(key[0], key[1], pos, pos, "+"),
        n_reads=len({s.read_id for s in members}),
        n_ins=kinds.count(SignatureKind.INS),
        n_split_ins=kinds.count(SignatureKind.SPLIT_INS),
        n_clip_left=sides.count(ClipSide.LEFT),
        n_clip_right=sides.count(ClipSide.RIGHT),
        insert_lengths=tuple(sorted(s.insert_length for s in insertions)),
        signature_ids=tuple(sorted(s.signature_id for s in members)),
    )


def _assign_clips(
    anchors: list[tuple[int, str, int]], clip_positions: list[int], window: int
) -> list[int | None]:
    """Subgroup index each clip joins, or None; O(log anchors) per clip.

    ``anchors`` are ``(position, first signature_id, subgroup index)`` sorted
    ascending, so for each position the first anchor wins same-position ties.
    The nearest position is next to the bisection point; on an equal
    distance the earlier (left) position wins.
    """
    positions = [a[0] for a in anchors]
    out: list[int | None] = []
    for x in clip_positions:
        right = bisect_left(positions, x)
        best: int | None = None
        best_dist = window + 1
        if right > 0:
            left = bisect_left(positions, positions[right - 1])
            best, best_dist = left, x - positions[left]
        if right < len(positions) and positions[right] - x < best_dist:
            best, best_dist = right, positions[right] - x
        out.append(anchors[best][2] if best is not None and best_dist <= window
                   else None)  # fmt: skip
    return out


def _cluster_key(key: tuple[str, str], sigs: list[Signature],
                 params: ClusterParams) -> list[Cluster]:  # fmt: skip
    window = params.window
    inserts = [s for s in sigs if s.kind is not SignatureKind.CLIP]
    subgroups = [
        sub
        for group in _single_linkage(inserts, window)
        for sub in _split_by_length(group, params.length_tolerance)
    ]
    # Sorted order encodes the tie rules: earlier position, then smaller
    # first signature_id.
    anchors = sorted(
        (_lower_median([s.locus.start for s in sub]),
         min(s.signature_id for s in sub), i)
        for i, sub in enumerate(subgroups)
    )  # fmt: skip
    clips = sorted((s for s in sigs if s.kind is SignatureKind.CLIP), key=_by_pos)
    attached: list[list[Signature]] = [[] for _ in subgroups]
    loose: list[Signature] = []
    targets = _assign_clips(anchors, [c.locus.start for c in clips], window)
    for clip, target in zip(clips, targets):
        if target is None:
            loose.append(clip)
        else:
            attached[target].append(clip)
    clusters = [_record(key, sub, attached[i]) for i, sub in enumerate(subgroups)]
    clusters += [_record(key, [], g) for g in _single_linkage(loose, window)]
    return clusters


def cluster_signatures(
    signatures: Iterable[Signature], params: ClusterParams
) -> list[Cluster]:
    """Cluster signatures into candidate insertion sites.

    Args:
        signatures: Signatures to cluster, consumed once, in any order.
        params: Window and length tolerance.

    Returns:
        Clusters ordered by reference, sequence, position and first
        signature id, with unique ``cluster_id``s.

    Raises:
        ValueError: If two signatures share a ``signature_id``.
    """
    by_key: dict[tuple[str, str], list[Signature]] = defaultdict(list)
    seen: set[str] = set()
    for s in signatures:
        if s.signature_id in seen:
            raise ValueError(f"duplicate signature_id {s.signature_id!r}")
        seen.add(s.signature_id)
        by_key[s.locus.key].append(s)

    clusters = [c for key in sorted(by_key)
                for c in _cluster_key(key, by_key[key], params)]  # fmt: skip
    clusters.sort(key=lambda c: (c.locus.reference_id, c.locus.sequence_id,
                                 c.locus.start, c.signature_ids[0]))  # fmt: skip

    out: list[Cluster] = []
    used: dict[str, int] = {}
    for c in clusters:
        base = f"{c.locus.reference_id}:{c.locus.sequence_id}:{c.locus.start}"
        n = used.get(base, 0)
        used[base] = n + 1
        out.append(replace(c, cluster_id=base if n == 0 else f"{base}.{n}"))
    return out
