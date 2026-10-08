"""Resolve junction-family evidence without treating arbitrary hit ties as votes.

The interface accepts (selected family, equally best families, compatible
geometry) per read. Empty alternatives mean the legacy single-family evidence.
Resolution is deterministic and uses junction evidence only, never truth or
supporting mates. ``support`` counts reads distinguishing one family; confidence
is its fraction of all informative reads, including ambiguous reads.
"""

from collections import Counter
from dataclasses import dataclass, field


@dataclass(frozen=True)
class TEFamilyEvidence:
    """A single primary family together with independently supported evidence."""

    primary: str
    support: dict[str, int]
    confidence: float
    status: str
    ambiguous_reads: int = 0
    resolution: str = "selected_votes"
    candidate_support: dict[str, int] = field(default_factory=dict)


def resolve_family(reads, *, fallback_primary: str | None = None) -> TEFamilyEvidence:
    """Use independent consensus only when sufficiently supported and compatible.

    At least two distinguishing reads, a strict majority of distinguishing
    reads, and compatibility with a strict majority of all informative reads
    are required to resolve tied-family evidence. Different trim geometries
    prevent resolution. Otherwise preserve deterministic selected-vote fallback
    and report ambiguity; missing/NA families do not vote. Callers deduplicating
    evidence may supply the pre-deduplication primary for unresolved ties.
    That label is a compatibility fallback, never an independent evidence vote.
    """
    selected, unique, candidates = Counter(), Counter(), Counter()
    ambiguous = total = 0
    geometry_ok = True
    for primary, alternatives, compatible in reads:
        if not primary or primary == "NA":
            continue
        families = {name for name in alternatives if name and name != "NA"} or {primary}
        if primary not in families:
            raise ValueError("Selected family must belong to its best-hit family set")
        total += 1
        selected[primary] += 1
        candidates.update(families)
        if len(families) == 1:
            unique[primary] += 1
        else:
            ambiguous += 1
            geometry_ok = geometry_ok and compatible
    if not total:
        return TEFamilyEvidence("NA", {}, 0.0, "unassigned", resolution="unassigned")
    primary = min(selected, key=lambda name: (-selected[name], name))
    if not ambiguous:
        count = selected[primary]
        status = (
            "unique"
            if len(selected) == 1
            else "dominant"
            if count * 2 > total
            else "ambiguous"
        )
        return TEFamilyEvidence(
            primary,
            dict(sorted(selected.items(), key=lambda r: (-r[1], r[0]))),
            count / total,
            status,
        )
    resolution = "unresolved_ties"
    if fallback_primary is not None:
        if fallback_primary not in selected:
            raise ValueError("Fallback family must occur in selected read evidence")
        primary = fallback_primary
    if unique:
        winner = min(unique, key=lambda name: (-unique[name], name))
        if (
            geometry_ok
            and unique[winner] >= 2
            and unique[winner] * 2 > sum(unique.values())
            and candidates[winner] * 2 > total
        ):
            primary = winner
            resolution = "unique_junction_consensus"
    support = dict(sorted(unique.items(), key=lambda r: (-r[1], r[0])))
    return TEFamilyEvidence(
        primary,
        support,
        unique[primary] / total,
        "ambiguous",
        ambiguous,
        resolution,
        dict(sorted(candidates.items(), key=lambda item: (-item[1], item[0]))),
    )
