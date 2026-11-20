"""Score long-read calls against a truth table of insertions (L3 gate).

**Scored calls** are the calls with ``tier == spanning`` that PASS (no
filter failed).

**Matching** is greedy one-to-one. Candidate pairs ``(truth, call)`` lie on
the same ``sequence_id``, within ``max_dist`` bp (``abs(call.locus.start -
truth position) <= max_dist``) and name the same family. Pairs are sorted by
``(distance, truth index, call index)`` and a pair is taken when neither its
truth site nor its call is used yet; ``tp`` is the number of pairs taken.

- ``recall = tp / n_truth``, ``precision = tp / n_pass`` and
  ``f1 = 2 * tp / (n_truth + n_pass)`` (the harmonic mean of the two), as
  4-decimal strings, empty when the denominator is 0.
- ``tp_position_only``: the same matching ignoring family (diagnostic).
- ``strand_correct``: matched pairs whose ``call.te_strand`` equals the truth
  ``strand``; ``strand_accuracy = strand_correct / tp``.
- ``n_filtered``: spanning calls that do not PASS.
- ``clip_only_pass``: clip-only calls that PASS. ``clip_only_hits``: truth
  sites left unmatched by the spanning match that a PASS clip-only call
  matches (same rule, family included, greedy).
- With ``--control``: ``control_pass`` (PASS spanning calls on the
  reference-only control) and ``control_clip_only_pass``; both are empty
  without it.

Matching is by sequence name only, so the truth table must name sequences of
the reference the calls were made against: the ``reference_id`` in each
calls file's metadata. Every call row must carry that ``reference_id``; a
row with another one raises ``ValueError``, as does a control whose
``reference_id`` differs from the sample's.

Every truth row is scored, whatever its ``biological_class`` (the small
panel is all homozygous). A truth ``strand`` other than ``+`` or ``-``
raises ``ValueError``.

``--per-family PATH`` writes ``family n_truth tp recall n_pass precision``
per family (``n_pass``: PASS spanning calls naming that family; ``tp``:
matched pairs of that family), one row for each family in the truth or in
the scored calls, sorted by family name. An existing ``PATH`` is refused
(``FileExistsError``).

Inputs:
    --truth     truth_events.tsv (columns include ``chrom``, ``position``
                (0-based insertion point), ``te_family`` and ``strand``)
    --calls     <sample>.calls.tsv.gz from ``relocaTE3 lr-call``
    --control   optional calls file from the reference-only control reads

Output: one TSV row on stdout (columns in ``COLUMNS``), preceded by the header
line with ``--header``.

Usage:
    score_calls.py --truth T --calls C [--control K] [--max-dist 50] \
        --label hifi_cov20x [--header] [--per-family PATH]
"""

from __future__ import annotations

import argparse
import bisect
import csv
import os
import sys
from collections import Counter
from collections.abc import Sequence
from pathlib import Path
from typing import NamedTuple

from RelocaTE3.longread.model import Call, CallTier
from RelocaTE3.longread.tables import read_calls

COLUMNS = (
    "label", "n_truth", "n_pass", "tp", "recall", "precision", "f1",
    "tp_position_only", "strand_correct", "strand_accuracy", "n_filtered",
    "clip_only_pass", "clip_only_hits", "control_pass",
    "control_clip_only_pass",
)  # fmt: skip

PER_FAMILY_COLUMNS = ("family", "n_truth", "tp", "recall", "n_pass", "precision")


class TruthSite(NamedTuple):
    """One truth insertion."""

    chrom: str
    position: int
    family: str
    strand: str


def read_truth(path: str | Path) -> list[TruthSite]:
    """Return every row of a truth table as a :class:`TruthSite`.

    Raises:
        ValueError: If a column is missing or a ``strand`` is not ``+``/``-``.
    """
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        need = {"chrom", "position", "te_family", "strand"}
        missing = need - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: truth table lacks columns {sorted(missing)}")
        sites: list[TruthSite] = []
        for row in reader:
            strand = row["strand"]
            if strand not in ("+", "-"):
                raise ValueError(
                    f"{path}: line {reader.line_num}: strand {strand!r} "
                    "is not '+' or '-'"
                )
            sites.append(
                TruthSite(row["chrom"], int(row["position"]), row["te_family"], strand)
            )
        return sites


def read_checked_calls(path: str | Path) -> tuple[str, list[Call]]:
    """Return the metadata ``reference_id`` and all calls in ``path``.

    Raises:
        ValueError: If a call's ``locus.reference_id`` differs from the
            file's metadata ``reference_id``.
    """
    meta, calls = read_calls(path)
    reference_id = meta["reference_id"]
    out: list[Call] = []
    for c in calls:
        if c.locus.reference_id != reference_id:
            raise ValueError(
                f"{path}: call {c.call_id} has reference_id "
                f"{c.locus.reference_id!r}, file metadata says {reference_id!r}"
            )
        out.append(c)
    return reference_id, out


def read_scored_calls(path: str | Path) -> tuple[list[Call], list[Call]]:
    """Return ``(all calls, scored calls)``; scored = PASS spanning."""
    _, calls = read_checked_calls(path)
    return calls, _pass_of(calls, CallTier.SPANNING)


def _pass_of(calls: Sequence[Call], tier: CallTier) -> list[Call]:
    return [c for c in calls if c.tier is tier and c.passed]


def greedy_match(
    truth: Sequence[TruthSite],
    calls: Sequence[Call],
    max_dist: int,
    use_family: bool = True,
) -> list[tuple[int, int]]:
    """Match truth sites to calls one-to-one, nearest pairs first.

    Returns:
        ``(truth index, call index)`` pairs in the order they were taken.
    """
    # Per sequence, (start, call index) sorted by start; a bisect window over
    # the starts gives the calls within max_dist of each truth site.
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for j, c in enumerate(calls):
        by_chrom.setdefault(c.locus.sequence_id, []).append((c.locus.start, j))
    starts: dict[str, list[int]] = {}
    for chrom, entries in by_chrom.items():
        entries.sort()
        starts[chrom] = [start for start, _ in entries]
    pairs: list[tuple[int, int, int]] = []
    for i, t in enumerate(truth):
        if t.chrom not in by_chrom:
            continue
        keys = starts[t.chrom]
        lo = bisect.bisect_left(keys, t.position - max_dist)
        hi = bisect.bisect_right(keys, t.position + max_dist)
        for start, j in by_chrom[t.chrom][lo:hi]:
            if not use_family or calls[j].family == t.family:
                pairs.append((abs(start - t.position), i, j))
    pairs.sort()
    used_truth: set[int] = set()
    used_calls: set[int] = set()
    taken: list[tuple[int, int]] = []
    for _, i, j in pairs:
        if i not in used_truth and j not in used_calls:
            used_truth.add(i)
            used_calls.add(j)
            taken.append((i, j))
    return taken


def _ratio(num: int, den: int) -> str:
    return f"{num / den:.4f}" if den else ""


def _write_per_family(
    path: str | Path,
    truth: Sequence[TruthSite],
    scored: Sequence[Call],
    matched: Sequence[tuple[int, int]],
) -> None:
    n_truth = Counter(t.family for t in truth)
    n_pass = Counter(c.family for c in scored)
    tp = Counter(truth[i].family for i, _ in matched)
    lines = ["\t".join(PER_FAMILY_COLUMNS)]
    for family in sorted(set(n_truth) | set(n_pass)):
        fields = (family, n_truth[family], tp[family], _ratio(tp[family], n_truth[family]),
                  n_pass[family], _ratio(tp[family], n_pass[family]))  # fmt: skip
        lines.append("\t".join(str(f) for f in fields))
    with open(path, "x") as fh:
        fh.write("\n".join(lines) + "\n")


def score(
    truth_path: str | Path,
    calls_path: str | Path,
    control_path: str | Path | None,
    max_dist: int,
    label: str = "",
    per_family: str | Path | None = None,
) -> dict[str, object]:
    """Score one sample's calls (and optional control) against the truth.

    Returns:
        A dict keyed by ``COLUMNS``; ratios are 4-decimal strings (``""``
        when the denominator is 0), control fields are ``""`` when
        ``control_path`` is None.

    Raises:
        FileExistsError: If ``per_family`` already exists.
        ValueError: On a bad truth row, a call row or control whose
            ``reference_id`` disagrees, or a negative ``max_dist``.
    """
    if max_dist < 0:
        raise ValueError(f"max_dist must be >= 0, got {max_dist}")
    if per_family is not None and os.path.lexists(per_family):
        raise FileExistsError(f"per-family output exists: {per_family}")
    truth = read_truth(truth_path)
    sample_ref, calls = read_checked_calls(calls_path)
    scored = _pass_of(calls, CallTier.SPANNING)
    control: list[Call] | None = None
    if control_path is not None:
        control_ref, control = read_checked_calls(control_path)
        if control_ref != sample_ref:
            raise ValueError(
                f"{control_path}: control reference_id {control_ref!r} differs "
                f"from the sample's {sample_ref!r} ({calls_path})"
            )
    n_truth, n_pass = len(truth), len(scored)

    matched = greedy_match(truth, scored, max_dist)
    tp = len(matched)
    strand_correct = sum(scored[j].te_strand == truth[i].strand for i, j in matched)
    matched_truth = {i for i, _ in matched}
    left = [t for i, t in enumerate(truth) if i not in matched_truth]
    clip_pass = _pass_of(calls, CallTier.CLIP_ONLY)
    row: dict[str, object] = {
        "label": label,
        "n_truth": n_truth,
        "n_pass": n_pass,
        "tp": tp,
        "recall": _ratio(tp, n_truth),
        "precision": _ratio(tp, n_pass),
        "f1": _ratio(2 * tp, n_truth + n_pass),
        "tp_position_only": len(greedy_match(truth, scored, max_dist, False)),
        "strand_correct": strand_correct,
        "strand_accuracy": _ratio(strand_correct, tp),
        "n_filtered": sum(c.tier is CallTier.SPANNING and not c.passed for c in calls),
        "clip_only_pass": len(clip_pass),
        "clip_only_hits": len(greedy_match(left, clip_pass, max_dist)),
        "control_pass": "",
        "control_clip_only_pass": "",
    }
    if control is not None:
        row["control_pass"] = len(_pass_of(control, CallTier.SPANNING))
        row["control_clip_only_pass"] = len(_pass_of(control, CallTier.CLIP_ONLY))
    if per_family is not None:
        _write_per_family(per_family, truth, scored, matched)
    return row


def main(argv: Sequence[str] | None = None) -> int:
    """Print one TSV row (and optionally the header) to stdout."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--truth", required=True, help="truth_events.tsv")
    parser.add_argument("--calls", required=True, help="sample calls.tsv.gz")
    parser.add_argument("--control", default=None,
                        help="reference-only control calls.tsv.gz")  # fmt: skip
    parser.add_argument("--max-dist", type=int, default=50,
                        help="largest |call position - truth position| (bp)")  # fmt: skip
    parser.add_argument("--label", required=True, help="row label")
    parser.add_argument("--header", action="store_true", help="print the header")
    parser.add_argument("--per-family", default=None,
                        help="write per-family recall/precision TSV here")  # fmt: skip
    args = parser.parse_args(argv)
    row = score(args.truth, args.calls, args.control, args.max_dist, args.label,
                args.per_family)  # fmt: skip
    if args.header:
        sys.stdout.write("\t".join(COLUMNS) + "\n")
    sys.stdout.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
