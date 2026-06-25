"""Shared helpers for the real-rice validation comparison scripts.

Both Stage A (raw non-reference calls) and Stage B (characterized/genotyped
calls) reuse the same TSV io, TE-family filter, greedy nearest-neighbor
position matcher, and matplotlib/venn loaders.
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

# Schemas published here so the normalize_* scripts share a single source of
# truth (and so importing them does not chain through compare_*.py).

NONREF_COLUMNS = [
    "sample",
    "chrom",
    "start",
    "end",
    "strand",
    "te_name",
    "tsd",
    "left_junction_reads",
    "right_junction_reads",
    "left_support_reads",
    "right_support_reads",
    "source_file",
]

CHAR_COLUMNS = [
    "sample",
    "chrom",
    "start",
    "end",
    "strand",
    "te_name",
    "tsd",
    "avg_flankers",
    "spanners",
    "status",
    "source_file",
]


def read_tsv(path: Path) -> list[dict]:
    """Read a tab-delimited file with a header into a list of dicts."""
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def filter_nonref(
    rows: Iterable[dict], te_family: str, min_junction_reads: int
) -> list[dict]:
    """Case-insensitive TE-family filter + minimum-junction-reads filter."""
    needle = te_family.lower() if te_family else ""
    out: list[dict] = []
    for r in rows:
        if needle and r["te_name"].lower() != needle:
            continue
        lj = int(r.get("left_junction_reads", "0") or 0)
        rj = int(r.get("right_junction_reads", "0") or 0)
        if lj + rj < min_junction_reads:
            continue
        out.append(r)
    return out


def filter_te_family(rows: Iterable[dict], te_family: str) -> list[dict]:
    """Case-insensitive TE-family filter only (used for characterized rows)."""
    needle = te_family.lower() if te_family else ""
    return [r for r in rows if not needle or r["te_name"].lower() == needle]


def match_sample(
    r2: list[dict], r3: list[dict], window: int
) -> tuple[list[tuple[dict, dict]], list[dict], list[dict]]:
    """Greedy nearest-neighbor match within ``window`` bp on (chrom, te_family).

    TE family names are compared case-insensitively because RelocaTE2 writes
    ``mPing`` while RelocaTE3 writes ``mping``.
    """
    by_key: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for rec in r3:
        by_key[(rec["chrom"], rec["te_name"].lower())].append(rec)
    for lst in by_key.values():
        lst.sort(key=lambda x: int(x["start"]))

    used_r3: set[int] = set()
    matched: list[tuple[dict, dict]] = []
    r2_only: list[dict] = []

    for rec2 in sorted(r2, key=lambda x: (x["chrom"], int(x["start"]))):
        candidates = by_key.get((rec2["chrom"], rec2["te_name"].lower()), [])
        s2 = int(rec2["start"])
        best_idx = -1
        best_dist = window + 1
        for i, rec3 in enumerate(candidates):
            if id(rec3) in used_r3:
                continue
            d = abs(int(rec3["start"]) - s2)
            if d <= window and d < best_dist:
                best_idx = i
                best_dist = d
        if best_idx >= 0:
            chosen = candidates[best_idx]
            used_r3.add(id(chosen))
            matched.append((rec2, chosen))
        else:
            r2_only.append(rec2)

    r3_only = [rec for rec in r3 if id(rec) not in used_r3]
    return matched, r2_only, r3_only


def write_rows(path: Path, rows: list[dict]) -> None:
    """Write a list of dicts as a TSV; empty list -> empty file."""
    if not rows:
        path.write_text("")
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def write_nonref_pairs(
    path: Path, matched: list[tuple[dict, dict]], r2_only: list[dict]
) -> None:
    """Write the Stage-A matched_calls.tsv (R2 rows + matched R3, or unmatched)."""
    cols = [
        "sample",
        "chrom",
        "te_name",
        "r2_start",
        "r2_end",
        "r2_strand",
        "r2_tsd",
        "r2_left_jr",
        "r2_right_jr",
        "r3_start",
        "r3_end",
        "r3_strand",
        "r3_tsd",
        "r3_left_jr",
        "r3_right_jr",
        "distance_bp",
        "matched",
    ]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for rec2, rec3 in matched:
            w.writerow(
                {
                    "sample": rec2["sample"],
                    "chrom": rec2["chrom"],
                    "te_name": rec2["te_name"],
                    "r2_start": rec2["start"],
                    "r2_end": rec2["end"],
                    "r2_strand": rec2["strand"],
                    "r2_tsd": rec2["tsd"],
                    "r2_left_jr": rec2["left_junction_reads"],
                    "r2_right_jr": rec2["right_junction_reads"],
                    "r3_start": rec3["start"],
                    "r3_end": rec3["end"],
                    "r3_strand": rec3["strand"],
                    "r3_tsd": rec3["tsd"],
                    "r3_left_jr": rec3["left_junction_reads"],
                    "r3_right_jr": rec3["right_junction_reads"],
                    "distance_bp": abs(int(rec3["start"]) - int(rec2["start"])),
                    "matched": 1,
                }
            )
        for rec2 in r2_only:
            w.writerow(
                {
                    "sample": rec2["sample"],
                    "chrom": rec2["chrom"],
                    "te_name": rec2["te_name"],
                    "r2_start": rec2["start"],
                    "r2_end": rec2["end"],
                    "r2_strand": rec2["strand"],
                    "r2_tsd": rec2["tsd"],
                    "r2_left_jr": rec2["left_junction_reads"],
                    "r2_right_jr": rec2["right_junction_reads"],
                    "r3_start": "",
                    "r3_end": "",
                    "r3_strand": "",
                    "r3_tsd": "",
                    "r3_left_jr": "",
                    "r3_right_jr": "",
                    "distance_bp": "",
                    "matched": 0,
                }
            )


def load_venn():
    """Return (plt, venn2) if matplotlib + matplotlib_venn import, else (None, None)."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib_venn import venn2
    except ImportError:
        print(
            "INFO: matplotlib_venn not installed; skipping venn diagrams "
            "(install with: pixi add matplotlib-venn)",
            file=sys.stderr,
        )
        return None, None
    return plt, venn2


def draw_one_venn(plt, venn2, shared, r2_only, r3_only, title, out_png) -> None:
    """Render a single 2-set Venn diagram and save to ``out_png``."""
    fig, ax = plt.subplots(figsize=(6, 5))
    venn2(
        subsets=(r2_only, r3_only, shared),
        set_labels=("RelocaTE2", "RelocaTE3"),
        ax=ax,
    )
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def load_pyplot():
    """Return matplotlib.pyplot if importable (Agg backend), else None.

    Used by Stage B plots that aren't venn diagrams.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    return plt
