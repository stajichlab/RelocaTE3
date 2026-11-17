# Long-read L2: alignment, signatures, clustering — implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** Add `lr-align`, `lr-signatures` and `lr-cluster` to RelocaTE3, built on
a multi-reference-ready data model, and show >= 95 % site recall at 20x on the
small Chr3 long-read panel.

**Architecture:** A new subpackage `src/RelocaTE3/longread/`. The BAM adapter
turns alignments into `PlacementRecord`s; pure functions turn those into
`Signature`s and `Cluster`s; versioned TSV.gz files join the stages. One
registration call in `src/RelocaTE3/cli.py` is the only edit to existing code.

**Tech Stack:** Python 3.12 (pixi env), pysam, minimap2 2.31, samtools, pytest.

Design: `plans/2026-11-17-long-read-l2-design.md` (approved 2026-11-17).

---

## Fixed facts

| Name | Value |
|---|---|
| Worktree | `/bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/RelocaTE3_jason/RelocaTE3-lr-l2`, branch `feat/long-read-l2` from main `5740931` |
| Genome | `tests/data/sim_genome/MSU7.Chr3_2M.fa`, one sequence `Chr3` |
| Panel | `/rhome/nmath020/wessler_bigdata/rice/make_simulated_genome/make_simulation_new/results/chr3_lr_small`; link it as `tests/data/longread_small` (gitignored) |
| Panel manifest | `<panel>/panel_manifest.tsv`: `sample platform coverage replicate reads control_reads` (paths relative to the panel) |
| Truth | `<panel>/truth_events.tsv`: `event_id chrom position te_family strand tsd tsd_length biological_class`. `position` is the 0-based insertion point on `Chr3` (TSD = `reference[position-L:position]`) |
| Inserts | 145 bp – 11 kb (median 962 bp), 14 families, all homozygous |
| Run tests | `pixi run --frozen python -m pytest -o addopts="" -q -p no:cacheprovider <files>` |
| Lint | `PRE_COMMIT_HOME=<scratch> pixi exec --spec pre-commit --spec python=3.11 pre-commit run --files <files>` |

Rules: never run heavy work on the login node (`echo $SLURM_JOB_ID` must be
non-empty, or use `sbatch`); never use bare `git stash`; existing short-read
modules are not edited except the one registration line in Task 9; every
commit ends with the session's attribution lines.

Conventions for all long-read code:
- Coordinates are 0-based, half-open. An insertion "at `pos`" lies between
  reference bases `pos-1` and `pos`, so its `Locus` is `start = end = pos`.
- Sequences in signatures are in **reference-forward orientation** (as
  stored in the BAM). Query coordinates are in **read orientation** (as
  sequenced), counting hard clips.
- Nothing is keyed by sequence name alone: use `(reference_id, sequence_id)`.

---

### Task 0: Worktree baseline

**Step 1:** `ln -s <panel> tests/data/longread_small`; confirm `git status`
does not list it.

**Step 2:** Run the full suite (compute node): record the counts. Expected
(from PR #72): about 757 passed / 49 skipped, 0 failed. If
`test_te_sam_header_keeps_a_complete_fai` fails on this fresh checkout, run
`pixi run --frozen samtools faidx tests/data/TE_lib/RiceTE.fa` and rerun (known
fresh-checkout ordering issue, not ours).

No commit.

---

### Task 1: Data model (`longread/model.py`)

**Files:** Create `src/RelocaTE3/longread/__init__.py` (docstring only),
`src/RelocaTE3/longread/model.py`, `tests/longread_model_test.py`.

**Step 1 (red): tests**

```python
"""Long-read data model: loci, references, placements, signatures, clusters."""

from __future__ import annotations

import dataclasses

import pytest

from RelocaTE3.longread.model import (
    ClipSide, Cluster, Locus, PlacementRecord, ReadType, Reference,
    Signature, SignatureKind,
)  # fmt: skip


def test_read_type_presets_cover_the_four_types():
    assert {t.value: t.minimap2_preset for t in ReadType} == {
        "ont": "map-ont", "ont-hq": "lr:hq", "hifi": "map-hifi", "clr": "map-pb",
    }  # fmt: skip


@pytest.mark.parametrize(
    ("read_type", "window"),
    [("hifi", 50), ("ont-hq", 50), ("ont", 100), ("clr", 100)],
)
def test_default_cluster_window_by_read_type(read_type, window):
    assert ReadType(read_type).default_cluster_window == window


def test_locus_is_frozen_and_validates():
    locus = Locus("ref", "Chr3", 10, 10, "+")
    with pytest.raises(dataclasses.FrozenInstanceError):
        locus.start = 5  # type: ignore[misc]
    with pytest.raises(ValueError, match="end"):
        Locus("ref", "Chr3", 10, 9, "+")
    with pytest.raises(ValueError, match="strand"):
        Locus("ref", "Chr3", 1, 2, "x")


def test_locus_key_includes_the_reference():
    assert Locus("a", "Chr3", 1, 1, "+").key != Locus("b", "Chr3", 1, 1, "+").key


def test_signature_kind_and_clip_side_values():
    assert {k.value for k in SignatureKind} == {"INS", "SPLIT_INS", "CLIP"}
    assert {s.value for s in ClipSide} == {"left", "right", ""}


def test_records_build():
    locus = Locus("ref", "Chr3", 100, 100, "+")
    rec = PlacementRecord(
        read_id="r1", read_length=1000, query_start=0, query_end=1000,
        locus=Locus("ref", "Chr3", 0, 900, "+"), is_supplementary=False,
        mapq=60, insertions=((100, 300, 400),), left_clip=0, right_clip=0,
        sequence="A" * 1000, has_supplementary=False,
    )  # fmt: skip
    sig = Signature(
        signature_id="s1", kind=SignatureKind.INS, locus=locus, insert_length=100,
        clip_side=ClipSide.NONE, read_id="r1", query_start=300, query_end=400,
        mapq=60, sequence="A" * 100, left_flank="", right_flank="",
    )  # fmt: skip
    cl = Cluster(
        cluster_id="ref:Chr3:100", locus=locus, n_reads=1, n_ins=1, n_split_ins=0,
        n_clip_left=0, n_clip_right=0, insert_lengths=(100,), signature_ids=("s1",),
    )  # fmt: skip
    assert rec.locus.sequence_id == sig.locus.sequence_id == cl.locus.sequence_id
    assert Reference("ref", "g.fa").index is None
```

Run: expect `ModuleNotFoundError: RelocaTE3.longread`.

**Step 2 (green): `model.py`**

```python
"""Long-read data model (milestone L2).

Rules from the long-read proposal ("Designing for multiple references"):
every position names its reference (:class:`Locus`); evidence is separated
from the alignment format (:class:`PlacementRecord` holds no pysam objects);
signatures keep read-side evidence (sequence, flanks, query coordinates).
Coordinates are 0-based half-open; an insertion at ``pos`` has
``start == end == pos``.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class ReadType(str, Enum):
    """Long-read technology, selecting the minimap2 preset."""

    ONT = "ont"
    ONT_HQ = "ont-hq"
    HIFI = "hifi"
    CLR = "clr"

    @property
    def minimap2_preset(self) -> str:
        return {"ont": "map-ont", "ont-hq": "lr:hq", "hifi": "map-hifi",
                "clr": "map-pb"}[self.value]  # fmt: skip

    @property
    def default_cluster_window(self) -> int:
        """Breakpoint join distance; noisier reads shift breakpoints more."""
        return 50 if self in (ReadType.HIFI, ReadType.ONT_HQ) else 100


class SignatureKind(str, Enum):
    INS = "INS"
    SPLIT_INS = "SPLIT_INS"
    CLIP = "CLIP"


class ClipSide(str, Enum):
    LEFT = "left"  # clip before the aligned part (alignment start)
    RIGHT = "right"  # clip after the aligned part (alignment end)
    NONE = ""


@dataclass(frozen=True)
class Locus:
    reference_id: str
    sequence_id: str
    start: int
    end: int
    strand: str  # "+" or "-"

    def __post_init__(self) -> None:
        if self.start < 0 or self.end < self.start:
            raise ValueError(f"bad locus: start={self.start} end={self.end}")
        if self.strand not in ("+", "-"):
            raise ValueError(f"bad locus strand {self.strand!r}")

    @property
    def key(self) -> tuple[str, str]:
        """Grouping key: never the sequence name alone."""
        return (self.reference_id, self.sequence_id)


@dataclass(frozen=True)
class Reference:
    id: str
    fasta: str
    index: str | None = None


@dataclass(frozen=True)
class PlacementRecord:
    """One primary or supplementary alignment, format-independent.

    ``query_start``/``query_end`` are in read orientation and count hard
    clips. ``insertions`` holds ``(reference_pos, stored_qstart, stored_qend)``
    for CIGAR insertions, with query offsets into ``sequence``.
    ``sequence`` is the stored (reference-forward) read sequence, or ``None``
    when hard clips make it incomplete. ``left_clip``/``right_clip`` are
    soft-clip lengths at the stored left/right ends.
    """

    read_id: str
    read_length: int
    query_start: int
    query_end: int
    locus: Locus
    is_supplementary: bool
    mapq: int
    insertions: tuple[tuple[int, int, int], ...]
    left_clip: int
    right_clip: int
    sequence: str | None
    has_supplementary: bool
    left_hard_clip: int = 0
    right_hard_clip: int = 0


@dataclass(frozen=True)
class Signature:
    signature_id: str
    kind: SignatureKind
    locus: Locus
    insert_length: int
    clip_side: ClipSide
    read_id: str
    query_start: int
    query_end: int
    mapq: int
    sequence: str
    left_flank: str
    right_flank: str


@dataclass(frozen=True)
class Cluster:
    cluster_id: str
    locus: Locus
    n_reads: int
    n_ins: int
    n_split_ins: int
    n_clip_left: int
    n_clip_right: int
    insert_lengths: tuple[int, ...]
    signature_ids: tuple[str, ...]
```

Run the test file: all pass. Lint. Commit
`longread: data model (Locus, Reference, PlacementRecord, Signature, Cluster)`.

---

### Task 2: Reference space (`longread/reference.py`)

**Files:** Create `src/RelocaTE3/longread/reference.py`,
`tests/longread_reference_test.py`.

**Tests (red first):**
- `LinearFastaReferenceSpace(Reference("ref", <tmp fasta>))` on a tmp FASTA
  with `>a` (10 bp) and `>b` (5 bp): `sequences()` returns
  `[("a", 10), ("b", 5)]` in file order; `fetch(Locus("ref","a",2,5,"+"))`
  returns the 3 bases; `"-"` strand returns the reverse complement; a locus
  with another `reference_id` raises `ValueError`; an unknown sequence raises
  `KeyError`.
- Missing FASTA raises `FileNotFoundError` naming the path.
- The `.fai` is created if missing (use `pysam.faidx`).

**Implementation:**

```python
"""Reference access behind one interface (linear FASTA first; graphs later)."""

from __future__ import annotations

from pathlib import Path
from typing import Protocol

import pysam

from RelocaTE3.longread.model import Locus, Reference

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


class ReferenceSpace(Protocol):
    reference: Reference

    def sequences(self) -> list[tuple[str, int]]: ...

    def fetch(self, locus: Locus) -> str: ...


class LinearFastaReferenceSpace:
    def __init__(self, reference: Reference) -> None:
        path = Path(reference.fasta)
        if not path.is_file():
            raise FileNotFoundError(f"reference FASTA not found: {path}")
        if not Path(str(path) + ".fai").is_file():
            pysam.faidx(str(path))
        self.reference = reference
        self._fasta = pysam.FastaFile(str(path))

    def sequences(self) -> list[tuple[str, int]]:
        return list(zip(self._fasta.references, self._fasta.lengths))

    def fetch(self, locus: Locus) -> str:
        if locus.reference_id != self.reference.id:
            raise ValueError(
                f"locus is on {locus.reference_id!r}, not {self.reference.id!r}"
            )
        if locus.sequence_id not in self._fasta.references:
            raise KeyError(locus.sequence_id)
        seq = self._fasta.fetch(locus.sequence_id, locus.start, locus.end).upper()
        return reverse_complement(seq) if locus.strand == "-" else seq
```

Commit `longread: ReferenceSpace with a linear FASTA implementation`.

---

### Task 3: Versioned tables (`longread/tables.py`)

**Files:** Create `src/RelocaTE3/longread/tables.py`,
`tests/longread_tables_test.py`.

Format (both tables): gzip text. Line 1 is the schema line, tab-separated:
`#relocate3:<schema>`, `v<version>`, then `key=value` metadata
(`read_type=…`, `reference_id=…`). Line 2 is the column header. Then rows.
Writing goes to `<path>.partial` and is renamed on success; an existing
`<path>` raises `FileExistsError`.

Schemas (version 1):

```python
SIGNATURE_COLUMNS = (
    "signature_id", "kind", "reference_id", "sequence_id", "start", "end",
    "strand", "insert_length", "clip_side", "read_id", "query_start",
    "query_end", "mapq", "sequence", "left_flank", "right_flank",
)
CLUSTER_COLUMNS = (
    "cluster_id", "reference_id", "sequence_id", "start", "end", "strand",
    "n_reads", "n_ins", "n_split_ins", "n_clip_left", "n_clip_right",
    "insert_lengths", "signature_ids",
)
```

`insert_lengths` and `signature_ids` are comma-joined. Cluster strand is `+`
(L2 does not resolve element strand).

API:

```python
def write_signatures(path, signatures, *, read_type: ReadType, reference_id: str) -> int
def read_signatures(path) -> tuple[dict[str, str], Iterator[Signature]]
def write_clusters(path, clusters, *, read_type: ReadType, reference_id: str) -> int
def read_clusters(path) -> tuple[dict[str, str], Iterator[Cluster]]
```

The returned dict is the metadata (`schema`, `version`, `read_type`,
`reference_id`). Readers raise `ValueError` with a clear message for: wrong
schema name (e.g. passing a clusters file to `read_signatures`), unknown
version, or a header row that differs from the schema.

**Tests (red first):** round-trip of 3 signatures (one per kind, including an
empty `sequence` and a `-` strand) and of 2 clusters; metadata comes back;
`.partial` absent after success; existing output refused; wrong-schema,
`v2`, and altered-header files each raise `ValueError` naming the problem; a
writer that raises mid-iteration leaves no `<path>` (only cleanup of
`.partial`).

Commit `longread: versioned TSV.gz tables for signatures and clusters`.

---

### Task 4: Stage timing (`longread/timing.py`)

**Files:** Create `src/RelocaTE3/longread/timing.py`,
`tests/longread_timing_test.py`.

```python
@contextmanager
def stage_timer(outdir: Path, stage: str, sample: str) -> Iterator[StageCounts]:
    """Append one row to <outdir>/timing.tsv when the block finishes.

    Columns: stage, sample, start_utc, wall_s, max_rss_mb (self and
    children, from resource.getrusage, KiB on Linux), records_in,
    records_out, relocate3_version, status ("ok" or the exception class).
    The header is written only when the file is new.
    """
```

`StageCounts` is a small mutable dataclass with `records_in` and
`records_out` the block fills in. On an exception the row is still written
with `status=<ExceptionClass>` and the exception re-raised.

**Tests:** two stages append two rows under one header; an exception inside
writes `status=RuntimeError` and propagates; `wall_s` is a float >= 0;
`relocate3_version` equals `RelocaTE3.__version__`.

Commit `longread: per-stage timing table`.

---

### Task 5: BAM adapter (`longread/placements.py`)

**Files:** Create `src/RelocaTE3/longread/placements.py`,
`tests/longread_placements_test.py`.

```python
def placement_from_segment(seg: pysam.AlignedSegment, reference_id: str,
                           min_ins_piece: int = 10) -> PlacementRecord
def iter_placements(bam_path, reference_id: str, *, min_mapq: int,
                    stats: dict[str, int] | None = None) -> Iterator[PlacementRecord]
```

`placement_from_segment` rules (CIGAR ops: M/=/X consume both, I query only,
D/N reference only, S query only, H neither):
- `left_hard_clip`/`right_hard_clip` = H lengths at the CIGAR ends;
  `left_clip`/`right_clip` = S lengths at the ends (stored orientation).
- `read_length = seg.infer_read_length()` (counts hard clips).
- Stored aligned span on the query (including hard clips):
  `qs = left_hard + left_soft`, `qe = qs + aligned query bases`.
  Read orientation: forward `(qs, qe)`; reverse `(read_length - qe,
  read_length - qs)`.
- `insertions`: every I op with length >= `min_ins_piece`, as
  `(reference_pos_before_which_it_lies, stored_qstart, stored_qend)` where
  the query offsets index `seg.query_sequence` (soft clips included, hard
  clips excluded).
- `sequence = seg.query_sequence` (stored orientation, soft clips included,
  hard-clipped bases absent); `None` only when the BAM has no SEQ (`*`).
  A full read sequence exists only when both hard clips are 0.
- `has_supplementary = seg.has_tag("SA")`.
- `locus = Locus(reference_id, seg.reference_name, seg.reference_start,
  seg.reference_end, "-" if seg.is_reverse else "+")`.

`iter_placements` skips unmapped and secondary records and records with MAPQ
< `min_mapq`; counts `records`, `skipped_mapq`, `secondary`, `unmapped`,
`hard_clipped` into `stats`.

**Tests (red first)**, building BAMs with pysam in `tmp_path` (header with
`Chr3` length 10000; write unsorted then `pysam.sort` + `pysam.index`):
1. forward read, CIGAR `100M200I100M`: one insertion at ref `start+100`,
   stored query `(100, 300)`; query interval `(0, 400)`; no clips.
2. reverse read, `50S300M`: `left_clip=50`; read-orientation interval
   `(0, 300)` for read length 350.
3. `20H300M80S` supplementary with SA tag: `sequence is None`,
   `left_hard_clip=20`, `right_clip=80`, read length 400,
   `query_start=20, query_end=320`, `has_supplementary` true.
4. small I ops below `min_ins_piece` are not listed.
5. `iter_placements` drops a secondary (flag 256), an unmapped record and a
   MAPQ-5 record when `min_mapq=20`, and `stats` counts them.

Commit `longread: BAM adapter producing PlacementRecords`.

---

### Task 6: Signatures (`longread/signatures.py`)

**Files:** Create `src/RelocaTE3/longread/signatures.py`,
`tests/longread_signatures_test.py`.

```python
@dataclass(frozen=True)
class SignatureParams:
    min_ins_len: int = 50
    ins_merge_dist: int = 20
    min_clip_len: int = 100
    split_max_ref_gap: int = 50
    flank: int = 200

def ins_signatures(rec, params) -> list[Signature]
def clip_signatures(rec, params) -> list[Signature]
def split_ins_signatures(recs: list[PlacementRecord], read_sequence: str | None,
                         params) -> list[Signature]
def collect_signatures(placements: Iterable[PlacementRecord], params,
                       stats: dict[str, int] | None = None) -> Iterator[Signature]
```

Rules:
- **INS:** sort `rec.insertions` by reference position; merge consecutive
  pieces whose reference positions differ by <= `ins_merge_dist`. A merged
  group spans stored query `[first.qstart, last.qend)`; its length is that
  span; keep groups with length >= `min_ins_len`. Locus: `start = end =`
  first piece's reference position, strand from the record. Sequence =
  `rec.sequence[qs:qe]`; flanks = up to `flank` bases of `rec.sequence`
  either side. If `rec.sequence` is None (no SEQ), sequence and flanks are
  `""` and the length is still exact.
- **CLIP:** a soft clip >= `min_clip_len`: LEFT at `locus.start`, RIGHT at
  `locus.end` (`start == end`). `insert_length` = clip length; sequence = the
  clipped bases; the inner flank = up to `flank` aligned bases next to the
  clip, the outer flank `""`.
- **SPLIT_INS:** for one read's records with the same `(reference_id,
  sequence_id)` and strand, order them by read-orientation `query_start`.
  For each A in order, pair it with the **first later B** that passes both
  tests below (records in between that fail are skipped; this handles reads
  whose TE body aligns to another TE copy on the same sequence and strand;
  review 2026-11-17 found 12/200 panel sites otherwise lost). For each pair
  (A, B): read gap = `B.query_start -
  A.query_end`; forward strand: ref gap = `B.locus.start - A.locus.end`,
  position = `A.locus.end`; reverse strand: ref gap = `A.locus.start -
  B.locus.end`, position = `B.locus.end`. Emit when `abs(ref gap) <=
  split_max_ref_gap` and read gap >= `min_ins_len`. Insert sequence is the
  read-orientation slice `[A.query_end, B.query_start)` of the read
  sequence, reverse-complemented on `-` so it is reference-forward;
  `insert_length` = read gap. The read sequence comes from the primary
  record (convert stored to read orientation: reverse-complement when the
  primary is reverse; use a record whose hard clips are both 0). If no full
  sequence exists, sequence/flanks `""`. `mapq` = `min(A.mapq, B.mapq)`.
- **Signature query coordinates** are in read orientation, counting hard
  clips. For a stored offset `q` (into `rec.sequence`) on a forward record:
  `left_hard_clip + q`; on a reverse record: `read_length - left_hard_clip -
  q`, swapping start and end so `query_start < query_end`. INS and CLIP
  signatures use this conversion; SPLIT_INS uses the records' read-orientation
  intervals directly.
- **collect_signatures:** stream placements; records with
  `has_supplementary` false yield INS + CLIP immediately; records with it
  true yield INS + CLIP immediately **and** are buffered by `read_id`
  (keeping one full sequence per read) for SPLIT_INS at the end of the
  stream. `signature_id` = `f"{read_id}:{kind}:{n}"` with `n` counting per
  read. `stats` gets `ins`, `split_ins`, `clip`, `buffered_reads`.

**Tests (red first)** use hand-built `PlacementRecord`s (no BAM):
1. spanning insertion: one INS, correct locus, sequence, flanks;
2. ONT-style split insertion: pieces 30 bp + 40 bp, 8 bp apart → one INS of
   78 bp; pieces 60 bp apart → not merged, each below 50 → none;
3. clip-only read: right clip 500 → one CLIP RIGHT at `locus.end`;
   clip 99 → none;
4. split read forward: A `(0,1000)` → ref `[1000,2000)`, B `(1800,3000)` →
   ref `[2010,3210)`: one SPLIT_INS at 2000, length 800;
5. same on reverse strand gives position `B.locus.end` and a
   reverse-complemented sequence;
6. ref gap 60 > 50 → no SPLIT_INS; different sequences → none;
7. nested insert: two separate inserts >= 50 bp, 500 bp apart in one read →
   two INS;
8. hard-clipped supplementary paired with a full primary: SPLIT_INS has
   the right sequence from the primary;
9. `collect_signatures` over a mixed stream: counts and ids are
   deterministic.

Commit `longread: INS, SPLIT_INS and CLIP signatures`.

---

### Task 7: Clustering (`longread/cluster.py`)

**Files:** Create `src/RelocaTE3/longread/cluster.py`,
`tests/longread_cluster_test.py`.

```python
@dataclass(frozen=True)
class ClusterParams:
    window: int
    length_tolerance: float = 0.30

def cluster_signatures(signatures: Iterable[Signature], params) -> list[Cluster]
```

Algorithm, per `locus.key` (sequences in sorted order), deterministic:
1. **Insertion groups:** sort INS + SPLIT_INS by `locus.start`; single
   linkage: a signature joins the current group when its position minus the
   previous one is <= `window`.
2. **Length split:** inside each group, sort by `insert_length`; start a new
   subgroup when a length exceeds the previous by more than
   `length_tolerance` of the previous (`l > prev * (1 + tol)`).
3. **Clips:** each CLIP joins the insertion subgroup whose median position
   is nearest and within `window` (ties: the earlier one). Unattached clips
   are single-linkage grouped by position within `window` into clip-only
   clusters.
4. **Cluster record:** locus `start = end =` median position (lower median)
   of its INS/SPLIT_INS members, or of its clips for clip-only clusters;
   counts per kind; left/right clip counts; `insert_lengths` of INS/SPLIT_INS
   members sorted; `signature_ids` sorted; `n_reads` = distinct `read_id`s.
   `cluster_id = f"{reference_id}:{sequence_id}:{pos}"`, with `.1`, `.2`
   suffixes on collisions, in position order.

Notes from the Task 6 review (panel, ONT-HQ 20x): minimap2 left-aligns
gaps, so INS and LEFT CLIP sit at `truth - TSD length`, while SPLIT_INS and
RIGHT CLIP sit at `truth`; all within the window. Cluster strand is written
as `+` (the tables module writes `locus.strand` as given). `n_reads` counts
distinct reads, never a sum of signature kinds (one read can give a CLIP and
a SPLIT_INS at one site).

**Tests (red first):** two reads at 1000 and 1030 (window 50) → one cluster,
`n_reads=2`; at 1000 and 1100 → two; chained 1000/1040/1080 → one (single
linkage); lengths 1000 and 1250 → one, 1000 and 1400 → two; a right clip at
1010 attaches; a lone clip at 5000 → clip-only cluster; one read with INS and
CLIP at the same site → `n_reads=1`; two sequences / two reference ids never
merge; output order and ids are stable across input order shuffles.

Commit `longread: cluster signatures into candidate sites`.

---

### Task 8: Alignment (`longread/align.py`)

**Files:** Create `src/RelocaTE3/longread/align.py`,
`tests/longread_align_test.py`.

```python
def minimap2_command(genome, reads: list[Path], read_type: ReadType, threads: int) -> list[str]
    # ["minimap2", "-a", "-x", preset, "-Y", "--MD", "--secondary=no", "-t", N, genome, *reads]

def align_reads(genome, reads, read_type, out_bam: Path, threads: int) -> Path
    # minimap2 | samtools sort -@ N -o <out>.partial.bam - ; index; rename both.
    # Refuses an existing out_bam. Missing minimap2/samtools -> RuntimeError naming the tool.
    # On a failed pipe: remove partials, raise RuntimeError with the last stderr lines.

def check_alignment_input(bam: Path, space: ReferenceSpace) -> None
    # SO:coordinate in @HD, an index (.bai/.csi/.crai), and @SQ names+lengths
    # equal to space.sequences() (same set; order free). Each failure is a
    # ValueError naming the samtools command that fixes it.
```

**Tests (red first):** `minimap2_command` per read type (exact list);
`check_alignment_input` rejects unsorted, unindexed and wrong-`@SQ` BAMs
(built with pysam) and accepts a good one. Integration (skips if minimap2
is not on PATH): a 20 kb random genome; a "sample" with a 1 kb insert at
10 000; three 8 kb reads spanning it, written as FASTQ; `align_reads` for
each of the four read types gives a sorted indexed BAM whose
`iter_placements` + `collect_signatures` yields an INS (or SPLIT_INS)
within 20 bp of 10 000 with length within 10 % of 1000.

Commit `longread: minimap2 alignment with per-read-type presets`.

---

### Task 9: Subcommands (`longread/cli.py`)

**Files:** Create `src/RelocaTE3/longread/cli.py`,
`tests/longread_cli_test.py`. Modify `src/RelocaTE3/cli.py`: one import and
one call `register_longread_commands(subparsers)` after the existing
subcommands (the only edit to existing code).

Subcommands (all: `--sample` (required), `-o/--outdir` (required),
`-v/--verbose`):
- `lr-align`: `--genome`, `--read-type {ont,ont-hq,hifi,clr}` (required),
  exactly one of `--reads FASTQ [FASTQ ...]` / `--bam BAM`, `--threads`
  (default 1), `--reference-id` (default `ref`). FASTQ → `align_reads` to
  `<outdir>/<sample>.lr.bam`. `--bam` → `check_alignment_input`, then
  symlink it (and its index) as `<outdir>/<sample>.lr.bam`. Timing stage
  `lr-align`.
- `lr-signatures`: `--bam`, `--read-type`, `--reference-id`,
  `--min-mapq 20`, `--min-ins-len 50`, `--ins-merge-dist 20`,
  `--min-clip-len 100`, `--split-max-ref-gap 50`. Writes
  `<outdir>/<sample>.signatures.tsv.gz`; logs the adapter and signature
  stats. Timing stage `lr-signatures`.
- `lr-cluster`: `--signatures`, `--cluster-window` (default: from the file's
  `read_type`), `--length-tolerance 0.30`. Writes
  `<outdir>/<sample>.clusters.tsv.gz` with the same metadata. Timing stage
  `lr-cluster`.

Errors are raised as exceptions with clear messages; `RelocaTE3.cli.main`
already turns them into exit status 1.

**Tests (red first):** `main(["lr-align", "--help"])` exits 0 for each of
the three; `--reads` and `--bam` together are rejected; a bad
`--read-type` is rejected; end-to-end on the Task 8 synthetic data (skip
without minimap2): the three commands in sequence produce the three files,
`timing.tsv` has three `ok` rows, and the clusters file has one cluster
within 20 bp of the insert. The existing `tests/main_test.py` (or the CLI
help tests) still pass unchanged.

Commit `longread: lr-align, lr-signatures and lr-cluster subcommands`.

---

### Task 10: Gate scorer and SLURM driver

**Files:** Create `validation/longread_l2/score_sites.py`,
`validation/longread_l2/run_gate.slurm`, `validation/longread_l2/README.md`,
`tests/longread_score_sites_test.py`. Add `/validation/longread_l2/output/`
to `.gitignore`.

`score_sites.py` (standard library + `RelocaTE3.longread.tables`):

```
score_sites.py --truth <panel>/truth_events.tsv --clusters A.clusters.tsv.gz
               [--control B.clusters.tsv.gz] --max-dist 50 --label hifi_cov20x
```

- Truth sites: `(chrom, position)` from `truth_events.tsv`.
- A truth site is **found** when any cluster on the same `sequence_id` has
  `abs(cluster.start - position) <= max_dist`. Recall = found / 200.
- **Unmatched clusters:** sample clusters not within `max_dist` of any truth
  site.
- **Control:** every cluster in the reference-only control is a false
  cluster; also report the **chance recall**: the recall computed from the
  control clusters against the same truth sites (how often a site is "found"
  by noise alone).
- Each metric is also reported for clusters with `n_reads >= 2`.
- Output: one TSV row per label to stdout (`label, n_truth, found, recall,
  found_min2, recall_min2, n_clusters, unmatched, control_clusters,
  chance_recall, chance_recall_min2`).

Tests: hand-made truth (3 sites) and cluster files written with
`write_clusters`; exact counts at the 50 bp boundary (50 found, 51 not);
`n_reads >= 2` columns; control chance recall.

`run_gate.slurm` (`#!/usr/bin/bash -l`, `-p epyc`, `--cpus-per-task=8`,
`--mem=32gb`, `--time=04:00:00`, `-o logs/longread_l2_gate.%j.log`,
`set -euo pipefail`, date and parameters echoed, inputs checked): reads
`tests/data/longread_small/panel_manifest.tsv`; for each row and for both
`reads` and `control_reads`, runs `relocaTE3 lr-align` (FASTQ,
`--read-type` = the row's platform), `lr-signatures`, `lr-cluster` into
`validation/longread_l2/output/<platform>_<sample>[_control]/`, skipping a
set whose clusters file exists; then runs `score_sites.py` for each sample
with its control into `validation/longread_l2/output/gate.tsv`; prints the
table and `GATE: PASS` when recall >= 0.95 for both 20x rows, else
`GATE: FAIL`. Run `bash -n`. Commands are invoked through
`pixi run --frozen relocaTE3 ...` from the repo root.

Commit `longread: L2 gate scorer and SLURM driver`.

---

### Task 11: Run the gate, document

**Step 1:** `mkdir -p logs && sbatch validation/longread_l2/run_gate.slurm`
from the worktree root. Wait for completion (no tight polling).

**Step 2:** Read `validation/longread_l2/output/gate.tsv`. Expected: recall
>= 0.95 at 20x for `ont-hq` and `hifi`. **If it fails, stop and report the
numbers; do not tune thresholds to pass.** Investigate missed sites
(insert length, family, coverage at the site) before changing anything.

**Step 3:** Write `docs/2026-11-17-long-read-l2.md` (date/time, purpose,
status, commands, the gate table, per-stage timing summary from the
`timing.tsv` files, false-cluster numbers, decisions, next steps: L3) and
`validation/longread_l2/README.md` (how to run, what each column means).

**Step 4:** Full suite + lint. Commit `docs: long-read L2 gate results`.

---

## Not in this plan

DEL/Split-DEL (L5); TE annotation, filters, calls (L3); consensus,
sequence-based cluster splitting, TSDs (L4); `--read-type` on
`run-all`/`run-batch`/sample sheet (L3); GAF input; benchmark long-read
manifests (with L3).
