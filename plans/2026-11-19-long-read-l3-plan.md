# Long-read L3: TE annotation, filters, first calls — implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** Add `lr-annotate` and `lr-call` to RelocaTE3 so the long-read path
produces TE insertion calls with a family, and show recall >= 0.95 and
precision >= 0.95 (within 50 bp, correct family) at 20x on ONT-HQ and HiFi
on the small Chr3 long-read panel.

**Architecture:** Two new modules in `src/RelocaTE3/longread/`:
`annotate.py` (one minimap2 run per sample: signature sequences vs the TE
library -> `TeHit` rows) and `call.py` (pure functions: per-signature TE
summary, per-read family, same-family merge of clusters, vote, tier,
filters -> `Call` rows). Two new versioned tables. L2 stages are unchanged
except a reference-id sidecar written by `lr-align` and checked by
`lr-signatures`.

**Tech Stack:** Python 3.12 (pixi env), pysam, minimap2 2.31, samtools, pytest.

Design: `plans/2026-11-19-long-read-l3-design.md` (approved 2026-11-19).

---

## Fixed facts

| Name | Value |
|---|---|
| Worktree | `/bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/RelocaTE3_jason/RelocaTE3-lr-l3`, branch `feat/long-read-l3` from main `30cbf44` |
| Genome | `tests/data/sim_genome/MSU7.Chr3_2M.fa`, one sequence `Chr3` |
| Panel | `tests/data/longread_small` -> `/rhome/nmath020/wessler_bigdata/rice/make_simulated_genome/make_simulation_new/results/chr3_lr_small` (symlink exists, gitignored) |
| TE library | `<panel>/config/te_library.fa`: 14 families, 145 bp (mGing) – 10,989 bp (spmlike). Panel inserts are exact library copies (no divergence) |
| Truth | `<panel>/truth_events.tsv`: `event_id chrom position te_family strand tsd tsd_length biological_class`; 200 rows, all homozygous |
| Run tests | `pixi run --frozen python -m pytest -o addopts="" -q -p no:cacheprovider <files>` |
| Lint | `PRE_COMMIT_HOME=<scratch> pixi exec --spec pre-commit --spec python=3.11 pre-commit run --files <files>` |

Rules: run heavy work only on a compute node (`echo $SLURM_JOB_ID` non-empty)
or via `sbatch`; never bare `git stash`; short-read modules are not edited;
one implementer per worktree at a time; every commit ends with
`Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

Conventions (from L2, unchanged):
- 0-based half-open; an insertion at `pos` has `Locus(start=end=pos)`.
- `Signature.sequence` is **reference-forward**. So a TE hit on the `+`
  strand of that sequence means the TE lies forward on the reference.
- Never key by sequence name alone: use `(reference_id, sequence_id)`.

### minimap2 settings for TE annotation (measured 2026-11-19)

Measured on compute node c23 with minimap2 2.31 against the panel library:
420 simulated inserts (each family, full copy and 300 bp fragment, both
strands, 1/5/10% errors) plus 50 random sequences.

| Flags | Short TEs (mGing, Gaijin) | Other families | Random (no hit wanted) |
|---|---|---|---|
| `-c -x map-hifi` | 0/30 | 268/360 | 50/50 |
| `-c -x lr:hq` | 41/60 (fails at 5–10%) | 339/360 | 50/50 |
| `-c -x map-ont` | 52/60 | 360/360 | 50/50 |
| `-c -x map-ont -k 11 -w 5` | **60/60** | **360/360** | **50/50** |

`-r 100,100` (chaining bandwidth) is added so a nested insert (mPing inside
Retro1) is reported as three pieces instead of one Retro1 alignment with a
430 bp gap; it changed no other result. A nested insert shorter than about
150 bp (mGing inside Bajie) is still absorbed: a known limit. Final:

```
TE_MINIMAP2_ARGS = ("-c", "-x", "map-ont", "-k", "11", "-w", "5", "-r", "100,100", "--secondary=no")
```

The same flags are used for every read type (the query is a short insert,
not a read; the library is tiny). Identity = PAF column 10 / column 11.
Measured identity: >= 0.98 at 1% error, >= 0.89 at 5%, median 0.91 (min
0.66) at 10%. Starting `min_te_identity` defaults: hifi 0.85, ont-hq 0.80,
ont 0.70, clr 0.70.

---

### Task 0: Worktree baseline

**Step 1:** Confirm `tests/data/longread_small` resolves and `git status` does
not list it.

**Step 2:** Run the full suite (compute node). Record counts in the task
report. Expected about 961 passed, 0 failed. If
`test_te_sam_header_keeps_a_complete_fai` fails on the fresh checkout, run
`pixi run --frozen samtools faidx tests/data/TE_lib/RiceTE.fa` and rerun
(known, unrelated).

No commit.

---

### Task 1: Model — read-type defaults, `TeHit`, `Call`, tiers, filter reasons

**Files:**
- Modify: `src/RelocaTE3/longread/model.py`
- Test: `tests/longread_model_test.py`

**Step 1: failing tests** (append):

```python
from RelocaTE3.longread.model import (
    NO_HIT, Call, CallTier, FilterReason, TeHit,
)  # fmt: skip


@pytest.mark.parametrize(
    "read_type,support,identity",
    [(ReadType.HIFI, 2, 0.85), (ReadType.ONT_HQ, 2, 0.80),
     (ReadType.ONT, 3, 0.70), (ReadType.CLR, 3, 0.70)],
)  # fmt: skip
def test_read_type_call_defaults(read_type, support, identity) -> None:
    assert read_type.default_min_support == support
    assert read_type.default_min_te_identity == identity


def test_tier_and_reason_text() -> None:
    assert f"{CallTier.CLIP_ONLY}" == "clip_only"
    assert str(CallTier.SPANNING) == "spanning"
    assert [r.value for r in FilterReason] == [
        "no_te_hit", "few_te_reads", "family_discordant", "low_support",
        "low_te_fraction",
    ]  # fmt: skip


def _hit(**kw) -> TeHit:
    base = dict(signature_id="s1", read_id="r1", family="mPing", te_strand="+",
                te_start=0, te_end=430, te_length=430, query_start=5,
                query_end=435, query_length=440, identity=0.99)  # fmt: skip
    base.update(kw)
    return TeHit(**base)


def test_te_hit_checks() -> None:
    assert _hit().aligned == 430
    _hit(family=NO_HIT, te_strand=".", te_start=0, te_end=0, te_length=0,
         query_start=0, query_end=0, identity=0.0)  # fmt: skip
    with pytest.raises(ValueError, match="strand"):
        _hit(te_strand=".")  # a real hit needs + or -
    with pytest.raises(ValueError, match="strand"):
        _hit(family=NO_HIT, te_strand="+")
    with pytest.raises(ValueError, match="query"):
        _hit(query_start=10, query_end=5)


def test_call_passed() -> None:
    loc = Locus("ref", "Chr3", 100, 100, "+")
    common = dict(call_id="c", locus=loc, tier=CallTier.SPANNING,
                  family="mPing", te_strand="+", family_fraction=1.0,
                  strand_fraction=1.0, te_read_fraction=1.0, n_reads=3,
                  n_ins=3, n_split_ins=0, n_clip_left=0, n_clip_right=0,
                  median_insert_length=435, median_te_coverage=0.99,
                  median_identity=0.99, nested=False, cluster_ids=("k",))  # fmt: skip
    assert Call(filters=(), **common).passed
    assert not Call(filters=(FilterReason.LOW_SUPPORT,), **common).passed
    with pytest.raises(ValueError, match="insertion point"):
        Call(filters=(), **{**common, "locus": Locus("ref", "Chr3", 1, 5, "+")})
```

**Step 2:** run `tests/longread_model_test.py`; expect ImportError.

**Step 3: implement** in `model.py`:

- In `ReadType` add:

```python
    @property
    def default_min_support(self) -> int:
        """Fewest supporting reads for a PASS call (noisier reads need more)."""
        return 2 if self in (ReadType.HIFI, ReadType.ONT_HQ) else 3

    @property
    def default_min_te_identity(self) -> float:
        """Lowest TE-hit identity used (provisional; see the L3 results note)."""
        return {"hifi": 0.85, "ont-hq": 0.80, "ont": 0.70, "clr": 0.70}[self.value]
```

- New enums `CallTier` (`SPANNING="spanning"`, `CLIP_ONLY="clip_only"`) and
  `FilterReason` (values in the order tested above), each with the same
  `__str__`/`__format__` pair the other enums have. Member order is the
  order reasons are written.
- `NO_HIT = "."`.
- Frozen dataclasses:

```python
@dataclass(frozen=True)
class TeHit:
    """One alignment of a signature's sequence to the TE library.

    A signature with no alignment gets one row with ``family == NO_HIT``,
    ``te_strand == "."`` and zero coordinates, so "looked, nothing found"
    differs from "not looked". ``te_strand`` is the TE's orientation on the
    reference-forward signature sequence. ``query_*`` are on that sequence;
    ``query_length`` is its length. ``identity`` is PAF matches / block length.
    """

    signature_id: str
    read_id: str
    family: str
    te_strand: str
    te_start: int
    te_end: int
    te_length: int
    query_start: int
    query_end: int
    query_length: int
    identity: float

    def __post_init__(self) -> None:
        """Reject a strand that does not fit the hit, or inverted spans."""
        allowed = (".",) if self.family == NO_HIT else ("+", "-")
        if self.te_strand not in allowed:
            raise ValueError(f"bad te_strand {self.te_strand!r} for family {self.family!r}")
        if not 0 <= self.query_start <= self.query_end <= self.query_length:
            raise ValueError(f"bad query span {self.query_start}-{self.query_end} "
                             f"of {self.query_length}")  # fmt: skip
        if not 0 <= self.te_start <= self.te_end <= self.te_length:
            raise ValueError(f"bad TE span {self.te_start}-{self.te_end} of {self.te_length}")

    @property
    def aligned(self) -> int:
        """Query bases covered by this hit."""
        return self.query_end - self.query_start


@dataclass(frozen=True)
class Call:
    """A TE insertion call made from one or more same-family clusters.

    ``locus`` is the insertion point (``start == end``; strand ``"+"``, as
    for clusters); the TE orientation is ``te_strand`` (``"."`` when no
    family). ``filters`` is empty for a PASS call, else the failed checks in
    :class:`FilterReason` order.
    """

    call_id: str
    locus: Locus
    tier: CallTier
    filters: tuple[FilterReason, ...]
    family: str
    te_strand: str
    family_fraction: float
    strand_fraction: float
    te_read_fraction: float
    n_reads: int
    n_ins: int
    n_split_ins: int
    n_clip_left: int
    n_clip_right: int
    median_insert_length: int
    median_te_coverage: float
    median_identity: float
    nested: bool
    cluster_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        """Reject a locus that is not an insertion point."""
        if self.locus.start != self.locus.end:
            raise ValueError(f"call locus must be an insertion point: {self.locus}")

    @property
    def passed(self) -> bool:
        """True when no filter failed."""
        return not self.filters
```

**Step 4:** tests pass; whole `tests/longread_*` still pass.

**Step 5:** lint the two files; commit `longread: model for TE hits and calls`.

---

### Task 2: Tables — `lr-te-hits` and `lr-calls`

**Files:**
- Modify: `src/RelocaTE3/longread/tables.py`
- Test: `tests/longread_tables_test.py`

**Step 1: failing tests** (append; reuse the file's existing imports/helpers):

- `write_te_hits` then `read_te_hits` returns equal rows and metadata
  `schema == "lr-te-hits"`, `read_type`, `reference_id`; include a no-hit
  row. Use identities with at most 4 decimals (e.g. `0.9876`) so they
  round-trip exactly.
- `write_calls`/`read_calls` round trip: one PASS call (`filters=()`) and
  one with `(FilterReason.NO_TE_HIT, FilterReason.LOW_SUPPORT)`, `nested`
  True and False, two `cluster_ids`. Floats with at most 4 decimals.
- The `filter` column text is `PASS` for the first and
  `no_te_hit,low_support` for the second (read the gz text directly).
- `read_calls` on a te-hits file raises `ValueError` matching `wrong schema`.
- `write_calls` with a `cluster_ids` item containing a comma raises
  `ValueError`.

**Step 2:** run; expect ImportError.

**Step 3: implement:**

```python
TE_HIT_SCHEMA = "lr-te-hits"
CALL_SCHEMA = "lr-calls"

TE_HIT_COLUMNS = (
    "signature_id", "read_id", "family", "te_strand", "te_start", "te_end",
    "te_length", "query_start", "query_end", "query_length", "identity",
)  # fmt: skip
CALL_COLUMNS = (
    "call_id", "reference_id", "sequence_id", "position", "tier", "filter",
    "family", "strand", "family_fraction", "strand_fraction",
    "te_read_fraction", "n_reads", "n_ins", "n_split_ins", "n_clip_left",
    "n_clip_right", "median_insert_length", "median_te_coverage",
    "median_identity", "nested", "cluster_ids",
)  # fmt: skip
PASS = "PASS"


def _fraction(x: float) -> str:
    return f"{x:.4f}"
```

- `_te_hit_fields`: fields in column order, `identity` via `_fraction`.
  `_te_hit_from`: inverse (`float(identity)`).
- `_call_fields`: check every `cluster_ids` item with `_check_list_item`;
  `filter` = `PASS` if no filters else `",".join(str(r) for r in filters)`;
  `nested` = `"yes"`/`"no"`; the call's `locus.strand` is not written
  (always `"+"`); `strand` column holds `te_strand`; floats via `_fraction`.
- `_call_from`: inverse. `filter == "PASS"` -> `()`, else
  `tuple(FilterReason(x) for x in text.split(","))`; `nested` must be
  `yes`/`no` (else `ValueError`); locus is
  `Locus(reference_id, sequence_id, position, position, "+")`.
- `write_te_hits`, `read_te_hits`, `write_calls`, `read_calls` with the same
  signatures and docstring shape as `write_clusters`/`read_clusters`.
- Update the module docstring list of schemas.

**Step 4:** tests pass. **Step 5:** lint; commit
`longread: lr-te-hits and lr-calls tables`.

---

### Task 3: `annotate.py` — minimap2 against the TE library

**Files:**
- Create: `src/RelocaTE3/longread/annotate.py`
- Test: `tests/longread_annotate_test.py`

**Step 1: failing tests.** Unit tests need no tools; integration tests are
`skipif(shutil.which("minimap2") is None)`. Build a small library in
`tmp_path` from random sequences with a fixed seed (`random.Random(7)`):
`mGing` 145 bp, `mPing` 430 bp, `Retro1` 1,079 bp. Helper `_sig(sig_id,
read_id, seq, kind=INS)` makes a `Signature` like
`tests/longread_cluster_test.py::_sig` but with `sequence=seq`. Helper
`_mutate(seq, rate, rng)` applies substitutions/indels at `rate`.

Unit:
- `te_minimap2_command(Path("lib.fa"), Path("q.fa"), 4)` ==
  `["minimap2", "-c", "-x", "map-ont", "-k", "11", "-w", "5", "-r",
  "100,100", "--secondary=no", "-t", "4", "lib.fa", "q.fa"]`;
  `threads=0` raises `ValueError`.
- `parse_paf_line` on a hand-written PAF line
  `"s1\t440\t5\t435\t-\tmPing\t430\t0\t430\t425\t431\t60\ttp:A:P\n"` with
  `read_of={"s1": "r1"}` returns a `TeHit` with `te_strand="-"`,
  `identity == round(425/431, 4)`; a line with `tp:A:S` returns `None`.

Integration (real minimap2):
1. a 145 bp mGing copy with 5% errors plus 5 random bases -> one hit,
   family `mGing`, strand `+`;
2. the reverse complement of mPing -> family `mPing`, strand `-`;
3. 500 random bases -> exactly one row, `family == NO_HIT`;
4. `Retro1[:500] + mPing + Retro1[500:]` -> hits of both `Retro1` and
   `mPing`, the mPing hit spanning query ~500–930;
5. a CLIP signature (`kind=CLIP`, `clip_side=RIGHT`) whose sequence is the
   first 300 bp of Retro1 -> family `Retro1`;
6. a signature with `sequence=""` -> one `NO_HIT` row, `query_length=0`,
   and minimap2 is not given it (count `stats["queries"]`);
7. output order follows input signature order; within a signature, hits
   are sorted by `(query_start, query_end, family)`;
8. duplicate `signature_id` raises `ValueError`;
9. zero signatures -> `[]` without running minimap2;
10. a library path that does not exist raises `FileNotFoundError`;
11. `stats` keys: `signatures`, `queries`, `with_hit`, `without_hit`, `hits`.

**Step 2:** run; expect ModuleNotFoundError.

**Step 3: implement** `annotate.py`:

```python
"""TE annotation of long-read signatures (milestone L3, step 4).

Every signature's sequence (inserted bases; clipped bases for CLIP) is
aligned to the TE library with one minimap2 run per sample. Flags were
measured on the small Chr3 panel (plans/2026-11-19-long-read-l3-plan.md):
the HiFi preset misses 145 bp TEs entirely; ``map-ont`` with k=11, w=5 hits
every family at 1-10 % error and nothing in random sequence; ``-r 100,100``
splits nested inserts. The same flags serve every read type, because the
query is a short insert, not a read.
"""

TE_MINIMAP2_ARGS = ("-c", "-x", "map-ont", "-k", "11", "-w", "5", "-r",
                    "100,100", "--secondary=no")  # fmt: skip


def te_minimap2_command(library: Path, queries: Path, threads: int) -> list[str]:
    if threads < 1:
        raise ValueError(f"threads must be at least 1, got {threads}")
    return ["minimap2", *TE_MINIMAP2_ARGS, "-t", str(threads), str(library), str(queries)]


def parse_paf_line(line: str, read_of: Mapping[str, str]) -> TeHit | None:
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


def annotate_signatures(
    signatures: Iterable[Signature], library: Path, *, threads: int,
    workdir: Path, stats: dict[str, int] | None = None,
) -> list[TeHit]:
    ...
```

`annotate_signatures` behaviour:
- `FileNotFoundError` if `library` is not a file; `RuntimeError` if
  `minimap2` is not on PATH (only checked when there is a query).
- Collect `(signature_id, read_id, len(sequence))` in input order;
  duplicate id -> `ValueError`.
- Write non-empty sequences to `<tmpdir>/queries.fa` (one record per
  signature, name = `signature_id`); `tmpdir` is a
  `tempfile.TemporaryDirectory(dir=workdir, prefix=".lr-annotate.")`, so
  nothing is left behind on success or failure.
- Run minimap2 with stdout to `<tmpdir>/hits.paf` and stderr to
  `<tmpdir>/minimap2.err`; non-zero exit -> `RuntimeError` naming the command
  and the last 20 lines of stderr.
- Parse every line, group by `signature_id`, sort each group by
  `(query_start, query_end, family)`. Emit per input signature: its hits,
  or one `_no_hit` row.
- Fill `stats`.

**Step 4:** tests pass. **Step 5:** lint; commit
`longread: annotate signatures against the TE library (minimap2)`.

---

### Task 4: `call.py` part 1 — signature and read TE summaries

**Files:**
- Create: `src/RelocaTE3/longread/call.py`
- Test: `tests/longread_call_test.py`

Expected values in tests come from the rules below, worked by hand, never
from running the code.

**Step 1: failing tests** for `summarize_signature(hits, min_identity)` and
`read_family(summaries)`:

- no hits (one `NO_HIT` row) -> `family == NO_HIT`, `coverage == 0.0`,
  `nested is False`;
- one mPing `+` hit covering 430 of 440 -> family `mPing`, strand `+`,
  coverage `430/440`;
- a hit below `min_identity` is ignored (alone -> `NO_HIT`);
- two families, larger aligned total wins; equal totals -> family name
  order (`"Gaijin" < "mPing"`), then strand `"+"` before `"-"`;
- coverage is the union of all kept hits: `[0,300)` + `[200,500)` of 500
  -> `1.0`;
- nested: Retro1 `[0,500)` + `[930,1509)`, mPing `[500,930)` -> best
  Retro1, `nested is True`; an mPing hit inside Retro1's span (overlap >
  10% of its length) -> `nested is False`; another-family hit shorter than
  50 bp never makes it nested;
- `read_family` sums `bases[(family, strand)]` over a read's signatures and
  returns the top `(family, strand)` (ties as above), or `None` when no
  signature has a kept hit.

**Step 2:** run; expect ModuleNotFoundError.

**Step 3: implement:**

```python
NESTED_MIN_BP = 50
NESTED_MAX_OVERLAP = 0.10


@dataclass(frozen=True)
class SignatureTe:
    """TE summary of one signature after the identity filter.

    ``bases`` maps (family, strand) to aligned query bases; ``identities``
    holds the identity of each kept hit of the best family.
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
    best_union = _union((h.query_start, h.query_end) for h in kept if h.family == family)
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
    total: Counter[tuple[str, str]] = Counter()
    for s in summaries:
        for key, n in s.bases:
            total[key] += n
    if not total:
        return None
    return sorted(total.items(), key=_rank)[0][0]
```

**Step 4:** tests pass. **Step 5:** lint; commit
`longread: per-signature and per-read TE summaries`.

---

### Task 5: `call.py` part 2 — merge, vote, tier, filters, ids

**Files:**
- Modify: `src/RelocaTE3/longread/call.py`
- Test: `tests/longread_call_test.py`

Rules (write them into the module docstring):

1. **Inputs:** clusters; for every signature id, a `SignatureFacts`
   (`read_id`, `kind`, `position`, `insert_length`, `clip_side`) and a
   `SignatureTe`. A cluster signature missing from either map ->
   `ValueError` naming it.
2. **Cluster family:** `read_family` per read of the cluster; the cluster
   family is the family held by most reads (ties: more summed bases of
   that family, then name). `NO_HIT` if no read has one.
3. **Merge:** walk clusters ordered by `(reference_id, sequence_id, start,
   cluster_id)`. A cluster with family F joins the open group of
   `(key, F)` if its start minus that group's last start is <=
   `merge_window`; otherwise it opens a new group for `(key, F)`. `NO_HIT`
   clusters are always their own group. Groups keep first-cluster order.
4. **Vote on the group's union of signatures:** `reads` = distinct read ids
   (first-seen order); per read `read_family`; `hit_reads` = reads with a
   family; `te_read_fraction = len(hit_reads)/n_reads`; top family as in
   rule 2 over `hit_reads`; `family_fraction = top count/len(hit_reads)`;
   strand = majority `te_strand` among reads whose family is the top family
   (tie -> `"+"`); `strand_fraction` = that count / those reads. With no
   hit reads: family `NO_HIT`, strand `"."`, all three fractions `0.0`.
5. **Tier:** `SPANNING` if the group has any INS or SPLIT_INS signature,
   else `CLIP_ONLY`.
6. **Position:** lower median of insert-signature positions (`SPANNING`),
   else of all signature positions.
7. **Medians:** `median_insert_length` = lower median of insert-signature
   `insert_length` (`SPANNING`), else `0`. `median_te_coverage` =
   `statistics.median` of `coverage` over insert signatures (`SPANNING`) or
   all signatures (`CLIP_ONLY`), rounded to 4 decimals.
   `median_identity` = median of the `identities` of every group signature
   whose family is the top family, rounded to 4 decimals, `0.0` if none.
   `nested` = at least half of `hit_reads` have a signature with
   `nested=True` (False with no hit reads).
8. **Filters,** each checked independently, written in enum order:
   `NO_TE_HIT` (no hit reads); `FEW_TE_READS` (hit reads > 0 and
   `te_read_fraction < min_te_reads_fraction`); `FAMILY_DISCORDANT` (hit
   reads > 0 and `family_fraction < min_family_fraction`); `LOW_SUPPORT`
   (`n_reads < min_support`); `LOW_TE_FRACTION` (`SPANNING` and
   `median_te_coverage < min_te_coverage`).
9. **Ids:** `"{reference_id}:{sequence_id}:{position}"`, with `.1`, `.2`
   ... for repeats in output order (as `cluster_signatures` does). Output
   ordered by `(reference_id, sequence_id, position, first cluster_id)`.

```python
@dataclass(frozen=True)
class CallParams:
    min_support: int
    min_te_identity: float
    merge_window: int
    min_te_reads_fraction: float = 0.80
    min_family_fraction: float = 0.80
    min_te_coverage: float = 0.80


@dataclass(frozen=True)
class SignatureFacts:
    read_id: str
    kind: SignatureKind
    position: int
    insert_length: int
    clip_side: ClipSide

    @classmethod
    def of(cls, s: Signature) -> SignatureFacts:
        return cls(s.read_id, s.kind, s.locus.start, s.insert_length, s.clip_side)


def call_clusters(
    clusters: Iterable[Cluster],
    facts: Mapping[str, SignatureFacts],
    te: Mapping[str, SignatureTe],
    params: CallParams,
) -> list[Call]:
    ...
```

(`min_te_identity` is applied when building `te` with
`summarize_signature`; `CallParams` carries it so the CLI logs one object.)

**Step 1: failing tests** (hand-built clusters; small helpers make
`SignatureFacts`/`SignatureTe` directly):

- 3 reads, all mPing `+`, INS -> one PASS `SPANNING` call, family `mPing`,
  fractions `1.0`, `n_reads 3`, position = lower median of the 3 positions.
- 5 reads: 4 mPing, 1 no hit -> `te_read_fraction 0.8`, PASS (0.8 is not
  < 0.8).
- 5 reads: 3 mPing, 2 Gaijin -> `FAMILY_DISCORDANT` (0.6).
- 5 reads: 3 hit, 2 none -> `FEW_TE_READS`.
- 1 read, `min_support=2` -> `LOW_SUPPORT`; with `min_support=1` PASS.
- coverage medians 0.5 -> `LOW_TE_FRACTION`; same values on a
  `CLIP_ONLY` cluster -> not filtered for coverage.
- no hits at all, 1 read, min 2 -> filters `(NO_TE_HIT, LOW_SUPPORT)` in
  that order; family `NO_HIT`, strand `"."`.
- clip-only cluster -> tier `CLIP_ONLY`, `median_insert_length 0`.
- merge: mPing clusters at 1000 and 1030 (window 50) -> one call with both
  `cluster_ids` and `n_reads` = distinct reads of both; at 1000 and 1060 ->
  two calls; mPing 1000, Gaijin 1010, mPing 1030 -> two calls (the mPing
  pair merges around the Gaijin one); clusters on `Chr1` and `Chr3` at the
  same position never merge; `NO_HIT` clusters at 1000 and 1010 -> two
  calls.
- a clip-only mPing cluster at 1040 merges into a spanning mPing call at
  1000 -> tier `SPANNING`, `n_clip_*` counts added.
- strand tie (1 `+`, 1 `-`) -> `"+"`, `strand_fraction 0.5`.
- duplicate positions -> ids `ref:Chr3:1000` and `ref:Chr3:1000.1`.
- a cluster naming a signature missing from `facts` or `te` ->
  `ValueError` matching the id.

**Step 2:** run; expect ImportError on `call_clusters`.

**Step 3:** implement to the rules. Keep it small: `_group_clusters`,
`_vote`, `_make_call`, `call_clusters`.

**Step 4:** tests pass; whole `tests/longread_*` passes.

**Step 5:** lint; commit `longread: calls from clusters (merge, vote, filters)`.

---

### Task 6: reference-id sidecar from `lr-align`, checked by `lr-signatures`

**Files:**
- Modify: `src/RelocaTE3/longread/cli.py`
- Test: `tests/longread_cli_test.py`

Behaviour:
- After a successful `lr-align` (both `--reads` and `--bam` branches),
  write `<outdir>/<sample>.lr.meta.tsv`: header `key\tvalue`, rows
  `reference_id`, `read_type`, `genome` (resolved path). Write
  `.partial` then rename; refuse an existing file (check it with the other
  outputs before aligning, so nothing is half-done).
- `lr-signatures`: sidecar path = the BAM path with its last suffix
  replaced by `.meta.tsv` (`x.lr.bam` -> `x.lr.meta.tsv`). If it exists and
  its `reference_id` differs from `--reference-id`, or its `read_type`
  differs from `--read-type`, raise `ValueError` naming both values and the
  file. If it is absent, log at INFO that no sidecar was found (a user BAM)
  and continue.

**Step 1: failing tests:**
- parser-level: write a fake sidecar beside a pysam-built BAM (use the
  helper in `tests/longread_align_test.py` or a minimal BAM) with
  `reference_id=other`; `main(["lr-signatures", ..., "--reference-id",
  "ref"])` returns 1 and the log names `other` and `ref`.
- read-type mismatch -> returns 1.
- no sidecar -> runs (returns 0).
- `needs_tools` end to end: `lr-align` writes the sidecar with the three
  keys; a second `lr-align` into the same outdir fails without touching the
  first sidecar.

**Step 2–4:** red, implement (`_write_meta`, `_check_meta` helpers in
`cli.py`), green.

**Step 5:** lint; commit `longread: reference-id sidecar between lr-align and lr-signatures`.

---

### Task 7: `lr-annotate` and `lr-call` subcommands

**Files:**
- Modify: `src/RelocaTE3/longread/cli.py`
- Test: `tests/longread_cli_test.py`

`lr-annotate`: `--sample`, `-o/--outdir`, `--signatures` (required),
`--te-library` (required), `--threads` (default 1). Reads signatures;
output `<sample>.te_hits.tsv.gz` with the signatures' `read_type` and
`reference_id`; workdir = outdir. Timing stage `lr-annotate`,
`records_in` = signatures, `records_out` = hit rows. Logs the minimap2
flags and the stats dict.

`lr-call`: `--sample`, `-o`, `--clusters`, `--signatures`, `--te-hits`
(all required); `--min-support` (default: read type), `--min-te-identity`
(default: read type), `--merge-window` (default: read type's cluster
window), `--min-te-reads-fraction`, `--min-family-fraction`,
`--min-te-coverage` (each default 0.80). All three inputs must share
`read_type` and `reference_id` (`ValueError` naming the files otherwise).
Reads signatures into `SignatureFacts` (no sequences kept); groups hits
by signature id and builds `SignatureTe` with `summarize_signature`;
signatures without any hit row -> `ValueError` ("not annotated").
Output `<sample>.calls.tsv.gz`. Timing stage `lr-call`, `records_in` =
clusters, `records_out` = calls. Logs params and counts by tier and
PASS/filtered.

Register both in `register_longread_commands` (help texts "Long reads,
4/5: annotate signatures against the TE library" and "Long reads, 5/5:
call TE insertions from clusters"); renumber the L2 help texts to `1/5`,
`2/5`, `3/5`; update the module docstring.

**Step 1: failing tests:**
- `--help` parametrize list gains `lr-annotate` and `lr-call`; top-level
  help lists them.
- missing `--signatures` file / missing `--te-library` -> return 1.
- `lr-call` with mismatched `reference_id` between clusters and
  signatures -> return 1, message names both files.
- `lr-call` refuses existing output; `lr-annotate` refuses existing output.
- `lr-call` on hand-written tables (no tools): 2 clusters, signatures and
  hits written with the Task 2 writers -> calls file with expected tiers
  and filters; timing row `lr-call` with `records_in 2`.
- `needs_tools` end to end on the synthetic 1 kb insert from
  `tests/longread_align_test.py`: run all five stages with a TE library made
  of that insert -> one PASS `SPANNING` call at `INS_AT ± 50` with the
  library family name.

**Step 2–4:** red, implement, green. Full `tests/longread_*` green.

**Step 5:** lint; commit `longread: lr-annotate and lr-call subcommands`.

---

### Task 8: panel end-to-end test

**Files:**
- Create: `tests/longread_small_panel_calls_test.py`

Uses the same panel lookup as `tests/longread_small_panel_test.py` (env var
`RELOCATE3_LR_SMALL_PANEL`, else `tests/data/longread_small`; skip with the
looked-at path when missing) and skips without minimap2/samtools. Runs the
five stages on the `hifi` `cov5x_rep1` sample (reads path from
`panel_manifest.tsv`) into `tmp_path`, threads 2. Asserts: the calls file
exists; at least 150 PASS `SPANNING` calls; at least 90% of PASS spanning
calls have a family that is one of the 14 library names; `timing.tsv` has
one `ok` row for each of the five stages. Keep it under 60 s; record the
measured time in the commit message.

Commit `tests: long-read small-panel calls end to end`.

---

### Task 9: gate scorer `validation/longread_l3/score_calls.py`

**Files:**
- Create: `validation/longread_l3/score_calls.py`
- Test: `tests/longread_score_calls_test.py`

Rules (module docstring):
- Scored calls = `tier == spanning` and PASS.
- **Match:** greedy one-to-one. Candidate pairs (truth, call) on the same
  `sequence_id` with `|position difference| <= max_dist` and equal family;
  sort by `(distance, truth index, call index)`; take a pair when both are
  unused. `tp` = pairs taken.
- `recall = tp/n_truth`, `precision = tp/n_pass`, `f1`; 4-decimal strings,
  empty when the denominator is 0.
- `tp_position_only`: the same matching ignoring family (diagnostic).
- `strand_correct`: matched pairs with `call.te_strand == truth strand`;
  `strand_accuracy = strand_correct/tp`.
- `n_filtered`: spanning calls that are not PASS.
- `clip_only_pass`: PASS clip-only calls. `clip_only_hits`: truth sites
  left unmatched by the spanning match that a PASS clip-only call matches
  (same rule, greedy).
- With `--control`: `control_pass` (PASS spanning calls on the control)
  and `control_clip_only_pass`.
- Every call row's `reference_id` must equal the file metadata (as L2's
  scorer).
- `--per-family PATH` writes `family n_truth tp recall n_pass precision`
  (n_pass = PASS spanning calls named that family; sorted by family).

```python
COLUMNS = (
    "label", "n_truth", "n_pass", "tp", "recall", "precision", "f1",
    "tp_position_only", "strand_correct", "strand_accuracy", "n_filtered",
    "clip_only_pass", "clip_only_hits", "control_pass",
    "control_clip_only_pass",
)  # fmt: skip
```

**Step 1: failing tests** (hand-built calls files via `write_calls`, truth
TSV as in `tests/longread_score_sites_test.py`): perfect calls; a wrong
family counts in `tp_position_only` but not `tp`; two calls near one truth
site -> one tp, precision 0.5; one call between two truth sites matches
the nearer only; filtered calls are not scored; a clip-only PASS call at an
unmatched site counts in `clip_only_hits`; strand mismatch lowers
`strand_accuracy`; control columns empty without `--control`; wrong
reference id raises; per-family file contents; CLI `--header` output.

**Step 2–4:** red, implement, green. **Step 5:** lint; commit
`validation: L3 call scorer`.

---

### Task 10: gate script and README

**Files:**
- Create: `validation/longread_l3/run_gate.slurm`, `validation/longread_l3/README.md`

Copy `validation/longread_l2/run_gate.slurm` and change:
- job name `rt3-lr-l3-gate`, log `logs/longread_l3_gate.%j.log`, default
  `OUT=validation/longread_l3/output`, `TE_LIBRARY=$PANEL/config/te_library.fa`
  (checked), scorer `validation/longread_l3/score_calls.py`,
  `MIN_PRECISION=0.95` beside `MIN_RECALL=0.95`.
- `run_set` adds `lr-annotate` (`--threads "$THREADS"`) and `lr-call`
  after `lr-cluster`, each skipped when its output exists. If an L2 output
  is rebuilt, later outputs must be rebuilt: refuse when a later output is
  older than an earlier one (`-nt` test) with a message to remove the set
  directory.
- Scoring writes `$OUT/gate.tsv` and `$OUT/<set>.per_family.tsv`.
- Timing summary includes the two new stages.
- Verdict: each gate platform's 20x row needs recall >= `MIN_RECALL` and
  precision >= `MIN_PRECISION`.

README: what the gate checks, how to run (`mkdir -p logs && sbatch
validation/longread_l3/run_gate.slurm`), how to read the log, overridable
variables, rerun behaviour.

Check: `bash -n validation/longread_l3/run_gate.slurm`; shellcheck if
available (`pixi exec --spec shellcheck shellcheck ...`).

Commit `validation: L3 gate script`.

---

### Task 11: run the gate (controller, not a subagent)

`mkdir -p logs && sbatch validation/longread_l3/run_gate.slurm`; wait for
it (Monitor on `sacct`). Read `gate.tsv`, per-family tables, timing.
If the gate fails, stop and report: the numbers, the per-family table, and
the filter reasons of missed sites. Do not tune thresholds silently; any
default change is a decision to bring to the maintainer.

---

### Task 12: results note, docs, PR

**Files:**
- Create: `docs/2026-11-xx-long-read-l3.md` (date of the gate run)
- Modify: `plans/2026-11-19-long-read-l3-design.md` and this plan (an
  "Implemented" pointer at the top), `docs/2026-11-18-long-read-l2.md`
  (mark the reference-id and split-cluster limits as addressed in L3)

The note: date/time, purpose, gate table, per-family summary, strand
accuracy, clip-only counts, control false positives, timing table, the
minimap2 flags table above and why, the identity defaults and that the
panel has no TE divergence (so they are untested against divergent
copies), known limits (nested inserts < ~150 bp absorbed; clip-only calls
not placed or scored until L4; read groups; per-signature memory in
`lr-call`), "Where the plan and the code differ", commands to rerun.

Full suite on a compute node; lint all changed files; push over HTTPS
(`git -c credential.helper= -c 'credential.helper=!gh auth git-credential'
push https://github.com/stajichlab/RelocaTE3.git feat/long-read-l3`); open a
PR against main. Merge only after the maintainer says so.
