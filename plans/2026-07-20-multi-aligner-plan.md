# RelocaTE3 Multi-Aligner Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Let the user pick the aligner for each stage — `--te-aligner` (TE-library search) and `--genome-aligner` (genome re-alignment) — with backends minimap2 (default), bwa, bwa-mem2, bowtie2, blat.

**Architecture:** An `AlignerBackend` ABC with `index` / `map_te_library` / `map_genome`, one subclass per aligner, a `BACKENDS` registry and `get_aligner()` factory. `Aligner` becomes a thin facade delegating to backends. Stage-specific tuning lives inside each backend. Consumers (`identify_TE_reads`, `align_to_genome`) dispatch via the factory.

**Tech Stack:** Python 3.10+, pysam, subprocess to aligner binaries, pixi (bioconda: minimap2/bwa/bwa-mem2/bowtie2/blat/samtools), pytest.

**Design:** `plans/2026-07-20-multi-aligner-design.md`.

---

## Ground rules

- Run everything under `pixi run`.
- The **contract** every `map_*` returns: coordinate-sorted + indexed BAM, mapped-only, read names preserved, `NM` tag present. A shared helper enforces the sorted/mapped/index tail.
- `MinimapBackend` must be **behavior-preserving** vs today's `align.py` — the acceptance test (`tests/acceptance_test.py`) is the gate.
- Each backend test is `skipif` its binary is absent (mirror `acceptance_test.py`'s tool check).

---

## Task 1: AlignerBackend base, registry, factory, shared helper, MinimapBackend

**Files:**
- Create: `src/RelocaTE3/aligners.py`
- Test: `tests/aligners_test.py`
- Modify: `src/RelocaTE3/align.py` (make `Aligner` delegate; keep public method names)

**Step 1: Write the failing contract test for minimap2**

```python
# tests/aligners_test.py
import os, shutil, tempfile, unittest
from pathlib import Path
import pysam
from RelocaTE3.aligners import get_aligner, BACKENDS

def _has(tool): return shutil.which(tool) is not None

def _write_ref(path, seqs):
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n{seq}\n")

def _assert_contract(test, bam):
    test.assertTrue(Path(bam).exists())
    test.assertTrue(Path(f"{bam}.bai").exists())
    with pysam.AlignmentFile(str(bam), "rb") as fh:
        test.assertEqual(fh.header["HD"].get("SO"), "coordinate")
        for rec in fh.fetch(until_eof=True):
            test.assertFalse(rec.is_unmapped)          # mapped-only
            test.assertTrue(rec.has_tag("NM"))         # NM present
            test.assertTrue(rec.query_name)            # names preserved

class TestMinimapBackend(unittest.TestCase):
    @unittest.skipUnless(_has("minimap2") and _has("samtools"), "minimap2/samtools required")
    def test_map_genome_contract(self):
        with tempfile.TemporaryDirectory() as d:
            ref = os.path.join(d, "ref.fa")
            _write_ref(ref, {"chr1": "ACGT"*200})
            fq = os.path.join(d, "r.fq")
            Path(fq).write_text("@r1\n%s\n+\n%s\n" % ("ACGT"*15, "I"*60))
            aln = get_aligner("minimap2", threads=1)
            aln.index(Path(ref))
            bam = aln.map_genome(Path(ref), [Path(fq)], Path(d)/"out.bam", paired=False, threads=1)
            _assert_contract(self, bam)
```

**Step 2: Run to verify it fails**

Run: `pixi run pytest tests/aligners_test.py -v`
Expected: FAIL — `ModuleNotFoundError: RelocaTE3.aligners`.

**Step 3: Implement `aligners.py` — base + helper + registry + MinimapBackend**

Move the minimap2 command logic out of `align.py`. Complete skeleton:

```python
# src/RelocaTE3/aligners.py
from __future__ import annotations
import os, subprocess, tempfile
from abc import ABC, abstractmethod
from pathlib import Path
import pysam
from RelocaTE3 import ReadLibrary

def _sam_to_sorted_mapped_bam(sam: str, out_bam: Path, threads: int) -> Path:
    """Shared contract tail: drop unmapped, coordinate-sort, index."""
    out_bam = Path(out_bam)
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    tmp = str(out_bam) + ".mapped.bam"
    subprocess.run(["samtools", "view", "-b", "-F", "0x4", "-o", tmp, sam], check=True)
    pysam.sort("-@", str(threads), "-o", str(out_bam), tmp)
    pysam.index(str(out_bam))
    os.unlink(tmp)
    return out_bam

class AlignerBackend(ABC):
    name: str = ""
    def __init__(self, threads: int = 1): self.threads = threads
    @abstractmethod
    def index(self, reference: Path, *, force: bool = False) -> None: ...
    @abstractmethod
    def map_te_library(self, reads, out_bam: Path, *, threads: int | None = None,
                       tmpdir: Path | None = None) -> Path: ...
    @abstractmethod
    def map_genome(self, reference: Path, fastqs: list[Path], out_bam: Path, *,
                   paired: bool, threads: int | None = None,
                   tmpdir: Path | None = None) -> Path: ...

class MinimapBackend(AlignerBackend):
    name = "minimap2"
    # ... move index_minimap / map_minimap_library / map_reads_to_genome bodies here,
    #     using the EXACT flags currently in align.py:
    #       TE stage:     -a -x sr -k 11 -w 5 -N 20 -p 0.5   (per-direction, keep multimappers)
    #       genome stage: -a -x sr -k 11 -w 5                (short-flank sensitive)
    #     and finish through _sam_to_sorted_mapped_bam.

BACKENDS: dict[str, type[AlignerBackend]] = {"minimap2": MinimapBackend}

def get_aligner(name: str, threads: int = 1) -> AlignerBackend:
    try:
        return BACKENDS[name](threads)
    except KeyError:
        raise ValueError(f"unknown aligner {name!r}; choices: {sorted(BACKENDS)}")
```

Keep `align.py::Aligner` working by delegating its public methods to
`MinimapBackend` (so no other module breaks yet). Do **not** change flags.

**Step 4: Run tests to verify PASS**

Run: `pixi run pytest tests/aligners_test.py -v`  → PASS
Run: `pixi run pytest tests/acceptance_test.py -v`  → PASS (behavior-preserving)

**Step 5: Commit**

```bash
git add src/RelocaTE3/aligners.py src/RelocaTE3/align.py tests/aligners_test.py
git commit -m "feat(aligners): AlignerBackend base + registry + MinimapBackend (extracted)"
```

---

## Task 2: BwaBackend

**Files:** Modify `src/RelocaTE3/aligners.py`; add test class to `tests/aligners_test.py`.

**Step 1: Failing test** — copy `TestMinimapBackend` as `TestBwaBackend`, `skipUnless(_has("bwa"))`, `get_aligner("bwa")`; add a `map_te_library` case with a 2-read `ReadLibrary`.

**Step 2: Run** → FAIL (`unknown aligner 'bwa'`).

**Step 3: Implement `BwaBackend`**

```python
class BwaBackend(AlignerBackend):
    name = "bwa"
    index_cmd = "bwa"        # overridden by bwa-mem2
    map_cmd = "bwa"
    def index(self, reference, *, force=False):
        if force or not Path(f"{reference}.bwt").exists():
            subprocess.run([self.index_cmd, "index", str(reference)], check=True)
    def _mem(self, reference, read_files, out_bam, threads, tmpdir):
        threads = threads or self.threads
        self.index(reference)
        sam = os.path.join(tmpdir, "bwa.sam")
        with open(sam, "w") as fh:
            subprocess.run([self.map_cmd, "mem", "-t", str(threads), str(reference),
                            *map(str, read_files)], stdout=fh, check=True)
        return _sam_to_sorted_mapped_bam(sam, out_bam, threads)
    def map_te_library(self, reads, out_bam, *, threads=None, tmpdir=None):
        # bwa mem reports secondaries by default; keep them for multi-copy TEs.
        # Run per direction like minimap to preserve read-name/side handling.
        ...  # index TE lib, mem each of reads.left()/reads.right(), return per-side BAMs
    def map_genome(self, reference, fastqs, out_bam, *, paired, threads=None, tmpdir=None):
        read_files = fastqs[:2] if paired else fastqs
        with _tmp(tmpdir) as td:
            return self._mem(reference, read_files, out_bam, threads, td)
```

(Provide a `_tmp` contextmanager helper wrapping `tempfile.TemporaryDirectory`.)
Note: `map_te_library`'s return shape must match `map_minimap_library` (a list of
per-side BAMs) since `librelocate` consumes per-side outputs — verify against the
minimap signature and mirror it.

Register: add `"bwa": BwaBackend` to `BACKENDS`.

**Step 4: Run** `pixi run pytest tests/aligners_test.py::TestBwaBackend -v` → PASS.

**Step 5: Commit** `feat(aligners): bwa-mem backend`.

---

## Task 3: BwaMem2Backend

**Files:** `src/RelocaTE3/aligners.py`, `tests/aligners_test.py`.

- **Test:** `TestBwaMem2Backend`, `skipUnless(_has("bwa-mem2"))`.
- **Implement:** subclass — `class BwaMem2Backend(BwaBackend): name="bwamem2"; index_cmd="bwa-mem2"; map_cmd="bwa-mem2"` and override `index` to check `f"{reference}.bwt.2bit.64"`. Register `"bwamem2"`.
- Run test → PASS. Commit `feat(aligners): bwa-mem2 backend`.

---

## Task 4: Bowtie2Backend

**Files:** `src/RelocaTE3/aligners.py`, `tests/aligners_test.py`.

- **Test:** `TestBowtie2Backend`, `skipUnless(_has("bowtie2") and _has("bowtie2-build"))`.
- **Implement:**
  ```python
  class Bowtie2Backend(AlignerBackend):
      name = "bowtie2"
      def index(self, reference, *, force=False):
          if force or not Path(f"{reference}.1.bt2").exists():
              subprocess.run(["bowtie2-build", "--quiet", str(reference), str(reference)], check=True)
      def _run(self, reference, args, out_bam, threads, tmpdir):
          sam = os.path.join(tmpdir, "bt2.sam")
          subprocess.run(["bowtie2", "-p", str(threads), "-x", str(reference), "-S", sam, *args], check=True)
          return _sam_to_sorted_mapped_bam(sam, out_bam, threads)
      # map_te_library: add "-k", "20" to keep multi-mappers for multi-copy TE families
      # map_genome: paired -> "-1 R1 -2 R2"; else "-U" comma-joined
  ```
  Register `"bowtie2"`.
- Run → PASS. Commit `feat(aligners): bowtie2 backend`.

---

## Task 5: BlatBackend (TE-search only, PSL->SAM)

**Files:** `src/RelocaTE3/aligners.py`, `tests/aligners_test.py`.

- **Test:** `TestBlatBackend`, `skipUnless(_has("blat"))`; assert `map_te_library` contract; assert `map_genome` raises `NotImplementedError`.
- **Implement:**
  ```python
  class BlatBackend(AlignerBackend):
      name = "blat"
      def index(self, reference, *, force=False): pass  # blat indexes at run time
      def map_te_library(self, reads, out_bam, *, threads=None, tmpdir=None):
          # blat is single-end and outputs PSL. For each side:
          #   blat <te_lib> <reads.fa/fq> -out=psl <psl>
          #   convert PSL -> SAM (pslToSam / a small internal converter), set NM.
          ...
      def map_genome(self, *a, **k):
          raise NotImplementedError("blat genome alignment is not supported; use --genome-aligner {minimap2,bwa,bwamem2,bowtie2}")
  ```
  Decide the PSL→SAM path: prefer UCSC `pslToSam`/`psl2sam.pl` if pinned, else a
  minimal internal converter that emits `NM` (if `NM` cannot be computed, omit it —
  `insertions.py::_mismatch_count` treats missing `NM` as passing). Document the
  choice in the module docstring.
- Run → PASS. Commit `feat(aligners): blat backend (TE-search, PSL->SAM)`.

---

## Task 6: Wire backends into the two stage functions

**Files:** Modify `src/RelocaTE3/librelocate.py` (`identify_TE_reads`), `src/RelocaTE3/genome_align.py` (`align_to_genome`). Tests: extend existing `tests/insertions_test.py`/pipeline tests if they exercise these; else a targeted test.

**Step 1: Failing test** — assert `identify_TE_reads(..., te_aligner="minimap2")` still works and that an unknown aligner raises `ValueError` (not `NotImplementedError`).

**Step 2: Implement**
- `librelocate.py`: rename `search_tool` param to `te_aligner` (keep `search_tool` as a back-compat kwarg alias); replace the `if "minimap" not in search_tool` guard with `backend = get_aligner(te_aligner, self.cpu_threads)` and call `backend.map_te_library(...)`.
- `genome_align.py`: `align_to_genome(..., genome_aligner="minimap2")` builds `backend = get_aligner(genome_aligner, threads)` and calls `backend.map_genome(...)` instead of the hardcoded `Aligner.map_reads_to_genome`.

**Step 3: Run** the pipeline/insertions tests → PASS.
**Step 4: Commit** `refactor: dispatch TE-search and genome stages through aligner backends`.

---

## Task 7: CLI flags + fix the two open todos

**Files:** Modify `src/RelocaTE3/cli.py`. Tests: `tests/main_test.py` (subprocess smoke tests).

**Step 1: Failing smoke tests**
- `relocaTE3 run --help` mentions `--te-aligner`; `relocaTE3 align-genome --help` mentions `--genome-aligner`.
- `relocaTE3 run --te-aligner bwa ...` on a tiny fixture produces output (skipif no bwa) — no runtime `NotImplementedError`.

**Step 2: Implement**
- `map`/`trim`/`run` menus: add `--te-aligner` (choices = `sorted(BACKENDS)`; default `minimap2`). Keep `--aligner` as a hidden/deprecated alias mapping to the same dest with a deprecation note.
- `align-genome` menu: add `--genome-aligner` (choices minus `blat`; default `minimap2`).
- `cmd_run`: forward `--min-match`→`len_cut_match`, `--min-trimmed`→`len_cut_trim`, `--mismatch`→`mismatch_allowance` into `identify_TE_reads` (closes `cmd-run-drops-trim-thresholds.md`) and pass `te_aligner=args.te_aligner`.
- `cmd_align_genome`: pass `genome_aligner=args.genome_aligner`.

**Step 3: Run** `pixi run pytest tests/main_test.py -v` → PASS.
**Step 4: Commit** `feat(cli): --te-aligner/--genome-aligner; forward trim thresholds; deprecate --aligner`.

---

## Task 8: Dependencies, docs, close todos

**Files:** `pixi.toml`; `README.md`; `docs/source/usage.rst`; `todo/run-aligner-bwa-unsupported.md`, `todo/cmd-run-drops-trim-thresholds.md`, `todo/TODO_REGISTRY.md`.

- `pixi.toml`: add `bwa`, `bwa-mem2`, `bowtie2`, `blat` (bioconda). `pixi install`.
- README + usage.rst: document `--te-aligner`/`--genome-aligner`, the backend list, and blat's TE-search-only limitation; update the migration table.
- Mark both todos resolved (or delete) and update `TODO_REGISTRY.md`.
- Run full suite: `pixi run test` → PASS (acceptance still minimap2-green).
- Commit `docs+deps: aligner backends (bwa/bwa-mem2/bowtie2/blat), close aligner todos`.

---

## Verification (end of plan)

- `pixi run pytest tests/aligners_test.py -v` — every installed backend passes the contract test.
- `pixi run pytest tests/acceptance_test.py` — unchanged recall/precision (minimap2 default).
- `relocaTE3 run --te-aligner bwa` and `relocaTE3 align-genome --genome-aligner bwa` run end to end on the small fixture.
- Then execute `plans/2026-07-20-multi-aligner-benchmark-plan.md` (in the relocate-benchmark repo) to compare aligners.
