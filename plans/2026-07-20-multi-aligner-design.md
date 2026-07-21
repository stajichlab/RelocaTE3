# Multi-Aligner Support — Design

Date: 2026-07-20
Status: approved (design); implementation plans to follow
Scope: RelocaTE3 (the tool) **and** the relocate-benchmark repo

## Goal

Let a RelocaTE3 user choose the aligner for each of the two alignment stages,
and let the benchmark run RelocaTE3 under multiple aligner combinations so we can
compare their effect on TE-insertion calling. All benchmark plots and the
dashboard label each variant `RelocaTE3-<te-aligner>/<genome-aligner>`.

## Decisions (locked)

- **Per-stage selection:** `--te-aligner` (TE-library search) and
  `--genome-aligner` (genome re-alignment) are independent.
- **Backends:** minimap2 (default), bwa-mem, bwa-mem2, bowtie2, blat.
- **Benchmark variants:** explicit `(te, genome)` combos, each a registered
  caller with a display label `RelocaTE3-<te>/<genome>`.
- **characterize full-reads BAM:** stays on minimap2 regardless of
  `--genome-aligner` (isolates the variable; it measures depth, not junctions).
- **Default-enabled benchmark combos:** `minimap2/minimap2`, `bwa/bwa`,
  `bowtie2/bwa`, `blat/bwa`; others opt-in.
- **Backend sequencing:** bwa → bwa-mem2 → bowtie2 → blat, each with a contract
  test before the benchmark wires it in.

## Current architecture (what exists today)

Two alignment stages, both hardwired to minimap2 inside `align.py::Aligner`:

1. **TE-library search** (trim/`run`/`map`): `librelocate.RelocaTE.identify_TE_reads(search_tool=)`
   → `Aligner.map_minimap_library`. Rejects any non-minimap2 `search_tool` with
   `NotImplementedError`. CLI `map`/`run` expose `--aligner {minimap2,bwa}` but
   bwa is a non-functional trap (`todo/run-aligner-bwa-unsupported.md`).
2. **Genome re-alignment** (`align-genome`): `genome_align.align_to_genome(aligner=)`
   → `Aligner.map_reads_to_genome`. No aligner flag exists on `align-genome`.
3. characterize builds a full-reads-to-genome BAM (`Aligner.map_library_to_genome`;
   in the benchmark the adapter builds it directly with minimap2).

Downstream **contract** every stage output must satisfy (already assumed today):
coordinate-sorted + indexed BAM, mapped reads only, **read names preserved**, and
an **`NM` tag** (`insertions.py::_mismatch_count` reads `NM` minus indel bases).
bwa/bwa-mem2/bowtie2 all emit `NM` and preserve names natively; blat does not emit
SAM at all (see below).

## A. Aligner abstraction (RelocaTE3)

New package `src/RelocaTE3/aligners/` (or `backends.py`) defining an
`AlignerBackend` base with exactly three operations:

```python
class AlignerBackend(ABC):
    name: str                      # registry key, e.g. "bwa"
    def index(self, reference: Path, *, force: bool = False) -> None: ...
    def map_te_library(self, reads: ReadLibrary, out_bam: Path, *,
                       threads: int, tmpdir: Path | None = None) -> Path: ...
    def map_genome(self, reference: Path, fastqs: list[Path], out_bam: Path, *,
                   paired: bool, threads: int, tmpdir: Path | None = None) -> Path: ...
```

- Each method returns the standard contract BAM (sorted, indexed, mapped-only,
  `NM`, names preserved). A shared helper (`_sam_to_sorted_mapped_bam`) does the
  common `samtools view -F 0x4 | sort | index` tail so backends only produce SAM.
- **Stage-specific tuning lives inside each backend**, not in a generic params
  bag. `map_te_library` keeps multi-mappers (minimap2 `-N 20 -p 0.5`; bowtie2
  `-k`; bwa reports secondaries) because multi-copy TE families need them;
  `map_genome` uses the short-flank-sensitive settings (minimap2 `-k 11 -w 5`).
  These reproduce today's tuned minimap2 behavior exactly for the minimap2 backend.

Concrete backends:
- `MinimapBackend` — extracted verbatim from the current `Aligner` methods (must
  be byte-for-byte behavior-preserving; the acceptance test guards this).
- `BwaBackend` — `bwa index`; `bwa mem` → SAM. Native PE for `map_genome`.
- `BwaMem2Backend` — subclass of `BwaBackend` overriding the binary + index
  command (`bwa-mem2 index` / `bwa-mem2 mem`); identical SAM contract.
- `Bowtie2Backend` — `bowtie2-build`; `bowtie2 -x`. `-k N` at the TE stage to keep
  multi-mappers; `--no-unal` or downstream `-F 0x4` for mapped-only. Native PE.
- `BlatBackend` — `blat` → PSL → SAM shim (`psl2sam` or a small internal
  converter), **single-end**, **`map_te_library` only**; `map_genome` raises
  `NotImplementedError("blat genome alignment is not supported")`.

Registry + factory:

```python
BACKENDS: dict[str, type[AlignerBackend]] = {
    "minimap2": MinimapBackend, "bwa": BwaBackend, "bwamem2": BwaMem2Backend,
    "bowtie2": Bowtie2Backend, "blat": BlatBackend,
}
def get_aligner(name: str, threads: int = 1) -> AlignerBackend: ...
```

`Aligner` is retained as a thin facade that delegates to `get_aligner`, so
existing call sites keep working while the logic moves into backends.

## B. Wiring & CLI

- `identify_TE_reads(te_aligner="minimap2")` → `get_aligner(te_aligner).map_te_library(...)`.
  Remove the minimap-only `NotImplementedError`.
- `align_to_genome(genome_aligner="minimap2")` → `get_aligner(genome_aligner).map_genome(...)`.
- CLI:
  - `map`, `trim`, `run`: add `--te-aligner {minimap2,bwa,bwamem2,bowtie2,blat}`
    (default minimap2). Keep `--aligner` as a **deprecated alias** for
    `--te-aligner` (fixes `run-aligner-bwa-unsupported.md`: the flag now works).
  - `align-genome`: add `--genome-aligner {minimap2,bwa,bwamem2,bowtie2}`
    (default minimap2). blat rejected at parse time for this stage.
  - While editing `cmd_run`, **fix `cmd-run-drops-trim-thresholds.md`**: forward
    `--min-match`/`--min-trimmed`/`--mismatch` to `identify_TE_reads`.
- `pixi.toml`: add `bwa`, `bwa-mem2`, `bowtie2`, `blat` (bioconda).

## C. Testing

- Per-backend **contract test** (`tests/aligners_test.py`), parametrized over
  available backends: build a tiny synthetic reference + a handful of reads, run
  `index` + `map_te_library` / `map_genome`, assert the output BAM is
  coordinate-sorted, indexed, mapped-only, carries `NM`, and preserves read
  names. Each backend `skipif` the binary is not on `PATH` (mirrors
  `acceptance_test.py`).
- `acceptance_test.py` stays pinned to minimap2 (its recall/precision floors are
  minimap2-tuned; other aligners will shift them — that's the benchmark's job).

## D. Benchmark repo changes (relocate-benchmark)

The scoring/plot layer is already caller-agnostic (`compare_callers.py` is
"N-caller-safe"; `make_report.R` colours by whichever callers appear and routes
names through `pretty_caller`; the dashboard derives callers dynamically). So the
work is registry + adapter plumbing + label rendering.

- **`config/benchmark.toml`:** register one entry per variant. Key is fs-safe
  (`relocate3-<te>-<genome>`); fields `adapter="callers/relocate3"`, `repo`,
  `tsd`, `te_aligner`, `genome_aligner`, `label="RelocaTE3-<te>/<genome>"`,
  `enabled`. Default-enabled: `relocate3-minimap2-minimap2`, `relocate3-bwa-bwa`,
  `relocate3-bowtie2-bwa`, `relocate3-blat-bwa`.
- **`pipeline/config_env.py`:** generalize `CALLER_ENV_MAP` so any caller whose
  adapter is `callers/relocate3` maps `te_aligner→RT3_TE_ALIGNER`,
  `genome_aligner→RT3_GENOME_ALIGNER` (plus existing `RT3_REPO`, `TSD_PATTERN`).
  Emit `label` too.
- **`pipeline/run_benchmark_array.sh`:** resolve the adapter dir from the config
  `adapter` field instead of `callers/$CALLER` (so all variants share one
  adapter); export `RT3_TE_ALIGNER`, `RT3_GENOME_ALIGNER`; resolve `normalize.py`
  from the adapter too.
- **`callers/relocate3/run.sh`:** replace `RT3_ALIGNER` with `RT3_TE_ALIGNER`
  (→ `run --te-aligner`) and `RT3_GENOME_ALIGNER` (→ `align-genome
  --genome-aligner`). Full-reads BAM stays minimap2.
- **Labels (single source of truth):** carry the config `label` into the combined
  reports (a `caller_label` column or a `caller→label` sidecar the plotters read);
  extend R `pretty_caller` and add a Python twin in the dashboard so every plot +
  the dashboard render `RelocaTE3-<te>/<genome>`. Widen the palette to N callers.
- **env:** pin `bwa`, `bwa-mem2`, `bowtie2`, `blat` in the relocate3 caller pixi
  env so benchmark runs are reproducible.
- **Compute:** enabled callers × 9 samples as one SLURM array; each variant runs
  the full pipeline (aligners produce different BAMs, so no cross-variant reuse).

## E. Forward-compatibility (Rust / Nextflow)

- Backends are **stateless subprocess wrappers with explicit file I/O** → each is
  a Nextflow process with the aligner as a parameter; the aligner choice never
  enters the Python parsing (trim/insertions) slated for Rust — the aligner binary
  does the compute, the wrapper is pure orchestration.
- Staged CLI is preserved, so every `(stage, aligner)` is independently
  invocable — the unit a Nextflow scatter needs.
- No shared mutable state between backends; deterministic outputs.

## Open risks / notes

- Per-aligner multi-mapping semantics at the TE stage are the highest-risk detail
  (multi-copy TE families). Each backend's `map_te_library` params need validation
  against the contract test and, ultimately, the benchmark recall numbers.
- blat's PSL→SAM shim must set `NM` (or the trim/insertions mismatch logic must
  tolerate its absence — `_mismatch_count` already returns `None` when `NM` is
  missing, treated as passing). Decide during implementation.
- The two open todos (`run-aligner-bwa-unsupported`, `cmd-run-drops-trim-thresholds`)
  are folded into this work and can be closed by it.
