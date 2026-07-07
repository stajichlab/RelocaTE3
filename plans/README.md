# Plans

This directory contains RelocaTE3 development plans, both long-lived roadmaps and dated execution plans for individual features.

## Roadmap documents (start here)

- **`PLAN.md`** — historical narrative and phase-by-phase status of every completed piece of R3 (Phase 0 → 4 done, Phase 5+ ongoing).
- **`PERFORMANCE.md`** — active parity work + simulated-data benchmarking. What to work on next before capability features.
- **`FEATURES.md`** — future capability roadmap: Rust acceleration, Nextflow orchestration, long-read support, pangenome-aware calling.

## Dated execution plans

Each `YYYY-MM-DD-<slug>.md` file is a specific implementation plan with failing-test → fix → measurement structure. They stay in-tree after the work lands as an audit trail. Recent examples:

- `2026-06-24-tsd-depth-inference.md` (done — depth-based TSD in the function path)
- `2026-06-25-tsd-class-path-port.md` (done — wildcard TSD for the class path)
- `2026-06-25-tsd-supporting-junction-port.md` (deferred — remaining ~40 sites)
- `2026-06-26-genotype-status-parity.md` (retired — dead-end targeting the wrong step)
- `2026-07-02-trim-recall-parity.md` (done — reverse-strand coord-frame fix)

## Conventions

- Roadmaps use uppercase file names (`PERFORMANCE.md`, `FEATURES.md`, `PLAN.md`).
- Execution plans use `YYYY-MM-DD-<slug>.md`.
- Retired / superseded plans stay in-tree as documentation, with a "supersedes:" or "retired:" note added at the top.
- Private scratch plans or local-only notes should not be committed.

## Writing new plans

When starting a feature marked "ready to work on" in `PERFORMANCE.md` or `FEATURES.md`, write a dated execution plan alongside it using the `superpowers:writing-plans` skill. The roadmap document sets priority and scope; the execution plan drives implementation.
