# Future Features Roadmap — Post-Parity Capabilities

> **For Claude / a future maintainer:** This is a **roadmap** document for capability work that lands *after* R3 has closed the parity gap and gained a simulated benchmark (see `plans/PERFORMANCE.md`). Do NOT attempt to execute this file task-by-task. When a feature is ready to work on, write a dated implementation plan alongside it (e.g. `plans/2026-08-XX-rust-hotspot.md`) with failing-test → fix → measurement structure, using `superpowers:writing-plans`.

## Ordering rule

Each feature below has real design surface. Do them **in order**. Skipping earlier features to reach later ones (e.g. "let's do pangenomes before Nextflow because it's exciting") ships shaky infrastructure. The ordering below is chosen so each feature's design work builds on the previous one:

1. **Rust acceleration** — profile, then rewrite hotspots. No new capabilities, just speed & memory.
2. **Nextflow orchestration** — turn R3 from a "run this CLI" tool into a "point at a sample sheet and go" pipeline.
3. **Long-read support (ONT + PacBio)** — new data modality, but the pipeline shape is unchanged.
4. **Pangenome-aware calling** — new reference model, changes what "insertion location" even means.

## Prerequisite for all four

The simulated benchmark panel from `plans/PERFORMANCE.md` Priority 1 must exist before starting any feature here. Every one of these introduces changes that could regress recall/precision; without ground truth you cannot tell whether a Rust port is correct, whether Nextflow is running the same pipeline as the CLI, whether long-read calls are real, or whether pangenome projection is honest.

---

## Feature 1: Rust acceleration

### What
Rewrite the CPU-bound / memory-heavy parts of R3 in Rust, exposed to Python via PyO3 so the rest of the pipeline stays Python. This is Phase 7 in `PLAN.md`.

### Motivation
- The tool already works. Any Rust work is pure performance, not correctness.
- Rice validation runtime is dominated by two loops that can be trivially parallelised in Rust: streaming BAM records for insertion clustering, and the per-record TE-alignment parser (`_parse_te_bam`).
- On real production datasets (hundreds of samples, whole-genome coverage) the Python loops become the bottleneck. On the current 10-sample rice validation they're fine.
- Long-term this makes R3 competitive with tools like McClintock2 and TEMP2 on runtime, not just accuracy.

### Which hotspots
Profile first — do NOT guess. Use `cProfile` on a single-sample run of `pixi run relocaTE3 run`, then sort by cumulative time. Almost certainly candidates:

- `src/RelocaTE3/insertions.py:_stream_clusters` and `_call_insertions` — iterates every BAM record, groups into clusters, does per-cluster work. Streamable, embarrassingly parallel per contig.
- `src/RelocaTE3/librelocate.py:_parse_te_bam` — same pattern, iterates the TE-library BAM.
- `src/RelocaTE3/insertions.py:_estimate_tsd_length_from_depth` — builds per-base depth pileup. Small arrays but called many times; a hot inner loop.

Not candidates (leave in Python):
- CLI plumbing, argument parsing, config loading.
- `characterize.py` — bcftools/samtools subprocess wrappers, I/O bound.
- Report / normalisation scripts.

### Design considerations
- **PyO3 or standalone crate?** PyO3 keeps the boundary clean and lets us swap Rust in without changing the Python CLI. A standalone `relocate3-core` crate under `RelocaTE3/rust/` with `maturin build` is the standard shape.
- **Data-model translation.** The Python `JunctionObservation`, `Insertion`, `TEReadAlignment` dataclasses need Rust equivalents (or serde-serialised shared structs). Keep the schema identical so tests stay valid across both implementations.
- **Feature flag.** Add a `--backend python|rust` CLI switch. Default to Python during transition. The simulated panel from `PERFORMANCE.md` must pass on BOTH backends before Rust becomes the default.
- **Build/distribution.** pixi already handles Python dependencies. Adding Rust means either shipping pre-built wheels (via `maturin` + CI matrix builds for linux-x86_64 / linux-arm64 / macOS) or requiring users to have a Rust toolchain locally. The former is much better UX.
- **Rayon for parallelism.** The per-cluster and per-contig work is embarrassingly parallel; `rayon::par_iter` gets you multi-threaded scaling with almost no code change vs a serial iterator.
- **Zero-copy pysam replacement.** `rust-htslib` reads BAMs natively without going through pysam. This avoids a Python ↔ Rust marshaling round-trip per record and is often the biggest single speedup.

### Milestones
1. **F1-M1** — profile. Produce a flame graph on a 10 × rice sample. Rank hotspots by cumulative time.
2. **F1-M2** — port `_stream_clusters` (or whatever profiling identifies as #1). Ship behind `--backend rust`. Simulated panel passes.
3. **F1-M3** — port `_parse_te_bam`. Simulated panel passes.
4. **F1-M4** — measure speedup on a large realistic dataset (100+ samples). If < 3 × on the wall-clock, revisit priorities before porting more.
5. **F1-M5** — make Rust the default. Python backend stays as a fallback / correctness reference.

### Open questions
- Do we bind to `rust-htslib` or wrap pysam? `rust-htslib` is faster but adds a build-time dep on htslib.
- How do we ship the Rust binary? maturin + PyPI wheels, or vendored via pixi?
- Rust rewrite of `_pair_breakpoints` — worth it? The function is small but called per cluster and has a nested loop. Profile first.

### Prior art
- pysam itself is C wrapping htslib.
- `mapq` (Rust) and `noodles` (Rust) are two current BAM libraries; `noodles` is pure Rust and has good BAM/CRAM/VCF support.
- McClintock2 is Python + external tools; TEMP2 is C. R3 in Rust would land between them in engineering effort and above both in runtime.

---

## Feature 2: Nextflow orchestration

### What
Wrap R3 in an nf-core-style Nextflow pipeline: sample sheet in, per-sample scatter, resource-aware execution on SLURM / AWS / GCP, containerised environment, standardised reporting. This is Phase 8 in `PLAN.md`.

### Motivation
- Users currently write ad-hoc SLURM scripts (see `validation/real_rice/run_relocate3.sh`). Every user reinvents the wheel.
- nf-core pipelines are the *de facto* standard for reproducible bioinformatics. Adopting it means R3 gets a documented, containerised workflow that survives the tool's own CLI changes.
- Scatter/gather is where R3's per-sample speed adds up — Nextflow does the SLURM job graph for you.

### Design considerations
- **Sample sheet schema.** CSV with columns `sample_id, r1_fq, r2_fq, te_library, reference_genome, [repeatmasker_out]`. Mirror nf-core conventions so users can bring their existing sheets.
- **Per-sample scatter.** Each sample runs `index-genome` → `trim` → `align-genome` → `find-insertions` → `characterize` in one job graph. `index-genome` is per-genome, so it's cached across samples pointing at the same reference.
- **Chromosome scatter (optional).** Per-chromosome parallelism inside `find-insertions` is possible (`--target chr1`, `--target chr2`, ...) if a user has one genome with many contigs. Currently a manual step; Nextflow can automate.
- **Container.** Ship a Docker/Singularity image with the pixi env baked in. That means minimap2, samtools, bcftools, bedtools, and the R3 wheel all pinned to the versions R3 was tested against.
- **Modules structure.** Follow nf-core patterns:
  - `modules/local/relocate3_trim/main.nf`
  - `modules/local/relocate3_align_genome/main.nf`
  - `modules/local/relocate3_find_insertions/main.nf`
  - `modules/local/relocate3_characterize/main.nf`
  - `subworkflows/local/one_sample.nf` chains them.
  - `workflows/relocate3.nf` iterates samples and does the reference-cached indexing.
- **Report aggregation.** After all samples finish, aggregate their `characTErized.txt` files into a single joint call set. `compare_char.py` already does this for the rice harness; the Nextflow version wraps that logic in a `report/` stage.
- **Config profiles.** `slurm`, `awsbatch`, `local`, `test` — the standard set. The `test` profile runs a tiny simulated dataset (from `PERFORMANCE.md` M1) end-to-end in < 5 min.
- **Version pinning.** `nextflow.config` pins the R3 wheel version so a pipeline run is fully reproducible.

### Milestones
1. **F2-M1** — one-sample subworkflow, local profile only, using the current CLI. Fast validation that the wrapping works.
2. **F2-M2** — multi-sample scatter with reference-cached indexing. SLURM profile added.
3. **F2-M3** — report aggregation (joint call set + a simple HTML summary).
4. **F2-M4** — nf-core lint pass; submit for community listing. Not required but signals quality.
5. **F2-M5** — deprecate `validation/real_rice/run_relocate3.sh` in favour of `nextflow run relocate3 --profile slurm`.

### Open questions
- nf-core has a formal template (`nf-core create`) that ships with lint, docs, CI. Do we bootstrap from that? Probably yes — pays for itself immediately.
- Do we need Nextflow DSL2 processes for each step, or is a monolithic per-sample process fine? DSL2 makes future changes (e.g. swapping the aligner) local; monolithic ships faster.
- Interaction with the Rust feature. If R3's CLI stays the same across backends, Nextflow is unaffected. Keep the `--backend` flag surfaced as an nf-core parameter.

### Prior art
- **nf-core/mag**, **nf-core/rnaseq**, **nf-core/sarek** — all good templates for a "one CLI tool, per-sample scatter" pattern.
- **McClintock2** ships its own pipeline glue; R3 + Nextflow would land in the same UX niche but with better container hygiene.

---

## Feature 3: Long-read support (ONT + PacBio)

### What
Extend R3 to accept long-read data (Oxford Nanopore, PacBio HiFi/CLR). Currently R3 assumes short reads: junction-based detection, paired-end supporting reads, a soft-clip signal. Long reads change the physics.

### Motivation
- Long reads span entire insertions (a 15 kb ONT read easily covers the flank-TE-flank triple). One read = one full genotype call.
- Nanopore is now cheap enough for population-scale TE surveys. R3 should be able to consume that data instead of forcing users to short-read the same samples again.
- Sim wrapper already has ONT and PacBio modules (`data_sim/simulate_data/tests/test_reads_ont.py`, `test_reads_pacbio.py`) — the simulated panel from `PERFORMANCE.md` can extend to long-read variants at low marginal cost.

### Design considerations
- **Read model fundamentally differs.** A long read can either:
  1. Span the insertion (full flank + TE + flank in one read) — genotype directly from the read.
  2. End inside the TE (partial coverage from one side) — closer to the short-read junction model.
  3. Contain multiple TE copies (nested / adjacent insertions in one read) — a case short-read R3 does not handle at all.

  R3's `trim` step assumes one TE match per read; that assumption breaks for case 3. Long-read trim needs to split reads at every TE boundary and emit *multiple* flanking segments per input read.

- **Aligner change.** minimap2 has `-ax map-ont` and `-ax map-pb` presets. The genome-align step can switch on read type. The TE-library step probably needs `-ax splice` or similar; test on the sim panel.

- **Error rates.** ONT R10 is ~99 % accurate now; older data / R9 is 90 %. PacBio HiFi is ~99.9 %. The `mismatch <= 2` gate in trim is a **short-read assumption**. Long-read trim needs a percentage-based mismatch tolerance (e.g. `mismatch / match_length < 0.05`).

- **TSD detection.** With a spanning read, R3 sees both the 5′ and 3′ flanks *and* the TSD *and* the TE in one alignment. Direct-observe the TSD from the read instead of inferring from junction breakpoints. This is *better* than short-read TSD calling but needs a new code path in `insertions.py`.

- **Genotyping.** Long-read genotyping is different: instead of counting flanking-vs-spanning reads, you count reads that contain the full insertion vs reads that span the empty locus. Adapt `characterize.py:_classify` for the long-read case rather than reusing thresholds designed around short-read pileup semantics. **This is the one place where "preserve R2 logic" no longer applies** — R2 has no long-read logic to preserve.

- **Mixed short + long input.** A user might want to combine short + long reads for the same sample (short for coverage, long for span). The pipeline shape has to handle both entering the same downstream calls without double-counting.

### Milestones
1. **F3-M1** — extend the simulated panel with ONT + PacBio variants at 5×, 15×, 30× coverage. This is `PERFORMANCE.md` M2 territory but with the sim wrapper's long-read modules.
2. **F3-M2** — new `--read-type long-ont|long-pb|short` flag on `trim`, `align-genome`. The `short` case is a strict alias for current behaviour so nothing breaks.
3. **F3-M3** — long-read trim: split reads at every TE boundary; percentage-based mismatch tolerance.
4. **F3-M4** — long-read insertion caller: direct-observation TSD detection from spanning reads.
5. **F3-M5** — long-read genotype classifier (spanning-count vs empty-locus-count instead of flanker-vs-spanner).
6. **F3-M6** — measure absolute recall/precision on the long-read panel. Target: better than the 10× short-read numbers on the same insertions (long reads should be a genuine win).

### Open questions
- Do we ship one CLI (`relocaTE3 run --read-type long-ont`) or separate binaries (`relocaTE3 run` + `relocaTE3 long-run`)? One CLI is cleaner but the internals diverge substantially; separating is more maintainable long-term.
- Should the long-read code path REPLACE `characterize.py` for long-read runs, or extend it? Extension keeps report format consistent but adds branching complexity; replacement is cleaner but two report formats.
- Reference-only vs reference-free long-read calling? Long reads support *reference-free* structural variant calling (assemble contigs, compare graphs); do we ship that mode or stick to reference-guided?

### Prior art
- **TLDR** (Ewing lab) — long-read TE caller. Reference for the "direct observation" TSD approach.
- **sniffles**, **cuteSV** — general long-read SV callers. Their read-model and genotype-classifier code is worth studying even though they're not TE-specific.
- **Iris** (McClintock2 module) — long-read TE detection.

---

## Feature 4: Pangenome-aware calling

### What
Two related capabilities:

**(a) Per-genome resolution.** Given a pangenome graph containing multiple assemblies (e.g. a rice pangenome with Nipponbare + IR64 + Zhenshan97 + N22 + ...), and short-read data from an individual, R3 reports *which specific assembly path* the observed TE sits on. For example: "the mping at chr1:X is present in Nipponbare and IR64 in this sample, absent in Zhenshan97 and N22."

**(b) Novel-insertion detection.** For reads whose implied TE insertion site does NOT match any assembly path in the graph, R3 reports a *pangenome-novel* insertion. This is the sample-specific TE polymorphism relative to the entire assembly panel, not just against one reference.

### Motivation
- Rice pangenomes are becoming standard (Rice Gene Index, RiceHRPD, etc.). Users increasingly want to call TE variants relative to a graph, not a linear reference — because a linear reference biases every downstream analysis.
- Both questions are useful, and both are hard for existing linear-reference tools to answer.
- Novel-insertion detection is a genuine research capability: it answers "does this sample carry a TE insertion that no assembled genome has seen?" which is directly relevant to breeding programmes.
- Sequential rice-genome / pan-cereal papers are already reporting graph-based genotyping — R3 needs this to stay relevant.

### Design considerations

- **Graph alignment.** R3 currently reads a linear BAM. Graph aligners produce GAM (`vg giraffe`), GAF (`GraphAligner`, `minigraph`), or projected BAMs (`vg surject`). Choosing which is the interchange format is a fundamental decision:

  - **Projected BAM** (via `vg surject`) — every alignment gets projected onto a chosen linear reference. Easiest integration: R3's downstream code changes little. Downside: information loss — path membership is not recorded in the projection.
  - **Native GAM/GAF** — full graph coordinates preserved. Downside: R3 needs a graph-reader (new dependency: `libbdsg-vg` bindings or `noodles-gaf`).

  Recommendation: start with projected BAMs plus a *sidecar TSV* recording per-read path membership from the graph aligner. Adds one column to R3's internal read model. Native GAM support is a follow-up.

- **Reference model for `find-insertions`.** Currently `find-insertions` clusters reads on `(chrom, pos)` in a linear reference. With a graph, a cluster's canonical position becomes `(node_id, offset)` or `(reference_path, chrom, pos)`. Pick one convention early and stick to it; the choice ripples through GFF output, characterize, and downstream reports.

- **Existing-TE catalogue on a graph.** RepeatMasker produces linear-reference annotations. For a pangenome, you either:
  1. RepeatMask every assembly separately and merge the annotations, keeping per-assembly provenance.
  2. Project a merged annotation onto graph coordinates.

  (1) is the correct answer — each assembly's TEs get their own row in the existing-TE table, tagged with the assembly of origin. R3's reference-TE filter then works per-assembly.

- **Per-genome resolution mechanism.** For a matched insertion, R3 checks which paths in the graph pass through the same node cluster. If a graph node is in paths {Nipponbare, IR64}, and R3 finds a mping insertion at that node, the sample carries mping-at-this-locus for those two assembly backgrounds. Requires the graph aligner to report path membership per alignment; both `vg` and `minigraph` do this.

- **Novel-insertion detection mechanism.** For an R3 insertion call, check whether the two flanking breakpoints ARE both reachable via existing paths in the graph. If they are but no path connects them (i.e., the connection is the TE itself), the insertion is *pangenome-novel*. This is basically "the graph has the flanks but not the joining edge with the TE in it." Formalising this test is the hard technical work.

- **TSD detection on a graph.** The current TSD-capture logic assumes a linear reference; on a graph, the "TSD" is the overlap between the left-flank path and the right-flank path at the insertion node. Mechanically similar but implemented differently.

- **Compute cost.** Graph alignment is 3-10 × slower than linear. Nextflow (Feature 2) becomes essential; Rust (Feature 1) helps but doesn't replace the aligner cost.

- **Downstream compatibility.** VCF is standard for TE polymorphisms in a linear world. There is no widely-adopted format for "TE insertion in path X of the pangenome graph." Options: output multi-sample VCF projected onto a reference path (loses graph info), or a new pangenome-native format (`.pan-vcf`?). Recommend both — VCF for compatibility, native TSV for graph-specific queries.

### Milestones
1. **F4-M1** — pick the graph aligner and interchange format. Sketch on one small pangenome (e.g. rice pan-cereal from OGE). Deliverable: a projected BAM with sidecar path-membership TSV.
2. **F4-M2** — teach `trim` and `align-genome` to consume the interchange format. No new insertion-calling capability yet; just plumbing.
3. **F4-M3** — per-genome resolution: for a matched insertion, output the list of assembly paths passing through the node.
4. **F4-M4** — novel-insertion detection: implement the "graph has flanks but not the joining edge" test.
5. **F4-M5** — simulated pangenome benchmark (extend `PERFORMANCE.md` panel with graph-aware truth: which assembly has each simulated insertion).
6. **F4-M6** — publish results. This is a paper-worthy capability at that point.

### Open questions
- Which pangenome tooling do we bind to? `vg` is the incumbent but complex; `minigraph` is faster but less mature; `pggb` gives you graph construction but leaves alignment as a follow-up.
- Do we require users to bring their own graph, or ship a construction workflow? Bringing-your-own is simpler but shifts the burden. Construction adds significant scope but makes R3 "the pangenome-TE tool."
- Read-level or population-level analysis? Read-level answers "which path does this read follow?"; population-level answers "which assembly background does this sample sit in for this locus?". The two questions overlap but aren't identical.

### Prior art
- **vg** (Sirén et al.) — the reference implementation of graph genome tooling. `vg call` does SV genotyping on graphs; TE-specific extensions would build on that machinery.
- **PanGenie** — pangenome genotyping tool; good design reference for population-scale calls.
- **minigraph-cactus** — graph construction from assemblies. Users producing a pangenome via this pipeline is the assumed input model.
- **Locityper**, **PanTax** — tools that answer "which strain / haplotype does this sample carry?" — analogous question, different biological domain.

---

## Cross-cutting design considerations

### Interaction of features

- **Rust + Nextflow:** the Rust backend must produce identical output to the Python backend at every subcommand boundary. Nextflow does not care which backend runs; it just needs deterministic outputs so job-caching works.
- **Long reads + Nextflow:** long-read pipelines have different resource profiles (fewer, larger jobs). The Nextflow config needs read-type-aware resource requests.
- **Pangenomes + long reads:** graph alignment of long reads is a research area unto itself. Do NOT bundle F3 and F4 into one milestone. F3 lands first with linear-reference long-read support; F4 lands later with linear-first-then-graph.

### Backwards compatibility

Every feature must ship with a `--legacy` or `--short-read-linear` mode that reproduces the current output byte-for-byte. That mode is the regression gate. Break parity in a new mode, never in the default mode until the new mode's simulated panel results justify it.

### Documentation cadence

For every feature that lands:
- Update `plans/PLAN.md` with a "DONE" bullet linking to the dated implementation plan.
- Update `docs/source/` if the CLI surface changed.
- Update the README's feature matrix table.
- Update this file (FEATURES.md) to move the feature into a "shipped" section.

---

## Cross-references

- Current parity & benchmarking roadmap: `plans/PERFORMANCE.md`
- Historical narrative & completed work: `plans/PLAN.md`
- Existing simulation infrastructure: `/bigdata/stajichlab/nmath020/github/github_tools/data_sim/simulate_data`
