# Proposed Aligner Combinations & Parameters per Benchmark

Created: 2026-08-04. Status: **proposal for later testing** (not implemented).

**Purpose:** Recommend TE-search × genome-placement aligner combinations and parameter settings to test against each benchmark panel (mPing, riceTElib multi-TE, riceTElib divergence). Grounded in the aligner algorithms and RelocaTE3's two-stage structure. Every non-default value below is a hypothesis to benchmark, not a committed change.

---

## 1. The two stages have opposite needs

RelocaTE3 aligns twice, and the right settings differ per stage:

- **Stage 1 — TE-library search** (`map` / `run --te-aligner`): find reads that contain TE sequence, including **junction reads** (part TE, part genomic flank). Needs **soft-clipping/local** alignment (so a partly-TE read still aligns), **mismatch tolerance** (critical for divergent copies), and **multi-mapping retention** (TE families are multi-copy). Over-sensitivity here mostly costs runtime, not precision, because the trim + genome stages filter.
- **Stage 2 — genome placement** (`align-genome --genome-aligner`): place the **trimmed flank** (often 10–30 bp) at its unique genomic locus. Needs **no long seed floor** (short flanks) and **accurate, near-unique** placement. Over-sensitivity here directly costs **precision** (spurious junctions → false insertions), as already observed (see §6).

`blat` is Stage-1 only (no SAM/genome mode). `bwaaln` is Stage-2 only (the parity genome backend).

## 2. Aligner algorithm summary + key knobs

Verified from installed `--help` (minimap2 2.31, bwa 0.7.19, bwa-mem2 2.2.1, bowtie2, blat). Defaults in brackets.

| Aligner | Algorithm | Soft-clip? | Multi-map knob | Mismatch/divergence knob | Seed/sensitivity knob |
|---|---|---|---|---|---|
| **minimap2 `-x sr`** | minimizer chain + DP | yes | `-N`[5, sr=20] `-p`[0.8, sr=0.5] | `-B` mismatch penalty [4; **sr=8**] `-A`[2] | `-k`[15; sr=21] `-w`[10; sr=11] `-s`[80] `-m`[40] |
| **bwa mem / bwa-mem2** | SMEM seed + SW | yes | `-a` (all) `-c`[500] | `-B`[4] `-A`[1] | `-k` seed floor [**19**] `-T` min score [30] `-r`[1.5] |
| **bwa aln** (`aln`+samse/sampe) | backtracking (BWT) | no (full read, gapped) | reports best; `samse`/`sampe` | `-n` max diff/err [**0.04**] `-M` mm penalty [3] | `-l` seed len [32] `-k` seed diff [2] `-o` gap opens [1] |
| **bowtie2 `--local`** | FM seed + SW extend | yes (`--local`) | `-k` [off] | `-N` seed mm 0/1 [0] `--mp`[6] `--ma`[2] | `-L` seed len [22] `-i`[S,1,0.75] `-D`[15] `-R`[2] `--score-min`[G,20,8] presets `--very-sensitive-local` |
| **blat** | tiled k-mer index + extend | yes (reports blocks) | reports all above `minScore` | `-minIdentity` [**90**] `-oneOff`[0] | `-tileSize`[11] `-minScore`[30] `-stepSize`[=tile] `-minMatch`[2] |

**Divergence-critical facts:**
- minimap2 `sr` sets **`-B 8`** (double the default mismatch penalty) — tuned for *low-divergence* short reads; it punishes divergent TE reads hard. Lowering `-B` is the main divergence lever for minimap2.
- bwa `aln -n` is an **edit-distance budget**: `0.04` ≈ ~4 % of read length. Raising it (e.g. `-n 0.10`) directly buys divergence tolerance; `-l`/`-k` govern the seed.
- bowtie2 seed `-N 0→1` and `--very-sensitive-local` (shorter, denser seeds) recover divergent reads.
- **blat `-minIdentity 90` caps usable divergence near ~10 %** regardless of the score/tile flags. This is *shared* by RelocaTE2 (it never sets minIdentity) and is why both callers collapse by 20 % divergence. Testing `-minIdentity 80/70` is the blat divergence lever.

## 3. Current RelocaTE3 settings (baseline to beat)

- **minimap2** TE stage: `-x sr -k 11 -w 5 -N 20 -p 0.5`; genome stage: `-x sr -k 11 -w 5` (primary only). (`align.py:104,293`)
- **bwa (mem)** TE stage: `mem -a`; genome uses **bwaaln** (`bwa aln` defaults + samse/sampe, RelocaTE2 parity). (`aligners.py`)
- **bowtie2** TE stage: `-k 20 --local` (default `--sensitive-local`). (`aligners.py:391`)
- **blat** TE stage: `-minScore=10 -tileSize=7` (RelocaTE2 parity, this branch). (`aligners.py:_blat_cmd`)

## 4. Benchmark challenges

- **mPing** — one 430 bp MITE, **0 % divergence**, 3 bp TSD; hard part is **short junction flanks** and **low-VAF somatic** reads, not TE identity. Lever: sensitive *genome* placement of short flanks.
- **riceTElib multi-TE** — 10 families, **multi-copy**, mixed TSD incl. **no-TSD Helitron**; hard part is **multi-mapping** in TE search + breadth. Lever: keep multi-mappers in Stage 1; precise Stage 2.
- **riceTElib divergence** — reads from copies **2–20 % diverged** from the canonical library. Hard part is **Stage-1 mismatch tolerance**. Lever: relax mismatch/identity knobs in TE search (and, as an experiment, blat `-minIdentity`).

## 5. Proposed combinations to test

Notation: **TE-search → genome-placement**. "P" = primary proposal, "A" = alternative.

### mPing (short-flank / low-VAF sensitivity)
- **P: minimap2(`-x sr -k11 -w5`) → bwaaln.** Keep the tuned short-seed TE search; swap genome placement to `bwa aln` (no seed floor) to catch 10–15 bp flanks minimap2/mem drop. This is the highest-value single change for mPing recall.
- **A: blat(`-minScore=10 -tileSize=7`) → bwaaln.** Full RelocaTE2-parity combo; the natural R2-vs-R3 head-to-head.
- **A: minimap2 → minimap2(`-k11 -w5`)** — current baseline, for reference.

### riceTElib multi-TE (multi-copy breadth)
- **P: bwa mem(`-a`) → bwaaln.** `-a` keeps every TE-copy hit in Stage 1; `bwa aln` places short flanks in Stage 2. Balanced multi-copy recall with RelocaTE2-style placement.
- **A: minimap2(`-x sr -k11 -w5 -N20 -p0.5`) → bwaaln.** Current TE search (already keeps 20 secondaries) with the short-flank genome swap.
- **A: bowtie2(`-k20 --local`) → bwaaln.** Retain the `--local` multimapper for junction reads; compare against mem/minimap2 for the awkward families (Helitron no-TSD, CACTA).
- Genome stage: keep **primary-only / no `-N20 -p0.5`** — §6 shows secondary-heavy genome placement tanks precision.

### riceTElib divergence (mismatch tolerance — the main experiment)
Hold genome placement fixed at **bwaaln** and sweep the **TE-search** knobs:
- **P: minimap2(`-x sr -k11 -w5 -B4 -N20 -p0.5`) → bwaaln.** Restore the default mismatch penalty (`-B 4` vs sr's 8) so divergent TE reads are not double-penalized. Expect the largest divergence-recall gain of the minimap2 options.
- **P: bwa aln TE-search(`-n 0.10 -l 20`) → bwaaln** *(needs a bwa-aln TE-search path; currently bwaaln is genome-only).* Raise the edit-distance budget and shorten the seed for divergent copies. **Flag:** implementing `bwa aln` as a Stage-1 backend is a prerequisite — scope before testing.
- **A: bowtie2(`--very-sensitive-local -N 1 -k 20`) → bwaaln.** Denser, mismatch-tolerant seeding.
- **A (experiment): blat(`-minScore=10 -tileSize=7 -minIdentity 80`) → bwaaln.** Lifts blat's ~10 % identity ceiling toward the 15–20 % levels; the only way blat reaches high divergence. Watch precision/runtime.

## 6. Cross-cutting notes & risks

- **Genome-stage over-sensitivity hurts precision.** Adding `--secondary=yes -N20 -p0.5` to the *genome* minimap2 step dropped precision 0.90→0.72 (`plans/2026-06-26-genotype-status-parity.md`). Keep Stage 2 primary-only; put sensitivity in Stage 1 + trim.
- **bwaaln for short flanks is likely a broad win** across all three panels (no `-k 19` seed floor); prioritize the `* → bwaaln` swaps first.
- **`minIdentity` and `-B` are shared limits, not bugs.** Testing them explores *beyond* RelocaTE2 — report them as "R3 exceeds R2" experiments, separate from parity.
- **Prerequisite gaps:** a `bwa aln` *TE-search* backend and per-backend parameter passthrough (e.g. blat `-minIdentity`, minimap2 `-B`) do not exist yet; several §5 proposals need small backend changes before they are testable. Scope these as their own tasks.
- **Test economically:** fix genome=bwaaln first, sweep TE-search per panel, then vary genome only for the winner. Reuse the divergence smoke slice (2 callers, div 0/20, 30×) to rank cheaply before full 324-task runs.

## 7. References
- Aligner backends: `src/RelocaTE3/aligners.py`; minimap2 driver: `src/RelocaTE3/align.py`
- Parity plan (blat params, bwaaln): `plans/2026-08-04-blat-bwaaln-rt2-parity.md`
- Trim-length recall gap: `plans/2026-07-02-trim-recall-parity.md`
- Genome-stage precision finding: `plans/2026-06-26-genotype-status-parity.md`
- Help sources: `minimap2 --help`, `bwa mem`/`bwa aln`, `bowtie2 --help`, `blat` (usage), all from the pinned benchmark envs.
