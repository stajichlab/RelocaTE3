# Release first; future architecture direction

Recorded: 2026-10-02, America/Los_Angeles (session began 01:12 PDT).
Project: RelocaTE3, `main`, after merged PR #49.
Status: maintainer direction and design notes, not implemented features.

## Immediate commitment

The first release should be a modern, organized replacement for RelocaTE2
without sacrificing its utility. Stabilizing short-read results takes priority
over Rust, Nextflow, long reads, and pangenomes. Avoid expanding the already
large uncommitted change set on `main` with speculative refactoring.

Accuracy, runtime, memory, and feature availability are separate requirements.
More pooled true positives does not establish noninferiority if particular
families, low-frequency insertions, or precision regress. Compare both tools
against truth and inspect paired gains/losses, not just agreement with R2.
Use the R2 source in `../references/RelocaTE2` to explain implementation
differences; reproduce the full relevant behavior rather than copying an
isolated coordinate or threshold rule.

Release evidence must cover mPing, riceTElib, and diverged riceTElib, including
coverage, TE group, somatic fraction, genotype/status, and TSD accuracy.
The 0%-divergence riceTElib panel overlaps the canonical panel; do not present
pooled totals as independent biological replicates or a universal guarantee.
Audit the R2 feature surface separately: paired/single-end and mapped-read
inputs, reference/shared/absence calls, and optional excision outputs cannot
be certified by the non-reference paired-end benchmark alone.

Before declaring the release ready: resolve the observed regressions, run the
frozen full pipeline and full tests on compute nodes, verify actual CLI paths,
and measure peak memory and runtime with the final code. Then organize focused
commits and provide the maintainer a reviewable push/merge handoff. A clean
checkpoint may precede a release, but must not be described as release-ready
while these gates remain open. The maintainer will push; no push is requested
from the agent here.

## Ordered future work

1. **Rust:** profile the stabilized Python implementation; replace measured
   hotspots behind small, tested module interfaces. Preserve a Python reference
   and differential tests. Batch/stream evidence across the language seam so
   per-read object conversion does not erase improvements. Select libraries
   and deployment mechanisms when implementation starts, not from old guesses.
2. **Nextflow:** schedule the tested CLI stages with explicit inputs, outputs,
   resource requests, versions, and resume behavior. Keep biological decisions
   inside library modules. Validate chunk/contig gathering, deterministic
   ordering, and duplicate support handling before increasing parallelism.
3. **Long reads:** accept ONT and PacBio evidence, including full insertion
   spans, partial junctions, empty-locus spans, and multiple/nested insertions
   within one molecule. Recover inserted sequence when supported; report
   partial sequences, coverage, and uncertainty rather than implying every
   insertion is reconstructable. Genotyping needs validated long-read evidence
   rules; a spanning read alone is not a complete genotype determination.
4. **Multiple references and pangenomes:** identify homologous insertion loci
   across assemblies/haplotypes and report which contain the TE. Add graph
   representations when validated. Distinguish assemblies, graph paths,
   haplotypes, and biological subgenomes; they are not interchangeable labels.
   Subgenome attribution needs supplied annotations and sufficient evidence.

Novelty is always relative to the supplied, assessable references. Failure to
map, a missing assembly segment, or a graph path ambiguity must remain unknown,
not become proof of TE absence or an insertion being evolutionarily new.
Shared TE family sequence alone does not establish a shared insertion locus:
flank homology and insertion structure must support that claim.

## Data-model principles to preserve

These are design requirements for future work, not a request to replace today's
dataclasses. Following the codebase-design skill, add a new seam only when a
second implementation or data modality makes the variation concrete.

| Concept | Information to retain when the modality requires it |
|---|---|
| Evidence identity | Sample, read, mate/molecule identity, platform, segment coordinates, source alignment |
| Junction or spanning evidence | Query interval, orientation, alignment quality, alternate placements, evidence type; no mandatory short-read mate |
| Reference placement | Reference/assembly identifier, contig or path, coordinate convention, interval uncertainty, alternative placements |
| Insertion allele | Locus identity separate from TE-family label; optional full/partial inserted sequence and sequence provenance |
| Family evidence | Primary label plus alternate families and distinguishing support; no conversion of ties into independent votes |
| Genotype and occupancy | Sample allele evidence separately from present/absent/unknown state in each reference, haplotype, or subgenome |

Count independent molecules separately from alignment segments and junctions.
A long read can yield two junctions or several insertions; an alternate mapping
is not another molecule. Evidence weights/confidence must have documented
meaning; the current family-support fraction is not a calibrated probability.
Do not erase alternate alignments or reference provenance merely to fit a
single-reference output table.

Use versioned, streamable intermediate schemas when these become necessary.
Keep current linear-reference output available; richer future evidence should
not require encoding graph paths into chromosome names or compound family
labels. Test through module interfaces and maintain small real-data fixtures.

## Relationship to existing plans

`plans/FEATURES.md` already lists the requested sequence. This note updates the
release gate and evidence requirements. Its older claims about easy speedups,
one-read genotypes, aligner presets, graph path membership, and novelty rules
are hypotheses to revisit, not implementation specifications or validated
capabilities. No new external tool or schema format is selected here.

Next action: complete the focused parity stabilization described in
`2026-10-02-family-replay-results.md`; keep future feature work deferred.
