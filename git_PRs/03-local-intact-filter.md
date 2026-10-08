# Reject false junctions supported by locally intact original reads

A short trimmed flank can map beside a reference TE even when its original
read aligns nearly completely elsewhere nearby. The previous full-read filter
missed this evidence because it required the original alignment to cross the
reported breakpoint.

Count nearly complete original-read alignments within the existing same-contig
500-base candidate window, preserving mate identity. Reject when either the
existing spanning-read rule meets its 30% threshold on both sides or the local
intact-read rule does so. The two evidence types are not pooled across sides.
Both are collected in one bounded fetch; there is no whole-BAM read-name map.

The 72-sample local experiment removed **28 false positives with no lost or
gained truth detections**: three canonical riceTElib FPs and 25 diverged-panel
FPs. It retained all five real mPing detections lost by the earlier global
experiment. The global rule remains experimental and is not enabled in
production; no family-specific exception was added.

| Panel | R2 TP / FP | Local-filter R3 TP / FP |
|---|---:|---:|
| mPing | 2820 / 0 | 2826 / 0 |
| riceTElib | 2494 / 18 | 2513 / 24 |
| Diverged riceTElib | 5888 / 77 | 5919 / 104 |

These are paired BLAT/bwa aln results under the same TSD-interval scoring
policy and 10-bp tolerance. Canonical/zero-divergence conditions overlap.
The comparison reuses alignments; it is not a fresh full-pipeline resource run.

The integrated production rule passed **331 full-suite tests, zero failures
or skips** (job 29419250). Controls include actual SINE read evidence, remote
and other-contig alignments, clipped reads, wrong mates, and separate threshold
semantics. Depends on the family-evidence PR. Remaining FP and detection
differences prevent a universal RelocaTE2 noninferiority claim.
