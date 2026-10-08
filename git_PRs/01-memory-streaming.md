# Reduce BLAT conversion and trimming memory

BLAT processing retained complete FASTA/SAM output buffers, and trimming could
retain both mates' TE-hit dictionaries at once. This change streams selected
FASTA and converted SAM records, releases each completed mate before parsing
the next, and drops the sequence dictionary before sorting.

The existing `psl_to_sam` list-returning interface remains available; production
writing uses the shared iterator. Profiling tools and regression tests cover
streaming behavior, subprocess failures, and record lifetime.

On one representative riceTElib 30x sample, measured peak workflow RSS fell
from **34.45 GiB to 19.10 GiB**, with unchanged trimming artifacts and normalized
calls. The final run took 11:50:55 versus the original 11:41:05; runs used
different nodes, so this establishes a memory improvement, not a speedup.

Validation: the memory replay passed 279 tests and byte-for-byte artifact/call
checks. The assembled checkpoint subsequently passed **331 tests with zero
failures or skips** (SLURM job 29419250). The complete suite is validation of
the assembled stack, not an assertion that every intermediate PR ran 331 tests.

This PR changes resource handling; accuracy improvements are in subsequent
PRs. The representative resource measurement does not establish whole-panel
memory/runtime parity with RelocaTE2.
