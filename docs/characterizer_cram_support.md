# Characterizer CRAM support

Date: 2026-06-19 14:56:29 PDT

## Purpose

Allow Step 7 characterization to consume sorted, indexed CRAM files in addition
to BAM files.

## Status

Implemented and covered by a unit test that creates a real CRAM file and checks
the resulting spanner count and genotype.

## Logic

- Detect CRAM input by a case-insensitive `.cram` suffix.
- Open CRAM with pysam mode `rc` and the supplied reference FASTA.
- Continue opening BAM with mode `rb`.
- Require `--genome-fasta` for CRAM so decoding uses an explicit,
  reproducible reference.
- Apply the same behavior to both characterizer APIs and both CLI front ends.

## Verification

```text
.pixi/envs/default/bin/pytest -q tests/characterize_test.py tests/main_test.py
8 passed

ruff check ...
All checks passed
```

Direct execution outside `pixi run` does not put minimap2 on `PATH`. A strict
Sphinx build rendered all pages but returned nonzero for existing warnings:
missing `docs/source/_static` and unavailable network-hosted intersphinx
inventories.

## Next steps

No further change is required for CRAM characterization. Callers must provide
the same reference FASTA used to encode the CRAM.
