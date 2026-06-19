# Test TE library path migration

Date: 2026-06-19 15:02:09 PDT

## Purpose

Update tests and real-rice validation configuration after moving test TE FASTA
files into `tests/data/TE_lib/`.

## Changes

- Unit and pipeline tests use `tests/data/TE_lib/mping.fa`.
- The acceptance test uses `tests/data/TE_lib/RiceTE.fa`.
- Real-rice validation uses the repository-owned
  `tests/data/TE_lib/mping.fa` rather than the external validation-data copy.
- The ignored local `validation/real_rice/config.toml` was updated along with
  the tracked `config.example.toml`.

## Verification

```text
pytest -q
35 passed

ruff check <changed test and validation Python files>
All checks passed
```

Both the example and local real-rice configurations resolve the mPing FASTA to
an existing file.
