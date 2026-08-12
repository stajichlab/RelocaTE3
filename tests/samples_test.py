"""Resolving a batch of samples from a sample sheet or a FASTQ directory.

RelocaTE2 took ``--fq_dir`` (a directory of paired FASTQs). RelocaTE3 required
one invocation per sample with explicit ``--left``/``--right``, so processing a
cohort meant hand-writing a loop. These helpers turn either input into an
explicit list of samples.

The sheet columns follow the schema already documented for the planned Nextflow
work (``plans/FEATURES.md``): ``sample_id, r1_fq, r2_fq`` plus optional
per-row ``te_library`` / ``reference_genome`` / ``repeatmasker`` overrides, so
sheets stay usable across both entry points.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from RelocaTE3.samples import discover_fastq_dir, read_sample_sheet


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("@r\nACGT\n+\nIIII\n")
    return path


# ---------------------------------------------------------------------------
# sample sheets
# ---------------------------------------------------------------------------


def test_reads_a_csv_sheet(tmp_path: Path):
    r1 = _touch(tmp_path / "a_R1.fq.gz")
    r2 = _touch(tmp_path / "a_R2.fq.gz")
    sheet = tmp_path / "samples.csv"
    sheet.write_text(f"sample_id,r1_fq,r2_fq\nA,{r1},{r2}\n")

    samples = read_sample_sheet(sheet)
    assert len(samples) == 1
    assert samples[0].name == "A"
    assert samples[0].r1 == str(r1)
    assert samples[0].r2 == str(r2)


def test_reads_a_tsv_sheet_and_skips_comments(tmp_path: Path):
    r1 = _touch(tmp_path / "b_1.fq")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "# cohort 2026-08\n"
        "sample_id\tr1_fq\tr2_fq\n"
        "\n"
        f"B\t{r1}\t\n"
    )
    samples = read_sample_sheet(sheet)
    assert len(samples) == 1
    assert samples[0].name == "B"
    assert samples[0].r2 is None, "empty r2 column means single-end"


def test_accepts_short_column_aliases(tmp_path: Path):
    r1 = _touch(tmp_path / "c_R1.fq")
    r2 = _touch(tmp_path / "c_R2.fq")
    sheet = tmp_path / "s.csv"
    sheet.write_text(f"sample,r1,r2\nC,{r1},{r2}\n")
    assert read_sample_sheet(sheet)[0].name == "C"


def test_per_row_overrides_are_carried(tmp_path: Path):
    r1 = _touch(tmp_path / "d_R1.fq")
    te = _touch(tmp_path / "other_TE.fa")
    sheet = tmp_path / "s.csv"
    sheet.write_text(f"sample_id,r1_fq,te_library\nD,{r1},{te}\n")
    sample = read_sample_sheet(sheet)[0]
    assert sample.te_library == str(te)
    assert sample.genome is None


def test_missing_required_column_is_a_clear_error(tmp_path: Path):
    sheet = tmp_path / "s.csv"
    sheet.write_text("sample_id,left\nA,x.fq\n")
    with pytest.raises(ValueError, match="r1_fq"):
        read_sample_sheet(sheet)


def test_missing_read_file_is_reported_with_the_sample(tmp_path: Path):
    sheet = tmp_path / "s.csv"
    sheet.write_text(f"sample_id,r1_fq\nA,{tmp_path / 'nope.fq'}\n")
    with pytest.raises(ValueError, match="A"):
        read_sample_sheet(sheet)


def test_duplicate_sample_ids_are_rejected(tmp_path: Path):
    r1 = _touch(tmp_path / "e_R1.fq")
    sheet = tmp_path / "s.csv"
    sheet.write_text(f"sample_id,r1_fq\nE,{r1}\nE,{r1}\n")
    with pytest.raises(ValueError, match="duplicate"):
        read_sample_sheet(sheet)


def test_empty_sheet_is_an_error(tmp_path: Path):
    sheet = tmp_path / "s.csv"
    sheet.write_text("sample_id,r1_fq\n")
    with pytest.raises(ValueError, match="no samples"):
        read_sample_sheet(sheet)


# ---------------------------------------------------------------------------
# --fq-dir discovery (the RelocaTE2 --fq_dir equivalent)
# ---------------------------------------------------------------------------


def test_discovers_r1_r2_pairs(tmp_path: Path):
    _touch(tmp_path / "sampleA_R1.fastq.gz")
    _touch(tmp_path / "sampleA_R2.fastq.gz")
    _touch(tmp_path / "sampleB_R1.fastq.gz")
    _touch(tmp_path / "sampleB_R2.fastq.gz")

    samples = discover_fastq_dir(tmp_path)
    assert [s.name for s in samples] == ["sampleA", "sampleB"], "sorted by name"
    assert all(s.r2 is not None for s in samples)
    assert samples[0].r1.endswith("sampleA_R1.fastq.gz")


def test_discovers_underscore_1_2_pairs(tmp_path: Path):
    _touch(tmp_path / "HEG4_1.fq.gz")
    _touch(tmp_path / "HEG4_2.fq.gz")
    samples = discover_fastq_dir(tmp_path)
    assert len(samples) == 1 and samples[0].name == "HEG4"
    assert samples[0].r2 is not None


def test_unpaired_file_becomes_a_single_end_sample(tmp_path: Path):
    _touch(tmp_path / "solo.fastq.gz")
    samples = discover_fastq_dir(tmp_path)
    assert len(samples) == 1
    assert samples[0].name == "solo" and samples[0].r2 is None


def test_orphaned_mate_is_an_error_not_a_silent_single_end(tmp_path: Path):
    """An R1 with no R2 is far more likely a mistake than a real SE sample."""
    _touch(tmp_path / "lonely_R1.fastq.gz")
    with pytest.raises(ValueError, match="lonely"):
        discover_fastq_dir(tmp_path)


def test_directory_without_fastqs_is_an_error(tmp_path: Path):
    (tmp_path / "notes.txt").write_text("hi")
    with pytest.raises(ValueError, match="no FASTQ"):
        discover_fastq_dir(tmp_path)
