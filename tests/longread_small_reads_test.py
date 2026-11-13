"""Tests for the pure parts of the small-panel read driver."""

from __future__ import annotations

import gzip
import hashlib
import sys
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "validation" / "longread_small"))

import simulate_reads as sr  # noqa: E402


def test_jobs_are_platform_by_coverage():
    jobs = sr.jobs(["ont-hq", "hifi"], [5, 20], replicates=1, seed=916)
    assert [(j.platform, j.sample) for j in jobs] == [
        ("ont-hq", "cov5x_rep1"),
        ("ont-hq", "cov20x_rep1"),
        ("hifi", "cov5x_rep1"),
        ("hifi", "cov20x_rep1"),
    ]
    assert len({j.seed for j in jobs} | {j.control_seed for j in jobs}) == 8


def test_manifest_matches_full_panel_columns(tmp_path):
    jobs = sr.jobs(["ont-hq"], [5], replicates=1, seed=916)
    sr.write_manifest(jobs, tmp_path / "panel_manifest.tsv")
    lines = (tmp_path / "panel_manifest.tsv").read_text().splitlines()
    assert lines[0].split("\t") == [
        "sample",
        "platform",
        "coverage",
        "replicate",
        "reads",
        "control_reads",
    ]
    assert lines[1].split("\t") == [
        "cov5x_rep1",
        "ont-hq",
        "5",
        "1",
        "reads/ont-hq/cov5x_rep1/reads.fastq.gz",
        "reads/ont-hq/cov5x_rep1_reference_control/reads.fastq.gz",
    ]


def test_needs_run_states(tmp_path):
    target = tmp_path / "reads" / "ont-hq" / "cov5x_rep1"
    assert sr.needs_run(target) is True  # absent
    target.mkdir(parents=True)
    assert sr.needs_run(target) is True  # empty directory
    (target / "reads.fastq.gz").write_bytes(b"x")
    with pytest.raises(RuntimeError, match="incomplete"):
        sr.needs_run(target)  # output without completion marker
    (target / ".complete").touch()
    assert sr.needs_run(target) is False


def test_incomplete_error_names_what_to_move_aside(tmp_path):
    target = tmp_path / "cov5x_rep1"
    target.mkdir()
    (target / "reads.fastq.gz.partial").write_bytes(b"x")
    with pytest.raises(RuntimeError) as caught:
        sr.needs_run(target)
    message = str(caught.value)
    assert str(target) in message
    assert "reads.fastq.gz.partial" in message
    assert f"mv {target} {target}.failed" in message


def test_validate_accepts_good_inputs(tmp_path):
    genome = tmp_path / "g.fa"
    genome.write_text(">Chr3\nACGT\n")
    sr.validate_inputs(["ont-hq", "hifi"], [2, 5], 1, [genome, genome])


@pytest.mark.parametrize(
    ("platforms", "coverages", "replicates", "match"),
    [
        (["ont-hq", "clr"], [5], 1, "unknown platform 'clr'"),
        (["ont-hq", "ont-hq"], [5], 1, "repeated platform 'ont-hq'"),
        (["hifi"], [5, 5], 1, "repeated coverage 5"),
        (["hifi"], [0], 1, "coverage must be >= 1, got 0"),
        (["hifi"], [5], 0, "replicates must be >= 1, got 0"),
    ],
)
def test_validate_rejects_bad_inputs(tmp_path, platforms, coverages, replicates, match):
    genome = tmp_path / "g.fa"
    genome.write_text(">Chr3\nACGT\n")
    with pytest.raises(ValueError, match=match):
        sr.validate_inputs(platforms, coverages, replicates, [genome])


def test_validate_rejects_missing_genome(tmp_path):
    missing = tmp_path / "absent.fa"
    with pytest.raises(FileNotFoundError, match=str(missing)):
        sr.validate_inputs(["hifi"], [5], 1, [missing])


def test_main_validates_before_simulating(tmp_path):
    genome = tmp_path / "g.fa"
    genome.write_text(">Chr3\nACGT\n")
    argv = [
        "--panel",
        str(tmp_path / "panel"),
        "--sample-genome",
        str(genome),
        "--reference",
        str(genome),
        "--platform",
        "clr",
    ]
    with pytest.raises(ValueError, match="unknown platform"):
        sr.main(argv)
    assert not (tmp_path / "panel").exists()


def test_pooled_gzip_is_byte_deterministic(tmp_path):
    outputs = []
    for name in ("a.fastq.gz", "b.fastq.gz"):
        path = tmp_path / name
        with sr.open_pooled(path) as handle:
            handle.write("@r1\nACGT\n+\nIIII\n")
        outputs.append(path.read_bytes())
        time.sleep(1.1)  # a gzip mtime would differ between the two files
    assert outputs[0] == outputs[1]
    assert gzip.decompress(outputs[0]) == b"@r1\nACGT\n+\nIIII\n"


def _fake_simulate(calls):
    """Stand-in for `simulate`: writes the files the real one writes."""

    def fake(platform_key, genome, coverage, seed, directory, name_prefix, *_):
        calls.append(directory)
        directory.mkdir(parents=True, exist_ok=True)
        with gzip.open(directory / "reads.fastq.gz", "wt") as handle:
            handle.write(f"@{name_prefix}:S1_1\nACGT\n+\nIIII\n")
        sr.write_simulation_record(
            directory, platform_key, coverage, seed, genome, 1, None
        )
        (directory / ".complete").touch()
        return 1

    return fake


def _argv(panel, sample, reference, *extra):
    return [
        "--panel",
        str(panel),
        "--sample-genome",
        str(sample),
        "--reference",
        str(reference),
        "--platform",
        "ont-hq",
        "--coverage",
        "5",
        *extra,
    ]


@pytest.fixture
def genomes(tmp_path):
    sample = tmp_path / "sample.fa"
    sample.write_text(">Chr3\nACGTACGT\n")
    reference = tmp_path / "reference.fa"
    reference.write_text(">Chr3\nACGT\n")
    return sample, reference


def test_simulation_record_holds_genome_identity(tmp_path, genomes):
    sample, _ = genomes
    sr.write_simulation_record(tmp_path, "hifi", 5, 916, sample, 12, 0.95)
    rows = (tmp_path / "simulation.tsv").read_text().splitlines()
    record = dict(zip(rows[0].split("\t"), rows[1].split("\t"), strict=True))
    assert record == {
        "platform": "hifi",
        "coverage": "5",
        "seed": "916",
        "genome": str(sample),
        "genome_bytes": str(sample.stat().st_size),
        "genome_sha256": hashlib.sha256(sample.read_bytes()).hexdigest(),
        "reads": "12",
        "ccs_yield": "0.9500",
    }


def test_matching_rerun_skips_and_writes_manifest(tmp_path, genomes, monkeypatch):
    calls = []
    monkeypatch.setattr(sr, "simulate", _fake_simulate(calls))
    panel = tmp_path / "panel"
    assert sr.main(_argv(panel, *genomes)) == 0
    assert len(calls) == 2
    (panel / "panel_manifest.tsv").unlink()
    assert sr.main(_argv(panel, *genomes)) == 0
    assert len(calls) == 2  # nothing simulated again
    assert (panel / "panel_manifest.tsv").is_file()


def test_rerun_with_changed_seed_is_refused(tmp_path, genomes, monkeypatch):
    monkeypatch.setattr(sr, "simulate", _fake_simulate([]))
    panel = tmp_path / "panel"
    sr.main(_argv(panel, *genomes))
    with pytest.raises(RuntimeError, match=r"cov5x_rep1.*seed"):
        sr.main(_argv(panel, *genomes, "--seed", "1000"))


def test_rerun_with_changed_genome_is_refused(tmp_path, genomes, monkeypatch):
    monkeypatch.setattr(sr, "simulate", _fake_simulate([]))
    panel = tmp_path / "panel"
    sr.main(_argv(panel, *genomes))
    sample, _ = genomes
    sample.write_text(">Chr3\nTTTTACGT\n")
    with pytest.raises(RuntimeError, match=r"cov5x_rep1.*genome_sha256"):
        sr.main(_argv(panel, *genomes))
    before = (panel / "reads/ont-hq/cov5x_rep1/reads.fastq.gz").read_bytes()
    assert gzip.decompress(before).startswith(b"@cov5x_rep1:sample")  # untouched


@pytest.mark.parametrize("content", [None, "", "platform\tseed\nhifi\t1\n"])
def test_complete_dir_without_usable_record_is_refused(tmp_path, genomes, content):
    target = tmp_path / "cov5x_rep1"
    target.mkdir()
    (target / ".complete").touch()
    if content is not None:
        (target / "simulation.tsv").write_text(content)
    with pytest.raises(RuntimeError, match=str(target)):
        sr.check_complete(target, "ont-hq", 5, 916, genomes[0])
