"""Tests for single-sample pipeline orchestration and the run CLI."""

from __future__ import annotations

import os
from pathlib import Path

import pysam

from RelocaTE3.cli import main
from RelocaTE3.pipeline import run_sample
from RelocaTE3.ReadLibrary import ReadLibrary

DATA = Path(__file__).parent / "data"
R1 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"
R2 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"
TELIB = DATA / "TE_lib" / "mping.fa"
GENOME = DATA / "sim_genome" / "MSU7.Chr3_2M.fa"


def test_run_sample_produces_outputs(tmp_path: Path):
    """run_sample trims, re-aligns, calls insertions, and returns the GFF path.

    Pinned to minimap2 rather than the shipped blat/bwaaln defaults: this covers
    pipeline plumbing, not aligner choice, and must keep running in the default
    pixi environment, which cannot carry blat.
    """
    reads = ReadLibrary([str(R1), str(R2)], "HEG4")
    gff = run_sample(
        reads, str(TELIB), str(GENOME), tmp_path, threads=2,
        te_aligner="minimap2", genome_aligner="minimap2",
    )
    assert (tmp_path / "te_containing" / "HEG4.read_repeat_name.txt").exists()
    left = tmp_path / "flanking" / "HEG4.left.flankingReads.fq"
    assert left.exists() and os.path.getsize(left) > 0

    # genome BAM produced and indexed, containing tagged junction reads
    bam = tmp_path / "genome_aln" / "HEG4.genome.bam"
    assert bam.exists() and Path(str(bam) + ".bai").exists()
    with pysam.AlignmentFile(str(bam), "rb") as bf:
        names = [r.query_name for r in bf.fetch() if not r.is_unmapped]
    assert any(":5" in n or ":3" in n for n in names)

    # insertion GFF returned and non-empty
    assert gff.exists() and gff == tmp_path / "results" / "HEG4.all_nonref_insert.gff"
    assert os.path.getsize(gff) > 0


def test_run_cli_generates_flanking_reads(tmp_path: Path):
    """The installed `run` subcommand (map + trim) emits flanking FASTQs.

    This is the semantic equivalent of the old cli.py FASTQ `trim` command;
    both drive RelocaTE.identify_TE_reads. Exercised through RelocaTE3.cli.main
    to also cover the shim path.
    """
    rc = main(
        [
            "run",
            "-l",
            str(R1),
            "-r",
            str(R2),
            "-T",
            str(TELIB),
            "-n",
            "HEG4",
            "-o",
            str(tmp_path),
            "--threads",
            "2",
            # minimap2, not the blat default: see test_run_sample_produces_outputs
            "--te-aligner",
            "minimap2",
        ]
    )
    assert rc == 0
    assert (tmp_path / "flanking" / "HEG4.right.flankingReads.fq").exists()


def test_cli_main_handles_installed_subcommand():
    """RelocaTE3.cli.main is the canonical parser and handles installed subcommands (index-genome exists only on it)."""
    assert main(["index-genome", "--help"]) == 0


def test_run_sample_forwards_aligner_choices(monkeypatch, tmp_path: Path):
    """run_sample must let callers pick the aligners.

    The defaults are blat/bwaaln (matching RelocaTE2), but blat is not in the
    default pixi environment. Without a passthrough, every caller of run_sample
    -- including the acceptance gate -- would hard-require blat with no way to
    ask for anything else.
    """
    import RelocaTE3.pipeline as pipeline_mod

    seen: dict[str, str] = {}

    def fake_identify(self, reads, outdir, **kwargs):
        seen["te_aligner"] = kwargs.get("te_aligner")
        (Path(outdir) / "te_containing").mkdir(parents=True, exist_ok=True)
        return 0

    def fake_align(reads, genome, outdir, **kwargs):
        seen["genome_aligner"] = kwargs.get("genome_aligner")
        return Path(outdir) / "g.bam", None

    monkeypatch.setattr(pipeline_mod.RelocaTE, "identify_TE_reads", fake_identify)
    monkeypatch.setattr(pipeline_mod, "align_to_genome", fake_align)
    monkeypatch.setattr(pipeline_mod, "find_insertions", lambda *a, **k: [])

    reads = ReadLibrary([str(R1), str(R2)], "HEG4")
    run_sample(
        reads, str(TELIB), str(GENOME), tmp_path,
        te_aligner="minimap2", genome_aligner="bowtie2",
    )
    assert seen == {"te_aligner": "minimap2", "genome_aligner": "bowtie2"}


def test_run_sample_defaults_match_relocate2(monkeypatch, tmp_path: Path):
    """Unspecified, run_sample uses the RelocaTE2-matching pair."""
    import RelocaTE3.pipeline as pipeline_mod

    seen: dict[str, str] = {}
    monkeypatch.setattr(
        pipeline_mod.RelocaTE, "identify_TE_reads",
        lambda self, reads, outdir, **kw: (
            seen.__setitem__("te_aligner", kw.get("te_aligner")),
            (Path(outdir) / "te_containing").mkdir(parents=True, exist_ok=True),
            0,
        )[-1],
    )
    monkeypatch.setattr(
        pipeline_mod, "align_to_genome",
        lambda reads, genome, outdir, **kw: (
            seen.__setitem__("genome_aligner", kw.get("genome_aligner")),
            (Path(outdir) / "g.bam", None),
        )[-1],
    )
    monkeypatch.setattr(pipeline_mod, "find_insertions", lambda *a, **k: [])

    run_sample(ReadLibrary([str(R1), str(R2)], "HEG4"), str(TELIB), str(GENOME), tmp_path)
    assert seen == {"te_aligner": "blat", "genome_aligner": "bwaaln"}
