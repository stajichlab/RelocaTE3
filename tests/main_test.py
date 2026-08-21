from __future__ import annotations

import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Generator

import pytest

from RelocaTE3 import __version__, cli
from RelocaTE3.cli import main


def mockreturn(args):
    """Stand-in command handler that records its Namespace and succeeds."""
    mockreturn.args = args
    return 0


def test_main_version(capsys: pytest.CaptureFixture):
    assert main(["--version"]) == main(["-V"])
    captured: str = capsys.readouterr().out
    captured = captured.split("\n")
    assert captured[0] == captured[1] == __version__


def test_main_map(monkeypatch: Generator):
    monkeypatch.setattr(cli, "cmd_map", mockreturn)
    assert main(["map", "-l", "r1.fq", "-T", "te.fa", "-n", "HEG4"]) == 0
    assert mockreturn.args.name == "HEG4"
    assert mockreturn.args.left == "r1.fq"
    assert mockreturn.args.te_library == "te.fa"


def test_main_map_verbose(monkeypatch: Generator, caplog: pytest.LogCaptureFixture):
    monkeypatch.setattr(cli, "cmd_map", mockreturn)
    with caplog.at_level(logging.DEBUG):
        main(["map", "-l", "r1.fq", "-T", "te.fa", "-n", "HEG4", "-v"])
    assert any(record.levelname == "DEBUG" for record in caplog.records)
    assert "Debug mode enabled." in caplog.text


def test_main_characterize(monkeypatch: Generator):
    monkeypatch.setattr(cli, "cmd_characterize", mockreturn)
    assert main(["characterize", "-s", "sites.txt", "-b", "a.bam", "b.bam"]) == 0
    assert mockreturn.args.sites_file == "sites.txt"
    assert mockreturn.args.bam == ["a.bam", "b.bam"]


def test_module_entry_point_version():
    """`python -m RelocaTE3 --version` runs and prints a version."""
    proc = subprocess.run(
        [sys.executable, "-m", "RelocaTE3", "--version"],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0
    assert proc.stdout.strip()


@pytest.mark.skipif(
    shutil.which("relocaTE3") is None,
    reason="console script not installed on PATH",
)
def test_console_script_version():
    """The installed `relocaTE3 --version` console script runs."""
    proc = subprocess.run(["relocaTE3", "--version"], capture_output=True, text=True)
    assert proc.returncode == 0
    assert proc.stdout.strip()


def test_build_parser_registers_subcommand():
    """build_parser() returns a parser whose subcommands dispatch to cmd_*."""
    from RelocaTE3.cli import build_parser, cmd_index_genome

    parsed = build_parser().parse_args(["index-genome", "-g", "ref.fa"])
    assert parsed.func is cmd_index_genome


def test_te_and_genome_opts_are_parsed_and_shlex_split():
    """--te-opts/--genome-opts reach the aligner backends as argument lists.

    The benchmark drives RelocaTE3 through this CLI, so the per-stage aligner
    option passthrough is only usable if these flags exist and split correctly.

    Values start with "-" (e.g. -minIdentity=80), which argparse would otherwise
    read as another flag, so the "--opt=value" form is required. Callers/adapters
    must use that form.
    """
    from RelocaTE3.cli import build_parser, split_aligner_opts

    parser = build_parser()

    te = parser.parse_args(
        ["run", "-l", "r1.fq", "-T", "te.fa", "-n", "s", "-o", "out",
         "--te-aligner", "blat", "--te-opts=-minIdentity=80"]
    )
    assert split_aligner_opts(te.te_opts) == ["-minIdentity=80"]

    # multi-token values must split on whitespace, not be passed as one arg
    mapped = parser.parse_args(
        ["map", "-l", "r1.fq", "-T", "te.fa", "-n", "s", "-o", "out",
         "--te-opts=-B 4"]
    )
    assert split_aligner_opts(mapped.te_opts) == ["-B", "4"]

    genome = parser.parse_args(
        ["align-genome", "-g", "ref.fa", "-f", "f.fq", "-n", "s", "-o", "out",
         "--genome-aligner", "bwaaln", "--genome-opts=-n 0.10"]
    )
    assert split_aligner_opts(genome.genome_opts) == ["-n", "0.10"]

    # absent flags default to no extra options
    assert split_aligner_opts(None) == []


# ---------------------------------------------------------------------------
# find-reference (steps 0/6): reference/shared insertion CALLS
# ---------------------------------------------------------------------------


def test_main_find_reference_parses_args(monkeypatch: Generator):
    """'find-reference' is registered and binds its arguments."""
    monkeypatch.setattr(cli, "cmd_find_reference", mockreturn)
    assert (
        main(
            [
                "find-reference",
                "-b", "g.bam",
                "-R", "read_repeat.txt",
                "--repeatmasker", "rm.out",
                "-n", "HEG4",
                "-o", "out",
            ]
        )
        == 0
    )
    assert mockreturn.args.bam == "g.bam"
    assert mockreturn.args.read_repeat == "read_repeat.txt"
    assert mockreturn.args.repeatmasker == "rm.out"
    assert mockreturn.args.name == "HEG4"


def test_find_reference_writes_ref_insert_outputs(tmp_path):
    """RelocaTE2 parity: emit all_ref_insert.{gff,txt} plus existingTE.bed.

    RelocaTE2 reports TEs present in the reference AND confirmed in the sample
    (all_ref_insert.*). RelocaTE3 could only do this from the library
    (pipeline.run_sample); without a CLI command a user migrating from
    RelocaTE2 silently loses that half of the result set.
    """
    import pysam

    # one intact reference TE, with junction reads at both of its boundaries
    rm = tmp_path / "rm.out"
    rm.write_text(
        "   SW  perc perc perc  query   position in query        matching  repeat\n"
        "score  div. del. ins.  sequence begin end (left)  repeat  class/family\n"
        "\n"
        "  500  1.0  0.0  0.0  Chr3  1000  1500  (100) +  mPing  DNA/MITE  1 500 (0)  1\n"
    )

    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "Chr3", "LN": 2_000_000}]}
    bam_path = tmp_path / "g.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for name, pos in (("a/1:end:5", 900), ("b/1:start:5", 1500)):
            rec = pysam.AlignedSegment(bam.header)
            rec.query_name = name
            rec.query_sequence = "A" * 100
            rec.flag = 0
            rec.reference_id = 0
            rec.reference_start = pos
            rec.mapping_quality = 60
            rec.cigartuples = [(0, 100)]
            bam.write(rec)
    pysam.index(str(bam_path))

    read_repeat = tmp_path / "read_repeat.txt"
    read_repeat.write_text("a/1:end:5\tmPing\t+\nb/1:start:5\tmPing\t+\n")

    outdir = tmp_path / "out"
    rc = main(
        [
            "find-reference",
            "-b", str(bam_path),
            "-R", str(read_repeat),
            "--repeatmasker", str(rm),
            "-n", "HEG4",
            "-o", str(outdir),
        ]
    )
    assert rc == 0

    gff = outdir / "results" / "HEG4.all_ref_insert.gff"
    txt = outdir / "results" / "HEG4.all_ref_insert.txt"
    assert gff.exists(), "reference/shared insertions GFF not written"
    assert txt.exists(), "reference/shared insertions table not written"
    assert (outdir / "existingTE.bed").exists()

    body = gff.read_text()
    assert "mPing" in body
    assert "Shared" in body


# ---------------------------------------------------------------------------
# run-all: one-command pipeline (must equal the staged workflow)
# ---------------------------------------------------------------------------


def test_main_run_all_parses_args(monkeypatch: Generator):
    """'run-all' is registered and binds the staged workflow's arguments."""
    monkeypatch.setattr(cli, "cmd_run_all", mockreturn)
    assert (
        main(
            [
                "run-all",
                "-1", "r1.fq", "-2", "r2.fq",
                "-T", "te.fa",
                "-g", "ref.fa",
                "-n", "HEG4",
                "-o", "out",
                "--threads", "4",
                "--te-aligner", "blat",
                "--genome-aligner", "bwaaln",
                "--tsd", "UNK",
                "--mismatch", "2",
                "--min-mapq", "1",
                "--require-both-junctions",
                "--repeatmasker", "rm.out",
                "--genotype",
            ]
        )
        == 0
    )
    a = mockreturn.args
    assert a.left == "r1.fq" and a.right == "r2.fq"
    assert a.te_library == "te.fa" and a.genome_fasta == "ref.fa"
    assert a.te_aligner == "blat" and a.genome_aligner == "bwaaln"
    assert a.tsd == "UNK" and a.require_both_junctions is True
    assert a.repeatmasker == "rm.out" and a.genotype is True


def _fake_stage_outputs(outdir, sample: str, aligner: str = "bwaaln"):
    """Lay down the intermediates the real steps would have produced.

    run-all checks between steps that each stage actually emitted something, so
    with the handlers stubbed out those artifacts have to be faked.
    """
    (outdir / "te_containing").mkdir(parents=True, exist_ok=True)
    (outdir / "te_containing" / f"{sample}.read_repeat_name.txt").write_text(
        "r1\tmPing\t+\n"
    )
    (outdir / "flanking").mkdir(parents=True, exist_ok=True)
    (outdir / "flanking" / f"{sample}.left.flankingReads.fq").write_text("@r1\nA\n+\nI\n")
    (outdir / f"{sample}.repeat.{aligner}.sorted.bam").write_bytes(b"")
    (outdir / "results").mkdir(parents=True, exist_ok=True)


def test_run_all_delegates_to_the_staged_commands(monkeypatch, tmp_path):
    """run-all must BE the staged workflow, not a parallel implementation.

    RelocaTE3 has two distinct insertion-calling code paths -- the CLI's
    InsertionFinder (what the benchmark validated) and pipeline.run_sample's
    module-level find_insertions. They are not equivalent, so a one-command mode
    built on the latter would silently produce different, unvalidated numbers.
    This pins run-all to dispatching the same cmd_* handlers the staged
    workflow uses, with the same arguments.
    """
    calls: list[tuple[str, object]] = []

    def rec(tag):
        def _inner(args):
            calls.append((tag, args))
            return 0
        return _inner

    monkeypatch.setattr(cli, "cmd_index_genome", rec("index-genome"))
    monkeypatch.setattr(cli, "cmd_run", rec("run"))
    monkeypatch.setattr(cli, "cmd_align_genome", rec("align-genome"))
    monkeypatch.setattr(cli, "cmd_find_insertions", rec("find-insertions"))
    monkeypatch.setattr(cli, "cmd_find_reference", rec("find-reference"))
    monkeypatch.setattr(cli, "cmd_characterize", rec("characterize"))

    _fake_stage_outputs(tmp_path, "HEG4")

    rc = main(
        [
            "run-all",
            "-1", "r1.fq", "-2", "r2.fq",
            "-T", "te.fa",
            "-g", "ref.fa",
            "-n", "HEG4",
            "-o", str(tmp_path),
            "--threads", "4",
            "--te-aligner", "blat",
            "--genome-aligner", "bwaaln",
            "--tsd", "UNK",
            "--te-name", "riceTElib",
            "--mismatch", "2",
            "--min-mapq", "1",
            "--require-both-junctions",
            "--repeatmasker", "rm.out",
        ]
    )
    assert rc == 0

    order = [tag for tag, _ in calls]
    assert order == [
        "index-genome",
        "run",
        "align-genome",
        "find-insertions",
        "find-reference",
    ], order

    by = dict(calls)
    # step 2-3 carries the TE-side settings
    assert by["run"].te_aligner == "blat"
    assert by["run"].mismatch_allowance == 2
    # step 4 carries the genome-side settings AND the mate-anchoring reads,
    # which is what closed most of the TSD gap vs RelocaTE2
    assert by["align-genome"].genome_aligner == "bwaaln"
    assert by["align-genome"].left == "r1.fq" and by["align-genome"].right == "r2.fq"
    # step 5 carries the calling thresholds validated in the benchmark
    fi = by["find-insertions"]
    assert fi.tsd == "UNK" and fi.te_name == "riceTElib"
    # min_mapq is passed explicitly as 1 above; this asserts pass-through, not
    # the default (which is 0 -- see relocate2_defaults_test.py).
    assert fi.require_both_junctions is True and fi.min_mapq == 1
    assert fi.reference_ins == "rm.out"


def test_run_all_skips_optional_stages_when_not_requested(monkeypatch, tmp_path):
    """No --repeatmasker => no find-reference; no --genotype => no characterize."""
    calls: list[str] = []

    def rec(tag):
        def _inner(args):
            calls.append(tag)
            return 0
        return _inner

    for name, tag in [
        ("cmd_index_genome", "index-genome"),
        ("cmd_run", "run"),
        ("cmd_align_genome", "align-genome"),
        ("cmd_find_insertions", "find-insertions"),
        ("cmd_find_reference", "find-reference"),
        ("cmd_characterize", "characterize"),
    ]:
        monkeypatch.setattr(cli, name, rec(tag))

    _fake_stage_outputs(tmp_path, "S", aligner="minimap")

    rc = main(
        ["run-all", "-1", "r1.fq", "-T", "te.fa", "-g", "ref.fa",
         "-n", "S", "-o", str(tmp_path)]
    )
    assert rc == 0
    assert "find-reference" not in calls
    assert "characterize" not in calls
    assert calls == ["index-genome", "run", "align-genome", "find-insertions"]


# ---------------------------------------------------------------------------
# run-batch: many samples from a sheet or a directory
# ---------------------------------------------------------------------------


def _sheet(tmp_path, n=2):
    """A sample sheet with n samples whose read files exist."""
    rows = ["sample_id,r1_fq,r2_fq"]
    for i in range(n):
        r1 = tmp_path / f"S{i}_R1.fq"
        r2 = tmp_path / f"S{i}_R2.fq"
        for f in (r1, r2):
            f.write_text("@r\nACGT\n+\nIIII\n")
        rows.append(f"S{i},{r1},{r2}")
    sheet = tmp_path / "samples.csv"
    sheet.write_text("\n".join(rows) + "\n")
    return sheet


def test_main_run_batch_parses_args(monkeypatch: Generator, tmp_path):
    monkeypatch.setattr(cli, "cmd_run_batch", mockreturn)
    sheet = _sheet(tmp_path)
    assert main(["run-batch", "--samples", str(sheet), "-T", "te.fa",
                 "-g", "ref.fa", "-o", str(tmp_path / "out")]) == 0
    assert mockreturn.args.samples == str(sheet)
    assert mockreturn.args.te_library == "te.fa"


def test_run_batch_runs_each_sample_in_its_own_subdir(monkeypatch, tmp_path):
    """One run-all per sample, isolated outputs, shared settings."""
    seen = []

    def fake_run_all(args):
        seen.append(args)
        Path(args.outdir, "results").mkdir(parents=True, exist_ok=True)
        Path(args.outdir, "results", "ALL.repeat.all_nonref_insert.txt").write_text("x\n")
        return 0

    monkeypatch.setattr(cli, "cmd_run_all", fake_run_all)
    sheet = _sheet(tmp_path, n=3)
    out = tmp_path / "out"

    rc = main(["run-batch", "--samples", str(sheet), "-T", "te.fa", "-g", "ref.fa",
               "-o", str(out), "--threads", "4", "--te-aligner", "blat",
               "--require-both-junctions"])
    assert rc == 0
    assert [a.name for a in seen] == ["S0", "S1", "S2"]
    assert [Path(a.outdir).name for a in seen] == ["S0", "S1", "S2"]
    # shared settings reach every sample
    assert all(a.te_aligner == "blat" for a in seen)
    assert all(a.require_both_junctions is True for a in seen)
    assert all(a.threads == 4 for a in seen)


def test_run_batch_per_row_overrides_beat_the_global_flags(monkeypatch, tmp_path):
    seen = []
    monkeypatch.setattr(cli, "cmd_run_all", lambda a: (seen.append(a), 0)[1])

    te = tmp_path / "special_TE.fa"
    te.write_text(">te\nACGT\n")
    r1 = tmp_path / "X_R1.fq"
    r1.write_text("@r\nACGT\n+\nIIII\n")
    sheet = tmp_path / "s.csv"
    sheet.write_text(f"sample_id,r1_fq,te_library\nX,{r1},{te}\n")

    main(["run-batch", "--samples", str(sheet), "-T", "global_TE.fa",
          "-g", "ref.fa", "-o", str(tmp_path / "out")])
    assert seen[0].te_library == str(te)


def test_run_batch_skips_completed_samples_and_force_reruns(monkeypatch, tmp_path):
    """Resumable: a finished sample is not redone unless --force."""
    calls = []
    monkeypatch.setattr(cli, "cmd_run_all", lambda a: (calls.append(a.name), 0)[1])

    sheet = _sheet(tmp_path, n=2)
    out = tmp_path / "out"
    done = out / "S0" / "results"
    done.mkdir(parents=True)
    (done / "ALL.repeat.all_nonref_insert.txt").write_text("already done\n")

    main(["run-batch", "--samples", str(sheet), "-T", "te.fa", "-g", "ref.fa",
          "-o", str(out)])
    assert calls == ["S1"], "completed sample should be skipped"

    calls.clear()
    main(["run-batch", "--samples", str(sheet), "-T", "te.fa", "-g", "ref.fa",
          "-o", str(out), "--force"])
    assert calls == ["S0", "S1"], "--force should redo everything"


def test_run_batch_reports_failures_and_exits_nonzero(monkeypatch, tmp_path):
    """A failing sample must not be silently swallowed."""
    attempted = []

    def flaky(args):
        attempted.append(args.name)
        if args.name == "S0":
            raise RuntimeError("boom")
        return 0

    monkeypatch.setattr(cli, "cmd_run_all", flaky)
    sheet = _sheet(tmp_path, n=3)
    rc = main(["run-batch", "--samples", str(sheet), "-T", "te.fa", "-g", "ref.fa",
               "-o", str(tmp_path / "out"), "--keep-going"])
    assert rc != 0, "batch with a failed sample must exit non-zero"
    assert attempted == ["S0", "S1", "S2"], (
        "--keep-going must carry on past the failure, not stop at it"
    )


def test_run_batch_stops_at_first_failure_by_default(monkeypatch, tmp_path):
    attempted = []

    def flaky(args):
        attempted.append(args.name)
        raise RuntimeError("boom")

    monkeypatch.setattr(cli, "cmd_run_all", flaky)
    sheet = _sheet(tmp_path, n=3)
    rc = main(["run-batch", "--samples", str(sheet), "-T", "te.fa", "-g", "ref.fa",
               "-o", str(tmp_path / "out")])
    assert rc != 0
    assert attempted == ["S0"], "should stop rather than burn the whole cohort"


def test_run_batch_accepts_a_fastq_directory(monkeypatch, tmp_path):
    """RelocaTE2's --fq_dir equivalent."""
    seen = []
    monkeypatch.setattr(cli, "cmd_run_all", lambda a: (seen.append(a.name), 0)[1])

    fq = tmp_path / "fq"
    fq.mkdir()
    for base in ("alpha", "beta"):
        for mate in ("R1", "R2"):
            (fq / f"{base}_{mate}.fastq.gz").write_text("@r\nACGT\n+\nIIII\n")

    rc = main(["run-batch", "--fq-dir", str(fq), "-T", "te.fa", "-g", "ref.fa",
               "-o", str(tmp_path / "out")])
    assert rc == 0
    assert seen == ["alpha", "beta"]


def test_run_batch_requires_exactly_one_input_source(tmp_path):
    rc = main(["run-batch", "-T", "te.fa", "-g", "ref.fa", "-o", str(tmp_path)])
    assert rc != 0, "neither --samples nor --fq-dir should be an error"
