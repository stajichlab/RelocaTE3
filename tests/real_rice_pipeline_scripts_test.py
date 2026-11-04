"""Real-rice validation through the Nextflow pipeline (M2 Task 5).

* ``validation/real_rice/make_samplesheet.py`` turns the real-rice config into
  the pipeline's sample sheet and prints the calling parameters, so the
  pipeline and ``run_relocate3.sh`` (run-all) read the same settings.
* ``validation/real_rice/run_pipeline.sh`` launches the pipeline; ``DRY_RUN=1``
  prints the Nextflow command instead of running it.
* ``validation/real_rice/compare_pipeline_to_rc2.sh`` compares the pipeline's
  per-sample output with the v0.2.0rc2 run-all output.

Everything here runs on a fake data tree: the example config is pointed at it
with ``RELOCATE3_VALIDATION_ROOT``, so the settings checked are the real ones.
"""

from __future__ import annotations

import csv
import io
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
HARNESS = ROOT / "validation" / "real_rice"
MAKE_SHEET = HARNESS / "make_samplesheet.py"
RUN_PIPELINE = HARNESS / "run_pipeline.sh"
COMPARE = HARNESS / "compare_pipeline_to_rc2.sh"
EXAMPLE_CONFIG = HARNESS / "config.example.toml"
SAMPLES = ["S_1", "S_2"]


def _fake_root(tmp_path: Path, samples=SAMPLES, index=".crai") -> Path:
    """A validation_root laid out like the real one, with empty files."""
    root = tmp_path / "real_rice"
    for rel in ("ref_genome/MSU_r7.fa", "repeatmasker/MSU_r7.fa.out"):
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        (root / rel).write_text("x\n")
    (root / "fq_input").mkdir()
    (root / "sample_file").mkdir()
    (root / "sample_file" / "samples.csv").write_text("Sample\n" + "".join(f"{s}\n" for s in samples))
    for s in samples:
        for r in ("R1", "R2"):
            (root / "fq_input" / f"{s}_trim_{r}.fastq.gz").write_text("")
        aln = root / "aln_input" / s
        aln.mkdir(parents=True)
        (aln / f"{s}.sorted.markdup.cram").write_text("")
        if index:
            (aln / f"{s}.sorted.markdup.cram{index}").write_text("")
    return root


def _env(root: Path, **extra: str) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items()
           if k not in ("DRY_RUN", "CONFIG", "RELOCATE3_VALIDATION_ROOT", "SLURM_JOB_ID")}
    env["RELOCATE3_VALIDATION_ROOT"] = str(root)
    env.update(extra)
    return env


def _sheet(root: Path, *args: str, config: Path = EXAMPLE_CONFIG) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(MAKE_SHEET), "--config", str(config), *args],
        cwd=ROOT, env=_env(root), capture_output=True, text=True,
    )


# ---------------------------------------------------------------- sample sheet

def test_sheet_has_absolute_paths_in_sample_list_order(tmp_path):
    root = _fake_root(tmp_path)
    out = tmp_path / "sheet.csv"
    proc = _sheet(root, "-o", str(out))
    assert proc.returncode == 0, proc.stderr
    rows = list(csv.DictReader(out.open()))
    assert [r["sample"] for r in rows] == SAMPLES
    assert list(rows[0]) == ["sample", "fastq_1", "fastq_2", "genotype_bam"]
    r = rows[1]
    assert r["fastq_1"] == str(root / "fq_input" / "S_2_trim_R1.fastq.gz")
    assert r["fastq_2"] == str(root / "fq_input" / "S_2_trim_R2.fastq.gz")
    assert r["genotype_bam"] == str(root / "aln_input" / "S_2" / "S_2.sorted.markdup.cram")
    for row in rows:
        for col in ("fastq_1", "fastq_2", "genotype_bam"):
            assert os.path.isabs(row[col])


def test_sheet_to_stdout(tmp_path):
    root = _fake_root(tmp_path)
    proc = _sheet(root, "-o", "-")
    assert proc.returncode == 0, proc.stderr
    rows = list(csv.DictReader(io.StringIO(proc.stdout)))
    assert [r["sample"] for r in rows] == SAMPLES


def test_missing_cram_index_is_a_clear_error(tmp_path):
    root = _fake_root(tmp_path)
    crai = root / "aln_input" / "S_2" / "S_2.sorted.markdup.cram.crai"
    crai.unlink()
    out = tmp_path / "sheet.csv"
    proc = _sheet(root, "-o", str(out))
    assert proc.returncode != 0
    assert "S_2" in proc.stderr and str(crai) in proc.stderr
    assert "index" in proc.stderr
    assert not out.exists()


def test_every_missing_input_is_reported(tmp_path):
    root = _fake_root(tmp_path)
    (root / "fq_input" / "S_1_trim_R2.fastq.gz").unlink()
    (root / "aln_input" / "S_2" / "S_2.sorted.markdup.cram").unlink()
    proc = _sheet(root, "-o", str(tmp_path / "sheet.csv"))
    assert proc.returncode != 0
    assert "S_1_trim_R2.fastq.gz" in proc.stderr
    assert "S_2.sorted.markdup.cram" in proc.stderr


@pytest.mark.parametrize("index", [".bai", ".csi"])
def test_bam_with_bai_or_csi_index_is_accepted(tmp_path, index):
    root = _fake_root(tmp_path, index=None)
    cfg_text = EXAMPLE_CONFIG.read_text().replace(".sorted.markdup.cram", ".sorted.markdup.bam")
    cfg = tmp_path / "validation" / "real_rice" / "config.toml"
    cfg.parent.mkdir(parents=True)
    cfg.write_text(cfg_text.replace('repository_root = "../.."', f'repository_root = "{ROOT}"'))
    for s in SAMPLES:
        aln = root / "aln_input" / s
        (aln / f"{s}.sorted.markdup.bam").write_text("")
        (aln / f"{s}.sorted.markdup.bam{index}").write_text("")
    proc = _sheet(root, "-o", "-", config=cfg)
    assert proc.returncode == 0, proc.stderr
    assert ".sorted.markdup.bam" in proc.stdout


def test_existing_sheet_is_kept_when_identical_and_refused_when_different(tmp_path):
    root = _fake_root(tmp_path)
    out = tmp_path / "sheet.csv"
    assert _sheet(root, "-o", str(out)).returncode == 0
    first = out.read_text()
    again = _sheet(root, "-o", str(out))
    assert again.returncode == 0, again.stderr
    assert out.read_text() == first
    out.write_text("sample,fastq_1\nother,x.fq\n")
    proc = _sheet(root, "-o", str(out))
    assert proc.returncode != 0
    assert "differs" in proc.stderr
    assert out.read_text() == "sample,fastq_1\nother,x.fq\n"
    assert _sheet(root, "-o", str(out), "--force").returncode == 0
    assert out.read_text() == first


# ------------------------------------------- run_relocate3.sh <-> pipeline

RUN_RELOCATE3 = HARNESS / "run_relocate3.sh"
# run-all options run_relocate3.sh passes that are not pipeline params: per
# sample (sample sheet columns), per run (--outdir) or task resources.
NOT_PARAMS = {"left", "right", "name", "outdir", "threads", "genotype_bam"}
# run-all argparse dest -> pipeline param (tests/nextflow_schema_test.py's
# RUN_ALL_DEST, inverted, plus the inputs and --genotype).
DEST_TO_PARAM = {
    "te_aligner": "te_aligner", "genome_aligner": "genome_aligner",
    "mismatch_allowance": "mismatch", "minimum_match_length": "min_match",
    "minimum_trimmed_length": "min_trimmed", "insert_size": "size",
    "min_mapq": "min_mapq", "tsd": "tsd", "require_both_junctions": "require_both_junctions",
    "target": "target", "te_name": "te_name", "repeatmasker": "repeatmasker",
    "te_opts": "te_opts", "genome_opts": "genome_opts", "distance": "distance",
    "excision": "excision", "bcftools": "bcftools",
    "te_library": "te_library", "genome_fasta": "genome", "genotype": "genotype",
}


def _run_all_dest(option: str) -> str:
    from RelocaTE3.cli import build_parser

    parser = build_parser()
    sub = next(a for a in parser._actions if a.__class__.__name__ == "_SubParsersAction")
    for action in sub.choices["run-all"]._actions:
        if option in action.option_strings:
            return action.dest
    raise AssertionError(f"run_relocate3.sh passes {option}, which run-all does not have")


def _run_all_call(script: str) -> list[tuple[str, str | None]]:
    """(option, value word) for every option of the ``relocaTE3 run-all`` call."""
    lines = script.splitlines()
    i = next(n for n, l in enumerate(lines) if l.strip().startswith("relocaTE3 run-all"))
    words = []
    for line in lines[i:]:
        words += shlex.split(line.strip().rstrip("\\"))
        if not line.rstrip().endswith("\\"):
            break
    expanded = []
    for w in words[2:]:
        m = re.fullmatch(r"\$\{(\w+)\[@\]\}", w)
        if m:  # an array of options, e.g. GENOTYPE_ARGS=(--genotype ...)
            decl = re.search(rf"^\s*{m.group(1)}=\((.+)\)\s*$", script, re.M)
            assert decl, f"no assignment of {m.group(1)} with options"
            expanded += shlex.split(decl.group(1))
        else:
            expanded.append(w)
    pairs: list[tuple[str, str | None]] = []
    for n, w in enumerate(expanded):
        if w.startswith("-"):
            nxt = expanded[n + 1] if n + 1 < len(expanded) else None
            pairs.append((w, nxt if nxt is not None and not nxt.startswith("-") else None))
    return pairs


def _config_source(script: str, word: str | None, option: str) -> tuple[str, str]:
    """The [section] key a value word comes from, via read_cfg's CFG array."""
    if word is None:
        # A switch: set when its guard variable is "1" (GENOTYPE_ARGS under
        # DO_CHARACTERIZE).
        m = re.search(r'if \[\[ "\$(\w+)" == "1" \]\]; then\n(?:.*\n){0,4}?\s*\w+=\(' + re.escape(option) + r"\b", script)
        assert m, f"{option}: no guard variable found"
        var = m.group(1)
    else:
        m = re.fullmatch(r"\$\{?(\w+)\}?", word)
        assert m, f"{option} {word}: not a plain variable"
        var = m.group(1)
    idx = re.search(rf'^{var}="\$\{{CFG\[(\d+)\]\}}"$', script, re.M)
    assert idx, f"{option}: {var} is not read from the config (CFG array)"
    body = script[script.index("read_cfg() {"):]
    prints = re.findall(r"^print\((.*)\)$", body, re.M)
    expr = prints[int(idx.group(1))]
    keys = re.findall(r"\['(\w+)'\]|\.get\('(\w+)'", expr)
    keys = [a or b for a, b in keys]
    assert len(keys) == 2, f"{option}: cannot read config key from {expr}"
    return keys[0], keys[1]


def _distinct_config(tmp_path: Path, root: Path) -> Path:
    """config.example.toml with a distinct value for every [relocate3] number,
    so a param mapped from the wrong key cannot match by accident."""
    text = EXAMPLE_CONFIG.read_text().replace('repository_root = "../.."', f'repository_root = "{ROOT}"')
    out, section, n = [], None, 0
    for line in text.splitlines():
        s = re.match(r"^\[(\w+)\]", line)
        if s:
            section = s.group(1)
        m = re.match(r"^(\w+)(\s*=\s*)(\d+)(.*)$", line)
        if section == "relocate3" and m and m.group(1) != "threads":
            n += 1
            line = f"{m.group(1)}{m.group(2)}{700 + n}{m.group(4)}"
        out.append(line)
    cfg = tmp_path / "config.toml"
    cfg.write_text("\n".join(out) + "\n")
    return cfg


def _drift(tmp_path: Path, script_path: Path) -> list[str]:
    """Settings run_relocate3.sh passes to run-all that --print-params does not
    map to the matching pipeline param with the value of the same config key,
    and printed params run_relocate3.sh does not pass."""
    sys.path.insert(0, str(HARNESS))
    from _config import load_config

    root = _fake_root(tmp_path)
    cfg_path = _distinct_config(tmp_path, root)
    proc = _sheet(root, "--print-params", config=cfg_path)
    assert proc.returncode == 0, proc.stderr
    printed = dict(line.split("\t", 1) for line in proc.stdout.splitlines())
    os.environ["RELOCATE3_VALIDATION_ROOT"] = str(root)
    try:
        cfg = load_config(cfg_path)
    finally:
        del os.environ["RELOCATE3_VALIDATION_ROOT"]

    script = script_path.read_text()
    problems, mapped = [], set()
    for option, word in _run_all_call(script):
        dest = _run_all_dest(option)
        if dest in NOT_PARAMS:
            continue
        param = DEST_TO_PARAM.get(dest)
        if param is None:
            problems.append(f"{option}: no pipeline param for run-all dest {dest}")
            continue
        mapped.add(param)
        section, key = _config_source(script, word, option)
        want = cfg[section][key]
        if isinstance(want, bool):
            want = "true" if want else "false"
        elif section == "paths":
            want = os.path.normpath(Path(cfg["_config_dir"]) / want)
        if param not in printed:
            problems.append(f"{option} ([{section}] {key}): --print-params has no {param}")
        elif printed[param] != str(want):
            problems.append(f"{option} ([{section}] {key} = {want}): --print-params {param} = {printed[param]}")
    for param in sorted(set(printed) - mapped):
        problems.append(f"--print-params {param}: run_relocate3.sh does not pass it")
    return problems


def test_every_run_all_setting_of_run_relocate3_is_a_pipeline_param(tmp_path):
    assert _drift(tmp_path, RUN_RELOCATE3) == []


def test_drift_guard_catches_a_new_run_all_flag(tmp_path):
    # The guard itself: a flag added to run_relocate3.sh but not to
    # make_samplesheet.py must be reported.
    copy = tmp_path / "run_relocate3.sh"
    copy.write_text(RUN_RELOCATE3.read_text().replace(
        '    --min-mapq "$MIN_MAPQ" \\\n',
        '    --min-mapq "$MIN_MAPQ" \\\n    --distance "$MISMATCH" \\\n',
    ))
    assert "--distance" in copy.read_text()
    problems = _drift(tmp_path, copy)
    assert any(p.startswith("--distance") for p in problems), problems


def test_print_samples(tmp_path):
    root = _fake_root(tmp_path)
    proc = _sheet(root, "--print-samples")
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.split() == SAMPLES


def test_default_config_prefers_config_toml(tmp_path):
    sys.path.insert(0, str(HARNESS))
    try:
        import make_samplesheet
    finally:
        sys.path.remove(str(HARNESS))
    (tmp_path / "config.example.toml").write_text("")
    assert make_samplesheet.default_config(tmp_path) == tmp_path / "config.example.toml"
    (tmp_path / "config.toml").write_text("")
    assert make_samplesheet.default_config(tmp_path) == tmp_path / "config.toml"


# ---------------------------------------------------------------- run_pipeline.sh

def _run_pipeline(root: Path, *args: str, dry: bool = True, **env: str) -> subprocess.CompletedProcess:
    extra = {"CONFIG": str(EXAMPLE_CONFIG), "PYTHON": sys.executable, **env}
    if dry:
        extra["DRY_RUN"] = "1"
    return subprocess.run(["bash", str(RUN_PIPELINE), *args], cwd=ROOT,
                          env=_env(root, **extra), capture_output=True, text=True)


def _nf_line(proc: subprocess.CompletedProcess) -> list[str]:
    lines = [l for l in proc.stdout.splitlines() if "nextflow run" in l]
    assert len(lines) == 1, proc.stdout
    return lines[0].split()


def _flag(words: list[str], flag: str) -> str:
    assert words.count(flag) == 1, f"{flag} must appear once: {words}"
    return words[words.index(flag) + 1]


def test_dry_run_prints_the_pipeline_command_with_the_config_settings(tmp_path):
    root = _fake_root(tmp_path)
    out = tmp_path / "nf_out"
    proc = _run_pipeline(root, "pr-7", str(out))
    assert proc.returncode == 0, proc.stderr
    words = _nf_line(proc)
    assert words[:8] == ["pixi", "run", "--manifest-path", str(ROOT / "pixi.toml"),
                         "--frozen", "-e", "nf", "nextflow"]
    assert words[8:10] == ["run", str(ROOT)]
    expected = {
        "-profile": "hpcc,apptainer",
        "--container_tag": "pr-7",
        "--input": str(out / "samplesheet.csv"),
        "--outdir": str(out),
        "-work-dir": str(tmp_path / "nf_out.nextflow" / "work"),
        "--te_aligner": "blat",
        "--genome_aligner": "bwaaln",
        "--mismatch": "2",
        "--min_match": "10",
        "--min_trimmed": "10",
        "--min_mapq": "0",
        "--size": "250",
        "--te_name": "mping",
        "--target": "ALL",
        "--tsd": "UNK",
        "--genotype": "true",
        "--genome": str(root / "ref_genome" / "MSU_r7.fa"),
        "--te_library": str(ROOT / "tests" / "data" / "TE_lib" / "mping.fa"),
        "--repeatmasker": str(root / "repeatmasker" / "MSU_r7.fa.out"),
    }
    for flag, value in expected.items():
        assert _flag(words, flag) == value, flag
    assert words.count("-resume") == 1
    # Launch directory (holds .nextflow/ history for -resume) is printed.
    assert re.search(rf"^launch dir +: {re.escape(str(tmp_path / 'nf_out.nextflow'))}$",
                     proc.stdout, re.M)
    # A dry run creates nothing.
    assert not out.exists() and not (tmp_path / "nf_out.nextflow").exists()


def test_dry_run_default_outdir_is_per_tag(tmp_path):
    root = _fake_root(tmp_path)
    proc = _run_pipeline(root, "sha-abc1234")
    assert proc.returncode == 0, proc.stderr
    words = _nf_line(proc)
    assert _flag(words, "--outdir") == str(ROOT / "results" / "nextflow-real-rice" / "sha-abc1234")
    assert not (ROOT / "results" / "nextflow-real-rice" / "sha-abc1234").exists()


def test_dry_run_prints_date_and_parameters(tmp_path):
    root = _fake_root(tmp_path)
    proc = _run_pipeline(root, "pr-7", str(tmp_path / "o"))
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.startswith("date")
    for key in ("container tag", "outdir", "work dir", "config", "apptainer cache", "size", "te_name"):
        assert key in proc.stdout, key


def test_dry_run_checks_the_inputs(tmp_path):
    root = _fake_root(tmp_path)
    (root / "aln_input" / "S_1" / "S_1.sorted.markdup.cram.crai").unlink()
    proc = _run_pipeline(root, "pr-7", str(tmp_path / "o"))
    assert proc.returncode != 0
    assert "S_1.sorted.markdup.cram.crai" in proc.stderr


@pytest.mark.parametrize("args", [(), ("a b",), ("pr-7", "out", "extra"), ("-x",)])
def test_bad_arguments_fail_with_usage(tmp_path, args):
    root = _fake_root(tmp_path)
    proc = _run_pipeline(root, *args)
    assert proc.returncode != 0
    assert "Usage" in proc.stderr


def test_help(tmp_path):
    proc = _run_pipeline(_fake_root(tmp_path), "--help", dry=False)
    assert proc.returncode == 0
    assert "Usage" in proc.stdout and "-resume" in proc.stdout


def test_outdir_from_another_container_tag_is_refused(tmp_path):
    root = _fake_root(tmp_path)
    out = tmp_path / "o"
    out.mkdir()
    (out / "container_tag.txt").write_text("pr-6\n")
    proc = _run_pipeline(root, "pr-7", str(out))
    assert proc.returncode != 0
    assert "pr-6" in proc.stderr


# ---------------------------------------------------------------- compare script

FILES = ["S.left.bam", "S.left.bam.bai", "existingTE.bed",
         "results/ALL.mping.all_nonref_insert.gff", "genome_aln/S.fullreads.genome.bam.bai"]


def _tree(base: Path, samples=SAMPLES) -> Path:
    for s in samples:
        for f in FILES:
            p = base / s / f.replace("S.", f"{s}.")
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(f"{s} {f}\n")
    return base


def _compare(nf: Path, rc2: Path, tmp_path: Path) -> subprocess.CompletedProcess:
    # samtools stub: `samtools view X` prints X (fake BAMs are text).
    stub = tmp_path / "bin"
    stub.mkdir(exist_ok=True)
    (stub / "samtools").write_text('#!/usr/bin/env bash\n[[ "$1" == view ]] && cat "$2"\n')
    (stub / "samtools").chmod(0o755)
    env = {k: v for k, v in os.environ.items() if k != "SLURM_JOB_ID"}
    env["PIXI_BIN"] = str(stub)
    return subprocess.run(["bash", str(COMPARE), str(nf), str(rc2)], cwd=ROOT,
                          env=env, capture_output=True, text=True)


def _with_sheet(nf: Path, samples=SAMPLES) -> Path:
    nf.mkdir(parents=True, exist_ok=True)
    (nf / "samplesheet.csv").write_text(
        "sample,fastq_1,fastq_2,genotype_bam\n" + "".join(f"{s},a,b,c\n" for s in samples))
    return nf


def test_compare_identical_trees_pass(tmp_path):
    rc2 = _tree(tmp_path / "rc2")
    nf = _tree(_with_sheet(tmp_path / "nf"))
    proc = _compare(nf, rc2, tmp_path)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "S_1 files=5 identical=3 differing=0 presence_only=2" in proc.stdout
    assert "S_2 files=5 identical=3 differing=0 presence_only=2" in proc.stdout
    assert proc.stdout.rstrip().endswith("RESULT: PASS")


def test_compare_names_a_differing_text_file(tmp_path):
    rc2 = _tree(tmp_path / "rc2")
    nf = _tree(_with_sheet(tmp_path / "nf"))
    (nf / "S_2" / "results" / "ALL.mping.all_nonref_insert.gff").write_text("changed\n")
    proc = _compare(nf, rc2, tmp_path)
    assert proc.returncode != 0
    assert "results/ALL.mping.all_nonref_insert.gff" in proc.stdout
    assert "S_2 files=5 identical=2 differing=1" in proc.stdout
    assert "S_1 files=5 identical=3 differing=0" in proc.stdout
    assert proc.stdout.rstrip().endswith("RESULT: FAIL")


def test_compare_names_a_differing_bam(tmp_path):
    rc2 = _tree(tmp_path / "rc2")
    nf = _tree(_with_sheet(tmp_path / "nf"))
    (nf / "S_1" / "S_1.left.bam").write_text("other records\n")
    proc = _compare(nf, rc2, tmp_path)
    assert proc.returncode != 0
    assert "S_1.left.bam" in proc.stdout
    assert proc.stdout.rstrip().endswith("RESULT: FAIL")


@pytest.mark.parametrize("missing", ["S_1.left.bam.bai", "existingTE.bed"])
def test_compare_missing_file_fails(tmp_path, missing):
    """A missing index counts too: indexes are compared by presence."""
    rc2 = _tree(tmp_path / "rc2")
    nf = _tree(_with_sheet(tmp_path / "nf"))
    (nf / "S_1" / missing).unlink()
    proc = _compare(nf, rc2, tmp_path)
    assert proc.returncode != 0
    assert f"only in rc2: ./{missing}" in proc.stdout
    assert proc.stdout.rstrip().endswith("RESULT: FAIL")


def test_compare_extra_file_fails(tmp_path):
    rc2 = _tree(tmp_path / "rc2")
    nf = _tree(_with_sheet(tmp_path / "nf"))
    (nf / "S_2" / "genome_aln" / "S_2.reads.genome.bam").write_text("x\n")
    proc = _compare(nf, rc2, tmp_path)
    assert proc.returncode != 0
    assert "only in pipeline: ./genome_aln/S_2.reads.genome.bam" in proc.stdout
    assert proc.stdout.rstrip().endswith("RESULT: FAIL")


def test_compare_missing_sample_dir_fails(tmp_path):
    rc2 = _tree(tmp_path / "rc2")
    nf = _tree(_with_sheet(tmp_path / "nf"), samples=["S_1"])
    proc = _compare(nf, rc2, tmp_path)
    assert proc.returncode != 0
    assert "S_2" in proc.stdout and "missing" in proc.stdout
    assert proc.stdout.rstrip().endswith("RESULT: FAIL")


def test_compare_rc2_sample_not_in_the_sheet_fails(tmp_path):
    rc2 = _tree(tmp_path / "rc2", samples=SAMPLES + ["S_3"])
    nf = _tree(_with_sheet(tmp_path / "nf"))
    proc = _compare(nf, rc2, tmp_path)
    assert proc.returncode != 0
    assert "S_3" in proc.stdout
    assert proc.stdout.rstrip().endswith("RESULT: FAIL")


def test_compare_needs_the_sample_sheet(tmp_path):
    rc2 = _tree(tmp_path / "rc2")
    nf = _tree(tmp_path / "nf")
    proc = _compare(nf, rc2, tmp_path)
    assert proc.returncode != 0
    assert "samplesheet.csv" in proc.stderr
