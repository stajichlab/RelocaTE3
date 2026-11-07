"""Real-rice validation modes (M5 Task 2).

* ``validation/real_rice/run_validation.sh --mode pipeline`` (the default) runs
  the Nextflow pipeline (``run_pipeline.sh``, pinned ``--container-tag``) and
  reports on its outdir; ``--mode run-all`` is the run-all SLURM array as before.
  Both feed the same normalize + compare commands; only the RelocaTE3 outdir
  they read differs. ``DRY_RUN=1`` prints the commands instead of running them.
* ``validation/real_rice/_config.py`` takes the RelocaTE3 outdir and the report
  dir from ``RELOCATE3_VALIDATION_OUTDIR`` / ``RELOCATE3_VALIDATION_REPORT_DIR``,
  so the report can point at a pipeline outdir without editing the config.
* ``validation/release/submit_real_rice.sh`` (release check) uses pipeline mode
  by default with a required ``--container-tag``; ``--mode run-all`` remains.

Fake data trees throughout (``RELOCATE3_VALIDATION_ROOT``), as in
``tests/real_rice_pipeline_scripts_test.py``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from tests.real_rice_pipeline_scripts_test import SAMPLES, _fake_root

ROOT = Path(__file__).resolve().parent.parent
HARNESS = ROOT / "validation" / "real_rice"
RUN_VALIDATION = HARNESS / "run_validation.sh"
SUBMIT = ROOT / "validation" / "release" / "submit_real_rice.sh"
EXAMPLE_CONFIG = HARNESS / "config.example.toml"
PYBIN = str(Path(sys.executable).parent)
CLEAR = ("DRY_RUN", "CONFIG", "CONTAINER_TAG", "RELOCATE3_VALIDATION_ROOT", "SLURM_JOB_ID",
         "RELOCATE3_VALIDATION_OUTDIR", "RELOCATE3_VALIDATION_REPORT_DIR", "ALLOW_MUTABLE_TAG")


def _config(tmp_path: Path) -> Path:
    """config.example.toml with the repository root and the outputs in tmp."""
    text = EXAMPLE_CONFIG.read_text().replace('repository_root = "../.."', f'repository_root = "{ROOT}"')
    cfg = tmp_path / "cfg" / "config.toml"
    cfg.parent.mkdir(parents=True)
    cfg.write_text(text)
    return cfg


def _guard(root: Path) -> Path:
    """sbatch/srun stubs that fail: a test must never reach the real scheduler."""
    bindir = root.parent / "guard_bin"
    if not bindir.exists():
        bindir.mkdir()
        for tool in ("sbatch", "srun", "scancel"):
            (bindir / tool).write_text(f"#!/usr/bin/env bash\necho 'test guard: {tool} called' >&2\nexit 99\n")
            (bindir / tool).chmod(0o755)
    return bindir


def _env(root: Path, **extra: str) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if k not in CLEAR}
    env["RELOCATE3_VALIDATION_ROOT"] = str(root)
    env["PYTHON"] = sys.executable
    env["PATH"] = f"{_guard(root)}:{PYBIN}:{env['PATH']}"
    env.update(extra)
    return env


def _validate(tmp_path: Path, *args: str, dry: bool = True, root: Path | None = None,
              **env: str) -> subprocess.CompletedProcess:
    root = root or _fake_root(tmp_path)
    cfg = tmp_path / "cfg" / "config.toml"
    if not cfg.exists():
        _config(tmp_path)
    extra = {"DRY_RUN": "1"} if dry else {}
    extra.update(env)
    return subprocess.run(["bash", str(RUN_VALIDATION), "--config", str(cfg), *args],
                          cwd=ROOT, env=_env(root, **extra), capture_output=True, text=True)


def _commands(proc: subprocess.CompletedProcess) -> list[str]:
    return [l[2:] for l in proc.stdout.splitlines() if l.startswith("+ ")]


def _report_commands(proc: subprocess.CompletedProcess) -> list[str]:
    return [c for c in _commands(proc) if c.startswith("python3 ")]


def _value(proc: subprocess.CompletedProcess, key: str) -> str:
    lines = [l for l in proc.stdout.splitlines() if l.startswith(f"[validate-rice] {key}:")]
    assert len(lines) == 1, proc.stdout
    return lines[0].split(":", 1)[1].strip()


# ------------------------------------------------------------- _config.py


def test_config_outdir_and_report_dir_overrides(tmp_path, monkeypatch):
    sys.path.insert(0, str(HARNESS))
    try:
        from _config import load_config
    finally:
        sys.path.remove(str(HARNESS))
    cfg = _config(tmp_path)
    monkeypatch.delenv("RELOCATE3_VALIDATION_OUTDIR", raising=False)
    monkeypatch.delenv("RELOCATE3_VALIDATION_REPORT_DIR", raising=False)
    plain = load_config(cfg)["paths"]
    assert plain["relocate3_outdir"] == str(cfg.parent / "results")
    assert plain["report_dir"] == str(cfg.parent / "report")

    monkeypatch.setenv("RELOCATE3_VALIDATION_OUTDIR", str(tmp_path / "nf_out"))
    monkeypatch.setenv("RELOCATE3_VALIDATION_REPORT_DIR", str(tmp_path / "rep"))
    paths = load_config(cfg)["paths"]
    assert paths["relocate3_outdir"] == str(tmp_path / "nf_out")
    assert paths["report_dir"] == str(tmp_path / "rep")
    assert paths["log_dir"] == plain["log_dir"]
    # A relative override is relative to the caller's directory, not the config's.
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("RELOCATE3_VALIDATION_OUTDIR", "rel/out")
    assert load_config(cfg)["paths"]["relocate3_outdir"] == str(tmp_path / "rel" / "out")


# ------------------------------------------------------------- run_validation.sh


def test_default_mode_is_pipeline_and_runs_run_pipeline(tmp_path):
    out = tmp_path / "nf_out"
    proc = _validate(tmp_path, "--container-tag", "sha-abc1237", "--pipeline-outdir", str(out))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert _value(proc, "mode") == "pipeline"
    sb = [c for c in _commands(proc) if c.startswith("sbatch ")]
    assert len(sb) == 1, proc.stdout
    words = sb[0].split()
    assert "--wait" in words
    assert words[-3:] == [str(HARNESS / "run_pipeline.sh"), "sha-abc1237", str(out)]
    assert "run_relocate3.sh" not in proc.stdout
    # run_pipeline.sh's own dry run checked the inputs and printed the command.
    assert "nextflow run" in proc.stdout and "--container_tag sha-abc1237" in proc.stdout
    # The report reads the pipeline outdir.
    assert _value(proc, "relocate3_outdir") == str(out)
    # A dry run creates nothing.
    assert not out.exists() and not (tmp_path / "cfg" / "logs").exists()


def test_run_all_mode_runs_the_array(tmp_path):
    proc = _validate(tmp_path, "--mode", "run-all")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert _value(proc, "mode") == "run-all"
    sb = [c for c in _commands(proc) if c.startswith("sbatch ")]
    assert len(sb) == 1, proc.stdout
    assert "--array=1-2" in sb[0] and "run_relocate3.sh" in sb[0]
    assert "run_pipeline.sh" not in proc.stdout
    assert _value(proc, "relocate3_outdir") == str(tmp_path / "cfg" / "results")


def test_both_modes_feed_the_same_normalize_and_compare_commands(tmp_path):
    root = _fake_root(tmp_path)
    pipe = _validate(tmp_path, "--container-tag", "sha-abc1237", "--pipeline-outdir", str(tmp_path / "o"), root=root)
    runall = _validate(tmp_path, "--mode", "run-all", root=root)
    assert pipe.returncode == 0 and runall.returncode == 0, pipe.stderr + runall.stderr
    cmds = _report_commands(pipe)
    assert cmds == _report_commands(runall)
    scripts = [c.split()[1].rsplit("/", 1)[-1] for c in cmds]
    assert scripts == [
        "normalize_relocate2_nonref.py", "normalize_relocate3_nonref.py", "compare_nonref.py",
        "normalize_relocate2_char.py", "normalize_relocate3_char.py", "compare_char.py",
    ]
    assert all(f"--config {tmp_path / 'cfg' / 'config.toml'}" in c for c in cmds)


def test_pipeline_mode_needs_a_container_tag(tmp_path):
    proc = _validate(tmp_path)
    assert proc.returncode != 0
    assert "--container-tag" in proc.stderr
    assert not _commands(proc)


def test_container_tag_from_the_environment(tmp_path):
    proc = _validate(tmp_path, "--pipeline-outdir", str(tmp_path / "o"), CONTAINER_TAG="sha-abc1237")
    assert proc.returncode == 0, proc.stderr
    assert "run_pipeline.sh sha-abc1237" in "\n".join(_commands(proc))


def test_default_pipeline_outdir_is_run_pipelines_default(tmp_path):
    proc = _validate(tmp_path, "--container-tag", "sha-abc1234")
    assert proc.returncode == 0, proc.stderr
    want = ROOT / "results" / "nextflow-real-rice" / "sha-abc1234"
    assert _value(proc, "relocate3_outdir") == str(want)
    assert not want.exists()


def test_mutable_tag_is_refused_before_submitting(tmp_path):
    proc = _validate(tmp_path, "--container-tag", "dev", "--pipeline-outdir", str(tmp_path / "o"))
    assert proc.returncode != 0
    assert "ALLOW_MUTABLE_TAG=1" in proc.stderr
    assert not [c for c in _commands(proc) if c.startswith("sbatch")]


def test_skip_run_in_pipeline_mode_needs_only_the_outdir(tmp_path):
    out = tmp_path / "nf_out"
    out.mkdir()
    proc = _validate(tmp_path, "--skip-run", "--pipeline-outdir", str(out))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert not [c for c in _commands(proc) if c.startswith("sbatch")]
    assert "run_pipeline.sh" not in proc.stdout
    assert _value(proc, "relocate3_outdir") == str(out)
    assert len(_report_commands(proc)) == 6


def test_skip_run_in_pipeline_mode_with_a_missing_outdir_fails(tmp_path):
    proc = _validate(tmp_path, "--skip-run", "--pipeline-outdir", str(tmp_path / "nope"))
    assert proc.returncode != 0
    assert "nope" in proc.stderr


def test_report_dir_override(tmp_path):
    out = tmp_path / "nf_out"
    out.mkdir()
    rep = tmp_path / "rep"
    proc = _validate(tmp_path, "--skip-run", "--pipeline-outdir", str(out), "--report-dir", str(rep))
    assert proc.returncode == 0, proc.stderr
    assert _value(proc, "report_dir") == str(rep)


@pytest.mark.parametrize("mode_args", [
    ("--mode", "run-all"),
    ("--container-tag", "sha-abc1237"),
])
@pytest.mark.parametrize("var", ["RELOCATE3_VALIDATION_OUTDIR", "RELOCATE3_VALIDATION_REPORT_DIR"])
def test_inherited_output_overrides_are_refused(tmp_path, mode_args, var):
    """An output override inherited from the caller's environment would silently
    redirect the run (exported through sbatch) and the report; the wrapper owns
    these variables (--pipeline-outdir / --report-dir) and refuses a caller's."""
    proc = _validate(tmp_path, *mode_args, **{var: str(tmp_path / "elsewhere")})
    assert proc.returncode != 0
    assert var in proc.stderr
    assert not _commands(proc)


def test_pipeline_mode_default_report_dir_is_per_tag(tmp_path):
    proc = _validate(tmp_path, "--container-tag", "sha-abc1237", "--pipeline-outdir", str(tmp_path / "o"))
    assert proc.returncode == 0, proc.stderr
    assert _value(proc, "report_dir") == str(tmp_path / "cfg" / "report" / "sha-abc1237")


def test_pipeline_skip_run_without_tag_reports_per_outdir(tmp_path):
    out = tmp_path / "nf_out" / "sha-abc1237"
    out.mkdir(parents=True)
    proc = _validate(tmp_path, "--skip-run", "--pipeline-outdir", str(out))
    assert proc.returncode == 0, proc.stderr
    assert _value(proc, "report_dir") == str(tmp_path / "cfg" / "report" / "sha-abc1237")


def test_run_all_mode_keeps_the_config_report_dir(tmp_path):
    proc = _validate(tmp_path, "--mode", "run-all")
    assert proc.returncode == 0, proc.stderr
    assert _value(proc, "report_dir") == str(tmp_path / "cfg" / "report")


@pytest.mark.parametrize("args", [
    ("--local",), ("--force",), ("S_1",), ("--mode", "other"),
])
def test_pipeline_mode_refuses_run_all_only_options(tmp_path, args):
    proc = _validate(tmp_path, "--container-tag", "sha-abc1237", *args)
    assert proc.returncode != 0
    assert not _commands(proc)


def test_pipeline_outdir_is_refused_in_run_all_mode(tmp_path):
    proc = _validate(tmp_path, "--mode", "run-all", "--pipeline-outdir", str(tmp_path / "o"))
    assert proc.returncode != 0
    assert "--pipeline-outdir" in proc.stderr


def _relocate2_gff(root: Path, sample: str, starts: list[int]) -> None:
    d = root / "relocate2_results" / sample / "repeat" / "results"
    d.mkdir(parents=True)
    (d / "ALL.all_nonref_insert.gff").write_text("".join(_gff(sample, s) for s in starts))


def _gff(sample: str, start: int) -> str:
    return (f"Chr1\tRelocaTE\t{sample}\t{start}\t{start + 2}\t.\t+\t.\t"
            "ID=x;Name=mPing;TSD=TAA;Note=Non-reference, not found in reference;"
            "Right_junction_reads=3;Left_junction_reads=4;Right_support_reads=0;Left_support_reads=0;\n")


def test_pipeline_mode_submits_waits_and_reports_on_the_pipeline_outdir(tmp_path):
    """Not a dry run: a stub sbatch records its call; the report then reads the
    (fake) pipeline outdir, not the config's relocate3_outdir."""
    root = _fake_root(tmp_path)
    out = tmp_path / "nf_out"
    for s, starts in (("S_1", [1000, 5000]), ("S_2", [2000])):
        _relocate2_gff(root, s, starts)
        res = out / s / "results"
        res.mkdir(parents=True)
        (res / "ALL.mping.all_nonref_insert.gff").write_text(_gff(s, starts[0]))
    stub = tmp_path / "stub"
    stub.mkdir()
    (stub / "sbatch").write_text(f'#!/usr/bin/env bash\nprintf "%s\\n" "$@" > {tmp_path}/sbatch.argv\n'
                                 f'echo "CONFIG=$CONFIG" >> {tmp_path}/sbatch.argv\n')
    (stub / "sbatch").chmod(0o755)
    rep = tmp_path / "rep"
    proc = _validate(tmp_path, "--container-tag", "sha-abc1237", "--pipeline-outdir", str(out),
                     "--report-dir", str(rep), "--skip-char", dry=False, root=root,
                     PATH=f"{stub}:{PYBIN}:{os.environ['PATH']}")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    argv = (tmp_path / "sbatch.argv").read_text().splitlines()
    assert "--wait" in argv
    assert argv[-4:-1] == [str(HARNESS / "run_pipeline.sh"), "sha-abc1237", str(out)]
    assert argv[-1] == f"CONFIG={tmp_path / 'cfg' / 'config.toml'}"
    r3 = (rep / "nonref" / "relocate3_calls.tsv").read_text().splitlines()
    assert len(r3) == 3 and all(str(out) in l for l in r3[1:])
    summary = (rep / "nonref" / "summary.txt").read_text()
    assert summary and "Stage A summary" in proc.stdout
    assert not (tmp_path / "cfg" / "report").exists()


# ------------------------------------------------------------- submit_real_rice.sh


def _clone(tmp_path: Path) -> Path:
    """Base checkout + fresh clone layout that submit_real_rice.sh expects."""
    base = tmp_path / "base"
    clone = base / "results" / "release-validation" / "v9.9.9rc1" / "clone"
    shutil.copytree(HARNESS, clone / "validation" / "real_rice",
                    ignore=shutil.ignore_patterns("__pycache__", "config.toml", "results", "report", "logs"))
    (clone / "tests" / "data" / "TE_lib").mkdir(parents=True)
    (clone / "tests" / "data" / "TE_lib" / "mping.fa").write_text(">mPing\nACGT\n")
    (clone / "main.nf").write_text("")
    subprocess.run(["git", "init", "-q", str(clone)], check=True)
    subprocess.run(["git", "-C", str(clone), "-c", "user.email=t@t", "-c", "user.name=t",
                    "commit", "-q", "--allow-empty", "-m", "x"], check=True)
    for b in (clone / ".pixi" / "envs" / "blat" / "bin" / "relocaTE3",
              clone / ".pixi" / "envs" / "blat" / "bin" / "python",
              base / ".pixi" / "envs" / "default" / "bin" / "python3"):
        b.parent.mkdir(parents=True, exist_ok=True)
        b.write_text("#!/usr/bin/env bash\n")
        b.chmod(0o755)
    return base


def _submit(tmp_path: Path, *args: str, dry: bool = True) -> tuple[subprocess.CompletedProcess, Path]:
    root = _fake_root(tmp_path)
    base = _clone(tmp_path)
    extra = {"DRY_RUN": "1"} if dry else {}
    proc = subprocess.run(["bash", str(SUBMIT), "v9.9.9rc1", *args], cwd=base,
                          env=_env(root, **extra), capture_output=True, text=True)
    return proc, base


def test_release_check_needs_a_container_tag_in_pipeline_mode(tmp_path):
    proc, base = _submit(tmp_path)
    assert proc.returncode != 0
    assert "--container-tag" in proc.stderr
    assert not (base / "results" / "release-validation" / "v9.9.9rc1" / "real_rice.submission").exists()


def test_release_check_pipeline_mode_chains_head_job_and_report(tmp_path):
    proc, base = _submit(tmp_path, "--container-tag", "0.3.0rc1")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    rr = base / "results" / "release-validation" / "v9.9.9rc1" / "clone" / "validation" / "real_rice"
    cmds = _commands(proc)
    assert len(cmds) == 2, proc.stdout
    head, report = cmds
    assert head.startswith("sbatch ") and "run_pipeline.sh 0.3.0rc1 " + str(rr / "results") in head
    assert "--array" not in head
    assert report.startswith("sbatch ") and "--dependency=afterok:" in report
    assert "run_validation.sh --mode pipeline --skip-run" in report
    assert f"--pipeline-outdir {rr / 'results'}" in report
    # The report goes to the tag's own report dir (report/<tag>).
    assert "--container-tag 0.3.0rc1" in report
    # A dry run submits and writes nothing.
    assert not (rr / "config.toml").exists()
    assert not (base / "results" / "release-validation" / "v9.9.9rc1" / "real_rice.submission").exists()


@pytest.mark.parametrize("tag", ["dev", "pr-12"])
def test_release_check_refuses_mutable_tags(tmp_path, tag):
    proc, _ = _submit(tmp_path, "--container-tag", tag)
    assert proc.returncode != 0
    assert tag in proc.stderr and not _commands(proc)


def test_release_check_run_all_mode_keeps_the_array(tmp_path):
    proc, _ = _submit(tmp_path, "--mode", "run-all", "--partition", "intel")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    head, report = _commands(proc)
    assert "--array=1-2" in head and "-p intel" in head and "run_relocate3.sh" in head
    assert "run_validation.sh --mode run-all --skip-run" in report
    assert "run_pipeline.sh" not in proc.stdout
