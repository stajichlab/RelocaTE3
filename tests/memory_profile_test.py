"""Lightweight tests for the profiling harness, not the bioinformatics pipeline."""

import importlib.util
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest


SPEC = importlib.util.spec_from_file_location(
    "memory_profile",
    Path(__file__).resolve().parents[1] / "scripts/profile_benchmark_memory.py",
)
profile = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(profile)


def test_call_comparison_ignores_order_but_preserves_duplicate_counts(tmp_path):
    first, second = tmp_path / "first.tsv", tmp_path / "second.tsv"
    first.write_text("chrom\tposition\nChr1\t2\nChr1\t1\n")
    second.write_text("position\tchrom\n1\tChr1\n2\tChr1\n")
    assert profile.table_rows(first) == profile.table_rows(second)
    second.write_text("position\tchrom\n1\tChr1\n2\tChr1\n2\tChr1\n")
    assert profile.table_rows(first) != profile.table_rows(second)


def test_process_session_filter():
    own = Path("/proc") / str(os.getpid())
    assert profile.read_process(own, -1) is None
    row = profile.read_process(own, os.getsid(0))
    assert row["pid"] == os.getpid()
    assert row["hwm_kib"] >= row["rss_kib"] > 0


def test_monitor_captures_parent_and_child(tmp_path):
    output = tmp_path / "profile"
    child = "import time; payload=bytearray(8*1024*1024); time.sleep(0.5)"
    parent = "import subprocess,sys; print('adapter marker', flush=True); subprocess.run([sys.executable, '-c', sys.argv[1]], check=True)"
    code = profile.monitor(
        [sys.executable, "-c", parent, child],
        tmp_path,
        os.environ.copy(),
        output,
        interval=0.02,
    )
    assert code == 0
    report = json.loads((output / "summary.json").read_text())
    assert len({row["pid"] for row in report["processes"]}) >= 2
    assert report["max_sampled_rss_sum_kib"] > 8 * 1024
    assert "adapter marker" in (output / "adapter.log").read_text()


def test_monitor_preserves_failure_and_refuses_overwrite(tmp_path):
    output = tmp_path / "profile"
    command = [sys.executable, "-c", "raise SystemExit(7)"]
    assert (
        profile.monitor(command, tmp_path, os.environ.copy(), output, interval=0.01)
        == 7
    )
    assert json.loads((output / "summary.json").read_text())["exit_code"] == 7
    with pytest.raises(FileExistsError):
        profile.monitor(command, tmp_path, os.environ.copy(), output)


def test_replay_refuses_login_node(monkeypatch, tmp_path):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="SLURM"):
        profile.run(SimpleNamespace(output=tmp_path))


def test_json_is_exclusive(tmp_path):
    path = tmp_path / "manifest.json"
    profile.write_json(path, {"first": True})
    with pytest.raises(FileExistsError):
        profile.write_json(path, {"first": False})
    assert json.loads(path.read_text()) == {"first": True}


@pytest.fixture
def trim_validator(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    import validate_trim_memory

    return validate_trim_memory


def test_trim_artifact_hashes_detect_empty_missing_extra_and_changed_files(
    tmp_path, trim_validator
):
    for name in ("flanking", "te_containing", "te_portions"):
        (tmp_path / name).mkdir()
    empty = tmp_path / "te_portions/empty.fa"
    empty.touch()
    first = trim_validator.artifact_hashes(tmp_path)
    assert set(first) == {"te_portions/empty.fa"}
    (tmp_path / "flanking/flank.fq").write_text("@read\nAC\n+\nII\n")
    second = trim_validator.artifact_hashes(tmp_path)
    assert first != second
    empty.write_text(">read\nA\n")
    assert trim_validator.artifact_hashes(tmp_path) != second


def test_trim_validation_refuses_login_node(tmp_path, monkeypatch, trim_validator):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="SLURM"):
        trim_validator.validate(tmp_path, tmp_path)


def test_snapshot_verification_detects_changes(tmp_path, trim_validator):
    path = tmp_path / "source.py"
    path.write_text("original")
    manifest = {"snapshot_sha256": {"source.py": trim_validator.digest(path)}}
    trim_validator.verify_snapshot(tmp_path, manifest)
    path.write_text("changed")
    with pytest.raises(ValueError, match="Snapshot changed"):
        trim_validator.verify_snapshot(tmp_path, manifest)
