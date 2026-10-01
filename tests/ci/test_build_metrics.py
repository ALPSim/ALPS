"""Performance reporting must preserve the command's actual result."""

import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPT = Path(__file__).resolve().parents[2] / ".github/scripts/build_metrics.py"
SPEC = importlib.util.spec_from_file_location("build_metrics", SCRIPT)
metrics = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(metrics)


@pytest.mark.parametrize("code", [0, 7])
def test_command_exit_status_and_output_are_preserved(tmp_path, code):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "run", "SDK build", sys.executable, "-c",
         f"print('compiler output'); raise SystemExit({code})"],
        env={**os.environ, "ALPS_CI_METRICS": str(tmp_path)},
        text=True, capture_output=True,
    )
    assert result.returncode == code
    assert "compiler output" in result.stdout
    row = json.loads((tmp_path / "phases.jsonl").read_text())
    assert row["exit_code"] == code
    assert row["seconds"] >= 0


def test_timeout_is_a_reported_failure(tmp_path):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "run", "--timeout", "0.1", "Native test",
         sys.executable, "-c", "import time; time.sleep(30)"],
        env={**os.environ, "ALPS_CI_METRICS": str(tmp_path)},
        text=True, capture_output=True,
    )
    assert result.returncode == 124
    assert json.loads((tmp_path / "phases.jsonl").read_text())["exit_code"] == 124


def test_summary_distinguishes_restored_cache_from_actual_hits(tmp_path):
    (tmp_path / "phases.jsonl").write_text(json.dumps({
        "phase": "SDK build", "seconds": 12.5, "exit_code": 0, "hits": 4, "misses": 8,
    }) + "\n")
    result = metrics.summarize(tmp_path, "previous-run")
    assert "restored: **yes**" in result
    assert "| SDK build | 12.5 | 0 | 4 | 8 |" in result
    assert "no (cold cache)" in metrics.summarize(tmp_path)


def test_actual_ccache_counters_are_reported_per_phase(tmp_path, monkeypatch):
    snapshots = iter([
        "direct_cache_hit\t20\npreprocessed_cache_hit\t3\ncache_miss\t10\n",
        "direct_cache_hit\t24\npreprocessed_cache_hit\t5\ncache_miss\t11\n",
    ])

    def execute(command, **kwargs):
        if command == ["ccache", "--print-stats"]:
            return subprocess.CompletedProcess(command, 0, next(snapshots))
        return subprocess.CompletedProcess(command, 0, "")

    monkeypatch.setattr(metrics.subprocess, "run", execute)
    monkeypatch.setattr(metrics, "DIRECTORY", tmp_path)
    assert metrics.run("SDK build", ["cmake", "--build", "build"]) == 0
    row = json.loads((tmp_path / "phases.jsonl").read_text())
    assert row["hits"] == 6
    assert row["misses"] == 1


def test_all_hit_run_does_not_save_a_duplicate_cache(tmp_path):
    assert not metrics.cache_changed(tmp_path)
    phases = tmp_path / "phases.jsonl"
    phases.write_text(json.dumps({"hits": 800, "misses": 0}) + "\n")
    assert not metrics.cache_changed(tmp_path)
    with phases.open("a") as stream:
        stream.write(json.dumps({"hits": 0, "misses": 1}) + "\n")
    assert metrics.cache_changed(tmp_path)
