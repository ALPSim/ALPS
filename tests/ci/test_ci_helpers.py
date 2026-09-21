"""Check matrix selection and summaries using the actual Actions interfaces."""

import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize(
    ("event", "ref", "tier", "count"),
    [
        ("pull_request", "refs/pull/1/merge", "auto", 4),
        ("merge_group", "refs/heads/gh-readonly-queue/master/test", "auto", 4),
        ("push", "refs/heads/master", "auto", 4),
        ("push", "refs/tags/v3.0.0", "auto", 14),
        ("schedule", "refs/heads/master", "auto", 14),
        ("workflow_dispatch", "refs/heads/feature", "full", 14),
        ("workflow_dispatch", "refs/heads/feature", "quick", 4),
    ],
)
def test_matrix_event_selection(tmp_path, event, ref, tier, count):
    output = tmp_path / "output"
    summary = tmp_path / "summary"
    subprocess.run(
        [sys.executable, str(ROOT / ".github/scripts/ci_matrix.py"), "--tier", tier],
        env={**os.environ, "GITHUB_EVENT_NAME": event, "GITHUB_REF": ref,
             "GITHUB_OUTPUT": str(output), "GITHUB_STEP_SUMMARY": str(summary)},
        check=True, capture_output=True, text=True,
    )
    key, value = output.read_text().strip().split("=", 1)
    assert key == "matrix"
    builds = json.loads(value)["include"]
    assert len(builds) == count
    assert len({build["id"] for build in builds}) == count
    assert all(build["boost_sha256"] for build in builds)
    assert any(build["extras"] for build in builds)
    assert any(build["os"] == "macos-15" for build in builds)
    assert any(build["boost"] == "1.76.0" for build in builds)
    if count == 14:
        assert any(build["extensive"] for build in builds)
        sanitized = [build for build in builds if build["sanitizer"]]
        assert sanitized and all(build["mpi"] == "OFF" for build in sanitized)
    assert f"{count} builds" in summary.read_text()


def load_helper(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / f".github/scripts/{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("problem", ["duplicate", "checksum", "empty"])
def test_reject_broken_matrix(problem):
    manifest = json.loads((ROOT / ".github/ci-matrix.json").read_text())
    if problem == "duplicate":
        manifest["builds"].append(manifest["builds"][0])
    elif problem == "checksum":
        manifest["boost"]["1.76.0"] = "not a checksum"
    else:
        manifest["builds"] = []
    with pytest.raises(ValueError):
        load_helper("ci_matrix").select_matrix(manifest, "full")


def test_mixed_junit_results(tmp_path):
    report = tmp_path / "results.xml"
    report.write_text(
        '<testsuites><testsuite><testcase name="pass"/>'
        '<testcase name="failure"><failure>bad result</failure></testcase>'
        '<testcase name="error"><error>setup failed</error></testcase>'
        '<testcase name="skip"><skipped/></testcase></testsuite></testsuites>'
    )
    result = load_helper("junit_summary").summarize([report])
    assert "| 1 | 2 | 1 |" in result


def test_missing_junit_is_explicit():
    result = load_helper("junit_summary").summarize([])
    assert "No test reports" in result
