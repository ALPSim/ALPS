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
        ("pull_request", "refs/pull/1/merge", "auto", 1),
        ("merge_group", "refs/heads/gh-readonly-queue/master/test", "auto", 1),
        ("push", "refs/heads/master", "auto", 4),
        ("push", "refs/tags/v3.0.0", "auto", 14),
        ("schedule", "refs/heads/master", "auto", 14),
        ("workflow_dispatch", "refs/heads/feature", "full", 14),
        ("workflow_dispatch", "refs/heads/feature", "quick", 1),
    ],
)
def test_matrix_event_selection(tmp_path, event, ref, tier, count):
    output = tmp_path / "output"
    summary = tmp_path / "summary"
    payload = tmp_path / "event.json"
    payload.write_text("{}")
    subprocess.run(
        [sys.executable, str(ROOT / ".github/scripts/ci_matrix.py"), "--tier", tier],
        env={**os.environ, "GITHUB_EVENT_NAME": event, "GITHUB_REF": ref,
             "GITHUB_EVENT_PATH": str(payload), "GITHUB_OUTPUT": str(output),
             "GITHUB_STEP_SUMMARY": str(summary)},
        check=True, capture_output=True, text=True,
    )
    values = dict(line.split("=", 1) for line in output.read_text().splitlines())
    assert values["developer"] == str(count == 14 or event in {"pull_request", "merge_group"}).lower()
    value = values["matrix"]
    builds = json.loads(value)["include"]
    assert len(builds) == count
    assert len({build["id"] for build in builds}) == count
    assert all(build["boost_sha256"] for build in builds)
    if count > 1:
        assert any(build["extras"] for build in builds)
        assert any(build["os"] == "macos-15" for build in builds)
        assert any(build["boost"] == "1.76.0" for build in builds)
    else:
        assert builds[0]["python"] and builds[0]["mpi"] == "OFF"
        assert not builds[0]["extras"]
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


@pytest.mark.parametrize(("tier", "only", "developer", "expected"), [
    ("full", "clang-sanitizers", "false", "false"),
    ("full", "clang-sanitizers", "true", "true"),
    ("quick", "", "true", "true"),
])
def test_targeted_build_and_developer_outputs(tmp_path, tier, only, developer, expected):
    output = tmp_path / "output"
    payload = tmp_path / "event.json"
    payload.write_text("{}")
    subprocess.run([
        sys.executable, str(ROOT / ".github/scripts/ci_matrix.py"),
        "--tier", tier, "--only", only,
    ], env={**os.environ, "GITHUB_EVENT_NAME": "workflow_dispatch",
            "GITHUB_REF": "refs/heads/feature", "GITHUB_EVENT_PATH": str(payload),
            "GITHUB_OUTPUT": str(output), "DEVELOPER_SETUP": developer},
        check=True, capture_output=True, text=True)
    values = dict(line.split("=", 1) for line in output.read_text().splitlines())
    builds = json.loads(values["matrix"])["include"]
    assert len(builds) == 1
    if only:
        assert builds[0]["id"] == only
    assert values["developer"] == expected


@pytest.mark.parametrize(("tier", "only"), [
    ("full", "unknown"), ("quick", "clang-sanitizers"),
])
def test_reject_unavailable_build_selection(tier, only):
    manifest = json.loads((ROOT / ".github/ci-matrix.json").read_text())
    with pytest.raises(ValueError, match="Empty"):
        load_helper("ci_matrix").select_matrix(manifest, tier, only)


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
