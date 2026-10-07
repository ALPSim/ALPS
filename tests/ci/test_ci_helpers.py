"""Validate matrix expansion and CI result summaries."""

import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]


def load_helper(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / f".github/scripts/{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_source_matrix_expands_every_manifest_entry():
    manifest = json.loads((ROOT / ".github/ci-matrix.json").read_text())
    builds = load_helper("ci_matrix").select_matrix(manifest)["include"]
    assert len(builds) == len(manifest["builds"])
    for entry, build in zip(manifest["builds"], builds):
        assert entry.items() <= build.items()
        assert build["boost_sha256"] == manifest["boost"][build["boost"]]
        assert build["dependency"] == f"boost-{build['boost']}"


def test_matrix_cli_writes_github_outputs(tmp_path):
    output = tmp_path / "output"
    summary = tmp_path / "summary"
    subprocess.run(
        [sys.executable, str(ROOT / ".github/scripts/ci_matrix.py")],
        env={**os.environ, "GITHUB_OUTPUT": str(output),
             "GITHUB_STEP_SUMMARY": str(summary)},
        check=True, capture_output=True, text=True,
    )
    values = dict(line.split("=", 1) for line in output.read_text().splitlines())
    builds = json.loads(values["matrix"])["include"]
    manifest = json.loads((ROOT / ".github/ci-matrix.json").read_text())
    assert len(builds) == len(manifest["builds"])
    assert len({build["id"] for build in builds}) == len(builds)
    assert str(len(builds)) in summary.read_text()


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
        load_helper("ci_matrix").select_matrix(manifest)


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
