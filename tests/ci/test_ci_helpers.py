"""Preserve upstream source coverage through the SDK/Python migration."""

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


def test_source_matrix_preserves_upstream_compatibility():
    manifest = json.loads((ROOT / ".github/ci-matrix.json").read_text())
    builds = load_helper("ci_matrix").select_matrix(manifest)["include"]
    # The upstream source matrix, minus only the intentionally unsupported
    # Python 3.9/3.10 rows. Compare combinations, not just totals or versions.
    expected = {
        (os_name, compiler, "3.14", "1.91.0", 17)
        for os_name, compiler in [
            ("ubuntu-22.04", "gcc-11"), ("ubuntu-22.04", "gcc-12"),
            ("ubuntu-24.04", "gcc-13"), ("ubuntu-24.04", "gcc-14"),
            ("ubuntu-24.04", "gcc-15"),
            ("ubuntu-22.04", "clang-14"), ("ubuntu-22.04", "clang-15"),
            *[("ubuntu-24.04", f"clang-{version}") for version in range(16, 23)],
            *[(os_name, "/usr/bin/clang") for os_name in
              ("macos-14", "macos-15", "macos-15-intel", "macos-26")],
            ("macos-15", "gcc-13"), ("macos-15", "gcc-14"),
        ]
    }
    expected |= {("ubuntu-24.04", "gcc-14", "3.14", f"1.{v}.0", 17)
                 for v in (76, 81, 86, 87, 88, 89, 90)}
    expected |= {("ubuntu-24.04", "gcc-14", v, "1.91.0", 17)
                 for v in ("3.11", "3.12", "3.13")}
    expected |= {("ubuntu-24.04", "gcc-14", "3.14", "1.91.0", v)
                 for v in (20, 23)}
    actual = {(b["os"], b["cc"], b["python"], b["boost"], b["standard"])
              for b in builds}
    assert len(builds) == len(actual) == len(expected) == 32
    assert actual == expected
    assert all(b["mpi"] == "ON" for b in builds)
    assert all(b["repository"] == "llvm" for b in builds
               if b["cc"] in {f"clang-{v}" for v in range(19, 23)})
    assert all(b["packages"] == b["cc"].replace("gcc-", "gcc@") for b in builds
               if b["os"].startswith("macos-") and b["cc"].startswith("gcc-"))


@pytest.mark.parametrize("event,paths", [
    ("pull_request", ["README.md"]),
    ("pull_request", ["python/pyalps/src/pyalps/tools.py"]),
    ("pull_request", ["src/alps/alea/src/alea/observable.C"]),
    ("push", []), ("workflow_dispatch", []),
])
def test_matrix_is_not_reduced_by_event_or_paths(tmp_path, event, paths):
    output = tmp_path / "output"
    summary = tmp_path / "summary"
    payload = tmp_path / "event.json"
    payload.write_text(json.dumps({"commits": [{"modified": paths}]}))
    subprocess.run(
        [sys.executable, str(ROOT / ".github/scripts/ci_matrix.py")],
        env={**os.environ, "GITHUB_EVENT_NAME": event,
             "GITHUB_EVENT_PATH": str(payload), "GITHUB_OUTPUT": str(output),
             "GITHUB_STEP_SUMMARY": str(summary)},
        check=True, capture_output=True, text=True,
    )
    values = dict(line.split("=", 1) for line in output.read_text().splitlines())
    builds = json.loads(values["matrix"])["include"]
    assert len(builds) == 32
    assert len({build["id"] for build in builds}) == 32
    assert "all 32 builds" in summary.read_text()


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
