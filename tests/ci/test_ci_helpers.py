"""Validate test reporting and installed-wheel execution."""

import importlib.util
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]


def load_helper(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / f".github/scripts/{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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


def test_junit_summary_includes_elapsed_time(tmp_path):
    report = tmp_path / "results.xml"
    report.write_text('<testsuite><testcase name="a" time="1.25"/>'
                      '<testcase name="b" time="0.5"/></testsuite>')
    result = load_helper("junit_summary").summarize([report])
    assert "Seconds" in result
    assert "| 2 | 0 | 0 | 1.75 |" in result


@pytest.mark.parametrize("source,code", [
    ("def test_result():\n    assert True\n", 0),
    ("def test_result():\n    assert False\n", 1),
    ("# no tests\n", 5),
])
@pytest.mark.parametrize("version", [(3, 11), (3, 12), (3, 14)])
def test_wheel_runner_preserves_test_failures_and_empty_collection(tmp_path, monkeypatch, source, code, version):
    import types
    root = tmp_path / "checkout"
    for suite in ("pyalps", "cmake"):
        (root / "tests" / suite).mkdir(parents=True)
    (root / "tests/pyalps/test_binding_surface.py").write_text(source)
    (root / "tests/pyalps/test_mapping_lifetimes.py").touch()
    (root / "tests/pyalps/test_wheel_payload.py").touch()
    helper = load_helper("run_wheel_tests")
    monkeypatch.setattr(helper, "sys", types.SimpleNamespace(version_info=version, executable=sys.executable))
    monkeypatch.setattr(helper, "ROOT", root)
    monkeypatch.setitem(sys.modules, "pyalps", types.SimpleNamespace(
        __file__=str(tmp_path / "installed/pyalps/__init__.py")))
    assert helper.main() == code
    reports = list((root / "_build/wheel-reports").glob("*.xml"))
    assert len(reports) == 1
    assert "test_binding_surface" in reports[0].read_text() or code == 5


def test_wheel_runner_rejects_source_tree_import(monkeypatch):
    import types
    helper = load_helper("run_wheel_tests")
    monkeypatch.setitem(sys.modules, "pyalps", types.SimpleNamespace(
        __file__=str(ROOT / "python/pyalps/src/pyalps/__init__.py")))
    with pytest.raises(RuntimeError, match="source-tree"):
        helper.main()
