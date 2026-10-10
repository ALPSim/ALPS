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


@pytest.mark.parametrize("source,code", [
    ("def test_result():\n    assert True\n", 0),
    ("def test_result():\n    assert False\n", 1),
    ("# no tests\n", 5),
])
@pytest.mark.parametrize("version", [(3, 11), (3, 12), (3, 14)])
def test_wheel_runner_preserves_test_failures_and_empty_collection(tmp_path, monkeypatch, source, code, version):
    import types
    root = tmp_path / "checkout"
    for suite in ("python/pyalps/tests", "tests/cmake"):
        (root / suite).mkdir(parents=True)
    (root / "python/pyalps/tests/test_binding_surface.py").write_text(source)
    (root / "python/pyalps/tests/test_mapping_lifetimes.py").touch()
    (root / "python/pyalps/tests/test_wheel_payload.py").touch()
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
