"""Guard unconditional CI coverage and the supported wheel/runtime matrix."""

import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / ".github/scripts/packaging_matrix.py"
SPEC = importlib.util.spec_from_file_location("packaging_matrix", SCRIPT)
packaging = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(packaging)


def test_packaging_builds_all_families_and_tests_every_supported_python():
    matrices = packaging.packaging_matrices()
    assert {p["family"] for p in matrices["wheel_matrix"]["plat"]} == {
        "manylinux", "musllinux", "macos"}
    built = {"cibw-wheels-" + p["family"] for p in matrices["wheel_matrix"]["plat"]}
    smoke = matrices["smoke_matrix"]
    assert {p["artifact"] for p in smoke["plat"]} <= built
    assert {p["os"] for p in smoke["plat"]} == {"ubuntu-24.04", "macos-15", "macos-26"}
    assert smoke["python"] == ["3.11", "3.12", "3.13", "3.14"]


@pytest.mark.parametrize("event", ["pull_request", "push", "workflow_dispatch"])
def test_packaging_outputs_are_complete_for_every_event(tmp_path, event):
    output = tmp_path / "output"
    subprocess.run([sys.executable, str(SCRIPT)], check=True, capture_output=True,
                   env={**os.environ, "GITHUB_EVENT_NAME": event, "GITHUB_OUTPUT": str(output)})
    values = {key: json.loads(value) for key, value in
              (line.split("=", 1) for line in output.read_text().splitlines())}
    assert values == packaging.packaging_matrices()


def test_workflows_do_not_filter_coverage():
    workflows = SCRIPT.parents[1] / "workflows"
    for name in ("build.yml", "build_wheels.yml"):
        text = (workflows / name).read_text()
        triggers = text.split("on:\n", 1)[1].split("permissions:\n", 1)[0]
        assert "  push:" in triggers and "  pull_request:" in triggers
        assert "tags:" in triggers and "v*" in triggers and "master" in triggers
        for reduction in ("paths:", "paths-ignore:", "schedule:", "merge_group:", "tier:"):
            assert reduction not in triggers
        assert "ci_policy" not in text
        assert "outputs.packaging" not in text
        assert "outputs.source" not in text
        assert "outputs.developer" not in text
    source = (workflows / "build.yml").read_text()
    assert "python-version: ${{ matrix.python }}" in source
    assert "python -m pytest tests/pyalps -q" in source
    assert "python -m pytest tests/pyalps tests/cmake -q" in source
    assert "if: matrix.python" not in source
    assert "llvm_version=\"${CC#clang-}\"" in source


def test_native_windows_is_manual_and_outside_release_gates():
    workflows = SCRIPT.parents[1] / 'workflows'
    native = (workflows / 'windows.yml').read_text()
    triggers = native.split('on:\n', 1)[1].split('permissions:\n', 1)[0]
    assert '  workflow_dispatch:' in triggers
    assert 'workflow_call:' not in triggers
    assert 'pull_request:' not in triggers
    assert 'push:' not in triggers
    assert 'schedule:' not in triggers
    for name in ('build.yml', 'build_wheels.yml'):
        required = (workflows / name).read_text().lower()
        assert 'windows' not in required, f'{name} must not depend on native Windows'
