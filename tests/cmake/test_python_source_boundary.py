"""Build the application bindings with only Python-owned sources available."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

SOURCE = Path(__file__).resolve().parents[2]
PACKAGE = SOURCE / "python/pyalps"
pytestmark = pytest.mark.skipif(not os.environ.get("ALPS_DIR"), reason="requires an installed SDK")


def test_standalone_python_sources(tmp_path):
    nanobind = pytest.importorskip("nanobind")
    if not (Path(os.environ["ALPS_DIR"]) / "ALPSApplicationTargets.cmake").is_file():
        pytest.skip("requires an SDK with solvers and applications")
    source = tmp_path / "package"
    # Deliberately provide no repository ancestors or vendored C++ sources.
    shutil.copytree(PACKAGE, source, ignore=shutil.ignore_patterns(
        "dist", "build", "_skbuild", "__pycache__", ".pytest_cache"))
    shutil.copy2(SOURCE / "ALPS_VERSION.txt", source)
    build = tmp_path / "build"
    subprocess.run([
        "cmake", "-S", str(source), "-B", str(build), "-G", "Ninja",
        "-DCMAKE_BUILD_TYPE=Release", "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DALPS_DIR=" + os.environ["ALPS_DIR"], "-DPython_EXECUTABLE=" + sys.executable,
        *json.loads(os.environ.get("ALPS_TEST_CMAKE_ARGS", "[]")),
    ], check=True)
    commands = json.loads((build / "compile_commands.json").read_text())
    nanobind_root = Path(nanobind.__file__).resolve().parent
    for command in commands:
        assert "Py_LIMITED_API=0x030C0000" in command["command"]
        compiled_source = Path(command["file"]).resolve()
        assert (compiled_source.is_relative_to(source)
                or compiled_source.is_relative_to(nanobind_root)), compiled_source
    # Link the three former source-boundary violations against the installed SDK.
    subprocess.run([
        "cmake", "--build", str(build), "--target", "maxent_c", "cthyb", "ctint",
        "--config", "Release", "--parallel", "2",
    ], check=True)
    # Core-only builds retain the same source boundary and need no solver targets.
    subprocess.run([
        "cmake", "-S", str(source), "-B", str(build),
        "-DPYALPS_BUILD_SOLVERS=OFF", "-DPYALPS_BUNDLE_APPLICATIONS=OFF",
    ], check=True)
