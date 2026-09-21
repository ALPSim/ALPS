"""Reject a mismatched SDK before compiling or labelling a Python wheel."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

SOURCE = Path(__file__).resolve().parents[2]
pytestmark = pytest.mark.skipif(not os.environ.get("ALPS_DIR"), reason="requires an installed SDK")


@pytest.mark.parametrize("advertised_version", ["99.0.0", None])
def test_python_rejects_mismatched_sdk(tmp_path, advertised_version):
    pytest.importorskip("nanobind")
    core = (SOURCE / "ALPS_VERSION.txt").read_text().strip()
    sdk = tmp_path / "sdk"
    sdk.mkdir()
    # The second case also catches configs whose version file and payload disagree.
    (sdk / "ALPSConfigVersion.cmake").write_text(
        f'set(PACKAGE_VERSION "{advertised_version or core}")\n'
        'if(PACKAGE_FIND_VERSION VERSION_EQUAL PACKAGE_VERSION)\n'
        '  set(PACKAGE_VERSION_COMPATIBLE TRUE)\n'
        '  set(PACKAGE_VERSION_EXACT TRUE)\n'
        'endif()\n', encoding="utf-8")
    (sdk / "ALPSConfig.cmake").write_text('set(ALPS_VERSION "99.0.0")\n', encoding="utf-8")
    environment = dict(os.environ)
    # Do not let find_package fall back to the real SDK after rejecting the mock.
    environment.pop("ALPS_DIR", None)
    result = subprocess.run([
        "cmake", "-S", str(SOURCE / "python/pyalps"), "-B", str(tmp_path / "build"),
        "-G", "Ninja", "-DALPS_DIR=" + str(sdk), "-DPython_EXECUTABLE=" + sys.executable,
        "-DPYALPS_BUILD_SOLVERS=OFF", "-DPYALPS_BUNDLE_APPLICATIONS=OFF",
        *json.loads(os.environ.get("ALPS_TEST_CMAKE_ARGS", "[]")),
    ], capture_output=True, text=True, env=environment)
    assert result.returncode != 0
    assert "99.0.0" in result.stdout + result.stderr
    assert core in result.stdout + result.stderr
