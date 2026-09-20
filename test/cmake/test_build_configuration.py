"""Configure-only contracts: exercise build choices without duplicate objects."""
import json
import os
from pathlib import Path
import subprocess

import pytest

SOURCE = Path(__file__).resolve().parents[2]
pytestmark = pytest.mark.skipif(
    not os.environ.get("ALPS_DIR"), reason="requires configured SDK dependencies")


def configure(build, *options):
    return subprocess.run([
        "cmake", "-S", str(SOURCE), "-B", str(build),
        *json.loads(os.environ.get("ALPS_TEST_CMAKE_ARGS", "[]")),
        "-DALPS_ENABLE_MPI=OFF", "-DBUILD_TESTING=OFF", *options,
    ], text=True, capture_output=True)


def test_applications_and_examples_respect_build_testing(tmp_path):
    result = configure(tmp_path, "-DALPS_BUILD_APPLICATIONS=ON",
                       "-DALPS_BUILD_EXAMPLES=ON")
    assert result.returncode == 0, result.stdout + result.stderr
    result = subprocess.run([
        "ctest", "--test-dir", str(tmp_path), "--show-only=json-v1",
    ], check=True, text=True, capture_output=True)
    assert json.loads(result.stdout)["tests"] == []


def test_missing_blas_is_a_configuration_error(tmp_path):
    result = configure(tmp_path, "-DALPS_BUILD_APPLICATIONS=ON",
                       "-DBLA_VENDOR=Generic",
                       "-DCMAKE_DISABLE_FIND_PACKAGE_BLAS=ON")
    assert result.returncode != 0
    assert "CMAKE_DISABLE_FIND_PACKAGE_BLAS" in result.stdout + result.stderr
