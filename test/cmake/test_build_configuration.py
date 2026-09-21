"""Configure-only contracts: exercise build choices without duplicate objects."""
import json
import os
from pathlib import Path
import shutil
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


def test_archive_requires_applications(tmp_path):
    result = configure(tmp_path, "-DALPS_BUILD_ARCHIVE=ON",
                       "-DALPS_BUILD_APPLICATIONS=OFF")
    assert result.returncode != 0
    assert "ALPS_BUILD_ARCHIVE requires ALPS_BUILD_APPLICATIONS=ON" in (
        result.stdout + result.stderr)


def test_embedded_in_source_build_is_rejected_before_project(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    # The guard must run before any project setup or dependency discovery.
    shutil.copy2(SOURCE / "CMakeLists.txt", source)
    (tmp_path / "CMakeLists.txt").write_text(
        "cmake_minimum_required(VERSION 3.27)\n"
        "project(parent LANGUAGES NONE)\n"
        'add_subdirectory(source "${CMAKE_CURRENT_SOURCE_DIR}/source")\n')
    result = subprocess.run([
        "cmake", "-S", str(tmp_path), "-B", str(tmp_path / "build"),
    ], text=True, capture_output=True)
    assert result.returncode != 0
    assert "Use an out-of-source build" in result.stdout + result.stderr


def test_tutorials_are_an_explicit_install_component(tmp_path):
    build = tmp_path / "build"
    install = tmp_path / "install"
    (tmp_path / "CMakeLists.txt").write_text(
        "cmake_minimum_required(VERSION 3.27)\n"
        "project(tutorial_install LANGUAGES NONE)\n"
        "set(CMAKE_INSTALL_DATADIR share)\n"
        f'add_subdirectory("{SOURCE.as_posix()}/tutorials" tutorials)\n')
    subprocess.run([
        "cmake", "-S", str(tmp_path), "-B", str(build),
        f"-DCMAKE_INSTALL_PREFIX={install}",
    ], check=True, capture_output=True, text=True)
    command = ["cmake", "--install", str(build)]
    subprocess.run(command, check=True, capture_output=True, text=True)
    tutorials = install / "share/alps/tutorials"
    assert not tutorials.exists()
    subprocess.run(command + ["--component", "tutorials"],
                   check=True, capture_output=True, text=True)
    assert (tutorials / "ngs/1_accumulator_only/CMakeLists.txt").is_file()
    assert not (tutorials / "ngs/5_export_python").exists()
