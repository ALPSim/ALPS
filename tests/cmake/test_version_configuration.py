"""The shared release version must be strict and participate in regeneration."""
from pathlib import Path
import shutil
import subprocess

import pytest


@pytest.fixture
def version_project(tmp_path):
    source = tmp_path / "source"
    (source / "cmake").mkdir(parents=True)
    shutil.copy2(Path(__file__).resolve().parents[2] / "cmake/ALPSVersion.cmake",
                 source / "cmake/ALPSVersion.cmake")
    (source / "CMakeLists.txt").write_text(
        'cmake_minimum_required(VERSION 4.3)\n'
        'include(cmake/ALPSVersion.cmake)\n'
        'project(version_contract VERSION ${ALPS_VERSION_CORE} LANGUAGES NONE)\n'
        'file(WRITE "${CMAKE_BINARY_DIR}/version.txt" "${PROJECT_VERSION}")\n',
        encoding="utf-8",
    )
    return source, tmp_path / "build"


def configure(source, build):
    return subprocess.run(
        ["cmake", "-S", str(source), "-B", str(build), "-G", "Ninja"],
        capture_output=True, text=True,
    )


def test_version_change_reconfigures_build(version_project):
    source, build = version_project
    version = source / "ALPS_VERSION.txt"
    version.write_text("3.0.0\n", encoding="utf-8")
    result = configure(source, build)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (build / "version.txt").read_text() == "3.0.0"

    version.write_text("3.1.7\n", encoding="utf-8")
    result = subprocess.run(["cmake", "--build", str(build)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (build / "version.txt").read_text() == "3.1.7"


@pytest.mark.parametrize("version", ["3.0.0\nextra\n", "03.0.0", "3.0.0-beta.1", ""])
def test_invalid_version_is_rejected(version_project, version):
    source, build = version_project
    (source / "ALPS_VERSION.txt").write_text(version, encoding="utf-8")
    result = configure(source, build)
    assert result.returncode != 0
    assert "must contain exactly MAJOR.MINOR.PATCH" in result.stdout + result.stderr
