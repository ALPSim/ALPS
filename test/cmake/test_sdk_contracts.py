"""Integration checks against the installed SDK selected by ALPS_DIR."""
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

SOURCE = Path(__file__).resolve().parents[2]
CONSUMER = Path(__file__).with_name("consumer")
pytestmark = pytest.mark.skipif(not os.environ.get("ALPS_DIR"), reason="requires an installed SDK")


def configure(build, *options, success=True):
    result = subprocess.run([
        "cmake", "-S", str(CONSUMER), "-B", str(build),
        "-DCMAKE_BUILD_TYPE=Release", "-DALPS_DIR=" + os.environ["ALPS_DIR"],
        *json.loads(os.environ.get("ALPS_TEST_CMAKE_ARGS", "[]")), *options,
    ], text=True, capture_output=True)
    if success:
        assert result.returncode == 0, result.stdout + result.stderr
    else:
        assert result.returncode != 0, result.stdout + result.stderr
    return result.stdout + result.stderr


def build_and_run(build):
    subprocess.run(["cmake", "--build", str(build), "--config", "Release", "--parallel", "2"], check=True)
    subprocess.run(["ctest", "--test-dir", str(build), "-C", "Release", "--output-on-failure"], check=True)


@pytest.mark.parametrize("standard", (17, 20))
def test_installed_sdk_preserves_parent_settings(tmp_path, standard):
    configure(tmp_path, f"-DCMAKE_CXX_STANDARD={standard}")
    build_and_run(tmp_path)


def test_embedded_defaults_and_mpi_isolation(tmp_path):
    # Configure only: no duplicate ALPS object files are needed for this contract.
    configure(tmp_path, f"-DALPS_SOURCE={SOURCE}")


def test_sdk_rejects_integer_abi_mismatch(tmp_path):
    output = configure(tmp_path, "-DBLA_SIZEOF_INTEGER=8", "-DEXPECT_ABI=ON", success=False)
    assert "requires BLA_SIZEOF_INTEGER=4" in output


def test_sdk_exports_installed_applications(tmp_path):
    if not (Path(os.environ["ALPS_DIR"]) / "ALPSApplicationTargets.cmake").is_file():
        pytest.skip("requires an SDK with applications")
    configure(tmp_path, "-DEXPECT_APPLICATIONS=ON")
    # Generator expressions resolve the selected configuration and suffix.
    paths = (tmp_path / "applications-Release.txt").read_text().splitlines()
    assert len(paths) == 18
    assert len(set(paths)) == 18
    assert all(Path(path).is_file() for path in paths)


def test_relocated_sdk(tmp_path):
    prefix = Path(os.environ["ALPS_DIR"]).resolve().parents[1]
    relocated = tmp_path / "relocated"
    # Only copy SDK artifacts, not application programs or dependency caches.
    libdir = os.environ.get("ALPS_TEST_INSTALL_LIBDIR", "lib")
    bindir = os.environ.get("ALPS_TEST_INSTALL_BINDIR", "bin")
    for directory in ("include", libdir, "share/alps"):
        shutil.copytree(prefix / directory, relocated / directory, symlinks=True)
    if os.name == "nt":
        (relocated / bindir).mkdir(parents=True)
        for library in (prefix / bindir).glob("*.dll"):
            shutil.copy2(library, relocated / bindir / library.name)
    build = tmp_path / "consumer"
    configure(build, f"-DALPS_DIR={relocated / 'share/alps'}")
    build_and_run(build)
