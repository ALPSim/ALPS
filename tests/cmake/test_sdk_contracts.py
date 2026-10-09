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


def build_and_run(build, *, environment=None):
    subprocess.run(["cmake", "--build", str(build), "--config", "Release", "--parallel", "2"], check=True)
    subprocess.run([shutil.which("ctest"), "--test-dir", str(build), "-C", "Release", "--output-on-failure"],
                   check=True, env=environment)


@pytest.mark.parametrize("standard", (17, 20))
def test_installed_sdk_preserves_parent_settings(tmp_path, standard):
    configure(tmp_path, f"-DCMAKE_CXX_STANDARD={standard}",
              "-DCMAKE_FIND_PACKAGE_PREFER_CONFIG=ON")
    build_and_run(tmp_path)


def test_embedded_defaults_and_mpi_isolation(tmp_path):
    # Configure only: no duplicate ALPS object files are needed for this contract.
    configure(tmp_path, f"-DALPS_SOURCE={SOURCE}",
              "-DBUILD_SHARED_LIBS=OFF",
              "-DCMAKE_POSITION_INDEPENDENT_CODE=OFF",
              *(f"-DCMAKE_{kind}_OUTPUT_DIRECTORY={tmp_path.as_posix()}/{kind.lower()}"
                for kind in ("RUNTIME", "LIBRARY", "ARCHIVE")))

def test_components_link_without_building_the_core_runtime(tmp_path):
    # Compile only the extracted components, including the static export modes.
    configure(tmp_path, f"-DALPS_SOURCE={SOURCE}", "-DBUILD_SHARED_LIBS=OFF",
              "-DALPS_BUILD_APPLICATIONS=ON", "-DEXPECT_SOLVERS=ON",
              "-DCMAKE_POSITION_INDEPENDENT_CODE=OFF",
              *(f"-DCMAKE_{kind}_OUTPUT_DIRECTORY={tmp_path.as_posix()}/{kind.lower()}"
                for kind in ("RUNTIME", "LIBRARY", "ARCHIVE")))
    subprocess.run([
        "cmake", "--build", str(tmp_path), "--config", "Release",
        "--target", "utilities_contract", "hdf5_contract", "params_contract", "osiris_contract",
        "xml_contract", "cli_contract", "numeric_contract", "numeric_io_contract", "numeric_xml_contract",
        "maxent_independent_contract", "maxent", "--parallel", "2",
    ], check=True)
    core_runtime = Path((tmp_path / "core-runtime-Release.txt").read_text().strip())
    assert not core_runtime.exists(), f"Extracted components built the core runtime: {core_runtime}"
    subprocess.run([
        "ctest", "--test-dir", str(tmp_path), "-C", "Release", "--output-on-failure",
        "-R", "^(utilities|hdf5|params|osiris|xml|cli|numeric|numeric_io|numeric_xml|maxent_independent)_contract$", "--no-tests=error",
    ], check=True)


@pytest.mark.parametrize("source", [
    "tutorials/00-examples",
])
@pytest.mark.parametrize("testing", ["ON", "OFF"])
def test_standalone_examples_respect_build_testing(tmp_path, source, testing):
    result = subprocess.run([
        "cmake", "-S", str(SOURCE / source), "-B", str(tmp_path),
        "-DALPS_DIR=" + os.environ["ALPS_DIR"],
        *json.loads(os.environ.get("ALPS_TEST_CMAKE_ARGS", "[]")),
        f"-DALPS_BUILD_TESTING={testing}",
    ], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    result = subprocess.run([
        "ctest", "--test-dir", str(tmp_path), "--show-only=json-v1",
    ], check=True, capture_output=True, text=True)
    assert bool(json.loads(result.stdout)["tests"]) == (testing == "ON")


def test_sdk_rejects_integer_abi_mismatch(tmp_path):
    output = configure(tmp_path, "-DBLA_SIZEOF_INTEGER=8", "-DEXPECT_ABI=ON", success=False)
    assert "requires BLA_SIZEOF_INTEGER=4" in output


@pytest.mark.parametrize("tutorial,arguments", [
    ("09-code/06-mcmain-c++", []),
    ("10-ngs/1_accumulator_only", ["10"]),
])
def test_tutorial_text_parameters(tmp_path, tutorial, arguments):
    """Build actual callers and verify text is not treated as an HDF5 archive."""
    h5py = pytest.importorskip("h5py")
    build = tmp_path / "build"
    subprocess.run([
        "cmake", "-S", str(SOURCE / "tutorials" / tutorial), "-B", str(build),
        "-DCMAKE_BUILD_TYPE=Release", "-DALPS_DIR=" + os.environ["ALPS_DIR"],
        *json.loads(os.environ.get("ALPS_TEST_CMAKE_ARGS", "[]")),
    ], check=True)
    subprocess.run(["cmake", "--build", str(build), "--config", "Release", "--parallel", "2"], check=True)
    parameters = tmp_path / "input.txt"
    parameters.write_text("L=4; T=2.0; THERMALIZATION=2; SWEEPS=16; SEED=42;\n")
    executable = build / ("Release/ising.exe" if os.name == "nt" else "ising")
    result = subprocess.run([str(executable), *arguments, str(parameters)],
                            cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    with h5py.File(tmp_path / "input.out.h5") as archive:
        # Legacy text parameters retain expression strings in their archives.
        assert int(archive["parameters/L"][()]) == 4
        assert float(archive["parameters/T"][()]) == 2.
        assert archive["simulation/results/Energy/count"][()] > 0


def test_sdk_preserves_missing_dependency_diagnostic(tmp_path):
    configure(tmp_path, "-DEXPECT_MISSING_DEPENDENCY=ON", "-DCMAKE_DISABLE_FIND_PACKAGE_Boost=ON")


def test_sdk_exports_installed_applications(tmp_path):
    if not (Path(os.environ["ALPS_DIR"]) / "ALPSApplicationTargets.cmake").is_file():
        pytest.skip("requires an SDK with applications")
    configure(tmp_path, "-DEXPECT_APPLICATIONS=ON")
    # Generator expressions resolve the selected configuration and suffix.
    paths = (tmp_path / "applications-Release.txt").read_text().splitlines()
    assert len(paths) == 23
    assert len(set(paths)) == 23
    assert all(Path(path).is_file() for path in paths)


def test_sdk_exports_solver_libraries(tmp_path):
    if not (Path(os.environ["ALPS_DIR"]) / "ALPSApplicationTargets.cmake").is_file():
        pytest.skip("requires an SDK with solvers")
    configure(tmp_path, "-DEXPECT_SOLVERS=ON")
    build_and_run(tmp_path)


def test_relocated_sdk(tmp_path):
    # Unix SDKs intentionally use externally installed dependencies; deployment
    # without those dependencies is tested against the repaired wheels in CI.
    prefix = Path(os.environ["ALPS_DIR"]).resolve().parents[1]
    relocated = tmp_path / "relocated"
    # Preserve the installed layout, including platforms that use lib64.
    # Application programs are unnecessary for this library-consumer check.
    bindir = os.environ.get("ALPS_TEST_INSTALL_BINDIR", "bin")
    shutil.copytree(prefix, relocated, symlinks=True,
                    ignore=lambda directory, names: [bindir] if Path(directory) == prefix else [])
    if os.name == "nt":
        (relocated / bindir).mkdir(parents=True)
        for library in (prefix / bindir).glob("*.dll"):
            shutil.copy2(library, relocated / bindir / library.name)
    build = tmp_path / "consumer"
    configure(build, f"-DALPS_DIR={relocated / 'share/alps'}")
    environment = os.environ.copy()
    for name in ("LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "DYLD_FALLBACK_LIBRARY_PATH",
                 "ALPS_ROOT", "ALPS_XML_PATH"):
        environment.pop(name, None)
    environment["PATH"] = (str(Path(os.environ["SystemRoot"]) / "System32")
                           if os.name == "nt" else "/usr/bin:/bin")
    build_and_run(build, environment=environment)
