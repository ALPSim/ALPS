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


@pytest.mark.parametrize("unavailable,shared", [("header", "OFF"), ("symbols", "ON")])
def test_utilities_without_platform_stacktrace(tmp_path, monkeypatch, unavailable, shared):
    # Reproduce missing execinfo headers and headers whose symbols need an
    # unavailable library, without relying on CI's wheel-only compiler flags.
    monkeypatch.delenv("CXXFLAGS", raising=False)
    header = "#error execinfo is unavailable\n"
    if unavailable == "symbols":
        header = ('extern "C" int alps_missing_backtrace(void **, int);\n'
                  'extern "C" char **alps_missing_backtrace_symbols(void *const *, int);\n'
                  '#define backtrace alps_missing_backtrace\n'
                  '#define backtrace_symbols alps_missing_backtrace_symbols\n')
    (tmp_path / "execinfo.h").write_text(header)
    (tmp_path / "main.cpp").write_text(
        '#include <alps/ngs/stacktrace.hpp>\n'
        'int main() { return !alps::ngs::stacktrace().empty(); }\n')
    (tmp_path / "CMakeLists.txt").write_text(
        'cmake_minimum_required(VERSION 3.27...4.3)\n'
        'project(stacktrace_contract LANGUAGES C CXX)\n'
        'set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")\n'
        # A parent using compile-only probes must not hide missing link symbols.
        'set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)\n'
        f'add_subdirectory("{SOURCE.as_posix()}" alps EXCLUDE_FROM_ALL)\n'
        'add_executable(stacktrace_contract main.cpp)\n'
        'target_link_libraries(stacktrace_contract PRIVATE ALPS::utilities)\n'
        'enable_testing()\n'
        'add_test(NAME stacktrace_contract COMMAND stacktrace_contract)\n')
    build = tmp_path / "build"
    result = subprocess.run([
        "cmake", "-S", str(tmp_path), "-B", str(build),
        *json.loads(os.environ.get("ALPS_TEST_CMAKE_ARGS", "[]")),
        "-DCMAKE_BUILD_TYPE=Release", f"-DBUILD_SHARED_LIBS={shared}",
        f'-DCMAKE_CXX_FLAGS=-I"{tmp_path.as_posix()}"',
    ], text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    build_and_run(build)


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


def test_installed_sdk_preserves_public_headers():
    prefix = Path(os.environ["ALPS_DIR"]).resolve().parents[1]
    baseline = json.loads(Path(__file__).with_name("preserved_public_headers.json").read_text())
    missing = [name for name in baseline["headers"] if not (prefix / "include" / name).is_file()]
    assert not missing, f"Public headers lost from {baseline['source_commit']}: {missing}"


@pytest.mark.parametrize("tutorial,target", [
    ("08-alpsize/01-cmake", "hello"),
    ("09-code/02-c++", "ising"),
])
def test_tutorials_build_against_exported_sdk(tmp_path, tutorial, target):
    """Build the actual tutorial consumers through their documented CMake path."""
    build = tmp_path / "build"
    subprocess.run([
        "cmake", "-S", str(SOURCE / "tutorials" / tutorial), "-B", str(build),
        "-DCMAKE_BUILD_TYPE=Release", "-DALPS_DIR=" + os.environ["ALPS_DIR"],
        *json.loads(os.environ.get("ALPS_TEST_CMAKE_ARGS", "[]")),
    ], check=True)
    subprocess.run(["cmake", "--build", str(build), "--config", "Release",
                    "--target", target, "--parallel", "2"], check=True)
    executable = build / (f"Release/{target}.exe" if os.name == "nt" else target)
    assert executable.is_file()
    if target == "hello":
        result = subprocess.run([str(executable)], check=True, capture_output=True, text=True)
        assert result.stdout.strip() == "hello, world"


@pytest.mark.skipif(os.name == "nt" or not shutil.which("make"), reason="requires Unix make")
def test_make_intro_builds_without_sdk(tmp_path):
    # This lesson is ordinary hello-world C++, preceding ALPS integration.
    source = SOURCE / "tutorials/08-alpsize/00-make"
    for name in ("hello.C", "Makefile"):
        shutil.copy2(source / name, tmp_path / name)
    environment = os.environ.copy()
    for name in ("ALPS_HOME", "ALPS_ROOT", "ALPS_DIR", "CMAKE_PREFIX_PATH",
                 "CPPFLAGS", "LDFLAGS", "LDLIBS"):
        environment.pop(name, None)
    subprocess.run(["make", "-C", str(tmp_path)], check=True, env=environment)
    result = subprocess.run([str(tmp_path / "hello")], check=True,
                            capture_output=True, text=True, env=environment)
    assert result.stdout.strip() == "hello, world"
