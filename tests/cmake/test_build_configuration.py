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
        "-DALPS_ENABLE_MPI=OFF", "-DALPS_BUILD_TESTING=OFF", *options,
    ], text=True, capture_output=True)


def test_applications_respect_build_testing(tmp_path):
    result = configure(tmp_path, "-DALPS_BUILD_APPLICATIONS=ON")
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


@pytest.mark.parametrize("option, diagnostic", [
    ("BLA_SIZEOF_INTEGER=8", "requires BLA_SIZEOF_INTEGER=4"),
    ("BLA_SIZEOF_INTEGER=ANY", "requires BLA_SIZEOF_INTEGER=4"),
    ("BIND_FORTRAN_LOWERCASE=ON", "lowercase underscore symbols"),
])
def test_unsupported_numerical_abi_is_rejected(tmp_path, option, diagnostic):
    result = configure(tmp_path, "-DALPS_BUILD_APPLICATIONS=OFF", f"-D{option}")
    assert result.returncode != 0
    assert diagnostic in result.stdout + result.stderr


def test_embedded_in_source_build_is_rejected_before_project(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    # The guard must run before any project setup or dependency discovery.
    shutil.copy2(SOURCE / "CMakeLists.txt", source)
    (tmp_path / "CMakeLists.txt").write_text(
        "cmake_minimum_required(VERSION 3.27...4.3)\n"
        "project(parent LANGUAGES NONE)\n"
        'add_subdirectory(source "${CMAKE_CURRENT_SOURCE_DIR}/source")\n')
    result = subprocess.run([
        "cmake", "-S", str(tmp_path), "-B", str(tmp_path / "build"),
    ], text=True, capture_output=True)
    assert result.returncode != 0
    assert "Use an out-of-source build" in result.stdout + result.stderr


def test_hdf5_runtime_paths_follow_imported_configurations(tmp_path):
    build = tmp_path / "build"
    # Dependencies extracted inside the build tree are omitted from CMake's
    # automatic install RPATH. Use separate directories to expose flattening.
    provider = build / "hdf5"
    (provider / "include").mkdir(parents=True)
    (provider / "hdf5-config.cmake").write_text('''
set(HDF5_VERSION 1.14.6)
set(HDF5_ENABLE_PARALLEL OFF)
set(HDF5_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/include")
add_library(hdf5::hdf5-shared SHARED IMPORTED)
set_target_properties(hdf5::hdf5-shared PROPERTIES
  IMPORTED_CONFIGURATIONS "DEBUG;RELEASE"
  INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIR}")
foreach(config IN ITEMS Debug Release)
  string(TOUPPER "${config}" upper)
  set(directory "${CMAKE_CURRENT_LIST_DIR}/${config}")
  file(MAKE_DIRECTORY "${directory}")
  set(library "${directory}/${CMAKE_SHARED_LIBRARY_PREFIX}hdf5${CMAKE_SHARED_LIBRARY_SUFFIX}")
  file(WRITE "${library}" "")
  set_target_properties(hdf5::hdf5-shared PROPERTIES
    IMPORTED_LOCATION_${upper} "${library}")
  if(WIN32)
    file(WRITE "${directory}/hdf5.lib" "")
    set_target_properties(hdf5::hdf5-shared PROPERTIES
      IMPORTED_IMPLIB_${upper} "${directory}/hdf5.lib")
  endif()
endforeach()
''')
    capture = tmp_path / "capture.cmake"
    capture.write_text(
        'include(CMakePackageConfigHelpers)\n'
        f'write_basic_package_version_file("{provider.as_posix()}/hdf5-config-version.cmake"\n'
        '  VERSION 1.14.6 COMPATIBILITY AnyNewerVersion)\n'
        'file(GENERATE OUTPUT "${CMAKE_BINARY_DIR}/runtime-$<CONFIG>.txt"\n'
        '  CONTENT "$<TARGET_GENEX_EVAL:alps_hdf5,$<TARGET_PROPERTY:alps_hdf5,INSTALL_RPATH>>;'
        '$<TARGET_RUNTIME_DLLS:alps_hdf5>" TARGET alps_hdf5)\n'
        'file(GENERATE OUTPUT "${CMAKE_BINARY_DIR}/program-runtime-$<CONFIG>.txt"\n'
        '  CONTENT "$<TARGET_GENEX_EVAL:dmrg,$<TARGET_PROPERTY:dmrg,INSTALL_RPATH>>"\n'
        '  TARGET dmrg)\n')
    result = subprocess.run([
        "cmake", "-S", str(SOURCE), "-B", str(build),
        *json.loads(os.environ.get("ALPS_TEST_CMAKE_ARGS", "[]")),
        "-G", "Ninja Multi-Config", "-DCMAKE_CONFIGURATION_TYPES=Debug;Release",
        "-DALPS_BUILD_TESTING=OFF", "-DALPS_BUILD_APPLICATIONS=ON", "-DALPS_ENABLE_MPI=OFF",
        "-DBUILD_SHARED_LIBS=ON", "-DHDF5_USE_STATIC_LIBRARIES=OFF",
        # This contract deliberately selects a fixture provider instead of
        # the caller's HDF5 wrapper.
        "-UHDF5_ROOT", "-DHDF5_NO_FIND_PACKAGE_CONFIG_FILE=OFF",
        "-DCMAKE_FIND_PACKAGE_PREFER_CONFIG=ON", f"-DHDF5_DIR={provider}",
        f"-DCMAKE_PROJECT_alps_INCLUDE={capture}",
    ], text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    for config, other in (("Debug", "Release"), ("Release", "Debug")):
        runtime = (build / f"runtime-{config}.txt").read_text().split(";")
        suffix = "/hdf5.dll" if os.name == "nt" else ""
        assert f"{provider.as_posix()}/{config}{suffix}" in runtime
        assert f"{provider.as_posix()}/{other}{suffix}" not in runtime
        if os.name != "nt":
            program_runtime = (build / f"program-runtime-{config}.txt").read_text().split(";")
            assert f"{provider.as_posix()}/{config}" in program_runtime
            assert f"{provider.as_posix()}/{other}" not in program_runtime


def test_tutorials_are_an_explicit_install_component(tmp_path):
    build = tmp_path / "build"
    install = tmp_path / "install"
    (tmp_path / "CMakeLists.txt").write_text(
        "cmake_minimum_required(VERSION 3.27...4.3)\n"
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
    assert (tutorials / "03-mc/01-autocorrelations/parm1a").is_file()
    assert (tutorials / "02-ed/06-fulldiag/parm6a").is_file()
    assert (tutorials / "10-ngs/1_accumulator_only/CMakeLists.txt").is_file()
    assert (tutorials / "README.md").is_file()
    assert (tutorials / "00-examples/README.md").is_file()
    assert (tutorials / "00-examples/CMakeLists.txt").is_file()
    assert (tutorials / "00-examples/alea/testfile.h5").is_file()
    assert (tutorials / "00-examples/parapack/exchange/params-ising").is_file()
    assert (tutorials / "00-examples/parapack/loop/params_disorder").is_file()
    assert (tutorials / "00-examples/parapack/wanglandau/params_learn").is_file()
    assert (tutorials / "00-examples/parapack/multiple/ising.op-2").is_file()
    assert (install / "share/alps/cmake/ALPSTesting.cmake").is_file()
    assert (install / "share/alps/cmake/run_test.cmake").is_file()
    assert not (tutorials / "10-ngs/5_export_python").exists()
