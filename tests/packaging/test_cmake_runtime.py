"""The wheel's CMake interface follows its metadata through relocation."""
import json
from pathlib import Path
import shutil
import subprocess

import pytest


CONFIG = (Path(__file__).resolve().parents[2] /
          "python/pyalps/_build_support/pyalpsConfig.cmake")


@pytest.fixture
def consumer(tmp_path):
    source = tmp_path / "consumer"
    source.mkdir()
    sdk = tmp_path / "sdk"
    sdk.mkdir()
    (sdk / "ALPSConfig.cmake").write_text(
        'set(ALPS_VERSION "3.0.0")\n'
        'set(ALPS_RUNTIME_TARGETS ALPS::alps ALPS::params ALPS::hdf5 ALPS::utilities ALPS::osiris ALPS::xml ALPS::cli)\n'
        'foreach(target ALPS::alps ALPS::params ALPS::hdf5 ALPS::utilities ALPS::osiris ALPS::xml ALPS::cli ALPS::headers Threads::Threads Boost::filesystem Boost::serialization HDF5::HDF5 Boost::regex Boost::program_options)\n'
        '  if(NOT TARGET ${target})\n'
        '    add_library(${target} INTERFACE IMPORTED)\n'
        '  endif()\n'
        'endforeach()\n'
        'set_property(TARGET ALPS::alps PROPERTY INTERFACE_LINK_LIBRARIES "ALPS::params;ALPS::hdf5;ALPS::utilities;ALPS::osiris;ALPS::xml;ALPS::cli;ALPS::headers")\n'
        'set_property(TARGET ALPS::params PROPERTY INTERFACE_LINK_LIBRARIES "ALPS::hdf5;Boost::serialization")\n'
        'set_property(TARGET ALPS::hdf5 PROPERTY INTERFACE_LINK_LIBRARIES "ALPS::utilities;HDF5::HDF5")\n'
        'set_property(TARGET ALPS::utilities PROPERTY INTERFACE_LINK_LIBRARIES "Boost::filesystem")\n'
        'set_property(TARGET ALPS::osiris PROPERTY INTERFACE_LINK_LIBRARIES "ALPS::headers;Boost::filesystem;Boost::serialization")\n'
        'set_property(TARGET ALPS::xml PROPERTY INTERFACE_LINK_LIBRARIES "ALPS::headers;Boost::filesystem;Boost::regex")\n'
        'set_property(TARGET ALPS::cli PROPERTY INTERFACE_LINK_LIBRARIES "ALPS::headers;ALPS::utilities;Boost::program_options")\n')
    (sdk / "src").mkdir()
    (sdk / "src/nb_abi.h").write_text("#  define NB_INTERNALS_VERSION 21\n")
    package = tmp_path / "original" / "pyalps"
    (package / "cmake").mkdir(parents=True)
    shutil.copy2(CONFIG, package / "cmake")
    (package / "lib").mkdir()
    (package / "include").mkdir()
    components = ("alps", "alps_params", "alps_hdf5", "alps_utilities", "alps_osiris", "alps_xml", "alps_cli")
    packaged_libraries = [{"path": f"lib/lib{component}-hashed.so.3"} for component in components]
    for library in packaged_libraries:
        (package / library["path"]).touch()
    (source / "CMakeLists.txt").write_text('''
cmake_minimum_required(VERSION 3.27)
project(wheel_runtime_contract LANGUAGES NONE)
find_package(pyalps CONFIG QUIET)
if(EXPECT_FAILURE)
  if(pyalps_FOUND OR TARGET pyalps::runtime OR TARGET pyalps::library0)
    message(FATAL_ERROR "Rejected package left a usable runtime target")
  endif()
  if(NOT pyalps_NOT_FOUND_MESSAGE MATCHES "${EXPECT_FAILURE}")
    message(FATAL_ERROR "Missing diagnostic: ${pyalps_NOT_FOUND_MESSAGE}")
  endif()
else()
  if(NOT pyalps_FOUND OR NOT TARGET pyalps::runtime)
    message(FATAL_ERROR "Runtime package was not found: ${pyalps_NOT_FOUND_MESSAGE}")
  endif()
  find_package(pyalps CONFIG REQUIRED) # Repeated discovery is harmless.
  get_target_property(headers pyalps::runtime INTERFACE_INCLUDE_DIRECTORIES)
  file(REAL_PATH "${pyalps_DIR}/../include" expected_headers)
  if(NOT headers STREQUAL expected_headers)
    message(FATAL_ERROR "Binding headers did not follow the relocated package")
  endif()
  get_target_property(links pyalps::runtime INTERFACE_LINK_LIBRARIES)
  set(components alps alps_params alps_hdf5 alps_utilities alps_osiris alps_xml alps_cli)
  set(index 0)
  foreach(component IN LISTS components)
    get_target_property(location pyalps::library${index} IMPORTED_LOCATION)
    file(REAL_PATH "${pyalps_DIR}/../lib/lib${component}-hashed.so.3" expected)
    if(NOT location STREQUAL expected OR NOT "pyalps::library${index}" IN_LIST links)
      message(FATAL_ERROR "Runtime component did not follow the relocated package: ${component}")
    endif()
    math(EXPR index "${index} + 1")
  endforeach()
  if(EXPECT_REPAIRED)
    get_target_property(location pyalps::library0 IMPORTED_LOCATION)
    file(REAL_PATH "${pyalps_DIR}/../lib/libalps-hashed.so.3" expected)
    if(NOT location STREQUAL expected)
      message(FATAL_ERROR "Runtime did not follow the relocated package: ${location}")
    endif()
  else()
    if("ALPS::utilities" IN_LIST links OR "ALPS::alps" IN_LIST links OR "ALPS::hdf5" IN_LIST links OR "ALPS::params" IN_LIST links OR "ALPS::osiris" IN_LIST links OR "ALPS::xml" IN_LIST links OR "ALPS::cli" IN_LIST links
       OR NOT "Boost::regex" IN_LIST links OR NOT "Boost::program_options" IN_LIST links
       OR NOT "Boost::filesystem" IN_LIST links OR NOT "HDF5::HDF5" IN_LIST links OR NOT "Boost::serialization" IN_LIST links)
      message(FATAL_ERROR "Developer extensions must use packaged components and SDK external dependencies")
    endif()
    get_target_property(location pyalps::library0 IMPORTED_LOCATION)
    if(NOT location MATCHES "/pyalps/lib/libalps-hashed.so.3$")
      message(FATAL_ERROR "Developer runtime must use the package-owned ALPS library")
    endif()
  endif()
endif()
''')

    def configure(*, repaired=True, version="3.0.0", libraries=None, failure="", abi="21"):
        metadata = {
            "nanobind": {"version": "2.15.0", "internals_abi": abi},
            "schema": 1, "alps_version": version, "repaired": repaired,
            "libraries": libraries if libraries is not None else packaged_libraries,
        }
        (package / "runtime.json").write_text(json.dumps(metadata))
        relocated = tmp_path / "relocated environment" / "pyalps"
        shutil.move(package, relocated)
        result = subprocess.run([
            "cmake", "-S", str(source), "-B", str(tmp_path / "build"), "-G", "Ninja",
            "-DCMAKE_SYSTEM_NAME=Linux", f"-DALPS_DIR={sdk}", f"-DNB_DIR={sdk}",
            f"-Dpyalps_DIR={relocated / 'cmake'}", f"-DEXPECT_REPAIRED={int(repaired)}",
            f"-DEXPECT_FAILURE={failure}",
        ], text=True, capture_output=True)
        assert result.returncode == 0, result.stdout + result.stderr

    return configure


def test_repaired_runtime_follows_wheel_relocation(consumer):
    consumer()


def test_developer_runtime_uses_sdk_targets(consumer):
    consumer(repaired=False)


def test_runtime_rejects_wrong_sdk_version(consumer):
    consumer(version="9.0.0", failure="ALPS SDK 9.0.0")


@pytest.mark.parametrize("libraries", [[], [{"path": "lib/missing.so"}], [{"path": "../../sdk/ALPSConfig.cmake"}]])
def test_runtime_rejects_missing_or_escaping_libraries(consumer, libraries):
    consumer(libraries=libraries, failure="no runtime libraries|Missing or invalid")


def test_runtime_rejects_incompatible_nanobind(consumer):
    consumer(abi="20", failure="does not match pyalps nanobind")
