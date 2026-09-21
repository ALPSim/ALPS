"""The wheel's CMake interface follows its metadata through relocation."""
import json
from pathlib import Path
import shutil
import subprocess

import pytest


CONFIG = (Path(__file__).resolve().parents[2] /
          "bindings/python/pyalps/src/pyalps/cmake/pyalpsConfig.cmake")


@pytest.fixture
def consumer(tmp_path):
    source = tmp_path / "consumer"
    source.mkdir()
    sdk = tmp_path / "sdk"
    sdk.mkdir()
    (sdk / "ALPSConfig.cmake").write_text(
        'set(ALPS_VERSION "3.0.0")\n'
        'foreach(target ALPS::alps ALPS::headers Threads::Threads)\n'
        '  if(NOT TARGET ${target})\n'
        '    add_library(${target} INTERFACE IMPORTED)\n'
        '  endif()\n'
        'endforeach()\n')
    package = tmp_path / "original" / "pyalps"
    (package / "cmake").mkdir(parents=True)
    shutil.copy2(CONFIG, package / "cmake")
    (package / "lib").mkdir()
    (package / "include").mkdir()
    (package / "lib/libalps-hashed.so.3").touch()
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
  if(EXPECT_REPAIRED)
    get_target_property(location pyalps::library0 IMPORTED_LOCATION)
    file(REAL_PATH "${pyalps_DIR}/../lib/libalps-hashed.so.3" expected)
    if(NOT location STREQUAL expected)
      message(FATAL_ERROR "Runtime did not follow the relocated package: ${location}")
    endif()
  else()
    get_target_property(links pyalps::runtime INTERFACE_LINK_LIBRARIES)
    if(NOT links STREQUAL "ALPS::alps")
      message(FATAL_ERROR "Developer runtime lost the SDK link interface")
    endif()
  endif()
endif()
''')

    def configure(*, repaired=True, version="3.0.0", libraries=None, failure=""):
        metadata = {
            "schema": 1, "alps_version": version, "repaired": repaired,
            "libraries": libraries if libraries is not None else [{"path": "lib/libalps-hashed.so.3"}],
        }
        (package / "runtime.json").write_text(json.dumps(metadata))
        relocated = tmp_path / "relocated environment" / "pyalps"
        shutil.move(package, relocated)
        result = subprocess.run([
            "cmake", "-S", str(source), "-B", str(tmp_path / "build"), "-G", "Ninja",
            "-DCMAKE_SYSTEM_NAME=Linux", f"-DALPS_DIR={sdk}",
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
