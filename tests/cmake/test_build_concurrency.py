"""DLL staging after parallel links must not write a shared directory concurrently."""

from pathlib import Path
import subprocess
import sys

import pytest


SOURCE = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("generator", ["Ninja", "Ninja Multi-Config"])
def test_windows_runtime_staging_does_not_overlap(tmp_path, generator):
    # Simulate Windows' exclusive file access on every test host. The private
    # build policy is independent of the compiler and can use native binaries.
    (tmp_path / "stage.py").write_text('''
from pathlib import Path
import sys
import time
lock = Path(sys.argv[1])
lock.mkdir()
try:
    time.sleep(0.3)
finally:
    lock.rmdir()
''')
    (tmp_path / "main.c").write_text("int main(void) { return 0; }\n")
    (tmp_path / "CMakeLists.txt").write_text(f'''
cmake_minimum_required(VERSION 3.27)
project(runtime_staging LANGUAGES C)
set(WIN32 TRUE)
set(CMAKE_JOB_POOLS parent_compile=2)
set(CMAKE_JOB_POOL_COMPILE parent_compile)
include("{SOURCE.as_posix()}/cmake/ALPSBuildConcurrency.cmake")
foreach(index RANGE 1 6)
  add_executable(program${{index}} main.c)
  add_custom_command(TARGET program${{index}} POST_BUILD
    COMMAND "{Path(sys.executable).as_posix()}" "${{CMAKE_SOURCE_DIR}}/stage.py"
      "${{CMAKE_BINARY_DIR}}/staging-lock"
    USES_TERMINAL VERBATIM)
endforeach()
''')
    build = tmp_path / "build"
    subprocess.run(["cmake", "-S", str(tmp_path), "-B", str(build),
                    "-G", generator, "-DCMAKE_BUILD_TYPE=Release"],
                   check=True, capture_output=True, text=True)
    result = subprocess.run(["cmake", "--build", str(build), "--config", "Release",
                             "--parallel", "4"], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
