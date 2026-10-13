# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#
# Deprecated compatibility for alps_target_link_pyalps(<target>
# [PYTHON_EXECUTABLE <python>]). It links <target> to pyalps::runtime from the
# installed pyalps package; new projects use that package directly, see
# python/pyalps/README.md#downstream-native-extensions.

include_guard(GLOBAL)
message(DEPRECATION "UsePyALPS.cmake is deprecated: find_package(pyalps) and link pyalps::runtime")

function(alps_target_link_pyalps target)
  if(NOT TARGET "${target}")
    message(FATAL_ERROR "alps_target_link_pyalps: '${target}' is not a CMake target")
  endif()
  cmake_parse_arguments(PYALPS "" "PYTHON_EXECUTABLE" "" ${ARGN})
  if(NOT PYALPS_PYTHON_EXECUTABLE)
    if(NOT Python_EXECUTABLE)
      message(FATAL_ERROR
        "alps_target_link_pyalps requires PYTHON_EXECUTABLE or a preceding "
        "find_package(Python ... Interpreter)")
    endif()
    set(PYALPS_PYTHON_EXECUTABLE "${Python_EXECUTABLE}")
  endif()
  if(NOT TARGET pyalps::runtime)
    # pyalps checks the consumer's nanobind ABI against its own.
    if(NOT NB_DIR)
      message(FATAL_ERROR "Call find_package(nanobind CONFIG REQUIRED) before alps_target_link_pyalps")
    endif()
    execute_process(
      COMMAND "${PYALPS_PYTHON_EXECUTABLE}" -m pyalps --cmake-dir
      OUTPUT_VARIABLE pyalps_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
    find_package(pyalps CONFIG REQUIRED)
  endif()
  target_link_libraries("${target}" PRIVATE pyalps::runtime)
endfunction()
