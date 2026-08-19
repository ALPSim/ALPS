# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT
#
# Link a downstream nanobind module to the same ALPS runtime as pyalps.
#
# Binary wheels relocate libalps and its non-system dependencies into a
# wheel-private directory. Linking a consumer module to a separately installed
# ALPS/HDF5 stack is unsafe: objects such as hdf5::archive carry handles that
# are valid only in the HDF5 image that created them. This helper discovers a
# repaired wheel's private runtime and links the target to those exact files.
# Source/developer installs without relocated libraries keep using the normal
# ALPSConfig.cmake library paths.

include_guard(GLOBAL)

function(alps_target_link_pyalps target)
  if(NOT TARGET "${target}")
    message(FATAL_ERROR
      "alps_target_link_pyalps: '${target}' is not a CMake target")
  endif()

  cmake_parse_arguments(PYALPS "" "PYTHON_EXECUTABLE" "" ${ARGN})
  if(NOT PYALPS_PYTHON_EXECUTABLE)
    if(Python_EXECUTABLE)
      set(PYALPS_PYTHON_EXECUTABLE "${Python_EXECUTABLE}")
    else()
      message(FATAL_ERROR
        "alps_target_link_pyalps requires PYTHON_EXECUTABLE or a preceding "
        "find_package(Python ... Interpreter)")
    endif()
  endif()

  set(_pyalps_link_libraries ${ALPS_LIBRARIES})
  list(TRANSFORM _pyalps_link_libraries REPLACE "^hdf5-shared$" "hdf5")

  set(_pyalps_runtime_paths ${ALPS_LIBRARY_DIRS})
  if(ALPS_HDF5_INCLUDE_DIR)
    get_filename_component(_pyalps_hdf5_prefix
      "${ALPS_HDF5_INCLUDE_DIR}" DIRECTORY)
    list(APPEND _pyalps_runtime_paths "${_pyalps_hdf5_prefix}/lib")
  endif()

  # Find the package without importing it. Importing an extension while CMake
  # configures would load its runtime only in this short-lived child process.
  execute_process(
    COMMAND "${PYALPS_PYTHON_EXECUTABLE}" -c
      "import importlib.util, pathlib; s=importlib.util.find_spec('pyalps'); print(pathlib.Path(next(iter(s.submodule_search_locations))).resolve() if s and s.submodule_search_locations else '')"
    OUTPUT_VARIABLE _pyalps_package_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _pyalps_location_result
    ERROR_QUIET)

  if(_pyalps_location_result EQUAL 0 AND _pyalps_package_dir)
    set(_pyalps_runtime_candidates
      "${_pyalps_package_dir}/.dylibs"
      "${_pyalps_package_dir}/../pyalps.libs")
    foreach(_candidate IN LISTS _pyalps_runtime_candidates)
      if(IS_DIRECTORY "${_candidate}")
        get_filename_component(_pyalps_private_runtime "${_candidate}" REALPATH)
        break()
      endif()
    endforeach()
  endif()

  if(_pyalps_private_runtime)
    set(_pyalps_private_libraries "")
    foreach(_library IN LISTS _pyalps_link_libraries)
      file(GLOB _matches LIST_DIRECTORIES FALSE
        "${_pyalps_private_runtime}/lib${_library}.so*"
        "${_pyalps_private_runtime}/lib${_library}-*.so*"
        "${_pyalps_private_runtime}/lib${_library}.dylib"
        "${_pyalps_private_runtime}/lib${_library}.*.dylib"
        "${_pyalps_private_runtime}/lib${_library}-*.dylib")
      list(REMOVE_DUPLICATES _matches)
      list(LENGTH _matches _match_count)
      if(NOT _match_count EQUAL 1)
        message(FATAL_ERROR
          "pyalps uses a private wheel runtime, but exactly one bundled "
          "${_library} library was expected in ${_pyalps_private_runtime}; "
          "found: ${_matches}")
      endif()
      list(GET _matches 0 _match)
      list(APPEND _pyalps_private_libraries "${_match}")

      # delocate gives copied dylibs collision-resistant /DLC install names.
      # Those names are intentionally not real paths, so rewrite this target's
      # references to @rpath and point that rpath at the wheel directory. This
      # also permits importing the consumer module before importing pyalps.
      if(APPLE)
        if(NOT CMAKE_OTOOL)
          find_program(CMAKE_OTOOL otool REQUIRED)
        endif()
        if(NOT CMAKE_INSTALL_NAME_TOOL)
          find_program(CMAKE_INSTALL_NAME_TOOL install_name_tool REQUIRED)
        endif()
        execute_process(
          COMMAND "${CMAKE_OTOOL}" -D "${_match}"
          OUTPUT_VARIABLE _install_names
          OUTPUT_STRIP_TRAILING_WHITESPACE
          COMMAND_ERROR_IS_FATAL ANY)
        string(REGEX MATCHALL "[^\r\n]+" _install_name_lines
          "${_install_names}")
        list(LENGTH _install_name_lines _install_name_line_count)
        if(_install_name_line_count LESS 2)
          message(FATAL_ERROR "Could not read the install name of ${_match}")
        endif()
        list(GET _install_name_lines 1 _install_name)
        string(STRIP "${_install_name}" _install_name)
        get_filename_component(_runtime_name "${_match}" NAME)
        add_custom_command(TARGET "${target}" POST_BUILD
          COMMAND "${CMAKE_INSTALL_NAME_TOOL}" -change
            "${_install_name}" "@rpath/${_runtime_name}"
            "$<TARGET_FILE:${target}>"
          VERBATIM)
      endif()
    endforeach()

    set(_pyalps_link_libraries ${_pyalps_private_libraries})
    set(_pyalps_runtime_paths "${_pyalps_private_runtime}")
    message(STATUS
      "${target}: using pyalps wheel runtime at ${_pyalps_private_runtime}")
  else()
    target_link_directories("${target}" PRIVATE ${_pyalps_runtime_paths})
  endif()

  target_link_libraries("${target}" PRIVATE ${_pyalps_link_libraries})
  set_property(TARGET "${target}" APPEND PROPERTY
    BUILD_RPATH ${_pyalps_runtime_paths})
  set_property(TARGET "${target}" APPEND PROPERTY
    INSTALL_RPATH ${_pyalps_runtime_paths})
endfunction()
