# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT
#
# Link a downstream nanobind module to the same ALPS runtime as pyalps.
#
# Binary wheels keep libalps and its non-system dependencies in
# wheel-private directories: the ALPS libraries in pyalps/lib, and the
# dependencies the repair tools relocated in pyalps.libs (auditwheel) or
# pyalps/.dylibs (delocate). Linking a consumer module to a separately installed
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
    # A repair-tool directory is what marks this install as a relocated wheel:
    # auditwheel writes <site-packages>/pyalps.libs, delocate writes
    # pyalps/.dylibs.
    set(_pyalps_repaired_dirs "")
    foreach(_candidate IN ITEMS
        "${_pyalps_package_dir}/.dylibs"
        "${_pyalps_package_dir}/../pyalps.libs")
      if(IS_DIRECTORY "${_candidate}")
        get_filename_component(_candidate "${_candidate}" REALPATH)
        list(APPEND _pyalps_repaired_dirs "${_candidate}")
      endif()
    endforeach()

    if(_pyalps_repaired_dirs)
      # The repaired runtime spans two directories. The ALPS libraries
      # themselves are bundled into pyalps/lib by the wheel build -- one copy,
      # shared with the programs in pyalps/bin -- while the repair tools
      # vendor the external dependencies (HDF5, LAPACK, ...) into their own
      # directory. Search both.
      #
      # pyalps/lib alone must not trigger this path: a developer install has
      # that directory too, without any of the vendored dependencies the loop
      # below insists on finding, and such an install is meant to keep using
      # the ordinary ALPSConfig.cmake library paths.
      if(IS_DIRECTORY "${_pyalps_package_dir}/lib")
        get_filename_component(_pyalps_bundled_libs
          "${_pyalps_package_dir}/lib" REALPATH)
        list(APPEND _pyalps_private_runtime "${_pyalps_bundled_libs}")
      endif()
      list(APPEND _pyalps_private_runtime ${_pyalps_repaired_dirs})
      list(REMOVE_DUPLICATES _pyalps_private_runtime)
    endif()
  endif()

  if(_pyalps_private_runtime)
    set(_pyalps_private_libraries "")
    foreach(_library IN LISTS _pyalps_link_libraries)
      # Installed ALPSConfig files may record dependencies as bare linker
      # names (alps), linker flags (-lalps), or absolute paths
      # (/usr/lib64/liblapack.so). Wheel repair tools rename all three forms
      # to a private file such as liblapack-<hash>.so. Normalize the original
      # entry to its library stem before looking up that repaired file.
      set(_library_stem "${_library}")
      if(IS_ABSOLUTE "${_library_stem}")
        get_filename_component(_library_stem "${_library_stem}" NAME)
      endif()
      string(REGEX REPLACE "^-l" "" _library_stem "${_library_stem}")
      string(REGEX REPLACE "^lib" "" _library_stem "${_library_stem}")
      string(REGEX REPLACE "\\.so(\\.[0-9]+)*$" "" _library_stem
        "${_library_stem}")
      string(REGEX REPLACE "(\\.[0-9]+)*\\.dylib$" "" _library_stem
        "${_library_stem}")
      string(REGEX REPLACE "\\.a$" "" _library_stem "${_library_stem}")
      if(_library_stem STREQUAL "hdf5-shared")
        set(_library_stem "hdf5")
      endif()

      set(_matches "")
      foreach(_runtime_dir IN LISTS _pyalps_private_runtime)
        file(GLOB _runtime_dir_matches LIST_DIRECTORIES FALSE
          "${_runtime_dir}/lib${_library_stem}.so*"
          "${_runtime_dir}/lib${_library_stem}-*.so*"
          "${_runtime_dir}/lib${_library_stem}.dylib"
          "${_runtime_dir}/lib${_library_stem}.*.dylib"
          "${_runtime_dir}/lib${_library_stem}-*.dylib")
        list(APPEND _matches ${_runtime_dir_matches})
      endforeach()
      list(REMOVE_DUPLICATES _matches)
      list(LENGTH _matches _match_count)
      if(NOT _match_count EQUAL 1)
        message(FATAL_ERROR
          "pyalps uses a private wheel runtime, but exactly one bundled "
          "${_library_stem} library (from '${_library}') was expected in "
          "${_pyalps_private_runtime}; found: ${_matches}")
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
        if(NOT _install_name STREQUAL "@rpath/${_runtime_name}")
          add_custom_command(TARGET "${target}" POST_BUILD
            COMMAND "${CMAKE_INSTALL_NAME_TOOL}" -change
              "${_install_name}" "@rpath/${_runtime_name}"
              "$<TARGET_FILE:${target}>"
            VERBATIM)
        endif()
      endif()
    endforeach()
    list(REMOVE_DUPLICATES _pyalps_private_libraries)

    set(_pyalps_link_libraries ${_pyalps_private_libraries})
    set(_pyalps_runtime_paths ${_pyalps_private_runtime})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
      # GNU DT_RUNPATH is searched only for direct dependencies. Wheel
      # libraries such as LAPACK can themselves depend on relocated runtime
      # libraries (for example auditwheel's libgfortran copy), so emit the
      # transitive DT_RPATH tag for this consumer module instead.
      target_link_options("${target}" PRIVATE "LINKER:--disable-new-dtags")
    endif()
    string(REPLACE ";" ", " _pyalps_runtime_report "${_pyalps_private_runtime}")
    message(STATUS
      "${target}: using pyalps wheel runtime at ${_pyalps_runtime_report}")
  else()
    target_link_directories("${target}" PRIVATE ${_pyalps_runtime_paths})
  endif()

  target_link_libraries("${target}" PRIVATE ${_pyalps_link_libraries})
  set_property(TARGET "${target}" APPEND PROPERTY
    BUILD_RPATH ${_pyalps_runtime_paths})
  set_property(TARGET "${target}" APPEND PROPERTY
    INSTALL_RPATH ${_pyalps_runtime_paths})
endfunction()
