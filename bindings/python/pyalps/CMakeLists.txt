# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.27)
# The installed SDK supplies one runtime directory. Per-module vcpkg copies
# would introduce duplicate DLL paths when CMake collects wheel dependencies.
set(VCPKG_APPLOCAL_DEPS OFF)
project(pyalps LANGUAGES C CXX)

option(PYALPS_BUILD_APPLICATIONS "Build optional ALPS application bindings" ON)
option(PYALPS_BUNDLE_APPLICATIONS "Bundle the ALPS application executables into pyalps/bin" ON)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(PYALPS_BUNDLE_APPLICATIONS)
  find_package(ALPS REQUIRED CONFIG COMPONENTS applications)
else()
  find_package(ALPS REQUIRED CONFIG)
endif()
get_target_property(_alps_library_type ALPS::alps TYPE)
if(NOT _alps_library_type STREQUAL "SHARED_LIBRARY")
  message(
    FATAL_ERROR "pyalps requires a shared ALPS SDK so every extension uses the same C++ runtime")
endif()
find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_VARIABLE _nanobind_cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
list(PREPEND CMAKE_PREFIX_PATH "${_nanobind_cmake_dir}")
find_package(nanobind 2.10 CONFIG REQUIRED)

get_filename_component(_repo_root "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
set(_alps_source_root "${_repo_root}")
if(NOT EXISTS "${_alps_source_root}/tool/maxent.cpp")
  set(_alps_source_root "${CMAKE_CURRENT_SOURCE_DIR}/_vendor")
endif()
set(_bindings "${CMAKE_CURRENT_SOURCE_DIR}/cpp")

# Policy (see README "Free-threading and stable-ABI policy"): do NOT
# add FREE_THREADED (libalps relies on the GIL as its lock around
# shared state) and do not add STABLE_ABI without revisiting the wheel
# matrix — per-version wheels are deliberate.
# NB_STATIC keeps nanobind's core private to each extension. On musllinux,
# MUSL_DYNAMIC_LIBCPP prevents nanobind from also embedding libstdc++ and
# libgcc into every module; auditwheel vendors those shared runtimes once.
function(pyalps_add_module target)
  nanobind_add_module(${target} NB_STATIC MUSL_DYNAMIC_LIBCPP ${ARGN})
  list(APPEND _pyalps_targets ${target})
  set(_pyalps_targets "${_pyalps_targets}" PARENT_SCOPE)
endfunction()

pyalps_add_module(pyalea_c "${_bindings}/pyalea.cpp")
pyalps_add_module(pymcdata_c "${_bindings}/pymcdata.cpp")
pyalps_add_module(pytools_c "${_bindings}/pytools.cpp")
pyalps_add_module(pyngsparams_c "${_bindings}/ngs/params.cpp")
pyalps_add_module(pyngshdf5_c "${_bindings}/ngs/hdf5.cpp")
pyalps_add_module(pyngsbase_c "${_bindings}/ngs/mcbase.cpp")
pyalps_add_module(pyngsobservable_c "${_bindings}/ngs/observable.cpp")
pyalps_add_module(pyngsobservables_c "${_bindings}/ngs/observables.cpp")
pyalps_add_module(pyngsresult_c "${_bindings}/ngs/result.cpp")
pyalps_add_module(pyngsresults_c "${_bindings}/ngs/results.cpp")
pyalps_add_module(pyngsapi_c "${_bindings}/ngs/api.cpp")
pyalps_add_module(pyngsrandom01_c "${_bindings}/ngs/random01.cpp")
pyalps_add_module(pyngsaccumulator_c "${_bindings}/ngs/accumulator.cpp")

if(PYALPS_BUILD_APPLICATIONS)
  if(NOT EXISTS "${_alps_source_root}/tool/maxent.cpp")
    message(
      FATAL_ERROR "PYALPS_BUILD_APPLICATIONS requires an ALPS source checkout or pyalps sdist. "
                  "Configure with -DPYALPS_BUILD_APPLICATIONS=OFF for the core-only package.")
  endif()

  set(_dmft "${_alps_source_root}/applications/dmft/qmc")

  include("${_alps_source_root}/tool/maxent_sources.cmake")
  include("${_dmft}/sources.cmake")
  pyalps_add_module(maxent_c "${_alps_source_root}/tool/maxent.cpp" ${ALPS_MAXENT_SOURCES})
  pyalps_add_module(cthyb "${_dmft}/hybridization/hybmain.cpp" ${ALPS_CTHYB_SOURCES})
  pyalps_add_module(ctint "${_dmft}/interaction_expansion2/main.cpp" ${ALPS_CTINT_SOURCES})

  foreach(_target IN ITEMS maxent_c cthyb ctint)
    target_compile_definitions(${_target} PRIVATE BUILD_PYTHON_MODULE)
    target_include_directories(${_target} PRIVATE "${_bindings}" "${_dmft}")
  endforeach()
endif()

# Where the extension modules look for libalps at run time.
#
# When the ALPS libraries ship inside the package (the default), the *only*
# route to libalps must be the in-package pyalps/lib copy -- the same one the
# bundled programs load. Leaving the SDK's own library directory on the list
# as well is not a harmless fallback: auditwheel/delocate resolve the first
# match, find libalps outside the wheel, and vendor a second copy into
# pyalps.libs (or pyalps/.dylibs), so the wheel carries the same 10 MB library
# twice. HDF5 stays on the list precisely because it *should* be vendored.
#
# With PYALPS_BUNDLE_APPLICATIONS=OFF nothing installs libalps into the
# package, so there the SDK directory is the only way for the repair tools to
# find it at all.
if(APPLE)
  set(_pyalps_in_package_libs "@loader_path/../lib")
else()
  set(_pyalps_in_package_libs "$ORIGIN/../lib")
endif()
# The SDK records the selected HDF5 target, including its configuration.
set(_pyalps_hdf5_runtime_dirs "$<TARGET_FILE_DIR:${ALPS_HDF5_TARGETS}>")
if(PYALPS_BUNDLE_APPLICATIONS)
  set(_pyalps_install_rpath "${_pyalps_in_package_libs};${_pyalps_hdf5_runtime_dirs}")
else()
  set(_pyalps_install_rpath "$<TARGET_FILE_DIR:ALPS::alps>;${_pyalps_hdf5_runtime_dirs}")
endif()

foreach(_target IN LISTS _pyalps_targets)
  target_link_libraries(${_target} PRIVATE ALPS::alps)
  set_target_properties(${_target} PROPERTIES INSTALL_RPATH "${_pyalps_install_rpath}")
endforeach()

set(_extension_dir "pyalps/_ext")
install(
  TARGETS ${_pyalps_targets} RUNTIME_DEPENDENCY_SET pyalps_runtime
  LIBRARY DESTINATION "${_extension_dir}"
  RUNTIME DESTINATION "${_extension_dir}")

install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/pyalps/"
  DESTINATION pyalps
  FILES_MATCHING
  PATTERN "*.py")

# Runtime fallback paths for pyalps.tools, used only when the package does
# not carry the resource in question. The XML/XSL library below is always
# bundled, so its fallback is always empty; the application binaries are
# bundled unless PYALPS_BUNDLE_APPLICATIONS is OFF, and only that
# configuration needs to point back at the SDK that built the package.
#
# Never fill these in for a bundled package: they would be absolute paths on
# the build machine, baked into a redistributable artifact.
set(PYALPS_ALPS_XML_FALLBACK "")
if(PYALPS_BUNDLE_APPLICATIONS)
  set(PYALPS_ALPS_BIN_FALLBACK "")
else()
  set(PYALPS_ALPS_BIN_FALLBACK "${ALPS_RUNTIME_LIBRARY_DIRS}")
endif()
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/pyalps/pyalps_config.py.in"
               "${CMAKE_CURRENT_BINARY_DIR}/pyalps_config.py" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pyalps_config.py" DESTINATION pyalps)

# Application locations come from the SDK export, including custom bin dirs.
if(PYALPS_BUNDLE_APPLICATIONS)
  if(NOT ALPS_applications_FOUND)
    message(
      FATAL_ERROR "Bundling applications requires an SDK built with ALPS_BUILD_APPLICATIONS=ON. "
                  "Use PYALPS_BUNDLE_APPLICATIONS=OFF for a bindings-only wheel.")
  endif()
  install(
    IMPORTED_RUNTIME_ARTIFACTS ${ALPS_APPLICATION_TARGETS}
    RUNTIME_DEPENDENCY_SET pyalps_runtime
    RUNTIME DESTINATION pyalps/bin)

  # Imported target metadata identifies the actual runtime and its SONAME on
  # ELF and macOS alike. Copy its bytes once, using the loader-visible name.
  if(NOT WIN32)
    install(
      FILES "$<TARGET_FILE:ALPS::alps>"
      DESTINATION pyalps/lib
      RENAME "$<TARGET_SONAME_FILE_NAME:ALPS::alps>")
  endif()
endif()

# Windows extensions and child processes share the DLLs next to the bundled
# executables. Python registers this directory before loading any extension.
if(WIN32)
  install(
    RUNTIME_DEPENDENCY_SET pyalps_runtime
    DIRECTORIES "${ALPS_RUNTIME_LIBRARY_DIRS}"
    PRE_EXCLUDE_REGEXES "api-ms-.*" "ext-ms-.*" "python[0-9]+.*\\.dll"
    POST_EXCLUDE_REGEXES
    ".*[/\\\\][Ww][Ii][Nn][Dd][Oo][Ww][Ss][/\\\\].*"
    RUNTIME DESTINATION pyalps/bin)
endif()
include("${_alps_source_root}/cmake/ALPSRuntime.cmake")
alps_install_compiler_runtime(pyalps/bin)

# Run after all runtime files have been installed. Wheel repair regenerates this
# manifest once auditwheel/delocate have finalized their paths and install names.
install(
  CODE "execute_process(COMMAND \"${Python_EXECUTABLE}\"
  \"${CMAKE_CURRENT_SOURCE_DIR}/_build_support/runtime_manifest.py\"
  --tree \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}\" --alps-version \"${ALPS_VERSION}\"
  COMMAND_ERROR_IS_FATAL ANY)")

# Reuse the installed SDK's XML and stylesheets.
if(NOT EXISTS "${ALPS_DATA_DIR}/xml/ALPS.xsl")
  message(FATAL_ERROR "The ALPS SDK is missing its XML resources: ${ALPS_DATA_DIR}/xml")
endif()
install(DIRECTORY "${ALPS_DATA_DIR}/xml/" DESTINATION pyalps/xml
  FILES_MATCHING PATTERN "*.xsl" PATTERN "*.xml")
