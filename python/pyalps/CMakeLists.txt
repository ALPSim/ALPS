# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 4.3)
# One dependency scan supplies the package's Windows DLL directory.
set(VCPKG_APPLOCAL_DEPS OFF)
project(pyalps LANGUAGES C CXX)

option(PYALPS_BUILD_SOLVERS "Build the MaxEnt, CT-HYB, and CT-INT bindings" ON)
option(PYALPS_BUNDLE_APPLICATIONS "Include the ALPS command-line programs" ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development.Module Development.SABIModule)
execute_process(
  COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/_build_support/alps_version.py" --core
  OUTPUT_VARIABLE PYALPS_VERSION_CORE
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
set(_alps_components "")
if(PYALPS_BUILD_SOLVERS)
  list(APPEND _alps_components solvers)
endif()
if(PYALPS_BUNDLE_APPLICATIONS)
  list(APPEND _alps_components applications)
endif()
find_package(ALPS ${PYALPS_VERSION_CORE} EXACT REQUIRED CONFIG COMPONENTS ${_alps_components})
if(NOT ALPS_VERSION VERSION_EQUAL PYALPS_VERSION_CORE)
  message(FATAL_ERROR "pyalps ${PYALPS_VERSION_CORE} requires a matching ALPS SDK; found ${ALPS_VERSION}")
endif()
get_target_property(_alps_library_type ALPS::alps TYPE)
if(NOT _alps_library_type STREQUAL "SHARED_LIBRARY")
  message(FATAL_ERROR "pyalps requires a shared ALPS SDK so extensions share one ALPS runtime")
endif()

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_VARIABLE nanobind_ROOT
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
find_package(nanobind 2.10...<3 CONFIG REQUIRED)
if(NOT Python_INTERPRETER_ID STREQUAL "Python" OR NB_FREE_THREADED)
  message(FATAL_ERROR "pyalps requires GIL-enabled CPython for its cp312-abi3 wheels")
endif()

if(APPLE)
  set(_loader_origin "@loader_path")
elseif(UNIX)
  set(_loader_origin "$ORIGIN")
endif()

# Every extension and the shared nanobind runtime use the CPython 3.12 stable ABI.
function(pyalps_add_module target source)
  nanobind_add_module(${target} NB_SHARED STABLE_ABI MUSL_DYNAMIC_LIBCPP ${source})
  target_link_libraries(${target} PRIVATE ALPS::alps ${ARGN})
  # CMake derives dependency locations from the resolved link interface.
  # Wheel repair removes external paths; the packaged runtime always comes first.
  set_target_properties(${target} PROPERTIES
    INSTALL_RPATH "${_loader_origin}/../lib"
    INSTALL_RPATH_USE_LINK_PATH TRUE)
  install(TARGETS ${target} RUNTIME_DEPENDENCY_SET pyalps_runtime
    LIBRARY DESTINATION pyalps/_ext
    RUNTIME DESTINATION pyalps/_ext)
endfunction()

pyalps_add_module(pyalea_c cpp/pyalea.cpp)
pyalps_add_module(pymcdata_c cpp/pymcdata.cpp)
pyalps_add_module(pytools_c cpp/pytools.cpp)
pyalps_add_module(pyngsparams_c cpp/params.cpp)
pyalps_add_module(pyngshdf5_c cpp/hdf5.cpp)
pyalps_add_module(pyngsbase_c cpp/mcbase.cpp)
pyalps_add_module(pyngsobservable_c cpp/observable.cpp)
pyalps_add_module(pyngsobservables_c cpp/observables.cpp)
pyalps_add_module(pyngsresult_c cpp/result.cpp)
pyalps_add_module(pyngsresults_c cpp/results.cpp)
pyalps_add_module(pyngsapi_c cpp/api.cpp)
pyalps_add_module(pyngsrandom01_c cpp/random01.cpp)
pyalps_add_module(pyngsaccumulator_c cpp/accumulator.cpp)
if(PYALPS_BUILD_SOLVERS)
  pyalps_add_module(maxent_c cpp/maxent_c.cpp ALPS::maxent)
  pyalps_add_module(cthyb cpp/cthyb.cpp ALPS::cthyb)
  pyalps_add_module(ctint cpp/ctint.cpp ALPS::ctint)
endif()

# A package-specific name avoids DLL/SONAME collisions with other nanobind users.
set_target_properties(nanobind-abi3 PROPERTIES OUTPUT_NAME pyalps_nanobind
  INSTALL_RPATH "${_loader_origin}" INSTALL_RPATH_USE_LINK_PATH TRUE)
install(TARGETS nanobind-abi3 RUNTIME_DEPENDENCY_SET pyalps_runtime
  LIBRARY DESTINATION pyalps/lib RUNTIME DESTINATION pyalps/bin)

# Required libraries are bundled whether or not command-line programs are wanted.
if(NOT WIN32)
  install(FILES "$<TARGET_FILE:ALPS::alps>" DESTINATION pyalps/lib
    RENAME "$<TARGET_SONAME_FILE_NAME:ALPS::alps>")
endif()
if(PYALPS_BUNDLE_APPLICATIONS)
  install(IMPORTED_RUNTIME_ARTIFACTS ${ALPS_APPLICATION_TARGETS}
    RUNTIME_DEPENDENCY_SET pyalps_runtime RUNTIME DESTINATION pyalps/bin)
endif()
if(WIN32)
  install(RUNTIME_DEPENDENCY_SET pyalps_runtime
    DIRECTORIES "${ALPS_RUNTIME_LIBRARY_DIRS}"
    PRE_EXCLUDE_REGEXES "api-ms-.*" "ext-ms-.*" "python[0-9]+.*\\.dll"
    POST_EXCLUDE_REGEXES ".*[/\\\\][Ww][Ii][Nn][Dd][Oo][Ww][Ss][/\\\\].*"
    RUNTIME DESTINATION pyalps/bin)
  set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP TRUE)
  include(InstallRequiredSystemLibraries)
  # The ARM64 redistributable directory also contains an ARM64EC-only DLL.
  if(CMAKE_CXX_COMPILER_ARCHITECTURE_ID STREQUAL "ARM64")
    list(FILTER CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS EXCLUDE REGEX "/vcruntime140_1d?\\.dll$")
  endif()
  install(PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS} DESTINATION pyalps/bin)
endif()

# scikit-build owns Python sources. CMake owns native build and runtime resources.
install(FILES _build_support/pyalpsConfig.cmake DESTINATION pyalps/cmake)
install(FILES cpp/export_simulation.hpp DESTINATION pyalps/include/pyalps)
if(NOT EXISTS "${ALPS_DATA_DIR}/xml/ALPS.xsl")
  message(FATAL_ERROR "The ALPS SDK is missing its XML resources: ${ALPS_DATA_DIR}/xml")
endif()
install(DIRECTORY "${ALPS_DATA_DIR}/xml/" DESTINATION pyalps/xml
  FILES_MATCHING PATTERN "*.xsl" PATTERN "*.xml")
install(FILES "${ALPS_DATA_DIR}/licenses/xdr/LICENSE" DESTINATION pyalps/licenses/xdr)

# Repair tools regenerate this manifest after changing library names and paths.
install(CODE "execute_process(COMMAND \"${Python_EXECUTABLE}\"
  \"${CMAKE_CURRENT_SOURCE_DIR}/_build_support/runtime_manifest.py\"
  --tree \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}\" --alps-version \"${ALPS_VERSION}\"
  COMMAND_ERROR_IS_FATAL ANY)")
