# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.27...4.3)
# One dependency scan supplies the package's Windows DLL directory.
set(VCPKG_APPLOCAL_DEPS OFF)
project(pyalps LANGUAGES C CXX)

option(PYALPS_BUILD_SOLVERS "Build the MaxEnt, CT-HYB, and CT-INT bindings" ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Python 3.11 REQUIRED COMPONENTS Interpreter Development.Module)
if(NOT Python_INTERPRETER_ID STREQUAL "Python" OR Python_SOABI MATCHES "[0-9]t")
  message(FATAL_ERROR "pyalps requires GIL-enabled CPython")
endif()
set(_pyalps_nanobind_runtime nanobind)
if(Python_VERSION VERSION_GREATER_EQUAL 3.12)
  find_package(Python REQUIRED COMPONENTS Development.SABIModule)
  set(_pyalps_nanobind_runtime nanobind-abi3)
endif()
execute_process(
  COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/_build_support/alps_version.py" --core
  OUTPUT_VARIABLE PYALPS_VERSION_CORE
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
# Entry points are resolved before CMake runs. Use the same environment setting
# for metadata and native payload, rather than a CMake-only option that could
# leave bindings-only wheels with launchers that overwrite SDK executables.
execute_process(
  COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/_build_support/alps_version.py"
    --bundle-applications
  OUTPUT_VARIABLE _pyalps_bundle_applications
  OUTPUT_STRIP_TRAILING_WHITESPACE
  COMMAND_ERROR_IS_FATAL ANY)
if(DEFINED PYALPS_BUNDLE_APPLICATIONS AND
   ((PYALPS_BUNDLE_APPLICATIONS AND NOT _pyalps_bundle_applications) OR
    (NOT PYALPS_BUNDLE_APPLICATIONS AND _pyalps_bundle_applications)))
  message(FATAL_ERROR
    "Set the environment variable PYALPS_BUNDLE_APPLICATIONS=OFF (or ON) for "
    "the entire Python build, instead of -DPYALPS_BUNDLE_APPLICATIONS. "
    "Wheel entry points and the CMake payload must use the same setting. "
    "Remove PYALPS_BUNDLE_APPLICATIONS from an existing CMake cache when switching modes.")
endif()
unset(PYALPS_BUNDLE_APPLICATIONS CACHE)
set(PYALPS_BUNDLE_APPLICATIONS "${_pyalps_bundle_applications}")

set(_alps_components utilities hdf5 params)
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
find_package(nanobind 2.15.0 EXACT CONFIG REQUIRED)
file(STRINGS "${NB_DIR}/src/nb_abi.h" _nb_abi_line
  REGEX "^#[ \t]*define[ \t]+NB_INTERNALS_VERSION[ \t]+[0-9]+")
string(REGEX MATCH "[0-9]+$" PYALPS_NANOBIND_ABI "${_nb_abi_line}")
if(NOT PYALPS_NANOBIND_ABI)
  message(FATAL_ERROR "Cannot determine the nanobind internals ABI")
endif()
if(APPLE)
  set(_loader_origin "@loader_path")
elseif(UNIX)
  set(_loader_origin "$ORIGIN")
endif()

# nanobind uses abi3 from Python 3.12 onward and the CPython ABI on 3.11.
function(pyalps_add_module target source)
  nanobind_add_module(${target} STABLE_ABI NB_SHARED MUSL_DYNAMIC_LIBCPP ${source})
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
set_target_properties(${_pyalps_nanobind_runtime} PROPERTIES OUTPUT_NAME pyalps_nanobind
  INSTALL_RPATH "${_loader_origin}" INSTALL_RPATH_USE_LINK_PATH TRUE)
install(TARGETS ${_pyalps_nanobind_runtime} RUNTIME_DEPENDENCY_SET pyalps_runtime
  LIBRARY DESTINATION pyalps/lib RUNTIME DESTINATION pyalps/bin)

# Required libraries are bundled whether or not command-line programs are wanted.
if(WIN32)
  # An unused import library need not appear in an extension's PE imports.
  # Install every component and include its dependencies in the DLL scan.
  install(IMPORTED_RUNTIME_ARTIFACTS ${ALPS_RUNTIME_TARGETS}
    RUNTIME_DEPENDENCY_SET pyalps_runtime RUNTIME DESTINATION pyalps/bin)
else()
  foreach(runtime IN LISTS ALPS_RUNTIME_TARGETS)
    install(FILES "$<TARGET_FILE:${runtime}>" DESTINATION pyalps/lib
      RENAME "$<TARGET_SONAME_FILE_NAME:${runtime}>")
  endforeach()
endif()
if(PYALPS_BUNDLE_APPLICATIONS)
  install(IMPORTED_RUNTIME_ARTIFACTS ${ALPS_APPLICATION_TARGETS}
    RUNTIME_DEPENDENCY_SET pyalps_runtime RUNTIME DESTINATION pyalps/bin)
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    find_program(PYALPS_PATCHELF patchelf REQUIRED)
    # Imported programs retain their SDK RPATH (which may name lib64 or a
    # custom libdir). Resolve our bundled runtime before wheel repair scans them.
    foreach(application IN LISTS ALPS_APPLICATION_TARGETS)
      install(CODE "execute_process(COMMAND \"${PYALPS_PATCHELF}\"
        --add-rpath \"\$ORIGIN/../lib\"
        \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/pyalps/bin/$<TARGET_FILE_NAME:${application}>\"
        COMMAND_ERROR_IS_FATAL ANY)")
    endforeach()
  endif()
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

# Bundle the catalog belonging to the linked SDK.
install(FILES "${ALPS_DATA_DIR}/CITATION.cff" "${ALPS_DATA_DIR}/CITATIONS.yaml"
              "${ALPS_DATA_DIR}/CITATION.md" DESTINATION pyalps/share/alps)

# scikit-build owns Python sources. CMake owns native build and runtime resources.
install(FILES _build_support/pyalpsConfig.cmake DESTINATION pyalps/cmake)
install(FILES cpp/export_simulation.hpp
  DESTINATION pyalps/include/pyalps)
if(NOT EXISTS "${ALPS_DATA_DIR}/xml/ALPS.xsl")
  message(FATAL_ERROR "The ALPS SDK is missing its XML resources: ${ALPS_DATA_DIR}/xml")
endif()
install(DIRECTORY "${ALPS_DATA_DIR}/xml/" DESTINATION pyalps/xml
  FILES_MATCHING PATTERN "*.xsl" PATTERN "*.xml")
install(FILES "${ALPS_DATA_DIR}/licenses/xdr/LICENSE" DESTINATION pyalps/licenses/xdr)

# Only bindings-only builds may fall back to their chosen SDK executables.
set(_pyalps_sdk_bin "")
if(NOT PYALPS_BUNDLE_APPLICATIONS)
  set(_pyalps_sdk_bin "${ALPS_BIN_DIR}")
endif()

# Repair tools regenerate this manifest after changing library names and paths.
install(CODE "execute_process(COMMAND \"${Python_EXECUTABLE}\"
  \"${CMAKE_CURRENT_SOURCE_DIR}/_build_support/runtime_manifest.py\"
  --tree \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}\" --alps-version \"${ALPS_VERSION}\"
  --sdk-bin \"${_pyalps_sdk_bin}\"
  --nanobind-version \"${nanobind_VERSION}\" --nanobind-abi \"${PYALPS_NANOBIND_ABI}\"
  COMMAND_ERROR_IS_FATAL ANY)")
