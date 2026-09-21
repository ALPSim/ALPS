#  Copyright Synge Todo and Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 4.3)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
  message(FATAL_ERROR "Use an out-of-source build: cmake -S . -B build")
endif()
include(cmake/ALPSVersion.cmake)
# These options must precede project(): vcpkg reads manifest features when the
# toolchain is loaded. Embedded builds default to just the library.
set(_alps_top_level OFF)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(_alps_top_level ON)
endif()
option(BUILD_TESTING "Build and register tests" ${_alps_top_level})
option(ALPS_ENABLE_MPI "Enable MPI parallelization" OFF)
option(ALPS_BUILD_APPLICATIONS "Build simulation applications and command-line tools"
       ${_alps_top_level})
if(BUILD_TESTING)
  list(APPEND VCPKG_MANIFEST_FEATURES tests)
endif()
if(ALPS_ENABLE_MPI)
  list(APPEND VCPKG_MANIFEST_FEATURES mpi)
endif()
project(
  alps
  VERSION ${ALPS_VERSION_CORE}
  LANGUAGES C CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(ALPS_ENABLE_OPENMP "Enable OpenMP parallelization" OFF)
option(ALPS_BUILD_EXAMPLES "Build examples" OFF)
option(ALPS_BUILD_EXTENSIVE_TESTS "Build expensive graph and HDF5 type-matrix tests" OFF)
set(ALPS_VERSION_PRERELEASE
    ""
    CACHE STRING "Optional prerelease version label")

# Language requirements belong to ALPS::headers and propagate to every consumer.
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT DEFINED CMAKE_POSITION_INDEPENDENT_CODE)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()
# Build mechanics are private to this directory tree, not SDK usage requirements.
add_compile_options("$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/bigobj;/permissive-;/Zc:__cplusplus>")
if(NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
endif()
foreach(kind IN ITEMS LIBRARY ARCHIVE)
  if(NOT DEFINED CMAKE_${kind}_OUTPUT_DIRECTORY)
    set(CMAKE_${kind}_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
  endif()
endforeach()

set(ALPS_VERSION ${PROJECT_VERSION})
if(ALPS_VERSION_PRERELEASE)
  string(APPEND ALPS_VERSION "-${ALPS_VERSION_PRERELEASE}")
endif()
string(TIMESTAMP ALPS_YEAR "%Y" UTC)

find_package(Threads REQUIRED)
set(ALPS_BOOST_COMPONENTS
    filesystem
    serialization
    program_options
    regex
    thread
    date_time
    chrono
    timer
    iostreams)
if(ALPS_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
  list(APPEND ALPS_BOOST_COMPONENTS mpi)
  set(ALPS_HAVE_MPI 1)
  set(ALPS_HAVE_BOOST_MPI 1)
endif()
find_package(Boost 1.76 CONFIG REQUIRED COMPONENTS ${ALPS_BOOST_COMPONENTS})

find_package(HDF5 REQUIRED COMPONENTS C)
# Prefer HDF5's exported targets: older FindHDF5 aggregate targets flatten
# Debug/Release import libraries into whichever configuration was found first.
set(ALPS_HDF5_TARGETS hdf5::hdf5)
if(TARGET hdf5::hdf5-shared AND NOT HDF5_USE_STATIC_LIBRARIES)
  set(ALPS_HDF5_TARGETS hdf5::hdf5-shared)
elseif(TARGET hdf5::hdf5-static)
  set(ALPS_HDF5_TARGETS hdf5::hdf5-static)
endif()
if(HDF5_IS_PARALLEL)
  if(NOT ALPS_ENABLE_MPI)
    message(FATAL_ERROR "Parallel HDF5 requires ALPS_ENABLE_MPI=ON")
  endif()
  set(ALPS_HAVE_HDF5_PARALLEL 1)
endif()
if(ALPS_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()

include(cmake/ALPSNumerics.cmake)

# The bundled XDR implementation is portable and preserves the legacy format.
# Using it consistently avoids incompatible platform RPC typedefs and ABIs.
configure_file(src/alps/config.h.in src/alps/config.h)
configure_file(src/alps/version.h.in src/alps/version.h)

file(RELATIVE_PATH _alps_bin_to_lib "${CMAKE_INSTALL_FULL_BINDIR}" "${CMAKE_INSTALL_FULL_LIBDIR}")
if(APPLE)
  set(ALPS_LOADER_ORIGIN "@loader_path")
elseif(UNIX)
  set(ALPS_LOADER_ORIGIN "$ORIGIN")
endif()
if(ALPS_LOADER_ORIGIN)
  list(APPEND CMAKE_INSTALL_RPATH "${ALPS_LOADER_ORIGIN}/${_alps_bin_to_lib}")
endif()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include(CTest)
include(cmake/ALPSTesting.cmake)
add_subdirectory(src/alps)
if(ALPS_BUILD_APPLICATIONS)
  add_subdirectory(src/tools)
  add_subdirectory(src/apps)
  get_target_property(ALPS_APPLICATION_TARGETS alps ALPS_APPLICATION_TARGETS)
endif()
add_subdirectory(tutorials)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()

configure_package_config_file(
  cmake/ALPSConfig.cmake.in "${PROJECT_BINARY_DIR}/cmake/ALPSConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATADIR}/alps
  NO_SET_AND_CHECK_MACRO
  PATH_VARS CMAKE_INSTALL_BINDIR CMAKE_INSTALL_DATADIR)
write_basic_package_version_file("${PROJECT_BINARY_DIR}/cmake/ALPSConfigVersion.cmake"
                                 COMPATIBILITY SameMinorVersion)
install(
  EXPORT ALPSTargets
  NAMESPACE ALPS::
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps
  COMPONENT build)
install(
  FILES "${PROJECT_BINARY_DIR}/cmake/ALPSConfig.cmake"
        "${PROJECT_BINARY_DIR}/cmake/ALPSConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps
  COMPONENT build)
install(
  FILES CHANGELOG.md CITATION.md CONTRIBUTING.md LICENSE.txt README.md
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps
  COMPONENT libraries)

# Windows resolves DLLs alongside executables. CMake follows the actual PE
# dependency graph, including numerical-library runtimes, at install time.
if(WIN32)
  install(
    RUNTIME_DEPENDENCY_SET alps_runtime
    DIRECTORIES "$<TARGET_FILE_DIR:alps>"
    PRE_EXCLUDE_REGEXES "api-ms-.*" "ext-ms-.*"
    POST_EXCLUDE_REGEXES
    ".*[/\\\\][Ww][Ii][Nn][Dd][Oo][Ww][Ss][/\\\\].*"
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT runtime)
endif()

# Use CMake's compiler-runtime discovery for redistributable installations.
set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP TRUE)
include(InstallRequiredSystemLibraries)
# The MSVC ARM64 redist also contains an ARM64EC-only exception runtime.
if(CMAKE_CXX_COMPILER_ARCHITECTURE_ID STREQUAL "ARM64")
  list(FILTER CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS EXCLUDE REGEX "/vcruntime140_1d?\\.dll$")
endif()
install(PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
  DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT runtime)
install(DIRECTORY ${CMAKE_INSTALL_SYSTEM_RUNTIME_DIRECTORIES}
  DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT runtime)
