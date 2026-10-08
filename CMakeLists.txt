#  Copyright Synge Todo and Matthias Troyer 2009 - 2010.
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the “Software”),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#  
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#  
#   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

cmake_minimum_required(VERSION 3.27...4.3)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
  message(FATAL_ERROR "Use an out-of-source build: cmake -S . -B build")
endif()
include(cmake/ALPSVersion.cmake)
# Embedded builds default to just the library.
string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}" _alps_top_level)
option(ALPS_BUILD_TESTING "Build and register ALPS tests" ${_alps_top_level})
option(ALPS_ENABLE_MPI "Enable MPI parallelization" OFF)
option(ALPS_BUILD_APPLICATIONS "Build simulation applications and command-line tools"
       ${_alps_top_level})
project(alps VERSION ${ALPS_VERSION_CORE} LANGUAGES C CXX)

# Removed legacy options are otherwise ignored with only CMake's generic
# "Manually-specified variables were not used" warning. See CHANGELOG.md.
foreach(legacy IN ITEMS
    "ALPS_NGS_USE_NEW_ALEA|the experimental backend was removed"
    "ALPS_BUILD_TESTS|use ALPS_BUILD_TESTING"
    "ALPS_BUILD_LIBS_ONLY|use ALPS_BUILD_APPLICATIONS=OFF"
    "ALPS_BUILD_EXAMPLES|build tutorials/00-examples against the installed SDK"
    "ALPS_BUILD_FORTRAN|ALPS::fortran is always built"
    "ALPS_INSTALL_HEADERS|headers are always installed"
    "ALPS_INCLUDE_TUTORIALS|use cmake --install <build> --component tutorials"
    "ALPS_ENABLE_OPENMP_WORKER|use ALPS_ENABLE_OPENMP"
    "ALPS_USE_SYSTEM_BOOST|an installed Boost >= 1.76 is always required"
    "Boost_SRC_DIR|bundled Boost builds were removed; set Boost_ROOT"
    "LAPACK_64_BIT|only LP64 BLAS/LAPACK is supported")
  string(REPLACE "|" ";" legacy "${legacy}")
  list(GET legacy 0 name)
  list(GET legacy 1 hint)
  if(DEFINED ${name})
    message(WARNING "${name} is no longer supported: ${hint}.")
  endif()
endforeach()
# ALPS 3.0 always cached ALPS_NGS_USE_NEW_ALEA and defaulted MPI to ON; a
# reused build directory keeps that value although the default is now OFF.
if(ALPS_ENABLE_MPI AND DEFINED CACHE{ALPS_NGS_USE_NEW_ALEA})
  message(WARNING "This build directory was configured by ALPS 3.0, whose "
          "ALPS_ENABLE_MPI default was ON; the cached value is still in effect. "
          "Pass -DALPS_ENABLE_MPI=OFF or use a fresh build directory if you do not need MPI.")
endif()
if(NOT ALPS_ENABLE_MPI)
  message(STATUS "MPI is disabled; pass -DALPS_ENABLE_MPI=ON to enable it")
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(cmake/ALPSCitations.cmake)

# A local default must not change later add_library() calls in a parent project.
if(NOT DEFINED BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS ON)
endif()
option(ALPS_NGS_OPENMPI_ULFM "Enable user-level failure mitigation with a compatible OpenMPI" OFF)
mark_as_advanced(ALPS_NGS_OPENMPI_ULFM)
option(ALPS_ENABLE_OPENMP "Enable OpenMP parallelization" OFF)
option(ALPS_BUILD_EXTENSIVE_TESTS "Build expensive graph and HDF5 type-matrix tests" OFF)
option(ALPS_ENABLE_SANITIZERS "Instrument ALPS and tests with address/undefined sanitizers" OFF)
if(ALPS_ENABLE_SANITIZERS)
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" OR MSVC)
    message(FATAL_ERROR "ALPS_ENABLE_SANITIZERS requires GCC or Clang with a Unix-style driver")
  endif()
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address,undefined)
endif()
set(ALPS_VERSION_PRERELEASE "" CACHE STRING "Optional prerelease version label")

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
set(ALPS_BOOST_COMPONENTS filesystem serialization program_options regex thread date_time chrono)
if(ALPS_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
  list(APPEND ALPS_BOOST_COMPONENTS mpi)
  set(ALPS_HAVE_MPI 1)
  set(ALPS_HAVE_BOOST_MPI 1)
endif()
find_package(Boost 1.76 CONFIG REQUIRED COMPONENTS ${ALPS_BOOST_COMPONENTS})

# FindHDF5 supplies HDF5::HDF5 and preserves the provider's configurations.
# Keep the manylinux_2_28 system package within the supported source-build range.
set(ALPS_HDF5_MIN_VERSION 1.10.5)
find_package(HDF5 ${ALPS_HDF5_MIN_VERSION} MODULE REQUIRED COMPONENTS C)
if(HDF5_IS_PARALLEL AND NOT ALPS_ENABLE_MPI)
  message(FATAL_ERROR "Parallel HDF5 requires ALPS_ENABLE_MPI=ON")
endif()
set(ALPS_HAVE_HDF5_PARALLEL "${HDF5_IS_PARALLEL}")
if(ALPS_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()

include(cmake/ALPSNumerics.cmake)

# The bundled XDR implementation is portable and preserves the legacy format.
# Using it consistently avoids incompatible platform RPC typedefs and ABIs.
configure_file(cmake/config/config.h.in generated/include/alps/config.h)
configure_file(cmake/config/version.h.in generated/include/alps/version.h)

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

if(ALPS_BUILD_TESTING)
  enable_testing()
endif()
include(cmake/ALPSTesting.cmake)
add_subdirectory(src/alps)
if(ALPS_BUILD_APPLICATIONS)
  if(ALPS_LOADER_ORIGIN)
    # Programs also link directly to the SDK's dependencies. ELF RUNPATH is
    # not transitive, and automatic RPATH discovery omits dependency archives
    # staged inside the source/build tree. Reuse the configuration-aware paths
    # already derived for libalps, preserving the programs' relative SDK path.
    get_target_property(_alps_dependency_rpaths alps INSTALL_RPATH)
    list(REMOVE_ITEM _alps_dependency_rpaths "${ALPS_LOADER_ORIGIN}")
    list(APPEND CMAKE_INSTALL_RPATH ${_alps_dependency_rpaths})
    list(REMOVE_DUPLICATES CMAKE_INSTALL_RPATH)
  endif()
  add_subdirectory(src/tools)
  add_subdirectory(src/apps)
  get_target_property(ALPS_APPLICATION_TARGETS alps ALPS_APPLICATION_TARGETS)
endif()
add_subdirectory(tutorials)
if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

get_property(ALPS_SDK_COMPONENTS GLOBAL PROPERTY ALPS_SDK_COMPONENTS)
get_property(ALPS_SDK_RUNTIME_TARGETS GLOBAL PROPERTY ALPS_SDK_RUNTIME_TARGETS)
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
  FILES CHANGELOG.md CONTRIBUTING.md LICENSE.txt README.md
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
