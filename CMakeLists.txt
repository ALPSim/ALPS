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

cmake_minimum_required(VERSION 3.22)
include(cmake/ALPSVersion.cmake)
# These options must precede project(): vcpkg reads manifest features when the
# toolchain is loaded. Embedded builds default to just the library.
set(_alps_top_level OFF)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(_alps_top_level ON)
endif()
option(BUILD_TESTING "Build and register tests" ${_alps_top_level})
option(ALPS_ENABLE_MPI "Enable MPI parallelization" ${_alps_top_level})
if(BUILD_TESTING)
  list(APPEND VCPKG_MANIFEST_FEATURES tests)
endif()
if(ALPS_ENABLE_MPI)
  list(APPEND VCPKG_MANIFEST_FEATURES mpi)
endif()
project(alps VERSION ${ALPS_VERSION_CORE} LANGUAGES C CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "Use an out-of-source build: cmake -S . -B build")
endif()
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(ALPS_ENABLE_OPENMP "Enable OpenMP parallelization" OFF)
option(ALPS_ENABLE_OPENMP_WORKER "Enable OpenMP worker support" OFF)
option(ALPS_BUILD_FORTRAN "Build Fortran interfaces and applications" OFF)
option(ALPS_BUILD_APPLICATIONS "Build simulation applications and command-line tools" ${PROJECT_IS_TOP_LEVEL})
option(ALPS_BUILD_EXAMPLES "Build examples" OFF)
option(ALPS_INCLUDE_TUTORIALS "Install tutorials" OFF)
option(ALPS_BUILD_DEVELOPER_TOOLS "Build maintenance tools" OFF)
option(BIND_FORTRAN_LOWERCASE "Fortran symbols have no trailing underscore" OFF)
set(ALPS_VERSION_PRERELEASE "" CACHE STRING "Optional prerelease version label")
set(ALPS_XML_PATH "${CMAKE_INSTALL_LIBDIR}" CACHE STRING "XML data parent directory")

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# Build mechanics are private to this directory tree, not SDK usage requirements.
add_compile_options("$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/bigobj;/permissive-;/Zc:__cplusplus>")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build configuration" FORCE)
endif()
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

set(ALPS_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(ALPS_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(ALPS_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(ALPS_VERSION ${PROJECT_VERSION})
if(ALPS_VERSION_PRERELEASE)
  string(APPEND ALPS_VERSION "-${ALPS_VERSION_PRERELEASE}")
endif()
set(ALPS_VERSION_STRING "ALPS Libraries version ${ALPS_VERSION}")
string(TIMESTAMP ALPS_YEAR "%Y" UTC)
set(ALPS_CONFIG_HOST unknown)
set(ALPS_CONFIG_USER unknown)
set(ALPS_PREFIX "${CMAKE_INSTALL_PREFIX}")
set(ALPS_XML_DIR "${CMAKE_INSTALL_PREFIX}/${ALPS_XML_PATH}/xml")

find_package(Threads REQUIRED)
set(ALPS_BOOST_COMPONENTS filesystem serialization program_options regex
  thread date_time chrono timer iostreams)
if(ALPS_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
  list(APPEND ALPS_BOOST_COMPONENTS mpi)
  set(ALPS_HAVE_MPI 1)
  set(ALPS_HAVE_BOOST_MPI 1)
endif()
find_package(Boost 1.76 CONFIG REQUIRED COMPONENTS ${ALPS_BOOST_COMPONENTS})
set(Boost_INCLUDE_DIR ${Boost_INCLUDE_DIRS})

if(ALPS_BUILD_FORTRAN)
  enable_language(Fortran)
  set(ALPS_HDF5_COMPONENTS C HL Fortran)
else()
  set(ALPS_HDF5_COMPONENTS C HL)
endif()
# HDF5's config version file requires a matching major release; check our
# lower bound separately so a newer major can still supply its own targets.
find_package(HDF5 REQUIRED COMPONENTS ${ALPS_HDF5_COMPONENTS})
if(HDF5_VERSION VERSION_LESS 1.8.2)
  message(FATAL_ERROR "ALPS requires HDF5 >= 1.8.2")
endif()
# Prefer HDF5's exported targets: older FindHDF5 aggregate targets flatten
# Debug/Release import libraries into whichever configuration was found first.
set(ALPS_HDF5_TARGETS HDF5::HDF5)
if(TARGET hdf5::hdf5-shared AND NOT HDF5_USE_STATIC_LIBRARIES)
  set(ALPS_HDF5_TARGETS hdf5::hdf5-shared)
  if(ALPS_BUILD_FORTRAN)
    list(APPEND ALPS_HDF5_TARGETS hdf5::hdf5_fortran-shared)
  endif()
elseif(TARGET hdf5::hdf5-static)
  set(ALPS_HDF5_TARGETS hdf5::hdf5-static)
  if(ALPS_BUILD_FORTRAN)
    list(APPEND ALPS_HDF5_TARGETS hdf5::hdf5_fortran-static)
  endif()
endif()
set(ALPS_HAVE_HDF5 1)
set(HDF5_INCLUDE_DIR ${HDF5_INCLUDE_DIRS})
if(HDF5_IS_PARALLEL)
  if(NOT ALPS_ENABLE_MPI)
    message(FATAL_ERROR "Parallel HDF5 requires ALPS_ENABLE_MPI=ON")
  endif()
  set(ALPS_HAVE_HDF5_PARALLEL 1)
endif()
if(ALPS_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()
if(ALPS_ENABLE_OPENMP_WORKER AND NOT ALPS_ENABLE_OPENMP)
  message(FATAL_ERROR "ALPS_ENABLE_OPENMP_WORKER requires ALPS_ENABLE_OPENMP=ON")
endif()

include(ALPSNumerics)
find_package(SQLite QUIET)
find_package(LPSolve QUIET)

include(CheckIncludeFile)
foreach(header IN ITEMS stdarg.h sys/stat.h sys/types.h unistd.h windows.h)
  string(TOUPPER "${header}" feature)
  string(REGEX REPLACE "[/.]" "_" feature "${feature}")
  check_include_file("${header}" "ALPS_HAVE_${feature}")
endforeach()
set(ALPS_HAVE_VALARRAY 1)
# The bundled XDR implementation is portable and preserves the legacy format.
# Using it consistently avoids incompatible platform RPC typedefs and ABIs.
configure_file(src/alps/config.h.in src/alps/config.h)
configure_file(src/alps/version.h.in src/alps/version.h)
configure_file(lib/xml/lattices.xml.in lib/xml/lattices.xml COPYONLY)
configure_file(lib/xml/models.xml.in lib/xml/models.xml COPYONLY)

if(APPLE)
  list(APPEND CMAKE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
  list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include(CTest)
include(add_alps_test)
add_subdirectory(src/alps)
if(ALPS_BUILD_APPLICATIONS)
  add_subdirectory(tool)
  add_subdirectory(applications)
endif()
if(ALPS_BUILD_EXAMPLES)
  add_subdirectory(example)
endif()
if(ALPS_INCLUDE_TUTORIALS)
  add_subdirectory(tutorials)
endif()
if(ALPS_BUILD_DEVELOPER_TOOLS)
  add_subdirectory(script)
endif()
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

install(DIRECTORY src/alps src/boost src/ietl src/mocasito bindings/boost
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT headers
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.ipp")
  install(FILES "${PROJECT_BINARY_DIR}/src/alps/config.h"
    "${PROJECT_BINARY_DIR}/src/alps/version.h"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/alps COMPONENT headers)
install(DIRECTORY lib/xml DESTINATION ${ALPS_XML_PATH} COMPONENT xml
  FILES_MATCHING PATTERN "*.xsl")
install(DIRECTORY "${PROJECT_BINARY_DIR}/lib/xml" DESTINATION ${ALPS_XML_PATH} COMPONENT xml)

configure_package_config_file(cmake/ALPSConfig.cmake.in
  "${PROJECT_BINARY_DIR}/cmake/ALPSConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATADIR}/alps
  PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR CMAKE_INSTALL_BINDIR)
write_basic_package_version_file("${PROJECT_BINARY_DIR}/cmake/ALPSConfigVersion.cmake"
  COMPATIBILITY SameMinorVersion)
install(EXPORT ALPSTargets NAMESPACE ALPS:: DESTINATION ${CMAKE_INSTALL_DATADIR}/alps COMPONENT build)
install(FILES "${PROJECT_BINARY_DIR}/cmake/ALPSConfig.cmake"
  "${PROJECT_BINARY_DIR}/cmake/ALPSConfigVersion.cmake"
  cmake/UseALPS.cmake cmake/UsePyALPS.cmake cmake/add_alps_test.cmake
  cmake/run_test.cmake cmake/run_test_mpi.cmake
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps COMPONENT build)
install(FILES CITATION.md LICENSE.txt README.md
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps COMPONENT libraries)

# Windows resolves DLLs alongside executables. CMake follows the actual PE
# dependency graph, including numerical-library runtimes, at install time.
if(WIN32)
  install(RUNTIME_DEPENDENCY_SET alps_runtime
    DIRECTORIES "$<TARGET_FILE_DIR:alps>"
    PRE_EXCLUDE_REGEXES "api-ms-.*" "ext-ms-.*"
    POST_EXCLUDE_REGEXES ".*[/\\\\][Ww][Ii][Nn][Dd][Oo][Ww][Ss][/\\\\].*"
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT runtime)
endif()

# Use CMake's compiler-runtime discovery for redistributable SDKs and wheels.
include(cmake/ALPSRuntime.cmake)
alps_install_compiler_runtime("${CMAKE_INSTALL_BINDIR}")
