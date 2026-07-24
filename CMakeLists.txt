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

cmake_minimum_required(VERSION 3.18.0)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt.")
endif()

if (NOT_BUILD_SHARED_LIBS)
	option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
else (NOT_BUILD_SHARED_LIBS)
	option(BUILD_SHARED_LIBS "Build shared libraries" ON)
endif (NOT_BUILD_SHARED_LIBS)
option(ALPS_INSTALL_DEPRECATED "Install deprecated headers for ALPS" OFF)
option(ALPS_BUILD_DEVELOPER_TOOLS "Build tools used by developers to maintain ALPS" OFF)
option(ALPS_ENABLE_OPENMP "Enable OpenMP parallelization" OFF)
option(ALPS_ENABLE_OPENMP_WORKER "Enable OpenMP worker support" OFF)
option(ALPS_BUILD_ON_CRAY "Build ALPS for cray" OFF)
option(LAPACK_64_BIT "Use 64-bit version of LAPACK" OFF)
option(ALPS_BUILD_FORTRAN "Build the ALPS fortran applications and tools" OFF)
option(ALPS_INCLUDE_TUTORIALS "Include the ALPS tutorials" ON)
option(BIND_FORTRAN_LOWERCASE "FORTRAN functions are compiled without a trailing underscore" OFF)
option(ALPS_LINK_BOOST_TEST "Link Boost.test pre-built library" OFF)
option(ALPS_INSTALL_BOOST_TEST "Install Boost Test framework" OFF)
option(ALPS_NGS_USE_NEW_ALEA "Use the ALPS ngs-alea instead of just alea" OFF)
option(ALPS_NGS_OPENMPI_ULFM "Userlevel failure mitigation is available" OFF)
option(ALPS_BUILD_LIBS_ONLY "Build only libraries" OFF)



mark_as_advanced(BIND_FORTRAN_LOWERCASE)
mark_as_advanced(ALPS_LINK_BOOST_TEST)
mark_as_advanced(ALPS_INSTALL_BOOST_TEST)
mark_as_advanced(ALPS_NGS_USE_NEW_ALEA)
mark_as_advanced(ALPS_NGS_OPENMPI_ULFM)
mark_as_advanced(ALPS_BUILD_LIBS_ONLY)

option(ALPS_USE_MKL_PARALLEL "Use parallel version of MKL" OFF)
mark_as_advanced(ALPS_USE_MKL_PARALLEL)
    
SET (APPLICATIONS_CAN_BE_BUILT ON)

set(ALPS_BOOST_LIBRARY_NAME "boost" CACHE STRING "name of the boost library")
set(ALPS_ENABLE_MPI ON CACHE BOOL "Enable MPI Parallelization")
set(ALPS_INSTALL_HEADERS  ON CACHE BOOL "Install headers for ALPS and all dependent libraries")
set(ALPS_BUILD_EXAMPLES ON CACHE BOOL "Build ALPS examples")
set(ALPS_BUILD_TESTS ON CACHE BOOL "Build ALPS tests")
set(ALPS_BUILD_APPLICATIONS ${APPLICATIONS_CAN_BE_BUILT} CACHE BOOL "Build ALPS applications")

mark_as_advanced(ALPS_BOOST_LIBRARY_NAME)

if(ALPS_BUILD_LIBS_ONLY)
  set(ALPS_BUILD_APPLICATIONS OFF)
endif()

if(ALPS_BUILD_FORTRAN)
  project(alps C CXX Fortran)
  if (APPLE)
   set (CMAKE_Fortran_RUNTIME_LIBRARIES "-lgcc_s.1")
  endif(APPLE)
else(ALPS_BUILD_FORTRAN)
  project(alps C CXX)
endif(ALPS_BUILD_FORTRAN)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# print extensive compiler information
execute_process(COMMAND ${CMAKE_CXX_COMPILER} --version COMMAND head -n1 OUTPUT_VARIABLE ALPS_COMPILER_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "Compiler version: ${ALPS_COMPILER_VERSION}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -DBOOST_NO_AUTO_PTR -DBOOST_FILESYSTEM_NO_CXX20_ATOMIC_REF -DBOOST_TIMER_ENABLE_DEPRECATED")


# set convenient warning flags for clang
STRING(REGEX MATCH "clang" IS_CLANG "${ALPS_COMPILER_VERSION}")
if(IS_CLANG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth-256 -Wall -Wno-unknown-pragmas -Wno-unused-function -Wno-overloaded-virtual -Wno-unneeded-internal-declaration -Wno-unused-variable -Wno-return-type-c-linkage")
endif(IS_CLANG)

# For gfortran 10 (see https://gcc.gnu.org/gcc-10/changes.html)
if(ALPS_BUILD_FORTRAN AND CMAKE_COMPILER_IS_GNUCC AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 9.99)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
endif(ALPS_BUILD_FORTRAN AND CMAKE_COMPILER_IS_GNUCC AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 9.99)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

######################################################################
# CMAKE_MODULE_PATH
######################################################################
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

######################################################################
# set default CMAKE_INSTALL_PREFIX
######################################################################

######################################################################
# Version information
######################################################################
set(ALPS_YEAR 2026)
set(ALPS_VERSION_MAJOR 2)
set(ALPS_VERSION_MINOR 3)
set(ALPS_VERSION_PATCH 3)
set(ALPS_VERSION_BUILD "")

if(ALPS_VERSION_BUILD)
  set(ALPS_VERSION "${ALPS_VERSION_MAJOR}.${ALPS_VERSION_MINOR}.${ALPS_VERSION_PATCH}-${ALPS_VERSION_BUILD}")
else(ALPS_VERSION_BUILD)
  set(ALPS_VERSION "${ALPS_VERSION_MAJOR}.${ALPS_VERSION_MINOR}.${ALPS_VERSION_PATCH}")
endif(ALPS_VERSION_BUILD)
set(ALPS_VERSION_STRING "ALPS Libraries version ${ALPS_VERSION}")
MESSAGE(STATUS "ALPS version: ${ALPS_VERSION}") 

set(ALPS_CONFIG_HOST unknown)
set(ALPS_CONFIG_USER unknown)

set(ALPS_PREFIX "${CMAKE_INSTALL_PREFIX}")
set(ALPS_SRCDIR "${CMAKE_SOURCE_DIR}")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(bindir "${CMAKE_INSTALL_PREFIX}/bin")

######################################################################
# C & C++ Headers
######################################################################
include(CheckIncludeFile)
CHECK_INCLUDE_FILE(rpc/rpc.h ALPS_HAVE_RPC_XDR_H)
CHECK_INCLUDE_FILE(stdarg.h ALPS_HAVE_STDARG_H)
CHECK_INCLUDE_FILE(sys/stat.h ALPS_HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILE(sys/systeminfo.h ALPS_HAVE_SYS_SYSTEMINFO_H)
CHECK_INCLUDE_FILE(sys/types.h ALPS_HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILE(unistd.h ALPS_HAVE_UNISTD_H)
CHECK_INCLUDE_FILE(windows.h ALPS_HAVE_WINDOWS_H)

include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(valarray ALPS_HAVE_VALARRAY)

######################################################################
# Packages
######################################################################

# MPI
if(ALPS_ENABLE_MPI)
  IF(NOT MPI_FOUND)
    find_package(MPIForALPS)
      message (STATUS "MPI compiler was ${MPI_COMPILER}")
  ENDIF(NOT MPI_FOUND)

  set(MPI_DEFINITIONS ${MPI_COMPILE_FLAGS})
  if(MPI_INCLUDE_DIR)
    set(MPI_DEFINITIONS "${MPI_COMPILE_FLAGS} -I${MPI_INCLUDE_DIR}")
  endif(MPI_INCLUDE_DIR)
 
  set(MPI_INCLUDE_DIR ${MPI_INCLUDE_PATH})
  IF(MPI_FOUND)
    include_directories(${MPI_INCLUDE_DIR})
    SET(ALPS_HAVE_MPI 1)
    SET(ALPS_HAVE_BOOST_MPI 1)
    set(BUILD_BOOST_MPI TRUE)
    foreach(arg ${MPI_LIBRARIES})
      set(MPI_LIBS "${MPI_LIBS} ${arg}")
    endforeach(arg ${MPI_LIBRARIES})
  ELSE(MPI_FOUND)
    set(BUILD_BOOST_MPI FALSE)
  ENDIF(MPI_FOUND)
else (ALPS_ENABLE_MPI)
    SET(ALPS_HAVE_MPI 0)
    SET(ALPS_HAVE_BOOST_MPI 0)
    set(BUILD_BOOST_MPI FALSE)
endif(ALPS_ENABLE_MPI)

# OpenMP
if(ALPS_ENABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    # Almost always OpenMP flags are same both for C and for Fortran.
    if(ALPS_BUILD_FORTRAN)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
    endif(ALPS_BUILD_FORTRAN)
  endif(OPENMP_FOUND)
endif(ALPS_ENABLE_OPENMP)

# Lapack
# For parallel MKL, this section should be executed after OpenMP check.
if (NOT LAPACK_FOUND)
  find_package(Lapack)
endif (NOT LAPACK_FOUND)
if (LAPACK_FOUND)
  set(ALPS_HAVE_LAPACK 1)
  set(ALPS_HAVE_BLAS 1)
endif(LAPACK_FOUND)
IF (REQUIRE_PTHREAD)
  SET(LAPACK_LIBRARY "${LAPACK_LIBRARY};${PTHREAD_LIBRARY}")
  SET(LAPACK_LIBRARIES "${LAPACK_LIBRARIES};${PTHREAD_LIBRARY}")
ENDIF (REQUIRE_PTHREAD)
if(MAC_VECLIB)
  set(LAPACK_LDFLAGS "-framework ${BLAS_LAPACK_MAC_FRAMEWORK}")
endif(MAC_VECLIB)
set(LAPACK_LINKER_FLAGS ${LAPACK_LDFLAGS})

# SQLite
find_package(SQLite)

# LPSolve
find_package(LPSolve)

# SZip
find_package(SZIP)
if(SZIP_FOUND)
  include_directories(${SZIP_INCLUDE_DIR})
  link_libraries(${SZIP_LIBRARIES})
endif(SZIP_FOUND)

# HDF5
if(ALPS_BUILD_FORTRAN)
  find_package(HDF5 COMPONENTS C HL Fortran)
else()
  find_package(HDF5 COMPONENTS C HL)
endif()
IF(HDF5_FOUND)
  SET(HAVE_LIBHDF5 1)
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  add_definitions(${HDF5_DEFINTIONS})
  SET(ALPS_HAVE_HDF5 1)
ELSE(HDF5_FOUND)
  MESSAGE(FATAL_ERROR "Require hdf5 1.8.2 or higher. Set HDF5_HOME")
ENDIF(HDF5_FOUND)


IF(HDF5_IS_PARALLEL)
  SET(ALPS_HAVE_HDF5_PARALLEL 1)
  MESSAGE(WARNING "parallel(MPI) hdf5 is detected. We will compile but ALPS does not use parallel HDF5. The standard version is preferred.")
  IF(NOT MPI_FOUND)
    MESSAGE(FATAL_ERROR "parallel(MPI) hdf5 needs MPI. Enable MPI or install serial HDF5 libraries.")
  ENDIF(NOT MPI_FOUND)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF(HDF5_IS_PARALLEL)

# Python bindings are built by the standalone scikit-build-core project in
# bindings/python/pyalps.  The C++ SDK deliberately has no Python dependency.

# Boost Libraries
find_package(BoostForALPS REQUIRED)
if(Boost_ROOT_DIR)
  message(STATUS "Adding Boost dir: ${Boost_ROOT_DIR} " )
  include_directories(${Boost_ROOT_DIR})
else(Boost_ROOT_DIR)
  message(STATUS "Adding Boost dir: ${Boost_INCLUDE_DIR} " )
  include_directories(${Boost_INCLUDE_DIR})
endif(Boost_ROOT_DIR)

# XML Path
set(ALPS_XML_PATH "lib" CACHE STRING "Directory to install ALPS XML directory in")
set(ALPS_XML_DIR "${CMAKE_INSTALL_PREFIX}/${ALPS_XML_PATH}/xml")

message (STATUS "ALPS XML dir is ${ALPS_XML_DIR}")


if (BIND_FORTRAN_LOWERCASE)
  add_definitions(-DBIND_FORTRAN_LOWERCASE)
endif (BIND_FORTRAN_LOWERCASE)
  
  
######################################################################
# test environment
######################################################################
include(add_alps_test)

######################################################################
# configuration files
######################################################################

set (BIND_FORTRAN_INTEGER_8 ${LAPACK_64_BIT})

configure_file(src/alps/config.h.in ${PROJECT_BINARY_DIR}/src/alps/config.h)
configure_file(src/alps/version.h.in ${PROJECT_BINARY_DIR}/src/alps/version.h)

configure_file(lib/xml/lattices.xml.in ${PROJECT_BINARY_DIR}/lib/xml/lattices.xml)
configure_file(lib/xml/models.xml.in ${PROJECT_BINARY_DIR}/lib/xml/models.xml)

######################################################################
# include directories
######################################################################
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(src)
include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/bindings)

# order the include directories so that directories which are in the
# source or build tree always come before directories outside the project.
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)

######################################################################
# RPATH setting
######################################################################
set(_alps_install_rpath)
if(CMAKE_INSTALL_RPATH)
  list(APPEND _alps_install_rpath ${CMAKE_INSTALL_RPATH})
endif()
if(APPLE)
  list(APPEND _alps_install_rpath "@loader_path/../lib")
else()
  list(APPEND _alps_install_rpath "$ORIGIN/../lib")
endif()
set(CMAKE_INSTALL_RPATH "${_alps_install_rpath}")
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_MACOSX_RPATH 1)

######################################################################
# setting up ALPSConfig.cmake and include.mk
######################################################################

set(ALPS_USE_FILE ${CMAKE_INSTALL_PREFIX}/share/alps/UseALPS.cmake)

set(Boost_INCLUDE_DIR_CONFIG ${Boost_INCLUDE_DIR})

set(ALPS_INCLUDE_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}/include)
list(APPEND ALPS_INCLUDE_DIRS_CONFIG ${Boost_INCLUDE_DIR_CONFIG})
set(ALPS_LIBRARY_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}/lib)
set(ALPS_RUNTIME_LIBRARY_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}/bin)

if(MPI_DEFINITIONS)
  set(ALPS_EXTRA_DEFINITIONS "${ALPS_EXTRA_DEFINITIONS} ${MPI_DEFINITIONS}")
endif(MPI_DEFINITIONS)
if(LAPACK_DEFINITIONS)
  set(ALPS_EXTRA_DEFINITIONS "${ALPS_EXTRA_DEFINITIONS} ${LAPACK_DEFINITIONS}")
endif(LAPACK_DEFINITIONS)
if(HDF5_DEFINITIONS)
  set(ALPS_EXTRA_DEFINITIONS "${ALPS_EXTRA_DEFINITIONS} ${HDF5_DEFINITIONS}")
endif(HDF5_DEFINITIONS)

if(MPI_INCLUDE_DIR)
  list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${MPI_INCLUDE_DIR})
  string(REPLACE ";" " -I" MPI_INCLUDE_DIR_FIXED "${MPI_INCLUDE_DIR}")
  set (MPI_INCLUDE_DIR_FIXED "-I${MPI_INCLUDE_DIR_FIXED}")
else(MPI_INCLUDE_DIR)
  set (MPI_INCLUDE_DIR_FIXED "")
endif(MPI_INCLUDE_DIR)
if(HDF5_INCLUDE_DIR)
  list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
endif(HDF5_INCLUDE_DIR)
if(Boost_INCLUDE_DIR_CONFIG)
  list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${Boost_INCLUDE_DIR_CONFIG})
endif(Boost_INCLUDE_DIR_CONFIG)

if(NOT Boost_FOUND)
  set(Boost_LIBRARIES ${ALPS_BOOST_LIBRARY_NAME})
  set(Boost_LIBS "-l${ALPS_BOOST_LIBRARY_NAME}")
else()
  string (REPLACE ";" " " Boost_LIBS "${Boost_LIBRARIES}") # variable for Makefile
endif(NOT Boost_FOUND)
set(ALPS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY "ALPS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY-NOTFOUND")
if(ALPS_INSTALL_BOOST_TEST)
  set(ALPS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY "-lboost_unit_test_framework -lboost_test_exec_monitor -lboost_prg_exec_monitor")
endif(ALPS_INSTALL_BOOST_TEST)

if(MPI_LINKER_FLAGS)
  set(ALPS_EXTRA_LINKER_FLAGS "${ALPS_EXTRA_LINKER_FLAGS} ${MPI_LINKER_FLAGS}")
endif(MPI_LINKER_FLAGS)
if(LAPACK_LINKER_FLAGS)
  set(ALPS_EXTRA_LINKER_FLAGS "${ALPS_EXTRA_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
endif(LAPACK_LINKER_FLAGS)

if(MPI_LIBRARIES)
  list(APPEND ALPS_EXTRA_LIBRARIES ${MPI_LIBRARIES})
endif(MPI_LIBRARIES)
if(LAPACK_LIBRARIES)
  list(APPEND ALPS_EXTRA_LIBRARIES ${LAPACK_LIBRARIES})
endif(LAPACK_LIBRARIES)
if(HDF5_LIBRARIES)
  list(APPEND ALPS_EXTRA_LIBRARIES ${HDF5_LIBRARIES})
endif(HDF5_LIBRARIES)
configure_file(cmake/ALPSConfig.cmake.in ${PROJECT_BINARY_DIR}/cmake/ALPSConfig.cmake @ONLY)
configure_file(cmake/include.mk.in ${PROJECT_BINARY_DIR}/cmake/include.mk)


# installation
######################################################################
include(InstallRequiredSystemLibraries)
if (ALPS_INSTALL_HEADERS)
  set(ALPS_HEADER_DIR "include")
  install(DIRECTORY src/alps src/boost src/ietl src/mocasito COMPONENT headers DESTINATION ${ALPS_HEADER_DIR}
          FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.ipp"
         )
  # ALPS' public headers include the vendored Boost.Numeric bindings.  Install
  # them as part of the C++ SDK so downstream consumers do not need the ALPS
  # source tree on their include path.
  install(DIRECTORY bindings/boost COMPONENT headers DESTINATION ${ALPS_HEADER_DIR}
          FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.ipp"
         )
  install(FILES ${PROJECT_BINARY_DIR}/src/alps/config.h ${PROJECT_BINARY_DIR}/src/alps/version.h DESTINATION ${ALPS_HEADER_DIR}/alps COMPONENT headers)
  if (ALPS_INSTALL_DEPRECATED)
    install(DIRECTORY legacy/alps COMPONENT deprecated DESTINATION ${ALPS_HEADER_DIR}
            FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
           )
  endif (ALPS_INSTALL_DEPRECATED)
endif(ALPS_INSTALL_HEADERS)

install(DIRECTORY lib/xml DESTINATION ${ALPS_XML_PATH} COMPONENT xml
        FILES_MATCHING PATTERN "*.xsl"
       )

install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/xml DESTINATION ${ALPS_XML_PATH} COMPONENT xml)

add_subdirectory(cmake)

install(FILES cmake/UseALPS.cmake
        ${PROJECT_BINARY_DIR}/cmake/ALPSConfig.cmake
        ${PROJECT_BINARY_DIR}/cmake/include.mk
        cmake/run_test.cmake
        cmake/run_test_mpi.cmake
        cmake/add_alps_test.cmake
        DESTINATION share/alps COMPONENT build)

install(FILES CITATION.md LICENSE.txt README.md
        DESTINATION share/alps COMPONENT libraries)
######################################################################
# libraries
######################################################################

add_subdirectory(src/boost)
add_subdirectory(src/alps)


######################################################################
# programs
######################################################################

if (NOT ALPS_BUILD_LIBS_ONLY)
  add_subdirectory(tool)
  if (ALPS_BUILD_EXAMPLES)
    add_subdirectory(example)
  endif (ALPS_BUILD_EXAMPLES)
endif (NOT ALPS_BUILD_LIBS_ONLY)


######################################################################
# applications
######################################################################

if (ALPS_BUILD_APPLICATIONS AND NOT ALPS_BUILD_LIBS_ONLY)
  add_subdirectory(applications)
endif (ALPS_BUILD_APPLICATIONS AND NOT ALPS_BUILD_LIBS_ONLY)

if (ALPS_INCLUDE_TUTORIALS AND NOT ALPS_BUILD_LIBS_ONLY)
  add_subdirectory(tutorials)
endif (ALPS_INCLUDE_TUTORIALS AND NOT ALPS_BUILD_LIBS_ONLY)

######################################################################
# developer tools
######################################################################

if (ALPS_BUILD_DEVELOPER_TOOLS)
  add_subdirectory(script)
endif (ALPS_BUILD_DEVELOPER_TOOLS)

######################################################################
# tests
######################################################################

if (ALPS_BUILD_TESTS)
  include(ALPSBuildname)
  enable_testing()
  include(CTest)
  add_subdirectory(test)
endif (ALPS_BUILD_TESTS)
