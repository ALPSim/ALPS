#  Copyright Synge Todo and Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# If you want to load a machine configuration,
# give the filename with full-path or the name of your machine:
#  $ cmake -DALPS_USE_MACHINE=/full/path/to/machine.cmake
# or
#  $ cmake -DALPS_USE_MACHINE=machine-name

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt.")
endif()

if (ALPS_USE_MACHINE)
  if (EXISTS ${ALPS_USE_MACHINE})
    # given full-path case
    set(ALPS_MACHINE_CONF ${ALPS_USE_MACHINE})
    message("-- Loading machine configuration: ${ALPS_MACHINE_CONF}")
    include(${ALPS_MACHINE_CONF})
  else (EXISTS ${ALPS_USE_MACHINE})
    # given machine-name case
    set(ALPS_MACHINE_CONF "${CMAKE_CURRENT_SOURCE_DIR}/config/${ALPS_USE_MACHINE}.cmake")
    if (EXISTS ${ALPS_MACHINE_CONF})
      message("-- Loading machine configuration: ${ALPS_MACHINE_CONF}")
      include(${ALPS_MACHINE_CONF})
    else (EXISTS ${ALPS_MACHINE_CONF})
      message(FATAL_ERROR "Can't find your macnine-configuration file")
    endif (EXISTS ${ALPS_MACHINE_CONF})
  endif (EXISTS ${ALPS_USE_MACHINE})
endif (ALPS_USE_MACHINE)

if (NOT_BUILD_SHARED_LIBS)
	option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
else (NOT_BUILD_SHARED_LIBS)
	option(BUILD_SHARED_LIBS "Build shared libraries" ON)
endif (NOT_BUILD_SHARED_LIBS)
option(ALPS_FOR_VISTRAILS "Build for use within VisTrails" OFF)
if (NOT_ALPS_BUILD_PYTHON)
	option(ALPS_BUILD_PYTHON "Build ALPS python extentions" OFF)
else (NOT_ALPS_BUILD_PYTHON)
	option(ALPS_BUILD_PYTHON "Build ALPS python extentions" ON)
endif (NOT_ALPS_BUILD_PYTHON)

option(ALPS_MACOS106 "Build ALPS for MacOS X 10.6" OFF)

option(ALPS_INSTALL_DEPRECATED "Install deprecated headers for ALPS" OFF)
option(ALPS_BUILD_DEVELOPER_TOOLS "Build tools used by developers to maintain ALPS" OFF)
option(ALPS_ENABLE_OPENMP "Enable OpenMP parallelization" OFF)
option(ALPS_ENABLE_OPENMP_WORKER "Enable OpenMP worker support" OFF)
option(ALPS_BUILD_ON_CRAY "Build ALPS for cray" OFF)
option(ALPS_BUILD_SOURCE "Build a source distribution" OFF)
option(LAPACK_64_BIT "Use 64-bit version of LAPACK" OFF)
option(ALPS_BUILD_FORTRAN "Build the ALPS fortran applications and tools" OFF)
option(ALPS_INCLUDE_TUTORIALS "Include the ALPS tutorials" ON)
option(ALPS_BUILD_TUTORIALS_PACKAGE "Build only a tutorials package" OFF)
option(BIND_FORTRAN_LOWERCASE "FORTRAN functions are compiled without a trailing underscore" OFF)
option(ALPS_VISTRAILS_32BIT "Force Vistrails packages to be built in 32 bit" OFF)
option(ALPS_LINK_BOOST_TEST "Link Boost.test pre-built library" OFF)
option(ALPS_INSTALL_BOOST_TEST "Install Boost Test framework" OFF)
option(ALPS_PACKAGE_LIBRARIES "Install all dependent libraries into package" OFF)
option(ALPS_DOWNLOAD_PACKAGES "Prepare an installer that downloads instead of containing the installation packages" OFF)
option(ALPS_PATCH_VISTRAILS "Patch VisTrails from VisTrails git repository" OFF)
option(ALPS_NGS_USE_NEW_ALEA "Use the ALPS ngs-alea instead of just alea" OFF)
option(ALPS_NGS_OPENMPI_ULFM "Userlevel failure mitigation is available" OFF)

mark_as_advanced(BIND_FORTRAN_LOWERCASE)
mark_as_advanced(ALPS_BUILD_TUTORIALS_PACKAGE)
mark_as_advanced(ALPS_LINK_BOOST_TEST)
mark_as_advanced(ALPS_INSTALL_BOOST_TEST)
mark_as_advanced(ALPS_PATCH_VISTRAILS)
mark_as_advanced(ALPS_NGS_USE_NEW_ALEA)
mark_as_advanced(ALPS_NGS_OPENMPI_ULFM)

option(ALPS_USE_MKL_PARALLEL "Use parallel version of MKL" OFF)
mark_as_advanced(ALPS_USE_MKL_PARALLEL)
    
if (ALPS_FOR_VISTRAILS)
  if (NOT WIN32 AND NOT APPLE)
    message(ERROR "Use ALPS_INSTALL_VISTRAILS_PACKAGES instead of ALPS_FOR_VISTRAILS on Linux")
  endif (NOT WIN32 AND NOT APPLE)
  set(NOT_ALPS_FOR_VISTRAILS OFF)
else (ALPS_FOR_VISTRAILS)
  set(NOT_ALPS_FOR_VISTRAILS ON)
endif (ALPS_FOR_VISTRAILS)

#IF(EXISTS ${PROJECT_SOURCE_DIR}/applications)
  SET (APPLICATIONS_CAN_BE_BUILT  ${NOT_ALPS_FOR_VISTRAILS})
#ELSE(EXISTS ${PROJECT_SOURCE_DIR}/applications)
#  SET (APPLICATIONS_CAN_BE_BUILT OFF)
#ENDIF(EXISTS ${PROJECT_SOURCE_DIR}/applications)

set(ALPS_BOOST_LIBRARY_NAME "boost" CACHE STRING "name of the boost library")
set(PYTHON_EXTENSION_DIR "lib" CACHE STRING "path for Python extensions")
set(ALPS_ENABLE_MPI ${NOT_ALPS_FOR_VISTRAILS} CACHE BOOL "Enable MPI Parallelization")
set(ALPS_INSTALL_HEADERS  ${NOT_ALPS_FOR_VISTRAILS} CACHE BOOL "Install headers for ALPS and all dependent libraries")
set(ALPS_BUILD_EXAMPLES ${NOT_ALPS_FOR_VISTRAILS} CACHE BOOL "Build ALPS examples")
set(ALPS_BUILD_TESTS ${NOT_ALPS_FOR_VISTRAILS} CACHE BOOL "Build ALPS tests")
set(ALPS_BUILD_DOC ${ALPS_BUILD_SOURCE} CACHE BOOL "Build ALPS documentation")
set(ALPS_BUILD_NEW_DOC ${ALPS_BUILD_SOURCE} CACHE BOOL "Build ALPS Sphinx documentation")
set(ALPS_INCLUDE_BOOST_SOURCE OFF CACHE BOOL "Include Boost sources in the source tarball")
set(ALPS_SYSTEM_NAME "" CACHE STRING "Override the package name")
set(ALPS_PACKAGE_VERSION "" CACHE STRING "Override the package version")
set(ALPS_BUILD_APPLICATIONS ${APPLICATIONS_CAN_BE_BUILT} CACHE BOOL "Build ALPS applications")
set(ALPS_USE_VISTRAILS ${ALPS_FOR_VISTRAILS} CACHE BOOL "Use the libraries deployed with VisTrails")
set(ALPS_INSTALL_VISTRAILS_PACKAGES ${ALPS_FOR_VISTRAILS} CACHE BOOL "Install the ALPS VisTrails extensions")

mark_as_advanced(ALPS_BOOST_LIBRARY_NAME)
mark_as_advanced(PYTHON_EXTENSION_DIR)

if(CMAKE_SIZEOF_VOID_P EQUAL 8 OR CMAKE_GENERATOR MATCHES Win64)
  if (NOT ALPS_VISTRAILS_32BIT)
    set (ALPS_64BIT ON)
  endif (NOT ALPS_VISTRAILS_32BIT)
endif(CMAKE_SIZEOF_VOID_P EQUAL 8 OR CMAKE_GENERATOR MATCHES Win64)

# maintenance building options
# 
# removed CMakeDetermineSystem since it is already called by
# the makefile generator, and cmake 2.8.10 tries to cache the
# configuration in /CMakeSystem.cmake

#include (CMakeDetermineSystem)
#include (CMakeDetermineFortranCompiler)
#message (STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
#include (CMakeTestFortranCompiler)
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
if(ALPS_BUILD_FORTRAN)
  project(alps C CXX Fortran)
  if (APPLE)
   set (CMAKE_Fortran_RUNTIME_LIBRARIES "-lgcc_s.1")
  endif(APPLE)
else(ALPS_BUILD_FORTRAN)
  project(alps C CXX)
endif(ALPS_BUILD_FORTRAN)

# print extensive compiler information
if(NOT WIN32)
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} --version COMMAND head -n1 OUTPUT_VARIABLE ALPS_COMPILER_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "Compiler version: ${ALPS_COMPILER_VERSION}")
endif(NOT WIN32)

# set convenient warning flags for clang
STRING(REGEX MATCH "clang" IS_CLANG "${ALPS_COMPILER_VERSION}")
if(IS_CLANG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth-256 -Wall -Wno-unknown-pragmas -Wno-unused-function -Wno-overloaded-virtual -Wno-unneeded-internal-declaration -Wno-unused-variable -Wno-return-type-c-linkage")
    if(NOT ${CMAKE_CXX_FLAGS} MATCHES "-(std=c\\+\\+11|stdlib=)")
        # compilation fails with -stdlib=libc++, which is new default on Mavericks
        message(STATUS "Forcing -stdlib=libstdc++. If you want to override this decision explicitly add -stdlib=libc++ or -std=c++11 to CMAKE_CXX_FLAGS.")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
    endif(NOT ${CMAKE_CXX_FLAGS} MATCHES "-(std=c\\+\\+11|stdlib=)")
endif(IS_CLANG)

IF(WIN32)
ADD_DEFINITIONS(/bigobj)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /bigobj")
SET(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} /bigobj")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /bigobj")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} /bigobj")
ENDIF(WIN32)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

######################################################################
# CMAKE_MODULE_PATH
######################################################################
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

######################################################################
# set default CMAKE_INSTALL_PREFIX
######################################################################

if (ALPS_BUILD_PYTHON OR WIN32 OR APPLE)
  find_package(Vistrails)
endif (ALPS_BUILD_PYTHON OR WIN32 OR APPLE)

if (ALPS_BUILD_SOURCE)
  set (CMAKE_INSTALL_PREFIX .)
else (ALPS_BUILD_SOURCE)
  if(WIN32 AND NOT UNIX)
    if (ALPS_FOR_VISTRAILS)
      if (ALPS_64BIT)
        set(CMAKE_INSTALL_PREFIX "$ENV{HOMEDRIVE}/Program Files/VisTrails" CACHE PATH "ALPS install prefix")
      else (ALPS_64BIT)
        set(CMAKE_INSTALL_PREFIX "$ENV{HOMEDRIVE}/Program Files (x86)/VisTrails" CACHE PATH "ALPS install prefix")
      endif(ALPS_64BIT)
      set(CMAKE_INSTALL_ALTERNATE_PREFIX "$ENV{HOMEDRIVE}/Program Files/VisTrails" CACHE PATH "ALPS install prefix")
    else (ALPS_FOR_VISTRAILS)
      if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        if (ALPS_64BIT)
          set(CMAKE_INSTALL_PREFIX "$ENV{HOMEDRIVE}/Program Files/ALPS" CACHE PATH "ALPS install prefix" FORCE)
        else (ALPS_64BIT)
          set(CMAKE_INSTALL_PREFIX "$ENV{HOMEDRIVE}/Program Files (x86)/VisTrails" CACHE PATH "ALPS install prefix")
        endif(ALPS_64BIT)
        set(CMAKE_INSTALL_PREFIX "$ENV{HOMEDRIVE}/Program Files (x86)/ALPS" CACHE PATH "ALPS install prefix" FORCE)
        set(CMAKE_INSTALL_ALTERNATE_PREFIX "$ENV{HOMEDRIVE}/Program Files/ALPS" CACHE PATH "ALPS install prefix")
      endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    endif (ALPS_FOR_VISTRAILS)
  else(WIN32 AND NOT UNIX)
    if (ALPS_MACOS106)
      set(CMAKE_OSX_DEPLOYMENT_TARGET "10.6")
    endif (ALPS_MACOS106)
    set(CPACK_OSX_PACKAGE_VERSION "10.5")
    if (ALPS_FOR_VISTRAILS)
      if (APPLE)
        if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
          set(CMAKE_INSTALL_PREFIX ${VISTRAILS_APP_DIR} CACHE PATH "ALPS install prefix" FORCE)
        endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        if (ALPS_VISTRAILS_32BIT)
          set(CMAKE_OSX_DEPLOYMENT_TARGET "10.5")
          set(CMAKE_OSX_ARCHITECTURES "i386")
          set(CPACK_OSX_PACKAGE_VERSION "10.5")
        endif (ALPS_VISTRAILS_32BIT)
      else (APPLE)
        # what happens to the above variables in this case?
      endif (APPLE)
    else (ALPS_FOR_VISTRAILS)
      if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(CMAKE_INSTALL_PREFIX "/opt/alps" CACHE PATH "ALPS install prefix" FORCE)
      endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    endif (ALPS_FOR_VISTRAILS)
  endif(WIN32 AND NOT UNIX)
endif(ALPS_BUILD_SOURCE)

######################################################################
# Version information
######################################################################
set(ALPS_YEAR 2013)
set(ALPS_VERSION_MAJOR 2)
set(ALPS_VERSION_MINOR 2)
set(ALPS_VERSION_PATCH b2)
set(ALPS_VERSION_BUILD "")

IF(NOT ALPS_VERSION_BUILD AND EXISTS ${PROJECT_SOURCE_DIR}/.svn) 
  FIND_PACKAGE(Subversion) 
  IF(Subversion_FOUND) 
    Subversion_WC_INFO(${PROJECT_SOURCE_DIR} ALPS) 
    set(ALPS_VERSION_BUILD "r${ALPS_WC_REVISION}")
  ENDIF(Subversion_FOUND) 
ENDIF(NOT ALPS_VERSION_BUILD AND EXISTS ${PROJECT_SOURCE_DIR}/.svn) 

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

if(MSVC)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS)
endif(MSVC)

######################################################################
# Packages
######################################################################

# Boost Libraries
find_package(BoostForALPS REQUIRED)
if(Boost_ROOT_DIR)
  message(STATUS "Adding Boost dir: ${Boost_ROOT_DIR} " )
  include_directories(${Boost_ROOT_DIR})
else(Boost_ROOT_DIR)
  message(STATUS "Adding Boost dir: ${Boost_INCLUDE_DIR} " )
  include_directories(${Boost_INCLUDE_DIR})
endif(Boost_ROOT_DIR)

# MPI
if(ALPS_ENABLE_MPI)
  IF(NOT MPI_FOUND)
    IF(ALPS_PACKAGE_LIBRARIES AND NOT ALPS_USE_VISTRAILS AND APPLE)
      set(MPI_COMPILER ${CMAKE_INSTALL_PREFIX}/bin/mpicxx)
      message (STATUS "MPI compiler is ${MPI_COMPILER}")
    ENDIF(ALPS_PACKAGE_LIBRARIES AND NOT ALPS_USE_VISTRAILS AND APPLE)
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
#add_definitions( -DBOOST_MPI_HOMOGENEOUS )
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
    if(ALPS_ENABLE_OPENMP_WORKER)
      # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DALPS_ENABLE_OPENMP_WORKER")
    endif(ALPS_ENABLE_OPENMP_WORKER)
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
if(MSVC)
  file(COPY ${LAPACK_LIBRARY} ${BLAS_LIBRARY} DESTINATION ${PROJECT_BINARY_DIR}/lib)
endif(MSVC)
IF(NOT WIN32)
  IF (REQUIRE_PTHREAD)
    SET(LAPACK_LIBRARY "${LAPACK_LIBRARY};${PTHREAD_LIBRARY}")
    SET(LAPACK_LIBRARIES "${LAPACK_LIBRARIES};${PTHREAD_LIBRARY}")
  ENDIF (REQUIRE_PTHREAD)
ENDIF(NOT WIN32)
if(MAC_VECLIB)
  set(LAPACK_LDFLAGS "-framework vecLib")
endif(MAC_VECLIB) 
set(LAPACK_LINKER_FLAGS ${LAPACK_LDFLAGS})

# SQLite
find_package(SQLite)
if(MSVC AND SQLite_FOUND)
  file(COPY ${SQLite_LIBRARIES} DESTINATION ${PROJECT_BINARY_DIR}/lib)
  file(COPY ${SQLite_DLLS} DESTINATION ${PROJECT_BINARY_DIR}/bin)
endif(MSVC AND SQLite_FOUND)

# FFTW
find_package(FFTW)

# LPSolve
find_package(LPSolve)

# Zlib
find_package(ZLIB)
if(ZLIB_FOUND AND NOT ALPS_USE_VISTRAILS)
  include_directories(${ZLIB_INCLUDE_DIR})
  link_libraries(${ZLIB_LIBRARIES})
  if(MSVC)
    file(COPY ${ZLIB_LIBRARIES} DESTINATION ${PROJECT_BINARY_DIR}/lib)
  endif(MSVC)
endif(ZLIB_FOUND AND NOT ALPS_USE_VISTRAILS)


# SZip
find_package(SZIP)
if(SZIP_FOUND AND NOT ALPS_USE_VISTRAILS)
  include_directories(${SZIP_INCLUDE_DIR})
  link_libraries(${SZIP_LIBRARIES})
  if(MSVC)
    file(COPY ${SZIP_LIBRARIES} DESTINATION ${PROJECT_BINARY_DIR}/lib)
  endif(MSVC)
endif(SZIP_FOUND AND NOT ALPS_USE_VISTRAILS)

if(MSVC AND NOT ALPS_FOR_VISTRAILS)
  file(COPY ${SZIP_DLLS} ${ZLIB_DLLS} ${PYTHON_DLLS} DESTINATION ${PROJECT_BINARY_DIR}/bin)
endif(MSVC AND NOT ALPS_FOR_VISTRAILS)
 
# HDF5
find_package(HDF5)
IF(HDF5_FOUND)
  SET(HAVE_LIBHDF5 1)  
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  add_definitions(${HDF5_DEFINTIONS})
  SET(ALPS_HAVE_HDF5 1)
  if(MSVC)
    file(COPY ${HDF5_hdf5_LIBRARY_RELEASE} ${HDF5_hdf5_hl_LIBRARY_RELEASE}
      DESTINATION ${PROJECT_BINARY_DIR}/lib)
    file(COPY ${HDF5_DLLS} DESTINATION ${PROJECT_BINARY_DIR}/bin)
  endif(MSVC)
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

# python


if (ALPS_BUILD_PYTHON)
  set(PYTHON_SCRIPTDIR "${CMAKE_INSTALL_PREFIX}/lib/python")
  if (ALPS_USE_VISTRAILS)
    #set up the installation for VisTrail's Python
    if (APPLE)
      SET(PYTHON_INCLUDE_DIRS ${VISTRAILS_APP_DIR}/${VISTRAILS_APP_NAME}/Contents/Frameworks/Python.framework/Headers CACHE PATH "Python include directory")
      SET(PYTHON_LIBRARY "-F${VISTRAILS_APP_DIR}/${VISTRAILS_APP_NAME}/Contents/Frameworks -framework Python" CACHE STRING "Python libraries")
      SET(PYTHON_NUMPY_INCLUDE_DIR ${VISTRAILS_APP_DIR}/${VISTRAILS_PYTHON_EXTENSION_DIR}/numpy/core/include CACHE PATH "Numpy include directory")
    endif(APPLE)
    if (WIN32)
      SET(PYTHON_LIBRARY "${VISTRAILS_APP_DIR}/${VISTRAILS_PYTHON_DIR}/libs/python27.lib" CACHE STRING "Python libraries")
      SET(PYTHON_INCLUDE_DIRS ${VISTRAILS_APP_DIR}/${VISTRAILS_PYTHON_DIR}/include CACHE PATH "Python include directory")
      SET(PYTHON_NUMPY_INCLUDE_DIR ${VISTRAILS_APP_DIR}/${VISTRAILS_PYTHON_EXTENSION_DIR}/numpy/core/include CACHE PATH "Numpy include directory")
    endif(WIN32)
    if(NOT APPLE AND NOT WIN32)
      set(ALPS_USE_VISTRAILS "FALSE")
      find_package(Python)
      set(ALPS_USE_VISTRAILS "TRUE")
    endif(NOT APPLE AND NOT WIN32)
    set(PYTHONLIBS_FOUND TRUE)
    set(PYTHON_FOUND TRUE)
  endif (ALPS_USE_VISTRAILS)
  find_package(Python)
endif (ALPS_BUILD_PYTHON)

if(WIN32 AND NOT UNIX AND ALPS_BUILD_PYTHON AND NOT ALPS_FOR_VISTRAILS)
  find_package(PyInstaller)
  include(UsePyInstaller)
endif(WIN32 AND NOT UNIX AND ALPS_BUILD_PYTHON AND NOT ALPS_FOR_VISTRAILS)

IF (PYTHONLIBS_FOUND AND ALPS_BUILD_PYTHON)
  include_directories(${PYTHON_NUMPY_INCLUDE_DIR})
    MESSAGE (STATUS "Numpy include in ${PYTHON_NUMPY_INCLUDE_DIR}")
  SET(ALPS_HAVE_PYTHON ON)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
  set(BUILD_BOOST_PYTHON TRUE)
ELSE (PYTHONLIBS_FOUND AND ALPS_BUILD_PYTHON)
  set(BUILD_BOOST_PYTHON OFF)
  SET(ALPS_HAVE_PYTHON OFF)
ENDIF (PYTHONLIBS_FOUND AND ALPS_BUILD_PYTHON)

set(ALPS_XML_PATH "lib" CACHE STRING "Directory to install ALPS XML directory in")
if (ALPS_FOR_VISTRAILS)
  set(ALPS_XML_DIR "${CMAKE_INSTALL_PREFIX}/${VISTRAILS_APP_NAME}/${VISTRAILS_LIB_DIR}/xml")
  if (WIN32 AND NOT UNIX)
    set(ALPS_XML_ALTERNATE_DIR "${CMAKE_INSTALL_ALTERNATE_PREFIX}/${VISTRAILS_APP_NAME}/${VISTRAILS_LIB_DIR}/xml")
  endif (WIN32 AND NOT UNIX)
else (ALPS_FOR_VISTRAILS)
  set(ALPS_XML_DIR "${CMAKE_INSTALL_PREFIX}/${ALPS_XML_PATH}/xml")
  if (WIN32 AND NOT UNIX)
    set(ALPS_XML_ALTERNATE_DIR "${CMAKE_INSTALL_ALTERNATE_PREFIX}/${ALPS_XML_PATH}/xml")
  endif (WIN32 AND NOT UNIX)
endif (ALPS_FOR_VISTRAILS)

message (STATUS "ALPS XML dir is ${ALPS_XML_DIR}")


if (BIND_FORTRAN_LOWERCASE)
  add_definitions(-DBIND_FORTRAN_LOWERCASE)
endif (BIND_FORTRAN_LOWERCASE)
  
  
######################################################################
# test environment
######################################################################
if (NOT ALPS_BUILD_SOURCE)
include(add_alps_test)
endif (NOT ALPS_BUILD_SOURCE)

######################################################################
# configuration files
######################################################################

set (BIND_FORTRAN_INTEGER_8 ${LAPACK_64_BIT})

configure_file(src/alps/config.h.in ${PROJECT_BINARY_DIR}/src/alps/config.h)
configure_file(src/alps/version.h.in ${PROJECT_BINARY_DIR}/src/alps/version.h)

configure_file(lib/xml/lattices.xml.in ${PROJECT_BINARY_DIR}/lib/xml/lattices.xml)
configure_file(lib/xml/models.xml.in ${PROJECT_BINARY_DIR}/lib/xml/models.xml)

configure_file(vistrails/alps/__init__.py.in ${PROJECT_BINARY_DIR}/vistrails/alps/__init__.py)

######################################################################
# include directories
######################################################################
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(src)

# order the include directories so that directories which are in the
# source or build tree always come before directories outside the project.
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)

######################################################################
# RPATH setting
######################################################################
if(APPLE)
  set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else(APPLE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif(APPLE)

######################################################################
# install dependent libraries under Windows, from applications
######################################################################

if(WIN32 AND NOT UNIX AND NOT ALPS_FOR_VISTRAILS AND NOT ALPS_BUILD_SOURCE)
  # Lapack
  if(LAPACK_LIBRARIES)
    foreach(file ${LAPACK_LIBRARIES})
      get_filename_component(LAPACK_PATH ${file} PATH)
    endforeach(file ${LAPACK_LIBRARIES})
    add_custom_target(liblapack ALL
      COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/lib
      COMMAND cmake -E copy_if_different ${LAPACK_PATH}/lapack.lib ${PROJECT_BINARY_DIR}/lib
      COMMAND cmake -E copy_if_different ${LAPACK_PATH}/blas.lib ${PROJECT_BINARY_DIR}/lib
      #COMMAND cmake -E copy_if_different ${LAPACK_PATH}/f2c.lib ${PROJECT_BINARY_DIR}/lib
    )
    install(FILES ${PROJECT_BINARY_DIR}/lib/lapack.lib ${PROJECT_BINARY_DIR}/lib/blas.lib
            DESTINATION lib COMPONENT dependentlibraries
#            ${PROJECT_BINARY_DIR}/lib/f2c.lib
    )
  endif(LAPACK_LIBRARIES)
  # SQLite
  if(SQLite_LIBRARIES)
    get_filename_component(SQLite_PATH ${SQLite_LIBRARIES} PATH)
    get_filename_component(SQLite_BASE ${SQLite_PATH} PATH)
    get_filename_component(SQLite_NAME ${SQLite_LIBRARIES} NAME_WE)
    add_custom_target(libsqlite ALL
      COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/lib
      COMMAND cmake -E copy_if_different ${SQLite_LIBRARIES} ${PROJECT_BINARY_DIR}/lib
      COMMAND cmake -E copy_if_different ${SQLite_BASE}/bin/${SQLite_NAME}.dll ${PROJECT_BINARY_DIR}/bin
    )
    install(FILES ${PROJECT_BINARY_DIR}/lib/${SQLite_NAME}.lib DESTINATION lib COMPONENT dependentlibraries)
    install(FILES ${PROJECT_BINARY_DIR}/bin/${SQLite_NAME}.dll DESTINATION bin COMPONENT dependentlibraries)
  endif(SQLite_LIBRARIES)
  # alps & boost
  # add_custom_target(libalps ALL
    # COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/bin
    # COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/lib
    # COMMAND cmake -E make_directory ${PROJECT_BINARY_DIR}/lib/xml
    # COMMAND cmake -E copy_if_different ${CMAKE_INSTALL_PREFIX}/bin/archive.exe ${PROJECT_BINARY_DIR}/bin
    # COMMAND cmake -E copy_if_different ${CMAKE_INSTALL_PREFIX}/bin/parameter2xml.exe ${PROJECT_BINARY_DIR}/bin
    # COMMAND cmake -E copy_if_different ${CMAKE_INSTALL_PREFIX}/bin/parameter2hdf5.exe ${PROJECT_BINARY_DIR}/bin
    # COMMAND cmake -E copy_if_different ${CMAKE_INSTALL_PREFIX}/bin/ising.exe ${PROJECT_BINARY_DIR}/bin
    # COMMAND cmake -E copy_directory ${CMAKE_INSTALL_PREFIX}/lib/xml ${PROJECT_BINARY_DIR}/lib/xml    
  # )
  if (ALPS_USE_VISTRAILS)
    install(FILES ${VISTRAILS_APP_DIR}/${VISTRAILS_APP_NAME}/${VISTRAILS_PYTHON_DIR}/python27.dll DESTINATION bin COMPONENT dependentlibraries)
    install(FILES ${VISTRAILS_APP_DIR}/${VISTRAILS_APP_NAME}/msvcr71.dll DESTINATION bin COMPONENT dependentlibraries)
    install(FILES ${VISTRAILS_APP_DIR}/${VISTRAILS_APP_NAME}/msvcr80.dll DESTINATION bin COMPONENT dependentlibraries)
    install(FILES ${VISTRAILS_APP_DIR}/${VISTRAILS_APP_NAME}/hdf5dll.dll DESTINATION bin COMPONENT dependentlibraries)
    install(FILES ${VISTRAILS_APP_DIR}/${VISTRAILS_APP_NAME}/szlibdll.dll DESTINATION bin COMPONENT dependentlibraries)
    install(FILES ${VISTRAILS_APP_DIR}/${VISTRAILS_APP_NAME}/zlib.dll DESTINATION bin COMPONENT dependentlibraries)
  endif (ALPS_USE_VISTRAILS)

endif(WIN32 AND NOT UNIX AND NOT ALPS_FOR_VISTRAILS AND NOT ALPS_BUILD_SOURCE)

######################################################################
# setting up ALPSConfig.cmake and include.mk
######################################################################

set(ALPS_USE_FILE ${CMAKE_INSTALL_PREFIX}/share/alps/UseALPS.cmake)

set(ALPS_INCLUDE_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}/include)
# if(NOT INSTALL_BOOST_HEADERS)
  list(APPEND ALPS_INCLUDE_DIRS_CONFIG ${Boost_INCLUDE_DIR})
# endif(NOT INSTALL_BOOST_HEADERS)
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
#if(SQLite_INCLUDE_DIR)
#  list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${SQLite_INCLUDE_DIR})
#endif(SQLite_INCLUDE_DIR)
if(HDF5_INCLUDE_DIR)
  list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
endif(HDF5_INCLUDE_DIR)
if(PYTHON_INCLUDE_DIRS)
  list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})
endif(PYTHON_INCLUDE_DIRS)
if(PYTHON_NUMPY_INCLUDE_DIR)
  list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${PYTHON_NUMPY_INCLUDE_DIR})
endif(PYTHON_NUMPY_INCLUDE_DIR)
# if(Boost_ROOT_DIR)
#   list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${Boost_ROOT_DIR}) # FindBoostSrc.cmake is setting Boost_INCLUDE_DIR=Boost_ROOT_DIR
# endif(Boost_ROOT_DIR)
if(Boost_INCLUDE_DIR)
  list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${Boost_INCLUDE_DIR})
endif(Boost_INCLUDE_DIR)

if(NOT Boost_FOUND)
	set(Boost_LIBRARIES ${ALPS_BOOST_LIBRARY_NAME})
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
#if(SQLite_LIBRARIES)
#  list(APPEND ALPS_EXTRA_LIBRARIES ${SQLite_LIBRARIES})
#endif(SQLite_LIBRARIES)
if(HDF5_LIBRARIES)
  list(APPEND ALPS_EXTRA_LIBRARIES ${HDF5_LIBRARIES})
endif(HDF5_LIBRARIES)
if(PYTHON_LIBRARY)
  list(APPEND ALPS_EXTRA_LIBRARIES ${PYTHON_LIBRARY})
endif(PYTHON_LIBRARY)

set (HDF5_LIBS "")
foreach (l ${HDF5_LIBRARIES})
  set (HDF5_LIBS "${HDF5_LIBS} ${l}")
  if (${l} STREQUAL "debug")
    set (HDF5_LIBS "")
  endif()
  if (${l} STREQUAL "optimized")
    set (HDF5_LIBS "")
  endif()
endforeach(l)

set (LAPACK_LIBS "")
foreach (l ${LAPACK_LIBRARY})
  set (LAPACK_LIBS "${LAPACK_LIBS} ${l}")
endforeach(l)
 if(PYTHON_FOUND)
  set(PYTHON_CPP_FLAGS "-I${PYTHON_INCLUDE_DIRS} -I${PYTHON_NUMPY_INCLUDE_DIR}")
endif(PYTHON_FOUND)
configure_file(config/ALPSConfig.cmake.in ${PROJECT_BINARY_DIR}/config/ALPSConfig.cmake @ONLY)
configure_file(config/include.mk.in ${PROJECT_BINARY_DIR}/config/include.mk)

set (HDF5_LIBRARIES ${HDF5_LIBRARIES_SAVE})


# installation
######################################################################
include(InstallRequiredSystemLibraries)

if (NOT ALPS_BUILD_TUTORIALS_PACKAGE)
if (ALPS_BUILD_SOURCE)
  FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.txt")
  INSTALL(FILES ${files} DESTINATION alps)  
  INSTALL(FILES script/hdf5.sh script/cmake.sh DESTINATION alps/script)  
  install(DIRECTORY config example legacy lib src test tool vistrails DESTINATION alps
          FILES_MATCHING PATTERN "*"  PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
  if (ALPS_INCLUDE_BOOST_SOURCE)
    install(DIRECTORY ${Boost_ROOT_DIR} DESTINATION . PATTERN "*" 
            PATTERN *fop-* EXCLUDE PATTERN "docbook-[dx]*" EXCLUDE)
  endif (ALPS_INCLUDE_BOOST_SOURCE)
else (ALPS_BUILD_SOURCE)
  if (ALPS_INSTALL_HEADERS)
    set(ALPS_HEADER_DIR "include")
    install(DIRECTORY src/alps src/boost src/ietl src/mocasito COMPONENT headers DESTINATION ${ALPS_HEADER_DIR}
            FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.ipp"
            PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
    install(FILES ${PROJECT_BINARY_DIR}/src/alps/config.h ${PROJECT_BINARY_DIR}/src/alps/version.h DESTINATION ${ALPS_HEADER_DIR}/alps COMPONENT headers)
    if (ALPS_INSTALL_DEPRECATED)
      install(DIRECTORY legacy/alps COMPONENT deprecated DESTINATION ${ALPS_HEADER_DIR}
              FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" 
              PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
    endif (ALPS_INSTALL_DEPRECATED)

    if(ALPS_PACKAGE_LIBRARIES)
       install(DIRECTORY ${Boost_ROOT_DIR}/boost DESTINATION ${ALPS_HEADER_DIR} COMPONENT dependentheaders
               FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.ipp"
               PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
       install(DIRECTORY ${PROJECT_BINARY_DIR}/deplibs/include
                         ${PROJECT_BINARY_DIR}/deplibs/etc
                         ${PROJECT_BINARY_DIR}/deplibs/share
               DESTINATION . COMPONENT dependentheaders)
       #install(DIRECTORY  DESTINATION . COMPONENT dependentheaders)
    endif(ALPS_PACKAGE_LIBRARIES)

  endif(ALPS_INSTALL_HEADERS)

  if (ALPS_FOR_VISTRAILS)
    install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/xml DESTINATION ${VISTRAILS_APP_NAME}/${VISTRAILS_LIB_DIR} COMPONENT vistrails)
  endif (ALPS_FOR_VISTRAILS)

  # if (APPLE AND NOT ALPS_FOR_VISTRAILS)
  #   configure_file(tool/postflight.sh.in ${PROJECT_BINARY_DIR}/tool/postflight.sh)
  #   file(WRITE ${PROJECT_BINARY_DIR}/tool/paths "${CMAKE_INSTALL_PREFIX}/bin")
  #   install(FILES ${PROJECT_BINARY_DIR}/tool/paths DESTINATION share/alps COMPONENT environment)
  # endif (APPLE AND NOT ALPS_FOR_VISTRAILS)

  if(ALPS_PACKAGE_LIBRARIES)
    install(DIRECTORY ${PROJECT_BINARY_DIR}/deplibs/lib DESTINATION . COMPONENT dependentlibraries)
    FILE(GLOB depfiles "${PROJECT_BINARY_DIR}/deplibs/bin/*")
    install(PROGRAMS ${depfiles} DESTINATION bin COMPONENT dependentlibraries)
    if(NOT ALPS_FOR_VISTRAILS AND MSVC)
      install(PROGRAMS ${HDF5_DLLS} ${PYTHON_DLLS} ${ZLIB_DLLS} ${SZIP_DLLS} DESTINATION bin COMPONENT dependentlibraries)
    endif(NOT ALPS_FOR_VISTRAILS AND MSVC)
  endif(ALPS_PACKAGE_LIBRARIES)

  if (ALPS_FOR_VISTRAILS)
    install(DIRECTORY lib/xml DESTINATION ${VISTRAILS_APP_NAME}/${VISTRAILS_LIB_DIR}  COMPONENT vistrails
          FILES_MATCHING PATTERN "*.xsl"
          PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)

    install(DIRECTORY lib/pyalps DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR} COMPONENT vistrails
          FILES_MATCHING PATTERN "*.py" PATTERN "*.pyc"
          PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE PATTERN "*.C" EXCLUDE)

      if (MSVC)
        install(PROGRAMS ${HDF5_DLLS} ${ZLIB_DLLS} ${SZIP_DLLS} DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps  COMPONENT vistrails)
      endif(MSVC)

    install(DIRECTORY vistrails/alps DESTINATION ${VISTRAILS_PACKAGE_DIR} COMPONENT vistrails
        FILES_MATCHING PATTERN "*.py"  PATTERN "*.xml" PATTERN "*.pyc"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
    install(FILES ${PROJECT_BINARY_DIR}/vistrails/alps/__init__.py DESTINATION ${VISTRAILS_PACKAGE_DIR}/alps COMPONENT vistrails)
  endif (ALPS_FOR_VISTRAILS)

  if (NOT ALPS_FOR_VISTRAILS)
    if(ALPS_INSTALL_VISTRAILS_PACKAGES)
      install(DIRECTORY vistrails/alps DESTINATION  ${VISTRAILS_APP_DIR}/${VISTRAILS_PACKAGE_DIR} COMPONENT vistrails
          FILES_MATCHING PATTERN "*.py"  PATTERN "*.xml" PATTERN "*.pyc"
          PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
      install(FILES ${PROJECT_BINARY_DIR}/vistrails/alps/__init__.py DESTINATION ${VISTRAILS_APP_DIR}/${VISTRAILS_PACKAGE_DIR}/alps COMPONENT vistrails)

      install(DIRECTORY lib/pyalps DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR} COMPONENT vistrails
              FILES_MATCHING PATTERN "*.py" PATTERN "*.pyc"
              PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE PATTERN "*.C" EXCLUDE)

    endif(ALPS_INSTALL_VISTRAILS_PACKAGES)

    install(DIRECTORY lib/xml DESTINATION ${ALPS_XML_PATH} COMPONENT xml
          FILES_MATCHING PATTERN "*.xsl"
          PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)

    install(DIRECTORY lib/pyalps DESTINATION lib  COMPONENT python
          FILES_MATCHING PATTERN "*.py" PATTERN "*.pyc"
          PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
          
    install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/xml DESTINATION ${ALPS_XML_PATH} COMPONENT xml)

    add_subdirectory(config)

    install(FILES config/UseALPS.cmake
            ${PROJECT_BINARY_DIR}/config/ALPSConfig.cmake
            ${PROJECT_BINARY_DIR}/config/include.mk
            config/run_test.cmake
            config/run_test_mpi.cmake
            config/add_alps_test.cmake
            DESTINATION share/alps COMPONENT build)
            
    install(FILES LICENSE.txt README.txt  LICENSE-applications.txt CITATIONS.txt
            DESTINATION share/alps COMPONENT libraries)
            
    if (VISTRAILS_FOUND AND (APPLE OR WIN32))
      message (STATUS "VisTrails found.")
      set(PYTHONEXEC "${CMAKE_INSTALL_PREFIX}/bin/vispython")
    else (VISTRAILS_FOUND AND (APPLE OR WIN32))
      set(PYTHONEXEC python)
    endif (VISTRAILS_FOUND AND (APPLE OR WIN32))

    install(FILES vistrails/latex/includevistrail.py vistrails/latex/alpsvistrails.sty 
            DESTINATION share/alps/latex COMPONENT vistrails)

  endif (NOT ALPS_FOR_VISTRAILS)
endif (ALPS_BUILD_SOURCE)
endif (NOT ALPS_BUILD_TUTORIALS_PACKAGE)

#
######################################################################
# patch VisTrailss
######################################################################

if(ALPS_PATCH_VISTRAILS)
  find_package(VistrailsGit)
  MESSAGE (STATUS "patching VisTrails from ${VISTRAILS_GIT_DIR}")
  if (APPLE)
    foreach (DIR api core db gui packages tests)
      install(DIRECTORY ${VISTRAILS_GIT_DIR}/vistrails  COMPONENT vistrails
              DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}
              PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE PATTERN "*.pyc" EXCLUDE)
    endforeach(DIR)
#    install(FILES ${VISTRAILS_GIT_DIR}/vistrails/vistrails.py  COMPONENT vistrails
#            DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}/../..)
  else(APPLE)
    if (WIN32) 
        install(DIRECTORY ${VISTRAILS_GIT_DIR}/vistrails  COMPONENT vistrails
                DESTINATION .
                PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE PATTERN "*.pyc" EXCLUDE)
    else (WIN32)
      message(ERROR "Cannot patch VisTrails on this system. Just install the latest Vistrails nightly snapshot.")
    endif (WIN32)
  endif(APPLE)
endif(ALPS_PATCH_VISTRAILS)


######################################################################
# documentation
######################################################################

if (ALPS_BUILD_DOC AND NOT ALPS_BUILD_TUTORIALS_PACKAGE)
add_subdirectory(doc)
endif (ALPS_BUILD_DOC AND NOT ALPS_BUILD_TUTORIALS_PACKAGE)

if (ALPS_BUILD_NEW_DOC AND NOT ALPS_BUILD_TUTORIALS_PACKAGE)
add_subdirectory(docs)
endif (ALPS_BUILD_NEW_DOC AND NOT ALPS_BUILD_TUTORIALS_PACKAGE)


######################################################################
# libraries
######################################################################

if (NOT ALPS_BUILD_TUTORIALS_PACKAGE)
  add_subdirectory(src/boost)
  if (NOT ALPS_BUILD_SOURCE)
    add_subdirectory(src/alps)
    add_subdirectory(lib/pyalps)
  endif (NOT ALPS_BUILD_SOURCE)
endif (NOT ALPS_BUILD_TUTORIALS_PACKAGE)


######################################################################
# programs
######################################################################

if (NOT ALPS_BUILD_TUTORIALS_PACKAGE)
  add_subdirectory(tool)
  if (ALPS_BUILD_EXAMPLES AND NOT ALPS_BUILD_SOURCE)
    add_subdirectory(example)
  endif (ALPS_BUILD_EXAMPLES AND NOT ALPS_BUILD_SOURCE)
endif (NOT ALPS_BUILD_TUTORIALS_PACKAGE)


######################################################################
# applications
######################################################################

if (ALPS_BUILD_APPLICATIONS AND NOT ALPS_BUILD_TUTORIALS_PACKAGE)
  add_subdirectory(applications)
endif (ALPS_BUILD_APPLICATIONS AND NOT ALPS_BUILD_TUTORIALS_PACKAGE)

if (NOT ALPS_FOR_VISTRAILS AND ALPS_INCLUDE_TUTORIALS)
  add_subdirectory(tutorials)
endif (NOT ALPS_FOR_VISTRAILS AND ALPS_INCLUDE_TUTORIALS)

######################################################################
# developer tools
######################################################################

if (ALPS_BUILD_DEVELOPER_TOOLS AND NOT ALPS_BUILD_SOURCE)
  add_subdirectory(script)
endif (ALPS_BUILD_DEVELOPER_TOOLS AND NOT ALPS_BUILD_SOURCE)

##########################################################################
# Binary packages                                                        #
##########################################################################

if (NOT WIN32)
  set (CPACK_SET_DESTDIR ON)
endif (NOT WIN32)


if (NOT ALPS_BUILD_TUTORIALS_PACKAGE)
 set(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "ALPS and Boost libraries")
 set(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "ALPS headers")
 set(CPACK_COMPONENT_DEPRECATED_DISPLAY_NAME "Deprecated ALPS headers")
 set(CPACK_COMPONENT_APPLICATIONS_DISPLAY_NAME "ALPS applications")
 set(CPACK_COMPONENT_XML_DISPLAY_NAME "XML files")
 set(CPACK_COMPONENT_PYTHON_DISPLAY_NAME "ALPS Python packages")
 set(CPACK_COMPONENT_TOOLS_DISPLAY_NAME "ALPS tools")
 set(CPACK_COMPONENT_DEPENDENTLIBRARIES_DISPLAY_NAME "Dependent libraries")
 set(CPACK_COMPONENT_DEPENDENTHEADERS_DISPLAY_NAME "Headers of dependent libraries")
 set(CPACK_COMPONENT_EXAMPLES_DISPLAY_NAME "ALPS examples")
 set(CPACK_COMPONENT_BUILD_DISPLAY_NAME "Build system support")
 set(CPACK_COMPONENT_DOCUMENTATION_DISPLAY_NAME "Documentation")
 set(CPACK_COMPONENT_TUTORIALS_DISPLAY_NAME "Tutorials")
 set(CPACK_COMPONENT_VISTRAILS_DISPLAY_NAME "ALPS VisTrails extensions")
 
 set(CPACK_COMPONENT_LIBRARIES_DESCRIPTION "The libraries needed by the ALPS applications")
 set(CPACK_COMPONENT_HEADERS_DESCRIPTION "The ALPS header files")
 set(CPACK_COMPONENT_DEPRECATED_DESCRIPTION "Deprecated ALPS header files")
 set(CPACK_COMPONENT_APPLICATIONS_DESCRIPTION "The ALPS application programs")
 set(CPACK_COMPONENT_XML_DESCRIPTION "The XML files and XSL stylefiles for ALPS input and evaluation")
 set(CPACK_COMPONENT_PYTHON_DESCRIPTION "The ALPS Python packages for data evaluation")
 set(CPACK_COMPONENT_TOOLS_DESCRIPTION "The ALPS tools")
 set(CPACK_COMPONENT_DEPENDENTLIBRARIES_DESCRIPTION "External libraries that ALPS depends on. They do not need to be installed if they already exist on the system and can be found by the ALPS applications")
 set(CPACK_COMPONENT_DEPENDENTHEADERS_DESCRIPTION "Headers of libraries that ALPS depends on")
 set(CPACK_COMPONENT_EXAMPLES_DESCRIPTION "ALPS examples")
 set(CPACK_COMPONENT_BUILD_DESCRIPTION "Include files for CMake and make to simplify the building of ALPS applications")
 set(CPACK_COMPONENT_DOCUMENTATION_DESCRIPTION "Documentation for the ALPS libraries")
 set(CPACK_COMPONENT_TUTORIALS_DESCRIPTION "Tutorials for the ALPS libraries and applications")
 set(CPACK_COMPONENT_VISTRAILS_DESCRIPTION "VisTrails extensions and patches to the latest VisTrails release needed for ALPS")

 set(CPACK_COMPONENT_LIBRARIES_GROUP "Runtime")
 set(CPACK_COMPONENT_APPLICATIONS_GROUP "Runtime")
 set(CPACK_COMPONENT_XML_GROUP "Runtime")
 set(CPACK_COMPONENT_PYTHON_GROUP "Runtime")
 set(CPACK_COMPONENT_TOOLS_GROUP "Runtime")
 set(CPACK_COMPONENT_DEPENDENTLIBRARIES_GROUP "Runtime")
 set(CPACK_COMPONENT_HEADERS_GROUP "Development")
 set(CPACK_COMPONENT_DEPRECATED_GROUP "Development")
 set(CPACK_COMPONENT_DEPENDENTHEADERS_GROUP "Development")
 set(CPACK_COMPONENT_BUILD_GROUP "Development")
 #set(CPACK_COMPONENT_EXAMPLES_GROUP "Examplecodes")

 #set(CPACK_COMPONENT_HEADERS_DISABLED ON)
 set(CPACK_COMPONENT_DEPRECATED_DISABLED ON)
 #set(CPACK_COMPONENT_DEPENDENTHEADERS_DISABLED ON)
 #set(CPACK_COMPONENT_BUILD_DISABLED ON)

 set(CPACK_COMPONENT_GROUP_DEVELOPMENT_DESCRIPTION "Components needed only to build applications based on the ALPS libraries")
 set(CPACK_COMPONENT_GROUP_RUNTIME_DESCRIPTION "Components needed to run the ALPS applications")
# set(CPACK_COMPONENT_GROUP_EXAMPLECODES_DESCRIPTION "Installs some example programs")
 set(CPACK_COMPONENT_APPLICATIONS_DEPENDS libraries)
 set(CPACK_COMPONENT_PYTHON_DEPENDS libraries)
 set(CPACK_COMPONENT_TOOLS_DEPENDS libraries)
 set(CPACK_COMPONENT_EXAMPLES_DEPENDS libraries)
 set(CPACK_COMPONENT_TUTORIALS_DEPENDS python applications libraries tools)
 
 # if (APPLE AND NOT ALPS_FOR_VISTRAILS)
 #  set(CPACK_COMPONENT_ENVIRONMENT_DISPLAY_NAME "ALPS PATH environment")
 #  set(CPACK_COMPONENT_ENVIRONMENT_DESCRIPTION "Set PATH variables to find ALPS applications in the Terminal")
 #  set(CPACK_POSTFLIGHT_SCRIPT ${PROJECT_BINARY_DIR}/tool/postflight.sh)
 # endif (APPLE AND NOT ALPS_FOR_VISTRAILS)
 
endif (NOT ALPS_BUILD_TUTORIALS_PACKAGE)
 
if (ALPS_FOR_VISTRAILS)
  set(CPACK_PACKAGE_NAME "alps-vistrails")
  set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY OFF)
  set(CPACK_NSIS_MODIFY_PATH OFF)
else(ALPS_FOR_VISTRAILS)
  if (ALPS_BUILD_TUTORIALS_PACKAGE)
    set(CPACK_PACKAGE_NAME "alps-tutorials")
  elseif (ALPS_BUILD_TUTORIALS_PACKAGE)
    if (BUILD_SHARED_LIBS)
      set(CPACK_PACKAGE_NAME "alps")
    elseif (BUILD_SHARED_LIBS)
      set(CPACK_PACKAGE_NAME "alps-static")
    endif (BUILD_SHARED_LIBS)

    set(CPACK_NSIS_MODIFY_PATH ON)
    if(ALPS_DOWNLOAD_PACKAGES)
      set(CPACK_PACKAGE_NAME "${CPACK_PACKAGE_NAME}-download")
    endif(ALPS_DOWNLOAD_PACKAGES)
  endif (ALPS_BUILD_TUTORIALS_PACKAGE)
endif(ALPS_FOR_VISTRAILS)

set(CPACK_PACKAGE_VENDOR "ALPS Collaboration")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ALPS ${ALPS_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/README.txt")
set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.txt")

if (ALPS_PACKAGE_LIBRARIES)
  set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README-package.txt")
  set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE-package.txt")
else (ALPS_PACKAGE_LIBRARIES)
  set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.txt")
  if (ALPS_BUILD_APPLCATIONS)
    set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE-applications.txt")
  else (ALPS_BUILD_APPLCATIONS)
    set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE.txt")
  endif (ALPS_BUILD_APPLCATIONS)
endif (ALPS_PACKAGE_LIBRARIES)

set(CPACK_RESOURCE_FILE_WELCOME "${PROJECT_SOURCE_DIR}/Welcome.txt")
set(CPACK_PACKAGE_VERSION "${ALPS_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${ALPS_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${ALPS_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_SUBMINOR "${ALPS_VERSION_SUBMINOR}")
if(ALPS_VERSION_BUILD)
  set(CPACK_PACKAGE_VERSION_PATCH "${ALPS_VERSION_PATCH}-${ALPS_VERSION_BUILD}")
else(ALPS_VERSION_BUILD)
  set(CPACK_PACKAGE_VERSION_PATCH "${ALPS_VERSION_PATCH}")
endif(ALPS_VERSION_BUILD)

if (ALPS_FOR_VISTRAILS)
  set(CPACK_PACKAGE_INSTALL_DIRECTORY "VisTrails")
else (ALPS_FOR_VISTRAILS)
  set(CPACK_PACKAGE_INSTALL_DIRECTORY "ALPS")
endif (ALPS_FOR_VISTRAILS)


set(CPACK_PACKAGING_INSTALL_PREFIX "/")
#set(CPACK_PACKAGE_DEFAULT_LOCATION "${CMAKE_INSTALL_PREFIX}")

# override package name manually
if (ALPS_SYSTEM_NAME)
  set (CPACK_SYSTEM_NAME ${ALPS_SYSTEM_NAME})
else (ALPS_SYSTEM_NAME)
  if(ALPS_BUILD_SOURCE)
    if (ALPS_INCLUDE_BOOST_SOURCE)
      set (CPACK_SYSTEM_NAME "src-with-boost")
    else (ALPS_INCLUDE_BOOST_SOURCE)
      set (CPACK_SYSTEM_NAME "src")
    endif (ALPS_INCLUDE_BOOST_SOURCE)
  endif(ALPS_BUILD_SOURCE)
endif (ALPS_SYSTEM_NAME)
if (ALPS_PACKAGE_VERSION)
  set (CPACK_PACKAGE_VERSION ${ALPS_PACKAGE_VERSION})
endif (ALPS_PACKAGE_VERSION)
  
if (ALPS_FOR_VISTRAILS)
set(CPACK_PACKAGE_DEFAULT_LOCATION "/")
else (ALPS_FOR_VISTRAILS)
set(CPACK_PACKAGE_DEFAULT_LOCATION ${CMAKE_INSTALL_PREFIX})
endif(ALPS_FOR_VISTRAILS)

set(CPACK_PACKAGE_RELOCATABLE ON)
# set(CPACK_MONOLITHIC_INSTALL ON)

# For Debian package
#  See config/debian/sid/*

# For RPM package
# set(CPACK_RPM_PACKAGE_SUMMARY "${CPACK_PACKAGE_DESCRIPTION}")
# set(CPACK_RPM_PACKAGE_NAME "${CPACK_PACKAGE_NAME}")
# set(CPACK_RPM_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION}")
set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
# set(CPACK_RPM_PACKAGE_RELEASE 1)
set(CPACK_RPM_PACKAGE_LICENCE "ALPS License")
set(CPACK_RPM_PACKAGE_GROUP "${CPACK_PACKAGE_VENDOR}")
# set(CPACK_RPM_PACKAGE_VENDOR "${CPACK_PACKAGE_VENDOR}")
# set(CPACK_RPM_PACKAGE_DESCRIPTION "${CPACK_PACKAGE_DESCRIPTION_FILE}")
# set(CPACK_RPM_PACKAGE_REQUIRES "")
# set(CPACK_RPM_SPEC_INSTALL_POST "")
# set(CPACK_RPM_SPEC_MORE_DEFINE "")
# set(CPACK_RPM_USER_SPECFILE "")
# set(CPACK_RPM_GENERATE_USER_BINARY_SPECFILE_TEMPLATE "")
# set(CPACK_RPM_PRE_INSTALL_SCRIPT_FILE "")
# set(CPACK_RPM_POST_INSTALL_SCRIPT_FILE "")
# set(CPACK_RPM_PRE_UNINSTALL_SCRIPT_FILE "")
# set(CPACK_RPM_POST_UNINSTALL_SCRIPT_FILE "")

if(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  # NOTE: No ALPS icon yet
#  set(CPACK_PACKAGE_ICON "${PROJECT_SOURCE_DIR}/tools/build/CMake\\\\InstallIcon.bmp")
#  set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\MyExecutable.exe")
#  set(CPACK_NSIS_DISPLAY_NAME "ALPS ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_SUBMINOR}${CPACK_PACKAGE_VERSION_PATCH}")
  set(CPACK_NSIS_DISPLAY_NAME "ALPS ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_SUBMINOR}")
  set(CPACK_NSIS_CONTACT "wistaria@comp-phys.org")
  set(CPACK_NSIS_HELP_LINK "http://alps.comp-phys.org")
  set(CPACK_NSIS_URL_INFO_ABOUT "http://alps.comp-phys.org")
  if (ALPS_64BIT)
    set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES64")
  endif (ALPS_64BIT)
  if (NOT ALPS_FOR_VISTRAILS)
    set(CPACK_NSIS_MODIFY_PATH ON)
  endif (NOT ALPS_FOR_VISTRAILS)
  
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${CPACK_NSIS_DISPLAY_NAME}")
endif(WIN32 AND NOT UNIX)

#set(CPACK_STRIP_FILES ON)
include(CPack)

if (ALPS_DOWNLOAD_PACKAGES)
  if(COMMAND cpack_configure_downloads)
    cpack_configure_downloads(
    "http://alps.comp-phys.org/static/software/releases/packages/"
      ALL ADD_REMOVE)
  endif()
endif(ALPS_DOWNLOAD_PACKAGES)

##########################################################################


######################################################################
# tests
######################################################################

if (ALPS_BUILD_TESTS AND NOT ALPS_BUILD_SOURCE)
  include(ALPSBuildname)
  enable_testing()
  include(CTest)
  add_subdirectory(test)
endif (ALPS_BUILD_TESTS AND NOT ALPS_BUILD_SOURCE)
