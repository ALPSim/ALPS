##########################################################################
# CMake Build Rules for ALPS                                             #
##########################################################################
# Copyright (C) 2008 Matthias Troyer <troyer@comp-phys.org>              #
# based on the Boost CMake files                                         #
# Copyright (C) 2007, 2008 Douglas Gregor <doug.gregor@gmail.com>        #
# Copyright (C) 2007 Troy Straszheim                                     #
#                                                                        #
# Distributed under the Boost Software License, Version 1.0.             #
# See accompanying file LICENSE_1_0.txt or copy at                       #
#   http://alps.comp-phys.org/LICENSE_1_0.txt                            #
##########################################################################
# Basic Usage:                                                           #
#                                                                        #
#   On Unix variants:                                                    #
#     ccmake ALPS_DIRECTORY                                              #
#                                                                        #
#     (c)onfigure options to your liking, then (g)enerate                #
#     makefiles. Use "make" to build, "make test" to test, "make         #
#     install" to install, and "make package" to build binary            #
#     packages.                                                          #
#                                                                        #
#   On Windows:                                                          #
#     run the CMake GNU, load the ALPS directory, and generate           #
#     project files or makefiles for your environment.                   #
#                                                                        #
# For more information about CMake, see http://www.cmake.org             #
##########################################################################
cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)
project(alps)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(ALPS_LIGHT "Build ALPS light components only" OFF)

######################################################################
# Version information
######################################################################
set(ALPS_VERSION_MAJOR 1)
set(ALPS_VERSION_MINOR 4)
set(ALPS_VERSION_SUBMINOR a1)
set(ALPS_VERSION "${ALPS_VERSION_MAJOR}.${ALPS_VERSION_MINOR}.${ALPS_VERSION_SUBMINOR}")

######################################################################
#define the paths for library and executable for out-of-soruce compilation
######################################################################
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------- Setup output Directories -------------------------
SET (LIBRARY_OUTPUT_PATH
  ${PROJECT_BINARY_DIR}/lib
  CACHE PATH
  "Directory for all Libraries"
  )

# --------- Setup the Executable output Directory -------------
SET (EXECUTABLE_OUTPUT_PATH
  ${PROJECT_BINARY_DIR}/bin
  CACHE PATH
  "Directory for all Executables."
  )

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif()
MESSAGE("Build type " ${CMAKE_BUILD_TYPE})

##########################################################################
# Local CMake modules                                                    #
##########################################################################
list(APPEND CMAKE_MODULE_PATH ${alps_SOURCE_DIR}/config/CMake)
include(CheckFunctionExists)
include(CheckIncludeFile)
include(CheckTypeSize)
include(FindPythonInterp)
include(BoostUtils)
include(ALPSCore)

add_custom_target(${PROJECT_NAME}-test)
include(ALPSTesting)

#include(BoostConfig)
#include(BoostBuildSlave)
#include(BoostDocs)

##########################################################################

######################################################################
# Parallelization options
######################################################################
option(ALPS_MPI "Enable/disable MPI" ON)
option(ALPS_OMP "Enable/disable OpenMP" OFF)

######################################################################
# Standard test
######################################################################
# CHECK_FUNCTION_EXISTS(sqrt HAVE_SQRT)

######################################################################
# Overwrite with environments if they differ from the set values
######################################################################
IF(NOT $ENV{ALPS_MPI} MATCHES ALPS_MPI)
  set(ALPS_MPI $ENV{ALPS_MPI}) 
ENDIF(NOT $ENV{ALPS_MPI} MATCHES ALPS_MPI)
IF(NOT $ENV{ALPS_OMP} MATCHES ALPS_OMP)
  set(ALPS_OMP $ENV{ALPS_OMP}) 
ENDIF(NOT $ENV{ALPS_OMP} MATCHES ALPS_OMP)

######################################################################
#list of options to determine libraries. First, try to use libraries available
######################################################################
set(INTEL_COMPILER 0)
set(ENABLE_OPENMP 0)
set(ENABLE_PHDF5 0 CACHE BOOL "Enable phdf5 for output")
set(HAVE_MPI 0)
set(HAVE_OOMPI 0)
set(HAVE_LIBFFTW 0)
set(HAVE_LIBXMLPP 0)
set(HAVE_LIBHDF5 0)
set(HAVE_LIBSPRNG 0)
set(HAVE_LIBBLITZ 0)
set(INTEL_COMPILER 0)
set(MAC_VECLIB 0)
set(HAVE_MKL 0)
set(HAVE_ACML 0)

######################################################
# Directory where customize cmake files reside
######################################################
set (PROJECT_CMAKE ${alps_SOURCE_DIR}/CMake)

######################################################
# Get the current SVN version
######################################################
FIND_PACKAGE(Subversion) 
IF(Subversion_FOUND) 
  Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Project) 
  MESSAGE("-- Current revision : ${Project_WC_REVISION} modified on ${Project_WC_LAST_CHANGED_DATE}") 
  set(ALPS_BRANCH ${Project_WC_REVISION})
  set(ALPS_LAST_CHANGED_DATE ${Project_WC_LAST_CHANGED_DATE})
ENDIF(Subversion_FOUND)

######################################################
#COMPILER choose one of the cmake files to customize the compiler options
#If nothing is chosen, default settings by cmake will  be used.
#If the automatic detection does not work, comment out everything 
#upto COMPILER MANUAL SELECTION and use one of the customized cmake file.
######################################################
include(ALPSChooseCompiler)

#COMPIELR MANUAL SELECTION
#If automatic detection does not work, comment out above and select one
#of these files
#-------------------------------------------------------------------
#include(GNUCompilers)
#include(IntelCompilers)
#include(IBMCompilers)
#include(KCCCompilers)

MESSAGE("-- CXX compiler flags: ${CMAKE_CXX_FLAGS}")
MESSAGE("-- CC  compiler flags: ${CMAKE_C_FLAGS}")

set(Boost_USE_STATIC_LIBS ${BUILD_SHARED_LIBS})
set(Boost_ADDITIONAL_VERSIONS "1.37" )
include(FindBoost)
#find_package( Boost 1.35.0 COMPONENTS date_time filesystem program_options serialization REQUIRED )

######################################################
# DO NOT MODIFY BELOW
######################################################

include(FindThreads)
include(FindLibXml2)
include(FindLapack)
#include(FindHDF5)
include(FindPkgConfig)
#include(FindFFTW)

IF(ALPS_MPI)
  include(addMPI)
  IF(FOUND_MPI)
    LINK_LIBRARIES(${MPI_LIBRARIES})
  ENDIF(FOUND_MPI)
ENDIF(ALPS_MPI)

IF(LIBXML2_FOUND)
  INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})
  LINK_LIBRARIES(${LIBXML2_LIBRARIES})
  set(HAVE_LIBXML2 1)
ENDIF(LIBXML2_FOUND)

IF(Boost_FOUND)
  set(HAVE_LIBBOOST 1)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
  LINK_LIBRARIES(${Boost_LIBRARIES})
#  LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
ELSE(Boost_FOUND)
  MESSAGE(FATAL_ERROR "Require boost 1.35.x or higher.")
ENDIF(Boost_FOUND)

IF(FOUND_FFTW)
  set(HAVE_LIBFFTW 1)
  INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIR})
  LINK_LIBRARIES(${FFTW_LIBRARY})
ENDIF(FOUND_FFTW)

IF(FOUND_HDF5)
  set(HAVE_LIBHDF5 1)
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  LINK_LIBRARIES(${HDF5_LIBRARY})
  #zlib
  include(FindZLIB)
  IF(ZLIB_FOUND)
    LINK_LIBRARIES(${ZLIB_LIBRARY})
  ENDIF(ZLIB_FOUND)
  #szlib
#  include(FindSZLIB)
#  IF(SZLIB_FOUND)
#    LINK_LIBRARIES(${SZLIB_LIBRARY})
#  ENDIF(SZLIB_FOUND)
ELSE(FOUND_HDF5)
#  MESSAGE(FATAL_ERROR "Require hdf5 1.6.4 or higher. set HDF5_HOME")
ENDIF(FOUND_HDF5)


##########################################################################
# Binary packages                                                        #
##########################################################################
if (ALPS_LIGHT)
  set(PACKNAME "ALPS-light")
else(ALPS_LIGHT)
  set(PACKNAME "ALPS")
endif(ALPS_LIGHT)  

set(CPACK_PACKAGE_NAME "${PACKNAME}")
set(CPACK_PACKAGE_VENDOR "alps.comp-phys.org")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PACKNAME} ${ALPS_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${alps_SOURCE_DIR}/README.txt")
set(CPACK_RESOURCE_FILE_README "${alps_SOURCE_DIR}/README.txt")
set(CPACK_RESOURCE_FILE_LICENSE "${alps_SOURCE_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_WELCOME "${alps_SOURCE_DIR}/Welcome.txt")
set(CPACK_PACKAGE_VERSION "${ALPS_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${ALPS_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${ALPS_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_SUBMINOR "${ALPS_VERSION_SUBMINOR}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "/usr/local")

if(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  # NOTE: No ALPS icon yet
#  set(CPACK_PACKAGE_ICON "${alps_SOURCE_DIR}/tools/build/CMake\\\\InstallIcon.bmp")
#  set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\MyExecutable.exe")
  set(CPACK_NSIS_DISPLAY_NAME "ALPS ${ALPS_VERSION_MAJOR}.${ALPS_VERSION_MINOR}.${ALPS_VERSION_SUBMINOR}")
  set(CPACK_NSIS_HELP_LINK "http:\\\\\\\\alps.comp-phys.org")
  set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\alps.comp-phys.org")
  set(CPACK_NSIS_CONTACT "comp-phys-alps-users@lists.comp-phys.org")
  set(CPACK_NSIS_MODIFY_PATH ON)
  
  # Encode the compiler name in the package 
  if (MSVC60)
    set(CPACK_PACKAGE_FILE_NAME "${PACKNAME}-${ALPS_VERSION}-vc6")
    set(CPACK_NSIS_DISPLAY_NAME "${CPACK_NSIS_DISPLAY_NAME} for Microsoft Visual C++ 6")
  elseif (MSVC70)
    set(CPACK_PACKAGE_FILE_NAME "${PACKNAME}-${ALPS_VERSION}-vc7")
    set(CPACK_NSIS_DISPLAY_NAME "${CPACK_NSIS_DISPLAY_NAME} for Microsoft Visual Studio 2002")
  elseif (MSVC71)
    set(CPACK_PACKAGE_FILE_NAME "${PACKNAME}-${ALPS_VERSION}-vc71")
    set(CPACK_NSIS_DISPLAY_NAME "${CPACK_NSIS_DISPLAY_NAME} for Microsoft Visual Studio 2003")
  elseif (MSVC80)
    set(CPACK_PACKAGE_FILE_NAME "${PACKNAME}-${ALPS_VERSION}-vc8")
    set(CPACK_NSIS_DISPLAY_NAME "${CPACK_NSIS_DISPLAY_NAME} for Microsoft Visual Studio 2005")    
  elseif (MSVC90)
    set(CPACK_PACKAGE_FILE_NAME "${PACKNAME}-${ALPS_VERSION}-vc9")
    set(CPACK_NSIS_DISPLAY_NAME "${CPACK_NSIS_DISPLAY_NAME} for Microsoft Visual Studio 2008")
  elseif (BORLAND)
    set(CPACK_PACKAGE_FILE_NAME "${PACKNAME}-${ALPS_VERSION}-borland")  
    set(CPACK_NSIS_DISPLAY_NAME "${CPACK_NSIS_DISPLAY_NAME} for Borland C++ Builder")    
  endif (MSVC60)
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${CPACK_NSIS_DISPLAY_NAME}")
endif(WIN32 AND NOT UNIX)
include(CPack)

if(COMMAND cpack_configure_downloads)
  cpack_configure_downloads(
  "http://alps.comp-phys.org/software/releases/"
    ALL ADD_REMOVE)
endif()
##########################################################################



###################################################
# Done searching libraries. Start building alps
###################################################

MESSAGE("-- Ready to build alps")

INCLUDE_DIRECTORIES(
  ${alps_SOURCE_DIR}/src
  ${alps_BINARY_DIR}/src
  )
  
ADD_DEFINITIONS(-DHAVE_CONFIG_H)



##########################################################################
# Building ALPS libraries                                               #
##########################################################################

add_subdirectory(src)
add_subdirectory(lib)
add_subdirectory(example)
add_subdirectory(test)
add_subdirectory(tool)
