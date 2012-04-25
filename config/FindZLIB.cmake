# - Find zlib
# Find the native ZLIB includes and library
#
#  ZLIB_INCLUDE_DIRS - where to find zlib.h, etc.
#  ZLIB_LIBRARIES    - List of libraries when using zlib.
#  ZLIB_DLLS         - List of DLLs when using zlib.
#  ZLIB_FOUND        - True if zlib found.

#=============================================================================
# Copyright 2001-2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)


IF(ALPS_PACKAGE_LIBRARIES)
  IF (UNIX AND NOT WIN32)
    MESSAGE(STATUS "Using ALPS-installed ZLIB")
    FIND_PATH(ZLIB_INCLUDE_DIR zlib.h.h /usr/include ${CMAKE_INSTALL_PREFIX}/include NO_DEFAULT_PATH)
    FIND_LIBRARY(ZLIB_LIBRARY z /usr/lib ${CMAKE_INSTALL_PREFIX}/lib NO_DEFAULT_PATH)
  ENDIF (UNIX AND NOT WIN32)
ENDIF(ALPS_PACKAGE_LIBRARIES)

IF (ZLIB_INCLUDE_DIR)
  # Already in cache, be silent
  SET(ZLIB_FIND_QUIETLY TRUE)
ENDIF (ZLIB_INCLUDE_DIR)

FIND_PATH(ZLIB_INCLUDE_DIR zlib.h PATHS "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/opt/include")

SET(ZLIB_NAMES z zlib)
FIND_LIBRARY(ZLIB_LIBRARY NAMES ${ZLIB_NAMES} PATHS "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/opt/lib")

FIND_FILE(ZLIB_DLL zlib.dll "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/opt/bin")

MARK_AS_ADVANCED( ZLIB_LIBRARY ZLIB_INCLUDE_DIR ZLIB_DLL )

# Per-recommendation
SET(ZLIB_INCLUDE_DIRS "${ZLIB_INCLUDE_DIR}")
SET(ZLIB_LIBRARIES    "${ZLIB_LIBRARY}")
SET(ZLIB_DLLS         "${ZLIB_DLL}")

# handle the QUIETLY and REQUIRED arguments and set ZLIB_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ZLIB DEFAULT_MSG ZLIB_LIBRARIES ZLIB_INCLUDE_DIRS)


IF(ALPS_USE_VISTRAILS AND WIN32 AND NOT UNIX)
  MESSAGE(STATUS "Using VisTrails zlib")
   SET(ZLIB_DLLS ${VISTRAILS_APP_DIR}/${VISTRAILS_APP_NAME}/zlib1.dll ${VISTRAILS_APP_DIR}/${VISTRAILS_APP_NAME}/zlib.dll)  
   SET(ZLIB_FOUND TRUE)  
ENDIF(ALPS_USE_VISTRAILS AND WIN32 AND NOT UNIX)
