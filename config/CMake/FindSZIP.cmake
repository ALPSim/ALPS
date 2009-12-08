# - Find szip
# Find the native SZIP includes and library
#
#  SZIP_INCLUDE_DIRS - where to find szip.h, etc.
#  SZIP_LIBRARIES    - List of libraries when using szip.
#  SZIP_FOUND        - True if szip found.

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

IF (SZIP_INCLUDE_DIR)
  # Already in cache, be silent
  SET(SZIP_FIND_QUIETLY TRUE)
ENDIF (SZIP_INCLUDE_DIR)

FIND_PATH(SZIP_INCLUDE_DIR szlib.h)

SET(SZIP_NAMES szip)
FIND_LIBRARY(SZIP_LIBRARY NAMES ${SZIP_NAMES} )
MARK_AS_ADVANCED( SZIP_LIBRARY SZIP_INCLUDE_DIR )

# Per-recommendation
SET(SZIP_INCLUDE_DIRS "${SZIP_INCLUDE_DIR}")
SET(SZIP_LIBRARIES    "${SZIP_LIBRARY}")

# handle the QUIETLY and REQUIRED arguments and set SZIP_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SZIP DEFAULT_MSG SZIP_LIBRARIES SZIP_INCLUDE_DIRS)

