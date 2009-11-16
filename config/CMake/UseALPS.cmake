#
# This module is provided as ALPS_USE_FILE by ALPSConfig.cmake.  It can
# be INCLUDEd in a project to load the needed compiler and linker
# settings to use ALPS.
#

IF(NOT ALPS_USE_FILE_INCLUDED)
  SET(ALPS_USE_FILE_INCLUDED 1)

  # Add include directories needed to use ALPS.
  INCLUDE_DIRECTORIES(${ALPS_INCLUDE_DIRS})

  # Add link directories needed to use ALPS.
  LINK_DIRECTORIES(${ALPS_LIBRARY_DIRS})

ENDIF(NOT ALPS_USE_FILE_INCLUDED)
