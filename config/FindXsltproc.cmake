#  Copyright Synge Todo and Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# This module looks for xsltproc and will define 
# XSLTPROC_FOUND, XSLTPROC and XSLTPROC_EXECUTABLE 


FIND_PROGRAM(XSLTPROC_EXECUTABLE
  NAMES 
  xsltproc
)

# for compatibility
SET(XSLTPROC ${XSLTPROC_EXECUTABLE})

# handle the QUIETLY and REQUIRED arguments and set XSLTPROC_FOUND to TRUE if 
# all listed variables are TRUE

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Xsltproc DEFAULT_MSG XSLTPROC_EXECUTABLE)

MARK_AS_ADVANCED( XSLTPROC_EXECUTABLE )

