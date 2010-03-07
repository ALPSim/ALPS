#  Copyright Synge Todo and Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# This module looks for fop and will define 
# FOP_FOUND, FOP and FOP_EXECUTABLE 


FIND_PROGRAM(FOP_EXECUTABLE
  NAMES fop fop.sh
  PATHS ${Boost_ROOT_DIR}/tools/boostbook
  PATH_SUFFIXES fop bin fop-0.94
)

# for compatibility
SET(FOP ${FOP_EXECUTABLE})

# handle the QUIETLY and REQUIRED arguments and set FOP_FOUND to TRUE if 
# all listed variables are TRUE

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Fop DEFAULT_MSG FOP_EXECUTABLE)

MARK_AS_ADVANCED( FOP_EXECUTABLE )

