#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

add_subdirectory(looper)
add_subdirectory(sse)
add_subdirectory(checksign)
add_subdirectory(worms)
add_subdirectory(qwl)

IF(EXISTS ${PROJECT_SOURCE_DIR}/applications/qmc/sse-loop)
  add_subdirectory(sse-loop)
ENDIF(EXISTS ${PROJECT_SOURCE_DIR}/applications/qmc/sse-loop)
