#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if (ALPS_BUILD_SOURCE)
if(ALPS_PYTHON_WHEEL)

else()
  install(FILES CMakeLists.txt qmc.h DESTINATION alps/applications/qmc)
  install(DIRECTORY checksign looper qwl sse4 worms dwa DESTINATION alps/applications/qmc
          FILES_MATCHING PATTERN "*"  PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
endif()
else (ALPS_BUILD_SOURCE)

  add_subdirectory(checksign)
  add_subdirectory(looper)
  add_subdirectory(qwl)
  #add_subdirectory(sse)
  add_subdirectory(sse4)
  add_subdirectory(worms)
  add_subdirectory(dwa)

endif (ALPS_BUILD_SOURCE)
