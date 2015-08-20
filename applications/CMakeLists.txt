#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

######################################################################
# applications
######################################################################

enable_testing()

if (ALPS_BUILD_SOURCE)
  install(FILES CMakeLists.txt DESTINATION alps/applications)
  install(DIRECTORY mc diag dmft dmrg DESTINATION alps/applications
          FILES_MATCHING PATTERN "*" PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
  add_subdirectory(qmc)
else (ALPS_BUILD_SOURCE)
  add_definitions(${LAPACK_DEFINITIONS})
  add_subdirectory(mc)
  add_subdirectory(qmc)
  add_subdirectory(diag)
  add_subdirectory(dmrg)
  add_subdirectory(dmft)
endif (ALPS_BUILD_SOURCE)


######################################################################
# GUI's
######################################################################

# add_subdirectory(gui/simulator)

