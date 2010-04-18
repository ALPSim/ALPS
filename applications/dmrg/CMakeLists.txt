#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if (ALPS_BUILD_SOURCE)
  install(FILES CMakeLists.txt DESTINATION applications/dmrg)
  install(DIRECTORY dmrg particle_in_box DESTINATION applications/dmrg
          FILES_MATCHING PATTERN "*"  PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
else (ALPS_BUILD_SOURCE)
  add_subdirectory(dmrg)
  add_subdirectory(particle_in_box)
endif (ALPS_BUILD_SOURCE)