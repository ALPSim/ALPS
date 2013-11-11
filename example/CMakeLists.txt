#  Copyright Matthias Troyer and Synge Todo 2010-2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)


if (ALPS_BUILD_SOURCE)

  INSTALL(FILES CMakeLists.txt DESTINATION alps/example)  
  install(DIRECTORY alea parapack scheduler sampling model ietl fortran DESTINATION alps/example
          FILES_MATCHING PATTERN "*"  PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)

else (ALPS_BUILD_SOURCE)

  enable_testing()

  add_subdirectory(alea)
  add_subdirectory(parapack)
  add_subdirectory(scheduler)
  add_subdirectory(ngs)
  add_subdirectory(sampling)
  add_subdirectory(model)
  add_subdirectory(ietl)
  add_subdirectory(hdf5)
  if (ALPS_BUILD_FORTRAN)
    add_subdirectory(fortran)
  endif (ALPS_BUILD_FORTRAN)

endif (ALPS_BUILD_SOURCE)
