#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if (ALPS_BUILD_EXAMPLES)
  add_executable(p2h5 p2h5.cpp)
  add_dependencies(p2h5 alps)
  target_link_libraries(p2h5 alps)
  if (MPI_FOUND)
    add_executable(ngs ising.cpp ngs.cpp)
    add_dependencies(ngs alps)
    target_link_libraries(ngs alps)
  endif(MPI_FOUND)
endif (ALPS_BUILD_EXAMPLES)
