#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if (ALPS_BUILD_NG AND ALPS_BUILD_EXAMPLES)
  if (MPI_FOUND)
    add_executable(ngs main.cpp)
    add_dependencies(ngs alps)
    target_link_libraries(ngs alps)
    
    add_executable(ngs2 ising2.cpp main2.cpp)
    add_dependencies(ngs2 alps)
    target_link_libraries(ngs2 alps)
  endif(MPI_FOUND)
endif (ALPS_BUILD_NG AND ALPS_BUILD_EXAMPLES)
