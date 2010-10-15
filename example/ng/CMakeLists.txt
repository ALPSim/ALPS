#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if (ALPS_BUILD_NG AND ALPS_BUILD_EXAMPLES)
  if (MPI_FOUND)
    add_executable(ng main.cpp)
    add_dependencies(ng alps)
    target_link_libraries(ng alps)
    
    add_executable(ngdeprecated ising2.cpp main2.cpp)
    add_dependencies(ngdeprecated alps)
    target_link_libraries(ngdeprecated alps)
  endif(MPI_FOUND)
endif (ALPS_BUILD_NG AND ALPS_BUILD_EXAMPLES)
