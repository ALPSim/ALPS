#  Copyright Synge Todo 2010-2012.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

enable_testing()

add_subdirectory(single)
add_subdirectory(ising)
add_subdirectory(heisenberg)
add_subdirectory(loop)
add_subdirectory(exchange)
add_subdirectory(wanglandau)

if(${ALPS_HAVE_MPI})
  add_subdirectory(multiple)
endif(${ALPS_HAVE_MPI})
