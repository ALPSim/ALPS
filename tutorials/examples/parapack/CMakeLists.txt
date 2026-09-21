#  Copyright Synge Todo 2010-2012.
#   SPDX-License-Identifier: MIT

add_subdirectory(single)
add_subdirectory(ising)
add_subdirectory(heisenberg)
add_subdirectory(loop)
add_subdirectory(exchange)
add_subdirectory(wanglandau)

if(ALPS_ENABLE_MPI)
  add_subdirectory(multiple)
endif()
