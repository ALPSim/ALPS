#  Copyright Synge Todo 2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})

add_executable(ising_fortran main.C ising_impl.f90)
target_link_libraries(ising_fortran alps alps_fortran)

add_alps_test(ising_fortran)
