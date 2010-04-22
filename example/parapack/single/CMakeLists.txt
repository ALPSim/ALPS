#  Copyright Synge Todo 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})

add_library(ising_single_core ising.C)
target_link_libraries(ising_single_core alps)

add_executable(ising_single main.C)
target_link_libraries(ising_single ising_single_core)

add_executable(ising_single_ns main_ns.C)
target_link_libraries(ising_single_ns ising_single_core)
add_alps_test(ising_single_ns ising_ns)
