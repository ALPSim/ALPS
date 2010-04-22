#  Copyright Synge Todo 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})

add_library(ising_multiple_core ising.C)
target_link_libraries(ising_multiple_core alps)

add_executable(ising_multiple main.C)
target_link_libraries(ising_multiple ising_multiple_core)

add_executable(ising_multiple_ns main_ns.C)
target_link_libraries(ising_multiple_ns ising_multiple_core)
