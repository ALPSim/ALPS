#  Copyright Synge Todo 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})

add_library(loop_single_core loop.C)
target_link_libraries(loop_single_core alps)

add_executable(loop_single main.C)
target_link_libraries(loop_single loop_single_core)

add_executable(loop_single_ns main_ns.C)
target_link_libraries(loop_single_ns loop_single_core)
add_alps_test(loop_single_ns loop_ns)
