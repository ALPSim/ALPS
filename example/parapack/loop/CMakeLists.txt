#  Copyright Synge Todo 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})

set(SOURCE_CORE loop.C)

if(BUILD_SHARED_LIBS AND NOT MSVC)
  add_library(loop_single_core ${SOURCE_CORE})
  target_link_libraries(loop_single_core alps)
  add_executable(loop_single main.C)
  target_link_libraries(loop_single loop_single_core)
  add_executable(loop_single_ns main_ns.C)
  target_link_libraries(loop_single_ns loop_single_core)
else(BUILD_SHARED_LIBS AND NOT MSVC)
  add_executable(loop_single main.C ${SOURCE_CORE})
  target_link_libraries(loop_single alps)
  add_executable(loop_single_ns main_ns.C ${SOURCE_CORE})
  target_link_libraries(loop_single_ns alps)
endif(BUILD_SHARED_LIBS AND NOT MSVC)

add_alps_test(loop_single_ns loop_ns)
