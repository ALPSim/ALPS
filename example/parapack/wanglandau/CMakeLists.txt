#  Copyright Synge Todo 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})

set(SOURCE_CORE wanglandau.C)

if(BUILD_SHARED_LIBS AND NOT MSVC)
  add_library(wanglandau_core ${SOURCE_CORE})
  target_link_libraries(wanglandau_core alps)
  add_executable(wanglandau main.C)
  target_link_libraries(wanglandau wanglandau_core)
  add_executable(wanglandau_ns main_ns.C)
  target_link_libraries(wanglandau_ns wanglandau_core)
else(BUILD_SHARED_LIBS AND NOT MSVC)
  add_executable(wanglandau main.C ${SOURCE_CORE})
  target_link_libraries(wanglandau alps)
  add_executable(wanglandau_ns main_ns.C ${SOURCE_CORE})
  target_link_libraries(wanglandau_ns alps)
endif(BUILD_SHARED_LIBS AND NOT MSVC)

add_alps_test(wanglandau_ns)
