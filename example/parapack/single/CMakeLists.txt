#  Copyright Synge Todo 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})

set(SOURCE_CORE ising.C)

if(MSVC)
  add_library(ising_single_core STATIC ${SOURCE_CORE})
  target_link_libraries(ising_single_core alps)
else(MSVC)
  if(BUILD_SHARED_LIBS)
    add_library(ising_single_core ${SOURCE_CORE})
    target_link_libraries(ising_single_core alps)
  endif(BUILD_SHARED_LIBS)
endif(MSVC)

if(MSVC OR BUILD_SHARED_LIBS)
  add_executable(ising_single main.C)
  target_link_libraries(ising_single ising_single_core)
  add_executable(ising_single_ns main_ns.C)
  target_link_libraries(ising_single_ns ising_single_core)
else(MSVC OR BUILD_SHARED_LIBS)
  add_executable(ising_single main.C ${SOURCE_CORE})
  target_link_libraries(ising_single alps)
  add_executable(ising_single_ns main_ns.C ${SOURCE_CORE})
  target_link_libraries(ising_single_ns alps)
endif(MSVC OR BUILD_SHARED_LIBS)

add_alps_test(ising_single_ns ising_ns)
