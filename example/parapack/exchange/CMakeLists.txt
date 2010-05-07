#  Copyright Synge Todo 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})

set(SOURCE_CORE ising.C loop.C)

if(MSVC)
  add_library(exchange_core STATIC ${SOURCE_CORE})
  target_link_libraries(exchange_core alps)
else(MSVC)
  if(BUILD_SHARED_LIBS)
    add_library(exchange_core ${SOURCE_CORE})
    target_link_libraries(exchange_core alps)
  endif(BUILD_SHARED_LIBS)
endif(MSVC)

if(MSVC OR BUILD_SHARED_LIBS)
  add_executable(exchange main.C)
  target_link_libraries(exchange exchange_core)
  add_executable(exchange_ns main_ns.C)
  target_link_libraries(exchange_ns exchange_core)
else(MSVC OR BUILD_SHARED_LIBS)
  add_executable(exchange main.C ${SOURCE_CORE})
  target_link_libraries(exchange alps)
  add_executable(exchange_ns main_ns.C ${SOURCE_CORE})
  target_link_libraries(exchange_ns alps)
endif(MSVC OR BUILD_SHARED_LIBS)

add_alps_test(exchange_ns)
