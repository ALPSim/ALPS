#  Copyright Synge Todo 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

set(SOURCE_CORE ising.C)

if(MSVC)
  add_library(ising_multiple_core STATIC ${SOURCE_CORE})
  target_link_libraries(ising_multiple_core alps)
else(MSVC)
  if(BUILD_SHARED_LIBS)
    add_library(ising_multiple_core ${SOURCE_CORE})
    target_link_libraries(ising_multiple_core alps)
  endif(BUILD_SHARED_LIBS)
endif(MSVC)

if(MSVC OR BUILD_SHARED_LIBS)
  add_executable(ising_multiple main.C)
  target_link_libraries(ising_multiple ising_multiple_core)
  add_executable(ising_multiple_ns main_ns.C)
  target_link_libraries(ising_multiple_ns ising_multiple_core)
else(MSVC OR BUILD_SHARED_LIBS)
  add_executable(ising_multiple main.C ${SOURCE_CORE})
  target_link_libraries(ising_multiple alps)
  add_executable(ising_multiple_ns main_ns.C ${SOURCE_CORE})
  target_link_libraries(ising_multiple_ns alps)
endif(MSVC OR BUILD_SHARED_LIBS)
