#  Copyright Synge Todo 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})

if(ALPS_BUILD_EXAMPLES)
add_library(ising_multiple_common ising.C)
target_link_libraries(ising_multiple_common alps)

add_executable(ising_multiple main.C)
target_link_libraries(ising_multiple ising_multiple_common)

add_executable(ising_multiple_ns main_ns.C)
target_link_libraries(ising_multiple_ns ising_multiple_common)
if(MSVC)
  get_target_property(EXE_NAME ising_multiple_ns LOCATION)
  add_custom_command(TARGET ising_multiple_ns POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${EXE_NAME} ${PROJECT_BINARY_DIR}/bin)
endif(MSVC)
endif(ALPS_BUILD_EXAMPLES)
