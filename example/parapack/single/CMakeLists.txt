#  Copyright Synge Todo 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})

if(ALPS_BUILD_EXAMPLES)
add_library(ising_single_core ising.C)
target_link_libraries(ising_single_core alps)

add_executable(ising_single main.C)
target_link_libraries(ising_single ising_single_core)

add_executable(ising_single_ns main_ns.C)
target_link_libraries(ising_single_ns ising_single_core)
if(MSVC)
  get_target_property(EXE_NAME ising_single_ns LOCATION)
  add_custom_command(TARGET ising_single_ns POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${EXE_NAME} ${PROJECT_BINARY_DIR}/bin)
endif(MSVC)

add_test(
  ising_single_ns 
  ${CMAKE_COMMAND}
    -Dcmd=ising_single_ns
    -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
    -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
    -Ddllexedir=${PROJECT_BINARY_DIR}/bin
    -Dinput=ising_ns.ip
    -Doutput=ising_ns.op
    -Dtimeout=600
    -P ${PROJECT_SOURCE_DIR}/config/run_test.cmake
)

endif(ALPS_BUILD_EXAMPLES)
