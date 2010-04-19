#  Copyright Synge Todo 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})

if(ALPS_BUILD_EXAMPLES)
add_library(exchange_common ising.C loop.C)
target_link_libraries(exchange_common alps)

add_executable(exchange main.C)
target_link_libraries(exchange exchange_common)

add_executable(exchange_ns main_ns.C)
target_link_libraries(exchange_ns exchange_common)
if(MSVC)
  get_target_property(EXE_NAME exchange_ns LOCATION)
  add_custom_command(TARGET exchange_ns POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${EXE_NAME} ${PROJECT_BINARY_DIR}/bin)
endif(MSVC)

add_test(
  exchange_ns 
  ${CMAKE_COMMAND}
    -Dcmd=exchange_ns
    -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
    -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
    -Ddllexedir=${PROJECT_BINARY_DIR}/bin
    -Dinput=exchange_ns.ip
    -Doutput=exchange_ns.op
    -Dtimeout=600
    -P ${PROJECT_SOURCE_DIR}/config/run_test.cmake
)

endif(ALPS_BUILD_EXAMPLES)
