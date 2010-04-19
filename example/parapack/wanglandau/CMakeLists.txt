#  Copyright Synge Todo 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})

if(ALPS_BUILD_EXAMPLES)
add_library(wanglandau_common wanglandau.C)
target_link_libraries(wanglandau_common alps)

add_executable(wanglandau main.C)
target_link_libraries(wanglandau wanglandau_common)

add_executable(wanglandau_ns main_ns.C)
target_link_libraries(wanglandau_ns wanglandau_common)
if(MSVC)
  get_target_property(EXE_NAME wanglandau_ns LOCATION)
  add_custom_command(TARGET wanglandau_ns POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${EXE_NAME} ${PROJECT_BINARY_DIR}/bin)
endif(MSVC)

add_test(
  wanglandau_ns 
  ${CMAKE_COMMAND}
    -Dcmd=wanglandau_ns
    -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
    -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
    -Ddllexedir=${PROJECT_BINARY_DIR}/bin
    -Dinput=wanglandau_ns.ip
    -Doutput=wanglandau_ns.op
    -Dtimeout=600
    -P ${PROJECT_SOURCE_DIR}/config/run_test.cmake
)

endif(ALPS_BUILD_EXAMPLES)
