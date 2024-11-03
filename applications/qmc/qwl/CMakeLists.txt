#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

add_executable(qwl qwl.C)
add_executable(qwl_evaluate qwl_evaluate.C)

target_link_libraries(qwl alps)
target_link_libraries(qwl_evaluate alps)
if(ALPS_PYTHON_WHEEL)
set_target_properties(qwl  PROPERTIES INSTALL_RPATH "${CMAKE_BINARY_DIR}/src/alps;${CMAKE_BINARY_DIR}/src/boost" )
set_target_properties(qwl_evaluate  PROPERTIES INSTALL_RPATH "${CMAKE_BINARY_DIR}/src/alps;${CMAKE_BINARY_DIR}/src/boost" )
install(TARGETS qwl qwl_evaluate RUNTIME DESTINATION pyalps/bin COMPONENT applications)
else()
install(TARGETS qwl qwl_evaluate RUNTIME DESTINATION bin COMPONENT applications)
endif()
