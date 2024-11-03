#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
 
add_library(worm_impl STATIC WModel.C WRun.C Wcheck.C Wdostep.C Winit.C Wmeas.C)
target_link_libraries(worm_impl alps)

add_executable(worm main.C)
add_executable(worm_evaluate evaluate.C)

target_link_libraries(worm worm_impl)
target_link_libraries(worm_evaluate worm_impl)
if(ALPS_PYTHON_WHEEL)
set_target_properties(worm  PROPERTIES INSTALL_RPATH "${CMAKE_BINARY_DIR}/src/alps;${CMAKE_BINARY_DIR}/src/boost" )
set_target_properties(worm_evaluate  PROPERTIES INSTALL_RPATH "${CMAKE_BINARY_DIR}/src/alps;${CMAKE_BINARY_DIR}/src/boost" )
install(TARGETS worm RUNTIME DESTINATION pyalps/bin COMPONENT applications)
install(TARGETS worm_evaluate RUNTIME DESTINATION pyalps/bin COMPONENT applications)
else()
install(TARGETS worm RUNTIME DESTINATION bin COMPONENT applications)
install(TARGETS worm_evaluate RUNTIME DESTINATION bin COMPONENT applications)
endif()