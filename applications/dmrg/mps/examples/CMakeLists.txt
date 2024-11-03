
add_executable(mps_apply_op mps_apply_op.cpp)
target_link_libraries(mps_apply_op ${MYAPP_LIBRARIES})
if(ALPS_PYTHON_WHEEL)
set_target_properties(mps_apply_op  PROPERTIES INSTALL_RPATH "${CMAKE_BINARY_DIR}/src/alps;${CMAKE_BINARY_DIR}/src/boost;${CMAKE_BINARY_DIR}/applications/dmrg/mps" )
install(TARGETS mps_apply_op EXPORT ALPSMPS-targets RUNTIME DESTINATION pyalps/bin COMPONENT applications)
else()
install(TARGETS mps_apply_op EXPORT ALPSMPS-targets RUNTIME DESTINATION bin COMPONENT applications)
endif()
