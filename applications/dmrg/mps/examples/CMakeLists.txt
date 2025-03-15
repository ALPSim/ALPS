
add_executable(mps_apply_op mps_apply_op.cpp)
target_link_libraries(mps_apply_op ${MYAPP_LIBRARIES})
if(ALPS_PYTHON_WHEEL)
set_target_properties(mps_apply_op  PROPERTIES INSTALL_RPATH "$ORIGIN/../lib;" )
install(TARGETS mps_apply_op EXPORT ALPSMPS-targets RUNTIME DESTINATION pyalps/bin COMPONENT applications)
else()
install(TARGETS mps_apply_op EXPORT ALPSMPS-targets RUNTIME DESTINATION bin COMPONENT applications)
endif()
