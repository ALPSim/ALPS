
add_executable(mps_apply_op mps_apply_op.cpp)
target_link_libraries(mps_apply_op ${MYAPP_LIBRARIES})

install(TARGETS mps_apply_op RUNTIME DESTINATION bin COMPONENT applications)
