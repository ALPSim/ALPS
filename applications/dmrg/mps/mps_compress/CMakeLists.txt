
add_executable(mps_compress main.cpp)
target_link_libraries(mps_compress ${MYAPP_LIBRARIES})

install(TARGETS mps_compress EXPORT ALPSMPS-targets RUNTIME DESTINATION bin COMPONENT applications)
