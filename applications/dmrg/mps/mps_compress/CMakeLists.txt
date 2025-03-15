
add_executable(mps_compress main.cpp)
target_link_libraries(mps_compress ${MYAPP_LIBRARIES})

if(ALPS_PYTHON_WHEEL)
set_target_properties(mps_compress  PROPERTIES INSTALL_RPATH "$ORIGIN/../lib;" )
install(TARGETS mps_compress EXPORT ALPSMPS-targets RUNTIME DESTINATION pyalps/bin COMPONENT applications)
else()
install(TARGETS mps_compress EXPORT ALPSMPS-targets RUNTIME DESTINATION bin COMPONENT applications)
endif()
