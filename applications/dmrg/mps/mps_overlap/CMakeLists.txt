set(MPS_OPTIM_SOURCES main.cpp run_sim.cpp)

add_executable(mps_overlap main.cpp)
target_link_libraries(mps_overlap ${MYAPP_LIBRARIES})

if(ALPS_PYTHON_WHEEL)
set_target_properties(mps_overlap  PROPERTIES INSTALL_RPATH "$ORIGIN/../lib;" )
install(TARGETS mps_overlap EXPORT ALPSMPS-targets RUNTIME DESTINATION pyalps/bin COMPONENT applications)
else()
install(TARGETS mps_overlap EXPORT ALPSMPS-targets RUNTIME DESTINATION bin COMPONENT applications)
endif()
