configure_symm_file("simulation_symm/simulation_tpl.cpp.in"
                    "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_{SYMM}.cpp" CMAKE_SYMM_GROUP)
set(MPS_EVOLVE_SOURCES main.cpp run_sim.cpp)
get_symm_files(APPEND MPS_EVOLVE_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_{SYMM}.cpp")

add_executable(mps_evolve ${MPS_EVOLVE_SOURCES})
target_link_libraries(mps_evolve ${MYAPP_LIBRARIES})
if(ALPS_PYTHON_WHEEL)
set_target_properties(mps_evolve  PROPERTIES INSTALL_RPATH "$ORIGIN/../lib;" )
install(TARGETS mps_evolve EXPORT ALPSMPS-targets RUNTIME DESTINATION pyalps/bin COMPONENT applications)
else()
install(TARGETS mps_evolve EXPORT ALPSMPS-targets RUNTIME DESTINATION bin COMPONENT applications)
endif()