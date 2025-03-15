configure_symm_file("simulation_symm/simulation_tpl.cpp.in"
                    "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_{SYMM}.cpp" CMAKE_SYMM_GROUP)
set(MPS_EIGENSTATE_RUN_SOURCES)
get_symm_files(APPEND MPS_EIGENSTATE_RUN_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_{SYMM}.cpp")
# TODO: change to add_library( OBJECT ) once CMake 2.8.8 is the minimum required
add_library(mps_eigenstate_run STATIC ${MPS_EIGENSTATE_RUN_SOURCES})
target_link_libraries(mps_eigenstate_run ${MYAPP_LIBRARIES})

add_executable(mps_optim mps_optim.cpp run_optim.cpp)
target_link_libraries(mps_optim mps_eigenstate_run ${MYAPP_LIBRARIES})

add_executable(mps_meas mps_meas.cpp run_measure.cpp)
target_link_libraries(mps_meas mps_eigenstate_run ${MYAPP_LIBRARIES})

if(ALPS_PYTHON_WHEEL)
  if(APPLE)
    set_target_properties(mps_optim PROPERTIES INSTALL_RPATH "@loader_path/../lib" )
    set_target_properties(mps_meas PROPERTIES INSTALL_RPATH "@loader_path/../lib" )
  else(APPLE)
    set_target_properties(mps_optim PROPERTIES INSTALL_RPATH "$ORIGIN/../lib" )
    set_target_properties(mps_meas PROPERTIES INSTALL_RPATH "$ORIGIN/../lib" )
  endif(APPLE)
  install(TARGETS mps_optim mps_meas EXPORT ALPSMPS-targets RUNTIME DESTINATION pyalps/bin COMPONENT applications)
else()
install(TARGETS mps_optim mps_meas EXPORT ALPSMPS-targets RUNTIME DESTINATION bin COMPONENT applications)
endif()