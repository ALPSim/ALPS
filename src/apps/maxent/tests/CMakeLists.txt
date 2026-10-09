add_executable(maxent_linear_grid_numeric linear_grid.cpp)
target_link_libraries(maxent_linear_grid_numeric PRIVATE ALPS::maxent)
alps_add_test(maxent_linear_grid_numeric)
set_property(TEST maxent_linear_grid_numeric APPEND PROPERTY LABELS "maxent;scientific;regression")
set_property(TEST maxent_linear_grid_numeric APPEND PROPERTY ENVIRONMENT_MODIFICATION
  "TMPDIR=set:${CMAKE_CURRENT_BINARY_DIR}/test-work/maxent_linear_grid_numeric")
