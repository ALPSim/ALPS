# Copyright Matthias Troyer and Synge Todo 2009-2013.
# SPDX-License-Identifier: MIT

add_executable(maxent_linear_grid_numeric linear_grid.cpp)
target_include_directories(maxent_linear_grid_numeric PRIVATE
  "$<TARGET_PROPERTY:alps_maxent,SOURCE_DIR>")
target_link_libraries(maxent_linear_grid_numeric PRIVATE ALPS::maxent)
alps_add_test(maxent_linear_grid_numeric)
