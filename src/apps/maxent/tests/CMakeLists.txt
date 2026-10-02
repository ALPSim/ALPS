# Copyright Matthias Troyer and Synge Todo 2009-2013.
# SPDX-License-Identifier: MIT

add_executable(maxent_linear_grid_numeric linear_grid.cpp)
target_include_directories(maxent_linear_grid_numeric PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_link_libraries(maxent_linear_grid_numeric PRIVATE ALPS::maxent)
alps_add_test(maxent_linear_grid_numeric)
set_property(TEST maxent_linear_grid_numeric PROPERTY LABELS maxent)
