# Copyright Matthias Troyer and Synge Todo 2009-2013.
# SPDX-License-Identifier: MIT

add_executable(maxent_linear_grid_numeric linear_grid.cpp)
target_include_directories(maxent_linear_grid_numeric PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_link_libraries(maxent_linear_grid_numeric PRIVATE ALPS::maxent)
alps_add_test(maxent_linear_grid_numeric)
set_property(TEST maxent_linear_grid_numeric PROPERTY LABELS maxent)

add_executable(maxent_reference reference.cpp)
target_link_libraries(maxent_reference PRIVATE ALPS::maxent)
add_test(NAME maxent_reference
  COMMAND maxent_reference "${CMAKE_CURRENT_SOURCE_DIR}/reference_cases.json")
set_tests_properties(maxent_reference PROPERTIES TIMEOUT 120 LABELS maxent)

find_package(Python 3.10 REQUIRED COMPONENTS Interpreter)
add_test(NAME maxent_cli_contract
  COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/cli_contract.py"
    "$<TARGET_FILE:maxent>")
set_tests_properties(maxent_cli_contract PROPERTIES TIMEOUT 60 LABELS "maxent;cli")
if(WIN32)
  foreach(directory IN LISTS ALPS_RUNTIME_LIBRARY_DIRS)
    set_property(TEST maxent_cli_contract maxent_reference APPEND PROPERTY ENVIRONMENT_MODIFICATION
      "PATH=path_list_prepend:${directory}")
  endforeach()
endif()
