# Copyright Matthias Troyer and Synge Todo 2009-2013.
# SPDX-License-Identifier: MIT

add_executable(maxent_linear_grid_numeric linear_grid.cpp)
target_include_directories(maxent_linear_grid_numeric PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_link_libraries(maxent_linear_grid_numeric PRIVATE ALPS::maxent)
alps_add_gtest(maxent_linear_grid_numeric COMPONENT maxent LABELS scientific TIMEOUT 120)

# Keep a native public-API/linkage probe alongside the shared scientific oracle.
add_executable(maxent_api_contract "${PROJECT_SOURCE_DIR}/tests/cmake/consumer/maxent_independent.cpp")
target_link_libraries(maxent_api_contract PRIVATE ALPS::maxent)
alps_add_test(maxent_api_contract)
set_property(TEST maxent_api_contract PROPERTY LABELS "maxent;integration")

find_package(Python 3.11 QUIET COMPONENTS Interpreter)
if(Python_Interpreter_FOUND)
  add_test(NAME maxent_cli_contract
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/cli_contract.py"
      "$<TARGET_FILE:maxent>")
  set_tests_properties(maxent_cli_contract PROPERTIES TIMEOUT 60 LABELS "maxent;cli")
  set(maxent_python_tests maxent_cli_contract)
  execute_process(COMMAND "${Python_EXECUTABLE}" -c "import numpy, h5py"
    RESULT_VARIABLE missing_reference_dependencies ERROR_QUIET)
  if(NOT missing_reference_dependencies)
    add_test(NAME maxent_reference
      COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/reference.py" "$<TARGET_FILE:maxent>")
    set_tests_properties(maxent_reference PROPERTIES TIMEOUT 120 LABELS "maxent;scientific"
      ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    list(APPEND maxent_python_tests maxent_reference)
  else()
    message(STATUS "Skipping MaxEnt reference test: NumPy and h5py are unavailable for ${Python_EXECUTABLE}")
  endif()
  if(WIN32)
    foreach(directory IN LISTS ALPS_RUNTIME_LIBRARY_DIRS)
      set_property(TEST ${maxent_python_tests} APPEND PROPERTY ENVIRONMENT_MODIFICATION
        "PATH=path_list_prepend:${directory}")
    endforeach()
  endif()
else()
  message(STATUS "Skipping MaxEnt Python tests: Python 3.11 or newer is unavailable")
endif()
