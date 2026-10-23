# SPDX-License-Identifier: MIT
add_executable(sparsediag_parameters sparsediag_parameters.cpp)
target_link_libraries(sparsediag_parameters PRIVATE ALPS::alps)
alps_add_gtest(sparsediag_parameters COMPONENT diag)

add_executable(fulldiag_parameters fulldiag_parameters.cpp)
target_link_libraries(fulldiag_parameters PRIVATE fulldiag_impl)
alps_add_gtest(fulldiag_parameters COMPONENT diag)

foreach(app sparsediag fulldiag)
  add_test(NAME diag.invalid_parameter_${app}
    COMMAND "${CMAKE_COMMAND}" "-DAPP=$<TARGET_FILE:${app}>"
      "-DXML_PATH=${PROJECT_SOURCE_DIR}/src/alps/resources"
      "-DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/invalid_parameter_${app}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/invalid_parameter.cmake")
  set_tests_properties(diag.invalid_parameter_${app} PROPERTIES LABELS "diag;cli"
    TIMEOUT 60 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()
