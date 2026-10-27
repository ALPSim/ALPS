# DMRG scratch-file and parameter tests, owned by the solver.

add_executable(dmrg_temporary_files dmrg_temporary_files.cpp)
target_include_directories(dmrg_temporary_files PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(dmrg_temporary_files PRIVATE ALPS::utilities)
alps_add_gtest(dmrg_temporary_files COMPONENT dmrg LABELS utilities compatibility)

add_executable(dmrg_signals dmrg_signals.cpp)
target_include_directories(dmrg_signals PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(dmrg_signals PRIVATE ALPS::alps)
alps_add_gtest(dmrg_signals COMPONENT dmrg)

if(TARGET dmrg AND TARGET parameter2xml)
  add_test(NAME dmrg_scratch_cleanup
    COMMAND ${CMAKE_COMMAND}
      "-Ddmrg=$<TARGET_FILE:dmrg>"
      "-Dparameter2xml=$<TARGET_FILE:parameter2xml>"
      "-Dxml_dir=${PROJECT_SOURCE_DIR}/src/alps/resources"
      "-Dtest_dir=${CMAKE_CURRENT_BINARY_DIR}/dmrg-scratch"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/dmrg_scratch_cleanup.cmake")
  set_tests_properties(dmrg_scratch_cleanup PROPERTIES LABELS "dmrg" TIMEOUT 120)
endif()

add_executable(dmrg_parameters dmrg_parameters.cpp)
target_include_directories(dmrg_parameters PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(dmrg_parameters PRIVATE ALPS::alps)
alps_add_gtest(dmrg_parameters COMPONENT dmrg)

if(TARGET dmrg)
  add_test(NAME dmrg.invalid_parameter
    COMMAND "${CMAKE_COMMAND}" "-DAPP=$<TARGET_FILE:dmrg>"
      "-DXML_PATH=${PROJECT_SOURCE_DIR}/src/alps/resources"
      "-DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/invalid_parameter"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/invalid_parameter.cmake")
  set_tests_properties(dmrg.invalid_parameter PROPERTIES LABELS "dmrg;cli"
    TIMEOUT 60 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

  add_test(NAME dmrg.same_seed
    COMMAND "${CMAKE_COMMAND}" "-DAPP=$<TARGET_FILE:dmrg>"
      "-DXML_PATH=${PROJECT_SOURCE_DIR}/src/alps/resources"
      "-DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/same_seed"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/same_seed.cmake")
  set_tests_properties(dmrg.same_seed PROPERTIES LABELS "dmrg"
    TIMEOUT 120 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()
