# DMRG scratch-file tests, owned by the solver.

add_executable(dmrg_temporary_files dmrg_temporary_files.cpp)
target_include_directories(dmrg_temporary_files PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(dmrg_temporary_files PRIVATE ALPS::utilities)
alps_add_gtest(dmrg_temporary_files COMPONENT dmrg LABELS utilities compatibility)

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
