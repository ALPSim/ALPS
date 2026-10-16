# DMRG scratch-file test of the application. The solver's own scratch
# directory is tested with the classic DMRG component.

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
