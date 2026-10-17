# DMRG scratch-file and SAVE_MPS tests of the application. The solver's own
# scratch directory and its MPS are tested with the classic DMRG component.

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

if(TARGET dmrg AND TARGET parameter2xml)
  add_executable(dmrg_mps_check dmrg_mps_check.cpp)
  target_link_libraries(dmrg_mps_check PRIVATE ALPS::alps)
  add_test(NAME dmrg_save_mps
    COMMAND ${CMAKE_COMMAND}
      "-Ddmrg=$<TARGET_FILE:dmrg>"
      "-Dparameter2xml=$<TARGET_FILE:parameter2xml>"
      "-Dchecker=$<TARGET_FILE:dmrg_mps_check>"
      "-Dxml_dir=${PROJECT_SOURCE_DIR}/src/alps/resources"
      "-Dtest_dir=${CMAKE_CURRENT_BINARY_DIR}/dmrg-save-mps"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/dmrg_save_mps.cmake")
  set_tests_properties(dmrg_save_mps PROPERTIES LABELS "dmrg" TIMEOUT 300)
endif()
