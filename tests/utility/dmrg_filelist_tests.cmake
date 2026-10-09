# DMRG scratch-file tests, included from CMakeLists.txt.

add_executable(dmrg_temporary_files dmrg_temporary_files.cpp)
target_include_directories(dmrg_temporary_files PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(dmrg_temporary_files alps)
add_test(NAME dmrg_temporary_files COMMAND dmrg_temporary_files)
set_tests_properties(dmrg_temporary_files PROPERTIES LABELS "utility;dmrg")

if(TARGET dmrg AND TARGET parameter2xml)
  add_test(NAME dmrg_scratch_cleanup
    COMMAND ${CMAKE_COMMAND}
      "-Ddmrg=$<TARGET_FILE:dmrg>"
      "-Dparameter2xml=$<TARGET_FILE:parameter2xml>"
      "-Dsource_dir=${PROJECT_SOURCE_DIR}"
      "-Dxml_dir=${PROJECT_BINARY_DIR}/src/alps/resources"
      "-Dtest_dir=${CMAKE_CURRENT_BINARY_DIR}/dmrg-scratch"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/dmrg_scratch_cleanup.cmake")
  set_tests_properties(dmrg_scratch_cleanup PROPERTIES LABELS "dmrg" TIMEOUT 120)
endif()
