# DMRG FileList tests, included from CMakeLists.txt.

add_executable(dmrg_temporary_files dmrg_temporary_files.cpp)
target_include_directories(dmrg_temporary_files PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(dmrg_temporary_files alps)
add_test(NAME dmrg_temporary_files COMMAND dmrg_temporary_files)
set_tests_properties(dmrg_temporary_files PROPERTIES LABELS "utility;dmrg")
