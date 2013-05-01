cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(alpsize NONE)

# enable C and C++ compilers
enable_language(C CXX)

# rule for generating wolff program
add_executable(wolff wolff.c)
target_link_libraries(wolff m)
