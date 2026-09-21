cmake_minimum_required(VERSION 4.3)
project(alpsize LANGUAGES C CXX)

# rule for generating wolff program
add_executable(wolff wolff.c)
if(UNIX)
  target_link_libraries(wolff PRIVATE m)
endif()
