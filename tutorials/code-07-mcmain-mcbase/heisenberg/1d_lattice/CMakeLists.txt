cmake_minimum_required(VERSION 4.3)
project(heisenberg LANGUAGES C CXX)

find_package(ALPS CONFIG REQUIRED)

# rule for generating single core heisenberg program
add_executable(single heisenberg.cpp single.cpp)
target_link_libraries(single PRIVATE ALPS::alps)
