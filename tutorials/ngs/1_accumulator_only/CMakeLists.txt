cmake_minimum_required(VERSION 3.27)
project(alps_ising_tutorial LANGUAGES C CXX)

find_package(ALPS CONFIG REQUIRED)
add_executable(ising ising.cpp main.cpp)
target_link_libraries(ising PRIVATE ALPS::alps)
