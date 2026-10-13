cmake_minimum_required(VERSION 3.27...4.3)
project(ising_tutorial LANGUAGES C CXX)

find_package(ALPS CONFIG REQUIRED)

add_executable(ising solution/ising.cpp)
target_link_libraries(ising PRIVATE ALPS::alps)

# The exercise contains deliberate placeholders; build it after filling them in.
add_executable(ising-exercise EXCLUDE_FROM_ALL ising-skeleton.cpp)
target_link_libraries(ising-exercise PRIVATE ALPS::alps)
