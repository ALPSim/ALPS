cmake_minimum_required(VERSION 3.27...4.3)
project(alpsize LANGUAGES C CXX)

find_package(ALPS CONFIG REQUIRED)

# rule for generating wolff program
add_executable(wolff wolff.C)
target_link_libraries(wolff PRIVATE ALPS::alps)
