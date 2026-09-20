cmake_minimum_required(VERSION 3.27)
project(alpsize LANGUAGES C CXX)

find_package(ALPS CONFIG REQUIRED)

# rule for generating 'hello world' program
add_executable(hello hello.C)
target_link_libraries(hello PRIVATE ALPS::alps)
