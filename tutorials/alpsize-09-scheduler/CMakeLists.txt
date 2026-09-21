cmake_minimum_required(VERSION 4.3)
project(alpsize LANGUAGES C CXX)

find_package(ALPS CONFIG REQUIRED)

# rule for generating 'hello world' program
add_executable(hello hello.C hello_worker.C)
target_link_libraries(hello PRIVATE ALPS::alps)

# rule for generating wolff program
add_executable(wolff wolff.C wolff_worker.C)
target_link_libraries(wolff PRIVATE ALPS::alps)
