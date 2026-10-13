cmake_minimum_required(VERSION 3.27...4.3)
project(alpsize LANGUAGES C CXX Fortran)

find_package(ALPS CONFIG REQUIRED COMPONENTS fortran)
find_package(OpenMP REQUIRED COMPONENTS Fortran)

# rule for generating 'hello world' program
add_executable(hello hello.C hello_impl.f90)
target_link_libraries(hello PRIVATE ALPS::fortran OpenMP::OpenMP_Fortran)
