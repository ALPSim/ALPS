cmake_minimum_required(VERSION 3.27)
project(alpsize LANGUAGES C CXX Fortran)

find_package(ALPS CONFIG REQUIRED COMPONENTS fortran)
find_package(OpenMP REQUIRED COMPONENTS Fortran)

# rule for generating ising program
add_executable(ising ising.C ising_impl.f90)
target_link_libraries(ising PRIVATE ALPS::fortran OpenMP::OpenMP_Fortran)
