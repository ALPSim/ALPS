#  Copyright Synge Todo 2011.
#   SPDX-License-Identifier: MIT

add_executable(ising_fortran main.C ising_impl.f90)
target_link_libraries(ising_fortran PRIVATE ALPS::alps ALPS::fortran)

alps_add_test(ising_fortran)
