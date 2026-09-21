#  Copyright Synge Todo 2011.
#   SPDX-License-Identifier: MIT

add_executable(hello_fortran main.C hello_impl.f90)
target_link_libraries(hello_fortran PRIVATE ALPS::alps ALPS::fortran)

alps_add_test(hello_fortran)
