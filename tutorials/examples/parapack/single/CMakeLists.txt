#  Copyright Synge Todo 2010.
#   SPDX-License-Identifier: MIT

add_executable(ising_single main.C ising.C)
target_link_libraries(ising_single PRIVATE ALPS::alps)

alps_add_test(
  ising_single
  INPUT ising
  OUTPUT ising)
