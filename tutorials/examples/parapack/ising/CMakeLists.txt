#  Copyright Synge Todo 2010-2012.
#   SPDX-License-Identifier: MIT

add_executable(ising ising.C)
target_link_libraries(ising PRIVATE ALPS::alps)
alps_add_test(
  ising
  INPUT ising
  OUTPUT ising)
