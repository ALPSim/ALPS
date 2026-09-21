#  Copyright Synge Todo 2010.
#   SPDX-License-Identifier: MIT

add_executable(loop_single main.C loop.C)
target_link_libraries(loop_single PRIVATE ALPS::alps)

alps_add_test(
  loop_single
  INPUT loop
  OUTPUT loop)
