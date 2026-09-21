#  Copyright Synge Todo 2010 - 2015.
#   SPDX-License-Identifier: MIT

add_executable(simplemc main.C evaluator.C ising.C xy.C heisenberg.C)
target_link_libraries(simplemc PRIVATE ALPS::alps)
alps_install_applications(simplemc)

alps_add_test(
  simplemc_ising
  TARGET simplemc
  INPUT ising
  OUTPUT ising)
alps_add_test(
  simplemc_xy
  TARGET simplemc
  INPUT xy
  OUTPUT xy)
alps_add_test(
  simplemc_heisenberg
  TARGET simplemc
  INPUT heisenberg
  OUTPUT heisenberg)
