#  Copyright Synge Todo 2010.
#   SPDX-License-Identifier: MIT

add_executable(ising_multiple main.C ising.C)
target_link_libraries(ising_multiple PRIVATE ALPS::alps)
