#  Copyright Synge Todo 2012.
#   SPDX-License-Identifier: MIT

add_executable(heisenberg heisenberg.C)
target_link_libraries(heisenberg PRIVATE ALPS::alps)
alps_add_test(heisenberg)
