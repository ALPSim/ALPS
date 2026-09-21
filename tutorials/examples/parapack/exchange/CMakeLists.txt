#  Copyright Synge Todo 2010.
#   SPDX-License-Identifier: MIT

add_executable(exchange main.C ising.C loop.C)
target_link_libraries(exchange PRIVATE ALPS::alps)

alps_add_test(exchange)
