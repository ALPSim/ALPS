#  Copyright Synge Todo 2010.
#   SPDX-License-Identifier: MIT

add_executable(wanglandau main.C wanglandau.C)
target_link_libraries(wanglandau PRIVATE ALPS::alps)

alps_add_test(wanglandau)
