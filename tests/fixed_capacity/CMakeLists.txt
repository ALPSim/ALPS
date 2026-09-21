#  Copyright Matthias Troyer and Synge Todo 2009 - 2010.
#   SPDX-License-Identifier: MIT

foreach(name fixed_capacity_deque fixed_capacity_traits fixed_capacity_vector test_deque
             test_vector)
  add_executable(${name} ${name}.C)

  target_link_libraries(${name} PRIVATE ALPS::alps)
  alps_add_test(${name})
  set_property(TEST ${name} PROPERTY LABELS fixed_capacity)
endforeach()
