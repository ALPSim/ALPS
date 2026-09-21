#  Copyright Matthias Troyer and Synge Todo 2009 - 2012.
#   SPDX-License-Identifier: MIT

foreach(name random_choice uniform_on_sphere_n)
  add_executable(${name} ${name}.C)

  target_link_libraries(${name} PRIVATE ALPS::alps)
  alps_add_test(${name})
  set_property(TEST ${name} PROPERTY LABELS random)
endforeach()
