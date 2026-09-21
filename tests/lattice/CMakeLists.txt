#  Copyright Matthias Troyer and Synge Todo 2009 - 2010.
#   SPDX-License-Identifier: MIT

foreach(
  name
  example1
  example2
  example3
  example4
  example5
  example6
  example7
  example8
  example9
  example10
  example11
  label)
  add_executable(lattice_${name} ${name}.C)

  target_link_libraries(lattice_${name} PRIVATE ALPS::alps)
  alps_add_test(
    lattice_${name}
    INPUT ${name}
    OUTPUT ${name})
  set_property(TEST lattice_${name} PROPERTY LABELS lattice)
endforeach()
foreach(name coloring parity)
  add_executable(${name} ${name}.C)

  target_link_libraries(${name} PRIVATE ALPS::alps)
  alps_add_test(${name})
  set_property(TEST ${name} PROPERTY LABELS lattice)
endforeach()
