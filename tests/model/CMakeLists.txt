#  Copyright Matthias Troyer 2009 - 2010.
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
  example12
  example13
  example14
  example15
  example16
  example17
  example18)
  add_executable(model_${name} ${name}.C)

  target_link_libraries(model_${name} PRIVATE ALPS::alps)
  alps_add_test(
    model_${name}
    INPUT ${name}
    OUTPUT ${name})
  set_property(TEST model_${name} PROPERTY LABELS model)
endforeach()
