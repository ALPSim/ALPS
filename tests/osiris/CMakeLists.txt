#  Copyright Matthias Troyer and Synge Todo 2009 - 2010.
#   SPDX-License-Identifier: MIT

configure_file(xdrdump2.input.cmake ${CMAKE_CURRENT_BINARY_DIR}/boostdump2.input)
configure_file(xdrdump2.input.cmake ${CMAKE_CURRENT_BINARY_DIR}/xdrdump2.input)

foreach(
  name
  os
  sizeof
  xdrdump
  xdrdump2
  boostdump
  boostdump2
  boostdump3
  boostdump4)
  add_executable(${name} ${name}.C)

  target_link_libraries(${name} PRIVATE ALPS::alps)
  alps_add_test(${name})
  set_property(TEST ${name} PROPERTY LABELS osiris)
endforeach()
