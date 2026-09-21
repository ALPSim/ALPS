#  Copyright Synge Todo 2012.
#   SPDX-License-Identifier: MIT

foreach(name vmusage bitops)
  add_executable(${name} ${name}.cpp)

  target_link_libraries(${name} PRIVATE ALPS::alps)
  alps_add_test(${name})
  set_property(TEST ${name} PROPERTY LABELS utility)
endforeach()
