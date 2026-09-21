#  Copyright Matthias Troyer 2009 - 2010
#            Lukas Gamper 2012
#   SPDX-License-Identifier: MIT

foreach(name default ordering stream assign not_found external)
  add_executable(param_${name} ${name}.cpp)

  target_link_libraries(param_${name} PRIVATE ALPS::alps)
  alps_add_test(param_${name})
endforeach()
