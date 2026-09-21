#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

#
# examples of using the model library for numeric and symbolic computations
#

set(EXAMPLES print_numeric print_symbolic print_numeric2 print_symbolic2 print_numeric3
             print_symbolic3)
set(INPUT parm_numeric parm_symbolic)

foreach(file ${INPUT})
  configure_file(${file}.in ${CMAKE_CURRENT_BINARY_DIR}/${file})
endforeach()

foreach(program ${EXAMPLES})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} PRIVATE ALPS::alps)
endforeach()
