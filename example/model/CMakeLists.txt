#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

#
# examples of using the model library for numeric and symbolic computations
# 

if (ALPS_BUILD_EXAMPLES)
  set (EXAMPLES print_numeric print_symbolic print_numeric2 print_symbolic2 print_numeric3 print_symbolic3)
  set (INPUT parm_numeric parm_symbolic)
  
  foreach(file ${INPUT})
    configure_file(${file}.in ${CMAKE_CURRENT_BINARY_DIR}/${file})
  endforeach(file)
     
  foreach(program ${EXAMPLES})
    add_executable(${program} ${program}.cpp)
    target_link_libraries(${program} alps)
  endforeach(program)
endif (ALPS_BUILD_EXAMPLES)

