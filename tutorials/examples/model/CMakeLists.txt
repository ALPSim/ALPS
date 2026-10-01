#  Copyright Matthias Troyer 2009 - 2010.
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the “Software”),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#  
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#  
#   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

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

