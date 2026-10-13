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

foreach(
  name
  example1
  example3
  example4
  example5
  example8
  example11
  example13
  example14
  example17
  example18)
  add_executable(model_${name} ${name}.C)

  target_link_libraries(model_${name} PRIVATE ALPS::alps)
  alps_add_gtest(model_${name} COMPONENT model LABELS compatibility)
endforeach()

# Related cases share compilation and linkage, but CTest still discovers and
# executes every GoogleTest case independently.
foreach(name site_basis basis_states bloch_basis_states)
  add_executable(model_${name} ${name}.cpp)
  target_link_libraries(model_${name} PRIVATE ALPS::alps)
  alps_add_gtest(model_${name} COMPONENT model LABELS compatibility)
endforeach()
