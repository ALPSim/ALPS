#  Copyright Synge Todo 2010 - 2015.
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

add_executable(simplemc main.C evaluator.C ising.C xy.C heisenberg.C)
target_link_libraries(simplemc alps)
if(ALPS_PYTHON_WHEEL)
  if(APPLE)
    set_target_properties(simplemc PROPERTIES INSTALL_RPATH "@loader_path/../lib" )
  else(APPLE)
    set_target_properties(simplemc PROPERTIES INSTALL_RPATH "$ORIGIN/../lib" )
  endif(APPLE)
  install(TARGETS simplemc RUNTIME DESTINATION pyalps/bin COMPONENT applications)
else()
  install(TARGETS simplemc RUNTIME DESTINATION bin COMPONENT applications)
endif()
enable_testing()
add_alps_test(simplemc_ising simplemc ising ising)
add_alps_test(simplemc_xy simplemc xy xy)
add_alps_test(simplemc_heisenberg simplemc heisenberg heisenberg)
