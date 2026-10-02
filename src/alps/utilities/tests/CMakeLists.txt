#  Copyright Synge Todo 2012.
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

foreach(name vmusage bitops)
  add_executable(${name} ${name}.cpp)

  target_link_libraries(${name} PRIVATE ALPS::utilities)
  alps_add_test(${name})
  set_property(TEST ${name} PROPERTY LABELS utility)
endforeach()

add_executable(ngs_hash ngs_hash.cpp)
target_link_libraries(ngs_hash PRIVATE ALPS::utilities)
alps_add_test(ngs_hash)

# This consumer has no aggregate or numerical include roots.
add_executable(foundation_headers foundation_headers.cpp)
target_link_libraries(foundation_headers PRIVATE ALPS::utilities)
alps_add_test(foundation_headers)
set_property(TEST foundation_headers PROPERTY LABELS utility architecture)
