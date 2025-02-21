#  Copyright Matthias Troyer and Synge Todo 2009 - 2012.
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
if(NOT ALPS_PYTHON_WHEEL AND NOT ALPS_BUILD_LIBS_ONLY)
add_subdirectory(accumulator)
add_subdirectory(alea)
add_subdirectory(fixed_capacity)
add_subdirectory(graph)
add_subdirectory(hdf5)
add_subdirectory(lattice)
add_subdirectory(model)
add_subdirectory(ngs)
add_subdirectory(numeric)
add_subdirectory(osiris)
add_subdirectory(parameter)
add_subdirectory(parapack)
add_subdirectory(parser)
add_subdirectory(pyalps)
add_subdirectory(random)
add_subdirectory(utility)
endif()