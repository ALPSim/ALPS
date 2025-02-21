#  Copyright Matthias Troyer and Synge Todo 2009 - 2010.
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

if(LAPACK_FOUND)
  include_directories(${PROJECT_BINARY_DIR}/src)
  include_directories(${PROJECT_SOURCE_DIR}/src)
  include_directories(${Boost_ROOT_DIR})
  add_definitions(-DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
  FOREACH (name matrix_unit_tests matrix_column_view matrix_transpose_view matrix_hdf5 matrix_deprecated_hdf5_format_test matrix_ublas_sparse_functions accumulate_if vector_functions vector_valarray_conversion matrix_algorithms real_tests)
    add_executable(${name} ${name}.C)
    add_dependencies(${name} alps)
    target_link_libraries(${name} alps ${BLAS_LIBRARY} ${LAPACK_LIBRARY})
    set_property(TARGET ${name} PROPERTY LABELS numeric)
    add_alps_test(${name})
    set_property(TEST ${name} PROPERTY LABELS numeric)
  ENDFOREACH(name)
endif(LAPACK_FOUND)

