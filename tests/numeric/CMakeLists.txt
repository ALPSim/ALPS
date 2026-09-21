#  Copyright Matthias Troyer and Synge Todo 2009 - 2010.
#   SPDX-License-Identifier: MIT

foreach(
  name
  matrix_unit_tests
  matrix_column_view
  matrix_transpose_view
  matrix_hdf5
  matrix_deprecated_hdf5_format_test
  matrix_ublas_sparse_functions
  accumulate_if
  vector_functions
  vector_valarray_conversion
  matrix_algorithms
  matrix_kron
  real_tests)
  add_executable(${name} ${name}.C)

  target_link_libraries(${name} PRIVATE ALPS::alps)
  target_compile_definitions(${name} PRIVATE ALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
  alps_add_test(${name})
  set_property(TEST ${name} PROPERTY LABELS numeric)
endforeach()

# This test reads a reference .h5 file from the source tree. That path is a
# property of this build, not of the installed library, so it is a private
# definition on the one target that needs it rather than a macro in the
# installed alps/version.h.
target_compile_definitions(matrix_deprecated_hdf5_format_test
                           PRIVATE ALPS_SRCDIR="${PROJECT_SOURCE_DIR}")
