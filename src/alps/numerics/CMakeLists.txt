add_library(alps_numerics INTERFACE)
add_library(ALPS::numerics ALIAS alps_numerics)
set_target_properties(alps_numerics PROPERTIES EXPORT_NAME numerics)
target_link_libraries(alps_numerics INTERFACE ALPS::configuration ALPS::containers ALPS::utilities
  ${ALPS_BLAS_TARGET} ${ALPS_LAPACK_TARGET})
# Numerics source ownership.
target_sources(alps_numerics INTERFACE
  FILE_SET numerics_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/multi_array.hpp
    include/alps/multi_array/functions.hpp
    include/alps/multi_array/operators.hpp
    include/alps/functional.h
    include/alps/utility/numeric_cast.hpp
    include/alps/ngs/numeric.hpp
    include/alps/ngs/numeric/array.hpp
    include/alps/ngs/numeric/detail.hpp
    include/alps/ngs/numeric/inf.hpp
    include/alps/ngs/numeric/multi_array.hpp
    include/alps/ngs/numeric/vector.hpp
    include/alps/numeric/abs2.hpp
    include/alps/numeric/accumulate_if.hpp
    include/alps/numeric/binomial.hpp
    include/alps/numeric/checked_divide.hpp
    include/alps/numeric/conj.hpp
    include/alps/numeric/diagonal_matrix.hpp
    include/alps/numeric/double2int.hpp
    include/alps/numeric/fourier.hpp
    include/alps/numeric/functional.hpp
    include/alps/numeric/imag.hpp
    include/alps/numeric/is_equal.hpp
    include/alps/numeric/is_negative.hpp
    include/alps/numeric/is_nonzero.hpp
    include/alps/numeric/is_positive.hpp
    include/alps/numeric/is_zero.hpp
    include/alps/numeric/isinf.hpp
    include/alps/numeric/isnan.hpp
    include/alps/numeric/matrix.hpp
    include/alps/numeric/matrix/algorithms.hpp
    include/alps/numeric/matrix/column_view.hpp
    include/alps/numeric/matrix/conj.hpp
    include/alps/numeric/matrix/detail/auto_deduce_multiply_return_type.hpp
    include/alps/numeric/matrix/detail/auto_deduce_plus_return_type.hpp
    include/alps/numeric/matrix/detail/blasmacros.hpp
    include/alps/numeric/matrix/detail/column_view_adaptor.hpp
    include/alps/numeric/matrix/detail/debug_output.hpp
    include/alps/numeric/matrix/detail/matrix_adaptor.hpp
    include/alps/numeric/matrix/detail/print_matrix.hpp
    include/alps/numeric/matrix/detail/print_vector.hpp
    include/alps/numeric/matrix/detail/transpose_view_adaptor.hpp
    include/alps/numeric/matrix/detail/vector_adaptor.hpp
    include/alps/numeric/matrix/entity.hpp
    include/alps/numeric/matrix/exchange_value_type.hpp
    include/alps/numeric/matrix/gemm.hpp
    include/alps/numeric/matrix/gemv.hpp
    include/alps/numeric/matrix/is_blas_dispatchable.hpp
    include/alps/numeric/matrix/matrix.hpp
    include/alps/numeric/matrix/matrix.ipp
    include/alps/numeric/matrix/matrix_concept_archetype.hpp
    include/alps/numeric/matrix/matrix_concept_check.hpp
    include/alps/numeric/matrix/matrix_element_iterator.hpp
    include/alps/numeric/matrix/matrix_interface.hpp
    include/alps/numeric/matrix/matrix_traits.hpp
    include/alps/numeric/matrix/operators/multiply.hpp
    include/alps/numeric/matrix/operators/multiply_matrix.hpp
    include/alps/numeric/matrix/operators/multiply_scalar.hpp
    include/alps/numeric/matrix/operators/op_assign.hpp
    include/alps/numeric/matrix/operators/op_assign_matrix.hpp
    include/alps/numeric/matrix/operators/op_assign_vector.hpp
    include/alps/numeric/matrix/operators/plus_minus.hpp
    include/alps/numeric/matrix/resizable_matrix_concept_check.hpp
    include/alps/numeric/matrix/resizable_matrix_interface.hpp
    include/alps/numeric/matrix/scalar_product.hpp
    include/alps/numeric/matrix/strided_iterator.hpp
    include/alps/numeric/matrix/transpose.hpp
    include/alps/numeric/matrix/transpose_view.hpp
    include/alps/numeric/matrix/ublas_sparse_functions.hpp
    include/alps/numeric/matrix/vector.hpp
    include/alps/numeric/matrix/vector_interface.hpp
    include/alps/numeric/matrix_as_vector.hpp
    include/alps/numeric/outer_product.hpp
    include/alps/numeric/polynomial.hpp
    include/alps/numeric/real.hpp
    include/alps/numeric/regression.hpp
    include/alps/numeric/round.hpp
    include/alps/numeric/scalar_product.hpp
    include/alps/numeric/sequence_comparisons.hpp
    include/alps/numeric/set_negative_0.hpp
    include/alps/numeric/special_functions.hpp
    include/alps/numeric/update_minmax.hpp
    include/alps/numeric/valarray_functions.hpp
    include/alps/numeric/vector_functions.hpp
    include/alps/numeric/vector_valarray_conversion.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
