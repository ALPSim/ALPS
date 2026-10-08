/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef ALPS_NUMERIC_MATRIX_GEMM_HPP
#define ALPS_NUMERIC_MATRIX_GEMM_HPP
#include <boost/numeric/bindings/blas/level3/gemm.hpp>
#include <alps/numeric/matrix/detail/debug_output.hpp>
#include <alps/numeric/matrix/is_blas_dispatchable.hpp>
#include <algorithm>
#include <cassert>


namespace alps {
namespace numeric {
    //
    // Default matrix matrix multiplication implementations
    // May be overlaoded for special Matrix types
    //

    template<typename MatrixA, typename MatrixB, typename MatrixC>
    void gemm(MatrixA const& a, MatrixB const& b, MatrixC& c, boost::mpl::false_)
    {
        using std::fill;
        assert( num_cols(a) == num_rows(b) );
        assert( num_rows(c) == num_rows(a) );
        assert( num_cols(c) == num_cols(b) );
        // Simple matrix matrix multiplication
        for(std::size_t j=0; j < num_cols(b); ++j)
        {
            std::pair<typename MatrixC::col_element_iterator,typename MatrixC::col_element_iterator> range = col(c,j);
            fill(range.first, range.second, typename MatrixC::value_type(0));
            for(std::size_t k=0; k < num_cols(a); ++k)
                for(std::size_t i=0; i < num_rows(a); ++i)
                    c(i,j) += a(i,k) * b(k,j);
        }
    }


    template<typename MatrixA, typename MatrixB, typename MatrixC>
    void gemm(MatrixA const& a, MatrixB const& b, MatrixC& c, boost::mpl::true_)
    {
        ALPS_NUMERIC_MATRIX_DEBUG_OUTPUT( "using blas gemm for " << typeid(a).name() << " " << typeid(b).name() << " -> " << typeid(c).name() );
        typedef typename MatrixA::value_type value_type;
        boost::numeric::bindings::blas::gemm(value_type(1),a,b,value_type(0),c);
    }

    // The classic gemm - as known from Fortran - writing the result to argument c
    template<typename MatrixA, typename MatrixB, typename MatrixC>
    void gemm(MatrixA const & a, MatrixB const & b, MatrixC & c)
    {
        assert( num_cols(a) == num_rows(b) );
        assert( num_rows(c) == num_rows(a) );
        assert( num_cols(c) == num_cols(b) );
        // TODO this check should also involve MatrixC
        gemm(a,b,c,is_blas_dispatchable<MatrixA,MatrixB>());
    }

} // end namespace numeric
} // end namespace alps

#endif // ALPS_NUMERIC_MATRIX_GEMM_HPP
