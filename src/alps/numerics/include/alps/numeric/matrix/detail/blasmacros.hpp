/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_MATRIX_BLASMACROS_HPP
#define ALPS_MATRIX_BLASMACROS_HPP

// provide overloads for types where blas can be used        

namespace alps {
    namespace numeric {

    #define ALPS_IMPLEMENT_FOR_REAL_BLAS_TYPES(F) F(float) F(double)

    #define ALPS_IMPLEMENT_FOR_COMPLEX_BLAS_TYPES(F) \
    F(std::complex<float>) \
    F(std::complex<double>)

    #define ALPS_IMPLEMENT_FOR_ALL_BLAS_TYPES(F) \
    ALPS_IMPLEMENT_FOR_REAL_BLAS_TYPES(F) \
    ALPS_IMPLEMENT_FOR_COMPLEX_BLAS_TYPES(F)
    } // namespave numeric
} // namespace alps

#endif // ALPS_MATRIX_BLASMACROS_HPP
