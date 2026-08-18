/*****************************************************************************
 *
 * ALPS DMFT Project - BLAS Compatibility headers
 *  BLAS headers for accessing BLAS from C++.
 *
 * Copyright (C) 2010 Matthias Troyer <gtroyer@ethz.ch>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <boost/numeric/bindings/blas.hpp>
#include<vector>

// provide overloads for types where blas can be used        

namespace blas{

#define IMPLEMENT_FOR_REAL_BLAS_TYPES(F) F(float) F(double)

#define IMPLEMENT_FOR_COMPLEX_BLAS_TYPES(F) \
F(std::complex<float>) \
F(std::complex<double>)

#define IMPLEMENT_FOR_ALL_BLAS_TYPES(F) \
IMPLEMENT_FOR_REAL_BLAS_TYPES(F) \
IMPLEMENT_FOR_COMPLEX_BLAS_TYPES(F) 
} // namespace
