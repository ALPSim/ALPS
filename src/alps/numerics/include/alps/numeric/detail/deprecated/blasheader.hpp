/*****************************************************************************
 *
 * ALPS DMFT Project - BLAS Compatibility headers
 *  BLAS headers for accessing BLAS from C++.
 *
 * Copyright (C) 2005 - 2009 by 
 *                              Emanuel Gull <gull@phys.columbia.edu>,
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef BOOST_NUMERIC_DETAIL_BLASHEADER_HPP
#define BOOST_NUMERIC_DETAIL_BLASHEADER_HPP

#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/lapack.hpp>

#include<complex>

namespace lapack{
extern "C" void vvexp(double * /* y */,const double * /* x */,const int * /* n */);
}
namespace acml{
extern "C" void vrda_exp(const int, double *, const double *);
}
namespace mkl{
extern "C" void vdExp(const int, const double *, double *);
extern "C" void vfExp(const int, const float *, float *);
extern "C" void vzmul_(const int*, const std::complex<double>*, const std::complex<double> *, std::complex<double> *);
}

#endif
