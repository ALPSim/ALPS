// Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/globals.hpp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk in this tree (Adrian Feiguin's dmtk).
// Adaptation: "dmtk/" include paths; dmtk::half_integer and dmtk/support.hpp
// replace the ALPS model and utility headers.

#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#ifndef MIN_VECTOR_SIZE
#define MIN_VECTOR_SIZE 1000
#endif
#ifndef MIN_MATRIX_SIZE
#define MIN_MATRIX_SIZE 100
#endif

#include "dmtk/vector.hpp"
#include "dmtk/matrix.hpp"

namespace dmtk
{

template<class T>
class DMTKglobals
{
  public:
    Matrix<T> m1;
    Matrix<T> m2;
    Matrix<T> m3;
    Matrix<T> m4;
    Vector<T> v1;
    Vector<T> v2;
    Vector<T> v3;
    Vector<T> v4;

    DMTKglobals() { m1 = m2 = m3 = m4 = Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE); v1 = v2 = v3 = v4 = Vector<T>(MIN_VECTOR_SIZE); }

    void swap(DMTKglobals& other) noexcept
    {
      m1.swap(other.m1); m2.swap(other.m2);
      m3.swap(other.m3); m4.swap(other.m4);
      v1.swap(other.v1); v2.swap(other.v2);
      v3.swap(other.v3); v4.swap(other.v4);
    }
};

static DMTKglobals<double> globals_double;
static DMTKglobals<complex<double> > globals_complex;

DMTKglobals<double> &
get_globals(const double&)
{
  return globals_double;
}

DMTKglobals<complex<double> >&
get_globals(const complex<double>&)
{
  return globals_complex;
}

} //namespace dmtk

#endif // __GLOBALS_H__
