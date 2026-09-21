/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef __DMCONJ_H__
#define __DMCONJ_H__

#include <complex>

// Borrowed from MTL (Matrix Template Library)

namespace dmtk {

// Bring the std::complex<T> overloads into this namespace so that
// dmtk::conj/real/imag also work for complex arguments (via ADL or
// explicit dmtk:: qualification), not just the dummy real-number
// overloads defined below.
using std::conj;
using std::real;
using std::imag;

// dummy conj function for real numbers
inline double conj(double a) {
  return a;
}
inline float conj(float a) {
  return a;
}
inline int conj(int a) {
  return a;
}
inline bool conj(bool a) {
  return a;
}

// dummy real and imag function for real numbers
inline double real(double a) {
  return a;
}
inline double imag(double) {
  return 0.0;
}

inline float real(float a) {
  return a;
}
inline float imag(float) {
  return 0.0;
}

} 

#endif /* __DMCONJ_H__ */
