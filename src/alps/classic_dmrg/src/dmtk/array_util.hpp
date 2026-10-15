// Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/array_util.hpp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk in this tree (Adrian Feiguin's dmtk).
// Adaptation: "dmtk/" include paths; dmtk::half_integer and dmtk/support.hpp
// replace the ALPS model and utility headers.

#ifndef __DMTK_ARRAY_UTIL_H__
#define __DMTK_ARRAY_UTIL_H__

#include <algorithm>
#include <complex>
#include "dmtk/conj.hpp"
#include <Eigen/Core>

namespace dmtk
{

using namespace std;


inline double
dot_product(const int& n,
            const double *v1, const int& stride1,
            const double *v2, const int& stride2)
{
  using EVec = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using StridedMap = Eigen::Map<const EVec, 0, Eigen::InnerStride<>>;
  StridedMap a(v1, n, Eigen::InnerStride<>(stride1));
  StridedMap b(v2, n, Eigen::InnerStride<>(stride2));
  return a.dot(b);
}

inline complex<double>
dot_product(const int& n,
            const complex<double> *v1, const int& stride1,
            const complex<double> *v2, const int& stride2)
{
  complex<double> sum(0.f,0.f);
  const complex<double> *pi = v1;
  const complex<double> *pj = v2;

//  sum = zdotc_(n, v1, stride1, v2, stride2);
//  return sum;

  if(stride1 == 1 && stride2 == 1)
    for(int k = 0; k < n; k++)
      sum += std::conj(*pi++)*(*pj++);
  else
    for(int k = 0; k < n; k++){
      sum += std::conj(*pi)*(*pj);
      pi += stride1;
      pj += stride2;
    }

  return sum;
}


// Pointer-to-pointer copy. Earlier versions hand-unrolled by 7 for
// 1990s-era compilers and added cblas_{d,z}copy specialisations for
// double and complex<double>; modern toolchains compile std::copy_n
// to memcpy for trivially-copyable T and auto-vectorise the rest, so
// the explicit unroll and the BLAS calls just obscured the intent.
template<class T>
inline void
array_copy(int n, const T* in, T* out)
{
  std::copy_n(in, n, out);
}

// Indexed copy for slice_iter / cslice_iter / gslice_iter — these
// proxy types support operator[] but are not random-access iterators
// in the std::iterator_traits sense, so std::copy_n cannot consume them.
template<class T>
inline void
array_copy(int n, T& in, T& out)
{
  for (int i = 0; i < n; ++i) out[i] = in[i];
}

template<class T1, class T2>
inline void
array_copy2(int n, const T1& in, T2& out)
{
  for (int i = 0; i < n; ++i) out[i] = in[i];
}

template<class T1, class T2>
inline void
array_copy2(int n, T1& in, T2& out)
{
  for (int i = 0; i < n; ++i) out[i] = in[i];
}

inline double quickran(long & idum)
{
    const int im = 134456;
    const int ia = 8121;
    const int ic = 28411;
    const double scale = 1.0 / im;
    idum = (idum*ia+ic)%im;
    return double(idum) * scale;
}

} // namespace dmtk

#endif // __DMTK_ARRAY_UTIL_H__
