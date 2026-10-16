// Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/util.hpp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk in this tree (Adrian Feiguin's dmtk).
// Adaptation: "dmtk/" include paths; dmtk::half_integer and dmtk/support.hpp
// replace the ALPS model and utility headers; the dense products call the
// platform BLAS (dmtk/blas.hpp) again, as the legacy engine did, instead of
// the rework's Eigen views.

#ifndef __DMTK_LA_UTIL_H__
#define __DMTK_LA_UTIL_H__

#include <algorithm>
#include <complex>
#include <cstddef>
#include "dmtk/conj.hpp"
#include "dmtk/blas.hpp"

namespace dmtk
{

using namespace std;


namespace detail {

// gemm-like operation: res := coefa * op1(a) * op2(b) + coefb * res,
// where op{1,2} is selected by the Fortran-style 'N'/'T'/'C' char, through
// the platform BLAS (dmtk/blas.hpp), as the legacy engine did. The alpsx
// rework had routed this through Eigen Map views, materialising op(a),
// op(b) and the product in temporaries; the platform BLAS is several times
// faster on the engine's block sizes and needs no copies. Empty extents are
// handled here because BLAS rejects a leading dimension of zero.
template <class T>
inline void gemm_op(char c1, char c2, int m, int n, int k, T coefa,
                    const T* a, int lda, const T* b, int ldb, T coefb,
                    T* res, int ldc) {
    if (m <= 0 || n <= 0) return;
    if (k <= 0) {
        for (int j = 0; j < n; ++j) {
            T* column = res + static_cast<std::ptrdiff_t>(j) * ldc;
            if (coefb == T(0))
                std::fill_n(column, m, T(0));
            else
                for (int i = 0; i < m; ++i) column[i] *= coefb;
        }
        return;
    }
    blas::gemm(c1, c2, m, n, k, coefa, a, blas::leading(lda), b, blas::leading(ldb),
               coefb, res, blas::leading(ldc));
}

// gemv-like operation: res := coefa * op(a) * v + coefb * res, where a is
// m × n with leading dimension lda.
template <class T>
inline void gemv_op(char c, int m, int n, T coefa,
                    const T* a, int lda, const T* v, T coefb, T* res) {
    int const out = (c == 'N') ? m : n;
    if (out <= 0) return;
    if (m <= 0 || n <= 0) {
        if (coefb == T(0))
            std::fill_n(res, out, T(0));
        else
            for (int i = 0; i < out; ++i) res[i] *= coefb;
        return;
    }
    blas::gemv(c, m, n, coefa, a, blas::leading(lda), v, coefb, res);
}

} // namespace detail

void
matrix_matrix_product(const char &c1, const char &c2,
                      const dmtk::Matrix<double> &a,
                      const dmtk::Matrix<double> &b,
                      dmtk::Matrix<double> &res,
                      const double &coefa = 1.f,
                      const double &coefb = 0.f)
{
  int m = c1 == 'N' ? a.rows() : a.cols();
  int n = c2 == 'N' ? b.cols() : b.rows();
  int k = c1 == 'N' ? a.cols() : a.rows();

  detail::gemm_op<double>(c1, c2, m, n, k, coefa,
                          a.array(), a.rows(),
                          b.array(), b.rows(),
                          coefb, res.array(), res.rows());
}

void
matrix_matrix_product(const char &c1, const char &c2,
                      const dmtk::Matrix<complex<double> >& a,
                      const dmtk::Matrix<complex<double> >& b,
                      dmtk::Matrix<complex<double> >& res,
                      const complex<double> &coefa = 1.f,
                      const complex<double> &coefb = 0.f)
{
  int m = c1 == 'N' ? a.rows() : a.cols();
  int n = c2 == 'N' ? b.cols() : b.rows();
  int k = c1 == 'N' ? a.cols() : a.rows();
  detail::gemm_op<complex<double> >(c1, c2, m, n, k, coefa,
                                    a.array(), a.rows(),
                                    b.array(), b.rows(),
                                    coefb, res.array(), res.rows());
}

void
matrix_matrix_product(const char &c1, const char &c2,
                      const int &m,
                      const int &n,
                      const int &k,
                      const double *a,
                      const int &lda,
                      const double *b,
                      const int &ldb,
                      double *res,
                      const int &ldc,
                      const double &coefa = 1.f,
                      const double &coefb = 0.f)
{
  detail::gemm_op<double>(c1, c2, m, n, k, coefa, a, lda, b, ldb,
                          coefb, res, ldc);
}

void
matrix_matrix_product(const char &c1, const char &c2,
                      const int &m,
                      const int &n,
                      const int &k,
                      const complex<double>* a,
                      const int &lda,
                      const complex<double>* b,
                      const int &ldb,
                      complex<double>* res,
                      const int &ldc,
                      const complex<double> &coefa = 1.f,
                      const complex<double> &coefb = 0.f)
{
  detail::gemm_op<complex<double> >(c1, c2, m, n, k, coefa, a, lda, b, ldb,
                                    coefb, res, ldc);
}


void
matrix_vector_product(const char& c,
                      const dmtk::Matrix<double> &a,
                      const dmtk::Vector<double> &v,
                      dmtk::Vector<double> &res,
                      const double &coefa = 1.f,
                      const double &coefb = 0.f)
{
  detail::gemv_op<double>(c, a.rows(), a.cols(), coefa,
                          a.array(), a.rows(),
                          v.array(), coefb, res.array());
}

void
matrix_vector_product(const char& c,
                      const dmtk::Matrix<complex<double> >&a,
                      const dmtk::Vector<complex<double> >&v,
                      dmtk::Vector<complex<double> >&res,
                      const complex<double> &coefa = 1.f,
                      const complex<double> &coefb = 0.f)
{
  detail::gemv_op<complex<double> >(c, a.rows(), a.cols(), coefa,
                                    a.array(), a.rows(),
                                    v.array(), coefb, res.array());
}

void
matrix_vector_product(const char& c,
                      const int& m,
                      const int& n,
                      const double *a,
                      const double *v,
                      double *res,
                      const double &coefa = 1.f,
                      const double &coefb = 0.f)
{
  detail::gemv_op<double>(c, m, n, coefa, a, m, v, coefb, res);
}

void
matrix_vector_product(const char& c,
                      const int& m,
                      const int& n,
                      const complex<double> *a,
                      const complex<double> *v,
                      complex<double> *res,
                      const complex<double> &coefa = 1.f,
                      const complex<double> &coefb = 0.f)
{
  detail::gemv_op<complex<double> >(c, m, n, coefa, a, m, v, coefb, res);
}

} // namespace dmtk

#endif // __DMTK_LA_UTIL_H__
