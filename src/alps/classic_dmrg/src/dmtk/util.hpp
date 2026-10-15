// Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/util.hpp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk in this tree (Adrian Feiguin's dmtk).
// Adaptation: "dmtk/" include paths; dmtk::half_integer and dmtk/support.hpp
// replace the ALPS model and utility headers.

#ifndef __DMTK_LA_UTIL_H__
#define __DMTK_LA_UTIL_H__

#include <complex>
#include "dmtk/conj.hpp"

#include <Eigen/Core>

namespace dmtk
{

using namespace std;


namespace detail {

// gemm-like operation: res := coefa * op1(a) * op2(b) + coefb * res,
// where op{1,2} is selected by the Fortran-style 'N'/'T'/'C' char.
//
// Replaces the old cblas_gemm dispatch with Eigen Map views over the
// raw column-major arrays. The product is materialised into a
// temporary EMat (cheap for the small DMRG blocks involved) so that
// aliasing between res and a/b stays safe.
template <class T>
inline void gemm_op(char c1, char c2, int m, int n, int k, T coefa,
                    const T* a, int lda, const T* b, int ldb, T coefb,
                    T* res, int ldc) {
    using EMat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using ConstMap = Eigen::Map<const EMat, 0, Eigen::OuterStride<>>;
    using Map = Eigen::Map<EMat, 0, Eigen::OuterStride<>>;

    // a is (lda × k') with k' = (c1 == 'N' ? k : m) physical cols.
    int a_rows = (c1 == 'N') ? m : k;
    int a_cols = (c1 == 'N') ? k : m;
    int b_rows = (c2 == 'N') ? k : n;
    int b_cols = (c2 == 'N') ? n : k;

    ConstMap A(a, a_rows, a_cols, Eigen::OuterStride<>(lda));
    ConstMap B(b, b_rows, b_cols, Eigen::OuterStride<>(ldb));
    Map      C(res, m, n, Eigen::OuterStride<>(ldc));

    auto op_a = [&]() -> EMat {
        if (c1 == 'C') return A.adjoint();
        if (c1 == 'T') return A.transpose();
        return A;
    };
    auto op_b = [&]() -> EMat {
        if (c2 == 'C') return B.adjoint();
        if (c2 == 'T') return B.transpose();
        return B;
    };

    EMat const product = op_a() * op_b();
    if (coefb == T(0))
        C = coefa * product;
    else
        C = coefa * product + coefb * C;
}

// gemv-like operation: res := coefa * op(a) * v + coefb * res.
template <class T>
inline void gemv_op(char c, int m, int n, T coefa,
                    const T* a, int lda, const T* v, T coefb, T* res) {
    using EMat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using EVec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    using ConstMatMap = Eigen::Map<const EMat, 0, Eigen::OuterStride<>>;
    using ConstVecMap = Eigen::Map<const EVec>;
    using VecMap = Eigen::Map<EVec>;

    ConstMatMap A(a, m, n, Eigen::OuterStride<>(lda));

    int vec_in_size  = (c == 'N') ? n : m;
    int vec_out_size = (c == 'N') ? m : n;
    ConstVecMap V(v, vec_in_size);
    VecMap      R(res, vec_out_size);

    EVec const product = [&]() -> EVec {
        if (c == 'C') return A.adjoint() * V;
        if (c == 'T') return A.transpose() * V;
        return A * V;
    }();
    if (coefb == T(0))
        R = coefa * product;
    else
        R = coefa * product + coefb * R;
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
