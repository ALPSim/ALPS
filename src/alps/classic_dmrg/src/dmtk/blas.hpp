// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Written for this tree: the four BLAS level-2/3 routines the dmtk engine
// uses for its dense block products, through the Fortran LP64 ABI that the
// ALPS build selects (lowercase, trailing underscore, 32-bit integers; see
// cmake/ALPSNumerics.cmake), as ALPS::ed declares its LAPACK routines. The
// legacy engine (src/apps/dmrg/dmrg/dmtk/util.h on master) reached the same
// routines through the Boost bindings' FORTRAN_ID; this header is private to
// ALPS::classic_dmrg, so consumers of the installed SDK need no BLAS headers.

#ifndef __DMTK_BLAS_H__
#define __DMTK_BLAS_H__

#include <algorithm>
#include <complex>

extern "C" {
void dgemm_(char const* transa, char const* transb, int const* m, int const* n, int const* k,
            double const* alpha, double const* a, int const* lda, double const* b, int const* ldb,
            double const* beta, double* c, int const* ldc);
void zgemm_(char const* transa, char const* transb, int const* m, int const* n, int const* k,
            std::complex<double> const* alpha, std::complex<double> const* a, int const* lda,
            std::complex<double> const* b, int const* ldb, std::complex<double> const* beta,
            std::complex<double>* c, int const* ldc);
void dgemv_(char const* trans, int const* m, int const* n, double const* alpha, double const* a,
            int const* lda, double const* x, int const* incx, double const* beta, double* y,
            int const* incy);
void zgemv_(char const* trans, int const* m, int const* n, std::complex<double> const* alpha,
            std::complex<double> const* a, int const* lda, std::complex<double> const* x,
            int const* incx, std::complex<double> const* beta, std::complex<double>* y,
            int const* incy);
}

namespace dmtk {
namespace blas {

inline void gemm(char transa, char transb, int m, int n, int k, double alpha, double const* a,
                 int lda, double const* b, int ldb, double beta, double* c, int ldc) {
    dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
inline void gemm(char transa, char transb, int m, int n, int k, std::complex<double> alpha,
                 std::complex<double> const* a, int lda, std::complex<double> const* b, int ldb,
                 std::complex<double> beta, std::complex<double>* c, int ldc) {
    zgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
inline void gemv(char trans, int m, int n, double alpha, double const* a, int lda, double const* x,
                 double beta, double* y) {
    int const one = 1;
    dgemv_(&trans, &m, &n, &alpha, a, &lda, x, &one, &beta, y, &one);
}
inline void gemv(char trans, int m, int n, std::complex<double> alpha, std::complex<double> const* a,
                 int lda, std::complex<double> const* x, std::complex<double> beta,
                 std::complex<double>* y) {
    int const one = 1;
    zgemv_(&trans, &m, &n, &alpha, a, &lda, x, &one, &beta, y, &one);
}

// A leading dimension below one is an argument error for BLAS (xerbla), which
// the engine's empty quantum-number blocks would otherwise trigger.
inline int leading(int ld) { return std::max(ld, 1); }

} // namespace blas
} // namespace dmtk

#endif // __DMTK_BLAS_H__
