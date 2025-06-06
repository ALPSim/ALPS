/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef __DMTK_LAPACK_INTERFACE_H__
#define __DMTK_LAPACK_INTERFACE_H__

#include <complex>
#include <alps/config.h> // needed to set up correct bindings
#include <boost/numeric/bindings/detail/config/fortran.hpp>
#include <boost/numeric/bindings/traits/type.h>

extern "C" {

void FORTRAN_ID(zheev)(
    const char &jobz,        // (input)
    const char &uplo,        // (input)
    const int &n,            // (input)
    std::complex<double> *a,    // a[n][lda] (input/output)
    const int &lda,            // (input)
    double *w,            // w[n] (output)
    std::complex<double> *work,    // work[lwork] (workspace/output)
    const int &lwork,        // (input)
    double *rwork,            // rwork[max(1, 3*n-2)] (workspace)
    int &info            // (output)
    );

void FORTRAN_ID(dsyev)(
    const char &jobz,        // (input)
    const char &uplo,        // (input)
    const int &n,            // (input)
    double *a,            // a[n][lda] (input/output)
    const int &lda,            // (input)
    double *w,            // w[n] (output)
    double *work,            // work[lwork] (workspace/output)
    const int &lwork,        // (input)
    int &info            // (output)
    );

void FORTRAN_ID(zhpev)(
    const char &jobz,        // (input)
    const char &uplo,        // (input)
    const int &n,            // (input)
    std::complex<double> *ap,    // ap[n*(n+1)/2] (input/output)
    double *w,            // w[n] (output)
    std::complex<double> *z,    // z[n][ldz] (output)
    const int &ldz,            // (input)
    std::complex<double> *work,    // work[max(1, 2*n-1)] (workspace)
    double *rwork,            // rwork[max(1, 3*n-2)] (workspace)
    int &info            // (output)
    );

void FORTRAN_ID(dstevd)(
    const char &jobz,        // (input)
    const int &n,            // (input)
    double *d,            // d[n] (input/output)
    double *e,            // e[n] (input/output)
    double *z,            // z[n][ldz] (output)
    const int &ldz,            // (input)
    double *work,            // work[?] (workspace/output)
    const int &lwork,        // (input)
    int *iwork,            // iwork[liwork] (workspace/output)
    const int &liwork,        // (input)
    int &info            // (output)
    );

double FORTRAN_ID(ddot)(
            const int&,
            const double*, const int&,
            const double*, const int&);

dcomplex_t FORTRAN_ID(zdotc)(
            const int&,
            const std::complex<double>*, const int&,
            const std::complex<double>*, const int&);

void FORTRAN_ID(dcopy)(const int &,
            const double *, const int &,
            double *, const int &);

void FORTRAN_ID(zcopy)(const int &,
            const std::complex<double> *, const int &,
            std::complex<double> *, const int &);

void FORTRAN_ID(dgemv)(const char&,
            const int&, const int&,
            const double&,
            const double*, const int&,
            const double*, const int&,
            const double&,
            double*, const int&);

void FORTRAN_ID(zgemv)(const char&,
            const int&, const int&,
            const std::complex<double>&,
            const std::complex<double>*, const int&,
            const std::complex<double>*, const int&,
            const std::complex<double>&,
            std::complex<double>*, const int&);

void FORTRAN_ID(dgemm)(const char&, const char&,
            const int&, const int&, const int&,
            const double&,
            const double*, const int&,
            const double*, const int& ,
            const double&,
            double *, const int&);

void FORTRAN_ID(zgemm)(const char&, const char&,
            const int&, const int&, const int&,
            const std::complex<double>&,
            const std::complex<double>*, const int&,
            const std::complex<double>*, const int&,
            const std::complex<double>&,
            std::complex<double>*, const int&);

void FORTRAN_ID(daxpy)(const int &n, const double& da,
            const double* dx, const int& incx,
            const double* dy, const int& incy);
}

#endif // __DMTK_LAPACK_INTERFACE_H__
